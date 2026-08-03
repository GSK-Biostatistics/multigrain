// graph_shortcut.cpp
// -----------------------------------------------------------------------------
// Serial and parallel (RcppParallel) implementations of the Bretz et al.
// (2009) shortcut algorithm for graph-based multiple testing.
//
// Both functions produce a rejection matrix only (0/1 entries) as a
// LogicalMatrix, matching the convention of apply_ctp / apply_ctp_parallel.
//
// Serial implementation (graph_shortcut):
//   1. Working graph stored as std::vector<double>.
//   2. Two reusable buffers (cur_G, new_G) ping-pong via std::swap (O(1),
//      exchanges internal pointers; no copy).
//   3. Raw pointer access to REAL(pvals), REAL(G), REAL(w), LOGICAL(h).
//   4. P-value row pre-extracted into a dense local buffer (column-major
//      stride avoided in the hot loop).
//   5. Early termination when all hypotheses are rejected.
//   6. Graph update skipped when only one hypothesis remains (unused work).
//   7. Division 1 / (1 - gir * gri) hoisted out of the inner j loop.
//
// Parallel implementation (graph_shortcut_parallel):
//   Parallelised across trials via RcppParallel::parallelFor. Thread-safety:
//   * pvals, G, init_a are read-only in the Worker.
//   * Each trial writes to a unique row of the output, so there are no
//     write-write races.
//   * Scratch buffers (cur_a, cur_G, new_G, cur_p) are allocated inside
//     operator(), i.e. once per chunk dispatch. TBB's thread-local allocator
//     caches make these allocations effectively free in practice.
// -----------------------------------------------------------------------------

// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <algorithm>

using namespace Rcpp;
using namespace RcppParallel;


// =============================================================================
// Serial implementation
// =============================================================================

//' Graphical shortcut algorithm (Bretz et al. 2009) -- fast serial backend
//'
//' Applies the sequentially rejective graphical test procedure to a matrix of
//' simulated p-values and returns only the rejection matrix (no adjusted
//' alphas or graph history). Matches the behaviour of gMCPLite::graphTest()
//' in the non-parametric, non-entangled, non-upscale case.
//'
//' @param pvals NumericMatrix (N x m) of simulated p-values.
//' @param alpha Scalar significance level.
//' @param w NumericVector (length m) of hypothesis weights.
//' @param G NumericMatrix (m x m) transition matrix. Diagonal must be 0,
//'   row sums must be <= 1.
//' @return LogicalMatrix (N x m) of rejection indicators.
//'
//' @seealso \code{\link{graph_shortcut_parallel}} for the multithreaded variant.
//' @name graph_shortcut
//' @noRd
// [[Rcpp::export]]
LogicalMatrix graph_shortcut(const NumericMatrix& pvals,
                             const double         alpha,
                             const NumericVector& w,
                             const NumericMatrix& G) {
  const int N = pvals.nrow();
  const int m = pvals.ncol();

  if (G.nrow() != m || G.ncol() != m)
    stop("G must be an m x m matrix matching ncol(pvals).");
  if (w.size() != m)
    stop("w must have length equal to ncol(pvals).");

  LogicalMatrix h(N, m);  // zero-initialised by R

  // Raw pointers to underlying column-major storage
  const double* p_ptr = REAL(pvals);
  const double* G_ptr = REAL(G);
  const double* w_ptr = REAL(w);
  int* h_ptr = LOGICAL(h);

  // Pre-allocated working buffers (reused across trials)
  std::vector<double> init_a(m);
  std::vector<double> cur_a(m);
  std::vector<double> cur_G(m * m);
  std::vector<double> new_G(m * m);
  std::vector<double> cur_p(m);

  // Compute initial local alphas once: a[i] = w[i] * alpha
  for (int i = 0; i < m; ++i) {
    init_a[i] = w_ptr[i] * alpha;
  }

  for (int set = 0; set < N; ++set) {
    // Copy initial state for this trial
    std::copy(init_a.begin(), init_a.end(), cur_a.begin());
    std::copy(G_ptr, G_ptr + m * m, cur_G.begin());

    // Extract this trial's p-value row once (column-major -> dense local copy)
    for (int i = 0; i < m; ++i) {
      cur_p[i] = p_ptr[set + i * N];
    }

    int sumrej = 0;

    while (true) {
      // Find first hypothesis meeting rejection condition
      int rej = -1;
      for (int i = 0; i < m; ++i) {
        if (cur_p[i] < cur_a[i]) {
          rej = i;
          h_ptr[set + i * N] = 1;   // TRUE in R's logical representation
          ++sumrej;
          break;
        }
      }
      if (rej == -1) break;      // no further rejection possible
      if (sumrej == m) break;    // all hypotheses rejected

      // Update local alphas: a[i] += a[rej] * G[rej, i]
      for (int i = 0; i < m; ++i) {
        cur_a[i] += cur_a[rej] * cur_G[rej + i * m];
      }

      // Update graph only if >1 hypothesis remains (otherwise unused work)
      if (sumrej < m - 1) {
        for (int i = 0; i < m; ++i) {
          const double gir = cur_G[i + rej * m];  // G[i, rej]
          const double gri = cur_G[rej + i * m];  // G[rej, i]
          const double denom = 1.0 - gir * gri;

          if (denom != 0.0) {
            const double inv_denom = 1.0 / denom;
            for (int j = 0; j < m; ++j) {
              if (i == j) {
                new_G[i + j * m] = 0.0;
              } else {
                const double grj = cur_G[rej + j * m];  // G[rej, j]
                const double gij = cur_G[i + j * m];    // G[i, j]
                new_G[i + j * m] = (gij + gir * grj) * inv_denom;
              }
            }
          } else {
            // Degenerate case (should not occur for valid graphs)
            for (int j = 0; j < m; ++j) {
              new_G[i + j * m] = 0.0;
            }
          }
        }

        // O(1) swap: exchanges internal pointers of the two vectors
        std::swap(cur_G, new_G);

        // Zero out the rejected row and column of the working graph
        for (int i = 0; i < m; ++i) {
          cur_G[rej + i * m] = 0.0;
          cur_G[i + rej * m] = 0.0;
        }
      }

      cur_a[rej] = 0.0;
    }
  }

  return h;
}


// =============================================================================
// Parallel implementation
// =============================================================================

// Worker: processes a chunk [begin, end) of trials.
struct GraphShortcutWorker : public Worker {
  const RMatrix<double>        pvals;    // N x m
  const RMatrix<double>        G_init;   // m x m
  const std::vector<double>&   init_a;   // length m (w * alpha, pre-computed)
  const int                    m;
  RMatrix<int>                 out;      // N x m (LogicalMatrix view)

  GraphShortcutWorker(const NumericMatrix&      p_,
                      const NumericMatrix&      G_,
                      const std::vector<double>& a_,
                      RMatrix<int>              out_)
    : pvals(p_), G_init(G_), init_a(a_), m(p_.ncol()), out(out_) {}

  void operator()(std::size_t begin, std::size_t end) {
    // Thread-local scratch (allocated once per chunk dispatch, reused across
    // the trials in this chunk)
    std::vector<double> cur_a(m);
    std::vector<double> cur_G(m * m);
    std::vector<double> new_G(m * m);
    std::vector<double> cur_p(m);

    for (std::size_t set = begin; set < end; ++set) {
      // Copy initial state
      std::copy(init_a.begin(), init_a.end(), cur_a.begin());
      for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
          cur_G[i + j * m] = G_init(i, j);
        }
      }

      // Extract this trial's p-value row once
      for (int i = 0; i < m; ++i) {
        cur_p[i] = pvals(set, i);
      }

      int sumrej = 0;

      while (true) {
        // Find first hypothesis meeting rejection condition
        int rej = -1;
        for (int i = 0; i < m; ++i) {
          if (cur_p[i] < cur_a[i]) {
            rej = i;
            out(set, i) = 1;     // TRUE in R's logical representation
            ++sumrej;
            break;
          }
        }
        if (rej == -1) break;      // no further rejection possible
        if (sumrej == m) break;    // all hypotheses rejected

        // Update local alphas: a[i] += a[rej] * G[rej, i]
        for (int i = 0; i < m; ++i) {
          cur_a[i] += cur_a[rej] * cur_G[rej + i * m];
        }

        // Update graph only if >1 hypothesis remains
        if (sumrej < m - 1) {
          for (int i = 0; i < m; ++i) {
            const double gir = cur_G[i + rej * m];    // G[i, rej]
            const double gri = cur_G[rej + i * m];    // G[rej, i]
            const double denom = 1.0 - gir * gri;

            if (denom != 0.0) {
              const double inv_denom = 1.0 / denom;
              for (int j = 0; j < m; ++j) {
                if (i == j) {
                  new_G[i + j * m] = 0.0;
                } else {
                  const double grj = cur_G[rej + j * m];  // G[rej, j]
                  const double gij = cur_G[i + j * m];    // G[i, j]
                  new_G[i + j * m] = (gij + gir * grj) * inv_denom;
                }
              }
            } else {
              // Degenerate case (should not occur for valid graphs)
              for (int j = 0; j < m; ++j) {
                new_G[i + j * m] = 0.0;
              }
            }
          }

          std::swap(cur_G, new_G);

          // Zero out the rejected row and column of the working graph
          for (int i = 0; i < m; ++i) {
            cur_G[rej + i * m] = 0.0;
            cur_G[i + rej * m] = 0.0;
          }
        }

        cur_a[rej] = 0.0;
      }
    }
  }
};


//' Graphical shortcut algorithm (Bretz et al. 2009) -- parallel backend
//'
//' Multithreaded implementation of the shortcut algorithm via RcppParallel,
//' parallelised across trials. Matches the behaviour of gMCPLite::graphTest()
//' in the non-parametric, non-entangled, non-upscale case.
//'
//' @param pvals NumericMatrix (N x m) of simulated p-values.
//' @param alpha Scalar significance level.
//' @param w NumericVector (length m) of hypothesis weights.
//' @param G NumericMatrix (m x m) transition matrix. Diagonal must be 0,
//'   row sums must be <= 1.
//' @param num_threads Number of parallel threads (>= 1). Default 1.
//' @param grain_size Chunk size for parallelFor. Default -1 (auto-tuned).
//'   When negative, grain size is chosen so each chunk does approximately
//'   TARGET_OPS worth of inner-loop operations. Per-trial work is O(m^3)
//'   worst case (up to m rejections, each triggering an m*m graph update);
//'   TARGET_OPS = 100000 gives ~195 trials/chunk at m=8 and ~800 trials/chunk
//'   at m=5, leaving thousands of chunks at N=1e6 for good load balancing.
//'
//' @return LogicalMatrix (N x m) of rejection indicators.
//'
//' @seealso \code{\link{graph_shortcut}} for the single-threaded variant.
//' @name graph_shortcut_parallel
//' @noRd
// [[Rcpp::export]]
LogicalMatrix graph_shortcut_parallel(const NumericMatrix& pvals,
                                       const double         alpha,
                                       const NumericVector& w,
                                       const NumericMatrix& G,
                                       int                  num_threads = 1,
                                       int                  grain_size  = -1) {

  const int N = pvals.nrow();
  const int m = pvals.ncol();

  if (G.nrow() != m || G.ncol() != m)
    stop("G must be an m x m matrix matching ncol(pvals).");
  if (w.size() != m)
    stop("w must have length equal to ncol(pvals).");
  if (num_threads < 1)
    stop("num_threads must be >= 1.");

  // Pre-compute initial local alphas once (shared, read-only across threads)
  std::vector<double> init_a(m);
  for (int i = 0; i < m; ++i) {
    init_a[i] = w[i] * alpha;
  }

  // Auto-tune grain size ------------------------------------------------------
  // Per-trial work is O(m^3) worst case: up to m rejections, each triggering
  // an m*m graph update. Typical trials do fewer rejections so actual cost
  // is lower -- m^3 is a conservative upper bound.
  //
  // TARGET_OPS is set higher than apply_ctp_parallel's 15000 because
  // graph_shortcut is algorithmically cheaper per trial (O(m^3) worst case
  // vs O(m * 2^m) for apply_ctp), so chunks need more trials to amortize
  // threading overhead.
  // -------------------------------------------------------------------------
  std::size_t effective_grain_size;
  if (grain_size < 0) {
    double ops_per_row = static_cast<double>(m) * m * m;
    const double TARGET_OPS = 100000.0;

    effective_grain_size = std::max(
      std::size_t(1),
      static_cast<std::size_t>(TARGET_OPS / ops_per_row)
    );

    // Ensure enough chunks for load balancing (>= 4 * num_threads)
    std::size_t min_chunks = static_cast<std::size_t>(num_threads) * 4;
    std::size_t max_grain_for_balance = std::max(
      std::size_t(1),
      N / min_chunks
    );
    effective_grain_size = std::min(effective_grain_size, max_grain_for_balance);

    // Hard bounds
    effective_grain_size = std::max(std::size_t(1), effective_grain_size);
    effective_grain_size = std::min(std::size_t(10000), effective_grain_size);

  } else {
    if (grain_size < 1)
      stop("grain_size must be >= 1 when specified.");
    effective_grain_size = static_cast<std::size_t>(grain_size);
  }

  LogicalMatrix ans(N, m);
  RMatrix<int> ansView(ans);
  GraphShortcutWorker worker(pvals, G, init_a, ansView);

  parallelFor(0, N, worker, effective_grain_size,
              static_cast<std::size_t>(num_threads));

  return ans;
}
