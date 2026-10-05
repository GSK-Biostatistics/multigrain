// graph_shortcut_gsd.cpp
// -----------------------------------------------------------------------------
// Serial and parallel (RcppParallel) implementations of the group sequential
// graphical procedure (Maurer and Bretz 2013, Algorithm 1) written as the
// fixed-sample shortcut cascade of graph_shortcut.cpp wrapped in an outer loop
// over analyses ("looks").
//
// The kernel never sees a boundary, a spending function or an information
// fraction: the input matrix already holds *repeated* (or, per hypothesis,
// sequential) p-values, so the group sequential rejection rule is the ordinary
// fixed-sample comparison against the hypothesis's current local level. See
// section 4.3 of dev/gsd_design_record.md.
//
// Input layout: pvals is N x (m * K); column (k - 1) * m + i holds hypothesis i
// at look k. That is the column-major reshape of the N x m x K transformed
// array, so R produces it with `dim <- c(N, m * K)` and no copy.
//
// Per trial: cur_a = w * alpha and cur_G = G are initialised once; both are
// carried across looks, which is exactly the graph state Maurer and Bretz carry
// into the next analysis. A rejected hypothesis stays rejected: its flag is set
// and it is skipped in every later scan.
//
// Outputs: rejected (N x m logical) and time (N x m integer, 0 = never).
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

//' Group sequential graphical shortcut algorithm -- fast serial backend
//'
//' Applies the sequentially rejective graphical test procedure of Maurer and
//' Bretz (2013) to a matrix of transformed (repeated or sequential) p-values
//' and returns the rejection matrix together with the analysis at which each
//' rejection was made.
//'
//' @param pvals NumericMatrix (N x m*K) of transformed p-values. Column
//'   (k - 1) * m + i is hypothesis i at look k.
//' @param alpha Scalar significance level.
//' @param w NumericVector (length m) of hypothesis weights.
//' @param G NumericMatrix (m x m) transition matrix. Diagonal must be 0,
//'   row sums must be <= 1.
//' @param K Number of analyses (>= 1). Must divide ncol(pvals).
//' @return A list with `rejected`, a LogicalMatrix (N x m) of rejection
//'   indicators, and `time`, an IntegerMatrix (N x m) of decision times
//'   (0 = never rejected).
//'
//' @seealso \code{\link{graph_shortcut_gsd_parallel}} for the multithreaded
//'   variant.
//' @name graph_shortcut_gsd
//' @noRd
// [[Rcpp::export]]
List graph_shortcut_gsd(const NumericMatrix& pvals,
                        const double         alpha,
                        const NumericVector& w,
                        const NumericMatrix& G,
                        const int            K) {
  const int N  = pvals.nrow();
  const int nc = pvals.ncol();

  if (K < 1)
    stop("K must be >= 1.");
  if (nc % K != 0)
    stop("ncol(pvals) must be a multiple of K.");

  const int m = nc / K;

  if (G.nrow() != m || G.ncol() != m)
    stop("G must be an m x m matrix matching ncol(pvals) / K.");
  if (w.size() != m)
    stop("w must have length equal to ncol(pvals) / K.");

  LogicalMatrix h(N, m);     // zero-initialised by R
  IntegerMatrix tau(N, m);   // zero-initialised by R; 0 = never rejected

  // Raw pointers to underlying column-major storage
  const double* p_ptr = REAL(pvals);
  const double* G_ptr = REAL(G);
  const double* w_ptr = REAL(w);
  int* h_ptr   = LOGICAL(h);
  int* tau_ptr = INTEGER(tau);

  // Pre-allocated working buffers (reused across trials)
  std::vector<double> init_a(m);
  std::vector<double> cur_a(m);
  std::vector<double> cur_G(m * m);
  std::vector<double> new_G(m * m);
  std::vector<double> cur_p(m);
  std::vector<int>    rej_flag(m);

  // Compute initial local alphas once: a[i] = w[i] * alpha
  for (int i = 0; i < m; ++i) {
    init_a[i] = w_ptr[i] * alpha;
  }

  for (int set = 0; set < N; ++set) {
    // Copy initial state for this trial; both are carried across looks
    std::copy(init_a.begin(), init_a.end(), cur_a.begin());
    std::copy(G_ptr, G_ptr + m * m, cur_G.begin());
    std::fill(rej_flag.begin(), rej_flag.end(), 0);

    int sumrej = 0;

    for (int k = 0; k < K; ++k) {
      // Extract this trial's p-value row for look k (column-major -> dense
      // local copy)
      const int offset = k * m;
      for (int i = 0; i < m; ++i) {
        cur_p[i] = p_ptr[set + (offset + i) * N];
      }

      while (true) {
        // Find first not-yet-rejected hypothesis meeting rejection condition
        int rej = -1;
        for (int i = 0; i < m; ++i) {
          if (!rej_flag[i] && cur_p[i] < cur_a[i]) {
            rej = i;
            rej_flag[i] = 1;
            h_ptr[set + i * N]   = 1;       // TRUE in R's logical representation
            tau_ptr[set + i * N] = k + 1;   // 1-based look index
            ++sumrej;
            break;
          }
        }
        if (rej == -1) break;      // no further rejection possible at this look
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

      if (sumrej == m) break;    // nothing left to test at later looks
    }
  }

  return List::create(Named("rejected") = h,
                      Named("time")     = tau);
}


// =============================================================================
// Parallel implementation
// =============================================================================

// Worker: processes a chunk [begin, end) of trials.
struct GraphShortcutGsdWorker : public Worker {
  const RMatrix<double>        pvals;    // N x m*K
  const RMatrix<double>        G_init;   // m x m
  const std::vector<double>&   init_a;   // length m (w * alpha, pre-computed)
  const int                    m;
  const int                    K;
  RMatrix<int>                 out;      // N x m (LogicalMatrix view)
  RMatrix<int>                 out_time; // N x m (IntegerMatrix view)

  GraphShortcutGsdWorker(const NumericMatrix&       p_,
                         const NumericMatrix&       G_,
                         const std::vector<double>& a_,
                         const int                  m_,
                         const int                  K_,
                         RMatrix<int>               out_,
                         RMatrix<int>               out_time_)
    : pvals(p_), G_init(G_), init_a(a_), m(m_), K(K_), out(out_),
      out_time(out_time_) {}

  void operator()(std::size_t begin, std::size_t end) {
    // Thread-local scratch (allocated once per chunk dispatch, reused across
    // the trials in this chunk)
    std::vector<double> cur_a(m);
    std::vector<double> cur_G(m * m);
    std::vector<double> new_G(m * m);
    std::vector<double> cur_p(m);
    std::vector<int>    rej_flag(m);

    for (std::size_t set = begin; set < end; ++set) {
      // Copy initial state
      std::copy(init_a.begin(), init_a.end(), cur_a.begin());
      for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
          cur_G[i + j * m] = G_init(i, j);
        }
      }
      std::fill(rej_flag.begin(), rej_flag.end(), 0);

      int sumrej = 0;

      for (int k = 0; k < K; ++k) {
        // Extract this trial's p-value row for look k
        const int offset = k * m;
        for (int i = 0; i < m; ++i) {
          cur_p[i] = pvals(set, offset + i);
        }

        while (true) {
          // Find first not-yet-rejected hypothesis meeting rejection condition
          int rej = -1;
          for (int i = 0; i < m; ++i) {
            if (!rej_flag[i] && cur_p[i] < cur_a[i]) {
              rej = i;
              rej_flag[i] = 1;
              out(set, i)      = 1;      // TRUE in R's logical representation
              out_time(set, i) = k + 1;  // 1-based look index
              ++sumrej;
              break;
            }
          }
          if (rej == -1) break;      // no further rejection at this look
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

        if (sumrej == m) break;    // nothing left to test at later looks
      }
    }
  }
};


//' Group sequential graphical shortcut algorithm -- parallel backend
//'
//' Multithreaded implementation of the group sequential shortcut algorithm via
//' RcppParallel, parallelised across trials. Bit-identical to
//' \code{graph_shortcut_gsd()}.
//'
//' @param pvals NumericMatrix (N x m*K) of transformed p-values. Column
//'   (k - 1) * m + i is hypothesis i at look k.
//' @param alpha Scalar significance level.
//' @param w NumericVector (length m) of hypothesis weights.
//' @param G NumericMatrix (m x m) transition matrix. Diagonal must be 0,
//'   row sums must be <= 1.
//' @param K Number of analyses (>= 1). Must divide ncol(pvals).
//' @param num_threads Number of parallel threads (>= 1). Default 1.
//' @param grain_size Chunk size for parallelFor. Default -1 (auto-tuned).
//'   When negative, grain size is chosen so each chunk does approximately
//'   TARGET_OPS worth of inner-loop operations. Per-trial work is O(K * m^3)
//'   worst case (K looks, each with up to m rejections triggering an m*m
//'   graph update).
//'
//' @return A list with `rejected`, a LogicalMatrix (N x m) of rejection
//'   indicators, and `time`, an IntegerMatrix (N x m) of decision times
//'   (0 = never rejected).
//'
//' @seealso \code{\link{graph_shortcut_gsd}} for the single-threaded variant.
//' @name graph_shortcut_gsd_parallel
//' @noRd
// [[Rcpp::export]]
List graph_shortcut_gsd_parallel(const NumericMatrix& pvals,
                                 const double         alpha,
                                 const NumericVector& w,
                                 const NumericMatrix& G,
                                 const int            K,
                                 int                  num_threads = 1,
                                 int                  grain_size  = -1) {

  const int N  = pvals.nrow();
  const int nc = pvals.ncol();

  if (K < 1)
    stop("K must be >= 1.");
  if (nc % K != 0)
    stop("ncol(pvals) must be a multiple of K.");

  const int m = nc / K;

  if (G.nrow() != m || G.ncol() != m)
    stop("G must be an m x m matrix matching ncol(pvals) / K.");
  if (w.size() != m)
    stop("w must have length equal to ncol(pvals) / K.");
  if (num_threads < 1)
    stop("num_threads must be >= 1.");

  // Pre-compute initial local alphas once (shared, read-only across threads)
  std::vector<double> init_a(m);
  for (int i = 0; i < m; ++i) {
    init_a[i] = w[i] * alpha;
  }

  // Auto-tune grain size ------------------------------------------------------
  // Per-trial work is O(K * m^3) worst case: K looks, each with up to m
  // rejections, each triggering an m*m graph update. Typical trials do fewer
  // rejections, so K * m^3 is a conservative upper bound. TARGET_OPS matches
  // the fixed-sample kernel.
  // -------------------------------------------------------------------------
  std::size_t effective_grain_size;
  if (grain_size < 0) {
    double ops_per_row = static_cast<double>(K) * m * m * m;
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
  IntegerMatrix ans_time(N, m);
  RMatrix<int> ansView(ans);
  RMatrix<int> ansTimeView(ans_time);
  GraphShortcutGsdWorker worker(pvals, G, init_a, m, K, ansView, ansTimeView);

  parallelFor(0, N, worker, effective_grain_size,
              static_cast<std::size_t>(num_threads));

  return List::create(Named("rejected") = ans,
                      Named("time")     = ans_time);
}
