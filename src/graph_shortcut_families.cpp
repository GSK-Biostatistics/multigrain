// graph_shortcut_families.cpp
// -----------------------------------------------------------------------------
// Serial and parallel (RcppParallel) shortcut for graphs whose nodes are
// families of hypotheses (Maurer and Bretz 2014, Algorithm 5.7). Within a
// family, hypotheses are tested with an alpha-consistent procedure at the
// family's local level; the family passes its weight on only when every member
// is rejected (the family is "emptied").
//
// Inputs are within-family adjusted p-values (see adjust_pvals_families.cpp)
// and the 0-based family of each hypothesis. Each trial runs in three steps:
//   1. Family p-values: q[f] is the largest adjusted p-value in family f.
//      Family f can be emptied at level a exactly when q[f] < a.
//   2. The Bretz et al. (2009) shortcut on (q, w, G) over the M families, with
//      the same arithmetic as graph_shortcut.cpp (strict `<`, lowest index
//      first). It records which families are emptied and the final local
//      levels.
//   3. Read-out: H_j is rejected if its family was emptied, or if its adjusted
//      p-value is below its family's final local level.
// With every family a singleton the result is identical to graph_shortcut().
//
// The per-trial work is in shortcut_families_trial(), shared by both kernels.
//
// Parallel implementation (graph_shortcut_families_parallel):
//   Parallelised across trials via RcppParallel::parallelFor. Thread-safety is
//   as for graph_shortcut_parallel:
//   * pvals, family, G and init_a are read-only in the Worker.
//   * Each trial writes to a unique row of the output.
//   * Scratch buffers are allocated inside operator(), once per chunk.
// -----------------------------------------------------------------------------

// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <algorithm>
#include <limits>

using namespace Rcpp;
using namespace RcppParallel;


// =============================================================================
// Shared per-trial procedure
// =============================================================================

// Scratch buffers for one trial, reused across the trials of a loop or chunk.
struct FamiliesScratch {
  std::vector<double> p;        // length m: this trial's adjusted p-values
  std::vector<double> q;        // length M: family p-values
  std::vector<double> a;        // length M: local levels
  std::vector<double> cur_G;    // M x M working graph
  std::vector<double> new_G;    // M x M update buffer
  std::vector<int>    emptied;  // length M: family emptied?
  std::vector<int>    rej;      // length m: H_j rejected?

  FamiliesScratch(const int m, const int M)
    : p(m), q(M), a(M), cur_G(M * M), new_G(M * M), emptied(M), rej(m) {}
};

// Runs one trial. On entry s.p holds the trial's adjusted p-values; on exit
// s.rej holds the rejection indicators of the m hypotheses.
static inline void shortcut_families_trial(FamiliesScratch&           s,
                                           const int*                 fam,
                                           const int                  m,
                                           const int                  M,
                                           const std::vector<double>& init_a,
                                           const double*              G_ptr) {
  // 1. Family p-values: the largest adjusted p-value in each family
  std::fill(s.q.begin(), s.q.end(), -std::numeric_limits<double>::infinity());
  for (int j = 0; j < m; ++j) {
    s.q[fam[j]] = std::max(s.q[fam[j]], s.p[j]);
  }

  // 2. Shortcut on the families (arithmetic as in graph_shortcut.cpp)
  std::copy(init_a.begin(), init_a.end(), s.a.begin());
  std::copy(G_ptr, G_ptr + M * M, s.cur_G.begin());
  std::fill(s.emptied.begin(), s.emptied.end(), 0);

  int sumrej = 0;

  while (true) {
    // Find first family that can be emptied
    int rej = -1;
    for (int i = 0; i < M; ++i) {
      if (s.q[i] < s.a[i]) {
        rej = i;
        s.emptied[i] = 1;
        ++sumrej;
        break;
      }
    }
    if (rej == -1) break;      // no further family can be emptied
    if (sumrej == M) break;    // all families emptied

    // Update local levels: a[i] += a[rej] * G[rej, i]
    for (int i = 0; i < M; ++i) {
      s.a[i] += s.a[rej] * s.cur_G[rej + i * M];
    }

    // Update graph only if >1 family remains
    if (sumrej < M - 1) {
      for (int i = 0; i < M; ++i) {
        const double gir = s.cur_G[i + rej * M];  // G[i, rej]
        const double gri = s.cur_G[rej + i * M];  // G[rej, i]
        const double denom = 1.0 - gir * gri;

        if (denom != 0.0) {
          const double inv_denom = 1.0 / denom;
          for (int j = 0; j < M; ++j) {
            if (i == j) {
              s.new_G[i + j * M] = 0.0;
            } else {
              const double grj = s.cur_G[rej + j * M];  // G[rej, j]
              const double gij = s.cur_G[i + j * M];    // G[i, j]
              s.new_G[i + j * M] = (gij + gir * grj) * inv_denom;
            }
          }
        } else {
          // Degenerate case (should not occur for valid graphs)
          for (int j = 0; j < M; ++j) {
            s.new_G[i + j * M] = 0.0;
          }
        }
      }

      std::swap(s.cur_G, s.new_G);

      // Zero out the emptied family's row and column of the working graph
      for (int i = 0; i < M; ++i) {
        s.cur_G[rej + i * M] = 0.0;
        s.cur_G[i + rej * M] = 0.0;
      }
    }

    s.a[rej] = 0.0;
  }

  // 3. Read-out. The emptied flag is needed because the loop sets the local
  //    level of an emptied family to 0.
  for (int j = 0; j < m; ++j) {
    const int f = fam[j];
    s.rej[j] = s.emptied[f] || s.p[j] < s.a[f];
  }
}

// Input checks shared by both kernels.
static void check_families_inputs(const int            m,
                                  const IntegerVector& family,
                                  const NumericVector& w,
                                  const NumericMatrix& G) {
  const int M = w.size();

  if (G.nrow() != M || G.ncol() != M)
    stop("G must be an M x M matrix matching length(w).");
  if (family.size() != m)
    stop("family must have length equal to ncol(pvals).");

  std::vector<int> family_size(M, 0);
  for (int j = 0; j < m; ++j) {
    if (family[j] < 0 || family[j] >= M)
      stop("family must take values in 0, ..., length(w) - 1.");
    ++family_size[family[j]];
  }
  for (int f = 0; f < M; ++f) {
    if (family_size[f] == 0)
      stop("Every family must have at least one hypothesis.");
  }
}


// =============================================================================
// Serial implementation
// =============================================================================

//' Graphical shortcut for families of hypotheses -- serial backend
//'
//' Applies Algorithm 5.7 of Maurer and Bretz (2014) with a Bonferroni test
//' across families to a matrix of within-family adjusted p-values and returns
//' the rejection matrix of the elementary hypotheses.
//'
//' @param pvals NumericMatrix (N x m) of within-family adjusted p-values.
//' @param family IntegerVector (length m): the 0-based family of each column.
//' @param alpha Scalar significance level.
//' @param w NumericVector (length M) of family weights.
//' @param G NumericMatrix (M x M) family transition matrix. Diagonal must be
//'   0, row sums must be <= 1.
//' @return LogicalMatrix (N x m) of rejection indicators.
//'
//' @seealso \code{\link{graph_shortcut_families_parallel}} for the
//'   multithreaded variant.
//' @name graph_shortcut_families
//' @noRd
// [[Rcpp::export]]
LogicalMatrix graph_shortcut_families(const NumericMatrix& pvals,
                                      const IntegerVector& family,
                                      const double         alpha,
                                      const NumericVector& w,
                                      const NumericMatrix& G) {
  const int N = pvals.nrow();
  const int m = pvals.ncol();
  const int M = w.size();

  check_families_inputs(m, family, w, G);

  LogicalMatrix h(N, m);  // zero-initialised by R

  const double* p_ptr = REAL(pvals);
  const double* G_ptr = REAL(G);
  const int* fam = INTEGER(family);
  int* h_ptr = LOGICAL(h);

  std::vector<double> init_a(M);
  for (int i = 0; i < M; ++i) {
    init_a[i] = w[i] * alpha;
  }

  FamiliesScratch s(m, M);

  for (int set = 0; set < N; ++set) {
    for (int j = 0; j < m; ++j) {
      s.p[j] = p_ptr[set + j * N];
    }

    shortcut_families_trial(s, fam, m, M, init_a, G_ptr);

    for (int j = 0; j < m; ++j) {
      h_ptr[set + j * N] = s.rej[j];
    }
  }

  return h;
}


// =============================================================================
// Parallel implementation
// =============================================================================

// Worker: processes a chunk [begin, end) of trials.
struct GraphShortcutFamiliesWorker : public Worker {
  const RMatrix<double>        pvals;   // N x m
  const std::vector<int>&      family;  // length m
  const std::vector<double>&   G;       // M x M, column-major
  const std::vector<double>&   init_a;  // length M (w * alpha)
  const int                    m;
  const int                    M;
  RMatrix<int>                 out;     // N x m (LogicalMatrix view)

  GraphShortcutFamiliesWorker(const NumericMatrix&       p_,
                              const std::vector<int>&    family_,
                              const std::vector<double>& G_,
                              const std::vector<double>& a_,
                              RMatrix<int>               out_)
    : pvals(p_), family(family_), G(G_), init_a(a_),
      m(p_.ncol()), M(a_.size()), out(out_) {}

  void operator()(std::size_t begin, std::size_t end) {
    FamiliesScratch s(m, M);

    for (std::size_t set = begin; set < end; ++set) {
      for (int j = 0; j < m; ++j) {
        s.p[j] = pvals(set, j);
      }

      shortcut_families_trial(s, family.data(), m, M, init_a, G.data());

      for (int j = 0; j < m; ++j) {
        out(set, j) = s.rej[j];
      }
    }
  }
};


//' Graphical shortcut for families of hypotheses -- parallel backend
//'
//' Multithreaded version of `graph_shortcut_families()` via RcppParallel,
//' parallelised across trials. Output is identical to the serial backend.
//'
//' @inheritParams graph_shortcut_families
//' @param num_threads Number of parallel threads (>= 1). Default 1.
//' @param grain_size Chunk size for parallelFor. Default -1 (auto-tuned as in
//'   `graph_shortcut_parallel()`, with per-trial work of order M^3 + m).
//'
//' @return LogicalMatrix (N x m) of rejection indicators.
//'
//' @seealso \code{\link{graph_shortcut_families}} for the single-threaded
//'   variant.
//' @name graph_shortcut_families_parallel
//' @noRd
// [[Rcpp::export]]
LogicalMatrix graph_shortcut_families_parallel(const NumericMatrix& pvals,
                                               const IntegerVector& family,
                                               const double         alpha,
                                               const NumericVector& w,
                                               const NumericMatrix& G,
                                               int num_threads = 1,
                                               int grain_size  = -1) {
  const int N = pvals.nrow();
  const int m = pvals.ncol();
  const int M = w.size();

  check_families_inputs(m, family, w, G);
  if (num_threads < 1)
    stop("num_threads must be >= 1.");

  // Read-only copies shared across threads
  std::vector<int> family_vec(family.begin(), family.end());
  std::vector<double> G_vec(G.begin(), G.end());
  std::vector<double> init_a(M);
  for (int i = 0; i < M; ++i) {
    init_a[i] = w[i] * alpha;
  }

  // Auto-tune grain size as in graph_shortcut_parallel(): up to M rejections,
  // each an M x M graph update, plus the O(m) family maximum and read-out.
  std::size_t effective_grain_size;
  if (grain_size < 0) {
    double ops_per_row = static_cast<double>(M) * M * M + m;
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
  GraphShortcutFamiliesWorker worker(pvals, family_vec, G_vec, init_a, ansView);

  parallelFor(0, N, worker, effective_grain_size,
              static_cast<std::size_t>(num_threads));

  return ans;
}
