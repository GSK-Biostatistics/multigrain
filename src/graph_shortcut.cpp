// graph_shortcut.cpp
// -----------------------------------------------------------------------------
// Serial and parallel (RcppParallel) implementations of the Bretz et al.
// (2009) shortcut algorithm for graph-based multiple testing.
//
// Both functions produce a rejection matrix only (0/1 entries) as a
// LogicalMatrix, matching the convention of apply_ctp / apply_ctp_parallel.
//
// Two kernels sit behind each exported function; the kernel is chosen at
// entry (see "Kernel selection" below):
//
// Table kernel:
//   The local critical values alpha * w_j(S) depend only on the set S of
//   unrejected hypotheses, so they are computed once per call for all
//   2^m - 1 non-empty S (O(2^m m^2) build, row-major table of 2^m * m
//   doubles indexed by the bitmask of S). Each trial then walks the table:
//   reject every j in S with p_j < c_j(S), remove them from S, repeat.
//   Each step is O(m) and there is no per-trial copy of the graph. All
//   hypotheses that qualify at S are rejected together; this reaches the
//   same final set because local weights only grow as S shrinks.
//
// Sequential kernel (the original implementation, kept as the fallback for
// large m and as the test oracle):
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
//   * pvals, G, init_a and the critical-value table are read-only in the
//     Workers. The table is built once on the main thread before
//     parallelFor.
//   * Each trial writes to a unique row of the output, so there are no
//     write-write races.
//   * Scratch buffers (cur_a, cur_G, new_G, cur_p) are allocated inside
//     operator(), i.e. once per chunk dispatch. TBB's thread-local allocator
//     caches make these allocations effectively free in practice.
//
// Kernel selection:
//   Two internal (undocumented) options are read once per call, on the main
//   thread, before any parallel work:
//   * multigrain.shortcut_kernel: "auto" (default), "sequential" or "table".
//   * multigrain.shortcut_table_max_m: whole number in [0, 20] (default 20);
//     caps the table kernel under "auto".
//   Under "auto" the table kernel is used when m <= max_m and
//   2^m <= TABLE_KAPPA * N. The rule deliberately ignores num_threads, so the
//   kernel (and therefore the output) never depends on the thread count.
// -----------------------------------------------------------------------------

// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>

using namespace Rcpp;
using namespace RcppParallel;


// =============================================================================
// Kernel selection
// =============================================================================

namespace {

// Largest m for which the table kernel may run: 2^20 * 20 doubles = 160 MB.
// Also keeps the 32-bit subset mask valid.
constexpr int TABLE_HARD_MAX_M = 20;

// "auto" uses the table kernel when 2^m <= TABLE_KAPPA * N. The table build
// costs O(2^m m^2) once per call; the sequential kernel costs
// O(N m^2 (1 + r)) for r rejections per trial, so the break-even ratio
// 2^m / N grows with r. Measured serially (dev/review/shortcut_table/):
// the table kernel was faster in every case up to 2^m / N = 1.64, and the
// sequential kernel first won (at low power) at 2.62. TABLE_KAPPA = 2 keeps
// "auto" no slower than the sequential kernel.
constexpr double TABLE_KAPPA = 2.0;

// Target inner-loop operations per parallelFor chunk for the table walk,
// with m^2 operations assumed per trial.
constexpr double TABLE_TARGET_OPS = 100000.0;

enum class ShortcutKernel { Auto, Sequential, Table };

struct ShortcutConfig {
  ShortcutKernel kernel;
  int            table_max_m;
};

// Read and validate the kernel options. Main thread only (calls the R API).
ShortcutConfig read_shortcut_config() {
  ShortcutConfig cfg = {ShortcutKernel::Auto, TABLE_HARD_MAX_M};

  SEXP kernel = Rf_GetOption1(Rf_install("multigrain.shortcut_kernel"));
  if (!Rf_isNull(kernel)) {
    std::string value;
    if (TYPEOF(kernel) == STRSXP && XLENGTH(kernel) == 1 &&
        STRING_ELT(kernel, 0) != NA_STRING) {
      value = CHAR(STRING_ELT(kernel, 0));
    }
    if (value == "auto") {
      cfg.kernel = ShortcutKernel::Auto;
    } else if (value == "sequential") {
      cfg.kernel = ShortcutKernel::Sequential;
    } else if (value == "table") {
      cfg.kernel = ShortcutKernel::Table;
    } else {
      stop("Option `multigrain.shortcut_kernel` must be one of \"auto\", "
           "\"sequential\" or \"table\".");
    }
  }

  SEXP max_m = Rf_GetOption1(Rf_install("multigrain.shortcut_table_max_m"));
  if (!Rf_isNull(max_m)) {
    const bool numeric = (TYPEOF(max_m) == INTSXP || TYPEOF(max_m) == REALSXP) &&
                         XLENGTH(max_m) == 1;
    const double value = numeric ? Rf_asReal(max_m) : NA_REAL;
    if (ISNAN(value) || value != std::floor(value) || value < 0 ||
        value > TABLE_HARD_MAX_M) {
      stop("Option `multigrain.shortcut_table_max_m` must be a whole number "
           "between 0 and %d.", TABLE_HARD_MAX_M);
    }
    cfg.table_max_m = static_cast<int>(value);
  }

  return cfg;
}

// TRUE when the table kernel should run for N trials of m hypotheses.
bool use_table_kernel(const int m, const int N, const ShortcutConfig& cfg) {
  if (m < 1) return false;

  switch (cfg.kernel) {
    case ShortcutKernel::Sequential:
      return false;
    case ShortcutKernel::Table:
      if (m > TABLE_HARD_MAX_M)
        stop("The table kernel supports at most %d hypotheses (m = %d).",
             TABLE_HARD_MAX_M, m);
      return true;
    case ShortcutKernel::Auto:
    default:
      if (m > cfg.table_max_m) return false;
      return std::ldexp(1.0, m) <= TABLE_KAPPA * static_cast<double>(N);
  }
}

// Grain size for parallelFor: user value when grain_size >= 1, otherwise
// auto-tuned so each chunk does about target_ops operations, with at least
// 4 * num_threads chunks for load balancing.
std::size_t choose_grain_size(const int    grain_size,
                              const double ops_per_row,
                              const double target_ops,
                              const int    N,
                              const int    num_threads) {
  std::size_t effective_grain_size;
  if (grain_size < 0) {
    effective_grain_size = std::max(
      std::size_t(1),
      static_cast<std::size_t>(target_ops / ops_per_row)
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
  return effective_grain_size;
}


// =============================================================================
// Table kernel
// =============================================================================

// Builds crit[S * m + i] = alpha * w_{i+1}(S) for every S (bit i of S is
// H_{i+1}; row 0, the empty set, is all zero).
//
// Depth-first over subsets, removing hypotheses in increasing index order so
// each subset is reached from exactly one parent. One (a, G) pair is kept per
// depth. The arithmetic is the sequential kernel's: local levels start at
// w[i] * alpha and are updated in alpha-space (a_i + a_k g_ki), with the same
// graph update and degenerate-denominator rule. A set reached by rejecting in
// increasing index order therefore gets bit-identical critical values.
class CritTableBuilder {
public:
  CritTableBuilder(const int m, double* crit)
    : m_(m), crit_(crit),
      as_(m + 1, std::vector<double>(m)),
      Gs_(m + 1, std::vector<double>(static_cast<std::size_t>(m) * m)) {}

  void build(const double* w, const double alpha, const double* G) {
    for (int i = 0; i < m_; ++i) {
      as_[0][i] = w[i] * alpha;
    }
    std::copy(G, G + m_ * m_, Gs_[0].begin());
    dfs((std::uint32_t(1) << m_) - 1, 0, -1);
  }

private:
  const int                        m_;
  double*                          crit_;
  std::vector<std::vector<double>> as_;  // local levels, one per depth
  std::vector<std::vector<double>> Gs_;  // graphs, one per depth

  // S contains every index above `last`, so each k > last is removable.
  void dfs(const std::uint32_t S, const int depth, const int last) {
    const std::vector<double>& a = as_[depth];
    const std::vector<double>& G = Gs_[depth];
    const int m = m_;

    double* row = crit_ + static_cast<std::size_t>(S) * m;
    for (int i = 0; i < m; ++i) {
      row[i] = a[i];
    }

    for (int k = last + 1; k < m; ++k) {
      const std::uint32_t child = S & ~(std::uint32_t(1) << k);
      if (child == 0) continue;

      // Local levels: a[i] += a[k] * G[k, i]; a[k] = 0
      std::vector<double>& ca = as_[depth + 1];
      for (int i = 0; i < m; ++i) {
        ca[i] = (i == k) ? 0.0 : a[i] + a[k] * G[k + i * m];
      }

      // Graph update, needed only when the child has children of its own
      if (k < m - 1) {
        std::vector<double>& cG = Gs_[depth + 1];
        for (int i = 0; i < m; ++i) {
          if (i == k) {
            for (int j = 0; j < m; ++j) cG[i + j * m] = 0.0;
            continue;
          }
          const double gik = G[i + k * m];  // G[i, k]
          const double gki = G[k + i * m];  // G[k, i]
          const double denom = 1.0 - gik * gki;

          if (denom != 0.0) {
            const double inv_denom = 1.0 / denom;
            for (int j = 0; j < m; ++j) {
              if (i == j || j == k) {
                cG[i + j * m] = 0.0;
              } else {
                cG[i + j * m] = (G[i + j * m] + gik * G[k + j * m]) * inv_denom;
              }
            }
          } else {
            // Degenerate case (should not occur for valid graphs)
            for (int j = 0; j < m; ++j) cG[i + j * m] = 0.0;
          }
        }
      }

      dfs(child, depth + 1, k);
    }
  }
};

void build_crit_table(const int m, const double alpha, const double* w,
                      const double* G, std::vector<double>& crit) {
  crit.assign((std::size_t(1) << m) * static_cast<std::size_t>(m), 0.0);
  CritTableBuilder builder(m, crit.data());
  builder.build(w, alpha, G);
}

// Walks the table for trials [begin, end). Strict "<", and only unrejected
// hypotheses are tested, so a hypothesis with p = 0 and level 0 is never
// rejected. h must be zero-initialised; only rejections are written.
void walk_crit_table(const double* p_ptr, const std::size_t N, const int m,
                     const double* crit, int* h_ptr,
                     const std::size_t begin, const std::size_t end) {
  const std::uint32_t full = (std::uint32_t(1) << m) - 1;
  std::vector<double> cur_p(m);

  for (std::size_t set = begin; set < end; ++set) {
    // Extract this trial's p-value row once
    for (int i = 0; i < m; ++i) {
      cur_p[i] = p_ptr[set + i * N];
    }

    std::uint32_t S = full;
    while (S) {
      const double* row = crit + static_cast<std::size_t>(S) * m;
      std::uint32_t rej = 0;
      for (int i = 0; i < m; ++i) {
        rej |= static_cast<std::uint32_t>(((S >> i) & 1u) & (cur_p[i] < row[i]))
               << i;
      }
      if (!rej) break;
      S &= ~rej;
    }

    const std::uint32_t R = full & ~S;
    for (int i = 0; i < m; ++i) {
      if ((R >> i) & 1u) h_ptr[set + i * N] = 1;  // TRUE
    }
  }
}

// Worker: walks the shared, read-only table for a chunk [begin, end) of
// trials. Raw pointers are taken on the main thread before parallelFor.
struct TableShortcutWorker : public Worker {
  const double*     p_ptr;  // N x m, column-major
  const std::size_t N;
  const int         m;
  const double*     crit;   // 2^m x m, row-major
  int*              h_ptr;  // N x m, column-major

  TableShortcutWorker(const double* p_, const std::size_t N_, const int m_,
                      const double* crit_, int* h_)
    : p_ptr(p_), N(N_), m(m_), crit(crit_), h_ptr(h_) {}

  void operator()(std::size_t begin, std::size_t end) {
    walk_crit_table(p_ptr, N, m, crit, h_ptr, begin, end);
  }
};

}  // namespace


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

  const ShortcutConfig cfg = read_shortcut_config();

  LogicalMatrix h(N, m);  // zero-initialised by R

  // Raw pointers to underlying column-major storage
  const double* p_ptr = REAL(pvals);
  const double* G_ptr = REAL(G);
  const double* w_ptr = REAL(w);
  int* h_ptr = LOGICAL(h);

  if (use_table_kernel(m, N, cfg)) {
    std::vector<double> crit;
    build_crit_table(m, alpha, w_ptr, G_ptr, crit);
    walk_crit_table(p_ptr, static_cast<std::size_t>(N), m, crit.data(), h_ptr,
                    0, static_cast<std::size_t>(N));
    return h;
  }

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
//'   a target number of inner-loop operations. For the sequential kernel,
//'   per-trial work is O(m^3) worst case (up to m rejections, each
//'   triggering an m*m graph update) and the target is 100000, giving ~195
//'   trials/chunk at m=8 and ~800 trials/chunk at m=5. For the table kernel,
//'   per-trial work is taken as m^2 and the target is TABLE_TARGET_OPS.
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

  const ShortcutConfig cfg = read_shortcut_config();

  if (use_table_kernel(m, N, cfg)) {
    const std::size_t effective_grain_size = choose_grain_size(
      grain_size, static_cast<double>(m) * m, TABLE_TARGET_OPS, N, num_threads
    );

    // Built once on the main thread; shared read-only by all workers
    std::vector<double> crit;
    build_crit_table(m, alpha, REAL(w), REAL(G), crit);

    LogicalMatrix ans(N, m);  // zero-initialised by R
    TableShortcutWorker worker(REAL(pvals), static_cast<std::size_t>(N), m,
                               crit.data(), LOGICAL(ans));

    parallelFor(0, N, worker, effective_grain_size,
                static_cast<std::size_t>(num_threads));

    return ans;
  }

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
  const std::size_t effective_grain_size = choose_grain_size(
    grain_size, static_cast<double>(m) * m * m, 100000.0, N, num_threads
  );

  LogicalMatrix ans(N, m);
  RMatrix<int> ansView(ans);
  GraphShortcutWorker worker(pvals, G, init_a, ansView);

  parallelFor(0, N, worker, effective_grain_size,
              static_cast<std::size_t>(num_threads));

  return ans;
}


// =============================================================================
// Internal helpers for tests and development scripts
// =============================================================================

//' Kernel chosen by graph_shortcut() under the current options
//'
//' @param n_trials Number of trials (rows of pvals).
//' @param m Number of hypotheses.
//' @return "table" or "sequential".
//' @name graph_shortcut_kernel
//' @noRd
// [[Rcpp::export]]
std::string graph_shortcut_kernel(const int n_trials, const int m) {
  const ShortcutConfig cfg = read_shortcut_config();
  return use_table_kernel(m, n_trials, cfg) ? "table" : "sequential";
}

//' Critical-value table used by the table kernel
//'
//' @param alpha Scalar significance level.
//' @param w NumericVector (length m) of hypothesis weights.
//' @param G NumericMatrix (m x m) transition matrix.
//' @return NumericMatrix (2^m x m). Row `mask + 1` holds alpha * w_j(S) for
//'   the set S whose bit j - 1 is H_j; 0 for hypotheses outside S. Row 1
//'   (the empty set) is all zero.
//' @name graph_shortcut_table
//' @noRd
// [[Rcpp::export]]
NumericMatrix graph_shortcut_table(const double         alpha,
                                   const NumericVector& w,
                                   const NumericMatrix& G) {
  const int m = w.size();

  if (G.nrow() != m || G.ncol() != m)
    stop("G must be an m x m matrix matching length(w).");
  if (m < 1 || m > TABLE_HARD_MAX_M)
    stop("The table kernel supports 1 to %d hypotheses (m = %d).",
         TABLE_HARD_MAX_M, m);

  std::vector<double> crit;
  build_crit_table(m, alpha, REAL(w), REAL(G), crit);

  const int n_sets = 1 << m;
  NumericMatrix out(n_sets, m);
  for (int S = 0; S < n_sets; ++S) {
    for (int i = 0; i < m; ++i) {
      out(S, i) = crit[static_cast<std::size_t>(S) * m + i];
    }
  }
  return out;
}
