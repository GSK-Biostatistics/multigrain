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
// Optional efficacy stopping rule: stop_rule is a 2^m x (K - 1) logical
// truth table. Row index (0-based) is the bitmask of the hypotheses rejected so
// far, with hypothesis i (0-based) at bit i; column k (0-based) is interim
// analysis k + 1. After the cascade at interim look k finishes, the trial stops
// if the entry for the cumulative rejection set is TRUE: later looks are not
// run and hypotheses not yet rejected keep decision time 0. Independently of
// any rule, a trial in which all m hypotheses are rejected stops at the look of
// the last rejection. A NULL rule, or one with zero columns, means no rule.
// Stopping can only remove rejections, so error control is unchanged.
//
// Outputs: rejected (N x m logical), time (N x m integer, 0 = never) and
// stop_look (length N integer): the 1-based look at which the trial stopped,
// either because the rule fired or because every hypothesis was rejected,
// and K otherwise (so always 1 when K = 1).
// -----------------------------------------------------------------------------

// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <algorithm>

using namespace Rcpp;
using namespace RcppParallel;


// Validate an optional stopping-rule truth table. Returns true and sets `rule`
// when there is a rule to apply; returns false (no rule) when stop_rule is
// NULL or has zero columns. A table with columns must be 2^m x (K - 1).
static bool resolve_stop_rule(const Nullable<LogicalMatrix>& stop_rule,
                              const int                      m,
                              const int                      K,
                              LogicalMatrix&                 rule) {
  if (stop_rule.isNull())
    return false;

  LogicalMatrix candidate(stop_rule.get());
  if (candidate.ncol() == 0)
    return false;

  // 1 << m overflows an int for m > 30; no R matrix has that many rows anyway
  if (m > 30 || candidate.nrow() != (1 << m))
    stop("stop_rule must have 2^m rows.");
  if (candidate.ncol() != K - 1)
    stop("stop_rule must have K - 1 columns.");

  rule = candidate;
  return true;
}


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
//' @param stop_rule Optional efficacy stopping rule: NULL (default) or a
//'   LogicalMatrix (2^m x (K - 1)) truth table. Row r + 1 holds the rejection
//'   set whose bitmask is r (hypothesis i at bit i - 1); column k holds interim
//'   analysis k. After the cascade at interim k the trial stops if the entry
//'   for the hypotheses rejected at analyses 1..k is TRUE. A matrix with zero
//'   columns means no rule.
//' @return A list with `rejected`, a LogicalMatrix (N x m) of rejection
//'   indicators, `time`, an IntegerMatrix (N x m) of decision times
//'   (0 = never rejected), and `stop_look`, an IntegerVector (length N) of the
//'   analysis at which each trial stopped. A trial stops at the first interim
//'   whose rule entry is TRUE, or, with or without a rule, at the analysis of
//'   its last rejection once all m hypotheses are rejected; otherwise
//'   `stop_look` is K. Later analyses are not run, so hypotheses not rejected
//'   by then keep decision time 0.
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
                        const int            K,
                        Rcpp::Nullable<Rcpp::LogicalMatrix> stop_rule =
                          R_NilValue) {
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

  LogicalMatrix rule;
  const bool has_rule = resolve_stop_rule(stop_rule, m, K, rule);
  const int  rule_nr  = has_rule ? rule.nrow() : 0;
  const int* rule_ptr = has_rule ? LOGICAL(rule) : nullptr;

  LogicalMatrix h(N, m);     // zero-initialised by R
  IntegerMatrix tau(N, m);   // zero-initialised by R; 0 = never rejected
  IntegerVector s(N);        // stopping look (1-based)

  // Raw pointers to underlying column-major storage
  const double* p_ptr = REAL(pvals);
  const double* G_ptr = REAL(G);
  const double* w_ptr = REAL(w);
  int* h_ptr   = LOGICAL(h);
  int* tau_ptr = INTEGER(tau);
  int* s_ptr   = INTEGER(s);

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
    s_ptr[set] = K;            // no stop: every look is run

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

      // All hypotheses rejected: nothing left to test, the trial stops here
      if (sumrej == m) {
        s_ptr[set] = k + 1;
        break;
      }

      // Efficacy stopping rule, evaluated at interim looks only
      if (has_rule && k < K - 1) {
        int mask = 0;
        for (int i = 0; i < m; ++i) {
          mask += rej_flag[i] << i;
        }
        if (rule_ptr[mask + k * rule_nr]) {
          s_ptr[set] = k + 1;
          break;
        }
      }
    }
  }

  return List::create(Named("rejected")  = h,
                      Named("time")      = tau,
                      Named("stop_look") = s);
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
  const RMatrix<int>           rule;     // 2^m x (K - 1), or 1 x 1 dummy
  const bool                   has_rule; // false: rule is never read
  RMatrix<int>                 out;      // N x m (LogicalMatrix view)
  RMatrix<int>                 out_time; // N x m (IntegerMatrix view)
  RVector<int>                 out_stop; // N (IntegerVector view)

  GraphShortcutGsdWorker(const NumericMatrix&       p_,
                         const NumericMatrix&       G_,
                         const std::vector<double>& a_,
                         const int                  m_,
                         const int                  K_,
                         const LogicalMatrix&       rule_,
                         const bool                 has_rule_,
                         RMatrix<int>               out_,
                         RMatrix<int>               out_time_,
                         RVector<int>               out_stop_)
    : pvals(p_), G_init(G_), init_a(a_), m(m_), K(K_), rule(rule_),
      has_rule(has_rule_), out(out_), out_time(out_time_),
      out_stop(out_stop_) {}

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
      out_stop[set] = K;         // no stop: every look is run

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

        // All hypotheses rejected: nothing left to test, the trial stops here
        if (sumrej == m) {
          out_stop[set] = k + 1;
          break;
        }

        // Efficacy stopping rule, evaluated at interim looks only
        if (has_rule && k < K - 1) {
          int mask = 0;
          for (int i = 0; i < m; ++i) {
            mask += rej_flag[i] << i;
          }
          if (rule(mask, k)) {
            out_stop[set] = k + 1;
            break;
          }
        }
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
//' @param stop_rule Optional efficacy stopping rule: NULL (default) or a
//'   LogicalMatrix (2^m x (K - 1)) truth table. Row r + 1 holds the rejection
//'   set whose bitmask is r (hypothesis i at bit i - 1); column k holds interim
//'   analysis k. After the cascade at interim k the trial stops if the entry
//'   for the hypotheses rejected at analyses 1..k is TRUE. A matrix with zero
//'   columns means no rule.
//'
//' @return A list with `rejected`, a LogicalMatrix (N x m) of rejection
//'   indicators, `time`, an IntegerMatrix (N x m) of decision times
//'   (0 = never rejected), and `stop_look`, an IntegerVector (length N) of the
//'   analysis at which each trial stopped. A trial stops at the first interim
//'   whose rule entry is TRUE, or, with or without a rule, at the analysis of
//'   its last rejection once all m hypotheses are rejected; otherwise
//'   `stop_look` is K. Later analyses are not run, so hypotheses not rejected
//'   by then keep decision time 0.
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
                                 int                  grain_size  = -1,
                                 Rcpp::Nullable<Rcpp::LogicalMatrix>
                                   stop_rule = R_NilValue) {

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

  // A 1 x 1 dummy stands in when there is no rule; the worker never reads it
  LogicalMatrix rule(1, 1);
  const bool has_rule = resolve_stop_rule(stop_rule, m, K, rule);

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
  IntegerVector ans_stop(N);
  RMatrix<int> ansView(ans);
  RMatrix<int> ansTimeView(ans_time);
  RVector<int> ansStopView(ans_stop);
  GraphShortcutGsdWorker worker(pvals, G, init_a, m, K, rule, has_rule,
                                ansView, ansTimeView, ansStopView);

  parallelFor(0, N, worker, effective_grain_size,
              static_cast<std::size_t>(num_threads));

  return List::create(Named("rejected")  = ans,
                      Named("time")      = ans_time,
                      Named("stop_look") = ans_stop);
}
