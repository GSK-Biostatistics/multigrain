// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;
using namespace RcppParallel;

// -----------------------------------------------------------------------------
// Internal parallel worker for closed testing procedure (CTP).
// This struct is invoked by RcppParallel::parallelFor.
//
// @param pvals NumericMatrix (N × m), p-values for simulated trials.
// @param W NumericMatrix (r × 2m), combined intersection membership and weights.
// @param alpha Significance level (double scalar).
// @param out IntegerMatrix (N × m), stores 0/1 rejections as integers.
// -----------------------------------------------------------------------------
  struct CTPWorker : public Worker {
    const RMatrix<double> pvals;      // N × m
    const RMatrix<double> W;          // r × 2m  (weighting_strategy)
    const double          alpha;
    const int             m, r;
    const double          threshold;  // 2^(m-1)
    RMatrix<int>          out;        // N × m   (result, 0/1)

    CTPWorker(const NumericMatrix&  p_,
              const NumericMatrix&  W_,
              double                alpha_,
              RMatrix<int>          out_)
    : pvals(p_), W(W_), alpha(alpha_),
    m(p_.ncol()), r(W_.nrow()),
    threshold(std::pow(2.0, m - 1)),
    out(out_) {}

    void operator()(std::size_t begin, std::size_t end) {
      for (std::size_t i = begin; i < end; ++i) {          // loop over rows
        std::vector<int> col_sum(m, 0);                    // counts per H_j

        for (int inter = 0; inter < r; ++inter) {          // loop over intersections
          bool rejected = false;
          for (int j = 0; j < m; ++j) {
            //  adjusted weight is in column (m + j)
            if (pvals(i, j) <= alpha * W(inter, m + j)) {
              rejected = true;
              break;                                       // early exit
            }
          }
          if (rejected) {
            // membership indicator is in column j
            for (int j = 0; j < m; ++j)
              col_sum[j] += (W(inter, j) > 0.0);
          }
        }
        // final decision for this trial
        for (int j = 0; j < m; ++j)
          out(i, j) = (col_sum[j] == threshold);
      }
    }
  };


//' Parallel implementation of the closed testing procedure (CTP)
//'
//' Applies the closed testing procedure with multithreading support for
//' evaluating rejection decisions across simulated trials.
//'
//' @param pvals NumericMatrix (N × m), simulated trial p-values.
//' @param alpha Significance level (double scalar).
//' @param weighting_strategy NumericMatrix (r × 2m), intersection memberships
//'   (cols 1:m) and adjusted local weights (cols (m+1):2m).
//' @param num_threads Number of parallel threads. Default is 1 (serial execution).
//'   On HPC systems, set explicitly based on resource allocation.
//' @param grain_size Chunk size for parallel processing. Default is -1
//'   (automatic tuning). When negative, grain size is computed as
//'   \code{15000 / (m * r)} to balance thread overhead with load distribution.
//'
//' @return LogicalMatrix (N × m) of rejection indicators.
//'
//' @details
//' This function uses RcppParallel for multithreaded evaluation of the closed
//' testing procedure. The grain size controls how work is distributed across
//' threads:
//' \itemize{
//'   \item Automatic mode (-1): Adapts to problem size (m, r, N)
//'   \item Manual mode (>0): User-specified chunk size
//' }
//'
//' Thread safety: Always set \code{num_threads} explicitly on shared systems.
//' The default of 1 ensures safe execution on all platforms.
//'
//' @keywords internal
//' @importFrom RcppParallel RcppParallelLibs
// [[Rcpp::export]]
LogicalMatrix apply_ctp_parallel(const NumericMatrix& pvals,
                                 double               alpha,
                                 const NumericMatrix& weighting_strategy,
                                 int                  num_threads = 1,   // Safe default: single thread
                                 int                  grain_size  = -1) { // -1 = auto

  const int N = pvals.nrow();
  const int m = pvals.ncol();
  const int r = weighting_strategy.nrow();

  if (weighting_strategy.ncol() != 2 * m)
    stop("weighting_strategy must have 2 * ncol(pvals) columns.");

  // Validate num_threads
  if (num_threads < 1) {
    stop("num_threads must be >= 1. Use RcppParallel::setThreadOptions() to configure.");
  }

  // Respect RcppParallel's global thread settings
  // Users should set threads explicitly via:
  // RcppParallel::setThreadOptions(numThreads = X)
  // or pass num_threads explicitly to this function

  // Auto-tune grain size based on problem complexity
  std::size_t effective_grain_size;
  if (grain_size < 0) {
    // Estimate operations per row: roughly m * r inner loop iterations
    double ops_per_row = static_cast<double>(m) * r;

    // Target ~10,000-20,000 basic operations per chunk
    const double TARGET_OPS = 15000.0;
    effective_grain_size = std::max(
      std::size_t(1),
      static_cast<std::size_t>(TARGET_OPS / ops_per_row)
    );

    // Ensure enough chunks for load balancing
    // Want at least 4x more chunks than threads
    std::size_t min_chunks = static_cast<std::size_t>(num_threads) * 4;
    std::size_t max_grain_for_balance = std::max(
      std::size_t(1),
      N / min_chunks
    );
    effective_grain_size = std::min(effective_grain_size, max_grain_for_balance);

    // Set reasonable bounds
    effective_grain_size = std::max(std::size_t(1), effective_grain_size);
    effective_grain_size = std::min(std::size_t(10000), effective_grain_size);

  } else {
    if (grain_size < 1) {
      stop("grain_size must be >= 1 when specified.");
    }
    effective_grain_size = static_cast<std::size_t>(grain_size);
  }

  LogicalMatrix ans(N, m);
  RMatrix<int> ansView(ans);
  CTPWorker worker(pvals, weighting_strategy, alpha, ansView);

  parallelFor(0, N, worker, effective_grain_size,
              static_cast<std::size_t>(num_threads));

  return ans;
}
