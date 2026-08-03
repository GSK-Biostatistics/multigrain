#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace Rcpp;


// -----------------------------------------------------------------------------
// Applies the closed testing procedure (CTP) for multiple testing correction
// to a matrix of simulated p-values.
//
// @param pvals NumericMatrix of size (N × m), each row representing a
//              simulated trial's m p-values.
// @param alpha Numeric scalar representing the global significance level
//              (e.g., 0.025).
// @param weighting_strategy NumericMatrix (r × 2m)
//        containing intersection hypothesis structure.
//        - First m columns: binary (0/1) membership of hypotheses
//          in intersection hypotheses.
//        - Last m columns: corresponding adjusted local weights.
//
// @return LogicalMatrix (N × m) indicating rejection (TRUE) or non-rejection
//         (FALSE) of each hypothesis.
// [[Rcpp::export]]
LogicalMatrix apply_ctp(const NumericMatrix &pvals,
                        const double         alpha,
                        const NumericMatrix &weighting_strategy) {
  const int N = pvals.nrow();          // # trials
  const int m = pvals.ncol();          // # elementary hypotheses
  const int r = weighting_strategy.nrow();

  if (weighting_strategy.ncol() != 2 * m)
    stop("weighting_strategy must have 2 * ncol(pvals) columns.");

  const double threshold = std::pow(2.0, m - 1);   // 2^(m-1)

  LogicalMatrix out(N, m);           // result

  // ---- main loop over trials -----------------------------------------------
  for (int i = 0; i < N; ++i) {
    // col_sum[j] counts how many rejected intersections include H_j
    std::vector<int> col_sum(m, 0);

    // ---- loop over all intersection hypotheses -----------------------------
    for (int inter = 0; inter < r; ++inter) {

      // --- Is intersection hypothesis rejected ------------------------------
      bool rejected = false;
      for (int j = 0; j < m; ++j) {
        double w = weighting_strategy(inter, j + m);      // adjusted weight
        if (pvals(i, j) <= alpha * w) {                   // test
          rejected = true;
          break;                                          // early exit
        }
      }

      // --- If rejected, update counts for elements in that intersection --
      if (rejected) {
        for (int j = 0; j < m; ++j)
          col_sum[j] += (weighting_strategy(inter, j) > 0.0);
      }
    }

    // ---- Final decision for p-value row -------------------------------------
    for (int j = 0; j < m; ++j)
      out(i, j) = (col_sum[j] == threshold);
  }

  return out;
}



// -----------------------------------------------------------------------------
// Closed testing using Z-statistics
//
// @param z_base NumericMatrix (N × m); row i holds m test statistics Z_{ij}.
// @param NCP m-length vector
// @param alpha double; FWER
// @param weighting_strategy NumericMatrix (r × 2m):
//        * cols [0, ..., m-1]     : 0/1 membership indicators
//        * cols [m, ..., 2m-1]    : adjusted local weights w_{ij} for each
//                                   intersection row i, hypothesis j.
//        These should sum (per row across members) to 1 but we do not assume.
//        Values may be zero (never reject that member in that intersection).
//
// @return LogicalMatrix (N × m) rejection indicators for elementary hypotheses.
//         TRUE if all intersections containing H_j are rejected.
// [[Rcpp::export]]
LogicalMatrix apply_ctp_z(const NumericMatrix &z,
                          const double         alpha,
                          const NumericMatrix &weighting_strategy) {

  const int N = z.nrow();
  const int m = z.ncol();
  const int ws_cols = weighting_strategy.ncol();

  if (ws_cols != 2 * m)
    stop("weighting_strategy must have 2 * ncol(z) columns.");

  const int r = weighting_strategy.nrow();
  const double threshold = std::pow(2.0, m - 1);  // intersections per H_j

  NumericMatrix zcrit(r, m);

  for (int i = 0; i < r; ++i) {
    for (int j = 0; j < m; ++j) {
      const double w = weighting_strategy(i, m + j);
      const double aw = alpha * w;

      if (aw <= 0.0) {
        zcrit(i, j) = R_PosInf;           // never reject via this member
      } else if (aw >= 1.0) {
        zcrit(i, j) = R_NegInf;           // always reject via this member
      } else {
        // upper-tail quantile: P(Z > q) = aw  <=> q = qnorm(aw, lower_tail = FALSE)
        zcrit(i, j) = R::qnorm(aw, 0.0, 1.0, /*lower_tail=*/false, /*log_p=*/false);
      }
    }
  }

  // Main loop
  LogicalMatrix out(N, m);

  for (int i = 0; i < N; ++i) {
    // counts of rejected intersections per hypothesis
    std::vector<int> col_sum(m, 0);

    for (int inter = 0; inter < r; ++inter) {

      // Intersection rejected if ANY member’s Z >= its crit
      bool rejected = false;
      for (int j = 0; j < m; ++j) {
        // If weight is zero, crit = +Inf so condition is FALSE;
        if (z(i, j) >= zcrit(inter, j)) {
          rejected = true;
          break;
        }
      }

      if (rejected) {
        // bump counts for all members present in this intersection
        for (int j = 0; j < m; ++j) {
          if (weighting_strategy(inter, j) > 0.0)
            col_sum[j] += 1;
        }
      }
    }

    // Final decision: all intersections involving H_j were rejected?
    for (int j = 0; j < m; ++j)
      out(i, j) = (col_sum[j] == threshold);
  }

  return out;
}


