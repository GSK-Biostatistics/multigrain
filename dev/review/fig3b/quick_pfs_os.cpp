#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List graph_2h_time(const NumericVector& w,
                         const NumericMatrix& pvals,
                         const NumericMatrix& os_bounds,
                         const double alpha = 0.05,
                         const double tol = 1e-12) {
  // Checks
  if (w.size() != 2) stop("w must have length 2: c(w_PFS, w_OS).");
  if (pvals.ncol() != 3) stop("pvals must have 3 columns: (PFS, OS_IA, OS_FA).");
  if (os_bounds.nrow() != 2 || os_bounds.ncol() != 2)
    stop("os_bounds must be 2x2: rows=(OS initial overall alpha, OS full overall alpha), cols=(IA, FA).");
  if (!(alpha > 0.0 && alpha < 1.0)) stop("alpha must be in (0,1).");

  double w1 = w[0], w2 = w[1]; // PFS, OS
  if (NumericVector::is_na(w1) || NumericVector::is_na(w2)) stop("w contains NA.");
  if (w1 < 0.0 || w2 < 0.0) stop("w must be nonnegative.");
  double wsum = w1 + w2;
  if (std::fabs(wsum - 1.0) > tol) stop("w1 + w2 must equal 1 (within tol).");

  const int N = pvals.nrow();

  // Outputs: columns are (PFS, OS)
  LogicalMatrix rejected(N, 2);
  IntegerMatrix time(N, 2);     // 0 none, 1 time1, 2 time2
  IntegerMatrix rej_type(N, 2); // 0 none, 1 outright, 2 required recycling (counterfactual)

  // PFS cutoffs (PFS p-value is from time 1, may be re-used at time 2)
  const double pfs_cut_init = alpha * w1;
  const double pfs_cut_full = alpha;

  // OS nominal cutoffs (user-supplied)
  const double os_ia_init = os_bounds(0, 0);
  const double os_fa_init = os_bounds(0, 1);
  const double os_ia_full = os_bounds(1, 0);
  const double os_fa_full = os_bounds(1, 1);

  for (int i = 0; i < N; ++i) {
    double p_pfs  = pvals(i, 0);
    double p_os1  = pvals(i, 1); // OS at IA (time 1)
    double p_os2  = pvals(i, 2); // OS at FA (time 2)

    bool ok_pfs = !NumericVector::is_na(p_pfs);
    bool ok_os1 = !NumericVector::is_na(p_os1);
    bool ok_os2 = !NumericVector::is_na(p_os2);

    // ----- TIME 1 under initial allocation -----
    bool r_pfs_init = ok_pfs && (p_pfs <= pfs_cut_init);
    bool r_os_init  = ok_os1 && (p_os1 <= os_ia_init);

    if (r_pfs_init && r_os_init) {
      rejected(i, 0) = true; time(i, 0) = 1; rej_type(i, 0) = 1; // PFS outright
      rejected(i, 1) = true; time(i, 1) = 1; rej_type(i, 1) = 1; // OS outright
      continue;
    }

    if (r_pfs_init && !r_os_init) {
      // PFS rejected at time 1 => recycle to OS (OS gets full alpha)
      rejected(i, 0) = true; time(i, 0) = 1; rej_type(i, 0) = 1;

      // Re-check OS at time 1 with full-alpha boundary
      bool r_os_full_t1 = ok_os1 && (p_os1 <= os_ia_full);
      if (r_os_full_t1) {
        rejected(i, 1) = true;
        time(i, 1) = 1;
        // counterfactual: would it have rejected without recycling?
        rej_type(i, 1) = (p_os1 <= os_ia_init) ? 1 : 2;
        continue;
      }

      // Otherwise, check OS at time 2 with full-alpha boundary
      bool r_os_full_t2 = ok_os2 && (p_os2 <= os_fa_full);
      if (r_os_full_t2) {
        rejected(i, 1) = true;
        time(i, 1) = 2;
        // counterfactual: would it have rejected at time 2 without recycling?
        rej_type(i, 1) = (p_os2 <= os_fa_init) ? 1 : 2;
      }
      continue;
    }

    if (!r_pfs_init && r_os_init) {
      // OS rejected at time 1 => recycle to PFS (PFS gets full alpha)
      rejected(i, 1) = true; time(i, 1) = 1; rej_type(i, 1) = 1;

      // Re-check PFS at time 1 using full alpha (same PFS p-value)
      bool r_pfs_full_t1 = ok_pfs && (p_pfs <= pfs_cut_full);
      if (r_pfs_full_t1) {
        rejected(i, 0) = true; time(i, 0) = 1; rej_type(i, 0) = 2;
      }
      continue;
    }

    // ----- TIME 2: neither rejected at time 1, so OS may reject at FA under initial alpha -----
    bool r_os_t2_init = ok_os2 && (p_os2 <= os_fa_init);
    if (r_os_t2_init) {
      rejected(i, 1) = true; time(i, 1) = 2; rej_type(i, 1) = 1; // OS outright at time 2

      // Now recycle to PFS and "retest" PFS using same PFS p-value
      bool r_pfs_full_t2 = ok_pfs && (p_pfs <= pfs_cut_full);
      if (r_pfs_full_t2) {
        rejected(i, 0) = true; time(i, 0) = 2; rej_type(i, 0) = 2;
      }
    }
  }

  return List::create(
    Named("rejected") = rejected,
    Named("time") = time,
    Named("rej_type") = rej_type,
    Named("colnames") = CharacterVector::create("PFS", "OS")
  );
}
