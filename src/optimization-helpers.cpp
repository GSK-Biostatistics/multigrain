#include <Rcpp.h>
using namespace Rcpp;


//' Graph violation score
//'
//' Computes scalar violation score used to penalise invalid graphs during
//' optimisation
//'
//' @param w numeric vector of weights (length m).
//' @param G numeric matrix (m x m) of edge weights.
//' @return scalar violation score (double).
//'
//' @noRd
// [[Rcpp::export]]
double graph_violation_score_cpp(NumericVector w, NumericMatrix G) {

   const R_xlen_t m = w.size();
   double v = 0.0;

   // bounds for w
   for (R_xlen_t i = 0; i < m; ++i) {
     double wi = w[i];
     if (Rcpp::NumericVector::is_na(wi)) continue; // ignore NA
     if (wi < 0.0)      v += -wi;        // amount below 0
     else if (wi > 1.0) v +=  (wi - 1.0);// amount above 1
   }

   // bounds for G
   std::vector<double> rowsum(m, 0.0);
   for (R_xlen_t i = 0; i < m; ++i) {
     for (R_xlen_t j = 0; j < m; ++j) {
       double gij = G(i, j);
       if (Rcpp::NumericVector::is_na(gij)) continue; // ignore NA
       if (gij < 0.0)      v += -gij;
       else if (gij > 1.0) v +=  (gij - 1.0);
       rowsum[i] += gij;
     }
   }

   // abs(sum(w) - 1)
   double wsum = 0.0;
   for (R_xlen_t i = 0; i < m; ++i) {
     double wi = w[i];
     if (!Rcpp::NumericVector::is_na(wi)) wsum += wi;
   }
   v += std::fabs(wsum - 1.0);

   // sum(abs(rowSums(G) - 1))
   for (R_xlen_t i = 0; i < m; ++i) {
     v += std::fabs(rowsum[i] - 1.0);
   }

   return v;
 }
