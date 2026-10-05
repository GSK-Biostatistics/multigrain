# cpp_code snapshots

    Code
      cat(ts$cpp_code)
    Output
      
      #include <Rcpp.h>
      using namespace Rcpp;
      
      #define std_min std::min
      
      static const double d_tab[] = {0.0, 1.0, 0.75};
      
      // [[Rcpp::export]]
      double powerFunc(IntegerMatrix t) {
          if (t.ncol() < 2) {
              stop("the decision-time matrix has %d column(s) but the trial success function refers to 2 hypotheses", t.ncol());
          }
          int n = t.nrow();
          double total = 0.0;
      
          for (int i = 0; i < n; i++) {
              total += (0.4 * d_tab[t(i, 0)] + 1.0 * d_tab[t(i, 1)]); // gain of simulated trial i from its decision times
          }
      
          return total / n; // Return the mean
      }

---

    Code
      cat(ts2$cpp_code)
    Output
      
      #include <Rcpp.h>
      using namespace Rcpp;
      
      #define std_min std::min
      
      
      // [[Rcpp::export]]
      double powerFunc(IntegerMatrix t) {
          if (t.ncol() < 2) {
              stop("the decision-time matrix has %d column(s) but the trial success function refers to 2 hypotheses", t.ncol());
          }
          int n = t.nrow();
          double total = 0.0;
      
          for (int i = 0; i < n; i++) {
              total += (double(double(t(i, 0)) == 1.0) * double(double(t(i, 1)) == 1.0) + 0.5 * (double(double(t(i, 0)) == 2.0) * double(double(t(i, 1)) == 2.0))); // gain of simulated trial i from its decision times
          }
      
          return total / n; // Return the mean
      }

# print and summary methods

    Code
      print(ts)
    Output
      <multigrain_trial_success_gsd/multigrain_trial_success>
      0.4 * d(t1) + d(t2)
      Analyses (K): 2
      Discount table d(t): 1.00, 0.75

---

    Code
      summary(ts)
    Output
      
      Trial success function (group sequential):
      0.4 * d(t1) + d(t2)
      Analyses (K): 2
      Discount table d(t): 1.00, 0.75

---

    Code
      print(no_k)
    Output
      <multigrain_trial_success_gsd/multigrain_trial_success>
      r1
      Analyses (K): not set

# verbose handling matches trial_success()

    Code
      trial_success_gsd(r1 + r2, verbose = "info")
    Message
      v Trial success function compiled and sourced successfully.
    Output
      <multigrain_trial_success_gsd/multigrain_trial_success>
      r1 + r2
      Analyses (K): not set

---

    Code
      trial_success_gsd(r1 + r2, verbose = "silent")
    Output
      <multigrain_trial_success_gsd/multigrain_trial_success>
      r1 + r2
      Analyses (K): not set

---

    Code
      trial_success_gsd(r1 + r2, verbose = TRUE)
    Message
      v Trial success function compiled and sourced successfully.
    Output
      <multigrain_trial_success_gsd/multigrain_trial_success>
      r1 + r2
      Analyses (K): not set

---

    Code
      trial_success_gsd(r1 + r2, verbose = FALSE)
    Output
      <multigrain_trial_success_gsd/multigrain_trial_success>
      r1 + r2
      Analyses (K): not set

---

    Code
      trial_success_gsd(r1 + r2, verbose = 2)
    Condition
      Error in `trial_success_gsd()`:
      ! `verbose` must be a single string, not the number 2.

