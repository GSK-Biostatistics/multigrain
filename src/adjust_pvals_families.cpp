// adjust_pvals_families.cpp
// -----------------------------------------------------------------------------
// Within-family adjusted p-values for graphs whose nodes are families of
// hypotheses (Maurer and Bretz 2014).
//
// Each row of `pvals` is one simulated trial. The members of each family are
// adjusted together, exactly as stats::p.adjust() adjusts that family's
// p-values on their own. Adjusted p-values do not depend on the graph, so they
// are computed once per sample, outside the objective function.
//
// Method codes (the position in `family_spec()`'s list of tests, 0-based):
//   0 = Hommel   (port of the O(k^2) loop in stats::p.adjust())
//   1 = Hochberg (step-up)
//   2 = Holm     (step-down)
//
// Implementation:
//   1. Members of each family are gathered once, in increasing column order.
//   2. Per row and family, the member p-values are sorted with a stable
//      insertion sort (ties keep column order, as order() does in p.adjust()),
//      adjusted in place, and written back to the members' columns.
//   3. The arithmetic follows p.adjust() term by term so that results are
//      identical, not just close.
// -----------------------------------------------------------------------------

#include <Rcpp.h>
#include <vector>
#include <algorithm>

using namespace Rcpp;


// =============================================================================
// Adjustment of one family's sorted p-values
// =============================================================================

// Each function takes the k p-values of one family sorted in increasing order
// and overwrites them with the adjusted p-values (still in sorted order).

// Holm: pmin(1, cummax((k - i + 1) * p))
static void adjust_holm(std::vector<double>& p, const int k) {
  double run_max = 0.0;
  for (int i = 0; i < k; ++i) {
    run_max = std::max(run_max, std::min(1.0, (k - i) * p[i]));
    p[i] = run_max;
  }
}

// Hochberg: pmin(1, cummin((k - i + 1) * p)), accumulated from the largest
// p-value down
static void adjust_hochberg(std::vector<double>& p, const int k) {
  double run_min = 1.0;
  for (int i = k - 1; i >= 0; --i) {
    run_min = std::min(run_min, (k - i) * p[i]);
    p[i] = run_min;
  }
}

// Hommel: port of the "hommel" branch of stats::p.adjust(). `q` and `pa` are
// scratch buffers of length >= k. As in p.adjust(), k = 2 uses Hochberg.
static void adjust_hommel(std::vector<double>& p,
                          const int            k,
                          std::vector<double>& q,
                          std::vector<double>& pa) {
  if (k == 2) {
    adjust_hochberg(p, k);
    return;
  }

  // q <- pa <- rep(min(k * p / i), k)
  double init = R_PosInf;
  for (int i = 0; i < k; ++i) {
    init = std::min(init, k * p[i] / (i + 1));
  }
  std::fill(q.begin(), q.begin() + k, init);
  std::fill(pa.begin(), pa.begin() + k, init);

  for (int s = k - 1; s >= 2; --s) {
    // 0-based: i1 = 0..(k - s), i2 = (k - s + 1)..(k - 1), divisors 2..s
    double q1 = R_PosInf;
    for (int i = k - s + 1; i < k; ++i) {
      q1 = std::min(q1, s * p[i] / (i - (k - s) + 1));
    }
    for (int i = 0; i <= k - s; ++i) {
      q[i] = std::min(s * p[i], q1);
    }
    for (int i = k - s + 1; i < k; ++i) {
      q[i] = q[k - s];
    }
    for (int i = 0; i < k; ++i) {
      pa[i] = std::max(pa[i], q[i]);
    }
  }

  for (int i = 0; i < k; ++i) {
    p[i] = std::max(pa[i], p[i]);
  }
}

static void adjust_sorted(std::vector<double>& p,
                          const int            k,
                          const int            method,
                          std::vector<double>& q,
                          std::vector<double>& pa) {
  if (k <= 1) return;  // a single p-value is its own adjusted p-value

  switch (method) {
  case 0: adjust_hommel(p, k, q, pa); break;
  case 1: adjust_hochberg(p, k);      break;
  case 2: adjust_holm(p, k);          break;
  }
}


// =============================================================================
// Row-wise adjustment
// =============================================================================

//' Within-family adjusted p-values
//'
//' Adjusts the p-values of each family separately, row by row, with the
//' family's within-family test. Identical to applying
//' `stats::p.adjust(method = )` to each family's columns of each row.
//'
//' @param pvals NumericMatrix (N x m) of raw p-values.
//' @param family IntegerVector (length m): the 0-based family of each column.
//' @param test IntegerVector (length M): the method code of each family
//'   (0 = Hommel, 1 = Hochberg, 2 = Holm).
//' @return NumericMatrix (N x m) of within-family adjusted p-values.
//'
//' @name adjust_pvals_families
//' @noRd
// [[Rcpp::export]]
NumericMatrix adjust_pvals_families(const NumericMatrix& pvals,
                                    const IntegerVector& family,
                                    const IntegerVector& test) {
  const int N = pvals.nrow();
  const int m = pvals.ncol();
  const int M = test.size();

  if (family.size() != m)
    stop("family must have length equal to ncol(pvals).");

  // Members of each family, in increasing column order
  std::vector< std::vector<int> > members(M);
  for (int j = 0; j < m; ++j) {
    if (family[j] < 0 || family[j] >= M)
      stop("family must take values in 0, ..., length(test) - 1.");
    members[family[j]].push_back(j);
  }
  for (int f = 0; f < M; ++f) {
    if (test[f] < 0 || test[f] > 2)
      stop("test must take values 0 (Hommel), 1 (Hochberg) or 2 (Holm).");
  }

  // Every column belongs to exactly one family, so every entry is written
  NumericMatrix out = no_init(N, m);
  const double* p_ptr = REAL(pvals);
  double* out_ptr = REAL(out);

  // Scratch buffers sized for the largest family
  std::size_t k_max = 0;
  for (int f = 0; f < M; ++f) {
    k_max = std::max(k_max, members[f].size());
  }
  std::vector<int>    ord(k_max);
  std::vector<double> raw(k_max);
  std::vector<double> sorted(k_max);
  std::vector<double> q(k_max);
  std::vector<double> pa(k_max);

  for (int f = 0; f < M; ++f) {
    const std::vector<int>& cols = members[f];
    const int k = cols.size();

    for (int set = 0; set < N; ++set) {
      for (int i = 0; i < k; ++i) {
        raw[i] = p_ptr[set + cols[i] * N];
      }

      // Stable insertion sort of member positions by p-value. Families are
      // small, and unlike std::stable_sort this allocates nothing per row.
      for (int i = 0; i < k; ++i) {
        int pos = i;
        while (pos > 0 && raw[ord[pos - 1]] > raw[i]) {
          ord[pos] = ord[pos - 1];
          --pos;
        }
        ord[pos] = i;
      }
      for (int i = 0; i < k; ++i) {
        sorted[i] = raw[ord[i]];
      }

      adjust_sorted(sorted, k, test[f], q, pa);

      for (int i = 0; i < k; ++i) {
        out_ptr[set + cols[ord[i]] * N] = sorted[i];
      }
    }
  }

  return out;
}
