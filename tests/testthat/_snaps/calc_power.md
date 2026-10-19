# calc_power_pvals stops suboptimal graphs as default

    Code
      calc_power_pvals(pvals = pvals, hyp_weight = fs@weights, trans_matrix = fs@m,
      custom_power = list(custom_obj = custom_objective), sum_to_one_constraint = TRUE)
    Condition
      Warning:
      One or more rows of `trans_matrix` do not sum to 1 within tolerance.
      Error:
      ! The supplied `hyp_weight` and `trans_matrix` do not build a valid graph.

# calc_power_pvals complains when users passes anything via dots

    Code
      calc_power_pvals(pvals = pvals, hyp_weight = bh@weights, trans_matrix = bh@m,
      alpha)
    Condition
      Error in `calc_power_pvals()`:
      ! `...` must be empty.
      x Problematic argument:
      * ..1 = alpha
      i Did you forget to name an argument?

# calc_power_pvals checks families dimensions

    Code
      calc_power_pvals(pvals, hyp_weight = c(0.5, 0.25, 0.25), trans_matrix = matrix(
        0.5, 3, 3) - diag(0.5, 3), families = fams)
    Condition
      Error:
      ! The graph has 3 nodes but `families` has 2 families.
      i With `families`, the graph has one node per family.

---

    Code
      calc_power_pvals(pvals[, 1:2], hyp_weight = c(0.5, 0.5), trans_matrix = matrix(
        c(0, 1, 1, 0), nrow = 2), families = fams)
    Condition
      Error:
      ! `pvals` has 2 columns but `families` covers 3 hypotheses.

---

    Code
      calc_power_pvals(pvals, hyp_weight = c(0.5, 0.5), trans_matrix = matrix(c(0, 1,
        1, 0), nrow = 2), families = list(1, 2:3))
    Condition
      Error in `calc_power_pvals()`:
      ! `families` must be a multigrain family specification object or `NULL`, not a list.

