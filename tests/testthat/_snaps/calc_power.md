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

