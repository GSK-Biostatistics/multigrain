# calc_power_pvals stops suboptimal graphs as default

    Code
      calc_power_pvals(alpha = alpha, pvals = pvals, hyp_weight = fs@weights,
      trans_matrix = fs@m, custom_power = list(custom_obj = custom_objective),
      sum_to_one_constraint = TRUE)
    Condition
      Warning:
      One or more rows of `trans_matrix` do not sum to 1 within tolerance.
      Error:
      ! The supplied `hyp_weight` and `trans_matrix` do not build a valid graph.

