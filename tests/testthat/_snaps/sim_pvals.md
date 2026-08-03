# simulate_pvalues complains when users pass anything via `...`

    Code
      simulate_pvalues(nominal_power, alpha_level)
    Condition
      Error in `simulate_pvalues()`:
      ! `...` must be empty.
      x Problematic argument:
      * ..1 = alpha_level
      i Did you forget to name an argument?

