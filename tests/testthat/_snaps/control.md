# multigrain_control

    Code
      multigrain_control()
    Output
      <multigrain_control>

# new_multigrain_control() (low-level constructor)

    Code
      new_multigrain_control()
    Output
      <multigrain_control>

# check_control() gives useful error

    Code
      check_control(1)
    Condition
      Error:
      ! `1` must be a multigrain control object, not the number 1.

# check_control() with correct input

    Code
      check_control(multigrain_ctrl)

# check_control() with allow_null

    Code
      check_control(NULL, allow_null = TRUE)

---

    Code
      check_control(NULL, allow_null = FALSE)
    Condition
      Error:
      ! `NULL` must be a multigrain control object, not `NULL`.

# multigrain_control print method

    Code
      control_prepare(multigrain_control(), pvals = test_pvals)
    Output
      <multigrain_control>
      local simulations: 100
      global simulations: 100
      local optimisation:
      * algorithm: "NLOPT_LN_COBYLA"
      * xtol_rel: 5e-08
      * xtol_abs: 5e-09
      * maxeval: 5000
      * print_level: 0
      global optimisation:
      * pcrossover: 0.2
      * pmutation: 0.8
      * maxiter: 1e+05
      * popSize: 200
      * run: 200
      * monitor: FALSE
      * optimArgs: 
          * method: "Nelder-Mead"
          * poptim: 0.2
          * pressel: 0.6

