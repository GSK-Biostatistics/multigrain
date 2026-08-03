# default_control

    Code
      print.default(default_control())
    Output
      $nsim_local
      integer(0)
      
      $nsim_global
      integer(0)
      
      $local_opt
      $local_opt$algorithm
      [1] "NLOPT_LN_COBYLA"
      
      $local_opt$xtol_rel
      [1] 5e-08
      
      $local_opt$xtol_abs
      [1] 5e-09
      
      $local_opt$maxeval
      [1] 5000
      
      $local_opt$print_level
      [1] 0
      
      
      $global_opt
      $global_opt$pcrossover
      [1] 0.2
      
      $global_opt$pmutation
      [1] 0.8
      
      $global_opt$maxiter
      [1] 1e+05
      
      $global_opt$popSize
      [1] 200
      
      $global_opt$run
      [1] 200
      
      $global_opt$monitor
      [1] FALSE
      
      $global_opt$optimArgs
      $global_opt$optimArgs$method
      [1] "Nelder-Mead"
      
      $global_opt$optimArgs$poptim
      [1] 0.2
      
      $global_opt$optimArgs$pressel
      [1] 0.6
      
      
      
      attr(,"class")
      [1] "multigrain_control"

---

    Code
      print.default(default_control())
    Output
      $nsim_local
      integer(0)
      
      $nsim_global
      integer(0)
      
      $local_opt
      $local_opt$algorithm
      [1] "NLOPT_LN_COBYLA"
      
      $local_opt$xtol_rel
      [1] 5e-08
      
      $local_opt$xtol_abs
      [1] 5e-09
      
      $local_opt$maxeval
      [1] 5000
      
      $local_opt$print_level
      [1] 0
      
      
      $global_opt
      $global_opt$pcrossover
      [1] 0.2
      
      $global_opt$pmutation
      [1] 0.8
      
      $global_opt$maxiter
      [1] 1e+05
      
      $global_opt$popSize
      [1] 200
      
      $global_opt$run
      [1] 200
      
      $global_opt$monitor
      [1] FALSE
      
      $global_opt$optimArgs
      $global_opt$optimArgs$method
      [1] "Nelder-Mead"
      
      $global_opt$optimArgs$poptim
      [1] 0.2
      
      $global_opt$optimArgs$pressel
      [1] 0.6
      
      
      
      attr(,"class")
      [1] "multigrain_control"

# control_prepare injects the expected defaults

    Code
      control_prepare(empty_ctrl, pvals = test_pvals)
    Output
      <multigrain_control>
      local simulations: 60
      global simulations: 60
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

