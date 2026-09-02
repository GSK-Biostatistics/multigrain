# graph_optimal_get_control

    Code
      graph_optimal_get_control(graph_custom_power)
    Output
      <multigrain_control>
      local simulations: 100000
      global simulations: 50000
      local optimisation:
      * algorithm: "NLOPT_LN_COBYLA"
      * xtol_rel: 5e-08
      * xtol_abs: 5e-09
      * maxeval: 5000
      * print_level: 1
      global optimisation:
      * pcrossover: 0.2
      * pmutation: 0.8
      * maxiter: 1e+05
      * popSize: 200
      * run: 200
      * monitor: TRUE
      * optimArgs: 
          * method: "Nelder-Mead"
          * poptim: 0.2
          * pressel: 0.6

---

    Code
      graph_optimal_get_control(graph_optimal_example)
    Output
      <multigrain_control>
      local simulations: 20000
      global simulations: 50000
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
      * run: 7
      * monitor: FALSE
      * optimArgs: 
          * method: "Nelder-Mead"
          * poptim: 0.2
          * pressel: 0.6

# summarise helpers

    Code
      summarise_solution_source(NULL)
    Output
      NULL

---

    Code
      summarise_power_object(NULL)
    Output
      NULL

---

    Code
      summarise_solution_source(list(opt_source = "GA"))
    Output
      
      Solution source:
      Global optimisation (genetic algorithm)

---

    Code
      summarise_solution_source(list(opt_source = "GA_minN"))
    Output
      
      Solution source:
      Global optimisation (genetic algorithm; min N)

---

    Code
      summarise_solution_source(list(opt_source = "foo"))
    Output
      
      Solution source:
      foo

# graph_optimal print & summary methods

    Code
      print(obj)
    Output
      <multigrain_graph_optimal>
      Optimal graph found (given user-defined constraints on graph and computational resources):
      
      Hypothesis weights:
       H1  H2  H3 
      0.5 0.3 0.2 
      
      Transition matrix:
           [,1] [,2] [,3]
      [1,]    1    0    0
      [2,]    0    1    0
      [3,]    0    0    1
      
      Trial success function:
      r1 || r2 || r3
      
      Value of trial success measure:
      0.85
      
      Solution source:
      Local optimisation (nloptr)

---

    Code
      summary(obj)
    Output
      Optimisation summary
      
      Optimisation results:
      
      Hypothesis weights:
       H1  H2  H3 
      0.5 0.3 0.2 
      
      Transition matrix:
           [,1] [,2] [,3]
      [1,]    1    0    0
      [2,]    0    1    0
      [3,]    0    0    1
      
      Trial success function:
      r1 || r2 || r3
      
      Value of trial success measure:
      0.85
      
      Power metrics:
      
      Solution source:
      Local optimisation (nloptr)

---

    Code
      print(graph_custom_power)
    Output
      <multigrain_graph_optimal>
      Optimal graph found (given user-defined constraints on graph and computational resources):
      
      Hypothesis weights:
          H1     H2     H3     H4 
      0.6172 0.3828 0.0000 0.0000 
      
      Transition matrix:
             H1     H2     H3     H4
      H1 0.0000 0.7409 0.2591 0.0000
      H2 0.7544 0.0000 0.0000 0.2456
      H3 0.0000 1.0000 0.0000 0.0000
      H4 1.0000 0.0000 0.0000 0.0000
      
      Trial success function:
      0.25 * (2 * (r1 && r2) + r1 * r3 + r2 * r4)
      
      Value of trial success measure:
      0.8197
      
      Solution source:
      Local optimisation (nloptr)

---

    Code
      summary(graph_custom_power)
    Output
      Optimisation summary
      
      Optimisation results:
      
      Hypothesis weights:
          H1     H2     H3     H4 
      0.6172 0.3828 0.0000 0.0000 
      
      Transition matrix:
             H1     H2     H3     H4
      H1 0.0000 0.7409 0.2591 0.0000
      H2 0.7544 0.0000 0.0000 0.2456
      H3 0.0000 1.0000 0.0000 0.0000
      H4 1.0000 0.0000 0.0000 0.0000
      
      Trial success function:
      0.25 * (2 * (r1 && r2) + r1 * r3 + r2 * r4)
      
      Value of trial success measure:
      0.8197
      
      Power metrics:
      Power for each hypothesis:
      [1] 0.9668 0.9036 0.7975 0.7258
      Expected number of rejections: 3.39
      Probability of at least one rejection: 0.9926
      Probability of rejecting all hypotheses: 0.6228
      
      Graph constraints:
      Constraints on hypothesis weights:
      H1 H2 H3 H4 
      NA NA  0  0 
      
      Constraints on transition matrix:
         H1 H2 H3 H4
      H1  0 NA NA  0
      H2 NA  0  0 NA
      H3  0  1  0  0
      H4  1  0  0  0
      
      Solution source:
      Local optimisation (nloptr)

---

    Code
      print(graph_optimal_example)
    Output
      <multigrain_graph_optimal>
      Optimal graph found (given user-defined constraints on graph and computational resources):
      
      Hypothesis weights:
      H1 H2 H3 H4 
       1  0  0  0 
      
      Transition matrix:
         H1 H2 H3 H4
      H1  0  1  0  0
      H2  0  0  1  0
      H3  0  0  0  1
      H4  1  0  0  0
      
      Trial success function:
      r1 && r2 && r3 && r4
      
      Value of trial success measure:
      0.4194
      
      Solution source:
      Local optimisation (nloptr)

---

    Code
      summary(graph_optimal_example)
    Output
      Optimisation summary
      
      Optimisation results:
      
      Hypothesis weights:
      H1 H2 H3 H4 
       1  0  0  0 
      
      Transition matrix:
         H1 H2 H3 H4
      H1  0  1  0  0
      H2  0  0  1  0
      H3  0  0  0  1
      H4  1  0  0  0
      
      Trial success function:
      r1 && r2 && r3 && r4
      
      Value of trial success measure:
      0.4194
      
      Power metrics:
      Power for each hypothesis:
      [1] 0.9005 0.7312 0.4674 0.4194
      Expected number of rejections: 2.52
      Probability of at least one rejection: 0.9005
      Probability of rejecting all hypotheses: 0.4194
      
      Graph constraints:
      Constraints on hypothesis weights:
      H1 H2 H3 H4 
      NA NA NA NA 
      
      Constraints on transition matrix:
         H1 H2 H3 H4
      H1  0 NA NA NA
      H2 NA  0 NA NA
      H3 NA NA  0 NA
      H4 NA NA NA  0
      
      Solution source:
      Local optimisation (nloptr)

