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
       0  1  0  0 
      
      Transition matrix:
             H1     H2     H3     H4
      H1 0.0000 0.8987 0.1013 0.0000
      H2 0.7679 0.0000 0.0000 0.2321
      H3 0.0000 1.0000 0.0000 0.0000
      H4 1.0000 0.0000 0.0000 0.0000
      
      Trial success function:
      0.25 * (2 * (r1 && r2) + r1 * r3 + r2 * r4)
      
      Value of trial success measure:
      0.6393

---

    Code
      summary(graph_optimal_example)
    Output
      Optimisation summary
      
      Optimisation results:
      
      Hypothesis weights:
      H1 H2 H3 H4 
       0  1  0  0 
      
      Transition matrix:
             H1     H2     H3     H4
      H1 0.0000 0.8987 0.1013 0.0000
      H2 0.7679 0.0000 0.0000 0.2321
      H3 0.0000 1.0000 0.0000 0.0000
      H4 1.0000 0.0000 0.0000 0.0000
      
      Trial success function:
      0.25 * (2 * (r1 && r2) + r1 * r3 + r2 * r4)
      
      Value of trial success measure:
      0.6393
      
      Power metrics:
      Power for each hypothesis:
      [1] 0.7253 0.7998 0.4419 0.6648
      Expected number of rejections: 2.63
      Probability of at least one rejection: 0.7998
      Probability of rejecting all hypotheses: 0.4139
      
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

# graph_optimal print & summary methods, group sequential

    Code
      print(obj)
    Output
      <multigrain_graph_optimal>
      Optimal graph found (given user-defined constraints on graph and computational resources):
      
      Hypothesis weights:
      [1] 0.5 0.5
      
      Transition matrix:
           [,1] [,2]
      [1,]    0    1
      [2,]    1    0
      
      Trial success function:
      d(t1) + d(t2)
      
      Stopping rule:
      Every interim: r1
      Hypotheses referenced: up to r1
      
      Value of trial success measure:
      1.227

---

    Code
      summary(obj)
    Output
      Optimisation summary
      
      Optimisation results:
      
      Hypothesis weights:
      [1] 0.5 0.5
      
      Transition matrix:
           [,1] [,2]
      [1,]    0    1
      [2,]    1    0
      
      Trial success function (group sequential):
      d(t1) + d(t2)
      Analyses (K): 2
      Discount table d(t): 1.00, 0.75
      
      Stopping rule (group sequential):
      Every interim: r1
      Hypotheses referenced: up to r1
      
      Value of trial success measure:
      1.227
      
      Power metrics:
      Power for each hypothesis:
      [1] 0.864 0.682
      Expected number of rejections: 1.55
      Probability of at least one rejection: 0.936
      Probability of rejecting all hypotheses: 0.61
      Power by analysis:
         analysis 1 analysis 2
      H1      0.162      0.864
      H2      0.108      0.682
      Decision-time distribution:
         never analysis 1 analysis 2
      H1 0.136      0.162      0.702
      H2 0.318      0.108      0.574
      Stopping-look distribution:
      analysis 1 analysis 2 
           0.162      0.838 
      Mean stopping look: 1.838

---

    Code
      summary(obj_no_rule)
    Output
      Optimisation summary
      
      Optimisation results:
      
      Hypothesis weights:
      [1] 0.5 0.5
      
      Transition matrix:
           [,1] [,2]
      [1,]    0    1
      [2,]    1    0
      
      Trial success function (group sequential):
      d(t1) + d(t2)
      Analyses (K): 2
      Discount table d(t): 1.00, 0.75
      
      Value of trial success measure:
      1.3065
      
      Power metrics:
      Power for each hypothesis:
      [1] 0.864 0.788
      Expected number of rejections: 1.65
      Probability of at least one rejection: 0.936
      Probability of rejecting all hypotheses: 0.716
      Power by analysis:
         analysis 1 analysis 2
      H1      0.162      0.864
      H2      0.108      0.788
      Decision-time distribution:
         never analysis 1 analysis 2
      H1 0.136      0.162      0.702
      H2 0.212      0.108      0.680
      Stopping-look distribution:
      analysis 1 analysis 2 
           0.036      0.964 
      Mean stopping look: 1.964

# summarise_power_object shows the group sequential fields

    Code
      summarise_power_object(power_object)
    Output
      
      Power metrics:
      Power by analysis:
         analysis 1 analysis 2
      H1     0.1235     0.8000
      H2     0.2000     0.6543
      Decision-time distribution:
          never analysis 1 analysis 2
      H1 0.2000     0.1235     0.6766
      H2 0.3457     0.2000     0.4543
      Stopping-look distribution:
      analysis 1 analysis 2 
            0.25       0.75 
      Mean stopping look: 1.75

