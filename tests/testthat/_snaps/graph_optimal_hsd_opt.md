# print() and summary() of a d = 1 result

    Code
      print(obj)
    Output
      <multigrain_graph_optimal>
      Optimal graph found (given user-defined constraints on graph and computational resources):
      
      Hypothesis weights:
       H1  H2  H3 
      0.6 0.4 0.0 
      
      Transition matrix:
          H1  H2  H3
      H1 0.0 0.5 0.5
      H2 1.0 0.0 0.0
      H3 0.5 0.5 0.0
      
      Trial success function:
      d(t1) + d(t2) + d(t3)
      
      Spending functions:
        H1: sfLDOF
        H2: sfHSD(param = -2.25)   [optimised: sec, bounds -4 to 1, step 0.25]
        H3: sfHSD(param = -2.25)   [optimised: sec]
      Gain profile in sec (graph re-optimised at each value):
        within 5e-4 of the maximum on [-3.25, -1.25]; fixed-graph plateau [-2.5, -2]
      
      Value of trial success measure:
      0.7123

---

    Code
      summary(obj)
    Output
      Optimisation summary
      
      Optimisation results:
      
      Hypothesis weights:
       H1  H2  H3 
      0.6 0.4 0.0 
      
      Transition matrix:
          H1  H2  H3
      H1 0.0 0.5 0.5
      H2 1.0 0.0 0.0
      H3 0.5 0.5 0.0
      
      Trial success function (group sequential):
      d(t1) + d(t2) + d(t3)
      Analyses (K): 2
      Discount table d(t): 1.00, 0.75
      
      Spending-parameter profile (sec):
       value   gain gain_fixed_graph   diff      se n_eval
       -4.00 0.7110           0.7001 0.0014 0.00030     41
       -3.75 0.7113           0.7033 0.0010 0.00024     42
       -3.50 0.7116           0.7061 0.0007 0.00020     43
       -3.25 0.7119           0.7083 0.0004 0.00016     44
       -3.00 0.7121           0.7101 0.0003 0.00014     45
       -2.75 0.7122           0.7113 0.0001 0.00012     46
       -2.50 0.7123           0.7121 0.0000 0.00010     47
       -2.25 0.7123           0.7123 0.0000 0.00000     48
       -2.00 0.7123           0.7121 0.0000 0.00010     49
       -1.75 0.7122           0.7113 0.0001 0.00012     50
       -1.50 0.7121           0.7101 0.0003 0.00014     51
       -1.25 0.7119           0.7083 0.0004 0.00016     52
       -1.00 0.7116           0.7061 0.0007 0.00020     53
       -0.75 0.7113           0.7033 0.0010 0.00024     54
       -0.50 0.7110           0.7001 0.0014 0.00030     55
       -0.25 0.7105           0.6963 0.0018 0.00036     56
        0.00 0.7101           0.6921 0.0023 0.00043     57
        0.25 0.7095           0.6873 0.0028 0.00050     58
        0.50 0.7089           0.6821 0.0034 0.00059     59
        0.75 0.7083           0.6763 0.0040 0.00068     60
        1.00 0.7076           0.6701 0.0048 0.00078     61
      GA visited 18 of 21 nodes
      
      Value of trial success measure:
      0.7123
      
      Power metrics:
      
      Graph constraints:
      Constraints on hypothesis weights:
      H1 H2 H3 
      NA NA NA 
      
      Constraints on transition matrix:
         H1 H2 H3
      H1  0 NA NA
      H2 NA  0 NA
      H3 NA NA  0

# print() and summary() of a d = 2 result

    Code
      print(obj)
    Output
      <multigrain_graph_optimal>
      Optimal graph found (given user-defined constraints on graph and computational resources):
      
      Hypothesis weights:
       H1  H2  H3  H4 
      0.5 0.2 0.2 0.1 
      
      Transition matrix:
             H1     H2     H3     H4
      H1 0.0000 0.3333 0.3333 0.3333
      H2 0.3333 0.0000 0.3333 0.3333
      H3 0.3333 0.3333 0.0000 0.3333
      H4 0.3333 0.3333 0.3333 0.0000
      
      Trial success function:
      d(t1) + d(t2) + d(t3) + d(t4)
      
      Spending functions:
        H1: sfLDOF
        H2: sfHSD(param = -2)   [optimised: gamma1, bounds -3 to -1, step 0.5]
        H3: sfHSD(param = -2.5)   [optimised: b, bounds -3 to -2, step 0.5]
        H4: sfHSD(param = -2.5)   [fixed at -2.5: one look at full information]
      Gain profile in gamma1 (graph re-optimised at each value):
        within 5e-4 of the maximum on [-2.5, -1.5]; fixed-graph plateau [-2, -1.5]
      Gain profile in b (graph re-optimised at each value):
        within 5e-4 of the maximum on [-2.5, -2]; fixed-graph plateau [-2.5, -2]
        profile did not converge in 3 rounds; the plateau is indicative
      
      Value of trial success measure:
      0.6543

---

    Code
      summary(obj)
    Output
      Optimisation summary
      
      Optimisation results:
      
      Hypothesis weights:
       H1  H2  H3  H4 
      0.5 0.2 0.2 0.1 
      
      Transition matrix:
             H1     H2     H3     H4
      H1 0.0000 0.3333 0.3333 0.3333
      H2 0.3333 0.0000 0.3333 0.3333
      H3 0.3333 0.3333 0.0000 0.3333
      H4 0.3333 0.3333 0.3333 0.0000
      
      Trial success function (group sequential):
      d(t1) + d(t2) + d(t3) + d(t4)
      Analyses (K): 2
      Discount table d(t): 1.0, 0.9
      
      Spending-parameter profile (gamma1):
       value   gain gain_fixed_graph  diff      se n_eval
        -3.0 0.6523           0.6483 2e-03 0.00039     41
        -2.5 0.6540           0.6523 3e-04 0.00014     42
        -2.0 0.6543           0.6543 0e+00 0.00000     43
        -1.5 0.6542           0.6539 1e-04 0.00011     44
        -1.0 0.6534           0.6513 9e-04 0.00023     45
      no global search
      
      Spending-parameter profile (b):
       value   gain gain_fixed_graph   diff      se n_eval
        -3.0 0.6536           0.6533  7e-04 0.00020     41
        -2.5 0.6543           0.6543  0e+00 0.00000     42
        -2.0 0.6545           0.6542 -2e-04 0.00013     43
      no global search
      
      Value of trial success measure:
      0.6543
      
      Power metrics:
      
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

