# validate_graph_constraint() works

    Code
      validate_graph_constraint(graph_constraint(c(NA, NA, 0, 0, 0), matrix(c(0, 0.8,
        0.2, 0, 0, NA, 0, 0, NA, 0, 0, 0.9, 0, 0.1, 0, 0.9, NA, NA, 0, 0, 1, 0, 0, 0,
        0), nrow = 5, byrow = TRUE)), diagnose = TRUE)
    Message
      
      -- Hypothesis weight constraint diagnosis: 
      v All hypothesis weights constraint values are less than or equal to 1.
      v All hypothesis weights constraint values are greater than or equal to 0.
      v The hypothesis weights constraint vector is correctly defined for optimising 2 weights.
      v The sum of the incomplete hypothesis weights constraint vector is less than 1.
      
      -- Transition matrix constraint diagnosis: 
      v All values on the transition matrix diagonal are 0.
      v The transition matrix is square.
      v All transition matrix values are less than or equal to 1.
      v All transition matrix values are greater than or equal to 0.
      v All transition matrix rows are either complete or have at least 2 values to optimise.
      v All complete transition matrix rows sum up to 1.
      v No incomplete transition matrix rows have a sum greater than 1.
      v No incomplete transition matrix rows have a sum equal to 1.
      
      -- Hypothesis weight and transition matrix constraints consistency: 
      v The hypothesis weights vector has 5 elements.
      v The transition matrix has 5 columns and 5 rows.
    Output
      <multigrain_graph_constraint>
      Constraints on hypothesis weights:
      H1 H2 H3 H4 H5 
      NA NA  0  0  0 
      
      Constraints on transition matrix:
          H1  H2  H3  H4 H5
      H1 0.0 0.8 0.2 0.0  0
      H2  NA 0.0 0.0  NA  0
      H3 0.0 0.9 0.0 0.1  0
      H4 0.9  NA  NA 0.0  0
      H5 1.0 0.0 0.0 0.0  0

---

    Code
      validate_graph_constraint(new_graph_constraint(c(NA, NA, 0, 0), matrix(c(0, 0.8,
        0.1, 0, 0, NA, 0, 0, NA, 0, 0, 0.9, 0, 0.2, 0, 1, 0, NA, 0, 0, 1.1, 0, 0, 0,
        0), nrow = 5, byrow = TRUE)))
    Condition
      Error:
      ! Some transition matrix values are greater than 1.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

# assert_hc_tc_consistency() works

    Code
      assert_hc_tc_consistency(new_graph_constraint(c(NA, NA, 0, 0), matrix(c(0, 0.8,
        0.1, 0, 0, NA, 0, 0, NA, 0, 0, 0.9, 0, 0.2, 0, 1, 0, NA, 0, 0, 1.1, 0, 0, 0,
        0), nrow = 5, byrow = TRUE)))
    Condition
      Error:
      ! The hypothesis weight and transition matrix constraints must be dimensionally consistent. The transition matrix constraint has 5 rows and 5 columns while the hypothesis weight vector contains 4 elements.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

# diagnose_hc_tc_consistency() works

    Code
      diagnose_hc_tc_consistency(graph_constraint(c(NA, NA, 0, 0, 0), matrix(c(0, 0.8,
        0.2, 0, 0, NA, 0, 0, NA, 0, 0, 0.9, 0, 0.1, 0, 0.8, NA, NA, 0, 0, 1, 0, 0, 0,
        0), nrow = 5, byrow = TRUE)))
    Message
      
      -- Hypothesis weight and transition matrix constraints consistency: 
      v The hypothesis weights vector has 5 elements.
      v The transition matrix has 5 columns and 5 rows.

---

    Code
      diagnose_hc_tc_consistency(new_graph_constraint(c(NA, NA, 0, 0), matrix(c(0,
        0.8, 0.2, 0, 0, NA, 0, 0, NA, 0, 0, 0.9, 0, 0.1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,
        0), nrow = 5, byrow = TRUE)))
    Message
      
      -- Hypothesis weight and transition matrix constraints consistency: 
      v The hypothesis weights vector has 4 elements.
      x The transition matrix must have 4 rows and 4 columns. It has 5 rows and 5 columns. 

---

    Code
      diagnose_hc_tc_consistency(new_graph_constraint(c(NA, NA, 0, 0, 0), matrix(c(0,
        0.8, 0.2, 0, 0, NA, 0, 0, NA, 0, 0, 0.9, 0, 0.1, 0, 1, 0, 0, 0, 0), nrow = 4,
      byrow = TRUE)))
    Message
      
      -- Hypothesis weight and transition matrix constraints consistency: 
      v The hypothesis weights vector has 5 elements.
      x The transition matrix must have 5 rows and 5 columns. It has 4 rows and 5 columns. 

# assert_graph_constraint_not_full() fails when no optimisable vals

    Code
      assert_graph_constraint_not_full(new_graph_constraint(c(1, 0, 0), rbind(c(0, 1,
        0), c(0.3, 0, 0.7), c(1, 0, 0))))
    Condition
      Error:
      ! No hypothesis or transition weights to optimise.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

# graph_constraint() fails when no optimisable values

    Code
      graph_constraint(c(1, 0, 0), rbind(c(0, 1, 0), c(0.3, 0, 0.7), c(1, 0, 0)),
      diagnose = TRUE)
    Message
      x The `multigrain_graph_constraint` object is complete. There are no hypothesis or transition weights to optimise.
      
      -- Hypothesis weight constraint diagnosis: 
      v All hypothesis weights constraint values are less than or equal to 1.
      v All hypothesis weights constraint values are greater than or equal to 0.
      i The hypothesis weights constraint vector is complete and no weights will be optimised.
      v The sum of the complete hypothesis weight constraint vector is 1.
      
      -- Transition matrix constraint diagnosis: 
      v All values on the transition matrix diagonal are 0.
      v The transition matrix is square.
      v All transition matrix values are less than or equal to 1.
      v All transition matrix values are greater than or equal to 0.
      v All transition matrix rows are either complete or have at least 2 values to optimise.
      v All complete transition matrix rows sum up to 1.
      
      -- Hypothesis weight and transition matrix constraints consistency: 
      v The hypothesis weights vector has 3 elements.
      v The transition matrix has 3 columns and 3 rows.
    Condition
      Error in `graph_constraint()`:
      ! No hypothesis or transition weights to optimise.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

# diagnose_graph_constraint_not_full() reports no optimisable vals

    Code
      diagnose_graph_constraint_not_full(new_graph_constraint(c(1, 0, 0), rbind(c(0,
        1, 0), c(0.3, 0, 0.7), c(1, 0, 0))))
    Message
      x The `multigrain_graph_constraint` object is complete. There are no hypothesis or transition weights to optimise.

