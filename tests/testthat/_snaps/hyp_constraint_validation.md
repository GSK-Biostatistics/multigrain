# assert_hyp_constraint_sum() works

    Code
      assert_hyp_constraint_sum(new_graph_constraint(c(0.8, 0.3, 0, 0, 0)))
    Condition
      Error:
      ! The sum of the hypothesis weight constraint vector cannot be greater than 1. It is 1.1.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

---

    Code
      assert_hyp_constraint_sum(new_graph_constraint(c(0.7, 0.2, 0, 0, 0)))
    Condition
      Error:
      ! The sum of a complete hypothesis weight constraint vector must be 1. It is 0.9.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

---

    Code
      assert_hyp_constraint_sum(new_graph_constraint(c(0.7, 0.3, NA, NA, 0)))
    Condition
      Error:
      ! The sum of an incomplete hypothesis weight constraint vector must not be equal to 1.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

# sum(hc) == 1 is tolerant of small floating point differences

    Code
      diagnose_hyp_constraint_sum(new_graph_constraint(b))
    Message
      v The sum of the complete hypothesis weight constraint vector is 1.

---

    Code
      graph_constraint(b, tolerance = 0)
    Condition
      Error in `graph_constraint()`:
      ! The sum of the hypothesis weight constraint vector cannot be greater than 1. It is 1.000000000005.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

---

    Code
      assert_hyp_constraint_sum(new_graph_constraint(b, tolerance = 0))
    Condition
      Error:
      ! The sum of the hypothesis weight constraint vector cannot be greater than 1. It is 1.000000000005.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

---

    Code
      diagnose_hyp_constraint_sum(new_graph_constraint(b, tolerance = 0))
    Message
      x The sum of the complete hypothesis weights constraint vector is 1.000000000005. It should be 1.

# assert_hyp_constraint_values() works

    Code
      assert_hyp_constraint_values(new_graph_constraint(c(NA, NA, 0, 0.2, 1.9)))
    Condition
      Error:
      ! Values in the hypothesis weight constraint vector cannot be greater than 1.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

---

    Code
      assert_hyp_constraint_values(new_graph_constraint(c(NA, NA, 0, 0.2, -0.1)))
    Condition
      Error:
      ! Values in the hypothesis weights constraint vector cannot be less than 0.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

# assert_hyp_constraint_na() works

    Code
      assert_hyp_constraint_na(new_graph_constraint(c(NA, 0, 0, 0.2, 0)))
    Condition
      Error:
      ! An incomplete hypothesis weight constraint vector cannot have a single optimisable (i.e. `NA`) value.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

# assert_hyp_constraint() works

    Code
      assert_hyp_constraint(new_graph_constraint(c(0.8, 0.3, 0, 0, 0)))
    Condition
      Error:
      ! The sum of the hypothesis weight constraint vector cannot be greater than 1. It is 1.1.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

---

    Code
      assert_hyp_constraint(new_graph_constraint(c(NA, NA, 0, 0.2, 1.9)))
    Condition
      Error:
      ! Values in the hypothesis weight constraint vector cannot be greater than 1.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

---

    Code
      assert_hyp_constraint(new_graph_constraint(c(NA, 0, 0, 0.2, 0)))
    Condition
      Error:
      ! An incomplete hypothesis weight constraint vector cannot have a single optimisable (i.e. `NA`) value.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

# diagnose_hyp_constraint_sum() works

    Code
      diagnose_hyp_constraint_sum(new_graph_constraint(c(NA, NA, 0, 0, 0)))
    Message
      v The sum of the incomplete hypothesis weights constraint vector is less than 1.

---

    Code
      diagnose_hyp_constraint_sum(new_graph_constraint(c(NA, NA, 0, 0.2, 0.9)))
    Message
      x The sum of the incomplete hypothesis weight constraint vector cannot be greater than 1. It is 1.1.

---

    Code
      diagnose_hyp_constraint_sum(new_graph_constraint(c(NA, NA, 0, 0.2, 0.8)))
    Message
      x The sum of the incomplete hypothesis weight constraint vector is 1.

---

    Code
      diagnose_hyp_constraint_sum(new_graph_constraint(c(0, 0, 0, 0.2, 0.5)))
    Message
      x The sum of the complete hypothesis weights constraint vector is 0.7. It should be 1.

---

    Code
      diagnose_hyp_constraint_sum(new_graph_constraint(c(0, 0.2, 0.1, 0.2, 0.5)))
    Message
      v The sum of the complete hypothesis weight constraint vector is 1.

# diagnose_hyp_constraint_gt1() works

    Code
      diagnose_hyp_constraint_gt1(new_graph_constraint(c(NA, NA, 0, 0, 0)))
    Message
      v All hypothesis weights constraint values are less than or equal to 1.

---

    Code
      diagnose_hyp_constraint_gt1(new_graph_constraint(c(NA, NA, 0, 0.2, 1.9)))
    Message
      x Some hypothesis weight constraint values are greater than 1.

# diagnose_hyp_constraint_lt0() works

    Code
      diagnose_hyp_constraint_lt0(new_graph_constraint(c(NA, NA, 0, 0, 0)))
    Message
      v All hypothesis weights constraint values are greater than or equal to 0.

---

    Code
      diagnose_hyp_constraint_lt0(new_graph_constraint(c(NA, NA, 0, 0.2, -0.1)))
    Message
      x Some hypothesis weight constraint values are less than 0.

# diagnose_hyp_constraint_na() works

    Code
      diagnose_hyp_constraint_na(new_graph_constraint(c(NA, NA, 0, 0, 0)))
    Message
      v The hypothesis weights constraint vector is correctly defined for optimising 2 weights.

---

    Code
      diagnose_hyp_constraint_na(new_graph_constraint(c(NA, 0, 0, 0.2, 0.7)))
    Message
      x A single hypothesis weight constraint cannot be optimised as weights must add up to 1.

---

    Code
      diagnose_hyp_constraint_na(new_graph_constraint(c(0.1, 0, 0, 0.2, 0.7)))
    Message
      i The hypothesis weights constraint vector is complete and no weights will be optimised.

# diagnose_hyp_constraint() works

    Code
      diagnose_hyp_constraint(new_graph_constraint(c(NA, NA, 0, 0, 0)))
    Message
      
      -- Hypothesis weight constraint diagnosis: 
      v All hypothesis weights constraint values are less than or equal to 1.
      v All hypothesis weights constraint values are greater than or equal to 0.
      v The hypothesis weights constraint vector is correctly defined for optimising 2 weights.
      v The sum of the incomplete hypothesis weights constraint vector is less than 1.

---

    Code
      diagnose_hyp_constraint(new_graph_constraint(c(NA, 0, 0.9, 1.1, -0.1)))
    Message
      
      -- Hypothesis weight constraint diagnosis: 
      x Some hypothesis weight constraint values are greater than 1.
      x Some hypothesis weight constraint values are less than 0.
      x A single hypothesis weight constraint cannot be optimised as weights must add up to 1.
      x The sum of the incomplete hypothesis weight constraint vector cannot be greater than 1. It is 1.9.

# diagnose_hyp_constraint() works with empty vector

    Code
      diagnose_hyp_constraint(new_graph_constraint())
    Message
      
      -- Hypothesis weight constraint diagnosis: 
      v All hypothesis weights constraint values are less than or equal to 1.
      v All hypothesis weights constraint values are greater than or equal to 0.
      i The hypothesis weights constraint vector is complete and no weights will be optimised.
      v The sum of the incomplete hypothesis weights constraint vector is less than 1.

# validate_graph_constraint() picks up hypothesis weight issues

    Code
      validate_graph_constraint(new_graph_constraint(c(0.8, 0.3, 0, 0, 0), tc))
    Condition
      Error:
      ! The sum of the hypothesis weight constraint vector cannot be greater than 1. It is 1.1.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

---

    Code
      validate_graph_constraint(new_graph_constraint(c(NA, NA, 0, 0.2, 1.9), tc))
    Condition
      Error:
      ! Values in the hypothesis weight constraint vector cannot be greater than 1.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

---

    Code
      validate_graph_constraint(new_graph_constraint(c(NA, 0, 0, 0.2, 0), tc))
    Condition
      Error:
      ! An incomplete hypothesis weight constraint vector cannot have a single optimisable (i.e. `NA`) value.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

---

    Code
      validate_graph_constraint(graph_constraint(c(NA, NA, 0, 0, 0), tc), diagnose = TRUE)
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
         H1 H2 H3 H4  H5
      H1  0  0 NA NA 0.5
      H2  0  0  1  0 0.0
      H3  0  0  0  1 0.0
      H4  0  0  0  0 1.0
      H5  1  0  0  0 0.0

---

    Code
      validate_graph_constraint(new_graph_constraint(c(NA, 0, 0.9, 1.1, -0.1), tc),
      diagnose = TRUE)
    Message
      
      -- Hypothesis weight constraint diagnosis: 
      x Some hypothesis weight constraint values are greater than 1.
      x Some hypothesis weight constraint values are less than 0.
      x A single hypothesis weight constraint cannot be optimised as weights must add up to 1.
      x The sum of the incomplete hypothesis weight constraint vector cannot be greater than 1. It is 1.9.
      
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
    Condition
      Error:
      ! Values in the hypothesis weight constraint vector cannot be greater than 1.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

# assert_hyp_constraint_values() with tolerance

    Code
      assert_hyp_constraint_values(new_graph_constraint(c(NA, NA, 0, 0, 1 + 1e-12)))

---

    Code
      assert_hyp_constraint_values(new_graph_constraint(c(NA, NA, 0, 0, 1 + 1e-12),
      tolerance = 0))
    Condition
      Error:
      ! Values in the hypothesis weight constraint vector cannot be greater than 1.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

# diagnose_hyp_constraint_gt1() with tolerance

    Code
      diagnose_hyp_constraint_gt1(new_graph_constraint(c(NA, NA, 0, 0.2, 1 + 1e-12)))
    Message
      v All hypothesis weights constraint values are less than or equal to 1.

---

    Code
      diagnose_hyp_constraint_gt1(new_graph_constraint(c(NA, NA, 0, 0.2, 1 + 1e-12),
      tolerance = 0))
    Message
      x Some hypothesis weight constraint values are greater than 1.

