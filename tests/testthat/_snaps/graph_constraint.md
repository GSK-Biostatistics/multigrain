# new_graph_constraint() works

    Code
      new_graph_constraint(hyp_constraint = c(NA, NA, 0, 0), trans_constraint = matrix(
        c(NA, NA, 0, 0, NA, NA, 0, 0, NA, NA, 0, 0, NA, NA, 0, 0), nrow = 4, byrow = TRUE),
      names = letters[1:4])
    Output
      <multigrain_graph_constraint>
      Constraints on hypothesis weights:
       a  b  c  d 
      NA NA  0  0 
      
      Constraints on transition matrix:
         a  b c d
      a NA NA 0 0
      b NA NA 0 0
      c NA NA 0 0
      d NA NA 0 0

---

    Code
      new_graph_constraint(hyp_constraint = c(NA, NA, 0, 0), trans_constraint = matrix(
        c(NA, NA, 0, 0, NA, NA, 0, 0, NA, NA, 0, 0, NA, NA, 0, 0), nrow = 4, byrow = TRUE))
    Output
      <multigrain_graph_constraint>
      Constraints on hypothesis weights:
      [1] NA NA  0  0
      
      Constraints on transition matrix:
           [,1] [,2] [,3] [,4]
      [1,]   NA   NA    0    0
      [2,]   NA   NA    0    0
      [3,]   NA   NA    0    0
      [4,]   NA   NA    0    0

# new_graph_constraint() complains with undesired inputs

    Code
      new_graph_constraint(hyp_constraint = "a")
    Condition
      Error in `new_graph_constraint()`:
      ! is.double(hyp_constraint) is not TRUE

---

    Code
      new_graph_constraint(hyp_constraint = c(NA, NA, 0, 0, 0), trans_constraint = "b")
    Condition
      Error in `new_graph_constraint()`:
      ! is.double(trans_constraint) is not TRUE

---

    Code
      new_graph_constraint(hyp_constraint = c(NA, NA, 0, 0, 0), trans_constraint = c(
        NA, NA, 0, 0, 0))
    Condition
      Error in `new_graph_constraint()`:
      ! is.matrix(trans_constraint) is not TRUE

---

    Code
      new_graph_constraint(hyp_constraint = c(NA, NA, 0, 0, 0), trans_constraint = matrix(
        c(NA, NA, 0, 0), nrow = 2, byrow = TRUE), names = c(NA, NA, 0, 0, 0))
    Condition
      Error in `new_graph_constraint()`:
      ! is.character(names) is not TRUE

---

    Code
      new_graph_constraint(hyp_constraint = c(NA, NA, 0, 0), trans_constraint = matrix(
        c(0, NA, NA, 0, NA, 0, NA, 0, NA, NA, 0, 0, NA, NA, 0, 0), nrow = 4, byrow = TRUE),
      tolerance = "a")
    Condition
      Error in `new_graph_constraint()`:
      ! is.double(tolerance) is not TRUE

# graph_constraint() when one of the inputs is NULL

    Code
      graph_constraint(hyp_constraint = c(NA, 0.4, NA))
    Output
      <multigrain_graph_constraint>
      Constraints on hypothesis weights:
       H1  H2  H3 
       NA 0.4  NA 
      
      Constraints on transition matrix:
         H1 H2 H3
      H1  0 NA NA
      H2 NA  0 NA
      H3 NA NA  0

---

    Code
      graph_constraint(trans_constraint = matrix(c(0, NA, NA, 0, NA, 0, NA, 0, NA, NA,
        0, 0, NA, NA, 0, 0), nrow = 4, byrow = TRUE))
    Output
      <multigrain_graph_constraint>
      Constraints on hypothesis weights:
      H1 H2 H3 H4 
      NA NA NA NA 
      
      Constraints on transition matrix:
         H1 H2 H3 H4
      H1  0 NA NA  0
      H2 NA  0 NA  0
      H3 NA NA  0  0
      H4 NA NA  0  0

# graph_constraint() errors with both inputs NULL

    Code
      graph_constraint()
    Condition
      Error in `graph_constraint()`:
      ! `hyp_constraint` and `trans_constraint` cannot both be `NULL` at the same time. At least one must be supplied.

# graph_constraint() errors trans_constraint not numeric

    Code
      graph_constraint(trans_constraint = matrix(letters[1:16], nrow = 4, byrow = TRUE))
    Condition
      Error in `graph_constraint()`:
      ! `trans_constraint` must be a double matrix, not a character matrix.

# graph_constraint complains when anything is passed via `...`

    Code
      graph_constraint(hyp_constraint = c(NA, NA, 0, 0), trans_constraint = matrix(c(
        0, NA, NA, 0, NA, 0, NA, 0, NA, NA, 0, 0, NA, NA, 0, 0), nrow = 4, byrow = TRUE),
      c("a", "b", "c", "d"))
    Condition
      Error in `graph_constraint()`:
      ! `...` must be empty.
      x Problematic argument:
      * ..1 = c("a", "b", "c", "d")
      i Did you forget to name an argument?

# graph_constraint() errors when tolerance not positive numeric

    Code
      graph_constraint(hyp_constraint = c(NA, NA, 0, 0), trans_constraint = matrix(c(
        0, NA, NA, 0, NA, 0, NA, 0, NA, NA, 0, 0, NA, NA, 0, 0), nrow = 4, byrow = TRUE),
      tolerance = "a")
    Condition
      Error in `graph_constraint()`:
      ! `tolerance` must be a number, not the string "a".

---

    Code
      graph_constraint(hyp_constraint = c(NA, NA, 0, 0), trans_constraint = matrix(c(
        0, NA, NA, 0, NA, 0, NA, 0, NA, NA, 0, 0, NA, NA, 0, 0), nrow = 4, byrow = TRUE),
      tolerance = -1e-06)
    Condition
      Error in `graph_constraint()`:
      ! `tolerance` must be a number larger than or equal to 0, not the number 0.

# graph_constraint_free() works

    Code
      graph_constraint_free(4)
    Output
      <multigrain_graph_constraint>
      Constraints on hypothesis weights:
      H1 H2 H3 H4 
      NA NA NA NA 
      
      Constraints on transition matrix:
         H1 H2 H3 H4
      H1  0 NA NA NA
      H2 NA  0 NA NA
      H3 NA NA  0 NA
      H4 NA NA NA  0

# graph_constraint_free() works with 2 hypotheses

    Code
      graph_constraint_free(2)
    Output
      <multigrain_graph_constraint>
      Constraints on hypothesis weights:
      H1 H2 
      NA NA 
      
      Constraints on transition matrix:
         H1 H2
      H1  0  1
      H2  1  0

# graph_constraint: users can update hyp_constraint

    Code
      gc$hyp_constraint <- "A"
    Condition
      Error in `graph_constraint()`:
      ! `hyp_constraint` must be a double or `NA`, not the string "A".

---

    Code
      gc[["hyp_constraint"]] <- c(1, NA, 0, 0)
    Condition
      Error in `graph_constraint()`:
      ! An incomplete hypothesis weight constraint vector cannot have a single optimisable (i.e. `NA`) value.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

---

    Code
      gc["hyp_constraint", tolerance = 0] <- c(1 + 1e-11, 0, 0, 0)
    Condition
      Error in `graph_constraint()`:
      ! Values in the hypothesis weight constraint vector cannot be greater than 1.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

---

    Code
      gc["hyp_constraint"] <- c(1 + 1e-12, 0, 0, 0)

---

    Code
      gc["hyp_constraint", tolerance = 0, diagnose = TRUE] <- c(1 + 1e-11, 0, 0, 0)
    Message
      
      -- Hypothesis weight constraint diagnosis: 
      x Some hypothesis weight constraint values are greater than 1.
      v All hypothesis weights constraint values are greater than or equal to 0.
      i The hypothesis weights constraint vector is complete and no weights will be optimised.
      x The sum of the complete hypothesis weights constraint vector is 1.00000000001. It should be 1.
      
      -- Transition matrix constraint diagnosis: 
      v All values on the transition matrix diagonal are 0.
      v The transition matrix is square.
      v All transition matrix values are less than or equal to 1.
      v All transition matrix values are greater than or equal to 0.
      v All transition matrix rows are either complete or have at least 2 values to optimise.
      v No incomplete transition matrix rows have a sum greater than 1.
      v No incomplete transition matrix rows have a sum equal to 1.
      
      -- Hypothesis weight and transition matrix constraints consistency: 
      v The hypothesis weights vector has 4 elements.
      v The transition matrix has 4 columns and 4 rows.
    Condition
      Error in `graph_constraint()`:
      ! Values in the hypothesis weight constraint vector cannot be greater than 1.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

---

    Code
      gc["hyp_constraint", tolerance = 1e-11, diagnose = TRUE] <- c(1 + 1e-12, 0, 0,
      0)
    Message
      
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
      v No incomplete transition matrix rows have a sum greater than 1.
      v No incomplete transition matrix rows have a sum equal to 1.
      
      -- Hypothesis weight and transition matrix constraints consistency: 
      v The hypothesis weights vector has 4 elements.
      v The transition matrix has 4 columns and 4 rows.

# graph_constraint: users can update trans_constraint

    Code
      gc$trans_constraint <- "A"
    Condition
      Error in `graph_constraint()`:
      ! `trans_constraint` must be a double matrix, not the string "A".

---

    Code
      gc[["trans_constraint"]] <- new_tc

---

    Code
      gc[["trans_constraint", tolerance = 0]] <- tolerance_tc
    Condition
      Error in `graph_constraint()`:
      ! Some transition matrix values are greater than 1.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

---

    Code
      gc[["trans_constraint", tolerance = 1e-11]] <- tolerance_tc

---

    Code
      gc[["trans_constraint", tolerance = 0, diagnose = TRUE]] <- tolerance_tc
    Message
      
      -- Hypothesis weight constraint diagnosis: 
      v All hypothesis weights constraint values are less than or equal to 1.
      v All hypothesis weights constraint values are greater than or equal to 0.
      v The hypothesis weights constraint vector is correctly defined for optimising 2 weights.
      v The sum of the incomplete hypothesis weights constraint vector is less than 1.
      
      -- Transition matrix constraint diagnosis: 
      v All values on the transition matrix diagonal are 0.
      v The transition matrix is square.
      x Some transition matrix values are greater than 1.
        * Row 1: [0, 0, 0, 1.00...] contains at least one value greater than 1:
        1.000000000001.
      v All transition matrix values are greater than or equal to 0.
      v All transition matrix rows are either complete or have at least 2 values to optimise.
      x At least one complete transition matrix row does not sum up to 1.
        * Row 1: [0, 0, 0, 1.00...] has a sum of 1.000000000001.
      v No incomplete transition matrix rows have a sum greater than 1.
      v No incomplete transition matrix rows have a sum equal to 1.
      
      -- Hypothesis weight and transition matrix constraints consistency: 
      v The hypothesis weights vector has 4 elements.
      v The transition matrix has 4 columns and 4 rows.
    Condition
      Error in `graph_constraint()`:
      ! Some transition matrix values are greater than 1.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

---

    Code
      gc[["trans_constraint", tolerance = 1e-11, diagnose = TRUE]] <- tolerance_tc
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
      v The hypothesis weights vector has 4 elements.
      v The transition matrix has 4 columns and 4 rows.

# graph_constraint: update incoming names are preferred

    Code
      gc[["trans_constraint"]] <- new_tc

# graph_constraint print and summary methods

    Code
      print(gc)
    Output
      <multigrain_graph_constraint>
      Constraints on hypothesis weights:
      H1 H2 H3 H4 
      NA NA  0  0 
      
      Constraints on transition matrix:
         H1 H2 H3 H4
      H1  0 NA NA  0
      H2 NA  0 NA  0
      H3 NA NA  0  0
      H4 NA NA  0  0

---

    Code
      summary(gc)
    Output
      
      Graph constraints:
      Constraints on hypothesis weights:
      H1 H2 H3 H4 
      NA NA  0  0 
      
      Constraints on transition matrix:
         H1 H2 H3 H4
      H1  0 NA NA  0
      H2 NA  0 NA  0
      H3 NA NA  0  0
      H4 NA NA  0  0

# set methods inherit the original tolerance if unspecified

    Code
      graph_constraint(hyp_constraint = c(NA, NA, 0, 0), tolerance = "a")
    Condition
      Error in `graph_constraint()`:
      ! `tolerance` must be a number, not the string "a".

---

    Code
      gc[["hyp_constraint", tolerance = 1e-05]] <- c(1 + 0.001, 0, 0, 0)
    Condition
      Error in `graph_constraint()`:
      ! Values in the hypothesis weight constraint vector cannot be greater than 1.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

---

    Code
      gc[["hyp_constraint", tolerance = 0]] <- c(1, 0, 0, 0)
    Condition
      Error in `graph_constraint()`:
      ! The sum of transition matrix constraint rows cannot be greater than 1.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

