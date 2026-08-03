# assert_trans_constr() works

    Code
      assert_trans_constr(new_graph_constraint(c(NA, NA, 0, 0), matrix(c(0, 0.8, 0.1,
        0, NA, 0, 0, NA, 0, 1, 0, 0, 1, 0, 0, 0), nrow = 4, byrow = TRUE)))
    Condition
      Error:
      ! The sum of complete transition matrix constraint rows must be 1.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

---

    Code
      assert_trans_constr(new_graph_constraint(c(NA, NA, 0, 0), matrix(c(0, 0.8, NA,
        NA, NA, 0.1, 0, NA, 0, 1, 0, 0, 1, 0, 0, 0), nrow = 4, byrow = TRUE)))
    Condition
      Error:
      ! All values on the transition matrix diagonal must be 0.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

---

    Code
      assert_trans_constr(new_graph_constraint(c(NA, NA, 0, 0, 0), matrix(c(0, 0.8,
        0.2, 0, 0, NA, 0, 0, NA, 0, 0, 0.9, 0, 0.1, 0, 1, 0, 0, 0, 0), nrow = 4,
      byrow = TRUE)))
    Condition
      Error:
      ! The transition matrix is not square.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

---

    Code
      assert_trans_constr(new_graph_constraint(c(NA, NA, 0, 0, 0), matrix(c(0, 0.8,
        0.2, 0, 0, NA, 0, 0, NA, 0, 0, 0.9, 0, 0.1, 0, 0, 0, 0, 0, 1.7, 1, 0, 0, 0, 0),
      nrow = 5, byrow = TRUE)))
    Condition
      Error:
      ! Some transition matrix values are greater than 1.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

---

    Code
      assert_trans_constr(new_graph_constraint(c(NA, NA, 0, 0, 0), matrix(c(0, 0.8,
        0.2, 0, 0, NA, 0, 0, NA, 0, 0, 0.9, 0, 0.1, 0, 0, 0, NA, 0, 0.9, 1, 0, 0, 0,
        0), nrow = 5, byrow = TRUE)))
    Condition
      Error:
      ! At least one incomplete transition matrix row has a single optimisable (i.e. `NA`) value.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

---

    Code
      assert_trans_constr(new_graph_constraint(c(NA, NA, 0, 0), matrix(c(0, 0.8, 0.1,
        0, NA, 0, 0, NA, 0, 1, 0, 0, 1, 0, 0, 0), nrow = 4, byrow = TRUE)))
    Condition
      Error:
      ! The sum of complete transition matrix constraint rows must be 1.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

# assert_trans_constr_diagonal() works

    Code
      assert_trans_constr_diagonal(new_graph_constraint(c(NA, NA, 0, 0), matrix(c(0,
        0.8, NA, NA, NA, 0.1, 0, NA, 0, 1, 0, 0, 1, 0, 0, 0), nrow = 4, byrow = TRUE)))
    Condition
      Error:
      ! All values on the transition matrix diagonal must be 0.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

---

    Code
      assert_trans_constr_diagonal(new_graph_constraint(c(NA, NA, 0, 0), matrix(c(0,
        0.8, NA, NA, NA, NA, 0, NA, 0, 1, 0, 0, 1, 0, 0, 0), nrow = 4, byrow = TRUE)))
    Condition
      Error:
      ! There must not be any `NA`s on the transition matrix diagonal. All values must be 0.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

# assert_trans_constr_square() works

    Code
      assert_trans_constr_square(new_graph_constraint(c(NA, NA, 0, 0, 0), matrix(c(0,
        0.8, 0.2, 0, 0, NA, 0, 0, NA, 0, 0, 0.9, 0, 0.1, 0, 1, 0, 0, 0, 0), nrow = 4,
      byrow = TRUE)))
    Condition
      Error:
      ! The transition matrix is not square.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

# assert_trans_constr_values() works

    Code
      assert_trans_constr_values(new_graph_constraint(c(NA, NA, 0, 0, 0), matrix(c(0,
        0.8, 0.2, 0, 0, NA, 0, 0, NA, 0, 0, 0.9, 0, 0.1, 0, 0, 0, 0, 0, 1.7, 1, 0, 0,
        0, 0), nrow = 5, byrow = TRUE)))
    Condition
      Error:
      ! Some transition matrix values are greater than 1.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

---

    Code
      assert_trans_constr_values(new_graph_constraint(c(NA, NA, 0, 0, 0), matrix(c(0,
        0.8, 0.2, 0, 0, NA, 0, 0, NA, 0, 0, 0.9, 0, 0.1, 0, 0, NA, NA, 0, -0.1, 1, 0,
        0, 0, 0), nrow = 5, byrow = TRUE)))
    Condition
      Error:
      ! Some transition matrix values are less than 0.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

# assert_trans_constr_row_na() works

    Code
      assert_trans_constr_row_na(new_graph_constraint(c(NA, NA, 0, 0, 0), matrix(c(0,
        0.8, 0.2, 0, 0, NA, 0, 0, NA, 0, 0, 0.9, 0, 0.1, 0, 0, 0, NA, 0, 0.9, 1, 0, 0,
        0, 0), nrow = 5, byrow = TRUE)))
    Condition
      Error:
      ! At least one incomplete transition matrix row has a single optimisable (i.e. `NA`) value.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

# assert_trans_constr_row_sums() works

    Code
      assert_trans_constr_row_sums(new_graph_constraint(c(NA, NA, 0, 0), matrix(c(0,
        0.8, 0.1, 0, NA, 0, 0, NA, 0, 1, 0, 0, 1, 0, 0, 0), nrow = 4, byrow = TRUE)))
    Condition
      Error:
      ! The sum of complete transition matrix constraint rows must be 1.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

---

    Code
      assert_trans_constr_row_sums(new_graph_constraint(c(NA, NA, 0, 0), matrix(c(0,
        0.8, 0.2, 0, 0.2, 0, 0, 0.9, 0, 1, 0, 0, 1, 0, 0, 0), nrow = 4, byrow = TRUE)))
    Condition
      Error:
      ! The sum of transition matrix constraint rows cannot be greater than 1.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

---

    Code
      assert_trans_constr_row_sums(new_graph_constraint(c(NA, NA, 0, 0, 0), matrix(c(
        0, 0.8, 0.2, 0, 0, 0.2, 0, NA, 0.9, NA, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0,
        0, 0), nrow = 5, byrow = TRUE)))
    Condition
      Error:
      ! The sum of transition matrix constraint rows cannot be greater than 1.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

---

    Code
      assert_trans_constr_row_sums(new_graph_constraint(c(NA, NA, 0, 0), matrix(c(0,
        0.8, 0.1, 0.1, NA, 0, NA, 1, 0, 1, 0, 0, 1, 0, 0, 0), nrow = 4, byrow = TRUE)))
    Condition
      Error:
      ! At least one incomplete transition matrix row has a sum equal to 1.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

# diagnose_trans_constr_diagonal() works

    Code
      diagnose_trans_constr_diagonal(graph_constraint(c(NA, NA, 0, 0, 0), matrix(c(0,
        0.8, 0.2, 0, 0, NA, 0, 0, NA, 0, 0, 0.9, 0, 0.1, 0, 0.8, NA, NA, 0, 0, 1, 0,
        0, 0, 0), nrow = 5, byrow = TRUE)))
    Message
      v All values on the transition matrix diagonal are 0.

---

    Code
      diagnose_trans_constr_diagonal(new_graph_constraint(c(NA, NA, 0, 0, 0), matrix(
        c(0, 0.8, 0.2, 0, 0, NA, 0, 0, NA, 0, 0, 0.9, 0, 0.1, 0, 1, 0, 0, 0, 0, 0.9,
          0, 0, 0, 0.1), nrow = 5, byrow = TRUE)))
    Message
      x Some values on the transition matrix diagonal are not 0.

---

    Code
      diagnose_trans_constr_diagonal(new_graph_constraint(c(NA, NA, 0, 0, 0), matrix(
        c(0, 0.8, 0.2, 0, 0, NA, 0, 0, NA, 0, 0, 0.9, 0, 0.1, 0, 0.8, 0, NA, NA, 0, 1,
          0, 0, 0, 0), nrow = 5, byrow = TRUE)))
    Message
      x All values on the transition matrix diagonal must be 0. Some values are `NA`.

# diagnose_trans_constr_square() works

    Code
      diagnose_trans_constr_square(graph_constraint(c(NA, NA, 0, 0, 0), matrix(c(0,
        0.8, 0.2, 0, 0, NA, 0, 0, NA, 0, 0, 0.9, 0, 0.1, 0, 0.8, NA, NA, 0, 0, 1, 0,
        0, 0, 0), nrow = 5, byrow = TRUE)))
    Message
      v The transition matrix is square.

---

    Code
      diagnose_trans_constr_square(new_graph_constraint(c(NA, NA, 0, 0, 0), matrix(c(
        0, 0.8, 0.2, 0, 0, NA, 0, 0, NA, 0, 0, 0.9, 0, 0.1, 0, 1, 0, 0, 0, 0), nrow = 4,
      byrow = TRUE)))
    Message
      x The transition matrix is not square.

# diagnose_trans_constr_row_sums() works

    Code
      diagnose_trans_constr_row_sums(new_graph_constraint(c(NA, NA, 0, 0, 0), matrix(
        c(0, 0.8, 0.1, 0, 0, NA, 0, 0, NA, 0, 0, 0.8, 0, 0.2, 0, 0.9, 0, 0.1, 0, 0, 1,
          0, 0, 0, 0), nrow = 5, byrow = TRUE)))
    Message
      x At least one complete transition matrix row does not sum up to 1.
        * Row 1: [0, 0.8, 0.1, 0, 0] has a sum of 0.9.
      v No incomplete transition matrix rows have a sum greater than 1.
      v No incomplete transition matrix rows have a sum equal to 1.

---

    Code
      diagnose_trans_constr_row_sums(new_graph_constraint(c(NA, NA, 0, 0, 0), matrix(
        c(0, 0.8, 0.2, 0, 0, NA, 0, 0, NA, 0, 0, 0.9, 0, 0.2, 0, 0.9, NA, NA, 0, 0, 1,
          0, 0, 0, 0), nrow = 5, byrow = TRUE)))
    Message
      x At least one complete transition matrix row does not sum up to 1.
        * Row 3: [0, 0.9, 0, 0.2, 0] has a sum of 1.1.
      v No incomplete transition matrix rows have a sum greater than 1.
      v No incomplete transition matrix rows have a sum equal to 1.

---

    Code
      diagnose_trans_constr_row_sums(new_graph_constraint(c(NA, NA, 0, 0, 0), matrix(
        c(0, 0.8, 0.2, 0, 0, NA, 0, 0, NA, 0, 0, 0.9, 0, 0.1, 0, 0.9, 0.2, NA, 0, NA,
          1, 0, 0, 0, 0), nrow = 5, byrow = TRUE)))
    Message
      v All complete transition matrix rows sum up to 1.
      x At least one incomplete transition matrix row has a sum greater than 1.
        * Row 4: [0.9, 0.2, NA, 0, NA] has a sum of 1.1.
      v No incomplete transition matrix rows have a sum equal to 1.

---

    Code
      diagnose_trans_constr_row_sums(new_graph_constraint(c(NA, NA, 0, 0, 0), matrix(
        c(0, 0.8, 0.2, 0, 0, NA, 0, 0, NA, 0, 0, 0.9, 0, 0.1, 0, 0.3, NA, NA, 0, 0.7,
          1, 0, 0, 0, 0), nrow = 5, byrow = TRUE)))
    Message
      v All complete transition matrix rows sum up to 1.
      v No incomplete transition matrix rows have a sum greater than 1.
      x At least one incomplete transition matrix row has a sum equal to 1.
        * Row 4: [0.3, NA, NA, 0, 0.7] has a sum of 1.

---

    Code
      diagnose_trans_constr_row_sums(graph_constraint(c(NA, NA, 0, 0, 0), matrix(c(0,
        0.8, 0.2, 0, 0, NA, 0, 0, NA, 0, 0, 0.9, 0, 0.1, 0, 0.3, NA, NA, 0, 0.6, 1, 0,
        0, 0, 0), nrow = 5, byrow = TRUE)))
    Message
      v All complete transition matrix rows sum up to 1.
      v No incomplete transition matrix rows have a sum greater than 1.
      v No incomplete transition matrix rows have a sum equal to 1.

# diagnose_trans_constr_gt1() works as expected

    Code
      diagnose_trans_constr_gt1(new_graph_constraint(c(NA, NA, 0, 0, 0), matrix(c(0,
        0.8, 0.2, 0, 0, NA, 0, 0, NA, 0, 0, 0.9, 0, 0.1, 0, 1.1, NA, NA, 0, 0, 1, 0,
        0, 0, 0), nrow = 5, byrow = TRUE)))
    Message
      x Some transition matrix values are greater than 1.
        * Row 4: [1.1, NA, NA, 0, 0] contains at least one value greater than 1: 1.1.

---

    Code
      diagnose_trans_constr_gt1(graph_constraint(c(NA, NA, 0, 0, 0), matrix(c(0, 0.8,
        0.2, 0, 0, NA, 0, 0, NA, 0, 0, 0.9, 0, 0.1, 0, 0.3, 0, NA, 0, NA, 1, 0, 0, 0,
        0), nrow = 5, byrow = TRUE)))
    Message
      v All transition matrix values are less than or equal to 1.

# diagnose_trans_constr_lt0() works

    Code
      diagnose_trans_constr_lt0(new_graph_constraint(c(NA, NA, 0, 0, 0), matrix(c(0,
        0.8, 0.2, 0, 0, NA, 0, 0, NA, 0, 0, 0.9, 0, 0.1, 0, -0.1, NA, NA, 0, 0, 1, 0,
        0, 0, 0), nrow = 5, byrow = TRUE)))
    Message
      x Some transition matrix values are less than 0.
        * Row 4: [-0.1, NA, NA, 0, 0] contains at least one value less than 0: -0.1.

---

    Code
      diagnose_trans_constr_lt0(graph_constraint(c(NA, NA, 0, 0, 0), matrix(c(0, 0.8,
        0.2, 0, 0, NA, 0, 0, NA, 0, 0, 0.9, 0, 0.1, 0, 0.3, NA, NA, 0, 0, 1, 0, 0, 0,
        0), nrow = 5, byrow = TRUE)))
    Message
      v All transition matrix values are greater than or equal to 0.

# diagnose_trans_constr_row_na() works

    Code
      diagnose_trans_constr_row_na(new_graph_constraint(c(NA, NA, 0, 0, 0), matrix(c(
        0, 0.8, NA, 0, 0, NA, 0, 0, NA, 0, 0, 0.8, 0, 0.2, 0, 0.9, 0, NA, 0, 0, 1, 0,
        0, 0, 0), nrow = 5, byrow = TRUE)))
    Message
      x At least one incomplete transition matrix row has a single optimisable (i.e. `NA`) value.
        * Row 1: [0, 0.8, NA, 0, 0] the optimisable value is effectively equal to
        0.2.
        * Row 4: [0.9, 0, NA, 0, 0] the optimisable value is effectively equal to
        0.1.

---

    Code
      diagnose_trans_constr_row_na(graph_constraint(c(NA, NA, 0, 0, 0), matrix(c(0,
        0.8, 0.2, 0, 0, NA, 0, 0, NA, NA, 0, 0.9, 0, 0.1, 0, 0.9, 0, NA, 0, NA, 1, 0,
        0, 0, 0), nrow = 5, byrow = TRUE)))
    Message
      v All transition matrix rows are either complete or have at least 2 values to optimise.

# diagnose_trans_constr() works

    Code
      diagnose_trans_constr(new_graph_constraint(c(NA, NA, 0, 0), matrix(c(0, 0.8,
        0.1, 0, NA, 0, 0, NA, 0, 1, 0, 0, 1, 0, 0, 0), nrow = 4, byrow = TRUE)))
    Message
      
      -- Transition matrix constraint diagnosis: 
      v All values on the transition matrix diagonal are 0.
      v The transition matrix is square.
      v All transition matrix values are less than or equal to 1.
      v All transition matrix values are greater than or equal to 0.
      v All transition matrix rows are either complete or have at least 2 values to optimise.
      x At least one complete transition matrix row does not sum up to 1.
        * Row 1: [0, 0.8, 0.1, 0] has a sum of 0.9.
      v No incomplete transition matrix rows have a sum greater than 1.
      v No incomplete transition matrix rows have a sum equal to 1.

---

    Code
      diagnose_trans_constr(new_graph_constraint(c(NA, NA, 0, 0), matrix(c(0, 0.8, NA,
        NA, NA, 0.1, 0, NA, 0, 1, 0, 0, 1, 0, 0, 0), nrow = 4, byrow = TRUE)))
    Message
      
      -- Transition matrix constraint diagnosis: 
      x Some values on the transition matrix diagonal are not 0.
      v The transition matrix is square.
      v All transition matrix values are less than or equal to 1.
      v All transition matrix values are greater than or equal to 0.
      v All transition matrix rows are either complete or have at least 2 values to optimise.
      v All complete transition matrix rows sum up to 1.
      v No incomplete transition matrix rows have a sum greater than 1.
      v No incomplete transition matrix rows have a sum equal to 1.

---

    Code
      diagnose_trans_constr(new_graph_constraint(c(NA, NA, 0, 0, 0), matrix(c(0, 0.8,
        0.2, 0, 0, NA, 0, 0, NA, 0, 0, 0.9, 0, 0.1, 0, 1, 0, 0, 0, 0), nrow = 4,
      byrow = TRUE)))
    Message
      
      -- Transition matrix constraint diagnosis: 
      v All values on the transition matrix diagonal are 0.
      x The transition matrix is not square.
      v All transition matrix values are less than or equal to 1.
      v All transition matrix values are greater than or equal to 0.
      v All transition matrix rows are either complete or have at least 2 values to optimise.
      v All complete transition matrix rows sum up to 1.
      v No incomplete transition matrix rows have a sum greater than 1.
      v No incomplete transition matrix rows have a sum equal to 1.

---

    Code
      diagnose_trans_constr(new_graph_constraint(c(NA, NA, 0, 0, 0), matrix(c(0, 0.8,
        0.2, 0, 0, NA, 0, 0, NA, 0, 0, 0.9, 0, 0.1, 0, 0, 0, 0, 0, 1.7, 1, 0, 0, 0, 0),
      nrow = 5, byrow = TRUE)))
    Message
      
      -- Transition matrix constraint diagnosis: 
      v All values on the transition matrix diagonal are 0.
      v The transition matrix is square.
      x Some transition matrix values are greater than 1.
        * Row 4: [0, 0, 0, 0, 1.7] contains at least one value greater than 1: 1.7.
      v All transition matrix values are greater than or equal to 0.
      v All transition matrix rows are either complete or have at least 2 values to optimise.
      x At least one complete transition matrix row does not sum up to 1.
        * Row 4: [0, 0, 0, 0, 1.7] has a sum of 1.7.
      v No incomplete transition matrix rows have a sum greater than 1.
      v No incomplete transition matrix rows have a sum equal to 1.

---

    Code
      diagnose_trans_constr(new_graph_constraint(c(NA, NA, 0, 0, 0), matrix(c(0, 0.8,
        0.2, 0, 0, NA, 0, 0, NA, 0, 0, 0.9, 0, 0.1, 0, 0, 0, NA, 0, 0.9, 1, 0, 0, 0,
        0), nrow = 5, byrow = TRUE)))
    Message
      
      -- Transition matrix constraint diagnosis: 
      v All values on the transition matrix diagonal are 0.
      v The transition matrix is square.
      v All transition matrix values are less than or equal to 1.
      v All transition matrix values are greater than or equal to 0.
      x At least one incomplete transition matrix row has a single optimisable (i.e. `NA`) value.
        * Row 4: [0, 0, NA, 0, 0.9] the optimisable value is effectively equal to
        0.1.
      v All complete transition matrix rows sum up to 1.
      v No incomplete transition matrix rows have a sum greater than 1.
      v No incomplete transition matrix rows have a sum equal to 1.

---

    Code
      diagnose_trans_constr(new_graph_constraint(c(NA, NA, 0, 0), matrix(c(0, 0.8,
        0.1, 0, NA, 0, 0, NA, 0, 1, 0, 0, 1, 0, 0, 0), nrow = 4, byrow = TRUE)))
    Message
      
      -- Transition matrix constraint diagnosis: 
      v All values on the transition matrix diagonal are 0.
      v The transition matrix is square.
      v All transition matrix values are less than or equal to 1.
      v All transition matrix values are greater than or equal to 0.
      v All transition matrix rows are either complete or have at least 2 values to optimise.
      x At least one complete transition matrix row does not sum up to 1.
        * Row 1: [0, 0.8, 0.1, 0] has a sum of 0.9.
      v No incomplete transition matrix rows have a sum greater than 1.
      v No incomplete transition matrix rows have a sum equal to 1.

---

    Code
      diagnose_trans_constr(new_graph_constraint(c(NA, NA, 0, 0), matrix(c(0, 0.8,
        0.1, 0.1, NA, 0, 0, NA, 0, 1, 0, 0, 1, 0, 0, 0), nrow = 4, byrow = TRUE)))
    Message
      
      -- Transition matrix constraint diagnosis: 
      v All values on the transition matrix diagonal are 0.
      v The transition matrix is square.
      v All transition matrix values are less than or equal to 1.
      v All transition matrix values are greater than or equal to 0.
      v All transition matrix rows are either complete or have at least 2 values to optimise.
      v All complete transition matrix rows sum up to 1.
      v No incomplete transition matrix rows have a sum greater than 1.
      v No incomplete transition matrix rows have a sum equal to 1.

---

    Code
      diagnose_trans_constr(new_graph_constraint(c(NA, NA, 0, 0, 0), matrix(c(0, 0.8,
        0.3, 0, 0, NA, 0, -0.1, NA, 1.2, 0, 0.9, NA, 0.1, 0, 1, 0, 0, NA, NA), nrow = 4,
      byrow = TRUE)))
    Message
      
      -- Transition matrix constraint diagnosis: 
      x All values on the transition matrix diagonal must be 0. Some values are `NA`.
      x The transition matrix is not square.
      x Some transition matrix values are greater than 1.
        * Row 2: [NA, 0, -0.1, NA, 1.2] contains at least one value greater than 1:
        1.2.
      x Some transition matrix values are less than 0.
        * Row 2: [NA, 0, -0.1, NA, 1.2] contains at least one value less than 0:
        -0.1.
      x At least one incomplete transition matrix row has a single optimisable (i.e. `NA`) value.
        * Row 3: [0, 0.9, NA, 0.1, 0] the optimisable value is effectively equal to
        0.
      x At least one complete transition matrix row does not sum up to 1.
        * Row 1: [0, 0.8, 0.3, 0, 0] has a sum of 1.1.
      x At least one incomplete transition matrix row has a sum greater than 1.
        * Row 2: [NA, 0, -0.1, NA, 1.2] has a sum of 1.1.
      x At least one incomplete transition matrix row has a sum equal to 1.
        * Row 3: [0, 0.9, NA, 0.1, 0] has a sum of 1.
        * Row 4: [1, 0, 0, NA, NA] has a sum of 1.

# validate_graph_constraint() picks up transition matrix issues

    Code
      validate_graph_constraint(new_graph_constraint(c(NA, NA, 0, 0), matrix(c(0, 0.8,
        0.1, 0, NA, 0, 0, NA, 0, 1, 0, 0, 1, 0, 0, 0), nrow = 4, byrow = TRUE)),
      diagnose = TRUE)
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
      x At least one complete transition matrix row does not sum up to 1.
        * Row 1: [0, 0.8, 0.1, 0] has a sum of 0.9.
      v No incomplete transition matrix rows have a sum greater than 1.
      v No incomplete transition matrix rows have a sum equal to 1.
      
      -- Hypothesis weight and transition matrix constraints consistency: 
      v The hypothesis weights vector has 4 elements.
      v The transition matrix has 4 columns and 4 rows.
    Condition
      Error:
      ! The sum of complete transition matrix constraint rows must be 1.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

---

    Code
      validate_graph_constraint(new_graph_constraint(c(NA, NA, 0, 0), matrix(c(0, 0.8,
        0.1, 0.1, NA, 0, 0, NA, 0, 1, 0, 0, 1, 0, 0, 0), nrow = 4, byrow = TRUE)))
    Output
      <multigrain_graph_constraint>
      Constraints on hypothesis weights:
      [1] NA NA  0  0
      
      Constraints on transition matrix:
           [,1] [,2] [,3] [,4]
      [1,]    0  0.8  0.1  0.1
      [2,]   NA  0.0  0.0   NA
      [3,]    0  1.0  0.0  0.0
      [4,]    1  0.0  0.0  0.0

# assert_trans_constr_values() with tolerance

    Code
      assert_trans_constr_values(new_graph_constraint(c(NA, NA, 0, 0, 0), matrix(c(0,
        0.8, 0.2, 0, 0, NA, 0, 0, NA, 0, 0, 0.9, 0, 0.1, 0, 0, 0, 0, 0, 1 + 1e-12, 1,
        0, 0, 0, 0), nrow = 5, byrow = TRUE)))

---

    Code
      assert_trans_constr_values(new_graph_constraint(c(NA, NA, 0, 0, 0), matrix(c(0,
        0.8, 0.2, 0, 0, NA, 0, 0, NA, 0, 0, 0.9, 0, 0.1, 0, 0, 0, 0, 0, 1 + 1e-12, 1,
        0, 0, 0, 0), nrow = 5, byrow = TRUE), tolerance = 0))
    Condition
      Error:
      ! Some transition matrix values are greater than 1.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

# trans_constraint row sum validation with small differences

    Code
      diagnose_trans_constr_row_sums(new_graph_constraint(trans_constraint = test_trans_constraint))
    Message
      v All complete transition matrix rows sum up to 1.

---

    Code
      diagnose_trans_constr_row_sums(new_graph_constraint(trans_constraint = test_trans_constraint,
        tolerance = 0))
    Message
      x At least one complete transition matrix row does not sum up to 1.
        * Row 2: [0.18..., 0, 0.34..., 0.25..., 0.21...] has a sum of 1.000000000004.
        * Row 3: [0.45..., 0.11..., 0, 0.39..., 0.03...] has a sum of 1.000000000004.
        * Row 4: [0.25..., 0.27..., 0.12..., 0, 0.33...] has a sum of 1.000000000004.
        * Row 5: [0.22..., 0.42..., 0.16..., 0.17..., 0] has a sum of 1.000000000004.

---

    Code
      graph_constraint(trans_constraint = test_trans_constraint, tolerance = 0)
    Condition
      Error in `graph_constraint()`:
      ! The sum of transition matrix constraint rows cannot be greater than 1.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

---

    Code
      assert_trans_constr_row_sums(new_graph_constraint(trans_constraint = test_trans_constraint,
        tolerance = 0))
    Condition
      Error:
      ! The sum of transition matrix constraint rows cannot be greater than 1.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

---

    Code
      diagnose_trans_constr_row_sums(new_graph_constraint(trans_constraint = test_trans_constraint,
        tolerance = 0))
    Message
      x At least one complete transition matrix row does not sum up to 1.
        * Row 2: [0.18..., 0, 0.34..., 0.25..., 0.21...] has a sum of 1.000000000004.
        * Row 3: [0.45..., 0.11..., 0, 0.39..., 0.03...] has a sum of 1.000000000004.
        * Row 4: [0.25..., 0.27..., 0.12..., 0, 0.33...] has a sum of 1.000000000004.
        * Row 5: [0.22..., 0.42..., 0.16..., 0.17..., 0] has a sum of 1.000000000004.

# incomplete trans_constraint rows with tolerance

    Code
      assert_trans_constr_row_sums(new_graph_constraint(trans_constraint = tc))
    Condition
      Error:
      ! At least one incomplete transition matrix row has a sum equal to 1.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

---

    Code
      diagnose_trans_constr_row_sums(new_graph_constraint(trans_constraint = tc))
    Message
      v All complete transition matrix rows sum up to 1.
      v No incomplete transition matrix rows have a sum greater than 1.
      x At least one incomplete transition matrix row has a sum equal to 1.
        * Row 1: [0, 0.07..., NA, 0.92..., NA] has a sum of 1.000000000002.
        * Row 3: [NA, NA, 0, 0.66..., 0.33...] has a sum of 1.000000000002.
        * Row 4: [0.02..., 0.97..., NA, 0, NA] has a sum of 1.000000000002.
        * Row 5: [0.85..., NA, 0.14..., NA, 0] has a sum of 1.000000000002.

---

    Code
      assert_trans_constr_row_sums(new_graph_constraint(trans_constraint = tc,
        tolerance = 0))
    Condition
      Error:
      ! The sum of transition matrix constraint rows cannot be greater than 1.
      i For a more detailed diagnosis run `graph_constraint()` with `diagnose = TRUE`.

# diagnose_trans_constr_gt1() with tolerance

    Code
      diagnose_trans_constr_gt1(new_graph_constraint(c(NA, NA, 0, 0, 0), matrix(c(0,
        0.8, 0.2, 0, 0, NA, 0, 0, NA, 0, 0, 0.9, 0, 0.1, 0, 1 + 1e-12, NA, NA, 0, 0,
        1, 0, 0, 0, 0), nrow = 5, byrow = TRUE)))
    Message
      v All transition matrix values are less than or equal to 1.

---

    Code
      diagnose_trans_constr_gt1(new_graph_constraint(c(NA, NA, 0, 0, 0), matrix(c(0,
        0.8, 0.2, 0, 0, NA, 0, 0, NA, 0, 0, 0.9, 0, 0.1, 0, 1 + 1e-12, NA, NA, 0, 0,
        1, 0, 0, 0, 0), nrow = 5, byrow = TRUE), tolerance = 0))
    Message
      x Some transition matrix values are greater than 1.
        * Row 4: [1.00..., NA, NA, 0, 0] contains at least one value greater than 1:
        1.000000000001.

