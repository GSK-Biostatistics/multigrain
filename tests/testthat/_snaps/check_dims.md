# check_length

    Code
      check_length("foo", 2)
    Condition
      Error:
      ! `"foo"` has length 1; expected 2.

---

    Code
      check_length(1:3, 2)
    Condition
      Error:
      ! `1:3` has length 3; expected 2.

# check_dim

    Code
      check_dim(matrix(c("foo", "bar")), 2)
    Condition
      Error:
      ! `matrix(c("foo", "bar"))` has dim (2 x 1); expected (2 x 2).

---

    Code
      check_dim(matrix(1:6, ncol = 2), 2)
    Condition
      Error:
      ! `matrix(1:6, ncol = 2)` has dim (3 x 2); expected (2 x 2).

