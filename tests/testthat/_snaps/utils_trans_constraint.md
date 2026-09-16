# print_offending_row

    Code
      print_offending_row(1, matrix(1:9, nrow = 3), type = "sum")
    Message
      * Row 1: [1, 4, 7] has a sum of 12.

---

    Code
      print_offending_row(1, matrix(1:9, nrow = 3), type = "value", ref = 0)
    Message
      * Row 1: [1, 4, 7] contains at least one value less than 0: .

---

    Code
      print_offending_row(1, matrix(1:9, nrow = 3), type = "value", ref = 1)
    Message
      * Row 1: [1, 4, 7] contains at least one value greater than 1: 4 and 7.

---

    Code
      print_offending_row(1, matrix(1:9, nrow = 3), type = "na")
    Message
      * Row 1: [1, 4, 7] the optimisable value is effectively equal to -11.

# offending_rows_bullets

    Code
      offending_rows_bullets(c(1, 2), matrix(1:9, nrow = 3), type = "sum")
    Message
          * Row 1: [1, 4, 7] has a sum of 12.
          * Row 2: [2, 5, 8] has a sum of 15.

