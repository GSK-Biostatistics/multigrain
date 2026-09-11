# bullets_with_header

    Code
      bullets_with_header("foo", list(x = 1, y = 2))
    Output
      foo
      * x: 1
      * y: 2

---

    Code
      bullets_with_header("foo", list())
    Output
      NULL

# normalise_sum complains when anything is passed via `...`

    Code
      normalise_sum(x, 1L)
    Condition
      Error in `normalise_sum()`:
      ! `...` must be empty.
      x Problematic argument:
      * ..1 = 1L
      i Did you forget to name an argument?

