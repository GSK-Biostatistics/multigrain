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

# calc_ncp

    Code
      calc_ncp(power = c(0.8, 0.9))
    Output
      [1] 2.801585 3.241516

---

    Code
      calc_ncp(power = c(0.8, 0.9), alpha = 0.01)
    Output
      [1] 3.167969 3.607899

