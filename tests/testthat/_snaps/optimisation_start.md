# .validate_start_graphs with other start graphs

    Code
      .validate_start_graphs(list(list(hyp_weight = c(0.1, 0.2, NA, NA),
      trans_matrix = trans_m)), m = 5)
    Condition
      Error:
      ! `start_graph[[1]]$hyp_weight` has length 4; expected 5.

---

    Code
      .validate_start_graphs(list(list(hyp_weight = hyp_w, trans_matrix = hyp_w)), m = 5)
    Condition
      Error:
      ! `start_graph[[1]]$trans_matrix` must be a double matrix or `NULL`, not a double vector.

---

    Code
      .validate_start_graphs(list(list(hyp_weight = hyp_w, trans_matrix = trans_m_4)),
      m = 5)
    Condition
      Error:
      ! `start_graph[[1]]$trans_matrix` has dim (4 x 4); expected (5 x 5).

