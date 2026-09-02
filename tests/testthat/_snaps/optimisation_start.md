# .build_start_matrix with default start graph

    Code
      .build_start_matrix(gc, start_graph = NULL)
    Output
           [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
      [1,]  0.2  0.2  0.2  0.2 0.25 0.25 0.25 0.25 0.25  0.25  0.25  0.25  0.25  0.25
      [2,]  1.0  0.0  0.0  0.0 1.00 0.00 0.00 0.00 1.00  0.00  0.00  0.00  1.00  0.00
           [,15] [,16] [,17] [,18] [,19]
      [1,]  0.25  0.25  0.25  0.25  0.25
      [2,]  0.00  0.00  0.25  0.25  0.25

---

    Code
      .build_start_matrix(gc, start_graph = list(list(hyp_weight = NULL,
        trans_matrix = NULL)))
    Output
           [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
      [1,]  0.2  0.2  0.2  0.2 0.25 0.25 0.25 0.25 0.25  0.25  0.25  0.25  0.25  0.25
      [2,]  1.0  0.0  0.0  0.0 1.00 0.00 0.00 0.00 1.00  0.00  0.00  0.00  1.00  0.00
           [,15] [,16] [,17] [,18] [,19]
      [1,]  0.25  0.25  0.25  0.25  0.25
      [2,]  0.00  0.00  0.25  0.25  0.25

# .validate_start_graphs with other start graphs

    Code
      .validate_start_graphs(list(list(hyp_weight = c(0.1, 0.2, NA, NA),
      trans_matrix = trans_m)), m = 5)
    Condition
      Error:
      ! graph_optimise: start_graph[[1]] has weight length 4; expected 5.

---

    Code
      .validate_start_graphs(list(list(hyp_weight = hyp_w, trans_matrix = hyp_w)), m = 5)
    Condition
      Error:
      ! graph_optimise: start_graph[[1]]$trans_matrix must be a matrix.

---

    Code
      .validate_start_graphs(list(list(hyp_weight = hyp_w, trans_matrix = trans_m_4)),
      m = 5)
    Condition
      Error:
      ! graph_optimise: start_graph[[1]]$trans_matrix has dim (4 x 4); expected (5 x 5).

