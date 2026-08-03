# Random (w, G) generator; row-sums of G are exactly 1 for m >= 2.
make_test_graph <- function(m, seed = NULL) {
    if (!is.null(seed)) {
        set.seed(seed)
    }
    gr <- graph_random(m)
    list(w = gr$hyp_weight, G = gr$trans_matrix)
}

# Strip dimnames and coerce to plain integer matrix for robust comparison
# between LogicalMatrix (graph_shortcut) and numeric-0/1 matrix (graphTest).
as_int_mat <- function(x) {
    m <- matrix(as.integer(x), nrow = nrow(x), ncol = ncol(x))
    dimnames(m) <- NULL
    m
}
