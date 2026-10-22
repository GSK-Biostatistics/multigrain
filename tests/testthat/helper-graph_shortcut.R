# Random (w, G) generator; row-sums of G are exactly 1 for m >= 2.
make_test_graph <- function(m, seed = NULL) {
    if (!is.null(seed)) {
        set.seed(seed)
    }
    gr <- graph_random(m)
    list(w = gr$hyp_weight, G = gr$trans_matrix)
}

# Random (w, G) on a grid of 1 / n: w and each row of G are multinomial
# counts divided by n, so critical values are round numbers that rounded
# p-values can hit exactly.
make_grid_graph <- function(m, seed = NULL, n = 20) {
    if (!is.null(seed)) {
        set.seed(seed)
    }
    w <- as.vector(stats::rmultinom(1, n, rep(1, m))) / n
    G <- matrix(0, m, m)
    for (i in seq_len(m)) {
        G[i, -i] <- as.vector(stats::rmultinom(1, n, rep(1, m - 1))) / n
    }
    list(w = w, G = G)
}

# Strip dimnames and coerce to plain integer matrix for robust comparison
# between LogicalMatrix (graph_shortcut) and numeric-0/1 matrix (graphTest).
as_int_mat <- function(x) {
    m <- matrix(as.integer(x), nrow = nrow(x), ncol = ncol(x))
    dimnames(m) <- NULL
    m
}

# Kernels behind graph_shortcut() / graph_shortcut_parallel(), selected with
# the internal option multigrain.shortcut_kernel
shortcut_kernels <- c("sequential", "table", "auto")

# Every kernel through every backend: serial, num_threads = 1 and 2.
shortcut_all_kernels <- function(pvals, alpha, w, G) {
    out <- list()
    for (k in shortcut_kernels) {
        withr::with_options(list(multigrain.shortcut_kernel = k), {
            out[[paste0(k, "_serial")]] <- graph_shortcut(pvals, alpha, w, G)
            out[[paste0(k, "_t1")]] <- graph_shortcut_parallel(
                pvals,
                alpha,
                w,
                G,
                num_threads = 1
            )
            out[[paste0(k, "_t2")]] <- graph_shortcut_parallel(
                pvals,
                alpha,
                w,
                G,
                num_threads = 2
            )
        })
    }
    out
}

# Every element of a named list of kernel results is identical to the first
expect_identical_kernels <- function(outs) {
    for (nm in names(outs)[-1]) {
        testthat::expect_identical(outs[[nm]], outs[[1]], label = nm)
    }
}
