test_that("graph_violation_score_cpp matches R for simple valid graph", {
    w <- c(0.25, 0.25, 0.25, 0.25)
    G <- matrix(0.25, 4, 4)
    diag(G) <- 0 # diag handled upstream; still ok
    # row sums not 1 now because diag=0: row sum = 0.75
    expect_identical(
        graph_violation_score_cpp(w, G),
        graph_violation_score_r(w, G),
        tolerance = sqrt(.Machine$double.eps)
    )
})

test_that("valid simplex gives zero violation when rows & weights sum 1", {
    w <- c(0.4, 0.3, 0.2, 0.1)
    G <- matrix(0, 4, 4)
    # distribute row mass off-diagonal to sum to 1
    # simple example: evenly across others
    for (i in 1:4) {
        G[i, -i] <- 1 / 3
    }
    expect_identical(
        sum(abs(rowSums(G) - 1)),
        0,
        tolerance = sqrt(.Machine$double.eps)
    )
    expect_identical(abs(sum(w) - 1), 0, tolerance = sqrt(.Machine$double.eps))

    v_cpp <- graph_violation_score_cpp(w, G)
    v_r <- graph_violation_score_r(w, G)
    expect_identical(v_cpp, v_r, tolerance = sqrt(.Machine$double.eps))
    expect_identical(v_cpp, 0, tolerance = sqrt(.Machine$double.eps))
})

test_that("out-of-bounds weights accumulate correctly", {
    w <- c(1.2, -0.1, 0.5, -0.6) # sums weirdly; >1 and <0
    G <- matrix(0, 4, 4)
    for (i in 1:4) {
        G[i, -i] <- 1 / 3
    }

    v_r <- graph_violation_score_r(w, G)
    v_cpp <- graph_violation_score_cpp(w, G)
    expect_identical(v_cpp, v_r, tolerance = sqrt(.Machine$double.eps))

    # manual: amount below 0 = 0.1 + 0.6; amount above 1 = 0.2
    manual_bounds <- (0.1 + 0.6) + 0.2
    # w sum = 1.0? let's check: 1.2 -0.1 +0.5 -0.6 = 1.0; simplex term 0
    expect_identical(manual_bounds, 0.9, tolerance = sqrt(.Machine$double.eps))
    # row simplex term 0 (rows sum 1)
    expect_identical(v_r, manual_bounds, tolerance = sqrt(.Machine$double.eps))
})

test_that("G bounds & row simplex contribute to violation", {
    w <- c(0.4, 0.3, 0.2, 0.1)
    G <- matrix(
        c(
            0,  1.2,  -0.1, 0,
            0,   0,    0.5, 1.0,
            0,   0,    0,   0,
            0,   0,    0,   0
        ),
        4,
        4,
        byrow = TRUE
    )

    v_r <- graph_violation_score_r(w, G)
    v_cpp <- graph_violation_score_cpp(w, G)
    expect_identical(v_cpp, v_r, tolerance = sqrt(.Machine$double.eps))
})

test_that("NAs ignored consistently", {
    w <- c(0.4, NA, 0.2, 0.4) # w sum (ignoring NA) = 1.0
    G <- matrix(NA_real_, 4, 4)
    G[1, 2] <- 0.5
    G[1, 3] <- 0.5
    G[2, 1] <- 1
    G[3, 4] <- 1
    G[4, 3] <- 1

    v_r <- graph_violation_score_r(w, G)
    v_cpp <- graph_violation_score_cpp(w, G)
    expect_identical(v_cpp, v_r, tolerance = sqrt(.Machine$double.eps))
})

### INPUTS ###

test_that("graph_shortcut errors on G and w dimension mismatch", {
    pvals <- matrix(stats::runif(30), nrow = 10, ncol = 3)
    w <- c(0.3, 0.4, 0.3)
    G <- matrix(
        c(
            0, 0.5, 0.5,
            0.5, 0,   0.5,
            0.5, 0.5, 0
        ),
        nrow = 3,
        byrow = TRUE
    )

    expect_error(graph_shortcut(pvals, 0.025, w, matrix(0, 4, 4)), "m x m")
    expect_error(graph_shortcut(pvals, 0.025, c(0.5, 0.5), G), "length")
})

test_that("graph_shortcut_parallel validates num_threads and grain_size", {
    pvals <- matrix(stats::runif(30), nrow = 10, ncol = 3)
    gr <- make_test_graph(3, seed = 1)

    expect_error(
        graph_shortcut_parallel(pvals, 0.025, gr$w, gr$G, num_threads = 0),
        "num_threads must be >= 1"
    )
    expect_error(
        graph_shortcut_parallel(
            pvals,
            0.025,
            gr$w,
            gr$G,
            num_threads = 1,
            grain_size = 0
        ),
        "grain_size must be >= 1"
    )
})

### OUTPUTS ###

test_that("graph_shortcut returns a LogicalMatrix with correct dims", {
    gr <- make_test_graph(5, seed = 2)
    out <- graph_shortcut(pvals_fixture[1:10, 1:5], 0.025, gr$w, gr$G)
    expect_true(is.matrix(out))
    expect_type(out, "logical")
    expect_identical(dim(out), c(10L, 5L))
})

test_that("graph_shortcut_parallel returns a LogicalMatrix with correct dims", {
    gr <- make_test_graph(5, seed = 2)
    out <- graph_shortcut_parallel(
        pvals_fixture[1:10, 1:5],
        0.025,
        gr$w,
        gr$G,
        num_threads = 2
    )
    expect_true(is.matrix(out))
    expect_type(out, "logical")
    expect_identical(dim(out), c(10L, 5L))
})

### CASE CHECK ###

test_that("all large p-values -> no rejections", {
    gr <- make_test_graph(4, seed = 11)
    pvals <- matrix(0.5, nrow = 20, ncol = 4)
    expect_false(any(graph_shortcut(pvals, 0.025, gr$w, gr$G)))
})

test_that("all tiny p-values -> everything rejected", {
    gr <- make_test_graph(4, seed = 12)
    pvals <- matrix(1e-10, nrow = 20, ncol = 4)
    expect_true(all(graph_shortcut(pvals, 0.025, gr$w, gr$G)))
})

test_that("Bonferroni (no recycling) == p < alpha/m rule element-wise", {
    m <- 4L
    w <- rep(1 / m, m)
    G <- matrix(0, m, m) # no edges -> no alpha recycling

    set.seed(42)
    pvals <- matrix(stats::runif(100 * m), nrow = 100, ncol = m)
    alpha <- 0.05

    out <- graph_shortcut(pvals, alpha, w, G)
    expect_identical(as_int_mat(out), as_int_mat(pvals < alpha / m))
})

test_that("fixed sequence: rejects contiguous prefix up to first failure", {
    w <- c(1, 0, 0)
    G <- matrix(
        c(
            0, 1, 0,
            0, 0, 1,
            0, 0, 0
        ),
        nrow = 3,
        byrow = TRUE
    )

    alpha <- 0.025

    pvals <- rbind(
        c(0.001, 0.02, 0.10), # H1, H2 reject; H3 does not
        c(0.10, 0.001, 0.001), # H1 fails -> nothing downstream
        c(0.001, 0.001, 0.001) # all reject
    )
    out <- graph_shortcut(pvals, alpha, w, G)
    expected <- matrix(
        c(
            TRUE,  TRUE,  FALSE,
            FALSE, FALSE, FALSE,
            TRUE,  TRUE,  TRUE
        ),
        nrow = 3,
        byrow = TRUE
    )
    expect_identical(as_int_mat(out), as_int_mat(expected))
})

test_that("m = 1 trivially reduces to p < alpha", {
    out <- graph_shortcut(
        matrix(c(0.001, 0.02, 0.03, 0.5), ncol = 1),
        0.025,
        1,
        matrix(0, 1, 1)
    )
    expect_identical(as.logical(out), c(TRUE, TRUE, FALSE, FALSE))
})

### VALIDATION WITH graphTest ###

test_that("graph_shortcut matches graphTest on improvedFallbackI (m = 3)", {
    skip_if_not_installed("gMCPLite")

    g <- gMCPLite::improvedFallbackI()
    w <- gMCPLite::getWeights(g)
    G <- gMCPLite::getMatrix(g)
    pvals <- pvals_fixture[, 1:3]
    alpha <- 0.025

    out_new <- graph_shortcut(pvals, alpha, w, G)
    out_par <- graph_shortcut_parallel(pvals, alpha, w, G, num_threads = 2)
    out_ref <- gMCPLite::graphTest(
        pvalues = pvals,
        weights = w,
        alpha = alpha,
        G = G
    )

    expect_identical(as_int_mat(out_new), as_int_mat(out_ref))
    expect_identical(out_par, out_new)
})

# nolint start: line_length_linter
test_that("graph_shortcut matches graphTest on improvedParallelGatekeeping (m = 4, eps)", {
    # nolint end
    skip_if_not_installed("gMCPLite")

    # substituteEps replaces symbolic 'eps' edges with a small numeric value
    # (default 1e-3). Tests the shortcut's graph update with near-zero edges.
    g <- gMCPLite::substituteEps(gMCPLite::improvedParallelGatekeeping())
    w <- gMCPLite::getWeights(g)
    G <- gMCPLite::getMatrix(g)
    pvals <- pvals_fixture[, 1:4]
    alpha <- 0.025

    out_new <- graph_shortcut(pvals, alpha, w, G)
    out_par <- graph_shortcut_parallel(pvals, alpha, w, G, num_threads = 2)
    out_ref <- gMCPLite::graphTest(
        pvalues = pvals,
        weights = w,
        alpha = alpha,
        G = G
    )

    expect_identical(as_int_mat(out_new), as_int_mat(out_ref))
    expect_identical(out_par, out_new)
})

test_that("graph_shortcut matches graphTest on BretzEtAl2011 (m = 6)", {
    skip_if_not_installed("gMCPLite")

    g <- gMCPLite::BretzEtAl2011()
    w <- gMCPLite::getWeights(g)
    G <- gMCPLite::getMatrix(g)
    pvals <- pvals_fixture[, 1:6]
    alpha <- 0.025

    out_new <- graph_shortcut(pvals, alpha, w, G)
    out_par <- graph_shortcut_parallel(pvals, alpha, w, G, num_threads = 2)
    out_ref <- gMCPLite::graphTest(
        pvalues = pvals,
        weights = w,
        alpha = alpha,
        G = G
    )

    expect_identical(as_int_mat(out_new), as_int_mat(out_ref))
    expect_identical(out_par, out_new)
})

### Equivalence to graphTest on random graphs ###

test_that("graph_shortcut matches graphTest on random graphs (m = 2..6)", {
    skip_if_not_installed("gMCPLite")

    for (m in 2:6) {
        gr <- make_test_graph(m, seed = 101L + m)
        pvals <- pvals_fixture[1:2000, seq_len(m)]
        alpha <- 0.025

        out_new <- graph_shortcut(pvals, alpha, gr$w, gr$G)
        out_ref <- gMCPLite::graphTest(
            pvalues = pvals,
            weights = gr$w,
            alpha = alpha,
            G = gr$G
        )
        expect_equal(
            as_int_mat(out_new),
            as_int_mat(out_ref),
            info = paste("m =", m)
        )
    }
})

### SERIAL equivalent to PARALLEL ###

test_that("graph_shortcut_parallel(num_threads = 1) == graph_shortcut", {
    for (m in c(3L, 5L, 6L)) {
        gr <- make_test_graph(m, seed = m * 7L)
        pvals <- pvals_fixture[1:5000, seq_len(m)]

        out_s <- graph_shortcut(pvals, 0.025, gr$w, gr$G)
        out_p <- graph_shortcut_parallel(
            pvals,
            0.025,
            gr$w,
            gr$G,
            num_threads = 1
        )
        expect_identical(out_s, out_p, info = paste("m =", m))
    }
})

test_that("graph_shortcut_parallel output is invariant to thread count", {
    m <- 6L
    gr <- make_test_graph(m, seed = 2024)
    pvals <- pvals_fixture[, seq_len(m)]
    alpha <- 0.025

    out_ref <- graph_shortcut(pvals, alpha, gr$w, gr$G)
    for (nt in c(1L, 2L, 4L, 8L)) {
        out_nt <- graph_shortcut_parallel(
            pvals,
            alpha,
            gr$w,
            gr$G,
            num_threads = nt
        )
        expect_identical(out_nt, out_ref, info = paste("num_threads =", nt))
    }
})

test_that("graph_shortcut_parallel output is invariant to grain size", {
    m <- 5L
    gr <- make_test_graph(m, seed = 99)
    pvals <- pvals_fixture[1:5000, seq_len(m)]

    out_auto <- graph_shortcut_parallel(
        pvals,
        0.025,
        gr$w,
        gr$G,
        num_threads = 2,
        grain_size = -1
    )
    for (gs in c(1L, 50L, 500L, 5000L)) {
        out_gs <- graph_shortcut_parallel(
            pvals,
            0.025,
            gr$w,
            gr$G,
            num_threads = 2,
            grain_size = gs
        )
        expect_identical(out_gs, out_auto, info = paste("grain_size =", gs))
    }
})

test_that("graph_shortcut_parallel is deterministic across repeated calls", {
    m <- 5L
    gr <- make_test_graph(m, seed = 1)
    pvals <- pvals_fixture[1:5000, seq_len(m)]

    runs <- replicate(
        5,
        graph_shortcut_parallel(pvals, 0.025, gr$w, gr$G, num_threads = 4),
        simplify = FALSE
    )
    for (i in 2:5) {
        expect_identical(runs[[i]], runs[[1]])
    }
})


### Equivalence to apply_ctp / apply_ctp_parallel ###
test_that("graph_shortcut matches apply_ctp on continuous p-values", {
    skip_if_not(
        exists("apply_ctp"),
        "apply_ctp not available (removed post-swap)"
    )

    for (m in c(3L, 5L, 6L)) {
        gr <- make_test_graph(m, seed = m)
        pvals <- pvals_fixture[1:2000, seq_len(m)]
        ws <- calc_local_weights(gr$w, gr$G)

        expect_identical(
            graph_shortcut(pvals, 0.025, gr$w, gr$G),
            apply_ctp(pvals, 0.025, ws),
            info = paste("m =", m)
        )
    }
})

test_that("graph_shortcut_parallel == apply_ctp_parallel on cont p-values", {
    skip_if_not(
        exists("apply_ctp_parallel"),
        "apply_ctp_parallel not available (removed post-swap)"
    )

    for (m in c(3L, 5L, 6L)) {
        gr <- make_test_graph(m, seed = m + 101L)
        pvals <- pvals_fixture[1:5000, seq_len(m)]
        ws <- calc_local_weights(gr$w, gr$G)

        expect_identical(
            graph_shortcut_parallel(pvals, 0.025, gr$w, gr$G, num_threads = 2),
            apply_ctp_parallel(pvals, 0.025, ws, num_threads = 2),
            info = paste("m =", m)
        )
    }
})

### TIE BREAK edge cases ###

test_that("graph_shortcut uses strict `<` at the alpha*w boundary", {
    w <- c(0.5, 0.5)
    G <- matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE)
    alpha <- 0.02

    # Row 1: p1 exactly on boundary -> NOT rejected under `<`
    # Row 2: p1 just below boundary -> rejected; alpha fully recycled to H2,
    #        p2 = 0.5 > new threshold so H2 not rejected
    # Row 3: both p's tiny -> both rejected
    pvals <- rbind(
        c(alpha * w[1], 0.5),
        c(alpha * w[1] - 1e-15, 0.5),
        c(1e-6, 1e-6)
    )
    expected <- matrix(
        c(
            FALSE, FALSE,
            TRUE,  FALSE,
            TRUE,  TRUE
        ),
        nrow = 3,
        byrow = TRUE
    )
    expect_identical(
        as_int_mat(graph_shortcut(pvals, alpha, w, G)),
        as_int_mat(expected)
    )
})

test_that("graph_shortcut and graphTest agree on exact-boundary p-values", {
    skip_if_not_installed("gMCPLite")

    w <- c(0.5, 0.5)
    G <- matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE)
    alpha <- 0.02
    pvals <- matrix(c(alpha * w[1], 0.5), nrow = 1) # exact-boundary row

    expect_identical(
        as_int_mat(graph_shortcut(pvals, alpha, w, G)),
        as_int_mat(
            gMCPLite::graphTest(
                pvalues = pvals,
                weights = w,
                alpha = alpha,
                G = G
            )
        )
    )
})


### Other edge cases ###

test_that("single-trial input (N = 1) is handled correctly", {
    gr <- make_test_graph(3L, seed = 1)
    out <- graph_shortcut(
        matrix(c(0.001, 0.5, 0.5), nrow = 1),
        0.025,
        gr$w,
        gr$G
    )
    expect_identical(dim(out), c(1L, 3L))
    expect_type(out, "logical")
})

test_that("zero-weight hypotheses can't reject without an active predecessor", {
    # H3, H4 carry no initial weight; alpha can only reach them by recycling
    # from H1, H2 respectively.
    w <- c(0.5, 0.5, 0, 0)
    G <- matrix(
        c(
            0, 0, 1, 0,
            0, 0, 0, 1,
            0, 0, 0, 0,
            0, 0, 0, 0
        ),
        nrow = 4,
        byrow = TRUE
    )

    out <- graph_shortcut(pvals_fixture[1:500, 1:4], 0.025, w, G)

    expect_false(any(!out[, 1] & out[, 3])) # H3 rejection implies H1 rejection
    expect_false(any(!out[, 2] & out[, 4])) # H4 rejection implies H2 rejection
})

test_that("weights summing below 1 (unused alpha) still match graphTest", {
    skip_if_not_installed("gMCPLite")

    w <- c(0.3, 0.3, 0.3) # sums to 0.9; some alpha never used
    G <- matrix(
        c(
            0,   0.5, 0.5,
            0.5, 0,   0.5,
            0.5, 0.5, 0
        ),
        nrow = 3,
        byrow = TRUE
    )
    pvals <- pvals_fixture[1:500, 1:3]

    expect_identical(
        as_int_mat(graph_shortcut(pvals, 0.025, w, G)),
        as_int_mat(
            gMCPLite::graphTest(
                pvalues = pvals,
                weights = w,
                alpha = 0.025,
                G = G
            )
        )
    )
})

test_that("read-only inputs are not mutated", {
    # Implementation uses raw REAL(G) / REAL(w) pointers; confirm no writes.
    gr <- make_test_graph(5L, seed = 3)
    pvals <- pvals_fixture[1:500, 1:5]

    w_before <- gr$w
    G_before <- gr$G
    p_before <- pvals

    invisible(
        graph_shortcut(
            pvals,
            0.025,
            gr$w,
            gr$G
        )
    )
    invisible(
        graph_shortcut_parallel(
            pvals,
            0.025,
            gr$w,
            gr$G,
            num_threads = 2
        )
    )

    expect_identical(gr$w, w_before)
    expect_identical(gr$G, G_before)
    expect_identical(pvals, p_before)
})

test_that("Validate calc_local_weights; gMCPLite using Bonferroni Holm graph", {
    bh_G5 <- gMCPLite::BonferroniHolm(5)
    g <- gMCPLite::getMatrix(bh_G5)
    w <- gMCPLite::getWeights(bh_G5)
    weighting_strategy <- gMCPLite::generateWeights(g, w)
    #fmt: skip
    weighting_strategy <- weighting_strategy[rev(seq_len(nrow(weighting_strategy))), ] # nolint
    dimnames(weighting_strategy) <- NULL

    expect_equal(
        calc_local_weights(w, g),
        weighting_strategy,
        tolerance = 100 * .Machine$double.eps
    )
})


test_that("Validate calc_local_weights; gMCPLite using fixed sequence", {
    fs4 <- gMCPLite::fixedSequence(4)
    g <- gMCPLite::getMatrix(fs4)
    w <- gMCPLite::getWeights(fs4)

    weighting_strategy <- gMCPLite::generateWeights(g, w)
    # fmt: skip
    weighting_strategy <- weighting_strategy[rev(seq_len(nrow(weighting_strategy))), ] # nolint
    dimnames(weighting_strategy) <- NULL

    expect_equal(
        calc_local_weights(w, g),
        weighting_strategy,
        tolerance = 100 * .Machine$double.eps
    )
})

# nolint start: line_length_linter
test_that("Validate calc_local_weights; gMCPLite + user-defined matrix and weights", {
    # nolint end
    # 3 hypothesis
    g <- matrix(
        c(
            0, 0.4, 0.6,
            0.3, 0, 0.7,
            1, 0, 0
        ),
        byrow = TRUE,
        nrow = 3
    )

    w <- c(0.3, 0.5, 0.2)
    weighting_strategy <- gMCPLite::generateWeights(g, w)
    # fmt: skip
    weighting_strategy <- weighting_strategy[rev(seq_len(nrow(weighting_strategy))), ] # nolint
    dimnames(weighting_strategy) <- NULL

    expect_equal(
        calc_local_weights(w, g),
        weighting_strategy,
        tolerance = 100 * .Machine$double.eps
    )

    # random matrix
    set.seed(112)
    w <- random_weights(7)
    g <- random_transitions(7)

    weighting_strategy <- gMCPLite::generateWeights(g, w)
    # fmt: skip
    weighting_strategy <- weighting_strategy[rev(seq_len(nrow(weighting_strategy))), ] # nolint
    dimnames(weighting_strategy) <- NULL

    expect_equal(
        calc_local_weights(w, g),
        weighting_strategy,
        tolerance = 100 * .Machine$double.eps
    )
})
