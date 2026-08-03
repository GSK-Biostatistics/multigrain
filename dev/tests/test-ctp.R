test_that("apply_ctp() rejects dimension or shape mismatch", {
    p <- matrix(0.5, 3, 2)
    ws_bad <- matrix(1, 3, 3) # not 2*m
    expect_error(
        apply_ctp(p, 0.025, ws_bad),
        "2 \\* ncol\\(pvals\\) columns"
    )
})

test_that("apply_ctp_z() rejects shape mismatch", {
    z <- matrix(0, 3, 2)
    ws_bad <- matrix(1, 3, 3)
    expect_error(
        apply_ctp_z(z, 0.025, ws_bad),
        "2 \\* ncol\\(z\\) columns"
    )
})

test_that("m = 1 closed test is just a single comparison (p & z)", {
    alpha <- 0.025
    p <- matrix(c(0.01, 0.03, 0.5), ncol = 1)
    z <- qnorm(1 - p) # one-sided upper tail
    w <- 1 # all alpha on H1
    ws <- matrix(1, nrow = 1, ncol = 2)
    got_p <- apply_ctp(p, alpha, ws)
    got_z <- apply_ctp_z(z, alpha, ws)
    true_result <- p <= alpha * w
    expect_identical(drop(got_p), drop(true_result))
    expect_identical(drop(got_z), drop(true_result))
})

test_that("m = 2 analytic check vs hand calculation", {
    alpha <- 0.025
    w <- c(0.6, 0.4) # initial weights
    ws <- calc_local_weights(w, G = matrix(c(0, 1, 1, 0), nrow = 2))
    # three trials
    p <- matrix(
        c(
            0.01, 0.02,
            0.03, 0.001,  # only H2 tiny
            0.5,  0.5     # none
        ),
        ncol = 2,
        byrow = TRUE
    )
    got <- apply_ctp(p, alpha, ws)

    rej1 <- p[, 1] <= alpha
    rej2 <- p[, 2] <= alpha
    rej12 <- p[, 1] <= alpha * 0.6 | p[, 2] <= alpha * 0.4

    want <- cbind(
        rej1 & rej12,
        rej2 & rej12
    )

    expect_identical(drop(got), drop(want))

    # Z version equivalence
    z <- stats::qnorm(1 - p) # upper-tail
    got_z <- apply_ctp_z(z, alpha, ws)
    expect_identical(got_z, got)
})


test_that("Z === p equivalence over random draws (upper-tail one-sided)", {
    set.seed(2)
    m <- 5
    n <- 1000
    alpha <- 0.025
    z <- matrix(rnorm(n * m), n, m)
    p <- 1 - stats::pnorm(z) # same tail convention as apply_ctp_z
    ws <- calc_local_weights(
        random_weights(m) |> t(),
        random_transitions(m)
    )
    got_p <- apply_ctp(p, alpha, ws)
    got_z <- apply_ctp_z(z, alpha, ws)
    expect_identical(got_z, got_p)
})


test_that("Compare rejection matrix with gMCPLite::graphTest", {
    set.seed(12)
    m <- 6
    n <- 10000
    alpha <- 0.025
    z <- matrix(rnorm(n * m), n, m) + 3
    p <- 1 - stats::pnorm(z) # same tail convention as apply_ctp_z
    w <- random_weights(m) |> t()
    G <- random_transitions(m)
    ws <- calc_local_weights(w, G)
    got_p <- apply_ctp(p, alpha, ws)
    got_z <- apply_ctp_z(z, alpha, ws)
    expect_identical(got_z, got_p)

    #gMCPLite comparison
    # gMCPLite::checkArgs uses strict `rowSums(G) > 1` with no tolerance,
    # which can fail on macOS ARM where FP non-associativity causes row
    # sums to land 1 ULP above 1. Clamp before passing.
    G_gmcp <- G
    rs <- rowSums(G_gmcp)
    overshoot <- rs > 1
    if (any(overshoot)) {
        G_gmcp[overshoot, ] <- G_gmcp[overshoot, , drop = FALSE] / rs[overshoot]
    }
    gMCP_result <- gMCPLite::graphTest(
        p,
        weights = w,
        alpha = 0.025,
        G = G_gmcp
    )
    gMCP_result <- apply(gMCP_result, 2, as.logical)
    dimnames(gMCP_result) <- NULL
    expect_identical(
        apply(gMCP_result, 2, as.logical) |> drop(),
        got_p
    )
})


test_that("Row order in weighting_strategy does not affect decisions", {
    set.seed(2)
    m <- 7
    n <- 5000
    alpha <- 0.025
    p <- matrix(runif(n * m), n, m)
    z <- qnorm(1 - p)
    w <- random_weights(m) |> t()
    G <- random_transitions(m)
    ws <- calc_local_weights(w, G)
    shuf <- ws[sample(nrow(ws)), , drop = FALSE]

    got_p <- apply_ctp(p, alpha, ws)
    got_p2 <- apply_ctp(p, alpha, shuf)
    got_z <- apply_ctp_z(z, alpha, ws)
    got_z2 <- apply_ctp_z(z, alpha, shuf)

    expect_identical(got_p, got_p2)
    expect_identical(got_z, got_z2)
})
