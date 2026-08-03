# Trial success functions to be used for tests
disjunctive_3m_power <- trial_success(r1 || r2 || r3, verbose = "silent")
conjunctive_4m_power <- trial_success(r1 && r2 && r3 && r4, verbose = "silent")
avg_6m_power <- trial_success(r1 + r2 + r3 + r4 + r5 + r6, verbose = "silent")

# simulate p-values for power calculations for rest of tests
power_vec <- c(0.9, 0.8, 0.6, 0.85, 0.85, 0.5)
ncp_vec <- rep(stats::qnorm(1 - 0.025), 6) - stats::qnorm(1 - power_vec)
corr_test <- matrix(
    0.2,
    nrow = 6,
    ncol = 6,
    byrow = TRUE
)
diag(corr_test) <- 1

pvals <- withr::with_seed(5, {
    sims <- mvtnorm::rmvnorm(2^20, mean = ncp_vec, sigma = corr_test)
    pvals <- stats::pnorm(sims, lower.tail = FALSE)
})


# create_obj_func --------------------------------------------------------

test_that("create_obj_func: 3-m disjunctive power: improvedFallbackI", {
    pvals_3m <- pvals[, 1:3]

    gMCP_3m <- gMCPLite::improvedFallbackI()
    G_3m <- gMCPLite::getMatrix(gMCP_3m)
    w_3m <- gMCPLite::getWeights(gMCP_3m)
    no_constraint_3m <- graph_constraint_free(3)
    x0_3m <- create_start_params(no_constraint_3m, w0 = w_3m, G0 = G_3m)
    out <- gMCPLite::graphTest(
        pvalues = pvals_3m,
        weights = w_3m,
        alpha = 0.025,
        G = G_3m
    )
    power_disj <- extract_power(out)

    ## Test objective function gets the same power
    obj_func_disj <- create_obj_func(
        3,
        power_criterion = disjunctive_3m_power$func,
        hyp_constraint = no_constraint_3m$hyp_constraint,
        trans_constraint = no_constraint_3m$trans_constraint,
        pvals = pvals_3m
    )
    test_power <- obj_func_disj(
        x = x0_3m
    )

    expect_type(test_power, "double")
    expect_gt(test_power, 0)

    # Check if the power objective function is identical
    expect_snapshot({
        cat("3-m disjunctive power: improvedFallbackI")
        cat(paste("multigrain power:", round(test_power, 4)))
        cat(paste("gMCP power:", round(power_disj$disj_power, 4)))
    })

    expect_identical(test_power, power_disj$disj_power)
})


test_that("create_obj_func: 4-m conjunctive pow: improvedParallelGatekeeping", {
    pvals_4m <- pvals[, 1:4]

    gMCP_4m <- gMCPLite::improvedParallelGatekeeping() |>
        gMCPLite::substituteEps()
    G_4m <- gMCPLite::getMatrix(gMCP_4m)
    w_4m <- gMCPLite::getWeights(gMCP_4m)
    no_constraint_4m <- graph_constraint_free(4)
    x0_4m <- create_start_params(no_constraint_4m, w0 = w_4m, G0 = G_4m)

    out <- gMCPLite::graphTest(
        pvalues = pvals[, 1:4],
        weights = w_4m,
        alpha = 0.025,
        G = G_4m
    )
    power_conj <- extract_power(out)

    obj_func_conj <- create_obj_func(
        4,
        power_criterion = conjunctive_4m_power$func,
        hyp_constraint = no_constraint_4m$hyp_constraint,
        trans_constraint = no_constraint_4m$trans_constraint,
        pvals = pvals[, 1:4]
    )
    test_power <- obj_func_conj(
        x = x0_4m
    )

    expect_type(test_power, "double")
    expect_gt(test_power, 0)
    # Check if the power objective function within 0.001 of gMCP version
    expect_snapshot({
        cat("4-m conjunctive power: improvedParallelGatekeeping")
        cat(paste("multigrain power:", round(test_power, 4)))
        cat(paste("gMCP power:", round(power_conj$conj_power, 4)))
    })

    expect_identical(test_power, power_conj$conj_power)
})


test_that("create_obj_func: 6-m average power: BretzEtAl2011", {
    pvals_6m <- pvals[, 1:6]
    gMCP_6m <- gMCPLite::BretzEtAl2011()
    G_6m <- gMCPLite::getMatrix(gMCP_6m)
    w_6m <- gMCPLite::getWeights(gMCP_6m)
    no_constraint_6m <- graph_constraint_free(6)
    x0_6m <- create_start_params(no_constraint_6m, w0 = w_6m, G0 = G_6m)

    out <- gMCPLite::graphTest(
        pvalues = pvals_6m,
        weights = w_6m,
        alpha = 0.025,
        G = G_6m
    )
    power_avg <- extract_power(out)

    obj_func_avg <- create_obj_func(
        6,
        power_criterion = avg_6m_power$func,
        hyp_constraint = no_constraint_6m$hyp_constraint,
        trans_constraint = no_constraint_6m$trans_constraint,
        pvals = pvals_6m
    )
    test_power <- obj_func_avg(
        x = x0_6m
    )

    expect_type(test_power, "double")
    expect_gt(test_power, 0)

    # Check if the power objective function within 0.001 of gMCP version
    expect_snapshot({
        cat("6-m average power: BretzEtAl2011")
        cat(paste("multigrain power:", round(test_power, 4)))
        cat(paste("gMCP power:", round(power_avg$exp_rejections, 4)))
    })

    expect_identical(test_power, power_avg$exp_rejections)
})

# Epsilon edges ----------------------------------------------------------

# Test objects for handling epsilons edges

# For edge < 1e-5, calculate power as if edge = 0:
gMCP_edge_nearzero <- gMCPLite::improvedParallelGatekeeping() |>
    gMCPLite::substituteEps(5e-6)
G_edge_nearzero <- gMCPLite::getMatrix(gMCP_edge_nearzero)
G_edge_nearzero[G_edge_nearzero < 1e-5] <- 0
w_4m <- gMCPLite::getWeights(gMCP_edge_nearzero)

out <- gMCPLite::graphTest(
    pvalues = pvals[, 1:4],
    weights = w_4m,
    alpha = 0.025,
    G = G_edge_nearzero
)
power_nearzero <- extract_power(out)

# For edge > 1e-5, do not adjust power:
gMCP_edge_noadjust <- gMCPLite::improvedParallelGatekeeping() |>
    gMCPLite::substituteEps(5e-5)
G_edge_noadjust <- gMCPLite::getMatrix(gMCP_edge_noadjust)
G_edge_noadjust[G_edge_noadjust < 1e-5] <- 0
w_4m <- gMCPLite::getWeights(gMCP_edge_noadjust)

out <- gMCPLite::graphTest(
    pvalues = pvals[, 1:4],
    weights = w_4m,
    alpha = 0.025,
    G = G_edge_noadjust
)

power_noadjust <- extract_power(out)

## Test objects for handling near-zero weights and transitions (epsilon)
test_that("create_obj_func handles near-zero transition weights", {
    # Use 4-m conjunctive power:
    # * pvals
    # * conjunctive_4m_power()
    pvals_4m <- pvals[, 1:4]

    no_constraint_4m <- graph_constraint_free(4)

    # Check that power for when edge < 1e-5; same as edge = 0
    x0_nearzero <- c(0.5, 0.5, 0, 0, 0.5, 0, 0.5, 5e-06, 0, 0, 5e-06)

    # Create objective function for 4-m
    obj_func_conj <- create_obj_func(
        4,
        power_criterion = conjunctive_4m_power$func,
        hyp_constraint = no_constraint_4m$hyp_constraint,
        trans_constraint = no_constraint_4m$trans_constraint,
        pvals = pvals[, 1:4]
    )

    package_power <- obj_func_conj(
        x = x0_nearzero
    )

    # Check if the power objective function within 0.001 of gMCP version
    expect_snapshot({
        cat(paste("multigrain power:", round(package_power, 4)))
        cat(paste("gMCP power:", round(power_nearzero$conj_power, 4)))
    })

    expect_identical(
        package_power,
        power_nearzero$conj_power
    )

    # Check: power for when edge > 1e-5; power unchanged (no edge adjustment)
    # Check: power for when edge < 1e-5; same as edge = 0
    x0_nearzero <- c(0.5, 0.5, 0, 0, 0.5, 0, 0.5, 5e-06, 0, 0, 5e-06)

    # Create objective function for 4-m
    package_power <- obj_func_conj(
        x = create_start_params(
            no_constraint_4m,
            w0 = w_4m,
            G0 = G_edge_noadjust
        )
    )

    # Check if the power objective function within 0.001 of gMCP version
    expect_snapshot({
        cat(paste("multigrain power:", round(package_power, 4)))
        cat(paste("gMCP power:", round(power_noadjust$conj_power, 4)))
    })

    expect_identical(
        package_power,
        power_noadjust$conj_power
    )
})


## ------ NAN guards in create_obj_func ------ ##
test_that("create_obj_func returns penalty (not error) when x contains NaN", {
    # Failure from bug fixed byPR #503
    # COBYLA proposes a parameter vector with NaN → recover_full_weights
    # propagates it → any(hyp_weight < 0) returns NA → if(NA) crashes.
    pvals_4m <- pvals[, 1:4]

    no_constraint_4m <- graph_constraint_free(4)

    obj_func <- create_obj_func(
        4,
        power_criterion = conjunctive_4m_power$func,
        hyp_constraint = no_constraint_4m$hyp_constraint,
        trans_constraint = no_constraint_4m$trans_constraint,
        pvals = pvals_4m
    )

    # Build a valid x0, then inject NaN at a weight-parameter position
    x0 <- create_start_params(no_constraint_4m)
    x_nan_weight <- x0
    x_nan_weight[1] <- NaN

    result <- obj_func(
        x = x_nan_weight
    )

    expect_type(result, "double")
    expect_length(result, 1)
    expect_false(is.nan(result))
    expect_identical(result, -1e6)
})

test_that("create_obj_func returns penalty when trans_matrix params are NaN", {
    pvals_4m <- pvals[, 1:4]

    no_constraint_4m <- graph_constraint_free(4)

    obj_func <- create_obj_func(
        4,
        power_criterion = conjunctive_4m_power$func,
        hyp_constraint = no_constraint_4m$hyp_constraint,
        trans_constraint = no_constraint_4m$trans_constraint,
        pvals = pvals_4m
    )

    # Inject NaN into a transition-matrix parameter position
    x0 <- create_start_params(no_constraint_4m)
    x_nan_trans <- x0
    x_nan_trans[length(x0)] <- NaN

    result <- obj_func(
        x = x_nan_trans
    )

    expect_type(result, "double")
    expect_length(result, 1)
    expect_false(is.nan(result))
    expect_identical(result, -1e6)
})

test_that("create_obj_func returns penalty when all params are NaN", {
    pvals_4m <- pvals[, 1:4]

    no_constraint_4m <- graph_constraint_free(4)

    obj_func <- create_obj_func(
        4,
        power_criterion = conjunctive_4m_power$func,
        hyp_constraint = no_constraint_4m$hyp_constraint,
        trans_constraint = no_constraint_4m$trans_constraint,
        pvals = pvals_4m
    )

    x0 <- create_start_params(no_constraint_4m)
    x_all_nan <- rep(NaN, length(x0))

    result <- obj_func(
        x = x_all_nan
    )

    expect_identical(result, -1e6)
})

test_that("create_obj_func returns penalty when params contain NA", {
    # NA (not NaN) should also be caught by the same guard
    pvals_4m <- pvals[, 1:4]

    no_constraint_4m <- graph_constraint_free(4)

    obj_func <- create_obj_func(
        4,
        power_criterion = conjunctive_4m_power$func,
        hyp_constraint = no_constraint_4m$hyp_constraint,
        trans_constraint = no_constraint_4m$trans_constraint,
        pvals = pvals_4m
    )

    x0 <- create_start_params(no_constraint_4m)
    x_na <- x0
    x_na[2] <- NA_real_

    result <- obj_func(
        x = x_na
    )

    expect_identical(result, -1e6)
})

# create_obj_func parallel -----------------------------------------------

test_that("create_obj_func: 3-m disjunctive power matches serial", {
    pvals_3m <- pvals[, 1:3]
    no_constraint_3m <- graph_constraint_free(3)

    gMCP_3m <- gMCPLite::improvedFallbackI()
    G_3m <- gMCPLite::getMatrix(gMCP_3m)
    w_3m <- gMCPLite::getWeights(gMCP_3m)
    x0_3m <- create_start_params(no_constraint_3m, w0 = w_3m, G0 = G_3m)

    # Serial version
    obj_func_serial <- create_obj_func(
        3,
        power_criterion = disjunctive_3m_power$func,
        hyp_constraint = no_constraint_3m$hyp_constraint,
        trans_constraint = no_constraint_3m$trans_constraint,
        pvals = pvals_3m
    )
    power_serial <- obj_func_serial(
        x = x0_3m
    )

    # Parallel version (1 thread - should match exactly)
    obj_func_parallel <- create_obj_func(
        3,
        disjunctive_3m_power$func,
        hyp_constraint = no_constraint_3m$hyp_constraint,
        trans_constraint = no_constraint_3m$trans_constraint,
        pvals = pvals_3m,
        num_threads = 1L
    )
    power_parallel_1t <- obj_func_parallel(
        x = x0_3m
    )

    # Parallel version (2 threads)
    obj_func_parallel_2t <- create_obj_func(
        3,
        disjunctive_3m_power$func,
        hyp_constraint = no_constraint_3m$hyp_constraint,
        trans_constraint = no_constraint_3m$trans_constraint,
        pvals = pvals_3m,
        num_threads = 2L
    )
    power_parallel_2t <- obj_func_parallel_2t(
        x = x0_3m
    )

    # All should be identical (deterministic operation)
    expect_identical(power_serial, power_parallel_1t)
    expect_identical(power_serial, power_parallel_2t)
    expect_type(power_parallel_1t, "double")
    expect_gt(power_parallel_1t, 0)
})

test_that("create_obj_func: 4-m conjunctive power matches serial", {
    pvals_4m <- pvals[, 1:4]
    no_constraint_4m <- graph_constraint_free(4)

    gMCP_4m <- gMCPLite::improvedParallelGatekeeping() |>
        gMCPLite::substituteEps()
    G_4m <- gMCPLite::getMatrix(gMCP_4m)
    w_4m <- gMCPLite::getWeights(gMCP_4m)
    x0_4m <- create_start_params(no_constraint_4m, w0 = w_4m, G0 = G_4m)

    # Serial version
    obj_func_serial <- create_obj_func(
        4,
        power_criterion = conjunctive_4m_power$func,
        hyp_constraint = no_constraint_4m$hyp_constraint,
        trans_constraint = no_constraint_4m$trans_constraint,
        pvals = pvals_4m
    )
    power_serial <- obj_func_serial(
        x = x0_4m
    )

    # Parallel versions
    obj_func_parallel_1t <- create_obj_func(
        4,
        conjunctive_4m_power$func,
        hyp_constraint = no_constraint_4m$hyp_constraint,
        trans_constraint = no_constraint_4m$trans_constraint,
        pvals = pvals_4m,
        num_threads = 1L
    )
    power_parallel_1t <- obj_func_parallel_1t(
        x = x0_4m
    )

    obj_func_parallel_2t <- create_obj_func(
        4,
        conjunctive_4m_power$func,
        hyp_constraint = no_constraint_4m$hyp_constraint,
        trans_constraint = no_constraint_4m$trans_constraint,
        pvals = pvals_4m,
        num_threads = 2L
    )
    power_parallel_2t <- obj_func_parallel_2t(
        x = x0_4m
    )

    expect_identical(power_serial, power_parallel_1t)
    expect_identical(power_serial, power_parallel_2t)
    expect_type(power_parallel_2t, "double")
    expect_gt(power_parallel_2t, 0)
})

test_that("create_obj_func: 6-m average power matches serial", {
    pvals_6m <- pvals[, 1:6]
    no_constraint_6m <- graph_constraint_free(6)

    gMCP_6m <- gMCPLite::BretzEtAl2011()
    G_6m <- gMCPLite::getMatrix(gMCP_6m)
    w_6m <- gMCPLite::getWeights(gMCP_6m)
    x0_6m <- create_start_params(no_constraint_6m, w0 = w_6m, G0 = G_6m)

    # Serial version
    obj_func_serial <- create_obj_func(
        6,
        power_criterion = avg_6m_power$func,
        hyp_constraint = no_constraint_6m$hyp_constraint,
        trans_constraint = no_constraint_6m$trans_constraint,
        pvals = pvals_6m
    )
    power_serial <- obj_func_serial(
        x = x0_6m
    )

    # Parallel version
    obj_func_parallel <- create_obj_func(
        6,
        avg_6m_power$func,
        hyp_constraint = no_constraint_6m$hyp_constraint,
        trans_constraint = no_constraint_6m$trans_constraint,
        pvals = pvals_6m,
        num_threads = 2L
    )

    power_parallel <- obj_func_parallel(
        x = x0_6m
    )

    expect_identical(power_serial, power_parallel)
    expect_type(power_parallel, "double")
    expect_gt(power_parallel, 0)
})

test_that("create_obj_func handles near-zero weights & transitions", {
    pvals_4m <- pvals[, 1:4]
    no_constraint_4m <- graph_constraint_free(4)

    # Test with near-zero edges
    x0_nearzero <- c(0.5, 0.5, 0, 0, 0.5, 0, 0.5, 5e-06, 0, 0, 5e-06)

    obj_func_serial <- create_obj_func(
        4,
        power_criterion = conjunctive_4m_power$func,
        hyp_constraint = no_constraint_4m$hyp_constraint,
        trans_constraint = no_constraint_4m$trans_constraint,
        pvals = pvals_4m
    )
    power_serial <- obj_func_serial(
        x = x0_nearzero
    )

    obj_func_parallel <- create_obj_func(
        4,
        conjunctive_4m_power$func,
        hyp_constraint = no_constraint_4m$hyp_constraint,
        trans_constraint = no_constraint_4m$trans_constraint,
        pvals = pvals_4m,
        num_threads = 2L
    )
    power_parallel <- obj_func_parallel(
        x = x0_nearzero
    )

    expect_identical(power_serial, power_parallel)
    expect_type(power_parallel, "double")
    expect_gt(power_parallel, 0)
})


# split_theta ------------------------------------------------------------

test_that("split_theta basic all-free 3x3 graph works", {
    # nolint start: commented_code_linter
    cw <- c(NA_real_, NA_real_, NA_real_) # 3 free weights -> w_len = 2
    cG <- matrix(NA_real_, 3, 3) # each row 2 free -> 1 each -> g_len = 3
    diag(cG) <- 0

    info <- .expected_encoded_lengths(cw, cG)
    expect_equal(info$w_len, 2L)
    expect_equal(info$g_len, 3L)

    theta <- c(0.3, 0.2, 0.6, 0.1, 1) # 2 + 3
    # nolint end

    parts <- split_theta(theta, cw)
    expect_length(parts$w_pars, info$w_len)
    expect_length(parts$g_pars, info$g_len)

    # reconstruct and check simplex
    w_full <- recover_full_weights(parts$w_pars, cw)
    expect_equal(sum(w_full), 1, tolerance = 1e-12)

    G_full <- recover_full_trans_matrix(parts$g_pars, cG)
    expect_equal(rowSums(G_full), rep(1, 3), tolerance = 1e-12)
})

test_that("split_theta handles fixed + free weight mix", {
    set.seed(47)

    # nolint start: commented_code_linter
    cw <- c(0.2, NA, 0.2, NA) # free_w = 2 -> w_len = 1
    cG <- matrix(NA, 4, 4)
    cG[3, ] <- c(0, 1, 0, 0) # no free row
    diag(cG) <- 0 # all free rows (6 NA each) -> (3-1)*3 = 6 g params
    # nolint end

    info <- .expected_encoded_lengths(cw, cG)
    expect_identical(info$w_len, 1L)
    expect_identical(info$g_len, 6L)

    rand_w <- runif(info$w_len)
    rand_w <- rand_w / sum(rand_w) / 5 # equal to 0.2
    rand_G <- runif(info$g_len)
    rand_G <- c(
        rand_G[1:2] / sum(rand_G[1:2]) / 1.25,
        rand_G[3:4] / sum(rand_G[3:4]) / 2,
        rand_G[5:6] / sum(rand_G[5:6]) / 3
    )
    theta <- c(rand_w, rand_G)
    parts <- split_theta(theta, cw)

    expect_length(parts$w_pars, info$w_len)
    expect_length(parts$g_pars, info$g_len)

    # round trip
    w_full <- recover_full_weights(parts$w_pars, cw)
    expect_equal(sum(w_full), 1, tolerance = 1e-12)
    expect_identical(w_full[!is.na(cw)], cw[!is.na(cw)])

    G_full <- recover_full_trans_matrix(parts$g_pars, cG)
    expect_equal(rowSums(G_full), rep(1, 4), tolerance = 1e-12)
    expect_equal(G_full[1, 4], 0.2, tolerance = 1e-12)
    expect_equal(G_full[2, 4], 0.5, tolerance = 1e-12)
    expect_equal(G_full[3, ], c(0, 1, 0, 0), tolerance = 1e-12)
    expect_equal(G_full[4, 3], 2 / 3, tolerance = 1e-12)
})

test_that("split_theta free_w == 1 encodes zero weight params", {
    set.seed(99)

    # nolint start: commented_code_linter
    cw <- c(0.6, NA, 0.2, 0.2) # free_w = 1 -> w_len = 0
    cG <- matrix(NA, 4, 4)
    cG[3, ] <- c(0, 1, 0, 0) # no free row
    diag(cG) <- 0 # all free rows (6 NA each) -> (3-1)*3 = 6 g params
    # nolint end

    info <- .expected_encoded_lengths(cw, cG)
    expect_identical(info$w_len, 0L)
    expect_identical(info$g_len, 6L)

    rand_G <- runif(info$g_len)
    rand_G <- c(
        rand_G[1:2] / sum(rand_G[1:2]) / 1.25,
        rand_G[3:4] / sum(rand_G[3:4]) / 2,
        rand_G[5:6] / sum(rand_G[5:6]) / 3
    )
    theta <- c(rand_G)
    parts <- split_theta(theta, cw)

    expect_length(parts$w_pars, 0L)
    expect_length(parts$g_pars, info$g_len)

    w_full <- recover_full_weights(parts$w_pars, cw)
    expect_equal(w_full, c(0.6, 0, 0.2, 0.2), tolerance = 1e-12)
    expect_identical(w_full[!is.na(cw)], cw[!is.na(cw)])
})

test_that("split_theta encodes zero transition weights when G is fixed", {
    set.seed(2399)

    # nolint start: commented_code_linter
    cw <- c(0.2, NA, 0.2, NA, NA) # free_w = 3 -> w_len = 2
    # nolint end
    cG <- random_transitions(5)

    info <- .expected_encoded_lengths(cw, cG)
    expect_identical(info$w_len, 2L)
    expect_identical(info$g_len, 0L)

    rand_w <- runif(info$w_len)
    rand_w <- rand_w / sum(rand_w) / 5 # equal to 0.2
    theta <- c(rand_w)
    parts <- split_theta(theta, cw)

    expect_length(parts$w_pars, 2L)
    expect_length(parts$g_pars, 0L)

    # round trip
    w_full <- recover_full_weights(parts$w_pars, cw)
    expect_equal(sum(w_full), 1, tolerance = 1e-12)
    expect_identical(w_full[!is.na(cw)], cw[!is.na(cw)])

    G_full <- recover_full_trans_matrix(parts$g_pars, cG)
    expect_equal(G_full, cG, tolerance = 1e-12)
})

test_that("split_theta respects row-wise g constraints", {
    # 4x4 example with row heterogeneity

    cG <- matrix(
        c(
            0, NA, NA, NA,
            0.5, 0, 0.5, 0,
            NA, NA,  0, 0.2,
            NA, 0, NA, 0
        ),
        nrow = 4,
        byrow = TRUE
    )

    cw <- c(NA, NA, NA, NA) # 4 free weights -> w_len = 3

    info <- .expected_encoded_lengths(cw, cG)
    expect_identical(info$w_len, 3L) # (4-1)
    expect_identical(info$g_len, 4L)

    theta <- c(c(0.1, 0.2, 0.5), c(0.1, 0.2, 0.5, 0.7))
    parts <- split_theta(theta, cw)

    expect_length(parts$w_pars, info$w_len)
    expect_length(parts$g_pars, info$g_len)

    # round trip
    w_full <- recover_full_weights(parts$w_pars, cw)
    expect_equal(sum(w_full), 1, tolerance = 1e-12)

    G_full <- recover_full_trans_matrix(parts$g_pars, cG)
    expect_equal(rowSums(G_full), rep(1, 4), tolerance = 1e-12)
    # row2 fixed row -> unchanged
    expect_equal(G_full[2, ], c(0.5, 0, 0.5, 0), tolerance = 1e-12)
})


# recover_full_weights ---------------------------------------------------

test_that("recover_full_weights() basic functionality", {
    # 3H no constraints
    no_constraints <- c(NA, NA, NA)
    x <- c(0.1, 0.4)

    expect_vector(recover_full_weights(x, no_constraints), 3)

    expect_identical(
        sum(recover_full_weights(x, no_constraints)),
        1
    )
    expect_true(all(recover_full_weights(x, no_constraints) >= 0))
    expect_true(all(recover_full_weights(x, no_constraints) <= 1))

    # 7H with 6 constraints
    constraint_7_6 <- c(0, 0, 0, 0, 0, 0, NA)
    x <- NULL
    expect_identical(
        recover_full_weights(x, constraint_7_6),
        c(0, 0, 0, 0, 0, 0, 1)
    )
})


# recover_full_trans_matrix ----------------------------------------------

test_that("recover_full_trans_matrix bug from issue #137 fix", {
    G_constr <- matrix(
        c(
            0,  1,  0,
            NA, 0, NA,
            NA, NA, 0
        ),
        nrow = 3,
        byrow = TRUE
    )

    x <- c(0.5, 0.5)
    expected <- matrix(
        c(
            0.0, 1.0, 0.0,
            0.5, 0.0, 0.5,
            0.5, 0.5, 0.0
        ),
        nrow = 3,
        byrow = TRUE
    )

    result <- recover_full_trans_matrix(x, G_constr)
    expect_identical(result, expected)
    expect_identical(rowSums(result), rep(1, 3))
})

test_that("recover_full_trans_matrix fully specified matrix passthrough", {
    G_constr <- matrix(
        c(
            0, 1, 0,
            0.5, 0, 0.5,
            0.2, 0.8, 0
        ),
        nrow = 3,
        byrow = TRUE
    )
    x <- numeric(0)
    result <- recover_full_trans_matrix(x, G_constr)
    expect_identical(result, G_constr)
})

test_that("recover_full_trans_matrix mixed rows", {
    G_constr <- matrix(
        c(
            0,  1,  0,  0,
            NA, 0, NA, NA,
            NA, NA, 0,  0,
            0,  0, NA, 0
        ),
        nrow = 4,
        byrow = TRUE
    )
    x <- c(0.2, 0.3, 0.4)
    result <- recover_full_trans_matrix(x, G_constr)
    expected <- matrix(
        c(
            0.0, 1.0, 0.0, 0.0,
            0.2, 0.0, 0.3, 0.5,
            0.4, 0.6, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0
        ),
        nrow = 4,
        byrow = TRUE
    )
    expect_identical(result, expected)
    expect_identical(rowSums(result), rep(1, 4))
})

test_that("recover_full_trans_matrix errors on wrong opt x length", {
    G_constr <- matrix(
        c(
            0, 1, 0,
            NA, 0, NA,
            NA, NA, 0
        ),
        nrow = 3,
        byrow = TRUE
    )
    expect_error(
        recover_full_trans_matrix(
            numeric(0),
            G_constr
        ),
        "does not match"
    )
})
