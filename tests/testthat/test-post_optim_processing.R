## param_to_solution tests
test_that("Correct param_to_solution() behaviour for 2 hypotheses", {
    gc <- graph_constraint(c(NA, NA), matrix(c(0, 1, 1, 0), nrow = 2))
    optim_param <- 0.25
    sol <- param_to_solution(optim_param, gc)
    expect_identical(sol$hyp_weight, c(0.25, 0.75))
    expect_identical(sol$trans_matrix, unname(gc$trans_constraint))
})

test_that("param_to_solution(process=TRUE) zeros weights below 1e-4 (m=2)", {
    # With hyp_constraint = c(NA, NA) and trans fully fixed,
    # optim_param is a single weight value for w[1]; w[2] = 1 - w[1].
    gc <- graph_constraint(c(NA, NA), matrix(c(0, 1, 1, 0), nrow = 2))
    # w[1] = 5e-5 (< 1e-4) _> should become 0 after processing
    sol <- param_to_solution(5e-5, gc, process = TRUE)
    expect_identical(sol$hyp_weight[1], 0)
    expect_equal(sum(sol$hyp_weight), 1, tolerance = 1e-14)
})

test_that("param_to_solution(process=TRUE) snaps weights to epsilon", {
    # snaps weights in (1e-4, 1e-3) to epsilon (m=2)
    gc <- graph_constraint(c(NA, NA), matrix(c(0, 1, 1, 0), nrow = 2))
    # w[1] = 5e-4 (between 1e-4 and 1e-3) -> should become 0.001
    sol <- param_to_solution(5e-4, gc, process = TRUE)
    expect_equal(sol$hyp_weight[1] - 0.001, 0, tolerance = 1e-6)
    expect_equal(sum(sol$hyp_weight) - 1, 0, tolerance = 1e-6)
})

test_that("param_to_solution(process=TRUE) doesn't change large weights(m=2)", {
    gc <- graph_constraint(c(NA, NA), matrix(c(0, 1, 1, 0), nrow = 2))
    # w[1] = 0.4 (well above 1e-3) -> should stay at 0.4
    sol_no <- param_to_solution(0.4, gc, process = FALSE)
    sol_yes <- param_to_solution(0.4, gc, process = TRUE)
    expect_identical(sol_no$hyp_weight, sol_yes$hyp_weight)
})

test_that("param_to_solution(process=TRUE) thresholds G matrix entries (m=3)", {
    # m=3, all free: hyp_constraint = c(NA, NA, NA)
    # trans_constraint: diag = 0, off-diag = NA
    gc <- graph_constraint_free(3)
    # For m=3 free: w_len = 2, g_len = 3 (one free param per row of G).
    # Row i: free cols are off-diag NAs. With 2 free off-diag per row and
    # 1 derived _> 1 param per row _> 3 G params total.
    # Total params = 2 + 3 = 5.
    #
    # Weights: w = c(0.5, 0.3, 0.2)  via  w_pars = c(0.5, 0.3)
    # G row 1 (cols 2,3 free): g_par = 5e-6 _> G[1,2]=5e-6, G[1,3]=1-5e-6
    # G row 2 (cols 1,3 free): g_par = 0.5 _> G[2,1]=0.5, G[2,3]=0.5
    # G row 3 (cols 1,2 free): g_par = 0.5 _> G[3,1]=0.5, G[3,2]=0.5
    optim_params <- c(0.5, 0.3, 5e-6, 0.5, 0.5)
    sol <- param_to_solution(optim_params, gc, process = TRUE)
    # G[1,2] was 5e-6 (< 1e-5) _> zeroed; row re-normalised
    expect_identical(sol$trans_matrix[1, 2], 0)
    # G[1,3] should absorb all row mass
    expect_equal(sol$trans_matrix[1, 3], 1, tolerance = 1e-14)
    # Diagonal stays 0
    expect_identical(sol$trans_matrix[1, 1], 0)
    # Each row sums to 1
    for (i in 1:3) {
        expect_equal(sum(sol$trans_matrix[i, ]), 1, tolerance = 1e-14)
    }
})

test_that("param_to_solution(process=TRUE) snaps G entries to epsilon", {
    # snaps G entries in (1e-5, 1e-3) to epsilon (m=3)
    gc <- graph_constraint_free(3)
    # G row 1: g_par = 5e-4 _> G[1,2]=5e-4, G[1,3]=1-5e-4
    # G[1,2] is in (1e-5, 1e-3) _> should snap to 0.001
    optim_params <- c(0.5, 0.3, 5e-4, 0.5, 0.5)
    sol <- param_to_solution(optim_params, gc, process = TRUE)
    expect_equal(sol$trans_matrix[1, 2] - 0.001, 0, tolerance = 1e-6)
})

test_that("param_to_solution(process=TRUE) works for m=4", {
    gc <- graph_constraint_free(4)
    # m=4 free: w_len = 3, G has 2 free params per row _> 8 G params.
    # Total is 3 + 8 = 11
    # Weights: w_pars = c(0.4, 0.3, 0.2) _> w = c(0.4, 0.3, 0.2, 0.1)
    # G: each row has 3 off-diag, 2 free params + 1 derived
    # Row 1: g = c(0.3, 0.3) _> G[1,] = c(0, 0.3, 0.3, 0.4)
    # Row 2: g = c(5e-5, 0.4) _> G[2,] = c(0, 0, 5e-5, 0.4)... actually
    w_pars <- c(0.4, 0.3, 0.2) # w = (0.4, 0.3, 0.2, 0.1)
    g_pars <- c(0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3)
    optim_params <- c(w_pars, g_pars)
    sol_proc <- param_to_solution(optim_params, gc, process = TRUE)
    sol_raw <- param_to_solution(optim_params, gc, process = FALSE)
    # All weights are well above 1e-3, so they should match
    expect_equal(sol_proc$hyp_weight, sol_raw$hyp_weight, tolerance = 1e-14)
    # Weights sum to 1
    expect_equal(sum(sol_proc$hyp_weight), 1, tolerance = 1e-14)
    # Diagonals are 0
    expect_identical(diag(sol_proc$trans_matrix), c(0, 0, 0, 0))
    # All rows sum to 1
    for (i in 1:4) {
        expect_equal(sum(sol_proc$trans_matrix[i, ]), 1, tolerance = 1e-14)
    }
})


test_that("param_to_solution(process=TRUE) boundary: weight exactly at 1e-4", {
    gc <- graph_constraint(c(NA, NA), matrix(c(0, 1, 1, 0), nrow = 2))
    # w[1] = 1e-4: condition is w < 1e-4 (strict), so 1e-4 does NOT get zeroed.
    # Then w > 1e-4 & w < 1e-3: 1e-4 is NOT > 1e-4, so it stays at 1e-4.
    sol <- param_to_solution(1e-4, gc, process = TRUE)
    expect_equal(sol$hyp_weight[1] - 1e-4, 0, tolerance = 1e-6)
})

test_that("param_to_solution(process=TRUE) boundary: weight exactly at 1e-3", {
    gc <- graph_constraint(c(NA, NA), matrix(c(0, 1, 1, 0), nrow = 2))
    # w[1] = 1e-3: condition w < 1e-4 is FALSE; w > 1e-4 & w < 1e-3 is FALSE
    # (1e-3 is NOT < 1e-3). So 1e-3 stays unchanged.
    sol <- param_to_solution(1e-3, gc, process = TRUE)
    expect_equal(sol$hyp_weight[1], 1e-3, tolerance = 1e-14)
})


# --- .free_recipients tests ---
test_that(".free_recipients returns correct indices", {
    expect_identical(.free_recipients(2L, integer(0), 4L), c(1L, 3L, 4L))
    expect_identical(.free_recipients(2L, c(1L, 4L), 4L), 3L)
    expect_identical(.free_recipients(1L, c(2L, 3L), 3L), integer(0))
})

test_that(".free_recipients excludes both drop and fixed", {
    # drop_idx overlaps with fixed_idx — still excluded just once
    expect_identical(.free_recipients(2L, c(2L, 3L), 4L), c(1L, 4L))
})


# --- .redistribute_mass tests ---
test_that(".redistribute_mass preserves total sum", {
    vec <- c(0.3, 0.5, 0.2)
    result <- .redistribute_mass(vec, drop_idx = 3L, fixed_idx = integer(0))
    expect_equal(sum(result), 1, tolerance = 1e-14)
    expect_identical(result[3], 0)
})

test_that(".redistribute_mass proportional redistribution", {
    vec <- c(0.6, 0.2, 0.2)
    result <- .redistribute_mass(vec, drop_idx = 3L, fixed_idx = integer(0))
    expect_identical(result[3], 0)
    # 0.6 and 0.2 in 3:1 ratio, should stay 3:1
    expect_equal(result[1] / result[2], 3, tolerance = 1e-12)
    expect_equal(sum(result), 1, tolerance = 1e-14)
})

test_that(".redistribute_mass uniform fallback when all recipients zero", {
    vec <- c(0.0, 0.0, 1.0, 0.0)
    result <- .redistribute_mass(vec, drop_idx = 3L, fixed_idx = integer(0))
    expect_identical(result[3], 0)
    # mass split equally among 1, 2, 4
    expect_equal(result[1], 1 / 3, tolerance = 1e-12)
    expect_equal(result[2], 1 / 3, tolerance = 1e-12)
    expect_equal(result[4], 1 / 3, tolerance = 1e-12)
})

test_that(".redistribute_mass respects fixed_idx", {
    vec <- c(0.4, 0.1, 0.3, 0.2)
    fixed <- c(1L, 4L)
    result <- .redistribute_mass(vec, drop_idx = 2L, fixed_idx = fixed)
    expect_identical(result[2], 0)
    # only index 3 is a free recipient
    expect_equal(result[3], 0.3 + 0.1, tolerance = 1e-12)
    expect_equal(sum(result), 1, tolerance = 1e-14)
})

test_that(".redistribute_mass with no free recipients returns vec unchanged", {
    vec <- c(0.5, 0.3, 0.2)
    # fixed = {1, 3}, drop = 2 -> no recipients
    result <- .redistribute_mass(vec, drop_idx = 2L, fixed_idx = c(1L, 3L))
    expect_identical(result[2], 0)
    # No recipients: fixed entries must not be rescaled.
    # Sum is 0.7 (not 1) — the downstream power check rejects if harmful.
    expect_equal(result[1], 0.5)
    expect_equal(result[3], 0.2)
})

test_that(".redistribute_mass works for transition rows", {
    row <- c(0, 0.6, 0.1, 0.3)
    # diagonal is index 1, treat as fixed
    result <- .redistribute_mass(row, drop_idx = 3L, fixed_idx = 1L)
    expect_identical(result[1], 0)
    expect_identical(result[3], 0)
    expect_equal(sum(result), 1, tolerance = 1e-14)
})


# --- .marginal_violated tests ---
test_that(".marginal_violated returns FALSE with no constraints", {
    expect_false(.marginal_violated(c(0.8, 0.7), integer(0), NULL))
})

test_that(".marginal_violated detects violation", {
    marg <- c(0.80, 0.60, 0.90)
    constr <- c(NA, 0.70, NA)
    idx <- which(!is.na(constr))
    expect_true(.marginal_violated(marg, idx, constr))
})

test_that(".marginal_violated passes when all satisfied", {
    marg <- c(0.80, 0.75, 0.90)
    constr <- c(NA, 0.70, 0.85)
    idx <- which(!is.na(constr))
    expect_false(.marginal_violated(marg, idx, constr))
})

test_that(".marginal_violated boundary: exactly equal passes", {
    marg <- 0.70
    constr <- 0.70
    idx <- 1L
    expect_false(.marginal_violated(marg, idx, constr))
})


# --- Setup shared across .try_prune() and prune_graph() E2E tests ---
conjunctive_4m_power <- trial_success(r1 && r2 && r3 && r4, verbose = "silent")
avg_6m_power <- trial_success(r1 + r2 + r3 + r4 + r5 + r6, verbose = "silent")

m <- 4
power_vector <- c(0.93, 0.91, 0.90, 0.85)
corr_mat <- matrix(0.2, nrow = m, ncol = m)
diag(corr_mat) <- 1

sims <- mvtnorm::rmvnorm(
    5e4,
    sigma = corr_mat,
    mean = calc_ncp(power_vector)
)
pvals <- stats::pnorm(sims, lower.tail = FALSE)

# 6-m setup for the constraint-respect test
power_vec_6 <- c(0.9, 0.8, 0.6, 0.85, 0.85, 0.5)
ncp_vec_6 <- qnorm(1 - 0.025) - qnorm(1 - power_vec_6)
corr_6 <- matrix(0.2, 6, 6)
diag(corr_6) <- 1
pvals_6m <- withr::with_seed(5, {
    s6 <- mvtnorm::rmvnorm(2^20, mean = ncp_vec_6, sigma = corr_6)
    stats::pnorm(s6, lower.tail = FALSE)
})

# --- .try_prune() unit tests ---

test_that(".try_prune() accepts when power improves and no constraints", {
    w <- c(0.5, 0.5, 0, 0)
    G <- rbind(
        c(0, 0, 0.5, 0.5),
        c(0, 0, 0.5, 0.5),
        c(0.5, 0, 0, 0.5),
        c(0, 0.5, 0.5, 0)
    )

    # Compute baseline power with a non-optimal graph
    baseline <- calc_power_pvals(
        pvals[, 1:4],
        hyp_weight = w,
        trans_matrix = G,
        custom_power = conjunctive_4m_power
    )

    # Concentrate weight on H1 — should improve conjunctive power
    w_better <- c(1, 0, 0, 0)
    result <- .try_prune(
        pvals = pvals[, 1:4],
        hyp_weight = w_better,
        trans_matrix = G,
        trial_success = conjunctive_4m_power,
        power_best = baseline$custom_power,
        constrained_idx = integer(0),
        power_constraint = NULL
    )

    expect_true(result$accepted)
    expect_gte(result$power_best, baseline$custom_power)
    expect_identical(result$hyp_weight, w_better)
    expect_identical(result$trans_matrix, G)
})

test_that(".try_prune() accepts when power exactly equals power_best", {
    w <- c(0.5, 0.5, 0, 0)
    G <- rbind(
        c(0, 0, 0.5, 0.5),
        c(0, 0, 0.5, 0.5),
        c(0.5, 0, 0, 0.5),
        c(0, 0.5, 0.5, 0)
    )

    # Compute the exact power with the same graph
    exact <- calc_power_pvals(
        pvals[, 1:4],
        hyp_weight = w,
        trans_matrix = G,
        custom_power = conjunctive_4m_power
    )

    # Pass the exact same graph with power_best == actual power
    result <- .try_prune(
        pvals = pvals[, 1:4],
        hyp_weight = w,
        trans_matrix = G,
        trial_success = conjunctive_4m_power,
        power_best = exact$custom_power,
        constrained_idx = integer(0),
        power_constraint = NULL
    )

    # Should accept because custom_power >= power_best (equal)
    expect_true(result$accepted)
    expect_identical(result$power_best, exact$custom_power)
})

test_that(".try_prune() rejects when power is lower than power_best", {
    w <- c(0.5, 0.5, 0, 0)
    G <- rbind(
        c(0, 0, 0.5, 0.5),
        c(0, 0, 0.5, 0.5),
        c(0.5, 0, 0, 0.5),
        c(0, 0.5, 0.5, 0)
    )

    # Set power_best artificially high so the candidate cannot match it
    result <- .try_prune(
        pvals = pvals[, 1:4],
        hyp_weight = w,
        trans_matrix = G,
        trial_success = conjunctive_4m_power,
        power_best = 1.0,
        constrained_idx = integer(0),
        power_constraint = NULL
    )

    expect_false(result$accepted)
    # power_best should remain unchanged at 1.0
    expect_identical(result$power_best, 1.0)
})

test_that(".try_prune() rejects when marginal constraint is violated", {
    w <- c(0.5, 0.5, 0, 0)
    G <- rbind(
        c(0, 0, 0.5, 0.5),
        c(0, 0, 0.5, 0.5),
        c(0.5, 0, 0, 0.5),
        c(0, 0.5, 0.5, 0)
    )

    # Compute baseline
    baseline <- calc_power_pvals(
        pvals[, 1:4],
        hyp_weight = w,
        trans_matrix = G,
        custom_power = conjunctive_4m_power
    )

    # Set an impossibly high marginal constraint on H4
    power_constr <- c(NA, NA, NA, 0.999)
    constrained_idx <- which(!is.na(power_constr))

    result <- .try_prune(
        pvals = pvals[, 1:4],
        hyp_weight = w,
        trans_matrix = G,
        trial_success = conjunctive_4m_power,
        power_best = 0,
        constrained_idx = constrained_idx,
        power_constraint = power_constr
    )

    # Should reject because H4 marginal power < 0.999
    expect_false(result$accepted)
    # power_best stays at the passed-in value when rejected
    expect_identical(result$power_best, 0)
})

test_that(".try_prune() returns correct list structure", {
    w <- c(0.5, 0.5, 0, 0)
    G <- rbind(
        c(0, 0, 0.5, 0.5),
        c(0, 0, 0.5, 0.5),
        c(0.5, 0, 0, 0.5),
        c(0, 0.5, 0.5, 0)
    )

    result <- .try_prune(
        pvals = pvals[, 1:4],
        hyp_weight = w,
        trans_matrix = G,
        trial_success = conjunctive_4m_power,
        power_best = 0,
        constrained_idx = integer(0),
        power_constraint = NULL
    )

    expect_named(
        result,
        c("hyp_weight", "trans_matrix", "power_best", "accepted")
    )
    expect_type(result$power_best, "double")
    expect_type(result$accepted, "logical")
    expect_identical(result$hyp_weight, w)
    expect_identical(result$trans_matrix, G)
})

test_that(".try_prune() rejects despite high power when constraint violated", {
    # Even when power improves, a violated marginal constraint causes rejection
    w_good <- c(1, 0, 0, 0)
    G <- rbind(
        c(0, 0, 0.5, 0.5),
        c(0, 0, 0.5, 0.5),
        c(0.5, 0, 0, 0.5),
        c(0, 0.5, 0.5, 0)
    )

    # Set a tight constraint on H2 — concentrating all weight on H1 means
    # H2 marginal power will be low
    power_constr <- c(NA, 0.99, NA, NA)
    constrained_idx <- which(!is.na(power_constr))

    result <- .try_prune(
        pvals = pvals[, 1:4],
        hyp_weight = w_good,
        trans_matrix = G,
        trial_success = conjunctive_4m_power,
        power_best = 0,
        constrained_idx = constrained_idx,
        power_constraint = power_constr
    )

    # Despite power_best = 0 (any power would improve), constraint violation
    # causes rejection
    expect_false(result$accepted)
})


# --- prune_hyp_weights() tests ---

test_that("`prune_hyp_weights()` prunes redundant hypothesis weights", {
    w <- c(0.5, 0.5, 0, 0)
    G <- rbind(
        c(0, 0, 0.5, 0.5),
        c(0, 0, 0.5, 0.5),
        c(0.001, 0, 0, 0.999),
        c(0, 0.001, 0.999, 0)
    )

    result <- prune_hyp_weights(
        pvals[, 1:4],
        gamma = 1,
        fixed_w = which(!is.na(graph_constraint_free(4)$hyp_constraint)),
        trial_success = conjunctive_4m_power,
        hyp_weight = w,
        trans_matrix = G
    )

    expect_identical(result$hyp_weight, c(1, 0, 0, 0))
    # trans_matrix is held fixed during weight pruning
    expect_identical(result$trans_matrix, G)
    expect_type(result$power_best, "double")
})

test_that("`prune_hyp_weights()` doesn't change hyp weights below gamma", {
    w <- c(0.5, 0.5, 0, 0)
    G <- rbind(
        c(0, 0, 0.5, 0.5),
        c(0, 0, 0.5, 0.5),
        c(0.001, 0, 0, 0.999),
        c(0, 0.001, 0.999, 0)
    )

    result <- prune_hyp_weights(
        pvals[, 1:4],
        gamma = 0.49,
        fixed_w = which(!is.na(graph_constraint_free(4)$hyp_constraint)),
        trial_success = conjunctive_4m_power,
        hyp_weight = w,
        trans_matrix = G
    )

    expect_equal(result$hyp_weight, c(0.5, 0.5, 0, 0))
})

test_that("`prune_hyp_weights()` respects graph_constraint fixed weights", {
    pvals_3m <- pvals_6m[, 1:3]

    gc <- graph_constraint(
        hyp_constraint = c(NA, NA, 0),
        trans_constraint = rbind(
            c(0, NA, NA),
            c(NA, 0, NA),
            c(0.5, 0.5, 0)
        )
    )

    w <- c(0.89, 0.11, 0.0)
    G <- rbind(
        c(0.0, 0.58, 0.42),
        c(0.36, 0.0, 0.64),
        c(0.5, 0.5, 0.0)
    )

    avg_power <- trial_success(1 / 3 * (r1 + r2 + r3), verbose = "silent")

    result <- prune_hyp_weights(
        pvals_3m,
        hyp_weight = w,
        trans_matrix = G,
        trial_success = avg_power,
        gamma = 1,
        fixed_w = which(!is.na(gc$hyp_constraint))
    )

    # H3 weight must remain 0 (fixed by constraint)
    expect_identical(result$hyp_weight[3], 0)
    # trans_matrix is unchanged
    expect_identical(result$trans_matrix, G)

    # Power must not decrease
    original_power <- calc_power_pvals(
        pvals_3m,
        hyp_weight = w,
        trans_matrix = G,
        custom_power = avg_power,
        sum_to_one_constraint = FALSE
    )
    expect_gte(result$power_best, original_power$custom_power)
})


# --- prune_edges() tests ---

test_that("`prune_edges()` prunes redundant edges", {
    w <- c(1, 0, 0, 0)
    G <- rbind(
        c(0, 0, 0.5, 0.5),
        c(0, 0, 0.5, 0.5),
        c(0.001, 0, 0, 0.999),
        c(0, 0.001, 0.999, 0)
    )

    baseline <- calc_power_pvals(
        pvals[, 1:4],
        hyp_weight = w,
        trans_matrix = G,
        custom_power = conjunctive_4m_power
    )

    result <- prune_edges(
        pvals[, 1:4],
        gamma = 1,
        fixed_edge = !is.na(graph_constraint_free(4)$trans_constraint),
        trial_success = conjunctive_4m_power,
        hyp_weight = w,
        trans_matrix = G,
        power_best = baseline$custom_power
    )

    expect_identical(
        result$trans_matrix,
        rbind(
            c(0, 0, 1, 0),
            c(0, 0, 1, 0),
            c(0, 0, 0, 1),
            c(0, 1, 0, 0)
        )
    )
    # hyp_weight is held fixed during edge pruning
    expect_identical(result$hyp_weight, w)
    expect_gte(result$power_best, baseline$custom_power)
})

test_that("`prune_edges()` doesn't change trans weights below gamma", {
    w <- c(1, 0, 0, 0)
    G <- rbind(
        c(0, 0, 0.5, 0.5),
        c(0, 0, 0.5, 0.5),
        c(0.01, 0, 0, 0.99),
        c(0, 0.01, 0.99, 0)
    )

    baseline <- calc_power_pvals(
        pvals[, 1:4],
        hyp_weight = w,
        trans_matrix = G,
        custom_power = conjunctive_4m_power
    )

    # Gamma below edge weight --> no pruning
    result <- prune_edges(
        pvals[, 1:4],
        gamma = 0.009,
        fixed_edge = !is.na(graph_constraint_free(4)$trans_constraint),
        trial_success = conjunctive_4m_power,
        hyp_weight = w,
        trans_matrix = G,
        power_best = baseline$custom_power
    )

    expect_identical(result$trans_matrix, G)
})

test_that("`prune_edges()` respects graph_constraint fixed edges", {
    pvals_3m <- pvals_6m[, 1:3]

    gc <- graph_constraint(
        hyp_constraint = c(NA, NA, 0),
        trans_constraint = rbind(
            c(0, NA, NA),
            c(NA, 0, NA),
            c(0.5, 0.5, 0)
        )
    )

    w <- c(1, 0, 0.0)
    G <- rbind(
        c(0.0, 0.58, 0.42),
        c(0.36, 0.0, 0.64),
        c(0.5, 0.5, 0.0)
    )

    avg_power <- trial_success(1 / 3 * (r1 + r2 + r3), verbose = "silent")

    baseline <- calc_power_pvals(
        pvals_3m,
        hyp_weight = w,
        trans_matrix = G,
        custom_power = avg_power,
        sum_to_one_constraint = FALSE
    )

    result <- prune_edges(
        pvals_3m,
        hyp_weight = w,
        trans_matrix = G,
        trial_success = avg_power,
        gamma = 1,
        fixed_edge = !is.na(gc$trans_constraint),
        power_best = baseline$custom_power
    )

    # Row 3 must remain exactly (0.5, 0.5, 0) — fixed by constraint
    expect_identical(result$trans_matrix[3, ], c(0.5, 0.5, 0.0))
    # Diagonals untouched
    expect_identical(result$trans_matrix[1, 1], 0)
    expect_identical(result$trans_matrix[2, 2], 0)
})


test_that("`prune_graph()` cleans redundant hyp or matrix weights", {
    w <- c(0.5, 0.5, 0, 0)
    G <- rbind(
        c(0, 0, 0.5, 0.5),
        c(0, 0, 0.5, 0.5),
        c(0.001, 0, 0, 0.999),
        c(0, 0.001, 0.999, 0)
    )

    clean <- prune_graph(
        pvals[, 1:4],
        gamma = 1,
        graph_constraint = graph_constraint_free(4),
        trial_success = conjunctive_4m_power,
        hyp_weight = w,
        trans_matrix = G,
        verbose = "silent"
    )

    expect_identical(clean$hyp_weight, c(1, 0, 0, 0))
    expect_identical(
        clean$trans_matrix,
        rbind(
            c(0, 0, 1, 0),
            c(0, 0, 1, 0),
            c(0, 0, 0, 1),
            c(0, 1, 0, 0)
        )
    )
})


test_that("`prune_graph()` doesn't change trans weights below gamma", {
    w <- c(1, 0, 0, 0)
    G <- rbind(
        c(0, 0, 0.5, 0.5),
        c(0, 0, 0.5, 0.5),
        c(0.01, 0, 0, 0.99),
        c(0, 0.01, 0.99, 0)
    )

    clean_1 <- prune_graph(
        pvals[, 1:4],
        gamma = 1,
        graph_constraint = graph_constraint_free(4),
        trial_success = conjunctive_4m_power,
        hyp_weight = w,
        trans_matrix = G,
        verbose = "silent"
    )

    expect_identical(
        clean_1$trans_matrix,
        rbind(
            c(0, 0, 1, 0),
            c(0, 0, 1, 0),
            c(0, 0, 0, 1),
            c(0, 1, 0, 0)
        )
    )

    # Gamma below edge weight --> no pruning
    clean_2 <- prune_graph(
        pvals[, 1:4],
        gamma = 0.009,
        graph_constraint = graph_constraint_free(4),
        trial_success = conjunctive_4m_power,
        hyp_weight = w,
        trans_matrix = G,
        verbose = "silent"
    )

    expect_identical(clean_2$trans_matrix, G)
})


test_that("`prune_graph()` doesn't change hyp weights below gamma", {
    w <- c(0.5, 0.5, 0, 0)
    G <- rbind(
        c(0, 0, 0.5, 0.5),
        c(0, 0, 0.5, 0.5),
        c(0.001, 0, 0, 0.999),
        c(0, 0.001, 0.999, 0)
    )

    clean <- prune_graph(
        pvals[, 1:4],
        gamma = 0.49,
        graph_constraint = graph_constraint_free(4),
        trial_success = conjunctive_4m_power,
        hyp_weight = w,
        trans_matrix = G,
        verbose = "silent"
    )

    expect_equal(clean$hyp_weight, c(0.5, 0.5, 0, 0))
})


test_that("`prune_graph()` respects graph_constraint fixed entries", {
    pvals_3m <- pvals_6m[, 1:3]

    gc <- graph_constraint(
        hyp_constraint = c(NA, NA, 0),
        trans_constraint = rbind(
            c(0, NA, NA),
            c(NA, 0, NA),
            c(0.5, 0.5, 0)
        )
    )

    w <- c(0.89, 0.11, 0.0)
    G <- rbind(
        c(0.0, 0.58, 0.42),
        c(0.36, 0.0, 0.64),
        c(0.5, 0.5, 0.0)
    )

    avg_power <- trial_success(1 / 3 * (r1 + r2 + r3), verbose = "silent")

    clean_gc <- prune_graph(
        pvals_3m,
        hyp_weight = w,
        trans_matrix = G,
        trial_success = avg_power,
        gamma = 1,
        graph_constraint = gc,
        verbose = "silent"
    )

    # H3 weight must remain 0
    expect_identical(clean_gc$hyp_weight[3], 0)
    # Row 3 must remain exactly (0.5, 0.5, 0)
    expect_identical(clean_gc$trans_matrix[3, ], c(0.5, 0.5, 0.0))
    # Diagonals untouched
    expect_identical(clean_gc$trans_matrix[1, 1], 0)
    expect_identical(clean_gc$trans_matrix[2, 2], 0)

    # Power must not decrease
    original_power <- calc_power_pvals(
        pvals_3m,
        hyp_weight = w,
        trans_matrix = G,
        custom_power = avg_power,
        sum_to_one_constraint = FALSE
    )
    pruned_power <- calc_power_pvals(
        pvals_3m,
        hyp_weight = clean_gc$hyp_weight,
        trans_matrix = clean_gc$trans_matrix,
        custom_power = avg_power,
        sum_to_one_constraint = FALSE
    )
    expect_gte(pruned_power$custom_power, original_power$custom_power)
})


## repair_graph tests

test_that("repair_graph is idempotent on a valid graph", {
    gc <- graph_constraint_free(3)
    w <- c(0.4, 0.3, 0.3)
    G <- matrix(
            c(
                0, 0.5, 0.5,
                0.3, 0, 0.7,
                0.6, 0.4, 0
            ),
            3,
            byrow = TRUE
        )
    result <- repair_graph(w, G, gc)
    expect_identical(result$hyp_weight, w)
    expect_identical(result$trans_matrix, G)
})

test_that("repair_graph clamps out-of-bound weights", {
    gc <- graph_constraint_free(3)
    w <- c(-0.001, 0.5, 0.501)
    G <- matrix(
            c(
                0, 0.5, 0.5,
                0.3, 0, 0.7,
                0.6, 0.4, 0
            ),
            3,
            byrow = TRUE
        )
    result <- repair_graph(w, G, gc)
    expect_true(all(result$hyp_weight >= 0))
    expect_true(all(result$hyp_weight <= 1))
    expect_equal(
        sum(result$hyp_weight),
        1,
        tolerance = sqrt(.Machine$double.eps)
    )
})

test_that("repair_graph normalises weights to sum to 1", {
    gc <- graph_constraint_free(3)
    w <- c(0.4, 0.3, 0.3002) # sum is 1.0002
    G <- matrix(
            c(
                0, 0.5, 0.5,
                0.3, 0, 0.7,
                0.6, 0.4, 0
            ),
            3,
            byrow = TRUE
        )
    result <- repair_graph(w, G, gc)
    expect_equal(
        sum(result$hyp_weight),
        1,
        tolerance = sqrt(.Machine$double.eps)
    )
})

test_that("repair_graph normalises each trans_matrix row", {
    gc <- graph_constraint_free(3)
    w <- c(0.4, 0.3, 0.3)
    G <- matrix(
            c(
                0, 0.5, 0.4998,
                0.3, 0, 0.7,
                0.6, 0.4, 0
            ),
            3,
            byrow = TRUE
        )
    result <- repair_graph(w, G, gc)
    for (i in 1:3) {
        expect_equal(
            sum(result$trans_matrix[i, ]),
            1,
            tolerance = sqrt(.Machine$double.eps)
        )
    }
})

test_that("repair_graph respects fixed constraints", {
    tc <- matrix(NA_real_, 3, 3)
    diag(tc) <- 0
    tc[1, 2] <- 0.3
    tc[1, 3] <- 0.7
    gc <- graph_constraint(
        hyp_constraint = c(0.5, NA, NA),
        trans_constraint = tc
    )
    w <- c(0.45, 0.3, 0.25) # first is fixed at 0.5
    G <- matrix(
            c(
                0, 0.5, 0.5,
                0.3, 0, 0.7,
                0.6, 0.4, 0
            ),
            3,
            byrow = TRUE
        )
    result <- repair_graph(w, G, gc)
    expect_equal(result$hyp_weight[1], 0.5)
    expect_equal(result$trans_matrix[1, 2], 0.3)
})

test_that("repair_graph forces diagonal to 0", {
    gc <- graph_constraint_free(3)
    w <- c(0.4, 0.3, 0.3)
    G <- matrix(
            c(
                1e-10, 0.5, 0.5,
                0.3, 1e-8, 0.7,
                0.6, 0.4, 1e-12
            ),
            3,
            byrow = TRUE
        )
    result <- repair_graph(w, G, gc)
    expect_identical(diag(result$trans_matrix), c(0, 0, 0))
})

test_that("repair_graph handles all-zero free weights", {
    tc <- matrix(NA_real_, 3, 3)
    diag(tc) <- 0
    gc <- graph_constraint(
        hyp_constraint = c(0.5, NA, NA),
        trans_constraint = tc
    )
    w <- c(0.5, 0, 0) # free weights are all zero
    G <- matrix(
            c(
                0, 0.5, 0.5,
                0.3, 0, 0.7,
                0.6, 0.4, 0
            ),
            3,
            byrow = TRUE
        )
    result <- repair_graph(w, G, gc)
    expect_identical(sum(result$hyp_weight), 1)
    expect_equal(result$hyp_weight[1], 0.5)
    expect_true(all(result$hyp_weight[2:3] > 0))
})

test_that("repair_graph handles all-zero free transition row", {
    tc <- matrix(NA_real_, 3, 3)
    diag(tc) <- 0
    gc <- graph_constraint(
        hyp_constraint = c(NA, NA, NA),
        trans_constraint = tc
    )
    w <- c(0.4, 0.3, 0.3)
    # Row 1: diagonal forced to 0, cols 2 and 3 are free but zero
    G <- matrix(
            c(
                0, 0, 0,
                0.3, 0, 0.7,
                0.6, 0.4, 0
            ),
            3,
            byrow = TRUE
        )
    result <- repair_graph(w, G, gc)
    expect_identical(sum(result$trans_matrix[1, ]), 1) # row sums to 1
    expect_true(all(result$trans_matrix[1, 2:3] > 0)) # free cols received mass
})


# --- graph_constraint tolerance flows to normalise_sum ---

# 3-hypothesis constraint with fixed weights summing to 1.005.
over_sum_gc <- function(tolerance) {
    tc <- matrix(NA_real_, 3, 3)
    diag(tc) <- 0
    new_graph_constraint(
        hyp_constraint = c(0.6, 0.405, NA_real_),
        trans_constraint = tc,
        names = c("H1", "H2", "H3"),
        tolerance = tolerance
    )
}

test_that("param_to_solution() uses the graph_constraint tolerance", {
    # 0 free weight params (both fixed) + 3 transition params.
    params <- c(0.5, 0.5, 0.5)

    expect_error(
        param_to_solution(params, over_sum_gc(sqrt(.Machine$double.eps)), TRUE),
        "exceeds target"
    )
    expect_no_error(
        param_to_solution(params, over_sum_gc(0.01), process = TRUE)
    )
})

test_that("repair_graph() uses the graph_constraint tolerance", {
    # Free weight (H3) is 0, so its mass must come from the fixed weights,
    w <- c(0.6, 0.405, 0)
    G <- rbind(c(0, 0.5, 0.5), c(0.5, 0, 0.5), c(0.5, 0.5, 0))

    expect_error(
        repair_graph(w, G, over_sum_gc(sqrt(.Machine$double.eps))),
        "exceeds target"
    )
    expect_no_error(repair_graph(w, G, over_sum_gc(0.01)))
})

test_that(".redistribute_mass() honours (and defaults) its tolerance", {
    # Fixed elements 1 and 4 sum to 1.005; dropping index 2 redistributes to 3.
    vec <- c(0.6, 0.10, 0.30, 0.405)

    # Default tolerance (sqrt(.Machine$double.eps)) rejects the over-sum.
    expect_error(
        .redistribute_mass(vec, drop_idx = 2L, fixed_idx = c(1L, 4L)),
        "exceeds target"
    )

    # A loose tolerance tolerates it.
    result <- .redistribute_mass(
        vec,
        drop_idx = 2L,
        fixed_idx = c(1L, 4L),
        tolerance = 0.01
    )
    expect_identical(result[2], 0)
})

test_that("prune_graph() forwards the graph_constraint tolerance", {
    gc <- graph_constraint(
        hyp_constraint = rep(NA_real_, 4),
        trans_constraint = trans_constraint_free(4),
        tolerance = 0.01
    )
    w <- c(0.5, 0.5, 0, 0)
    G <- rbind(
        c(0, 0, 0.5, 0.5),
        c(0, 0, 0.5, 0.5),
        c(0.001, 0, 0, 0.999),
        c(0, 0.001, 0.999, 0)
    )

    tolerances <- new.env()
    tolerances$values <- numeric(0)
    local_mocked_bindings(
        normalise_sum = function(
            x,
            ...,
            fixed_idx = integer(0),
            target = 1,
            tolerance = sqrt(.Machine$double.eps)
        ) {
            tolerances$values <- c(tolerances$values, tolerance)
            x
        }
    )

    prune_graph(
        pvals[, 1:4],
        hyp_weight = w,
        trans_matrix = G,
        trial_success = conjunctive_4m_power,
        graph_constraint = gc,
        gamma = 1,
        verbose = "silent"
    )

    expect_gt(length(tolerances$values), 0)
    expect_true(all(tolerances$values == 0.01))
})
