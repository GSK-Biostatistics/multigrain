# --------------------------------------------------------------------------
# Tests for feasibility-preserving pruning (prune_graph_n and helpers).
# A mocked trial-success object avoids Rcpp compilation; rejection matrices
# come from the precompiled graph_shortcut().
# --------------------------------------------------------------------------

mock_ts <- function(f, m = 3L) {
    structure(
        list(func = f, m = m, objective = "mock"),
        class = "multigrain_trial_success"
    )
}

disj_fn <- function(rej) mean(rowSums(rej) > 0)

# deterministic fixture: CRN-style p-values with non-trivial power, a graph
# with one small weight, and the reference quantities for prune decisions
make_prune_setup <- function() {
    pvals <- withr::with_seed(
        1,
        1 - pnorm(matrix(rnorm(400 * 3), ncol = 3) + 0.2 * sqrt(150))
    )
    w <- c(0.5, 0.45, 0.05)
    g <- matrix(
        c(
            0, 0.5, 0.5,
            0.5, 0, 0.5,
            0.5, 0.5, 0
        ),
        nrow = 3,
        byrow = TRUE
    )
    alpha <- 0.025

    rej_base <- graph_shortcut(pvals, alpha, w, g)
    # candidate after pruning the small weight w3
    w_pruned <- .redistribute_mass(w, drop_idx = 3L, fixed_idx = integer(0))
    rej_pruned <- graph_shortcut(pvals, alpha, w_pruned, g)

    list(
        pvals = pvals,
        w = w,
        g = g,
        alpha = alpha,
        gc = graph_constraint_free(3),
        lp_base = colMeans(rej_base),
        lp_pruned = colMeans(rej_pruned),
        ts_base = disj_fn(rej_base),
        ts_pruned = disj_fn(rej_pruned)
    )
}


## ------ .try_prune_n ------ ##

test_that(".try_prune_n: target semantics", {
    s <- make_prune_setup()

    # target below the achieved trial success: accepted
    res_ok <- .try_prune_n(
        pvals = s$pvals,
        hyp_weight = s$w,
        trans_matrix = s$g,
        alpha = s$alpha,
        ts_func = disj_fn,
        target = s$ts_base - 0.05,
        constrained_idx = integer(0),
        power_constraint = NULL
    )
    expect_true(res_ok$accepted)
    expect_identical(res_ok$ts_value, s$ts_base)
    expect_identical(res_ok$local_power, s$lp_base)

    # target above the achieved trial success: rejected
    res_bad <- .try_prune_n(
        pvals = s$pvals,
        hyp_weight = s$w,
        trans_matrix = s$g,
        alpha = s$alpha,
        ts_func = disj_fn,
        target = s$ts_base + 0.05,
        constrained_idx = integer(0),
        power_constraint = NULL
    )
    expect_false(res_bad$accepted)
})

test_that(".try_prune_n: marginal floors and NULL ts_func", {
    s <- make_prune_setup()
    floors <- c(NA, NA, s$lp_base[3] + 0.05) # H3 floor unattainable

    res <- .try_prune_n(
        pvals = s$pvals,
        hyp_weight = s$w,
        trans_matrix = s$g,
        alpha = s$alpha,
        ts_func = NULL,
        target = NULL,
        constrained_idx = 3L,
        power_constraint = floors
    )

    expect_false(res$accepted)
    expect_identical(res$ts_value, NA_real_)
})

test_that(".try_prune_n: num_threads does not change the decision", {
    s <- make_prune_setup()

    args <- list(
        pvals = s$pvals,
        hyp_weight = s$w,
        trans_matrix = s$g,
        alpha = s$alpha,
        ts_func = disj_fn,
        target = s$ts_base - 0.05,
        constrained_idx = integer(0),
        power_constraint = NULL
    )
    res1 <- do.call(.try_prune_n, c(args, num_threads = 1L))
    res2 <- do.call(.try_prune_n, c(args, num_threads = 2L))

    expect_identical(res1, res2)
})


## ------ .prune_hyp_weights_n ------ ##

test_that(".prune_hyp_weights_n: accepts a prune when floors survive", {
    s <- make_prune_setup()
    floors <- c(NA, NA, s$lp_pruned[3] - 0.02) # satisfied even after prune

    out <- .prune_hyp_weights_n(
        pvals = s$pvals,
        hyp_weight = s$w,
        trans_matrix = s$g,
        ts_func = NULL,
        fixed_w = integer(0),
        alpha = s$alpha,
        gamma = 0.1, # only the small w3 is a candidate
        target = NULL,
        constrained_idx = 3L,
        power_constraint = floors
    )

    expect_identical(out$hyp_weight[3], 0)
    expect_equal(sum(out$hyp_weight), 1, tolerance = 1e-12)
    expect_identical(out$trans_matrix, s$g)
})

test_that(".prune_hyp_weights_n: rejects a prune that breaks a floor", {
    s <- make_prune_setup()
    # floor sits between the pruned and unpruned marginal power of H3
    skip_if(
        s$lp_pruned[3] >= s$lp_base[3],
        "fixture did not separate the pruned and unpruned power"
    )
    floors <- c(NA, NA, (s$lp_pruned[3] + s$lp_base[3]) / 2)

    out <- .prune_hyp_weights_n(
        pvals = s$pvals,
        hyp_weight = s$w,
        trans_matrix = s$g,
        ts_func = NULL,
        fixed_w = integer(0),
        alpha = s$alpha,
        gamma = 0.1,
        target = NULL,
        constrained_idx = 3L,
        power_constraint = floors
    )

    expect_identical(out$hyp_weight, s$w)
})

test_that(".prune_hyp_weights_n: fixed weights are never pruned", {
    s <- make_prune_setup()
    floors <- c(0.1, NA, NA) # easily satisfied

    out <- .prune_hyp_weights_n(
        pvals = s$pvals,
        hyp_weight = s$w,
        trans_matrix = s$g,
        ts_func = NULL,
        fixed_w = 3L, # w3 fixed by the constraint
        alpha = s$alpha,
        gamma = 0.1,
        target = NULL,
        constrained_idx = 1L,
        power_constraint = floors
    )

    expect_identical(out$hyp_weight, s$w)
})


## ------ .prune_edges_n ------ ##

test_that(".prune_edges_n: prunes free edges while floors survive", {
    s <- make_prune_setup()
    fixed_edge <- !is.na(s$gc$trans_constraint) # only the diagonal is fixed
    floors <- c(0.05, NA, NA) # very low: every prune stays feasible

    out <- .prune_edges_n(
        pvals = s$pvals,
        hyp_weight = s$w,
        trans_matrix = s$g,
        ts_func = NULL,
        fixed_edge = fixed_edge,
        alpha = s$alpha,
        gamma = 1,
        target = NULL,
        constrained_idx = 1L,
        power_constraint = floors
    )

    expect_identical(out$hyp_weight, s$w)
    # some edge was pruned, rows still sum to 1, diagonal untouched
    expect_gt(sum(out$trans_matrix == 0), sum(s$g == 0))
    expect_equal(
        unname(rowSums(out$trans_matrix)),
        rep(1, 3),
        tolerance = 1e-12
    )
    expect_identical(unname(diag(out$trans_matrix)), rep(0, 3))
})

test_that(".prune_edges_n: fixed edges are never pruned", {
    s <- make_prune_setup()
    fixed_edge <- !is.na(s$gc$trans_constraint)
    fixed_edge[1, 2] <- TRUE # pin one off-diagonal edge
    floors <- c(0.05, NA, NA)

    out <- .prune_edges_n(
        pvals = s$pvals,
        hyp_weight = s$w,
        trans_matrix = s$g,
        ts_func = NULL,
        fixed_edge = fixed_edge,
        alpha = s$alpha,
        gamma = 1,
        target = NULL,
        constrained_idx = 1L,
        power_constraint = floors
    )

    expect_identical(out$trans_matrix[1, 2], s$g[1, 2])
})


## ------ prune_graph_n ------ ##

test_that("prune_graph_n: guards against vacuous pruning", {
    s <- make_prune_setup()

    expect_error(
        prune_graph_n(
            pvals = s$pvals,
            hyp_weight = s$w,
            trans_matrix = s$g,
            graph_constraint = s$gc
        ),
        "at least one active constraint"
    )
    expect_error(
        prune_graph_n(
            pvals = s$pvals,
            hyp_weight = s$w,
            trans_matrix = s$g,
            graph_constraint = s$gc,
            target = NULL,
            power_constraint = c(NA, NA, NA)
        ),
        "at least one active constraint"
    )
    expect_error(
        prune_graph_n(
            pvals = s$pvals,
            hyp_weight = s$w,
            trans_matrix = s$g,
            graph_constraint = s$gc,
            target = 0.8,
            trial_success = NULL
        ),
        "does not provide a function"
    )
})

test_that("prune_graph_n: prunes while preserving feasibility", {
    s <- make_prune_setup()
    ts <- mock_ts(disj_fn)
    target <- s$ts_pruned - 0.05 # feasible even after pruning w3

    out <- prune_graph_n(
        pvals = s$pvals,
        hyp_weight = s$w,
        trans_matrix = s$g,
        graph_constraint = s$gc,
        trial_success = ts,
        target = target,
        alpha = s$alpha
    )

    expect_named(
        out,
        c(
            "hyp_weight",
            "trans_matrix",
            "trial_success_value",
            "local_power",
            "feasible"
        )
    )
    expect_true(out$feasible)
    expect_identical(out$hyp_weight[3], 0) # the small weight was pruned
    expect_equal(sum(out$hyp_weight), 1, tolerance = 1e-12)
    expect_equal(
        unname(rowSums(out$trans_matrix)),
        rep(1, 3),
        tolerance = 1e-12
    )
    expect_gte(out$trial_success_value, target)
    expect_length(out$local_power, 3L)
})

test_that("prune_graph_n: returned graph always meets the target", {
    s <- make_prune_setup()
    ts <- mock_ts(disj_fn)
    # tight target: only prunes that keep trial success at the current level
    # can be accepted, and the result must still report feasibility
    target <- s$ts_base

    out <- prune_graph_n(
        pvals = s$pvals,
        hyp_weight = s$w,
        trans_matrix = s$g,
        graph_constraint = s$gc,
        trial_success = ts,
        target = target,
        alpha = s$alpha
    )

    expect_true(out$feasible)
    expect_gte(out$trial_success_value, target)
})

test_that("prune_graph_n: marginal-only mode and num_threads parity", {
    s <- make_prune_setup()
    floors <- c(s$lp_pruned[1] - 0.05, NA, NA)

    out1 <- prune_graph_n(
        pvals = s$pvals,
        hyp_weight = s$w,
        trans_matrix = s$g,
        graph_constraint = s$gc,
        power_constraint = floors,
        alpha = s$alpha,
        num_threads = 1L
    )
    out2 <- prune_graph_n(
        pvals = s$pvals,
        hyp_weight = s$w,
        trans_matrix = s$g,
        graph_constraint = s$gc,
        power_constraint = floors,
        alpha = s$alpha,
        num_threads = 2L
    )

    expect_identical(out1, out2)
    expect_true(out1$feasible)
    expect_identical(out1$trial_success_value, NA_real_)
    expect_gte(out1$local_power[1], floors[1])
})
