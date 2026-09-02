# --------------------------------------------------------------------------
# Unit tests for graph_optimise_n() internals (R/optimization_n.R).
# These tests avoid trial_success() so no Rcpp compilation is triggered;
# plain R functions / mocked objects stand in for compiled criteria.
# Integration tests live in test-optimization_n-integration.R.
# --------------------------------------------------------------------------

## ------ Shared helpers ------ ##

make_test_setup <- function(m = 3L, nsim = 500L) {
    m <- as.integer(m)
    nsim <- as.integer(nsim)
    constraint <- graph_constraint_free(m)
    z_null <- withr::with_seed(42, matrix(rnorm(nsim * m), ncol = m))
    list(
        gc = constraint,
        z_null = z_null,
        effect_size = rep(0.2, m),
        alpha = 0.025,
        m = m,
        nsim = nsim
    )
}

# fixed-sequence graph for m = 3: test H1, pass to H2, then H3
fixed_seq_graph_3m <- function() {
    g <- matrix(0, nrow = 3, ncol = 3)
    g[1, 2] <- 1
    g[2, 3] <- 1
    list(hyp_weight = c(1, 0, 0), trans_matrix = g)
}


## ------ .make_pvalue_generator_n ------ ##

test_that(".make_pvalue_generator_n: correct dimensions and n attribute", {
    s <- make_test_setup()
    gen <- .make_pvalue_generator_n(s$z_null, s$effect_size)
    pvals <- gen(100L)

    expect_identical(dim(pvals), c(s$nsim, s$m))
    expect_identical(attr(pvals, "n"), 100L)
    expect_true(all(pvals >= 0 & pvals <= 1))
})

test_that(".make_pvalue_generator_n: larger n gives smaller p-values", {
    s <- make_test_setup()
    gen <- .make_pvalue_generator_n(s$z_null, s$effect_size)

    expect_lt(mean(gen(200L)), mean(gen(50L)))
})

test_that(".make_pvalue_generator_n: CRN determinism across generators", {
    s <- make_test_setup()
    gen1 <- .make_pvalue_generator_n(s$z_null, s$effect_size)
    gen2 <- .make_pvalue_generator_n(s$z_null, s$effect_size)

    expect_identical(gen1(100L), gen2(100L))
})

test_that(".make_pvalue_generator_n: matches the manual formula", {
    s <- make_test_setup(m = 2, nsim = 100)
    gen <- .make_pvalue_generator_n(s$z_null, s$effect_size)
    pvals <- gen(100L)

    expected <- 1 -
        pnorm(
            sweep(s$z_null, 2, s$effect_size * sqrt(100), "+")
        )

    expect_equal(pvals, expected, tolerance = 1e-12, ignore_attr = TRUE)
})

test_that(".make_pvalue_generator_n: FIFO cache eviction", {
    s <- make_test_setup(m = 2, nsim = 50)
    gen <- .make_pvalue_generator_n(s$z_null, s$effect_size, cache_size = 2L)
    info <- attr(gen, "cache_info")

    gen(10L)
    gen(20L)
    expect_identical(info()$keys, c("10", "20"))

    # a third n evicts the oldest entry (FIFO)
    gen(30L)
    expect_identical(info()$keys, c("20", "30"))
    expect_identical(info()$size, 2L)
    expect_identical(info()$max_size, 2L)

    # a cache hit does not reorder (FIFO, not LRU)
    gen(20L)
    expect_identical(info()$keys, c("20", "30"))

    # re-requesting the evicted n recomputes and evicts again
    gen(10L)
    expect_identical(info()$keys, c("30", "10"))
})

test_that(".make_pvalue_generator_n: cache_size = 0 disables caching", {
    s <- make_test_setup(m = 2, nsim = 50)
    gen <- .make_pvalue_generator_n(s$z_null, s$effect_size, cache_size = 0L)
    info <- attr(gen, "cache_info")

    p1 <- gen(10L)
    p2 <- gen(10L)
    expect_identical(info()$size, 0L)
    expect_identical(p1, p2) # recomputed but identical (same CRN)
})

test_that(".make_pvalue_generator_n: input validation", {
    s <- make_test_setup(m = 2, nsim = 50)

    expect_error(
        .make_pvalue_generator_n(c(1, 2, 3), s$effect_size),
        "must be a matrix"
    )
    expect_error(
        .make_pvalue_generator_n(s$z_null, "a"),
        "must be numeric"
    )
    expect_error(
        .make_pvalue_generator_n(s$z_null, rep(0.2, 5)),
        "must match"
    )
    expect_error(
        .make_pvalue_generator_n(s$z_null, s$effect_size, cache_size = -1L),
        "non-negative"
    )

    gen <- .make_pvalue_generator_n(s$z_null, s$effect_size)
    expect_error(gen(c(10, 20)), "single positive number")
    expect_error(gen(-5), "single positive number")
})


## ------ .binary_search_n ------ ##

# build a deterministic monotone predicate: feasible iff n >= threshold
make_threshold_predicate <- function(threshold) {
    function(n) {
        feasible <- n >= threshold
        list(
            feasible = feasible,
            best_idx = 1L,
            f_best = if (feasible) 0.5 else -0.5,
            best_theta = matrix(as.numeric(n), nrow = 1)
        )
    }
}

test_that(".binary_search_n: finds the exact feasibility boundary", {
    result <- .binary_search_n(
        make_threshold_predicate(100L),
        n_range = c(2L, 200L)
    )

    expect_identical(result$n, 100L)
    expect_identical(result$best_idx, 1L)
    expect_identical(result$f_best, 0.5)
    expect_identical(result$best_theta, matrix(100, nrow = 1))
    expect_s3_class(result$history, "data.frame")
    expect_named(result$history, c("n", "feasible", "best_idx", "f_best"))
    # first probe is always the upper bound
    expect_identical(result$history$n[1], 200L)
    expect_true(result$history$feasible[1])
})

test_that(".binary_search_n: everything feasible returns the lower bound", {
    result <- .binary_search_n(
        make_threshold_predicate(0L),
        n_range = c(2L, 1000L)
    )

    expect_identical(result$n, 2L)
})

test_that(".binary_search_n: degenerate range returns immediately", {
    result <- .binary_search_n(
        make_threshold_predicate(0L),
        n_range = c(2L, 2L)
    )

    expect_identical(result$n, 2L)
    expect_identical(nrow(result$history), 1L)
})

test_that(".binary_search_n: aborts when the upper bound is infeasible", {
    expect_error(
        .binary_search_n(
            make_threshold_predicate(300L),
            n_range = c(2L, 200L)
        ),
        "No feasible seed graph"
    )
})

test_that(".binary_search_n: guards against non-monotone feasibility", {
    # stateful predicate: the upper bound is feasible only on its first
    # evaluation, so the final re-check of `hi` must catch the inconsistency
    pred_state <- new.env(parent = emptyenv())
    assign("calls", 0L, envir = pred_state)
    pred <- function(n) {
        calls <- get("calls", envir = pred_state) + 1L
        assign("calls", calls, envir = pred_state)
        feasible <- n >= 100L && calls <= 1L
        list(
            feasible = feasible,
            best_idx = 1L,
            f_best = if (feasible) 0.5 else -0.5,
            best_theta = matrix(as.numeric(n), nrow = 1)
        )
    }

    expect_error(
        .binary_search_n(pred, n_range = c(2L, 100L)),
        "not monotone"
    )
})


## ------ .create_fitness_min_n ------ ##

# decode an encoded parameter vector with the same helpers the closure uses
decode_params <- function(x, gc) {
    th <- split_theta(as.numeric(x), gc$hyp_constraint)
    list(
        w = recover_full_weights(th$w_pars, gc$hyp_constraint),
        g = recover_full_trans_matrix(th$g_pars, gc$trans_constraint)
    )
}

disj_power_fn <- function(rej) mean(rowSums(rej) > 0)
conj_power_fn <- function(rej) mean(rowSums(rej) == ncol(rej))

test_that(".create_fitness_min_n: requires an active constraint at creation", {
    s <- make_test_setup()

    expect_error(
        .create_fitness_min_n(
            alpha = s$alpha,
            hyp_constraint = s$gc$hyp_constraint,
            trans_constraint = s$gc$trans_constraint
        ),
        "At least one of"
    )
    expect_error(
        .create_fitness_min_n(
            alpha = s$alpha,
            hyp_constraint = s$gc$hyp_constraint,
            trans_constraint = s$gc$trans_constraint,
            target = 0.8
        ),
        "not a function"
    )
    # all-NA floors do not count as an active constraint
    expect_error(
        .create_fitness_min_n(
            alpha = s$alpha,
            hyp_constraint = s$gc$hyp_constraint,
            trans_constraint = s$gc$trans_constraint,
            local_power_target = c(NA, NA, NA)
        ),
        "At least one of"
    )
})

test_that(".create_fitness_min_n: -1-shifted boundary penalties are exact", {
    s <- make_test_setup()
    gen <- .make_pvalue_generator_n(s$z_null, s$effect_size)
    pvals <- gen(100L)

    fit <- .create_fitness_min_n(
        alpha = s$alpha,
        hyp_constraint = s$gc$hyp_constraint,
        trans_constraint = s$gc$trans_constraint,
        power_criterion = disj_power_fn,
        target = 0.8
    )

    x0 <- as.numeric(create_start_params(s$gc))

    # corrupt a weight parameter: recovered weights contain negatives
    x_bad_w <- x0
    x_bad_w[1] <- 5
    dec <- decode_params(x_bad_w, s$gc)
    expect_true(any(dec$w < 0))
    expect_identical(
        fit(x_bad_w, pvals),
        -1 + sum(dec$w[dec$w < 0])
    )

    # corrupt a transition parameter: recovered row contains a negative
    x_bad_g <- x0
    x_bad_g[3] <- 1.4
    dec_g <- decode_params(x_bad_g, s$gc)
    expect_false(any(dec_g$w < 0))
    expect_true(any(dec_g$g < 0))
    expect_identical(
        fit(x_bad_g, pvals),
        -1 + sum(dec_g$g[dec_g$g < 0])
    )

    # NA / NaN parameters return the hard penalty
    x_nan <- x0
    x_nan[1] <- NaN
    expect_identical(fit(x_nan, pvals), -1e6)
})

test_that(".create_fitness_min_n: invalid always ranks below valid", {
    s <- make_test_setup()
    gen <- .make_pvalue_generator_n(s$z_null, s$effect_size)
    pvals <- gen(10L) # tiny n: deeply infeasible for a 0.999 target

    fit <- .create_fitness_min_n(
        alpha = s$alpha,
        hyp_constraint = s$gc$hyp_constraint,
        trans_constraint = s$gc$trans_constraint,
        power_criterion = conj_power_fn,
        target = 0.999
    )

    x0 <- as.numeric(create_start_params(s$gc))
    f_valid <- fit(x0, pvals)

    x_bad <- x0
    x_bad[1] <- 1.01 # mildly invalid encoding
    f_invalid <- fit(x_bad, pvals)

    # valid graphs are bounded below by -1; invalid ones bounded above by -1
    expect_gte(f_valid, -1)
    expect_lte(f_invalid, -1)
    expect_lt(f_invalid, f_valid)
})

test_that(".create_fitness_min_n: fitness equals min(slacks) exactly", {
    s <- make_test_setup(m = 3, nsim = 500)
    gen <- .make_pvalue_generator_n(s$z_null, s$effect_size)
    pvals <- gen(300L)

    x0 <- as.numeric(create_start_params(s$gc))
    dec <- decode_params(x0, s$gc)
    rej <- graph_shortcut(pvals, s$alpha, dec$w, dec$g)
    local_power <- colMeans(rej)

    # trial-success only
    fit_ts <- .create_fitness_min_n(
        alpha = s$alpha,
        hyp_constraint = s$gc$hyp_constraint,
        trans_constraint = s$gc$trans_constraint,
        power_criterion = disj_power_fn,
        target = 0.3
    )
    expect_identical(fit_ts(x0, pvals), disj_power_fn(rej) - 0.3)

    # marginal floors only (NA = no floor on H2)
    floors <- c(0.5, NA, 0.4)
    fit_marg <- .create_fitness_min_n(
        alpha = s$alpha,
        hyp_constraint = s$gc$hyp_constraint,
        trans_constraint = s$gc$trans_constraint,
        local_power_target = floors
    )
    expect_identical(
        fit_marg(x0, pvals),
        min(local_power[c(1, 3)] - c(0.5, 0.4))
    )

    # combined: the minimum over all active slacks (no trial-success bonus)
    fit_both <- .create_fitness_min_n(
        alpha = s$alpha,
        hyp_constraint = s$gc$hyp_constraint,
        trans_constraint = s$gc$trans_constraint,
        power_criterion = disj_power_fn,
        target = 0.3,
        local_power_target = floors
    )
    expect_identical(
        fit_both(x0, pvals),
        min(
            disj_power_fn(rej) - 0.3,
            local_power[c(1, 3)] - c(0.5, 0.4)
        )
    )
})

test_that(".create_fitness_min_n: sign encodes feasibility", {
    s <- make_test_setup(m = 3, nsim = 500)
    gen <- .make_pvalue_generator_n(s$z_null, s$effect_size)
    x0 <- as.numeric(create_start_params(s$gc))

    fit_easy <- .create_fitness_min_n(
        alpha = s$alpha,
        hyp_constraint = s$gc$hyp_constraint,
        trans_constraint = s$gc$trans_constraint,
        power_criterion = disj_power_fn,
        target = 0.3
    )
    # large n: easily feasible for a disjunctive 0.3 target
    expect_gt(fit_easy(x0, gen(500L)), 0)

    fit_hard <- .create_fitness_min_n(
        alpha = s$alpha,
        hyp_constraint = s$gc$hyp_constraint,
        trans_constraint = s$gc$trans_constraint,
        power_criterion = conj_power_fn,
        target = 0.999
    )
    # tiny n: infeasible for a 0.999 conjunctive target
    expect_lt(fit_hard(x0, gen(10L)), 0)
})

test_that(".create_fitness_min_n: num_threads does not change the fitness", {
    s <- make_test_setup(m = 3, nsim = 500)
    gen <- .make_pvalue_generator_n(s$z_null, s$effect_size)
    pvals <- gen(200L)
    x0 <- as.numeric(create_start_params(s$gc))

    args <- list(
        alpha = s$alpha,
        hyp_constraint = s$gc$hyp_constraint,
        trans_constraint = s$gc$trans_constraint,
        power_criterion = disj_power_fn,
        target = 0.5
    )
    fit_serial <- do.call(.create_fitness_min_n, c(args, num_threads = 1L))
    fit_parallel <- do.call(.create_fitness_min_n, c(args, num_threads = 2L))

    expect_identical(fit_serial(x0, pvals), fit_parallel(x0, pvals))
})


## ------ .try_step_down_n / .execute_stepdown_n ------ ##

test_that(".try_step_down_n and .execute_stepdown_n: synthetic fitness", {
    s <- make_test_setup(m = 2, nsim = 10)
    gen <- .make_pvalue_generator_n(s$z_null, s$effect_size)

    # feasible iff n >= 80, judged purely from the p-value matrix attribute
    fit <- function(x, pvals) {
        attr(pvals, "n") - 79.5
    }

    expect_true(.try_step_down_n(numeric(0), gen, fit, n_current = 81L))
    expect_false(.try_step_down_n(numeric(0), gen, fit, n_current = 80L))
    # at the floor no step-down is attempted
    expect_false(
        .try_step_down_n(numeric(0), gen, fit, n_current = 80L, n_min = 80L)
    )

    expect_identical(
        .execute_stepdown_n(100L, numeric(0), gen, fit),
        80L
    )
    # n_min floor is respected even when lower n would be feasible
    expect_identical(
        .execute_stepdown_n(100L, numeric(0), gen, fit, n_min = 90L),
        90L
    )
    # an already-minimal n is returned unchanged
    expect_identical(
        .execute_stepdown_n(80L, numeric(0), gen, fit),
        80L
    )
})

## ------ .run_cauchy_mutation_block ------ ##

test_that(".run_cauchy_mutation_block: structure and frozen p-values", {
    s <- make_test_setup(m = 2, nsim = 20)
    gen <- .make_pvalue_generator_n(s$z_null, s$effect_size)

    # record the n attribute of every p-value matrix the fitness sees
    fit_state <- new.env(parent = emptyenv())
    assign("seen_n", integer(0), envir = fit_state)
    fit <- function(x, pvals) {
        assign(
            "seen_n",
            c(get("seen_n", envir = fit_state), attr(pvals, "n")),
            envir = fit_state
        )
        -sum((x - 0.7)^2)
    }

    out <- withr::with_seed(
        1,
        .run_cauchy_mutation_block(
            n = 50L,
            x0 = matrix(runif(12), ncol = 3),
            get_pvals = gen,
            fitness_func = fit,
            block_config = list(
                kappa = 10L,
                mu = 4L,
                lambda = 6L,
                poptim = 0.5
            )
        )
    )

    expect_named(
        out,
        c(
            "best_x",
            "best_f",
            "parents",
            "f_par",
            "history",
            "n",
            "gens_run",
            "local_searches"
        )
    )
    expect_identical(out$gens_run, 10L)
    expect_identical(out$n, 50L)
    expect_true(all(out$best_x >= 0 & out$best_x <= 1))
    expect_gte(out$local_searches, 0L)
    # the block evaluates against a single frozen p-value matrix
    expect_identical(unique(get("seen_n", envir = fit_state)), 50L)
})


## ------ .phase2_cauchy_stepdown (mocked blocks) ------ ##

# block mock factory: best_f follows `f_seq` (recycled), gens_run = kappa
make_block_mock <- function(f_seq, kappa = 5L, mu = 4L) {
    state <- new.env(parent = emptyenv())
    assign("counter", 0L, envir = state)
    assign("captured", list(), envir = state)
    fn <- function(
        n,
        x0,
        get_pvals,
        fitness_func,
        block_config = list(),
        verbose = FALSE,
        projector = .proj_box01
    ) {
        counter <- get("counter", envir = state) + 1L
        assign("counter", counter, envir = state)
        captured <- get("captured", envir = state)
        captured[[counter]] <- x0
        assign("captured", captured, envir = state)
        d <- ncol(x0)
        f_val <- f_seq[(counter - 1L) %% length(f_seq) + 1L]
        list(
            best_x = rep(0.5, d),
            best_f = f_val,
            parents = matrix(seq_len(mu * d) / (mu * d), nrow = mu, ncol = d),
            f_par = rev(seq_len(mu)), # row 1 is the best parent
            history = data.frame(
                gen = seq_len(kappa),
                f_best = f_val,
                f_mean = f_val,
                f_median = f_val
            ),
            n = n,
            gens_run = kappa,
            local_searches = 0L
        )
    }
    list(
        fn = fn,
        calls = function() get("counter", envir = state),
        captured = function() get("captured", envir = state)
    )
}

dummy_get_pvals <- function(n) matrix(0.5, nrow = 2, ncol = 2)
never_feasible <- function(x, pvals) -1

test_that(".phase2_cauchy_stepdown: steady improvement runs all blocks", {
    mock <- make_block_mock(f_seq = 1:7, kappa = 5L)
    local_mocked_bindings(.run_cauchy_mutation_block = mock$fn)

    out <- .phase2_cauchy_stepdown(
        n_init = 100L,
        x0 = matrix(0.5, nrow = 4, ncol = 2),
        get_pvals = dummy_get_pvals,
        fitness_func = never_feasible, # step-down never accepted
        control_n = list(
            kappa = 5L,
            patience = 100L,
            max_blocks = 7L,
            mu = 4L
        )
    )

    expect_identical(mock$calls(), 7L)
    expect_identical(out$total_gens, 35L)
    expect_identical(out$n_final, 100L)
    expect_identical(out$best_f, 7L)
    expect_identical(nrow(out$history), 7L)
    expect_named(
        out$history,
        c("block", "n", "gens", "f_best_block", "f_best_global", "step_down")
    )
    expect_identical(out$step_downs, integer(0))
})

test_that(".phase2_cauchy_stepdown: patience is counted in generations", {
    # constant best_f: block 1 improves from -Inf, then stagnation
    mock <- make_block_mock(f_seq = 1, kappa = 5L)
    local_mocked_bindings(.run_cauchy_mutation_block = mock$fn)

    out <- .phase2_cauchy_stepdown(
        n_init = 100L,
        x0 = matrix(0.5, nrow = 4, ncol = 2),
        get_pvals = dummy_get_pvals,
        fitness_func = never_feasible,
        control_n = list(
            kappa = 5L,
            patience = 10L,
            max_blocks = 100L,
            mu = 4L
        )
    )

    # blocks 2 and 3 accumulate 5 + 5 = 10 >= patience generations
    expect_identical(mock$calls(), 3L)
    expect_identical(out$total_gens, 15L)
})

test_that(".phase2_cauchy_stepdown: step-down reseeds the population", {
    mock <- make_block_mock(f_seq = 1, kappa = 5L, mu = 4L)
    step_state <- new.env(parent = emptyenv())
    assign("calls", 0L, envir = step_state)
    mock_step <- function(
        n_current,
        candidate,
        get_pvals,
        fitness_func,
        n_min = 2L,
        verbose = FALSE
    ) {
        calls <- get("calls", envir = step_state) + 1L
        assign("calls", calls, envir = step_state)
        if (calls == 1L) n_current - 5L else n_current
    }
    local_mocked_bindings(
        .run_cauchy_mutation_block = mock$fn,
        .execute_stepdown_n = mock_step
    )

    out <- .phase2_cauchy_stepdown(
        n_init = 100L,
        x0 = matrix(0.5, nrow = 4, ncol = 2),
        get_pvals = dummy_get_pvals,
        fitness_func = never_feasible,
        control_n = list(
            kappa = 5L,
            patience = 10L,
            max_blocks = 100L,
            mu = 4L
        )
    )

    expect_identical(out$n_final, 95L)
    expect_identical(out$step_downs, 95L)
    # history records n after the step-down and never increases
    expect_identical(out$history$n[1], 95L)
    expect_true(all(diff(out$history$n) <= 0))
    expect_true(out$history$step_down[1])

    # the reseeded population for block 2: global best first, mu rows total
    pops <- mock$captured()
    expect_gte(length(pops), 2L)
    expect_identical(nrow(pops[[2]]), 4L)
    expect_identical(pops[[2]][1, ], rep(0.5, 2))
})

test_that(".phase2_cauchy_stepdown: check_stepdown = FALSE skips step-down", {
    mock <- make_block_mock(f_seq = 1, kappa = 5L)
    step_state <- new.env(parent = emptyenv())
    assign("called", FALSE, envir = step_state)
    mock_step <- function(...) {
        assign("called", TRUE, envir = step_state)
        2L
    }
    local_mocked_bindings(
        .run_cauchy_mutation_block = mock$fn,
        .execute_stepdown_n = mock_step
    )

    out <- .phase2_cauchy_stepdown(
        n_init = 100L,
        x0 = matrix(0.5, nrow = 4, ncol = 2),
        get_pvals = dummy_get_pvals,
        fitness_func = never_feasible,
        control_n = list(
            kappa = 5L,
            patience = 10L,
            max_blocks = 100L,
            mu = 4L,
            check_stepdown = FALSE
        )
    )

    expect_false(get("called", envir = step_state))
    expect_identical(out$n_final, 100L)
})


## ------ graph_optimise_n: validation ------ ##

mock_trial_success <- function(f, m = 3L, objective = "r1 || r2 || r3") {
    structure(
        list(func = f, m = m, objective = objective),
        class = "multigrain_trial_success"
    )
}

test_that("graph_optimise_n: exported alongside its US alias", {
    expect_type(graph_optimise_n, "closure")
    expect_type(graph_optimize_n, "closure")
    expect_identical(graph_optimize_n, graph_optimise_n)
})

test_that("graph_optimise_n: validates inputs", {
    gc <- graph_constraint_free(3)
    ts <- mock_trial_success(disj_power_fn)

    # neither target nor local_power_target
    expect_error(
        graph_optimise_n(
            graph_constraint = gc,
            effect_size = rep(0.2, 3),
            trial_success = ts
        ),
        "target.*local_power_target"
    )

    # target without a trial_success object
    expect_error(
        graph_optimise_n(
            graph_constraint = gc,
            effect_size = rep(0.2, 3),
            target = 0.8
        ),
        "requires.*trial_success"
    )

    # trial_success dimension mismatch
    expect_error(
        graph_optimise_n(
            graph_constraint = gc,
            effect_size = rep(0.2, 3),
            trial_success = mock_trial_success(disj_power_fn, m = 4L),
            target = 0.8
        ),
        "defined for 4"
    )

    # effect_size length mismatch
    expect_error(
        graph_optimise_n(
            graph_constraint = gc,
            effect_size = c(0.1, 0.2),
            trial_success = ts,
            target = 0.8
        ),
        "effect_size"
    )

    # local_power_target length mismatch
    expect_error(
        graph_optimise_n(
            graph_constraint = gc,
            effect_size = rep(0.2, 3),
            local_power_target = c(0.8, NA)
        ),
        "local_power_target"
    )

    # floors must lie strictly inside (0, 1)
    expect_error(
        graph_optimise_n(
            graph_constraint = gc,
            effect_size = rep(0.2, 3),
            local_power_target = c(1.2, NA, NA)
        ),
        "strictly between"
    )

    # all-NA floors do not count as an active criterion
    expect_error(
        graph_optimise_n(
            graph_constraint = gc,
            effect_size = rep(0.2, 3),
            local_power_target = c(NA, NA, NA)
        ),
        "target.*local_power_target"
    )

    # non-square sigma
    expect_error(
        graph_optimise_n(
            graph_constraint = gc,
            effect_size = rep(0.2, 3),
            trial_success = ts,
            target = 0.8,
            sigma = matrix(0, nrow = 3, ncol = 2)
        ),
        "square"
    )

    # n_range must be ordered with n_min >= 2
    expect_error(
        graph_optimise_n(
            graph_constraint = gc,
            effect_size = rep(0.2, 3),
            trial_success = ts,
            target = 0.8,
            n_range = c(500L, 100L)
        ),
        "n_range"
    )

    # a marginal floor on a non-positive effect size is unattainable
    expect_error(
        graph_optimise_n(
            graph_constraint = gc,
            effect_size = c(0.2, -0.1, 0.2),
            local_power_target = c(NA, 0.8, NA)
        ),
        "non-positive"
    )
})

test_that("graph_optimise_n: warns for non-positive effects without floors", {
    gc <- graph_constraint_free(3)
    ts <- mock_trial_success(conj_power_fn)

    # the warning fires first; the run then aborts in phase 1 because a
    # conjunctive 0.99 target is infeasible at n_max = 2 with negative effect
    expect_warning(
        expect_error(
            graph_optimise_n(
                graph_constraint = gc,
                effect_size = c(0.2, -0.1, 0.2),
                trial_success = ts,
                target = 0.99,
                n_range = c(2L, 2L),
                nsim = 100L,
                seed = 1,
                verbose = FALSE
            ),
            "No feasible seed graph"
        ),
        "monotone"
    )
})


## ------ graph_optimise_n: mocked end-to-end smoke test ------ ##

test_that("graph_optimise_n: wiring works end-to-end with a mocked criterion", {
    skip_on_cran()

    gc <- graph_constraint_free(3)
    ts <- mock_trial_success(disj_power_fn)

    ctrl <- multigrain_control() |>
        control_global(
            popSize = 8L,
            generations_per_block = 3L,
            max_blocks = 2L,
            run = 6L
        )

    result <- graph_optimise_n(
        graph_constraint = gc,
        effect_size = rep(0.25, 3),
        trial_success = ts,
        target = 0.5,
        n_range = c(20L, 200L),
        nsim = 200L,
        num_threads = 1L,
        control = ctrl,
        seed = 7,
        verbose = FALSE
    )

    expect_s3_class(result, "multigrain_graph_optimal")
    expect_identical(result$N, result$n_final)
    expect_gte(result$n_final, 20L)
    expect_lte(result$n_final, 200L)
    expect_lte(result$n_final, result$n_init)
    expect_equal(sum(result$hyp_weight), 1, tolerance = 1e-12)
    expect_equal(
        unname(rowSums(result$trans_matrix)),
        rep(1, 3),
        tolerance = 1e-12
    )
    expect_identical(result$solution$opt_source, "sample_size")
    expect_true(result$solution$graph_valid[["sample_size"]])
    expect_gte(result$power$trial_success, 0.5 - 0.1)
    expect_s3_class(result$phase1$history, "data.frame")
    expect_s3_class(result$phase2$history, "data.frame")

    # reproducibility: same seed, same answer; RNG state restored
    rng_before <- .Random.seed
    result2 <- graph_optimise_n(
        graph_constraint = gc,
        effect_size = rep(0.25, 3),
        trial_success = ts,
        target = 0.5,
        n_range = c(20L, 200L),
        nsim = 200L,
        num_threads = 1L,
        control = ctrl,
        seed = 7,
        verbose = FALSE
    )
    expect_identical(rng_before, .Random.seed)
    expect_identical(result$n_final, result2$n_final)
    expect_identical(result$hyp_weight, result2$hyp_weight)
    expect_identical(result$trans_matrix, result2$trans_matrix)
})


## ------ control_prepare_n / default_control_n ------ ##

test_that("control_prepare_n: injects sample-size defaults", {
    ctrl <- control_prepare_n(multigrain_control())

    expect_identical(ctrl$global_opt$popSize, 100L)
    expect_identical(ctrl$global_opt$generations_per_block, 25L)
    expect_identical(ctrl$global_opt$max_blocks, 100L)
    expect_identical(ctrl$global_opt$run, 200L)
    expect_equal(ctrl$global_opt$cauchy_loc, 0.5)
    expect_true(ctrl$global_opt$check_stepdown)
    expect_equal(ctrl$local_opt$poptim, 0.2)
    expect_equal(ctrl$local_opt$pressel, 0.6)
    expect_identical(ctrl$local_opt$nm_maxit, 150L)
    # lambda is resolved at the point of use, not a default
    expect_null(ctrl$global_opt$lambda)
})

test_that("control_prepare_n: user values survive the merge", {
    ctrl <- multigrain_control() |>
        control_global(popSize = 50L, generations_per_block = 10L) |>
        control_local(poptim = 0.4) |>
        control_prepare_n()

    expect_identical(ctrl$global_opt$popSize, 50L)
    expect_identical(ctrl$global_opt$generations_per_block, 10L)
    expect_identical(ctrl$global_opt$max_blocks, 100L) # default preserved
    expect_equal(ctrl$local_opt$poptim, 0.4)
    expect_equal(ctrl$local_opt$pressel, 0.6) # default preserved
})

test_that("control_prepare_n: warns about and drops nsim controls", {
    ctrl <- multigrain_control() |>
        control_nsim_local(1000)

    expect_warning(
        prepared <- control_prepare_n(ctrl),
        "ignored by"
    )
    expect_null(prepared$nsim_local)
    expect_null(prepared$nsim_global)
})


test_that("step-down known answer: fixed-sequence conjunctive power", {
    skip_on_cran()

    # For a fixed-sequence test of m = 3 independent hypotheses with common
    # effect size 0.15, conjunctive power at sample size n is
    # pnorm(0.15 * sqrt(n) - qnorm(1 - 0.025))^3; solve for power = 0.8.
    pow <- function(n) pnorm(0.15 * sqrt(n) - qnorm(0.975))^3
    n_ref <- ceiling(
        uniroot(function(n) pow(n) - 0.8, c(50, 2000), tol = 1e-6)$root
    )

    graph <- fixed_seq_graph_3m()
    z_null <- withr::with_seed(
        99,
        rqmvnorm_qr(1000, mean = rep(0, 3), sigma = diag(3))
    )
    gen <- .make_pvalue_generator_n(z_null, rep(0.15, 3))

    fit <- function(x, pvals) {
        rej <- graph_shortcut(
            pvals,
            0.025,
            graph$hyp_weight,
            graph$trans_matrix
        )
        mean(rej[, 1] & rej[, 2] & rej[, 3]) - 0.8
    }

    n_final <- .execute_stepdown_n(
        n_current = n_ref + 30L,
        candidate = numeric(0),
        get_pvals = gen,
        fitness_func = fit
    )

    # QMC at nsim = 1e3: allow +-10 of the analytic boundary
    expect_lte(abs(n_final - n_ref), 10L)
})
