# --------------------------------------------------------------------------
# Tests for the Cauchy-mutation evolutionary search engine
# (R/cauchy_evolution.R) used by graph_optimise_n()
# --------------------------------------------------------------------------

## ------ .rtrunc_cauchy01 ------ ##

test_that(".rtrunc_cauchy01: draws stay inside the acceptance band", {
    draws <- withr::with_seed(
        42,
        .rtrunc_cauchy01(200, loc = 0.5, scale = 1, band = 0.4)
    )

    expect_length(draws, 200)
    expect_true(all(draws >= 0.3 & draws <= 0.7))
})

test_that(".rtrunc_cauchy01: respects the [0, 1] box for off-centre bands", {
    draws <- withr::with_seed(
        42,
        .rtrunc_cauchy01(100, loc = 0.9, scale = 1, band = 0.5)
    )

    expect_true(all(draws >= 0.65 & draws <= 1))
})

test_that(".rtrunc_cauchy01: deterministic under a fixed seed", {
    d1 <- withr::with_seed(1, .rtrunc_cauchy01(50, loc = 0.5))
    d2 <- withr::with_seed(1, .rtrunc_cauchy01(50, loc = 0.5))

    expect_identical(d1, d2)
})

test_that(".rtrunc_cauchy01: aborts on an empty acceptance band", {
    expect_error(
        .rtrunc_cauchy01(10, loc = 2, scale = 1, band = 1),
        "Empty acceptance band"
    )
})


## ------ .proj_box01 ------ ##

test_that(".proj_box01 clamps to the unit box", {
    expect_identical(
        .proj_box01(c(-0.5, 0, 0.3, 1, 1.7)),
        c(0, 0, 0.3, 1, 1)
    )
})


## ------ .rank_selection_prob / .select_idx_rank ------ ##

test_that(".rank_selection_prob: hand-computed case", {
    # fitness c(3, 1, 2) with q = 0.5: ranks are 1, 3, 2, so unnormalised
    # probabilities are 0.5, 0.125, 0.25 -> normalised c(4, 1, 2) / 7
    prob <- .rank_selection_prob(c(3, 1, 2), pressel = 0.5)

    expect_equal(prob, c(4, 1, 2) / 7, tolerance = 1e-12)
})

test_that(".rank_selection_prob: sums to one and ranks best highest", {
    fitness <- c(0.2, -1.5, 3.7, 0.9)
    prob <- .rank_selection_prob(fitness, pressel = 0.6)

    expect_equal(sum(prob), 1, tolerance = 1e-12)
    expect_identical(which.max(prob), which.max(fitness))
    # geometric decay: consecutive ranks differ by a factor of (1 - q)
    ord <- order(fitness, decreasing = TRUE)
    ratios <- prob[ord][-1] / prob[ord][-length(prob)]
    expect_equal(ratios, rep(1 - 0.6, 3), tolerance = 1e-12)
})

test_that(".rank_selection_prob: ties broken by first occurrence", {
    prob <- .rank_selection_prob(c(2, 2, 1), pressel = 0.5)

    expect_equal(prob, c(0.5, 0.25, 0.125) / 0.875, tolerance = 1e-12)
})

test_that(".rank_selection_prob: handles -Inf fitness and non-finite q", {
    prob <- .rank_selection_prob(c(1, -Inf), pressel = 0.5)
    expect_equal(sum(prob), 1, tolerance = 1e-12)
    expect_gt(prob[1], prob[2])

    # non-finite pressel resets to 0.25 rather than erroring
    prob_na <- .rank_selection_prob(c(1, 2), pressel = NaN)
    expect_equal(sum(prob_na), 1, tolerance = 1e-12)
})

test_that(".rank_selection_prob: all-equal fitness gives valid distribution", {
    prob <- .rank_selection_prob(rep(5, 4), pressel = 0.5)

    expect_equal(sum(prob), 1, tolerance = 1e-12)
    expect_true(all(prob > 0))
})

test_that(".select_idx_rank returns a valid index", {
    idx <- withr::with_seed(1, .select_idx_rank(c(1, 5, 2), pressel = 0.6))

    expect_true(idx %in% 1:3)
})


## ------ .apply_local_search_nm ------ ##

test_that(".apply_local_search_nm: poptim = 0 never attempts a search", {
    parents <- matrix(runif(10), nrow = 5)
    f_par <- rep(0, 5)

    res <- .apply_local_search_nm(
        parents,
        f_par,
        eval_f = function(x) 0,
        poptim = 0
    )

    expect_false(res$improved)
    expect_null(res$idx)
})

test_that(".apply_local_search_nm: improves a convex objective", {
    eval_f <- function(x) -((x[1] - 0.7)^2 + (x[2] - 0.3)^2)

    parents <- withr::with_seed(7, matrix(runif(10), nrow = 5))
    f_par <- apply(parents, 1, eval_f)

    res <- withr::with_seed(
        7,
        .apply_local_search_nm(parents, f_par, eval_f, poptim = 1)
    )

    expect_true(res$idx %in% 1:5)
    expect_true(res$improved)
    expect_gt(res$f_new, f_par[res$idx])
    # converges near the optimum at (0.7, 0.3)
    expect_lt(abs(res$x_new[1] - 0.7), 0.05)
    expect_lt(abs(res$x_new[2] - 0.3), 0.05)
})

test_that(".apply_local_search_nm: improved result stays in the box", {
    eval_f <- function(x) sum(x) # optimum at the (1, 1) corner

    parents <- withr::with_seed(7, matrix(runif(10), nrow = 5))
    f_par <- apply(parents, 1, eval_f)

    res <- withr::with_seed(
        7,
        .apply_local_search_nm(parents, f_par, eval_f, poptim = 1)
    )

    if (res$improved) {
        expect_true(all(res$x_new >= 0 & res$x_new <= 1))
    }
    expect_true(res$idx %in% 1:5)
})


## ------ .cauchy_evolution_core ------ ##

test_that(".cauchy_evolution_core: converges on the sphere benchmark", {
    # Sphere on [-5, 5]^d affine-mapped to the unit box; optimum at x = 0.5
    d <- 4
    f <- function(x) {
        y <- -5 + 10 * x
        -sum(y^2)
    }
    x0 <- withr::with_seed(1, matrix(runif(20 * d), ncol = d))

    out <- withr::with_seed(
        123,
        .cauchy_evolution_core(
            f,
            x0,
            mu = 20,
            lambda = 30,
            max_gens = 100L,
            patience = 200L
        )
    )

    expect_named(
        out,
        c("best_x", "best_f", "parents", "f_par", "history")
    )
    expect_length(out$best_x, d)
    expect_true(all(out$best_x >= 0 & out$best_x <= 1))
    expect_gt(out$best_f, -2)
    expect_identical(dim(out$parents), c(20L, 4L))
    expect_length(out$f_par, 20L)
    expect_named(out$history, c("gen", "f_best", "f_mean", "f_median"))
})

test_that(".cauchy_evolution_core: record best is monotone with keep_elite", {
    f <- function(x) -sum((x - 0.5)^2)
    x0 <- withr::with_seed(1, matrix(runif(10 * 3), ncol = 3))

    out <- withr::with_seed(
        2,
        .cauchy_evolution_core(
            f,
            x0,
            mu = 10,
            lambda = 15,
            max_gens = 30L,
            keep_elite = TRUE
        )
    )

    expect_true(all(diff(out$history$f_best) >= 0))
    expect_identical(out$best_f, max(out$history$f_best))
})

test_that(".cauchy_evolution_core: population top-up and truncation", {
    f <- function(x) sum(x)

    # top-up: 2 seed rows but mu = 8
    out_up <- withr::with_seed(
        3,
        .cauchy_evolution_core(
            f,
            matrix(runif(2 * 3), ncol = 3),
            mu = 8,
            lambda = 10,
            max_gens = 2L
        )
    )
    expect_identical(nrow(out_up$parents), 8L)

    # truncation: 12 seed rows but mu = 5
    out_down <- withr::with_seed(
        3,
        .cauchy_evolution_core(
            f,
            matrix(runif(12 * 3), ncol = 3),
            mu = 5,
            lambda = 10,
            max_gens = 2L
        )
    )
    expect_identical(nrow(out_down$parents), 5L)
})

test_that(".cauchy_evolution_core: callback injection replaces + re-evals", {
    d <- 3
    f <- function(x) -sum((x - 0.5)^2) # optimum 0 at x = 0.5

    cb_state <- new.env(parent = emptyenv())
    assign("injected", FALSE, envir = cb_state)
    cb <- function(state) {
        if (state$gen == 1L && !get("injected", envir = cb_state)) {
            assign("injected", TRUE, envir = cb_state)
            return(list(
                replace_idx = 1L,
                X = matrix(rep(0.5, d), nrow = 1)
            ))
        }
        NULL
    }

    out <- withr::with_seed(
        4,
        .cauchy_evolution_core(
            f,
            matrix(runif(5 * d), ncol = d),
            mu = 5,
            lambda = 8,
            max_gens = 3L,
            callback = cb
        )
    )

    expect_true(get("injected", envir = cb_state))
    expect_identical(out$best_f, 0)
    expect_identical(out$best_x, rep(0.5, d))
})

test_that(".cauchy_evolution_core: callback abort stops immediately", {
    f <- function(x) sum(x)
    cb <- function(state) {
        if (state$gen == 3L) {
            return(list(abort = TRUE))
        }
        NULL
    }

    out <- withr::with_seed(
        5,
        .cauchy_evolution_core(
            f,
            matrix(runif(5 * 2), ncol = 2),
            mu = 5,
            lambda = 8,
            max_gens = 50L,
            callback = cb
        )
    )

    expect_identical(nrow(out$history), 3L)
})

test_that(".cauchy_evolution_core: patience stops a stagnant run", {
    f <- function(x) 0 # constant: never improves on the initial best

    out <- withr::with_seed(
        6,
        .cauchy_evolution_core(
            f,
            matrix(runif(5 * 2), ncol = 2),
            mu = 5,
            lambda = 8,
            max_gens = 100L,
            patience = 5L
        )
    )

    expect_identical(nrow(out$history), 5L)
})

test_that(".cauchy_evolution_core: non-finite objective values are dominated", {
    # NaN over half the domain must not propagate or error
    f <- function(x) {
        if (x[1] > 0.5) {
            return(NaN)
        }
        sum(x)
    }

    out <- withr::with_seed(
        8,
        .cauchy_evolution_core(
            f,
            matrix(runif(10 * 2), ncol = 2),
            mu = 10,
            lambda = 15,
            max_gens = 10L
        )
    )

    expect_true(is.finite(out$best_f))
    expect_lte(out$best_x[1], 0.5)
})
