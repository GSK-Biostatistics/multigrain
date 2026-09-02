# Sample-size optimisation: find the minimum per-arm sample size $n$ for
# which an optimised graphical MTP meets a trial-success target and/or
# marginal power floors. Phase 1 bisects over $n$ with seed graphs; phase 2
# runs Cauchy-mutation evolutionary search blocks with greedy step-down.
# The search engine itself lives in R/cauchy_evolution.R.


#' Memoized CRN p-value generator
#'
#' Returns a function `get_pvals(n)` producing the \eqn{nsim \times m} matrix
#' of p-values at per-arm sample size \eqn{n} from the common-random-number
#' baseline draws:
#' \deqn{p(n) = 1 - \Phi(Z_{null} + \delta \sqrt{n})}
#' where \eqn{\delta} is the standardised `effect_size` vector. Because the
#' same `z_null` underpins every \eqn{n}, feasibility comparisons across
#' sample sizes share their Monte Carlo noise (variance reduction for the
#' step-down decisions).
#'
#' Results are cached per integer \eqn{n} with FIFO eviction; the cache holds
#' at most `cache_size` matrices. `attr(get_pvals, "cache_info")()` exposes
#' the cache state for testing.
#'
#' @param z_null Numeric \eqn{nsim \times m} matrix of standard MVN draws.
#' @param effect_size Numeric length-\eqn{m} vector of standardised effect
#'   sizes (the non-centrality at sample size \eqn{n} is
#'   \eqn{\delta_j \sqrt{n}}).
#' @param cache_size Integer; maximum number of cached sample sizes
#'   (default 3). Set to 0 to disable caching.
#'
#' @return A function `function(n)` returning the p-value matrix with
#'   attribute `"n"` recording the sample size used.
#' @noRd
.check_pvalue_generator_args <- function(z_null, effect_size, cache_size) {
    if (!is.matrix(z_null)) {
        cli::cli_abort("{.arg z_null} must be a matrix.")
    }
    if (!is.numeric(effect_size)) {
        cli::cli_abort("{.arg effect_size} must be numeric.")
    }
    if (ncol(z_null) != length(effect_size)) {
        cli::cli_abort(
            "Number of columns in {.arg z_null} ({ncol(z_null)}) must match \\
            the length of {.arg effect_size} ({length(effect_size)})."
        )
    }

    cache_size <- as.integer(cache_size)
    if (is.na(cache_size) || cache_size < 0L) {
        cli::cli_abort("{.arg cache_size} must be non-negative.")
    }

    cache_size
}

.make_pvalue_generator_n <- function(z_null, effect_size, cache_size = 3L) {
    cache_size <- .check_pvalue_generator_args(
        z_null,
        effect_size,
        cache_size
    )

    # FIFO insertion order lives in the cache environment under ".keys",
    # which cannot collide with the integer-string keys of cached matrices
    cache <- new.env(parent = emptyenv())
    assign(".keys", character(0), envir = cache)

    get_pvals <- function(n) {
        if (length(n) != 1L || !is.numeric(n) || !is.finite(n) || n < 1) {
            cli::cli_abort("{.arg n} must be a single positive number.")
        }

        n_int <- as.integer(n)
        key <- as.character(n_int)

        if (cache_size > 0L && exists(key, envir = cache, inherits = FALSE)) {
            return(get(key, envir = cache, inherits = FALSE))
        }

        z_shifted <- sweep(z_null, 2, effect_size * sqrt(n_int), "+")
        pvals <- 1 - stats::pnorm(z_shifted)
        attr(pvals, "n") <- n_int

        if (cache_size > 0L) {
            assign(key, pvals, envir = cache)
            keys <- c(get(".keys", envir = cache), key)
            if (length(keys) > cache_size) {
                rm(list = keys[1L], envir = cache)
                keys <- keys[-1L]
            }
            assign(".keys", keys, envir = cache)
        }

        pvals
    }

    attr(get_pvals, "cache_info") <- function() {
        keys <- get(".keys", envir = cache)
        list(
            size = length(keys),
            keys = keys,
            max_size = cache_size
        )
    }

    get_pvals
}


#' Integer bisection for the smallest feasible sample size
#'
#' Bisects over integer \eqn{n \in} `n_range` using a caller-supplied
#' feasibility predicate, assuming feasibility is monotone non-decreasing in
#' \eqn{n}. The upper bound is evaluated first and the search aborts with an
#' actionable message if even `n_range[2]` is infeasible.
#'
#' @param is_feasible_at_n Function `function(n)` returning a list with at
#'   least `feasible` (logical), `best_idx` (integer), `f_best` (numeric),
#'   and `best_theta` (the best seed's parameter row).
#' @param n_range Length-2 integer vector `c(n_min, n_max)`.
#' @param verbose Logical; print progress messages.
#' @param call Caller environment for error reporting.
#'
#' @return A list with `n` (smallest feasible sample size found),
#'   `best_idx`, `f_best`, `best_theta` (all evaluated at `n`), and
#'   `history` (data frame with columns `n`, `feasible`, `best_idx`,
#'   `f_best` for every bisection probe).
#' @noRd
.binary_search_n <- function(
    is_feasible_at_n,
    n_range,
    verbose = FALSE,
    call = rlang::caller_env()
) {
    n_min <- as.integer(n_range[1L])
    n_max <- as.integer(n_range[2L])

    top <- is_feasible_at_n(n_max)
    if (!top$feasible) {
        cli::cli_abort(
            c(
                "No feasible seed graph at {.code n = {n_max}} (best fitness \\
                {round(top$f_best, 4)}).",
                i = "Increase {.code n_range[2]}, relax {.arg target} / \\
                {.arg local_power_target}, or supply a better \\
                {.arg start_graph}."
            ),
            call = call
        )
    }
    if (verbose) {
        cli::cli_inform("Feasible at n = {n_max}; starting bisection.")
    }

    lo <- n_min - 1L # guard: feasibility unknown below n_min
    hi <- n_max # known feasible

    probe_log <- data.frame(
        n = n_max,
        feasible = TRUE,
        best_idx = top$best_idx,
        f_best = top$f_best
    )

    while (lo + 1L < hi) {
        mid <- (lo + hi) %/% 2L
        probe <- is_feasible_at_n(mid)
        probe_log <- rbind(
            probe_log,
            data.frame(
                n = mid,
                feasible = probe$feasible,
                best_idx = probe$best_idx,
                f_best = probe$f_best
            )
        )
        if (verbose) {
            cli::cli_inform(
                "n = {mid} {if (probe$feasible) 'feasible' else \\
                'infeasible'}: f_best = {round(probe$f_best, 5)} \\
                (seed {probe$best_idx})"
            )
        }
        if (probe$feasible) {
            hi <- mid
        } else {
            lo <- mid
        }
    }

    final <- is_feasible_at_n(hi)
    if (!final$feasible) {
        cli::cli_abort(
            "Bisection ended infeasible at {.code n = {hi}}; feasibility is \\
            not monotone in {.code n}. Check that {.arg effect_size} is \\
            positive for all targeted hypotheses.",
            call = call
        )
    }

    list(
        n = hi,
        best_idx = final$best_idx,
        f_best = final$f_best,
        best_theta = final$best_theta,
        history = probe_log
    )
}


#' Signed feasibility fitness for sample-size optimisation
#'
#' Builds a closure `f(x, pvals)` returning a signed feasibility score for an
#' encoded candidate graph evaluated against a p-value matrix:
#' \deqn{f = \min(\text{slacks})}
#' where the slack vector collects `trial_success - target` (when `target` is
#' active) and `local_power_j - local_power_target_j` for every non-`NA`
#' marginal floor. The score is positive iff every active constraint is
#' satisfied, zero on the feasibility boundary, and negative otherwise —
#' feasibility everywhere in the sample-size path is `f > 0`.
#'
#' Boundary checks mirror `create_obj_func()` but each non-`NA` penalty is
#' shifted by \eqn{-1}: valid graphs score in \eqn{[-1, 1]} here (slacks are
#' differences of probabilities), whereas in `graph_optimise()` valid
#' objectives are non-negative, so the unshifted penalties never collide with
#' them. Without the shift, a mildly invalid encoding (penalty \eqn{\approx}
#' -0.01) would outrank a valid but deeply infeasible graph
#' (\eqn{\approx} -0.9). The shift keeps every invalid encoding strictly
#' below every valid graph while preserving the gradient toward validity.
#'
#' @param alpha Numeric scalar; overall one-sided significance level.
#' @param hyp_constraint Numeric vector; `NA` marks free weights.
#' @param trans_constraint Numeric matrix; `NA` marks free transitions.
#' @param power_criterion Function mapping a rejection matrix to a scalar
#'   trial-success measure (`trial_success$func`); required iff `target` is
#'   supplied.
#' @param target Numeric scalar trial-success target, or `NULL`.
#' @param local_power_target Numeric length-\eqn{m} vector of marginal power
#'   floors (`NA` = no floor), or `NULL`.
#' @param num_threads Integer; threads for `graph_shortcut_parallel()`
#'   (serial `graph_shortcut()` is used when `num_threads < 2`).
#'
#' @return A function `function(x, pvals)` returning a numeric scalar.
#' @noRd
.boundary_penalty <- function(hyp_weight, trans_matrix) {
    if (anyNA(hyp_weight) || anyNA(trans_matrix)) {
        return(-1e6)
    }
    if (any(hyp_weight < 0)) {
        return(-1 + sum(hyp_weight[hyp_weight < 0]))
    }
    if (any(trans_matrix < 0)) {
        return(-1 + sum(trans_matrix[trans_matrix < 0]))
    }
    if (any(hyp_weight > 1)) {
        return(-1 - sum(hyp_weight[hyp_weight > 1]))
    }
    if (any(trans_matrix > 1)) {
        return(-1 - sum(trans_matrix[trans_matrix > 1]))
    }
    NULL
}

.create_fitness_min_n <- function(
    alpha,
    hyp_constraint,
    trans_constraint,
    power_criterion = NULL,
    target = NULL,
    local_power_target = NULL,
    num_threads = 1L
) {
    force(alpha)
    force(hyp_constraint)
    force(trans_constraint)
    force(power_criterion)
    force(target)
    force(num_threads)

    has_ts <- !is.null(target)
    marg_idx <- if (is.null(local_power_target)) {
        integer(0)
    } else {
        which(!is.na(local_power_target))
    }
    has_marg <- length(marg_idx) > 0L
    marg_floors <- if (has_marg) local_power_target[marg_idx] else numeric(0)

    # Internal programming errors: entry-point validation must catch these
    # before the closure is ever constructed
    if (!has_ts && !has_marg) {
        cli::cli_abort(
            "At least one of {.arg target} or {.arg local_power_target} \\
            must be active."
        )
    }
    if (has_ts && !is.function(power_criterion)) {
        cli::cli_abort(
            "{.arg target} is set but {.arg power_criterion} is not a \\
            function."
        )
    }

    use_parallel <- num_threads >= 2L

    function(x, pvals) {
        theta <- split_theta(x, hyp_constraint)
        hyp_weight <- recover_full_weights(theta$w_pars, hyp_constraint)
        trans_matrix <- recover_full_trans_matrix(
            theta$g_pars,
            trans_constraint
        )

        # create_obj_func() boundary penalties shifted by -1 (see header)
        penalty <- .boundary_penalty(hyp_weight, trans_matrix)
        if (!is.null(penalty)) {
            return(penalty)
        }

        hyp_weight[hyp_weight < 1e-4] <- 0
        trans_matrix[trans_matrix < 1e-5] <- 0

        rej_matrix <- if (use_parallel) {
            graph_shortcut_parallel(
                pvals = pvals,
                alpha = alpha,
                w = hyp_weight,
                G = trans_matrix,
                num_threads = num_threads,
                grain_size = -1L
            )
        } else {
            graph_shortcut(
                pvals = pvals,
                alpha = alpha,
                w = hyp_weight,
                G = trans_matrix
            )
        }

        slacks <- c(
            if (has_ts) power_criterion(rej_matrix) - target,
            if (has_marg) colMeans(rej_matrix)[marg_idx] - marg_floors
        )

        val <- min(slacks)
        if (!is.finite(val)) {
            return(-1e6)
        }
        val
    }
}


#' Test feasibility one sample size down
#'
#' Returns `TRUE` when `candidate` remains feasible (fitness > 0) at
#' `n_current - 1`.
#'
#' @param candidate Numeric parameter vector to test.
#' @param get_pvals P-value generator from `.make_pvalue_generator_n()`.
#' @param fitness_func Function `function(x, pvals)` returning the signed
#'   feasibility fitness.
#' @param n_current Current sample size.
#' @param n_min Smallest admissible sample size (the floor of `n_range`).
#'
#' @return Logical scalar.
#' @noRd
.try_step_down_n <- function(
    candidate,
    get_pvals,
    fitness_func,
    n_current,
    n_min = 2L
) {
    if (n_current <= n_min) {
        return(FALSE)
    }

    pvals_down <- get_pvals(n_current - 1L)
    fitness_func(candidate, pvals_down) > 0
}


#' Greedy step-down loop
#'
#' Repeatedly decrements the sample size while `candidate` stays feasible,
#' stopping at the first infeasible decrement or at `n_min`.
#'
#' @inheritParams .try_step_down_n
#' @param verbose Logical; print each accepted step-down.
#'
#' @return The final (smallest feasible) sample size as an integer.
#' @noRd
.execute_stepdown_n <- function(
    n_current,
    candidate,
    get_pvals,
    fitness_func,
    n_min = 2L,
    verbose = FALSE
) {
    n_new <- as.integer(n_current)
    while (n_new > n_min) {
        stepped <- .try_step_down_n(
            candidate,
            get_pvals,
            fitness_func,
            n_new,
            n_min
        )
        if (stepped) {
            n_new <- n_new - 1L
            if (verbose) {
                cli::cli_inform("Step down: {n_new + 1L} -> {n_new}")
            }
        } else {
            break
        }
    }
    n_new
}


#' Run one Cauchy-mutation search block at fixed sample size
#'
#' Freezes the p-value matrix at `n` and runs `block_config$kappa`
#' generations of `.cauchy_evolution_core()` with a Nelder-Mead callback.
#' `patience` is set to `kappa + 1` so a block never stops early — global
#' stagnation is handled across blocks by `.phase2_cauchy_stepdown()`.
#'
#' @param n Integer; the fixed sample size for this block.
#' @param x0 Numeric matrix of seeds for the block population.
#' @param get_pvals P-value generator from `.make_pvalue_generator_n()`.
#' @param fitness_func Function `function(x, pvals)` to maximise.
#' @param block_config List of knobs: `kappa`, `mu`, `lambda`, `poptim`,
#'   `pressel`, `nm_ctrl`, `cauchy_loc`, `cauchy_scale`, `cauchy_band`.
#' @param verbose Logical; print progress messages.
#' @param projector Box projector for the Nelder-Mead callback.
#'
#' @return A list with `best_x`, `best_f`, `parents`, `f_par`, `history`,
#'   `n`, `gens_run`, and `local_searches` (count of accepted Nelder-Mead
#'   improvements).
#' @noRd
.run_cauchy_mutation_block <- function(
    n,
    x0,
    get_pvals,
    fitness_func,
    block_config = list(),
    verbose = FALSE,
    projector = .proj_box01
) {
    block_gens <- block_config$kappa %||% 25L
    mu <- block_config$mu %||% max(10L, nrow(x0))
    lambda <- block_config$lambda %||% as.integer(ceiling(1.5 * mu))
    poptim <- block_config$poptim %||% 0.2
    pressel <- block_config$pressel %||% 0.6
    nm_ctrl <- block_config$nm_ctrl %||% list(maxit = 150L, reltol = 1e-8)

    pvals <- get_pvals(n)
    f_at_n <- function(x) fitness_func(x, pvals)

    # the success count lives in an environment so the callback can update
    # it without `<<-`
    nm_state <- new.env(parent = emptyenv())
    assign("successes", 0L, envir = nm_state)
    cb <- function(state) {
        nm_step <- .apply_local_search_nm(
            parents = state$parents,
            f_par = state$f_par,
            eval_f = state$eval_f,
            poptim = poptim,
            pressel = pressel,
            nm_ctrl = nm_ctrl,
            projector = projector,
            verbose = verbose
        )
        if (nm_step$improved) {
            assign(
                "successes",
                get("successes", envir = nm_state) + 1L,
                envir = nm_state
            )
            return(list(
                replace_idx = nm_step$idx,
                X = matrix(nm_step$x_new, nrow = 1)
            ))
        }
        NULL
    }

    if (verbose) {
        cli::cli_inform("  Block at n = {n}: {block_gens} generation{?s}")
    }

    # NOTE: unlike the December prototype, the Cauchy mutation parameters are
    # actually forwarded to the engine here (the prototype collected them in
    # block_config but never passed them on)
    es <- .cauchy_evolution_core(
        f = f_at_n,
        x0 = x0,
        mu = mu,
        lambda = lambda,
        max_gens = block_gens,
        patience = block_gens + 1L, # no early stop inside a block
        cauchy_loc = block_config$cauchy_loc %||% 0.5,
        cauchy_scale = block_config$cauchy_scale %||% 1.0,
        cauchy_band = block_config$cauchy_band %||% 1.0,
        callback = cb,
        verbose = FALSE
    )

    list(
        best_x = es$best_x,
        best_f = es$best_f,
        parents = es$parents,
        f_par = es$f_par,
        history = es$history,
        n = n,
        gens_run = nrow(es$history),
        local_searches = get("successes", envir = nm_state)
    )
}


#' Reseed the population after an accepted step-down
#'
#' Builds the next-block population: the global best, the top `mu - 1`
#' parents from the last block, then fresh uniform rows to fill `mu`.
#' @noRd
.reseed_after_stepdown <- function(global_best_x, block_result, mu_size, d) {
    ord <- order(block_result$f_par, decreasing = TRUE)
    keep <- max(0L, min(mu_size - 1L, length(ord)))
    kept <- if (keep > 0L) {
        block_result$parents[ord[seq_len(keep)], , drop = FALSE]
    } else {
        matrix(numeric(0), nrow = 0, ncol = d)
    }
    n_fill <- max(0L, mu_size - 1L - nrow(kept))
    fresh <- if (n_fill > 0L) {
        matrix(stats::runif(n_fill * d), ncol = d)
    } else {
        matrix(numeric(0), nrow = 0, ncol = d)
    }
    rbind(matrix(global_best_x, nrow = 1), kept, fresh)
}


#' Attempt a single opportunistic step-down after a block
#'
#' @return A list with the (possibly reduced) `n_current`, `stepped_down`
#'   (logical), `step_down` (the accepted `n`, or `NA`), and the reseeded
#'   `current_pop`.
#' @noRd
.maybe_step_down_block <- function(
    n_current,
    n_min,
    check_stepdown,
    global_best_x,
    global_best_f,
    block_result,
    get_pvals,
    fitness_func,
    block_config,
    current_pop,
    verbose = FALSE
) {
    result <- list(
        n_current = n_current,
        stepped_down = FALSE,
        step_down = NA_integer_,
        current_pop = current_pop
    )

    ok <- check_stepdown && n_current > n_min && is.finite(global_best_f)
    if (!ok) {
        return(result)
    }

    n_new <- .execute_stepdown_n(
        n_current = n_current,
        candidate = global_best_x,
        get_pvals = get_pvals,
        fitness_func = fitness_func,
        n_min = n_min,
        verbose = verbose
    )

    if (!(is.finite(n_new) && n_new < n_current)) {
        return(result)
    }

    if (verbose) {
        cli::cli_inform("  Step-down accepted: n {n_current} -> {n_new}")
    }

    result$n_current <- as.integer(n_new)
    result$stepped_down <- TRUE
    result$step_down <- as.integer(n_new)
    result$current_pop <- .reseed_after_stepdown(
        global_best_x,
        block_result,
        block_config$mu,
        ncol(current_pop)
    )
    result
}


#' Phase 2: Cauchy-mutation blocks with opportunistic step-down
#'
#' Orchestrates the main optimisation loop: runs fixed-\eqn{n} search blocks,
#' tracks the global best candidate, and after each block greedily steps the
#' sample size down while the global best stays feasible. After a step-down
#' the population is reseeded with the global best, the top `mu - 1` parents
#' from the last block, and fresh uniform rows.
#'
#' Stops when `max_blocks` is reached or after `patience` generations
#' (accumulated across blocks) without improving the global best fitness.
#'
#' @param n_init Integer; starting sample size from phase 1.
#' @param x0 Numeric matrix of initial population seeds.
#' @param get_pvals P-value generator from `.make_pvalue_generator_n()`.
#' @param fitness_func Function `function(x, pvals)` to maximise.
#' @param control_n List of knobs: `kappa` (generations per block),
#'   `patience` (global, in generations), `max_blocks`, `check_stepdown`,
#'   `mu`, `lambda`, `poptim`, `pressel`, `nm_ctrl`, `cauchy_loc`,
#'   `cauchy_scale`, `cauchy_band`.
#' @param n_min Smallest admissible sample size (floor of `n_range`).
#' @param verbose Logical; print progress messages.
#' @param projector Box projector for the Nelder-Mead callback.
#'
#' @return A list with `n_final`, `best_x`, `best_f`, `history` (block-level
#'   data frame with columns `block`, `n`, `gens`, `f_best_block`,
#'   `f_best_global`, `step_down`), `total_gens`, and `step_downs` (integer
#'   vector of accepted step-down targets).
#' @noRd
.phase2_cauchy_stepdown <- function(
    n_init,
    x0,
    get_pvals,
    fitness_func,
    control_n = list(),
    n_min = 2L,
    verbose = FALSE,
    projector = .proj_box01
) {
    stopifnot(
        is.numeric(n_init),
        length(n_init) == 1L,
        is.finite(n_init),
        is.matrix(x0),
        is.numeric(x0),
        is.function(get_pvals),
        is.function(fitness_func)
    )

    block_gens <- as.integer(control_n$kappa %||% 25L)
    patience_gen <- as.integer(control_n$patience %||% 200L)
    max_blocks <- as.integer(control_n$max_blocks %||% 100L)
    check_stepdown <- isTRUE(control_n$check_stepdown %||% TRUE)

    block_config <- list(
        kappa = block_gens,
        mu = control_n$mu %||% 100L,
        lambda = control_n$lambda,
        poptim = control_n$poptim %||% 0.2,
        pressel = control_n$pressel %||% 0.6,
        nm_ctrl = control_n$nm_ctrl %||% list(maxit = 150L, reltol = 1e-8),
        cauchy_loc = control_n$cauchy_loc %||% 0.5,
        cauchy_scale = control_n$cauchy_scale %||% 1.0,
        cauchy_band = control_n$cauchy_band %||% 1.0
    )

    n_current <- as.integer(n_init)
    current_pop <- x0

    global_best_x <- NULL
    global_best_f <- -Inf
    global_gens_since_improve <- 0L
    total_gens <- 0L
    block_id <- 0L
    step_downs <- integer(0)

    hist_df <- data.frame(
        block = integer(0),
        n = integer(0),
        gens = integer(0),
        f_best_block = numeric(0),
        f_best_global = numeric(0),
        step_down = logical(0)
    )

    if (verbose) {
        cli::cli_inform("Phase 2: start at n = {n_current}")
    }

    while (block_id < max_blocks) {
        block_id <- block_id + 1L
        if (verbose) {
            cli::cli_inform("=== Block {block_id} @ n = {n_current} ===")
        }

        block_result <- .run_cauchy_mutation_block(
            n = n_current,
            x0 = current_pop,
            get_pvals = get_pvals,
            fitness_func = fitness_func,
            block_config = block_config,
            verbose = verbose,
            projector = projector
        )

        gens_done <- as.integer(block_result$gens_run %||% block_gens)
        total_gens <- total_gens + gens_done

        improved <- is.finite(block_result$best_f) &&
            (block_result$best_f > global_best_f)
        if (improved) {
            global_best_f <- block_result$best_f
            global_best_x <- block_result$best_x
            global_gens_since_improve <- 0L
            if (verbose) {
                cli::cli_inform(
                    "  New global best f = {round(global_best_f, 6)}"
                )
            }
        } else {
            global_gens_since_improve <- global_gens_since_improve + gens_done
        }

        stepped_down <- FALSE
        step_res <- .maybe_step_down_block(
            n_current = n_current,
            n_min = n_min,
            check_stepdown = check_stepdown,
            global_best_x = global_best_x,
            global_best_f = global_best_f,
            block_result = block_result,
            get_pvals = get_pvals,
            fitness_func = fitness_func,
            block_config = block_config,
            current_pop = current_pop,
            verbose = verbose
        )
        stepped_down <- step_res$stepped_down
        if (stepped_down) {
            step_downs <- c(step_downs, step_res$step_down)
            n_current <- step_res$n_current
            current_pop <- step_res$current_pop
        }

        # n is recorded after any step-down
        hist_df <- rbind(
            hist_df,
            data.frame(
                block = block_id,
                n = n_current,
                gens = gens_done,
                f_best_block = block_result$best_f,
                f_best_global = global_best_f,
                step_down = stepped_down
            )
        )

        if (global_gens_since_improve >= patience_gen) {
            if (verbose) {
                cli::cli_inform(
                    "Stopping: global patience exhausted \\
                    ({global_gens_since_improve} generation{?s} without \\
                    improvement)."
                )
            }
            break
        }

        if (!stepped_down) {
            current_pop <- block_result$parents
        }
    }

    if (block_id >= max_blocks && verbose) {
        cli::cli_inform("Stopping: reached max_blocks = {max_blocks}.")
    }

    list(
        n_final = n_current,
        best_x = global_best_x,
        best_f = global_best_f,
        history = hist_df,
        total_gens = total_gens,
        step_downs = step_downs
    )
}


#' Validate the feasibility criteria (`target` and marginal floors)
#'
#' Errors when neither criterion is supplied, when `target` lacks a matching
#' `trial_success`, or when a marginal floor is misspecified.
#' @noRd
.check_feasibility_criteria_n <- function(
    target,
    trial_success,
    local_power_target,
    m,
    call = rlang::caller_env()
) {
    has_floors <- !is.null(local_power_target) &&
        !all(is.na(local_power_target))
    if (is.null(target) && !has_floors) {
        cli::cli_abort(
            "Provide at least one feasibility criterion: {.arg target} \\
            and/or {.arg local_power_target}.",
            call = call
        )
    }

    if (!is.null(target)) {
        rlang::check_number_decimal(target, min = 0, max = 1, call = call)
        if (!is_trial_success(trial_success)) {
            cli::cli_abort(
                "{.arg target} requires {.arg trial_success} to be a \\
                {.cls multigrain_trial_success} object; see \\
                {.fn trial_success}.",
                call = call
            )
        }
        if (trial_success$m != m) {
            cli::cli_abort(
                "{.arg trial_success} is defined for {trial_success$m} \\
                hypothes{?is/es} but {.arg graph_constraint} defines {m}.",
                call = call
            )
        }
    }

    if (!is.null(local_power_target)) {
        check_double(local_power_target, allow_na = TRUE, call = call)
        if (length(local_power_target) != m) {
            cli::cli_abort(
                "{.arg local_power_target} has length \\
                {length(local_power_target)} but {.arg graph_constraint} \\
                defines {m} hypothes{?is/es}. Use {.code NA} for \\
                hypotheses without a floor.",
                call = call
            )
        }
        floor_vals <- local_power_target[!is.na(local_power_target)]
        if (any(floor_vals <= 0 | floor_vals >= 1)) {
            cli::cli_abort(
                "Non-{.code NA} entries of {.arg local_power_target} must \\
                lie strictly between 0 and 1.",
                call = call
            )
        }
    }

    invisible(NULL)
}


#' Validate `n_range` and resolve the integer bounds
#'
#' @return A list with integer `n_min` and `n_max`.
#' @noRd
.resolve_n_range_n <- function(n_range, call = rlang::caller_env()) {
    if (length(n_range) != 2L || anyNA(n_range) || !is.numeric(n_range)) {
        cli::cli_abort(
            "{.arg n_range} must be a length-2 integer vector \\
            {.code c(n_min, n_max)}.",
            call = call
        )
    }
    n_min <- as.integer(n_range[1L])
    n_max <- as.integer(n_range[2L])
    if (n_min < 2L || n_max < n_min) {
        cli::cli_abort(
            "{.arg n_range} must satisfy {.code 2 <= n_range[1] <= \\
            n_range[2]}.",
            call = call
        )
    }
    list(n_min = n_min, n_max = n_max)
}


#' Validate the inputs of `graph_optimise_n()`
#'
#' Performs all entry-point validation and resolves derived quantities.
#' Errors are reported against the caller (`graph_optimise_n()`).
#'
#' @inheritParams graph_optimise_n
#' @param call Caller environment for error reporting.
#'
#' @return A list with `m`, `sigma` (defaulted to the identity when the
#'   input was `NULL`), `sigma_defaulted`, `n_min`, and `n_max`.
#' @noRd
.check_args_graph_optimise_n <- function(
    graph_constraint,
    effect_size,
    trial_success,
    target,
    local_power_target,
    alpha,
    sigma,
    n_range,
    nsim,
    start_graph,
    num_threads,
    control,
    seed,
    verbose,
    trace,
    call = rlang::caller_env()
) {
    check_graph_constraint(graph_constraint, call = call)
    m <- graph_constraint_get_m(graph_constraint)

    check_double(effect_size, call = call)
    if (length(effect_size) != m) {
        cli::cli_abort(
            "{.arg effect_size} has length {length(effect_size)} but \\
            {.arg graph_constraint} defines {m} hypothes{?is/es}.",
            call = call
        )
    }

    check_trial_success(trial_success, allow_null = TRUE, call = call)

    .check_feasibility_criteria_n(
        target = target,
        trial_success = trial_success,
        local_power_target = local_power_target,
        m = m,
        call = call
    )

    rlang::check_number_decimal(alpha, min = 0, max = 1, call = call)

    sigma_defaulted <- is.null(sigma)
    if (sigma_defaulted) {
        sigma <- diag(m)
    } else {
        check_double_matrix(sigma, call = call)
        if (nrow(sigma) != ncol(sigma) || nrow(sigma) != m) {
            cli::cli_abort(
                "{.arg sigma} must be a square {m} x {m} matrix.",
                call = call
            )
        }
    }

    n_bounds <- .resolve_n_range_n(n_range, call = call)
    n_min <- n_bounds$n_min
    n_max <- n_bounds$n_max

    rlang::check_number_whole(nsim, min = 1, call = call)
    rlang::check_number_whole(num_threads, min = 1, call = call)
    check_control(control, call = call)
    if (!is.null(seed)) {
        rlang::check_number_whole(seed, call = call)
    }
    check_logical(verbose, allow_na = FALSE, call = call)
    check_logical(trace, allow_na = FALSE, call = call)
    .validate_start_graphs(start_graph, m, where = "graph_optimise_n")

    # Bisection and step-down assume feasibility is monotone in n; that can
    # fail for non-positive effect sizes (a marginal floor on such a
    # hypothesis may even be unattainable)
    nonpos <- which(effect_size <= 0)
    if (length(nonpos) > 0L) {
        floored_nonpos <- !is.null(local_power_target) &&
            !all(is.na(local_power_target[nonpos]))
        if (floored_nonpos) {
            cli::cli_abort(
                "{.arg local_power_target} sets a floor for \\
                hypothes{?is/es} {nonpos} with non-positive \\
                {.arg effect_size}; increasing {.code n} cannot meet it.",
                call = call
            )
        }
        cli::cli_warn(
            "{.arg effect_size} has non-positive entr{?y/ies} (position{?s} \\
            {nonpos}); the search assumes feasibility is monotone in \\
            {.code n}."
        )
    }

    list(
        m = m,
        sigma = sigma,
        sigma_defaulted = sigma_defaulted,
        n_min = n_min,
        n_max = n_max
    )
}


#' Optimise the sample size of a graph-based multiple testing procedure
#'
#' `graph_optimise_n()` finds the smallest per-arm sample size \eqn{n} for
#' which an optimised graphical multiple testing procedure satisfies a trial
#' success target (`target`) and/or per-hypothesis marginal power floors
#' (`local_power_target`). The graph weights and transition matrix are
#' optimised jointly with the search over \eqn{n}.
#'
#' @details
#' The optimisation proceeds in two phases over a single set of
#' common-random-number (CRN) draws \eqn{Z_{null} \sim MVN(0, \Sigma)}
#' generated with randomised quasi-Monte Carlo (Sobol' sequences). P-values
#' at sample size \eqn{n} are \eqn{p(n) = 1 - \Phi(Z_{null} + \delta
#' \sqrt{n})}, where \eqn{\delta} is `effect_size`; sharing `z_null` across
#' all \eqn{n} removes Monte Carlo noise from comparisons between adjacent
#' sample sizes.
#'
#' **Phase 1** performs an integer bisection over `n_range`, evaluating a
#' bank of seed graphs (defaults plus `start_graph`) at each candidate
#' \eqn{n}; the smallest \eqn{n} with at least one feasible seed becomes the
#' phase-2 starting point. **Phase 2** runs blocks of a Cauchy-mutation
#' evolutionary search with probabilistic Nelder-Mead refinement; after each
#' block the sample size is greedily stepped down while the best graph found
#' so far remains feasible. Finally, redundant weights and edges are pruned
#' whenever removal preserves feasibility at the final sample size.
#'
#' Feasibility throughout is judged by the signed score
#' \eqn{\min(\text{slacks})}, where the slacks are `trial_success - target`
#' and `local_power - local_power_target` (non-`NA` entries only); a graph
#' is feasible when all slacks are positive.
#'
#' Bisection and step-down assume feasibility is monotone non-decreasing in
#' \eqn{n}, which requires positive `effect_size` entries for all targeted
#' hypotheses.
#'
#' Search hyperparameters are set through [multigrain_control()]:
#' [control_global()] carries the evolutionary search knobs (`popSize`,
#' `generations_per_block`, `max_blocks`, `run` (patience in generations),
#' `lambda`, `cauchy_loc`, `cauchy_scale`, `cauchy_band`, `check_stepdown`)
#' and [control_local()] the Nelder-Mead knobs (`poptim`, `pressel`,
#' `nm_maxit`, `nm_reltol`). Fields meaningful only to [graph_optimise()]
#' (e.g. `pcrossover`) are ignored on this path, as are
#' [control_nsim_local()] / [control_nsim_global()] — use `nsim` instead.
#'
#' @param graph_constraint A `multigrain_graph_constraint` object containing
#'   constraints on the graph's weights and transition matrix. Created with
#'   [graph_constraint()].
#' @param effect_size (numeric) Length-\eqn{m} vector of standardised effect
#'   sizes \eqn{\delta}: the non-centrality of test statistic \eqn{j} at
#'   per-arm sample size \eqn{n} is \eqn{\delta_j \sqrt{n}}. [calc_ncp()]
#'   can help derive values from nominal power at a reference sample size.
#' @param trial_success A `trial_success` object defining the trial success
#'   measure. Required when `target` is supplied; otherwise optional (used
#'   for reporting only). Created with [trial_success()].
#' @param target (numeric scalar) Feasibility threshold for the trial
#'   success measure (e.g. `0.8`). At least one of `target` or
#'   `local_power_target` must be supplied.
#' @param local_power_target (numeric) Length-\eqn{m} vector of marginal
#'   power floors; use `NA` for hypotheses without a floor.
#' @param alpha (numeric scalar) Overall one-sided significance level.
#'   Default is 0.025.
#' @param sigma (numeric matrix) \eqn{m \times m} correlation matrix of the
#'   test statistics. Defaults to the identity (independence) when `NULL`.
#' @param n_range (integer) Length-2 vector `c(n_min, n_max)` bounding the
#'   sample-size search. The optimisation aborts if no seed graph is
#'   feasible at `n_max`.
#' @param nsim (integer) Number of CRN simulations used for every power
#'   evaluation. Default is `50000`.
#' @param start_graph Optional list of starting graphs, each a list with
#'   elements `hyp_weight` and `trans_matrix`; these supplement the default
#'   seed bank in both phases.
#' @param num_threads (integer) Number of threads for parallel execution of
#'   the shortcut algorithm. On shared systems (HPC clusters, login nodes),
#'   always explicitly set `num_threads` based on your resource allocation.
#'   Default is `1L` (serial execution).
#' @param control Optional. A `multigrain_control` object; see Details for
#'   the fields used by the sample-size path.
#' @param seed Optional integer. When supplied, the run (QMC digital shift
#'   and evolutionary search) is reproducible and the caller's random number
#'   generator state is restored on exit. When `NULL`, the current RNG
#'   state is consumed; wrap the call in [withr::with_seed()] for
#'   reproducibility without an explicit argument.
#' @param verbose (logical) Print progress milestones. Default is `TRUE`.
#' @param trace (logical) Print detailed per-probe and per-block search
#'   output. Default is `FALSE`.
#'
#' @return An object of class `multigrain_graph_optimal` with the elements
#'   of [graph_optimise()]'s return value that apply (`hyp_weight`,
#'   `trans_matrix`, `constraints`, `trial_success`, `power`, `solution`,
#'   `control`, `start_graph`) plus sample-size-specific elements:
#'     * `N`, `n_final`: The selected (minimal feasible) sample size.
#'     * `n_init`: The phase-1 bisection result the search started from.
#'     * `fitness`: Signed feasibility of the optimised graph at `n_final`
#'       (positive means all targets met, before pruning).
#'     * `target`, `local_power_target`: The supplied feasibility floors.
#'     * `phase1`: Bisection summary (`n`, `best_idx`, `best_fitness`,
#'       `history`).
#'     * `phase2`: Search summary (`history`, `total_gens`, `step_downs`).
#'
#' @export
#' @examples
#' \dontrun{
#' # Effect sizes giving 80% / 75% / 70% nominal power at n = 500
#' delta <- calc_ncp(c(0.80, 0.75, 0.70), alpha = 0.025) / sqrt(500)
#'
#' ts <- trial_success(r1 || r2 || r3)
#'
#' result <- graph_optimise_n(
#'   graph_constraint = graph_constraint_free(3),
#'   effect_size = delta,
#'   trial_success = ts,
#'   target = 0.9,
#'   n_range = c(50L, 1000L),
#'   nsim = 20000L,
#'   num_threads = 4L,
#'   seed = 42
#' )
#' result$n_final
#' }
graph_optimise_n <- function(
    graph_constraint,
    effect_size,
    trial_success = NULL,
    target = NULL,
    local_power_target = NULL,
    alpha = 0.025,
    sigma = NULL,
    n_range = c(2L, 1000L),
    nsim = 50000L,
    start_graph = NULL,
    num_threads = 1L,
    control = multigrain_control(),
    seed = NULL,
    verbose = TRUE,
    trace = FALSE
) {
    checked <- .check_args_graph_optimise_n(
        graph_constraint = graph_constraint,
        effect_size = effect_size,
        trial_success = trial_success,
        target = target,
        local_power_target = local_power_target,
        alpha = alpha,
        sigma = sigma,
        n_range = n_range,
        nsim = nsim,
        start_graph = start_graph,
        num_threads = num_threads,
        control = control,
        seed = seed,
        verbose = verbose,
        trace = trace
    )
    m <- checked$m
    sigma <- checked$sigma
    n_min <- checked$n_min
    n_max <- checked$n_max

    control <- control_prepare_n(control)

    if (!is.null(seed)) {
        # local_preserve_seed() restores the caller's RNG state on exit;
        # withr::local_seed() alone would leave the stream advanced
        withr::local_preserve_seed()
        set.seed(seed)
    }

    if (checked$sigma_defaulted && verbose) {
        cli::cli_inform(
            "No {.arg sigma} supplied; assuming independent test statistics."
        )
    }

    # ---- Common random numbers (shared across all sample sizes) ----
    z_null <- rqmvnorm_qr(
        n = as.integer(nsim),
        mean = rep(0, m),
        sigma = sigma
    )

    ts_func <- if (!is.null(trial_success)) trial_success$func else NULL

    fitness_func <- .create_fitness_min_n(
        alpha = alpha,
        hyp_constraint = graph_constraint$hyp_constraint,
        trans_constraint = graph_constraint$trans_constraint,
        power_criterion = if (!is.null(target)) ts_func else NULL,
        target = target,
        local_power_target = local_power_target,
        num_threads = num_threads
    )

    x0 <- .build_start_matrix(graph_constraint, start_graph)
    get_pvals <- .make_pvalue_generator_n(z_null, effect_size)

    # ---- Phase 1: bisection over n ----
    if (verbose) {
        cli::cli_progress_step(
            "Phase 1: bisection of n in [{n_min}, {n_max}]"
        )
    }

    eval_seeds_at_n <- function(n) {
        pvals_n <- get_pvals(n)
        fvals <- vapply(
            seq_len(nrow(x0)),
            function(i) fitness_func(x0[i, ], pvals_n),
            numeric(1)
        )
        fvals[!is.finite(fvals)] <- -Inf
        feas_idx <- which(fvals > 0)
        best_idx <- if (length(feas_idx) > 0L) {
            feas_idx[which.max(fvals[feas_idx])]
        } else {
            which.max(fvals)
        }
        list(
            feasible = length(feas_idx) > 0L,
            best_idx = best_idx,
            f_best = fvals[best_idx],
            best_theta = x0[best_idx, , drop = FALSE]
        )
    }

    ph1 <- .binary_search_n(
        eval_seeds_at_n,
        n_range = c(n_min, n_max),
        verbose = trace
    )

    if (verbose) {
        cli::cli_inform(
            "Phase 1 complete: minimal feasible n among seeds = {ph1$n} \\
            (seed #{ph1$best_idx})."
        )
    }

    # ---- Phase 2: Cauchy-mutation search with step-down ----
    control_n <- list(
        kappa = control$global_opt$generations_per_block,
        patience = control$global_opt$run,
        max_blocks = control$global_opt$max_blocks,
        check_stepdown = control$global_opt$check_stepdown,
        mu = control$global_opt$popSize,
        lambda = control$global_opt$lambda %||%
            as.integer(ceiling(1.5 * control$global_opt$popSize)),
        poptim = control$local_opt$poptim,
        pressel = control$local_opt$pressel,
        nm_ctrl = list(
            maxit = control$local_opt$nm_maxit,
            reltol = control$local_opt$nm_reltol,
            trace = if (trace) 1L else 0L
        ),
        cauchy_loc = control$global_opt$cauchy_loc,
        cauchy_scale = control$global_opt$cauchy_scale,
        cauchy_band = control$global_opt$cauchy_band
    )

    if (verbose) {
        cli::cli_progress_step(
            "Phase 2: Cauchy-mutation search with step-down checks"
        )
    }

    ph2 <- .phase2_cauchy_stepdown(
        n_init = ph1$n,
        x0 = rbind(ph1$best_theta, x0),
        get_pvals = get_pvals,
        fitness_func = fitness_func,
        control_n = control_n,
        n_min = n_min,
        verbose = trace
    )

    if (is.null(ph2$best_x) || !is.finite(ph2$best_f)) {
        cli::cli_abort(
            "Phase 2 did not produce a valid solution; increase \\
            {.code control_global(run = , max_blocks = )} or relax \\
            {.arg target} / {.arg local_power_target}."
        )
    }

    # ---- Decode, re-verify feasibility, prune, and report ----
    pvals_final <- get_pvals(ph2$n_final)
    constrained_idx <- if (is.null(local_power_target)) {
        integer(0)
    } else {
        which(!is.na(local_power_target))
    }

    verify_feasibility <- function(graph) {
        .try_prune_n(
            pvals = pvals_final,
            hyp_weight = graph$hyp_weight,
            trans_matrix = graph$trans_matrix,
            alpha = alpha,
            ts_func = if (!is.null(target)) ts_func else NULL,
            target = target,
            constrained_idx = constrained_idx,
            power_constraint = local_power_target,
            num_threads = num_threads
        )$accepted
    }

    sol <- param_to_solution(ph2$best_x, graph_constraint, process = TRUE)
    sol <- repair_graph(sol$hyp_weight, sol$trans_matrix, graph_constraint)
    sol_feasible <- verify_feasibility(sol)

    if (!sol_feasible) {
        # epsilon-edge processing can break a tightly feasible graph; fall
        # back to the unprocessed decode the optimiser actually evaluated
        sol_raw <- param_to_solution(
            ph2$best_x,
            graph_constraint,
            process = FALSE
        )
        sol <- repair_graph(
            sol_raw$hyp_weight,
            sol_raw$trans_matrix,
            graph_constraint
        )
        sol_feasible <- verify_feasibility(sol)
    }

    clean_graph <- prune_graph_n(
        pvals = pvals_final,
        hyp_weight = sol$hyp_weight,
        trans_matrix = sol$trans_matrix,
        graph_constraint = graph_constraint,
        trial_success = if (!is.null(target)) trial_success else NULL,
        target = target,
        power_constraint = local_power_target,
        alpha = alpha,
        gamma = 1,
        num_threads = num_threads,
        verbose = verbose
    )

    if (verbose) {
        cli::cli_progress_step(
            "Evaluating trial success of pruned graph at n = {ph2$n_final}"
        )
    }

    final_power <- calc_power_pvals(
        pvals = pvals_final,
        alpha = alpha,
        hyp_weight = clean_graph$hyp_weight,
        trans_matrix = clean_graph$trans_matrix,
        custom_power = if (!is.null(trial_success)) {
            list(trial_success = trial_success)
        } else {
            NULL
        }
    )

    # Apply hypothesis names
    hyp_names <- names(graph_constraint$hyp_constraint)
    names(clean_graph$hyp_weight) <- hyp_names
    dimnames(clean_graph$trans_matrix) <- list(hyp_names, hyp_names)

    # Shared elements are built through the graph_optimal() constructor so the
    # object matches graph_optimise(); the sample-size-specific elements are
    # added afterwards as graph_optimal() does not carry them.
    out <- graph_optimal(
        hyp_weight = clean_graph$hyp_weight,
        trans_matrix = clean_graph$trans_matrix,
        constraints = graph_constraint,
        trial_success = trial_success,
        power = final_power,
        solution = list(
            opt_source = "sample_size",
            graph_valid = c(sample_size = isTRUE(clean_graph$feasible))
        ),
        control = control,
        start_graph = start_graph
    )

    out$N <- ph2$n_final
    out$n_init <- ph1$n
    out$n_final <- ph2$n_final
    out$fitness <- ph2$best_f
    out$target <- target
    out$local_power_target <- local_power_target
    out$phase1 <- list(
        n = ph1$n,
        best_idx = ph1$best_idx,
        best_fitness = ph1$f_best,
        history = ph1$history
    )
    out$phase2 <- list(
        history = ph2$history,
        total_gens = ph2$total_gens,
        step_downs = ph2$step_downs
    )
    out$opt_settings <- list(
        n_range = c(n_min, n_max),
        N_max = n_max,
        nsim = as.integer(nsim),
        alpha = alpha,
        effect_size = effect_size,
        sigma = sigma,
        num_threads = num_threads,
        seed = seed
    )

    out
}

#' @export
#' @rdname graph_optimise_n
#' @usage NULL
graph_optimize_n <- graph_optimise_n


# ==========================================================================
# Sample-size support functions
#
# These helpers back graph_optimise_n(): the sample-size control defaults and
# their injection, the quasi-random CRN generator, and the
# feasibility-preserving pruning pipeline. They are kept in this file (rather
# than in R/control_prepare.R, R/sim_pvals.R and R/post_optim_processing.R) so
# the sample-size feature stays self-contained.
# ==========================================================================


#' Create the default control object for sample-size optimisation
#'
#' `graph_optimise_n()` runs a Cauchy-mutation evolutionary search rather
#' than `GA::ga()` + COBYLA, so its control fields differ from
#' `default_control()`: `global_opt` carries the evolutionary search knobs
#' and `local_opt` the Nelder-Mead callback knobs. `lambda` (offspring count)
#' is deliberately absent — it defaults to `ceiling(1.5 * popSize)` at the
#' point of use.
#'
#' @returns a [multigrain_control] object with sample-size defaults
#'
#' @noRd
default_control_n <- function() {
    global_opt <- list(
        popSize = 100L,
        generations_per_block = 25L,
        max_blocks = 100L,
        run = 200L,
        cauchy_loc = 0.5,
        cauchy_scale = 1.0,
        cauchy_band = 1.0,
        check_stepdown = TRUE
    )

    local_opt <- list(
        poptim = 0.2,
        pressel = 0.6,
        nm_maxit = 150L,
        nm_reltol = 1e-8
    )

    new_multigrain_control(
        local_opt = local_opt,
        global_opt = global_opt
    )
}


#' Prepare the control object for sample-size optimisation
#'
#' Injects the sample-size defaults from `default_control_n()` underneath
#' any user-modified values. `nsim_local` / `nsim_global` have no meaning on
#' this path — a single common-random-number draw (sized by the `nsim`
#' argument of `graph_optimise_n()`) is shared across all sample sizes, so
#' per-stage subsampling would break the variance-reduction guarantee of the
#' step-down comparisons. Set values trigger a warning and are dropped.
#'
#' @inheritParams control_nsim_local
#' @inheritParams rlang::args_error_context
#'
#' @returns a modified [multigrain_control]
#'
#' @noRd
control_prepare_n <- function(ctrl, call = rlang::caller_env()) {
    check_control(ctrl, call = call)

    if (!is.null(ctrl$nsim_local) || !is.null(ctrl$nsim_global)) {
        cli::cli_warn(c(
            "{.fn control_nsim_local} and {.fn control_nsim_global} are \\
            ignored by {.fn graph_optimise_n}.",
            i = "Use the {.arg nsim} argument instead; a single \\
            common-random-number draw is shared across all sample sizes."
        ))
        ctrl$nsim_local <- NULL
        ctrl$nsim_global <- NULL
    }

    default_ctrl <- default_control_n()

    ctrl$global_opt <- purrr::list_modify(
        default_ctrl$global_opt,
        !!!ctrl$global_opt
    )

    ctrl$local_opt <- purrr::list_modify(
        default_ctrl$local_opt,
        !!!ctrl$local_opt
    )

    ctrl
}


#' Quasi-random multivariate normal draws (Sobol' + \eqn{\Phi^{-1}})
#'
#' Generates `n` draws from \eqn{N_p(\mu, \Sigma)} using a randomised
#' low-discrepancy Sobol' sequence mapped through the standard-normal quantile
#' function, then linearly transformed to impose the target covariance.
#' Used by `graph_optimise_n()` to generate the common-random-number (CRN)
#' baseline draws \eqn{Z_{null}} shared across all candidate sample sizes;
#' randomised QMC gives lower Monte Carlo error than `mvtnorm::rmvnorm()` at
#' the feasibility boundary for the same number of simulations.
#'
#' A pivoted Cholesky factorisation of \eqn{\Sigma} is attempted first; if
#' \eqn{\Sigma} is numerically rank-deficient the function falls back to an
#' SVD-based square root.
#'
#' @param n Integer. Number of draws to return.
#' @param mean Numeric vector of length \eqn{p}. Mean vector \eqn{\mu}.
#' @param sigma Symmetric, positive semi-definite \eqn{p \times p} covariance
#'   matrix \eqn{\Sigma}.
#' @param round_pow2 Logical; if `TRUE`, rounds `n` up to the next power of two
#'   (at least \eqn{2^{10}}) before generation and truncates back to `n` rows,
#'   mimicking embedded-net workflows. Defaults to `FALSE` (exact-`n` draws).
#' @param randomize One of `"digital.shift"` (default, recommended for unbiased
#'   RQMC estimators) or `"none"`.
#'
#' @return An \eqn{n \times p} numeric matrix. Column names are propagated
#'   from `mean` if present.
#' @noRd
rqmvnorm_qr <- function(
    n,
    mean,
    sigma,
    round_pow2 = FALSE,
    randomize = c("digital.shift", "none")
) {
    randomize <- rlang::arg_match(randomize)

    if (!is.matrix(sigma)) {
        cli::cli_abort("{.arg sigma} must be a matrix.")
    }
    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps))) {
        cli::cli_abort("{.arg sigma} must be symmetric.")
    }
    if (length(mean) != nrow(sigma)) {
        cli::cli_abort(
            "{.arg mean} (length {length(mean)}) and {.arg sigma} \\
            ({nrow(sigma)} x {ncol(sigma)}) have non-conforming sizes."
        )
    }
    p <- length(mean)

    if (round_pow2) {
        log2n <- max(ceiling(log2(n)), 10L)
        if (log2n > 30L) {
            cli::cli_abort(
                "Requested {.arg n} rounds above 2^30; refusing to generate \\
                that many points."
            )
        }
        n_eff <- 2^log2n
    } else {
        n_eff <- as.integer(n)
    }

    # Factor sigma: pivoted Cholesky, else SVD-based square root. Pivoted
    # chol warns on rank-deficient input; the SVD fallback handles that case,
    # so the warning is suppressed rather than surfaced to the user.
    use_chol <- TRUE
    chol_factor <- try(
        suppressWarnings(chol(sigma, pivot = TRUE)),
        silent = TRUE
    )
    if (inherits(chol_factor, "try-error")) {
        use_chol <- FALSE
    } else {
        rnk <- attr(chol_factor, "rank")
        if (!is.null(rnk) && rnk < p) {
            use_chol <- FALSE
        }
    }

    if (!use_chol) {
        # Eigendecomposition, not SVD: for a symmetric indefinite matrix the
        # singular values are |eigenvalues|, so an SVD-based check can never
        # detect indefiniteness.
        ev <- eigen(sigma, symmetric = TRUE)
        tol <- sqrt(.Machine$double.eps) * max(abs(ev$values))
        if (any(ev$values < -tol)) {
            cli::cli_abort(
                "{.arg sigma} has negative eigenvalues (beyond numerical \\
                tolerance)."
            )
        }
        # sqrt_sigma = diag(sqrt(lambda)) V^T, so t(sqrt_sigma) %*% sqrt_sigma
        # recovers sigma
        sqrt_sigma <- sqrt(pmax(ev$values, 0)) * t(ev$vectors)
    }

    # Sobol' points in (0, 1)^p; NULL randomize resolves to "none" in qrng
    u <- qrng::sobol(
        n = n_eff,
        d = p,
        randomize = if (randomize == "none") NULL else "digital.shift"
    )
    z <- stats::qnorm(u)

    if (use_chol) {
        piv <- attr(chol_factor, "pivot")
        x <- (z %*% chol_factor)[, order(piv), drop = FALSE]
    } else {
        x <- z %*% sqrt_sigma
    }

    x <- sweep(x, 2, mean, "+")
    if (!is.null(names(mean))) {
        colnames(x) <- names(mean)
    }

    if (n_eff != n) {
        x <- x[seq_len(n), , drop = FALSE]
    }
    x
}


#' Evaluate a candidate prune against the feasibility floors
#'
#' Evaluates the candidate graph with `graph_shortcut()` (or the parallel
#' variant) directly rather than through `calc_power_pvals()`, both to honour
#' `num_threads` and because candidates produced by `.redistribute_mass()`
#' are valid by construction, so re-validation per candidate is unnecessary.
#'
#' @param pvals Numeric matrix of CRN p-values at `n_final`.
#' @param hyp_weight,trans_matrix Candidate graph.
#' @param alpha Significance level.
#' @param ts_func Compiled trial-success function, or `NULL`.
#' @param target Trial-success floor, or `NULL`.
#' @param constrained_idx Indices of hypotheses with marginal floors.
#' @param power_constraint Numeric vector of marginal floors (`NA` = none).
#' @param num_threads Threads for `graph_shortcut_parallel()`.
#'
#' @return A list with `accepted` (logical: all active floors met),
#'   `ts_value` (`NA` when no `ts_func`), and `local_power`.
#' @noRd
.try_prune_n <- function(
    pvals,
    hyp_weight,
    trans_matrix,
    alpha,
    ts_func,
    target,
    constrained_idx,
    power_constraint,
    num_threads = 1L
) {
    rej <- if (num_threads >= 2L) {
        graph_shortcut_parallel(
            pvals = pvals,
            alpha = alpha,
            w = hyp_weight,
            G = trans_matrix,
            num_threads = num_threads,
            grain_size = -1L
        )
    } else {
        graph_shortcut(
            pvals = pvals,
            alpha = alpha,
            w = hyp_weight,
            G = trans_matrix
        )
    }

    local_power <- colMeans(rej)
    ts_value <- if (is.function(ts_func)) ts_func(rej) else NA_real_

    accepted <- (is.null(target) || ts_value >= target) &&
        !.marginal_violated(local_power, constrained_idx, power_constraint)

    list(
        accepted = accepted,
        ts_value = ts_value,
        local_power = local_power
    )
}


#' Prune hypothesis weights while preserving feasibility
#'
#' Unlike `prune_hyp_weights()`, every free weight below `gamma` is a
#' candidate (the fixed-design routine restricts candidates to
#' power-constrained hypotheses when floors are present); in the sample-size
#' context any free small weight may be pruned provided all floors survive.
#'
#' @inheritParams .try_prune_n
#' @param fixed_w Indices of weights fixed by the graph constraint.
#' @param gamma Threshold below which a weight is considered for removal.
#'
#' @return A list with `hyp_weight` and `trans_matrix` (unchanged).
#' @noRd
.prune_hyp_weights_n <- function(
    pvals,
    hyp_weight,
    trans_matrix,
    ts_func,
    fixed_w,
    alpha,
    gamma,
    target,
    constrained_idx,
    power_constraint,
    num_threads = 1L
) {
    m <- length(hyp_weight)

    for (i in m:1) {
        if (hyp_weight[i] <= 0 || hyp_weight[i] >= gamma) {
            next
        }
        if (i %in% fixed_w) {
            next
        }

        w_candidate <- .redistribute_mass(
            hyp_weight,
            drop_idx = i,
            fixed_idx = fixed_w
        )

        result <- .try_prune_n(
            pvals = pvals,
            hyp_weight = w_candidate,
            trans_matrix = trans_matrix,
            alpha = alpha,
            ts_func = ts_func,
            target = target,
            constrained_idx = constrained_idx,
            power_constraint = power_constraint,
            num_threads = num_threads
        )

        if (result$accepted) {
            hyp_weight <- w_candidate
        }
    }

    list(hyp_weight = hyp_weight, trans_matrix = trans_matrix)
}


#' Prune transition edges while preserving feasibility
#'
#' @inheritParams .prune_hyp_weights_n
#' @param fixed_edge Logical matrix; `TRUE` where `trans_constraint` is
#'   non-`NA` (the diagonal is always fixed).
#'
#' @return A list with `hyp_weight` (unchanged) and `trans_matrix`.
#' @noRd
.prune_edges_n <- function(
    pvals,
    hyp_weight,
    trans_matrix,
    ts_func,
    fixed_edge,
    alpha,
    gamma,
    target,
    constrained_idx,
    power_constraint,
    num_threads = 1L
) {
    g_best <- trans_matrix
    m <- length(hyp_weight)

    for (j in m:1) {
        for (i in 1:m) {
            if (i == j) {
                next
            }
            if (g_best[i, j] <= 0 || g_best[i, j] >= gamma) {
                next
            }
            if (fixed_edge[i, j]) {
                next
            }

            g_candidate <- g_best
            g_candidate[i, ] <- .redistribute_mass(
                g_best[i, ],
                drop_idx = j,
                fixed_idx = which(fixed_edge[i, ])
            )

            result <- .try_prune_n(
                pvals = pvals,
                hyp_weight = hyp_weight,
                trans_matrix = g_candidate,
                alpha = alpha,
                ts_func = ts_func,
                target = target,
                constrained_idx = constrained_idx,
                power_constraint = power_constraint,
                num_threads = num_threads
            )

            if (result$accepted) {
                g_best <- g_candidate
            }
        }
    }

    list(hyp_weight = hyp_weight, trans_matrix = g_best)
}


#' Prune a graph while preserving sample-size feasibility
#'
#' Sets small hypothesis weights and transition edges to zero (one at a
#' time, with mass redistributed to free recipients) whenever the pruned
#' graph remains feasible at `n_final`: trial success at least `target` (if
#' supplied) and every non-`NA` entry of `power_constraint` still met.
#' Because candidates are evaluated against the same CRN p-values used by
#' the optimiser, accepted prunes cannot break the reported feasibility.
#'
#' @param pvals Numeric matrix of CRN p-values at the final sample size.
#' @param hyp_weight,trans_matrix Graph to prune.
#' @param graph_constraint A `multigrain_graph_constraint` object.
#' @param trial_success A `multigrain_trial_success` object, or `NULL` when
#'   only marginal floors are active.
#' @param target Trial-success floor, or `NULL`.
#' @param power_constraint Numeric length-\eqn{m} vector of marginal power
#'   floors (`NA` = no floor), or `NULL`.
#' @param alpha Significance level.
#' @param gamma Threshold below which parameters are considered for removal.
#' @param num_threads Threads for `graph_shortcut_parallel()`.
#' @param verbose Logical; print a progress step.
#'
#' @return A list with `hyp_weight`, `trans_matrix`, `trial_success_value`
#'   (`NA` when no `trial_success`), `local_power`, and `feasible` (logical;
#'   whether the returned graph meets all active floors).
#' @noRd
prune_graph_n <- function(
    pvals,
    hyp_weight,
    trans_matrix,
    graph_constraint,
    trial_success = NULL,
    target = NULL,
    power_constraint = NULL,
    alpha = 0.025,
    gamma = 1,
    num_threads = 1L,
    verbose = FALSE
) {
    constrained_idx <- if (is.null(power_constraint)) {
        integer(0)
    } else {
        which(!is.na(power_constraint))
    }

    # without an active floor every prune would be vacuously "feasible"
    if (is.null(target) && length(constrained_idx) == 0L) {
        cli::cli_abort(
            "{.fn prune_graph_n} requires at least one active constraint: \\
            {.arg target} and/or a non-{.code NA} entry of \\
            {.arg power_constraint}."
        )
    }

    ts_func <- if (!is.null(trial_success)) trial_success$func else NULL
    if (!is.null(target) && !is.function(ts_func)) {
        cli::cli_abort(
            "{.arg target} is set but {.arg trial_success} does not provide \\
            a function."
        )
    }

    if (verbose) {
        cli::cli_progress_step(
            "Pruning redundant weights and edges (feasibility-preserving)"
        )
    }

    fixed_w <- which(!is.na(graph_constraint$hyp_constraint))
    fixed_edge <- !is.na(graph_constraint$trans_constraint)

    pruned_w <- .prune_hyp_weights_n(
        pvals = pvals,
        hyp_weight = hyp_weight,
        trans_matrix = trans_matrix,
        ts_func = ts_func,
        fixed_w = fixed_w,
        alpha = alpha,
        gamma = gamma,
        target = target,
        constrained_idx = constrained_idx,
        power_constraint = power_constraint,
        num_threads = num_threads
    )

    pruned <- .prune_edges_n(
        pvals = pvals,
        hyp_weight = pruned_w$hyp_weight,
        trans_matrix = pruned_w$trans_matrix,
        ts_func = ts_func,
        fixed_edge = fixed_edge,
        alpha = alpha,
        gamma = gamma,
        target = target,
        constrained_idx = constrained_idx,
        power_constraint = power_constraint,
        num_threads = num_threads
    )

    final <- .try_prune_n(
        pvals = pvals,
        hyp_weight = pruned$hyp_weight,
        trans_matrix = pruned$trans_matrix,
        alpha = alpha,
        ts_func = ts_func,
        target = target,
        constrained_idx = constrained_idx,
        power_constraint = power_constraint,
        num_threads = num_threads
    )

    list(
        hyp_weight = pruned$hyp_weight,
        trans_matrix = pruned$trans_matrix,
        trial_success_value = final$ts_value,
        local_power = final$local_power,
        feasible = final$accepted
    )
}
