# Problem-agnostic (mu + lambda) evolutionary search engine with
# band-truncated Cauchy mutation and an optional Nelder-Mead callback.
# Used by the sample-size optimisation path (`graph_optimise_n()`).


#' Band-truncated Cauchy draws on \eqn{[0, 1]}
#'
#' Draws from \eqn{Cauchy(loc, scale)} by rejection sampling until each draw
#' lands in \eqn{[loc - band/2, loc + band/2] \cap [0, 1]}. Heavy-tailed but
#' bounded, so mutated genes always remain valid unit-interval values.
#'
#' @param n Integer; number of accepted draws to return.
#' @param loc,scale Location and scale of the Cauchy distribution.
#' @param band Width of the acceptance band centred on `loc`.
#'
#' @return Numeric vector of length `n` with values in \eqn{[0, 1]}.
#' @noRd
.rtrunc_cauchy01 <- function(n, loc = 0.5, scale = 1, band = 1) {
    a <- max(0, loc - band / 2)
    b <- min(1, loc + band / 2)
    if (b < a) {
        cli::cli_abort(
            "Empty acceptance band: {.arg loc} = {loc} with {.arg band} = \\
            {band} leaves no mass in [0, 1]."
        )
    }

    out <- numeric(n)
    i <- 1L
    attempts <- 0L
    # Guard against degenerate loc/scale/band combinations where the
    # acceptance probability is vanishingly small (would hang otherwise).
    max_attempts <- 1e6L
    while (i <= n) {
        attempts <- attempts + 1L
        if (attempts > max_attempts) {
            cli::cli_abort(
                "Rejection sampling failed to accept {n} draw{?s} within \\
                {max_attempts} attempts; check {.arg cauchy_loc}, \\
                {.arg cauchy_scale}, and {.arg cauchy_band}."
            )
        }
        y <- stats::rcauchy(1L, location = loc, scale = scale)
        if (y >= a && y <= b && y >= 0 && y <= 1) {
            out[i] <- y
            i <- i + 1L
        }
    }
    out
}


#' Clamp a vector to the \eqn{[0, 1]} box
#' @noRd
.proj_box01 <- function(x) {
    pmin(1, pmax(0, x))
}


#' Rank-geometric selection probabilities
#'
#' Pure-R equivalent of `GA::ga()`-style rank selection: candidate with
#' descending-rank \eqn{r} (1 = best, ties broken by first occurrence) gets
#' probability proportional to \eqn{q (1 - q)^{r - 1}}, normalised to sum to
#' one. Falls back to a uniform distribution when all probabilities underflow
#' to zero.
#'
#' @param fitness Numeric vector of fitness values (larger is better).
#' @param pressel Selection pressure \eqn{q \in (0, 1)}; clamped to
#'   \eqn{[\epsilon, 1 - \epsilon]}, non-finite values reset to 0.25.
#'
#' @return Numeric vector of probabilities summing to one.
#' @noRd
.rank_selection_prob <- function(fitness, pressel = 0.5) {
    n <- length(fitness)
    eps <- sqrt(.Machine$double.eps)
    press <- if (is.finite(pressel)) {
        min(max(pressel, eps), 1 - eps)
    } else {
        0.25
    }

    r <- rank(-fitness, ties.method = "first")
    prob <- exp(log(press) + (r - 1) * log(1 - press))
    prob[!is.finite(prob)] <- 0

    total <- sum(prob)
    if (total > 0) {
        prob / total
    } else {
        rep(1 / n, n)
    }
}


#' Sample one index by rank-geometric selection
#' @noRd
.select_idx_rank <- function(fitness, pressel = 0.5) {
    sample.int(
        length(fitness),
        1L,
        prob = .rank_selection_prob(fitness, pressel)
    )
}


#' Nelder-Mead local search over a selected parent
#'
#' With probability `poptim`, selects a parent by rank-geometric selection and
#' refines it with `stats::optim(method = "Nelder-Mead")`, maximising `eval_f`
#' via negation. Candidate points are passed through `projector` so the search
#' stays on the \eqn{[0, 1]} box.
#'
#' @param parents Numeric \eqn{\mu \times d} matrix of current parents.
#' @param f_par Numeric vector of their fitness values (length \eqn{\mu}).
#' @param eval_f Function `function(x)` returning a scalar to maximise.
#' @param poptim Probability of attempting the local search at all.
#' @param pressel Selection pressure for choosing the starting parent.
#' @param nm_ctrl Control list passed to `stats::optim()`.
#' @param projector Function used to project candidate points into the box.
#' @param verbose Logical; print progress messages.
#'
#' @return A list with `improved` (logical), `idx` (selected parent index or
#'   `NULL` if no attempt was made), and — when `improved` — `x_new`/`f_new`.
#' @noRd
.apply_local_search_nm <- function(
    parents,
    f_par,
    eval_f,
    poptim = 0.2,
    pressel = 0.6,
    nm_ctrl = list(maxit = 150L, reltol = 1e-8, trace = 0),
    projector = .proj_box01,
    verbose = FALSE
) {
    stopifnot(is.matrix(parents), length(f_par) == nrow(parents))
    if (stats::runif(1) > poptim) {
        return(list(improved = FALSE, idx = NULL, x_new = NULL, f_new = NULL))
    }

    i <- .select_idx_rank(f_par, pressel = pressel)
    x0 <- parents[i, ]

    # Maximise fitness by minimising its negation; constraints via projection
    obj <- function(x) {
        -eval_f(projector(x))
    }

    opt <- try(
        stats::optim(
            par = x0,
            fn = obj,
            method = "Nelder-Mead",
            control = nm_ctrl
        ),
        silent = TRUE
    )

    if (inherits(opt, "try-error") || !is.list(opt)) {
        if (verbose) {
            cli::cli_inform("  Nelder-Mead failed.")
        }
        return(list(improved = FALSE, idx = i, x_new = NULL, f_new = NULL))
    }

    x_new <- projector(opt$par)
    f_new <- eval_f(x_new)
    improved <- is.finite(f_new) && (f_new > f_par[i])

    if (verbose) {
        if (improved) {
            cli::cli_inform(
                "  Nelder-Mead improved: {round(f_par[i], 6)} -> \\
                {round(f_new, 6)}"
            )
        } else {
            cli::cli_inform("  Nelder-Mead did not improve.")
        }
    }

    list(
        improved = improved,
        idx = i,
        x_new = if (improved) x_new else NULL,
        f_new = if (improved) f_new else NULL
    )
}


#' Build one generation of offspring by crossover and Cauchy mutation
#'
#' One-point crossover of random parent pairs (degenerating to cloning when
#' \eqn{d = 1}), followed by band-truncated Cauchy mutation of
#' `total_mutations` random (offspring, gene) entries and clamping to the
#' unit box.
#'
#' @param parents Numeric \eqn{\mu \times d} parent matrix.
#' @param lambda Number of offspring to produce.
#' @param total_mutations Number of (offspring, gene) entries to mutate.
#' @param cauchy_loc,cauchy_scale,cauchy_band Mutation parameters; see
#'   `.rtrunc_cauchy01()`.
#'
#' @return Numeric \eqn{\lambda \times d} offspring matrix.
#' @noRd
.make_offspring <- function(
    parents,
    lambda,
    total_mutations,
    cauchy_loc = 0.5,
    cauchy_scale = 1.0,
    cauchy_band = 1.0
) {
    mu <- nrow(parents)
    d <- ncol(parents)

    crossover_one_point <- function(p1, p2) {
        if (d == 1L) {
            return(p1)
        }
        crosspoint <- sample.int(d, 1L) - 1L # crosspoint in {0, ..., d-1}
        if (crosspoint == 0L) {
            return(p2)
        }
        c(p1[1:crosspoint], p2[(crosspoint + 1L):d])
    }

    offspr <- matrix(NA_real_, nrow = lambda, ncol = d)
    for (k in seq_len(lambda)) {
        i1 <- sample.int(mu, 1L)
        i2 <- if (mu > 1L) {
            sample(setdiff(seq_len(mu), i1), 1L)
        } else {
            i1
        }
        offspr[k, ] <- crossover_one_point(parents[i1, ], parents[i2, ])
    }

    if (total_mutations > 0L) {
        idx_off <- sample.int(lambda, total_mutations, replace = TRUE)
        idx_gene <- sample.int(d, total_mutations, replace = TRUE)
        draws <- .rtrunc_cauchy01(
            total_mutations,
            loc = cauchy_loc,
            scale = cauchy_scale,
            band = cauchy_band
        )
        for (midx in seq_len(total_mutations)) {
            offspr[idx_off[midx], idx_gene[midx]] <- draws[midx]
        }
        offspr[] <- .proj_box01(offspr)
    }

    offspr
}


#' Invoke the per-generation callback and apply its instructions
#'
#' Supports the callback contract of `.cauchy_evolution_core()`: a `NULL`
#' return means no action, `list(replace_idx, X)` injects rows (which are
#' re-evaluated), and `list(abort = TRUE)` requests an immediate stop.
#'
#' @param callback The user callback, or `NULL`.
#' @param gen Current generation number.
#' @param parents,f_par Current population and fitness values.
#' @param f The objective function.
#' @param eval_rows Row-wise evaluator used to score injected rows.
#'
#' @return A list with possibly-updated `parents` and `f_par`, and `abort`
#'   (logical).
#' @noRd
.run_evolution_callback <- function(
    callback,
    gen,
    parents,
    f_par,
    f,
    eval_rows
) {
    if (!is.function(callback)) {
        return(list(parents = parents, f_par = f_par, abort = FALSE))
    }

    cb <- callback(list(
        gen = gen,
        parents = parents,
        f_par = f_par,
        best_x = parents[which.max(f_par), ],
        best_f = max(f_par),
        eval_f = f,
        box = "unit"
    ))

    if (is.list(cb) && !is.null(cb$replace_idx) && !is.null(cb$X)) {
        idx <- cb$replace_idx
        inj <- cb$X
        if (!is.matrix(inj)) {
            inj <- as.matrix(inj)
        }
        stopifnot(length(idx) == nrow(inj))
        parents[idx, ] <- inj
        f_par[idx] <- eval_rows(inj)
    }

    list(
        parents = parents,
        f_par = f_par,
        abort = is.list(cb) && isTRUE(cb$abort)
    )
}


#' Batch objective evaluator for the evolutionary search
#'
#' Returns a `function(rows)` that evaluates each row of a matrix through
#' `f`, mapping non-scalar or non-finite returns to `-Inf` (dominated).
#' @noRd
.make_eval_rows <- function(f) {
    function(rows) {
        if (!is.matrix(rows)) {
            rows <- as.matrix(rows)
        }
        out <- numeric(nrow(rows))
        for (i in seq_len(nrow(rows))) {
            val <- f(rows[i, ])
            if (length(val) != 1L || !is.finite(val)) {
                val <- -Inf
            }
            out[i] <- val
        }
        out
    }
}


#' Resolve population sizing and defaults for the evolutionary search
#'
#' Coerces `x0`, resolves `mu` / `lambda` / `total_mutations_per_gen`
#' defaults, clamps seeds into \eqn{[0, 1]}, and tops up or truncates the
#' seed matrix to `mu` rows.
#'
#' @return A list with `x0`, `mu`, `lambda`, `d`, and
#'   `total_mutations_per_gen`.
#' @noRd
.init_evolution_pop <- function(x0, mu, lambda, total_mutations_per_gen) {
    if (!is.matrix(x0)) {
        x0 <- as.matrix(x0)
    }
    storage.mode(x0) <- "double"

    d <- ncol(x0)
    stopifnot(d >= 1L)

    if (is.null(mu) || length(mu) == 0L) {
        mu <- max(10L, nrow(x0))
    }
    if (is.null(lambda) || length(lambda) == 0L) {
        lambda <- as.integer(ceiling(1.5 * mu))
    }
    stopifnot(is.finite(mu), mu >= 1L, is.finite(lambda), lambda >= 1L)

    # Clamp seeds into [0,1] while preserving dimensions
    x0[] <- .proj_box01(x0)

    # Population sizing: top up with U(0,1) rows or truncate to mu seeds
    if (nrow(x0) < mu) {
        extra <- mu - nrow(x0)
        x0 <- rbind(x0, matrix(stats::runif(extra * d), nrow = extra))
    } else if (nrow(x0) > mu) {
        x0 <- x0[seq_len(mu), , drop = FALSE]
    }

    if (is.null(total_mutations_per_gen)) {
        total_mutations_per_gen <- max(1L, floor(lambda * d / 10))
    }

    list(
        x0 = x0,
        mu = mu,
        lambda = lambda,
        d = d,
        total_mutations_per_gen = total_mutations_per_gen
    )
}


#' \eqn{(\mu + \lambda)} elitist survivor selection
#'
#' Keeps the `mu` best of the combined parent + offspring pool; when
#' `keep_elite` is `TRUE` the best parent is guaranteed to survive.
#'
#' @return A list with the surviving `parents` and their fitness `f_par`.
#' @noRd
.select_survivors <- function(parents, offspr, f_par, f_off, mu, keep_elite) {
    combined <- rbind(parents, offspr)
    f_comb <- c(f_par, f_off)
    ord <- order(f_comb, decreasing = TRUE)
    keep <- ord[seq_len(mu)]

    if (keep_elite) {
        # Parents occupy rows 1..mu of `combined`, so the best parent's
        # global index equals its index within f_par
        elite_glob <- which.max(f_par)
        if (!any(keep == elite_glob)) {
            keep[which.min(f_comb[keep])] <- elite_glob
        }
    }

    list(
        parents = combined[keep, , drop = FALSE],
        f_par = f_comb[keep]
    )
}


#' Cauchy-mutation evolutionary search on the unit box
#'
#' \eqn{(\mu + \lambda)} evolutionary search that maximises a deterministic
#' black-box objective \eqn{f(x)} over \eqn{[0, 1]^d}, using one-point
#' crossover, band-truncated Cauchy mutation, elitist selection, early
#' stopping on stagnation, and an optional per-generation `callback` that can
#' inject replacements (e.g. Nelder-Mead refinements) or abort the run.
#'
#' Non-finite or non-scalar returns from `f` are treated as `-Inf`
#' (dominated). RNG state is the caller's responsibility; wrap calls in
#' `withr::with_seed()` for reproducibility.
#'
#' @param f Function `function(x)` returning a `numeric(1)` to maximise.
#' @param x0 Numeric matrix of initial seeds (rows = candidates, columns =
#'   dimensions); values are clamped into \eqn{[0, 1]}. Topped up with
#'   \eqn{U(0, 1)} rows (or truncated) to reach `mu` rows.
#' @param mu,lambda Integers; number of parents and offspring. Defaults are
#'   `mu = max(10, nrow(x0))` and `lambda = ceiling(1.5 * mu)`.
#' @param max_gens Integer; maximum number of generations.
#' @param patience Integer; early stop after this many consecutive
#'   generations without improving the record best.
#' @param total_mutations_per_gen `NULL` or integer; number of (offspring,
#'   gene) entries mutated per generation. Defaults to
#'   `max(1, floor(lambda * d / 10))` (about 10% of all genes).
#' @param cauchy_loc,cauchy_scale,cauchy_band Mutation distribution
#'   parameters; see `.rtrunc_cauchy01()`.
#' @param keep_elite Logical; if `TRUE` the best parent is guaranteed to
#'   survive each generation, making the record best monotone.
#' @param verbose Logical; print light progress messages.
#' @param callback Optional function called once per generation with a list
#'   `(gen, parents, f_par, best_x, best_f, eval_f, box)`. It may return
#'   `NULL` (no action), `list(replace_idx, X)` to inject rows (which are
#'   re-evaluated), or `list(abort = TRUE)` to stop immediately.
#'
#' @return A list with `best_x`, `best_f`, `parents` (\eqn{\mu \times d}),
#'   `f_par` (length \eqn{\mu}), and `history` (data frame with columns
#'   `gen`, `f_best`, `f_mean`, `f_median`).
#' @noRd
.cauchy_evolution_core <- function(
    f,
    x0,
    mu = NULL,
    lambda = NULL,
    max_gens = 200L,
    patience = 100L,
    total_mutations_per_gen = NULL,
    cauchy_loc = 0.5,
    cauchy_scale = 1.0,
    cauchy_band = 1.0,
    keep_elite = TRUE,
    verbose = FALSE,
    callback = NULL
) {
    stopifnot(is.function(f))
    init <- .init_evolution_pop(x0, mu, lambda, total_mutations_per_gen)
    x0 <- init$x0
    mu <- init$mu
    lambda <- init$lambda
    total_mutations_per_gen <- init$total_mutations_per_gen

    # Batch evaluation; non-scalar / non-finite values become -Inf (dominated)
    eval_rows <- .make_eval_rows(f)

    parents <- x0
    f_par <- eval_rows(parents)

    best_idx <- which.max(f_par)
    best_x <- parents[best_idx, ]
    best_f <- f_par[best_idx]

    if (verbose) {
        cli::cli_inform(
            "Cauchy evolution init: mu = {mu}, lambda = {lambda}, \\
            d = {init$d}, \\
            best = {round(best_f, 6)}"
        )
    }

    hist_df <- data.frame(
        gen = integer(0),
        f_best = numeric(0),
        f_mean = numeric(0),
        f_median = numeric(0)
    )

    gens_since_improve <- 0L

    for (gen in seq_len(max_gens)) {
        # (1)-(2) Crossover and mutation
        offspr <- .make_offspring(
            parents,
            lambda = lambda,
            total_mutations = total_mutations_per_gen,
            cauchy_loc = cauchy_loc,
            cauchy_scale = cauchy_scale,
            cauchy_band = cauchy_band
        )

        # (3) Evaluate offspring
        f_off <- eval_rows(offspr)

        # (4) (mu + lambda) elitist selection
        sel <- .select_survivors(
            parents,
            offspr,
            f_par,
            f_off,
            mu,
            keep_elite
        )
        parents <- sel$parents
        f_par <- sel$f_par

        # (5) Optional callback for local search / injection / abort
        cb_state <- .run_evolution_callback(
            callback,
            gen = gen,
            parents = parents,
            f_par = f_par,
            f = f,
            eval_rows = eval_rows
        )
        parents <- cb_state$parents
        f_par <- cb_state$f_par
        if (cb_state$abort) {
            hist_df <- rbind(
                hist_df,
                data.frame(
                    gen = gen,
                    f_best = max(f_par),
                    f_mean = mean(f_par),
                    f_median = stats::median(f_par)
                )
            )
            break
        }

        # (6) Update record best and patience counter
        gen_best <- max(f_par)
        gen_best_x <- parents[which.max(f_par), ]
        if (gen_best > best_f) {
            best_f <- gen_best
            best_x <- gen_best_x
            gens_since_improve <- 0L
        } else {
            gens_since_improve <- gens_since_improve + 1L
        }

        # (7) Append a history row
        hist_df <- rbind(
            hist_df,
            data.frame(
                gen = gen,
                f_best = gen_best,
                f_mean = mean(f_par),
                f_median = stats::median(f_par)
            )
        )

        if (verbose && (gen %% 10L == 0L || gen == 1L)) {
            cli::cli_inform(
                "gen = {gen}, best = {round(gen_best, 6)}, mean = \\
                {round(mean(f_par), 6)}"
            )
        }

        # (8) Early stopping on stagnation
        if (gens_since_improve >= patience) {
            if (verbose) {
                cli::cli_inform(
                    "Early stop after {patience} generation{?s} without \\
                    improvement."
                )
            }
            break
        }
    }

    list(
        best_x = best_x,
        best_f = best_f,
        parents = parents,
        f_par = f_par,
        history = hist_df
    )
}
