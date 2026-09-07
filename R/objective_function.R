#' Lexicographic score: feasibility, then fewer edges, then higher gain
#'
#' Ranks a candidate graph against a floor on the trial success measure. A
#' feasible graph (`u >= threshold`) scores `u - edge_price * n_edges`;
#' because `edge_price`
#' exceeds the largest possible difference in `u` between two feasible graphs,
#' one fewer edge always beats any gain in `u`, and among graphs with the same
#' edge count the higher `u` wins. Every infeasible graph is pushed below every
#' feasible one but still rises with `u`, so an optimiser that strays outside
#' the feasible region is pulled back towards the boundary.
#'
#' @param u (numeric scalar) Trial success measure of the candidate.
#' @param n_edges (integer scalar) Number of non-zero free transition entries.
#' @param threshold (numeric scalar) Floor on `u`.
#' @param edge_price (numeric scalar) Per-edge price, written `D` in the
#'   design record; `(u_max - threshold) + 1`.
#' @param n_free (integer scalar) Number of free transition entries.
#'
#' @returns A numeric scalar to be maximised.
#' @noRd
.lexico <- function(u, n_edges, threshold, edge_price, n_free) {
    if (u >= threshold) {
        u - edge_price * n_edges
    } else {
        u - edge_price * (n_free + 1)
    }
}


#' Exact range of the gain function over all rejection patterns
#'
#' Evaluates the compiled trial success function on each of the `2^m` possible
#' rejection patterns. Used to price an edge in [.lexico()] without assuming
#' anything about the scale of the user's objective.
#'
#' @param trial_success A `multigrain_trial_success` object.
#'
#' @returns A numeric vector with named elements `min` and `max`.
#' @noRd
.trial_success_range <- function(trial_success) {
    m <- trial_success$m
    patterns <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), m)))
    vals <- vapply(
        seq_len(nrow(patterns)),
        function(i) trial_success$func(patterns[i, , drop = FALSE]),
        numeric(1)
    )
    c(min = min(vals), max = max(vals))
}


#' Scorer for a decoded graph on a whole p-value sample
#'
#' `choose_graph()` and the final acceptance test compare graphs on the full
#' sample, and must use the same ranking the search used. This returns a
#' function of `(u, trans_matrix)` giving that ranking: the identity on `u`
#' when no threshold is in force, so the comparison is bit-identical to the
#' trial-success comparison made before this existed, and [.lexico()]
#' otherwise. The threshold is recomputed here from the reference graph on
#' `pvals`, so it belongs to this sample rather than to a subsample.
#'
#' @param pvals (numeric matrix) The sample to score on.
#' @param alpha (numeric scalar) Overall one-sided significance level.
#' @param trial_success A `multigrain_trial_success` object.
#' @param trans_constraint (numeric matrix) Transition matrix constraints.
#' @param objective_args (list) The `gain_tolerance`, `ref_graph` and
#'   `u_range` arguments of [create_obj_func()], or an empty list.
#'
#' @returns A function with signature `function(u, trans_matrix)`.
#' @noRd
.make_lexico_scorer <- function(
    pvals,
    alpha,
    trial_success,
    trans_constraint,
    objective_args
) {
    if (length(objective_args) == 0L) {
        return(function(u, trans_matrix) u)
    }

    free_mask <- is.na(trans_constraint)
    n_free <- sum(free_mask)
    ref <- objective_args$ref_graph
    u_ref <- trial_success$func(
        graph_shortcut(
            pvals = pvals,
            alpha = alpha,
            w = ref$hyp_weight,
            G = ref$trans_matrix
        )
    )
    threshold <- (1 - objective_args$gain_tolerance) * u_ref
    edge_price <- (objective_args$u_range[["max"]] - threshold) + 1

    function(u, trans_matrix) {
        .lexico(
            u = u,
            n_edges = sum(trans_matrix[free_mask] != 0),
            threshold = threshold,
            edge_price = edge_price,
            n_free = n_free
        )
    }
}


#' Create a parallelised objective function used for optimisation
#'
#' Constructs a closure that captures all data needed for fitness evaluation.
#' The returned function has signature `function(x)` — it takes only the
#' encoded parameter vector — so that optimisers (`GA::ga`, `nloptr::nloptr`)
#' do not need to forward large objects like `pvals` via `...`.
#'
#' @details When `gain_tolerance` is supplied the closure returns the
#'   lexicographic score of [.lexico()] rather than the raw trial success
#'   measure, so that the optimisers minimise the edge count subject to a floor
#'   on trial success. The floor is computed by this closure, from `ref_graph`,
#'   on the p-values this closure captured, so the reference and the candidates
#'   are always judged on the same simulated trials.
#'
#'   With the default `gain_tolerance = NULL` every returned value is
#'   arithmetically the value returned before these arguments existed
#'   (`penalty_base` is exactly `0`), and no extra random numbers are drawn.
#'
#' @param m Number of hypotheses.
#' @param power_criterion Function accessed from [trial_success()] used to
#'   calculate trial success measure.
#' @param hyp_constraint (numeric) Vector of hypothesis weight constraints.
#' @param trans_constraint (numeric matrix) Transition matrix constraints.
#' @param alpha (numeric scalar) Overall one-sided significance level. Default
#'   is 0.025.
#' @param pvals (numeric matrix) Matrix of p-values (`nsim × m`). Each row
#'   is a simulated trial; each column corresponds to a hypothesis.
#' @param gain_tolerance (numeric scalar) Fraction of the reference graph's
#'   trial success that may be traded for a smaller edge count. `NULL` (the
#'   default) disables the edge-count objective entirely.
#' @param ref_graph (list) The reference graph, with elements `hyp_weight` and
#'   `trans_matrix`. Required when `gain_tolerance` is supplied.
#' @param u_range (numeric) Named vector with elements `min` and `max` giving
#'   the exact range of `power_criterion` over all rejection patterns, as
#'   returned by [.trial_success_range()]. Required when `gain_tolerance` is
#'   supplied.
#' @inheritParams graph_optimise
#'
#' @returns A function with signature `function(x)` that evaluates the trial
#'   success measure for a given encoded parameter vector `x`. All other
#'   inputs (`alpha`, `pvals`, constraints, thread count) are captured in the
#'   closure at creation time.
#' @noRd
create_obj_func <- function(
    m,
    power_criterion,
    hyp_constraint,
    trans_constraint,
    pvals,
    alpha = 0.025,
    num_threads = 1L,
    gain_tolerance = NULL,
    ref_graph = NULL,
    u_range = NULL
) {
    force(power_criterion)
    force(hyp_constraint)
    force(trans_constraint)
    force(alpha)
    force(pvals)
    force(num_threads)

    use_parallel <- num_threads >= 2L

    shortcut <- function(w, G) {
        if (use_parallel) {
            graph_shortcut_parallel(
                pvals = pvals,
                alpha = alpha,
                w = w,
                G = G,
                num_threads = num_threads,
                grain_size = -1L
            )
        } else {
            graph_shortcut(
                pvals = pvals,
                alpha = alpha,
                w = w,
                G = G
            )
        }
    }

    enabled <- !is.null(gain_tolerance)
    penalty_base <- 0
    free_mask <- NULL
    n_free <- NULL
    threshold <- NULL
    edge_price <- NULL

    if (enabled) {
        free_mask <- is.na(trans_constraint)
        n_free <- sum(free_mask)
        u_ref <- power_criterion(
            shortcut(ref_graph$hyp_weight, ref_graph$trans_matrix)
        )
        threshold <- (1 - gain_tolerance) * u_ref
        edge_price <- (u_range[["max"]] - threshold) + 1
        penalty_base <- u_range[["min"]] - edge_price * (n_free + 2)
    }

    function(x) {
        theta <- split_theta(x, hyp_constraint)
        hyp_weight <- recover_full_weights(theta$w_pars, hyp_constraint)
        trans_matrix <- recover_full_trans_matrix(
            theta$g_pars,
            trans_constraint
        )

        if (anyNA(hyp_weight) || anyNA(trans_matrix)) {
            return(penalty_base - 1e6)
        }
        if (any(hyp_weight < 0)) {
            return(penalty_base + sum(hyp_weight[hyp_weight < 0]))
        }
        if (any(trans_matrix < 0)) {
            return(penalty_base + sum(trans_matrix[trans_matrix < 0]))
        }
        if (any(hyp_weight > 1)) {
            return(penalty_base - sum(hyp_weight[hyp_weight > 1]))
        }
        if (any(trans_matrix > 1)) {
            return(penalty_base - sum(trans_matrix[trans_matrix > 1]))
        }

        hyp_weight[hyp_weight < 1e-4] <- 0
        trans_matrix[trans_matrix < 1e-5] <- 0

        rej_matrix <- shortcut(hyp_weight, trans_matrix)

        u <- power_criterion(rej_matrix)

        if (!enabled) {
            return(u)
        }

        .lexico(
            u = u,
            n_edges = sum(trans_matrix[free_mask] != 0),
            threshold = threshold,
            edge_price = edge_price,
            n_free = n_free
        )
    }
}

#' Split encoded parameter vector into hypothesis weights, and transition
#' matrix.
#'
#' @param theta (numeric) Stacked sample-size-optimisation vector.
#' @param hyp_constraint (numeric) Vector; NA marks free weights, otherwise
#'     fixed.
#'
#' @returns list(w_pars = numeric(), g_pars = numeric()).
#' @noRd
split_theta <- function(theta, hyp_constraint) {
    stopifnot(is.numeric(theta), length(theta) >= 1L)

    free_w <- sum(is.na(hyp_constraint))
    w_len <- max(free_w - 1L, 0L)

    if (w_len) {
        w_pars <- theta[seq_len(w_len)]
        g_pars <- theta[-seq_len(w_len)]
    } else {
        w_pars <- numeric(0)
        g_pars <- theta
    }

    list(w_pars = w_pars, g_pars = g_pars)
}


#' Recover full weights from starting values or solution
#'
#' @param x A vector containing the free hypothesis weight parameters being
#'   optimised by [GA::ga()] or [nloptr::nloptr()].
#' @param hyp_constraint A vector containing constraints on hypothesis weights.
#'
#' @returns A numeric vector with the same length as `hyp_constraint`.
#' @noRd
recover_full_weights <- function(x, hyp_constraint) {
    hyp_constraint_full <- numeric(length = length(hyp_constraint))
    free_vars_indices <- which(is.na(hyp_constraint))
    hyp_constraint_full[!is.na(hyp_constraint)] <- hyp_constraint[
        !is.na(hyp_constraint)
    ]

    hyp_constraint_full[free_vars_indices] <- c(x, NA)

    hyp_constraint_full[free_vars_indices[length(free_vars_indices)]] <-
        1 - sum(hyp_constraint_full, na.rm = TRUE)

    hyp_constraint_full
}

#' Recover full transition matrix from starting values or solution
#'
#' @param x A vector containing the free transition matrix weight parameters
#'   being optimised by [GA::ga()] or [nloptr::nloptr()].
#' @param trans_constraint A vector containing constraints on transition matrix.
#'
#' @returns A matrix with the same dimensions as `trans_constraint`.
#' @noRd
recover_full_trans_matrix <- function(x, trans_constraint) {
    stopifnot(is.matrix(trans_constraint))
    m <- nrow(trans_constraint)
    n <- ncol(trans_constraint)

    G <- matrix(0, nrow = m, ncol = n)

    # Identify free elements
    free_idx_list <- apply(
        trans_constraint,
        1,
        function(x) which(is.na(x)),
        simplify = FALSE
    )
    fixed_idx_list <- apply(
        trans_constraint,
        1,
        function(x) which(!is.na(x)),
        simplify = FALSE
    )

    # How many params do we expect (k-1 per row w/ k free)
    n_params_needed <- sum(pmax(lengths(free_idx_list) - 1L, 0L))
    if (length(x) != n_params_needed) {
        stop(
            "Length of 'x' (",
            length(x),
            ") does not match number of required G parameters (",
            n_params_needed,
            ").",
            call. = FALSE
        )
    }

    x_pos <- 1L

    for (i in seq_len(m)) {
        fixed_cols <- fixed_idx_list[[i]]
        if (length(fixed_cols)) {
            G[i, fixed_cols] <- trans_constraint[i, fixed_cols]
        }

        free_cols <- free_idx_list[[i]]
        k <- length(free_cols)

        if (k == 0L) {
            next
        } else if (k == 1L) {
            G[i, free_cols] <- 1 - sum(G[i, ])
        } else {
            take <- k - 1L
            G[i, free_cols[seq_len(take)]] <- x[x_pos:(x_pos + take - 1L)]
            x_pos <- x_pos + take

            G[i, free_cols[k]] <- 1 - sum(G[i, ])
        }
    }

    G
}
