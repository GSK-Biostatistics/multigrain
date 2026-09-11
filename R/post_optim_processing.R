# (Optimised) parameters to full solution (option to process to select epsilon
# edges)
param_to_solution <- function(
    optimised_params,
    graph_constraint,
    process = FALSE
) {
    hyp_constraint <- graph_constraint$hyp_constraint
    trans_constraint <- graph_constraint$trans_constraint
    tolerance <- graph_constraint_get_tolerance(graph_constraint)

    theta <- split_theta(optimised_params, hyp_constraint)
    w_sol <- recover_full_weights(theta$w_pars, hyp_constraint)
    G_sol <- recover_full_trans_matrix(theta$g_pars, trans_constraint)

    if (process) {
        # Convert <1e-4 w & G parameters to zero and epsilon edges
        w_sol[w_sol < 1e-4] <- 0
        w_sol[w_sol > 1e-4 & w_sol < 1e-3] <- 0.001
        fixed_w <- which(!is.na(hyp_constraint))
        w_sol <- normalise_sum(
            w_sol,
            fixed_idx = fixed_w,
            tolerance = tolerance
        )

        G_sol[G_sol < 1e-5] <- 0
        G_sol[G_sol > 1e-5 & G_sol < 1e-3] <- 0.001
        m <- nrow(G_sol)
        for (row_i in seq_len(m)) {
            fixed_in_row <- c(row_i, which(!is.na(trans_constraint[row_i, ])))
            G_sol[row_i, ] <- normalise_sum(
                G_sol[row_i, ],
                fixed_idx = fixed_in_row,
                tolerance = tolerance
            )
        }
    }

    output <- list(
        hyp_weight = w_sol,
        trans_matrix = G_sol
    )

    output
}


#' Repair a graph that may have small constraint violations
#'
#' Projects `hyp_weight` and `trans_matrix` onto the feasible region.
#' Called immediately after `param_to_solution()` inside
#' `.graph_optimise_local()` to fix tiny boundary violations caused by
#' COBYLA's tolerance-based termination.
#'
#' @param hyp_weight  (numeric) Hypothesis weight vector (length $m$).
#' @param trans_matrix (numeric) $m \times m$ transition matrix.
#' @param graph_constraint A `multigrain_graph_constraint` object.
#'
#' @returns A list with `hyp_weight` and `trans_matrix`, both guaranteed to
#'   satisfy:
#'   - All elements in \eqn{[0, 1]}
#'   - `sum(hyp_weight) == 1` (within machine precision)
#'   - Each row of `trans_matrix` sums to 1 (within machine precision)
#'   - Diagonal of `trans_matrix` is 0
#'   - Fixed elements from `multigrain_graph_constraint` are respected
#'
#' @noRd
repair_graph <- function(hyp_weight, trans_matrix, graph_constraint) {
    hc <- graph_constraint$hyp_constraint
    tc <- graph_constraint$trans_constraint
    m <- graph_constraint_get_m(graph_constraint)
    tolerance <- graph_constraint_get_tolerance(graph_constraint)

    # 1. Clamp all values to [0, 1]
    hyp_weight <- pmin(pmax(hyp_weight, 0), 1)
    trans_matrix <- pmin(pmax(trans_matrix, 0), 1)

    # 2. Force diagonal to 0
    diag(trans_matrix) <- 0

    # 3. Pin fixed elements from graph_constraint
    fixed_w <- which(!is.na(hc))
    if (length(fixed_w) > 0L) {
        hyp_weight[fixed_w] <- hc[fixed_w]
    }
    fixed_g <- which(!is.na(tc))
    if (length(fixed_g) > 0L) {
        trans_matrix[fixed_g] <- tc[fixed_g]
    }

    # 4. Normalise hyp_weight to sum to 1
    free_w <- which(is.na(hc))
    if (length(free_w) > 0L && sum(hyp_weight[free_w]) == 0) {
        # All free weights are zero — distribute remaining mass uniformly
        remaining <- 1 - sum(hyp_weight[fixed_w])
        if (remaining > 0) {
            hyp_weight[free_w] <- remaining / length(free_w)
        }
    }
    hyp_weight <- normalise_sum(
        hyp_weight,
        fixed_idx = fixed_w,
        tolerance = tolerance
    )

    # 5. Normalise each row of trans_matrix to sum to 1
    for (i in seq_len(m)) {
        fixed_in_row <- c(i, which(!is.na(tc[i, ])))
        free_cols <- setdiff(seq_len(m), fixed_in_row)

        if (length(free_cols) > 0L && sum(trans_matrix[i, free_cols]) == 0) {
            remaining <- 1 - sum(trans_matrix[i, fixed_in_row])
            if (remaining > 0) {
                trans_matrix[i, free_cols] <- remaining / length(free_cols)
            }
        }
        trans_matrix[i, ] <- normalise_sum(
            trans_matrix[i, ],
            fixed_idx = fixed_in_row,
            tolerance = tolerance
        )
    }

    list(hyp_weight = hyp_weight, trans_matrix = trans_matrix)
}


#' Find hypothesis weight indices or transition matrix row
#' indices eligible to receive redistributed mass
#'
#' @param drop_idx Integer; the index being zeroed out.
#' @param fixed_idx Integer vector; indices that are fixed by constraints.
#' @param m Integer; total number of elements.
#'
#' @returns Integer vector of recipient indices (may be length 0).
#' @noRd
.free_recipients <- function(drop_idx, fixed_idx, m) {
    setdiff(seq_len(m), c(fixed_idx, drop_idx))
}


#' Redistribute mass from a dropped element to free recipients
#'
#' Zeros out `vec[drop_idx]` and redistributes its mass proportionally among
#' the free recipients. If all recipients are zero, mass is split uniformly.
#'
#' @param vec Numeric vector (e.g. hypothesis weights or a transition row).
#' @param drop_idx Integer scalar; position to zero out.
#' @param fixed_idx Integer vector; positions that must not change.
#' @param tolerance numeric >= 0. Sum-to-one tolerance forwarded to
#'   [normalise_sum()]. Defaults to `sqrt(.Machine$double.eps)`.
#'
#' @returns Numeric vector of same length, with `vec[drop_idx] == 0`.
#' @noRd
.redistribute_mass <- function(
    vec,
    drop_idx,
    fixed_idx,
    tolerance = sqrt(.Machine$double.eps)
) {
    m <- length(vec)
    freed <- vec[drop_idx]
    vec[drop_idx] <- 0

    recipients <- .free_recipients(drop_idx, fixed_idx, m)
    if (length(recipients) == 0L) {
        return(vec)
    }

    s <- sum(vec[recipients])
    if (s > 0) {
        vec[recipients] <- vec[recipients] +
            vec[recipients] * (freed / s)
    } else {
        vec[recipients] <- freed / length(recipients)
    }

    normalise_sum(vec, fixed_idx = fixed_idx, tolerance = tolerance)
}


#' Check whether marginal power constraints are violated
#'
#' @param marg_power Numeric vector of marginal powers (length m).
#' @param constrained_idx Integer vector of indices with constraints.
#' @param power_constraint Numeric vector of required marginal powers.
#'
#' @returns Logical scalar; `TRUE` if any constraint is violated.
#' @noRd
.marginal_violated <- function(marg_power, constrained_idx, power_constraint) {
    if (length(constrained_idx) == 0L) {
        return(FALSE)
    }
    any(marg_power[constrained_idx] < power_constraint[constrained_idx])
}


#' Evaluate a candidate graph and accept if it improves on the current best
#'
#' Computes trial-success power for the candidate graph. Returns the candidate
#' values if power does not decrease and marginal constraints are satisfied;
#' otherwise returns the current best unchanged.
#'
#' @param pvals Numeric matrix of p-values (n.sim × m).
#' @param hyp_weight Candidate hypothesis weight vector.
#' @param trans_matrix Candidate transition matrix.
#' @param alpha Significance level. Default is 0.025.
#' @param trial_success A `multigrain_trial_success` object.
#' @param power_best Current best trial-success power (scalar).
#' @param constrained_idx Integer vector; which hypotheses carry marginal
#'   constraints.
#' @param power_constraint Numeric vector of marginal thresholds.
#'
#' @returns A list with `hyp_weight`, `trans_matrix`, `power_best`, and
#'   `accepted` (logical).
#' @noRd
.try_prune <- function(
    pvals,
    hyp_weight,
    trans_matrix,
    trial_success,
    power_best,
    constrained_idx,
    power_constraint,
    alpha = 0.025
) {
    power_all <- calc_power_pvals(
        pvals,
        hyp_weight = hyp_weight,
        trans_matrix = trans_matrix,
        alpha = alpha,
        custom_power = trial_success
    )

    violated <- .marginal_violated(
        power_all$local_power,
        constrained_idx,
        power_constraint
    )

    output <- list(
        hyp_weight = hyp_weight,
        trans_matrix = trans_matrix,
        power_best = power_best,
        accepted = FALSE
    )

    if (!violated && power_all$custom_power >= power_best) {
        output$power_best <- power_all$custom_power
        output$accepted <- TRUE
    }

    output
}


#' Prune hypothesis weights below a threshold if doing so does not decrease
#' estimated power or violate marginal constraints.
#'
#' Each weight \eqn{w_i < \gamma} is tentatively set to zero (with the
#' remaining weights renormalised), and the removal is accepted if the
#' estimated objective (gain) function does not decrease.
#'
#' @param pvals (numeric) Numeric matrix of p-values (n.sim x m).
#' @param hyp_weight (numeric) Vector of hypothesis weights (length m).
#' @param trans_matrix (numeric) m × m transition matrix (held fixed during
#'   weight pruning).
#' @param trial_success A `multigrain_trial_success` object defining the custom
#'   power.
#' @param fixed_w Integer vector; indices of hypothesis weights that are
#'   fixed by the graph constraint (i.e. `which(!is.na(hyp_constraint))`)
#' @param alpha (numeric) Overall one-sided significance level.
#'   Default is 0.025.
#' @param gamma Threshold below which a weight should be considered for removal.
#' @param power_constraint Optional numeric vector of length m; only entries
#'   that are non-`NA` impose a marginal-power requirement. Defaults to `NULL`
#'   (no marginal constraints).
#' @param tolerance numeric >= 0. Sum-to-one tolerance forwarded to
#'   [normalise_sum()] during mass redistribution. Defaults to
#'   `sqrt(.Machine$double.eps)`.
#'
#' @returns A list with elements `hyp_weight` (pruned weights), `trans_matrix`
#'   (unchanged), and `power_best` (the best power achieved).
#'
#' @noRd
prune_hyp_weights <- function(
    pvals,
    hyp_weight,
    trans_matrix,
    trial_success,
    fixed_w,
    alpha = 0.025,
    gamma = 1,
    power_constraint = NULL,
    tolerance = sqrt(.Machine$double.eps)
) {
    power_all <- calc_power_pvals(
        pvals,
        hyp_weight = hyp_weight,
        trans_matrix = trans_matrix,
        alpha = alpha,
        custom_power = trial_success
    )
    power_best <- power_all$custom_power

    m <- length(hyp_weight)

    constrained_idx <- if (is.null(power_constraint)) {
        integer(0)
    } else {
        which(!is.na(power_constraint))
    }

    for (i in m:1) {
        if (hyp_weight[i] <= 0 || hyp_weight[i] >= gamma) {
            next
        }
        if (i %in% fixed_w) {
            next
        }
        if (length(constrained_idx) > 0L && !(i %in% constrained_idx)) {
            next
        }

        w_candidate <- .redistribute_mass(
            hyp_weight,
            drop_idx = i,
            fixed_idx = fixed_w,
            tolerance = tolerance
        )

        result <- .try_prune(
            pvals = pvals,
            hyp_weight = w_candidate,
            trans_matrix = trans_matrix,
            alpha = alpha,
            trial_success = trial_success,
            power_best = power_best,
            constrained_idx = constrained_idx,
            power_constraint = power_constraint
        )

        if (result$accepted) {
            hyp_weight <- result$hyp_weight
            power_best <- result$power_best
        }
    }

    list(
        hyp_weight = hyp_weight,
        trans_matrix = trans_matrix,
        power_best = power_best
    )
}


#' Prune transition-matrix edges below a threshold if doing so does not
#' decrease estimated power or violate marginal constraints.
#'
#' Each edge \eqn{g_{ij} < \gamma} is tentatively set to zero (with the
#' remaining outgoing edges renormalised), and the removal is accepted if the
#' estimated objective (gain) function does not decrease.
#'
#' @param pvals (numeric) Numeric matrix of p-values (n.sim x m).
#' @param hyp_weight (numeric) Vector of hypothesis weights (held fixed during
#'   edge pruning).
#' @param trans_matrix (numeric) m × m transition matrix.
#' @param trial_success A `multigrain_trial_success` object defining the custom
#'   power.
#' @param graph_constraint A `multigrain_graph_constraint` object.
#' @param power_best (numeric) Current best trial-success power.
#' @param alpha (numeric) Overall one-sided significance level.
#'   Default is 0.025.
#' @param gamma Threshold below which an edge should be considered for removal.
#' @param power_constraint Optional numeric vector of length m; only entries
#'   that are non-`NA` impose a marginal-power requirement. Defaults to `NULL`
#'   (no marginal constraints).
#' @param tolerance numeric >= 0. Sum-to-one tolerance forwarded to
#'   [normalise_sum()] during mass redistribution. Defaults to
#'   `sqrt(.Machine$double.eps)`.
#'
#' @returns A list with elements `hyp_weight` (unchanged), `trans_matrix`
#'   (pruned edges), and `power_best` (the best power achieved).
#'
#' @noRd
prune_edges <- function(
    pvals,
    hyp_weight,
    trans_matrix,
    trial_success,
    fixed_edge,
    power_best,
    alpha = 0.025,
    gamma = 1,
    power_constraint = NULL,
    tolerance = sqrt(.Machine$double.eps)
) {
    G_best <- trans_matrix

    m <- length(hyp_weight)

    constrained_idx <- if (is.null(power_constraint)) {
        integer(0)
    } else {
        which(!is.na(power_constraint))
    }

    for (j in m:1) {
        for (i in 1:m) {
            if (i == j) {
                next
            }
            if (G_best[i, j] <= 0 || G_best[i, j] >= gamma) {
                next
            }
            if (fixed_edge[i, j]) {
                next
            }

            # fixed_idx for row i: diagonal is always fixed in trans_constraint,
            # so fixed_edge[i, i] == TRUE and the diagonal is protected.
            fixed_in_row <- which(fixed_edge[i, ])

            G_candidate <- G_best
            G_candidate[i, ] <- .redistribute_mass(
                G_best[i, ],
                drop_idx = j,
                fixed_idx = fixed_in_row,
                tolerance = tolerance
            )

            result <- .try_prune(
                pvals = pvals,
                hyp_weight = hyp_weight,
                trans_matrix = G_candidate,
                alpha = alpha,
                trial_success = trial_success,
                power_best = power_best,
                constrained_idx = constrained_idx,
                power_constraint = power_constraint
            )

            if (result$accepted) {
                G_best <- result$trans_matrix
                power_best <- result$power_best
            }
        }
    }

    list(
        hyp_weight = hyp_weight,
        trans_matrix = G_best,
        power_best = power_best
    )
}


#' Remove edges best-first, subject to a floor on the trial success measure
#'
#' Where `prune_edges()` sweeps the matrix in a fixed index order and accepts
#' any removal that does not lower the trial success measure, this evaluates
#' every remaining removable edge on the full sample, keeps the candidates that
#' stay at or above `threshold` and reduce the edge count, removes the one with
#' the highest trial success, and repeats. Spending the budget on the cheapest
#' edge first is what stops a single expensive removal from consuming it.
#'
#' Candidates that would leave a row unable to sum to one -- because the
#' dropped entry has no free recipient -- are skipped, as are candidates that
#' do not actually reduce the edge count (the uniform fallback in
#' `.redistribute_mass()` can raise it).
#'
#' @param pvals (numeric) Numeric matrix of p-values (n.sim x m).
#' @param hyp_weight (numeric) Vector of hypothesis weights (length m); not
#'   modified here.
#' @param trans_matrix (numeric) m x m transition matrix to prune.
#' @param trial_success A `multigrain_trial_success` object.
#' @param fixed_edge Logical matrix (m x m); `TRUE` where `trans_constraint`
#'   is non-`NA`.
#' @param threshold (numeric) Floor on the trial success measure.
#' @param alpha (numeric) Overall one-sided significance level.
#' @param tolerance Tolerance passed to `.redistribute_mass()` and used for the
#'   row-sum check.
#'
#' @returns A list with `hyp_weight`, `trans_matrix`, `power_best`,
#'   `prune_loss` (the exact loss across the accepted removals, negative when
#'   removing edges raised the measure) and `n_removed`.
#'
#' @noRd
.prune_edges_best_first <- function(
    pvals,
    hyp_weight,
    trans_matrix,
    trial_success,
    fixed_edge,
    threshold,
    alpha = 0.025,
    tolerance = sqrt(.Machine$double.eps)
) {
    u_of <- function(G) {
        calc_power_pvals(
            pvals,
            hyp_weight = hyp_weight,
            trans_matrix = G,
            alpha = alpha,
            custom_power = trial_success
        )$custom_power
    }

    G_best <- trans_matrix
    u_best <- u_of(G_best)
    loss <- 0
    n_removed <- 0L

    repeat {
        n_edges <- sum(G_best != 0)
        cand <- which(G_best != 0 & !fixed_edge, arr.ind = TRUE)
        best <- NULL

        for (k in seq_len(nrow(cand))) {
            i <- cand[k, 1L]
            j <- cand[k, 2L]
            G_try <- G_best
            G_try[i, ] <- .redistribute_mass(
                G_best[i, ],
                drop_idx = j,
                fixed_idx = which(fixed_edge[i, ]),
                tolerance = tolerance
            )
            if (sum(G_try != 0) >= n_edges) {
                next
            }
            if (abs(sum(G_try[i, ]) - 1) > tolerance) {
                next
            }
            u_try <- u_of(G_try)
            if (u_try >= threshold && (is.null(best) || u_try > best$u)) {
                best <- list(G = G_try, u = u_try)
            }
        }

        if (is.null(best)) {
            break
        }

        loss <- loss + (u_best - best$u)
        u_best <- best$u
        G_best <- best$G
        n_removed <- n_removed + 1L
    }

    list(
        hyp_weight = hyp_weight,
        trans_matrix = G_best,
        power_best = u_best,
        prune_loss = loss,
        n_removed = n_removed
    )
}


#' Prune graph hypothesis and transition weights according to p-value
#' distribution and trial success measure.
#'
#' @param pvals (numeric) Numeric matrix of p-values (n.sim x m).
#' @param hyp_weight (numeric) Vector of initial hypothesis weights (length m).
#' @param trans_matrix (numeric) m × m transition matrix.
#' @param trial_success A `multigrain_trial_success` object defining the custom
#'   power.
#' @param alpha (numeric) Overall one-sided significance level.
#'   Default is 0.025.
#' @param gamma Threshold below which a weight or edge should be considered for
#'   removing.
#' @param fixed_edge Logical matrix (m × m); `TRUE` where `trans_constraint`
#'   is non-`NA` (i.e. fixed by the constraint).
#' @param power_constraint Optional numeric vector of length m; only entries
#'   that are non-`NA` impose a marginal-power requirement. Defaults to `NULL`
#'   (no marginal constraints).
#' @param threshold Optional floor on the trial success measure. When supplied,
#'   edges are removed best-first by `.prune_edges_best_first()` and a removal
#'   is accepted even if it lowers the trial success, provided the result stays
#'   at or above `threshold`. `NULL` (the default) keeps the fixed-order
#'   `prune_edges()`, which never accepts a removal that lowers it.
#'
#' @returns A list with elements `hyp_weight` and `trans_matrix`, the pruned
#' graph, and `prune_loss`, the exact loss of trial success across the accepted
#' edge removals (zero unless `threshold` was supplied, and negative when the
#' removals raised the measure).
#'
#' @noRd
prune_graph <- function(
    pvals,
    hyp_weight,
    trans_matrix,
    trial_success,
    graph_constraint,
    alpha = 0.025,
    gamma = 1,
    power_constraint = NULL,
    threshold = NULL,
    verbose = c("info", "detail", "silent")
) {
    verbose <- rlang::arg_match(verbose)

    if (verbose != "silent") {
        cli::cli_progress_step("Pruning redundant weights and edges")
    }

    # Derive constraint metadata once
    fixed_w <- which(!is.na(graph_constraint$hyp_constraint))
    fixed_edge <- !is.na(graph_constraint$trans_constraint)
    tolerance <- graph_constraint_get_tolerance(graph_constraint)

    pruned_weights <- prune_hyp_weights(
        pvals = pvals,
        hyp_weight = hyp_weight,
        trans_matrix = trans_matrix,
        trial_success = trial_success,
        fixed_w = fixed_w,
        alpha = alpha,
        gamma = gamma,
        power_constraint = power_constraint,
        tolerance = tolerance
    )

    pruned <- if (is.null(threshold)) {
        prune_edges(
            pvals = pvals,
            hyp_weight = pruned_weights$hyp_weight,
            trans_matrix = pruned_weights$trans_matrix,
            trial_success = trial_success,
            fixed_edge = fixed_edge,
            alpha = alpha,
            gamma = gamma,
            power_best = pruned_weights$power_best,
            power_constraint = power_constraint,
            tolerance = tolerance
        )
    } else {
        .prune_edges_best_first(
            pvals = pvals,
            hyp_weight = pruned_weights$hyp_weight,
            trans_matrix = pruned_weights$trans_matrix,
            trial_success = trial_success,
            fixed_edge = fixed_edge,
            threshold = threshold,
            alpha = alpha,
            tolerance = tolerance
        )
    }

    list(
        hyp_weight = pruned$hyp_weight,
        trans_matrix = pruned$trans_matrix,
        prune_loss = pruned$prune_loss %||% 0
    )
}
