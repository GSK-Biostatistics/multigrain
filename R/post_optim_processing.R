# (Optimised) parameters to full solution (option to process to select epsilon
# edges)
param_to_solution <- function(
    optimised_params,
    graph_constraint,
    process = FALSE
) {
    hyp_constraint <- graph_constraint$hyp_constraint
    trans_constraint <- graph_constraint$trans_constraint

    theta <- split_theta(optimised_params, hyp_constraint)
    w_sol <- recover_full_weights(theta$w_pars, hyp_constraint)
    G_sol <- recover_full_trans_matrix(theta$g_pars, trans_constraint)

    if (process) {
        # Convert <1e-4 w & G parameters to zero and epsilon edges
        w_sol[w_sol < 1e-4] <- 0
        w_sol[w_sol > 1e-4 & w_sol < 1e-3] <- 0.001
        fixed_w <- which(!is.na(hyp_constraint))
        w_sol <- normalise_sum(w_sol, fixed_idx = fixed_w)

        G_sol[G_sol < 1e-5] <- 0
        G_sol[G_sol > 1e-5 & G_sol < 1e-3] <- 0.001
        m <- nrow(G_sol)
        for (row_i in seq_len(m)) {
            fixed_in_row <- c(row_i, which(!is.na(trans_constraint[row_i, ])))
            G_sol[row_i, ] <- normalise_sum(
                G_sol[row_i, ],
                fixed_idx = fixed_in_row
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
#' @return A list with `hyp_weight` and `trans_matrix`, both guaranteed to
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
    hyp_weight <- normalise_sum(hyp_weight, fixed_idx = fixed_w)

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
            fixed_idx = fixed_in_row
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
#' @return Integer vector of recipient indices (may be length 0).
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
#'
#' @return Numeric vector of same length, with `vec[drop_idx] == 0`.
#' @noRd
.redistribute_mass <- function(vec, drop_idx, fixed_idx) {
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

    normalise_sum(vec, fixed_idx = fixed_idx)
}


#' Check whether marginal power constraints are violated
#'
#' @param marg_power Numeric vector of marginal powers (length m).
#' @param constrained_idx Integer vector of indices with constraints.
#' @param power_constraint Numeric vector of required marginal powers.
#'
#' @return Logical scalar; `TRUE` if any constraint is violated.
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
#' @param alpha Significance level.
#' @param trial_success A `multigrain_trial_success` object.
#' @param power_best Current best trial-success power (scalar).
#' @param constrained_idx Integer vector; which hypotheses carry marginal
#'   constraints.
#' @param power_constraint Numeric vector of marginal thresholds.
#'
#' @return A list with `hyp_weight`, `trans_matrix`, `power_best`, and
#'   `accepted` (logical).
#' @noRd
.try_prune <- function(
    pvals,
    hyp_weight,
    trans_matrix,
    alpha,
    trial_success,
    power_best,
    constrained_idx,
    power_constraint
) {
    power_all <- calc_power_pvals(
        pvals,
        hyp_weight = hyp_weight,
        alpha = alpha,
        trans_matrix = trans_matrix,
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
#' @param gamma Threshold below which a weight should be considered for removal.
#' @param power_constraint Optional numeric vector of length m; only entries
#'   that are non-`NA` impose a marginal-power requirement. Defaults to `NULL`
#'   (no marginal constraints).
#'
#' @return A list with elements `hyp_weight` (pruned weights), `trans_matrix`
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
    power_constraint = NULL
) {
    power_all <- calc_power_pvals(
        pvals,
        hyp_weight = hyp_weight,
        alpha = alpha,
        trans_matrix = trans_matrix,
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
            fixed_idx = fixed_w
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
#' @param gamma Threshold below which an edge should be considered for removal.
#' @param power_constraint Optional numeric vector of length m; only entries
#'   that are non-`NA` impose a marginal-power requirement. Defaults to `NULL`
#'   (no marginal constraints).
#'
#' @return A list with elements `hyp_weight` (unchanged), `trans_matrix`
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
    power_constraint = NULL
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
                fixed_idx = fixed_in_row
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


#' Prune graph hypothesis and transition weights according to p-value
#' distribution and trial success measure.
#'
#' @param pvals (numeric) Numeric matrix of p-values (n.sim x m).
#' @param hyp_weight (numeric) Vector of initial hypothesis weights (length m).
#' @param trans_matrix (numeric) m × m transition matrix.
#' @param trial_success A `multigrain_trial_success` object defining the custom
#'   power.
#' @param alpha (numeric) Overall one-sided significance level.
#' @param gamma Threshold below which a weight or edge should be considered for
#'   removing.
#' @param fixed_edge Logical matrix (m × m); `TRUE` where `trans_constraint`
#'   is non-`NA` (i.e. fixed by the constraint).
#' @param power_constraint Optional numeric vector of length m; only entries
#'   that are non-`NA` impose a marginal-power requirement. Defaults to `NULL`
#'   (no marginal constraints).
#'
#' @return A list with elements `hyp_weight` and `trans_matrix`, the pruned
#' graph.
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
    verbose = FALSE
) {
    if (verbose) {
        cli::cli_progress_step("Pruning redundant weights and edges")
    }

    # Derive constraint metadata once
    fixed_w <- which(!is.na(graph_constraint$hyp_constraint))
    fixed_edge <- !is.na(graph_constraint$trans_constraint)

    pruned_weights <- prune_hyp_weights(
        pvals = pvals,
        hyp_weight = hyp_weight,
        trans_matrix = trans_matrix,
        trial_success = trial_success,
        fixed_w = fixed_w,
        alpha = alpha,
        gamma = gamma,
        power_constraint = power_constraint
    )

    pruned <- prune_edges(
        pvals = pvals,
        hyp_weight = pruned_weights$hyp_weight,
        trans_matrix = pruned_weights$trans_matrix,
        trial_success = trial_success,
        fixed_edge = fixed_edge,
        alpha = alpha,
        gamma = gamma,
        power_best = pruned_weights$power_best,
        power_constraint = power_constraint
    )

    list(
        hyp_weight = pruned$hyp_weight,
        trans_matrix = pruned$trans_matrix
    )
}
