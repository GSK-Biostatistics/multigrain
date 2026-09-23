# Group sequential twins of the greedy pruning helpers in
# `R/post_optim_processing.R`. The loops, the acceptance rule and the mass
# redistribution are those of the fixed-sample versions; only the evaluation
# changes, from `calc_power_pvals()` to `calc_power_pvals_gsd()`.
# `.redistribute_mass()` and `.marginal_violated()` are reused unchanged.


#' Evaluate a candidate group sequential graph and accept if it improves on
#' the current best
#'
#' Computes trial-success power for the candidate graph. Returns the candidate
#' values if power does not decrease and marginal constraints are satisfied;
#' otherwise returns the current best unchanged.
#'
#' @param pvals A `multigrain_pvals_gsd` object.
#' @param hyp_weight Candidate hypothesis weight vector.
#' @param trans_matrix Candidate transition matrix.
#' @param alpha Significance level. Default is 0.025.
#' @param trial_success A `multigrain_trial_success_gsd` object.
#' @param power_best Current best trial-success power (scalar).
#' @param constrained_idx Integer vector; which hypotheses carry marginal
#'   constraints.
#' @param power_constraint Numeric vector of marginal thresholds.
#'
#' @returns A list with `hyp_weight`, `trans_matrix`, `power_best`, and
#'   `accepted` (logical).
#' @noRd
.try_prune_gsd <- function(
    pvals,
    hyp_weight,
    trans_matrix,
    trial_success,
    power_best,
    constrained_idx,
    power_constraint,
    alpha = 0.025
) {
    power_all <- calc_power_pvals_gsd(
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
#' the estimated group sequential gain or violate marginal constraints.
#'
#' @param pvals A `multigrain_pvals_gsd` object.
#' @param hyp_weight (numeric) Vector of hypothesis weights (length m).
#' @param trans_matrix (numeric) m × m transition matrix (held fixed during
#'   weight pruning).
#' @param trial_success A `multigrain_trial_success_gsd` object defining the
#'   custom power.
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
prune_hyp_weights_gsd <- function(
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
    power_all <- calc_power_pvals_gsd(
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

        result <- .try_prune_gsd(
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
#' decrease the estimated group sequential gain or violate marginal
#' constraints.
#'
#' @param pvals A `multigrain_pvals_gsd` object.
#' @param hyp_weight (numeric) Vector of hypothesis weights (held fixed during
#'   edge pruning).
#' @param trans_matrix (numeric) m × m transition matrix.
#' @param trial_success A `multigrain_trial_success_gsd` object defining the
#'   custom power.
#' @param fixed_edge Logical matrix (m × m); `TRUE` where `trans_constraint`
#'   is non-`NA` (i.e. fixed by the constraint).
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
prune_edges_gsd <- function(
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

            result <- .try_prune_gsd(
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


#' Prune graph hypothesis and transition weights according to a group
#' sequential p-value object and trial success measure.
#'
#' @param pvals A `multigrain_pvals_gsd` object.
#' @param hyp_weight (numeric) Vector of initial hypothesis weights (length m).
#' @param trans_matrix (numeric) m × m transition matrix.
#' @param trial_success A `multigrain_trial_success_gsd` object defining the
#'   custom power.
#' @param graph_constraint A `multigrain_graph_constraint` object.
#' @param alpha (numeric) Overall one-sided significance level.
#'   Default is 0.025.
#' @param gamma Threshold below which a weight or edge should be considered for
#'   removing.
#' @param power_constraint Optional numeric vector of length m; only entries
#'   that are non-`NA` impose a marginal-power requirement. Defaults to `NULL`
#'   (no marginal constraints).
#' @inheritParams graph_optimise_gsd
#'
#' @returns A list with elements `hyp_weight` and `trans_matrix`, the pruned
#' graph.
#'
#' @noRd
prune_graph_gsd <- function(
    pvals,
    hyp_weight,
    trans_matrix,
    trial_success,
    graph_constraint,
    alpha = 0.025,
    gamma = 1,
    power_constraint = NULL,
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

    pruned_weights <- prune_hyp_weights_gsd(
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

    pruned <- prune_edges_gsd(
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

    list(
        hyp_weight = pruned$hyp_weight,
        trans_matrix = pruned$trans_matrix
    )
}
