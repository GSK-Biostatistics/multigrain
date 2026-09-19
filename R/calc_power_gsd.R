#' Calculate power for a group sequential graph-based multiple test procedure
#'
#' Group sequential counterpart of [calc_power_pvals()]. Applies a graph to an
#' object of repeated p-values built by [transform_pvalues_gsd()] and reports
#' power by hypothesis and by analysis, together with the distribution of the
#' analysis at which each hypothesis was declared rejected.
#'
#' @details `calc_power_pvals_gsd()` runs the group sequential kernel once and
#'   summarises both of its outputs: the rejection indicators and the decision
#'   times (the analysis at which the procedure declared each hypothesis
#'   rejected, `0` for never). Every quantity is on the analysis-index scale;
#'   calendar times are not stored by the package.
#'
#' @param pvals A `multigrain_pvals_gsd` object of repeated (or, per
#'   hypothesis, sequential) p-values, created with [transform_pvalues_gsd()].
#' @inheritParams is_graph_valid
#' @inheritParams rlang::args_dots_empty
#' @param alpha A single numeric value representing the overall one-sided
#'   significance level to test at, or `NULL` (the default) to use the level
#'   of `pvals`. A supplied value must be greater than `0` and at most the
#'   level the repeated p-values were built up to.
#' @param custom_power A list of user-defined power measures. Alternatively a
#'   single measure may be given on its own. Three kinds are accepted, and
#'   each receives a different matrix:
#'
#'   - a `multigrain_trial_success_gsd` object from [trial_success_gsd()] is
#'     evaluated by its compiled function on the integer matrix of decision
#'     times (simulations by hypotheses, `0` for never rejected);
#'   - a fixed-sample `multigrain_trial_success` object from [trial_success()]
#'     is evaluated on the logical matrix of rejections, so a rejection-only
#'     gain reports the same number here as [calc_power_pvals()] would;
#'   - a plain function is called once per simulated trial with that trial's
#'     integer vector of decision times, and the results are averaged.
#'
#'   If the list has no names, the measures are reported as `"func1"`,
#'   `"func2"`, and so on.
#' @param sum_to_one_constraint A logical value controlling whether to allow
#'   graphs where transition matrix rows are not constrained to sum to one,
#'   for example in a fixed sequence. Defaults to `TRUE`.
#' @inheritParams rlang::args_error_context
#'
#' @returns A list containing:
#'   * `local_power`: the proportion of simulations in which each hypothesis
#'     is rejected at any analysis.
#'   * `local_power_by_analysis`: an `m` by `K` matrix whose entry
#'     \eqn{(i, k)} is the proportion of simulations in which hypothesis
#'     \eqn{i} was rejected at analysis \eqn{k} or earlier. Its last column is
#'     `local_power`.
#'   * `exp_rejections`: the expected number of rejections.
#'   * `disj_power`: the probability of rejecting at least one hypothesis.
#'   * `conj_power`: the probability of rejecting all hypotheses.
#'   * `mean_decision_look`: for each hypothesis, the mean analysis index over
#'     the simulations in which it was rejected, and `NA` when it never was.
#'   * `time_distribution`: an `m` by `K + 1` matrix of proportions, with a
#'     `"never"` column followed by one column per analysis. Its rows sum to
#'     one, and one minus its `"never"` column equals `local_power` up to
#'     floating-point rounding.
#'   * `"..."`: one entry per measure in `custom_power`.
#'
#' @seealso [calc_power_pvals()] for the fixed-sample version and
#'   [graph_optimise_gsd()] for the optimiser.
#'
#' @export
#' @examples
#' # two hypotheses, an interim at 60% information and a final analysis
#' corr <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#'
#' set.seed(1)
#' raw <- simulate_pvalues_gsd(
#'     power_nominal = c(0.9, 0.8),
#'     corr_matrix = corr,
#'     info_frac = c(0.6, 1),
#'     nsim = 500
#' )
#'
#' pvals <- transform_pvalues_gsd(
#'     raw,
#'     spending = gsDesign::sfLDOF,
#'     grid_size = 128L
#' )
#'
#' # a rejection at the final analysis is worth 75% of one at the interim
#' gain <- trial_success_gsd(
#'     d(t1) + d(t2),
#'     d = c(1, 0.75),
#'     verbose = "silent"
#' )
#'
#' result <- calc_power_pvals_gsd(
#'     pvals,
#'     hyp_weight = c(0.5, 0.5),
#'     trans_matrix = matrix(c(0, 1, 1, 0), nrow = 2),
#'     custom_power = list(gain = gain)
#' )
#'
#' result$local_power_by_analysis
#' result$time_distribution
#' result$gain
calc_power_pvals_gsd <- function(
    pvals,
    hyp_weight,
    trans_matrix,
    ...,
    alpha = NULL,
    custom_power = NULL,
    sum_to_one_constraint = TRUE,
    call = rlang::caller_env()
) {
    check_pvals_gsd(pvals)
    check_double(hyp_weight)
    check_double_matrix(trans_matrix)
    rlang::check_dots_empty()
    alpha <- .gsd_check_alpha(alpha, pvals)
    check_logical(sum_to_one_constraint)

    if (
        !is_graph_valid(
            hyp_weight = hyp_weight,
            trans_matrix = trans_matrix,
            sum_to_one_constraint = sum_to_one_constraint
        )
    ) {
        cli::cli_abort(
            "The supplied {.arg hyp_weight} and {.arg trans_matrix} do not \\
            build a valid graph.",
            call = call
        )
    }

    custom_power <- .auto_name_custom_power_gsd(custom_power)
    # Every compiled gain is checked, group sequential or fixed-sample: both
    # index their matrix by hypothesis without a bounds check, so a gain over
    # more hypotheses than the graph has reads past the row. A fixed-sample
    # object carries no `K`, so `.gsd_check_gain_dims()` skips that branch.
    for (nm in names(custom_power)) {
        if (is_trial_success(custom_power[[nm]])) {
            .gsd_check_gain_dims(
                custom_power[[nm]],
                pvals,
                arg = paste0("custom_power$", nm),
                call = call
            )
        }
    }

    res <- graph_shortcut_gsd(
        pvals = .gsd_kernel_matrix(pvals, call = call),
        alpha = alpha,
        w = hyp_weight,
        G = trans_matrix,
        K = pvals$K
    )

    # Cannot happen with the kernel of design record section 4.3, but a
    # decision time above K would be an out-of-bounds read in a gain that
    # carries discount tables (section 10 item 6).
    if (max(res$time) > pvals$K) {
        cli::cli_abort(
            "The group sequential kernel reported a decision time of \\
            {max(res$time)} for a design with K = {pvals$K}.",
            call = call
        )
    }

    rej_counts <- rowSums(res$rejected)
    c(
        list(
            local_power = colMeans(res$rejected),
            local_power_by_analysis = .gsd_power_by_analysis(
                res$time,
                m = pvals$m,
                n_look = pvals$K
            ),
            exp_rejections = sum(res$rejected) / nrow(res$rejected),
            disj_power = mean(rej_counts > 0L),
            conj_power = mean(rej_counts == ncol(res$rejected)),
            mean_decision_look = .gsd_mean_decision_look(res$time),
            time_distribution = .gsd_time_distribution(
                res$time,
                m = pvals$m,
                n_look = pvals$K
            )
        ),
        .eval_custom_power_gsd(
            custom_power,
            res$rejected,
            res$time,
            call = call
        )
    )
}


# Cumulative local power by analysis: entry (i, k) is the proportion of trials
# in which hypothesis i had been rejected by analysis k.
.gsd_power_by_analysis <- function(time, m, n_look) {
    out <- matrix(NA_real_, nrow = m, ncol = n_look)

    for (i in seq_len(m)) {
        tau <- time[, i]
        for (k in seq_len(n_look)) {
            out[i, k] <- mean(tau >= 1L & tau <= k)
        }
    }

    dimnames(out) <- list(
        .gsd_hyp_labels(m),
        paste("analysis", seq_len(n_look))
    )
    out
}


# Mean analysis index over the trials in which a hypothesis was rejected;
# NA when it never was.
.gsd_mean_decision_look <- function(time) {
    vapply(
        seq_len(ncol(time)),
        function(i) {
            tau <- time[, i]
            rejected <- tau > 0L
            if (!any(rejected)) {
                return(NA_real_)
            }
            mean(tau[rejected])
        },
        numeric(1L)
    )
}


# Proportions of decision times at "never" and at each analysis; rows sum to 1.
.gsd_time_distribution <- function(time, m, n_look) {
    out <- matrix(NA_real_, nrow = m, ncol = n_look + 1L)

    for (i in seq_len(m)) {
        tau <- time[, i]
        for (k in seq.int(0L, n_look)) {
            out[i, k + 1L] <- mean(tau == k)
        }
    }

    dimnames(out) <- list(
        .gsd_hyp_labels(m),
        c("never", paste("analysis", seq_len(n_look)))
    )
    out
}


# Normalise custom_power to a named list; assign "funcN" to unnamed positions.
# Validates each element is a function, a fixed-sample trial success object or
# a group sequential one. Returns an empty list for NULL input so callers need
# not guard against NULL.
.auto_name_custom_power_gsd <- function(x, call = rlang::caller_env()) {
    if (is.null(x)) {
        return(list())
    }
    if (is.function(x) || is_trial_success(x)) {
        return(list(custom_power = x))
    }
    if (!is.list(x)) {
        cli::cli_abort(
            "{.arg custom_power} must be a function, a \\
            {.cls multigrain_trial_success_gsd} object, a \\
            {.cls multigrain_trial_success} object, or a list of these.",
            call = call
        )
    }
    for (i in seq_along(x)) {
        item <- x[[i]]
        if (!is.function(item) && !is_trial_success(item)) {
            cli::cli_abort(
                "Each element of {.arg custom_power} must be a function, a \\
                {.cls multigrain_trial_success_gsd} object or a \\
                {.cls multigrain_trial_success} object, \\
                not {.obj_type_friendly item}.",
                call = call
            )
        }
    }
    nms <- names(x)
    if (is.null(nms)) {
        nms <- character(length(x))
    }
    blank <- which(nms == "" | is.na(nms))
    nms[blank] <- sprintf("func%d", blank)
    names(x) <- nms
    x
}


# Evaluate each entry in a normalised custom_power list. The matrix an entry
# receives depends on its class, never on inheritance: a group sequential gain
# scores the decision times, a fixed-sample gain the rejection indicators, and
# a plain R function is called once per trial on its row of decision times
# (design record, section 4.7 [Rev 2026-09-18]).
.eval_custom_power_gsd <- function(
    custom_power,
    rejected,
    time,
    call = rlang::caller_env()
) {
    lapply(custom_power, function(item) {
        if (is_trial_success_gsd(item)) {
            item$func(time)
        } else if (is_trial_success(item)) {
            item$func(rejected)
        } else if (is.function(item)) {
            mean(
                vapply(
                    seq_len(nrow(time)),
                    \(i) item(time[i, ]),
                    FUN.VALUE = numeric(1L)
                )
            )
        } else {
            cli::cli_abort(
                "Each element of {.arg custom_power} must be a function, a \\
                {.cls multigrain_trial_success_gsd} object or a \\
                {.cls multigrain_trial_success} object.",
                call = call
            )
        }
    })
}
