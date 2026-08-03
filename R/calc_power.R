#' Calculate power for a graph-based multiple test procedure using p-values
#'
#' Calculate multiplicity-adjusted (marginal) power values for a graph-based
#' multiple test procedure applied to a matrix of p-values. Provide both the
#' `hyp_weight` vector and transition matrix `trans_matrix` for the graph you
#' wish to evaluate.
#'
#' @details `calc_power_pvals()` calculates several power metrics:
#'   * _local power_ (probability to reject each individual hypothesis),
#'   * _disjunctive power_ (the probability to reject at least one hypothesis),
#'   * _conjunctive power_ (the probability to reject all hypotheses), and
#'   * the expected number of rejections.
#' Optionally, you can specify user-defined trial success criteria via
#' `custom_power`.
#'
#' @inheritParams graph_optimise pvals
#' @inheritParams is_graph_valid hyp_weight trans_matrix
#' @inheritParams rlang::args_dots_empty
#' @inheritParams graph_optimise alpha
#' @param custom_power A list of user-defined power functions. Alternatively, a
#'   single function or `multigrain_trial_success` object can be provided.
#'   There are two ways to specify this:
#'
#'   - Anonymous functions can be provided to specify the success criteria.
#'     Functions must take one simulation's logical vector of results as an
#'     input, and return a scalar. For example, if one is interested in the
#'     power to reject hypotheses 1 and 3 one could specify:
#'     `f = function(x) {x[1] && x[3]}`. If the power of rejecting hypotheses
#'     1 and 2 is also of interest one would use an (optionally named) list:
#'     `f = list(
#'         power1and3 = function(x) {x[1] && x[3]},
#'         power1and2 = function(x) {x[1] && x[2]}
#'       )`.
#'     If the list has no names, the functions will be referenced as `"func1"`,
#'     `"func2"`, etc. in the output. The user can also provide a
#'     `multigrain_trial_success` object instead (resulting in a faster
#'     calculation).
#'
#'   - Instead of anonymous functions, one can pass a `multigrain_trial_success`
#'     object (or within a list as with the anonymous functions). See
#'     [trial_success()] and the examples.
#' @param sum_to_one_constraint A logical value controlling whether to allow
#'   graphs where transition matrix rows are not constrained to sum to one,
#'   for example in a fixed sequence. Defaults to `TRUE`.
#' @inheritParams rlang::args_error_context call
#'
#' @returns A list containing:
#'   * `local_power`: The local power for each hypothesis: the proportion of
#'     simulations in which each hypothesis is rejected.
#'   * `exp_rejections`: The expected number of rejections across all
#'     simulations.
#'   * `disj_power`: Disjunctive power: the probability of rejecting at least
#'     one hypothesis.
#'   * `conj_power`: Conjunctive power: the probability of rejecting all
#'     hypotheses.
#'   * `"..."`: Additional results as specified by user-defined success
#'     functions in `custom_power`.
#'
#' @export
#' @examples
#'
#' # First we simulate our p-value distribution
#'
#' power_nominal <- c(0.90, 0.87, 0.73)
#' corr_matrix <- matrix(
#'   c(
#'     1, 0.2, 0.2,
#'     0.2, 1, 0.2,
#'     0.2, 0.2, 1
#'   ),
#'   nrow = 3,
#'   byrow = TRUE
#' )
#'
#' pvals <- simulate_pvalues(
#'   power_nominal = power_nominal,
#'   corr_matrix = corr_matrix
#' )
#'
#' # Second we construct our graph (in this case a fixed sequence)
#'
#' hyp_weights <- c(1, 0, 0)
#' trans_matrix <- matrix(
#'   c(
#'     0, 1, 0,
#'     0, 0, 1,
#'     0, 0, 0
#'   ),
#'   nrow = 3,
#'   byrow = TRUE
#' )
#'
#' # Third, we construct a list of metrics we wish to evaluate our graph with
#'
#' trial_success_measure <- trial_success((r1 && r2) || r3)
#' power_metrics <- list(
#'   trial_success = trial_success_measure,
#'   average_power = function(x) {
#'     x[1] + x[2] + x[3]
#'   }
#' )
#'
#' # Finally we calculate the power of graph conditional on p-value distribution
#'
#' result <- calc_power_pvals(
#'   pvals = pvals,
#'   hyp_weight = hyp_weights,
#'   trans_matrix = trans_matrix,
#'   custom_power = power_metrics,
#'   sum_to_one_constraint = FALSE # As third row of graph does not sum to 1
#' )
#'
#' result
calc_power_pvals <- function(
    pvals,
    hyp_weight,
    trans_matrix,
    ...,
    alpha = 0.025,
    custom_power = NULL,
    sum_to_one_constraint = TRUE,
    call = rlang::caller_env()
) {
    check_double_matrix(pvals)
    check_double(hyp_weight)
    check_double_matrix(trans_matrix)
    rlang::check_dots_empty()
    rlang::check_number_decimal(alpha, min = 0)
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

    custom_power <- .auto_name_custom_power(custom_power)

    rej_mat <- graph_shortcut(
        pvals = pvals,
        alpha = alpha,
        w = hyp_weight,
        G = trans_matrix
    )

    rej_counts <- rowSums(rej_mat)
    c(
        list(
            local_power = colMeans(rej_mat),
            exp_rejections = sum(rej_mat) / nrow(rej_mat),
            disj_power = mean(rej_counts > 0L),
            conj_power = mean(rej_counts == ncol(rej_mat))
        ),
        .eval_custom_power(custom_power, rej_mat)
    )
}


# Normalise custom_power to a named list; assign "funcN" to unnamed positions.
# Validates each element is a function or multigrain_trial_success.
# Returns an empty list for NULL input so callers need not guard against NULL.
.auto_name_custom_power <- function(x, call = rlang::caller_env()) {
    if (is.null(x)) {
        return(list())
    }
    if (is.function(x) || is_trial_success(x)) {
        return(list(custom_power = x))
    }
    if (!is.list(x)) {
        cli::cli_abort(
            "{.arg custom_power} must be a function, a \\
            {.cls multigrain_trial_success} object, or a list of these.",
            call = call
        )
    }
    for (i in seq_along(x)) {
        item <- x[[i]]
        if (!is.function(item) && !is_trial_success(item)) {
            cli::cli_abort(
                "Each element of {.arg custom_power} must be a function or a \\
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


# Evaluate each entry in a normalised custom_power list against `rej_mat`.
.eval_custom_power <- function(
    custom_power,
    rej_mat,
    call = rlang::caller_env()
) {
    lapply(custom_power, function(item) {
        if (is_trial_success(item)) {
            item$func(rej_mat)
        } else if (is.function(item)) {
            mean(
                vapply(
                    seq_len(nrow(rej_mat)),
                    \(i) item(rej_mat[i, ]),
                    FUN.VALUE = numeric(1L)
                )
            )
        } else {
            cli::cli_abort(
                "Each element of {.arg custom_power} must be a function or a \\
                {.cls multigrain_trial_success} object.",
                call = call
            )
        }
    })
}
