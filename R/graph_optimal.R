new_graph_optimal <- function(
    hyp_weight = double(),
    trans_matrix = matrix(double()),
    constraints = NULL,
    trial_success = NULL,
    power = NULL,
    solution = NULL,
    global_search = logical(),
    control = NULL,
    global_output = NULL,
    local_output = NULL,
    start_graph = NULL
) {
    structure(
        list(
            hyp_weight = hyp_weight,
            trans_matrix = trans_matrix,
            constraints = constraints,
            trial_success = trial_success,
            power = power,
            solution = solution,
            global_search = global_search,
            control = control,
            global_output = global_output,
            local_output = local_output,
            start_graph = start_graph
        ),
        class = "multigrain_graph_optimal"
    )
}

graph_optimal <- function(
    hyp_weight,
    trans_matrix,
    constraints = NULL,
    trial_success = NULL,
    power = NULL,
    solution = NULL,
    global_search = NULL,
    control = NULL,
    global_output = NULL,
    local_output = NULL,
    start_graph = NULL,
    stopping_rule = NULL,
    spending = NULL
) {
    check_double(hyp_weight)
    check_double_matrix(trans_matrix)
    check_graph_constraint(constraints, allow_null = TRUE)
    check_trial_success(trial_success, allow_null = TRUE)

    if (!is.null(solution)) {
        check_character(solution$opt_source)
        check_logical(solution$graph_valid)
    }

    check_logical(global_search, allow_null = TRUE)
    check_control(control, allow_null = TRUE)
    check_ga(global_output, allow_null = TRUE)
    check_nloptr(local_output, allow_null = TRUE)
    check_stopping_rule_gsd(stopping_rule, allow_null = TRUE)
    valid_spending <- is.list(spending) && rlang::is_named(spending)
    if (!is.null(spending) && !valid_spending) {
        cli::cli_abort("{.arg spending} must be a named list or {.code NULL}.")
    }

    # remove some elements from the GA and nloptr outputs as they bloat the
    # graph_optimal object

    # remove the @call slot from the global output. GA outputs are S4 objects
    # and the @call slot is recorded as an attribute
    attr(global_output, "call") <- NULL

    # remove the nloptr environment and eval_f since they're a bit chunky
    local_output$nloptr_environment <- NULL
    local_output$eval_f <- NULL

    out <- new_graph_optimal(
        hyp_weight = hyp_weight,
        trans_matrix = trans_matrix,
        constraints = constraints,
        trial_success = trial_success,
        power = power,
        solution = solution,
        global_search = global_search,
        control = control,
        global_output = global_output,
        local_output = local_output,
        start_graph = start_graph
    )

    # a group sequential result records its stopping rule; without one (and
    # for every fixed-sample result) the object is unchanged
    if (!is.null(stopping_rule)) {
        out$stopping_rule <- stopping_rule
    }

    # an optimised-spending result records its spending fields (build plan,
    # reconciliation 10); the fixed path never passes `spending`
    if (!is.null(spending)) {
        out[names(spending)] <- spending
    }

    out
}


#' Get the optimisation control
#'
#' Retrieve the settings used for the graph optimisation.
#'
#' @param graph_optimal A `multigrain_graph_optimal` object.
#'
#' @returns The `multigrain_control` object holding the optimisation settings.
#'
#' @export
#' @examples
#' # graph_optimal_example is an example optimised graph
#' graph_optimal_get_control(graph_optimal_example)
graph_optimal_get_control <- function(graph_optimal) {
    check_graph_optimal(graph_optimal)
    graph_optimal$control
}

summarise_power_object <- function(power_object) {
    if (is.null(power_object)) {
        return()
    }

    if (!is.null(power_object$trial_success)) {
        cli::cat_line(
            cli::style_underline(
                "\nValue of trial success measure"
            ),
            ":"
        )

        cli::cat_line(
            round(power_object$trial_success, 4)
        )
    }

    cli::cat_line(
        cli::style_underline("\nPower metrics"),
        ":"
    )

    # exact match: `$` would partially match `local_power_by_analysis`
    if (!is.null(power_object[["local_power"]])) {
        cli::cat_line("Power for each hypothesis:")
        print(round(power_object[["local_power"]], 4))
    }

    if (!is.null(power_object$exp_rejections)) {
        cli::cat_line(
            "Expected number of rejections: ",
            round(power_object$exp_rejections, 2)
        )
    }

    if (!is.null(power_object$disj_power)) {
        cli::cat_line(
            "Probability of at least one rejection: ",
            round(power_object$disj_power, 4)
        )
    }

    if (!is.null(power_object$conj_power)) {
        cli::cat_line(
            "Probability of rejecting all hypotheses: ",
            round(power_object$conj_power, 4)
        )
    }

    # group sequential fields; a fixed-sample power list has none of them
    if (!is.null(power_object$local_power_by_analysis)) {
        cli::cat_line("Power by analysis:")
        print(round(power_object$local_power_by_analysis, 4))
    }

    if (!is.null(power_object$time_distribution)) {
        cli::cat_line("Decision-time distribution:")
        print(round(power_object$time_distribution, 4))
    }

    if (!is.null(power_object$stop_distribution)) {
        cli::cat_line("Stopping-look distribution:")
        print(round(power_object$stop_distribution, 4))
        cli::cat_line(
            "Mean stopping look: ",
            round(power_object$mean_stop_look, 4)
        )
    }
}

summarise_solution_source <- function(solution_object) {
    if (is.null(solution_object)) {
        return()
    }

    cli::cat_line(
        cli::style_underline("\nSolution source"),
        ":"
    )

    src <- solution_object$opt_source
    if (identical(src, "local")) {
        cli::cat_line("Local optimisation (nloptr)")
    } else if (
        identical(src, "global") ||
            identical(src, "GA") ||
            identical(src, "GA_minN")
    ) {
        lab <- if (identical(src, "GA_minN")) {
            "Global optimisation (genetic algorithm; min N)"
        } else {
            "Global optimisation (genetic algorithm)"
        }
        cli::cat_line(lab)
    } else if (!is.null(src)) {
        cli::cat_line(src)
    }
}

# print method -----------------------------------------------------------

#' @export
print.multigrain_graph_optimal <- function(x, ...) {
    if (is.null(x)) {
        return()
    }

    cli::cat_line(cli::format_inline("{.cls {class(x)}}"))

    cli::cat_line(
        "Optimal graph found (given user-defined constraints ",
        "on graph and computational resources):"
    )

    cli::cat_line("\nHypothesis weights:")
    print(round(x$hyp_weight, 4))

    cli::cat_line("\nTransition matrix:")
    print(round(x$trans_matrix, 4))

    cli::cat_line("\nTrial success function:")
    cli::cat_line(x$trial_success$objective)

    if (!is.null(x$spending_param)) {
        .hsd_cat_print(x)
    }

    if (!is.null(x$stopping_rule)) {
        cli::cat_line("\nStopping rule:")
        .gsd_rule_cat_detail(x$stopping_rule)
    }

    if (!is.null(x$power) && !is.null(x$power$trial_success)) {
        cli::cat_line("\nValue of trial success measure:")
        cli::cat_line(format(round(x$power$trial_success, 4)))
    }

    invisible(x)
}

# summary method ---------------------------------------------------------

#' @export
summary.multigrain_graph_optimal <- function(object, ...) {
    if (is.null(object)) {
        return()
    }

    cli::cat_line(cli::style_bold("Optimisation summary"))

    cli::cat_line(
        cli::style_underline("\nOptimisation results"),
        ":"
    )

    cli::cat_line("\nHypothesis weights:")
    print(round(object$hyp_weight, 4))

    cli::cat_line("\nTransition matrix:")
    print(round(object$trans_matrix, 4))

    summary(object$trial_success)

    if (!is.null(object$spending_param)) {
        .hsd_cat_summary(object)
    }

    if (!is.null(object$stopping_rule)) {
        summary(object$stopping_rule)
    }

    summarise_power_object(object$power)

    summary(object$constraints)

    summarise_solution_source(object$solution)

    invisible(object)
}

# check helpers ----------------------------------------------------------

is_graph_optimal <- function(x) {
    inherits(x, "multigrain_graph_optimal")
}

check_graph_optimal <- function(
    graph_opt,
    arg = rlang::caller_arg(graph_opt),
    call = rlang::caller_env()
) {
    if (!missing(graph_opt) && is_graph_optimal(graph_opt)) {
        return(invisible(NULL))
    }

    rlang::stop_input_type(
        graph_opt,
        "a multigrain graph optimal object",
        allow_null = FALSE,
        arg = arg,
        call = call
    )
}
