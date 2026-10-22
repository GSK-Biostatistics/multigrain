checkArgs <- function(
    alpha,
    graph_constraint,
    power_obj,
    pvals,
    power_nominal,
    corr_matrix,
    startGraph,
    maxeval,
    popSize,
    global_search,
    max_generations,
    max_stagnation,
    power.nsim,
    power.nsim.global
) {
    # Check alpha
    if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
        stop(
            "`alpha` must be a numeric value between 0 and 1 (exclusive).",
            call. = FALSE
        )
    }

    # Check graph_constraint
    if (
        !inherits(graph_constraint, "graph_constraint") &&
            identical(graph_constraint, "none")
    ) {
        stop(
            "Provide a `graph_constraint` or specify `graph_constraint == 'none'`",
            call. = FALSE
        )
    }

    # Check power_obj
    if (!inherits(power_obj, "trial_success")) {
        stop(
            "Must provide a power objective function `power_obj`. Create one using `trial_success()`.",
            call. = FALSE
        )
    }

    # Check pvals
    if (!is.null(pvals) && !is.matrix(pvals)) {
        stop(
            "`pvals` must be a matrix if provided.",
            call. = FALSE
        )
    }

    # Validate correlation matrix
    if (!is.null(corr_matrix)) {
        if (!is.matrix(corr_matrix)) {
            stop(
                "`corr_matrix` must be a matrix.",
                call. = FALSE
            )
        }
        if (nrow(corr_matrix) != ncol(corr_matrix)) {
            stop(
                "`corr_matrix` must be square.",
                call. = FALSE
            )
        }
        if (any(diag(corr_matrix) != 1)) {
            stop(
                "The diagonal elements of `corr_matrix` must all be 1.",
                call. = FALSE
            )
        }
        if (any(corr_matrix < -1 | corr_matrix > 1)) {
            stop(
                "All elements of `corr_matrix` must be in the range [-1, 1].",
                call. = FALSE
            )
        }
    }

    # Validate startGraph
    if (!is.list(startGraph)) {
        stop(
            "`startGraph` must be a list of lists, each list containing `hyp_weight` and `trans_matrix`.",
            call. = FALSE
        )
    }

    if (
        !identical(
            startGraph,
            list(list(hyp_weight = NULL, trans_matrix = NULL))
        )
    ) {
        for (i in seq_along(startGraph)) {
            graph <- startGraph[[i]]
            if (
                !is.list(graph) ||
                    !all(c("hyp_weight", "trans_matrix") %in% names(graph))
            ) {
                stop(
                    "Each element of `startGraph` must be a list containing 'hyp_weight' and 'trans_matrix'.",
                    call. = FALSE
                )
            }

            if (!is_graph_valid(graph$hyp_weight, graph$trans_matrix)) {
                stop(
                    sprintf(
                        "The graph at index %d in `startGraph` is not valid.",
                        i
                    ),
                    call. = FALSE
                )
            }
        }
    }

    # Check global_search
    if (!is.logical(global_search) || length(global_search) != 1) {
        stop(
            "`global_search` must be a logical value.",
            call. = FALSE
        )
    }

    # Check maxeval, popSize, max_generations, max_stagnation, power.nsim,
    # power.nsim.global, print_level
    args_list <- list(
        maxeval,
        popSize,
        max_generations,
        max_stagnation,
        power.nsim,
        power.nsim.global
    )

    for (arg in args_list) {
        if (
            !is.null(arg) && (!is.numeric(arg) || length(arg) != 1 || arg <= 0)
        ) {
            stop(
                paste0(
                    "`maxeval`,
                    `popSize`,
                    `max_generations`,
                    `max_stagnation`,
                    `power.nsim`,
                    `power.nsim.global`",
                    "must be positive integers."
                ),
                call. = FALSE
            )
        }
    }
}
