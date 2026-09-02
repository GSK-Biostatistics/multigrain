#' Optimise graph-based multiple testing procedures
#'
#' Optimise the hypothesis weights and transition matrix of a graph-based
#' multiple testing procedure using [nloptr::nloptr()] local optimisation and,
#' optionally, a genetic algorithm using [GA::ga()] for global optimisation.
#'
#' The output is a graph where a specified objective function
#' - a *trial success measure* - is maximised under given constraints on the
#' graph structure, conditional on a p-value distribution supplied.
#'
#' @param pvals A numeric matrix of p-values. Each row represents a simulated
#'   trial, and each column corresponds to a hypothesis.
#' @param graph_constraint A `multigrain_graph_constraint` object containing
#'   constraints on the graph's weights and transition matrix. Created with
#'   [graph_constraint()].
#' @param trial_success A `multigrain_trial_success` object defining the trial
#'   success measure (power objective) function to maximize. Created with
#'   [trial_success()].
#' @param alpha Numeric, the overall one-sided significance level. Default is
#'   0.025.
#' @param start_graph Optional. Initial list of graphs suggested. Each graph is
#'   defined as a list containing starting weight vector `hyp_weight` and
#'   transition matrix `trans_matrix`. If `NULL`, default starting values are
#'   generated (currently with a Bonferroni-Holm graph).
#' @param power_nsim_local `r lifecycle::badge("deprecated")` `power_nsim_local`
#'   is no longer supported. You can set the number of simulations for the local
#'   optimisation with [control_nsim_local()].
#' @param power_nsim_global `r lifecycle::badge("deprecated")`
#'   `power_nsim_global` is no longer supported. You can set the number of
#'   simulations for the global optimisation with [control_nsim_global()].
#' @param global_search A logical indicate whether to perform a global
#'   optimisation before the local optimisation. Defaults to `TRUE`.
#' @param num_threads Number of threads to use for parallel execution of the
#'   shortcut algorithm. On shared systems (HPC clusters, login nodes), always
#'   explicitly set `num_threads` based on your resource allocation. Default is
#'   `1L` (serial execution).
#' @param control An optional `multigrain_control` object can be used to set
#'   various graph optimisation parameters. Created with [multigrain_control()].
#' @param verbose A logical controlling how much information to print (defaults
#'   to `TRUE`):
#'     * `FALSE`: no output
#'     * `TRUE`: show milestones
#' @param trace A logical controlling how much detailed optimisation information
#'   to print (defaults to `FALSE`):
#'     * `FALSE`: no detailed optimisation messages.
#'     * `TRUE`: all global and local optimisation messages.
#'   In addition to `trace`, [control_local()] and [control_global()] can be
#'   used to control the detail of information coming from the global and local
#'   optimisations separately.
#'
#' @return An object of class `multigrain_graph_optimal` containing the
#'   following elements:
#'     * `hyp_weight`: Optimised hypothesis weights (numeric vector).
#'     * `trans_matrix`: Optimised transition matrix (numeric matrix).
#'     * `constraints`: List of constraints used in the optimisation for weights
#'     and transition matrix.
#'     * `trial_success`: The trial success function used in the optimisation.
#'     * `power`: Power metrics for the optimised graph.
#'     * `solution`: A list containing:
#'       * `opt_source`: Source of the optimal solution (`local` or `global`).
#'       * `graph_valid`: Named logical vector indicating validity of the local
#'       and global solutions (`c("local" = TRUE/FALSE,
#'       "global" = TRUE/FALSE)`).
#'     * `global_search`: `TRUE` or `FALSE` indicating whether a global
#'     optimisation was performed.
#'     * `control`: A modified [multigrain_control()] object used. The values
#'     passed on by the user are complemented with contextual defaults.
#'     * `global_output`: Output from the genetic algorithm if global
#'     optimisation was performed.
#'     * `local_output`: Output from the NLOPT optimisation.
#'     * `start_graph`: Initial starting values used in the optimisation.
#'
#' @export
#' @examples
#'
#' # Generate test data
#' pvals <- simulate_pvalues(
#'   power_nominal = c(0.9, 0.85, 0.8, 0.75),
#'   alpha = 0.025,
#'   corr_matrix = diag(4),
#'   nsim = 5000
#' )
#'
#' # Create trial success function
#' ts <- trial_success(r1 + r2 + r3 + r4)
#'
#' \donttest{
#' # Optimise graph
#' result <- graph_optimise(
#'   pvals = pvals,
#'   graph_constraint = graph_constraint_free(4),
#'   trial_success = ts,
#'   alpha = 0.025,
#'   num_threads = 2
#' )
#' }
graph_optimise <- function(
    pvals,
    graph_constraint,
    trial_success,
    alpha = 0.025,
    start_graph = list(
        list(
            hyp_weight = NULL,
            trans_matrix = NULL
        )
    ),
    power_nsim_local = lifecycle::deprecated(),
    power_nsim_global = lifecycle::deprecated(),
    global_search = TRUE,
    num_threads = 1L,
    control = multigrain_control(),
    verbose = TRUE,
    trace = FALSE
) {
    check_double_matrix(pvals)
    check_graph_constraint(graph_constraint)
    check_trial_success(trial_success)
    rlang::check_number_decimal(alpha, min = 0, max = 1)
    rlang::check_number_whole(num_threads, min = 1)
    check_control(control)
    check_logical(global_search, allow_na = FALSE)
    check_logical(verbose, allow_na = FALSE)
    check_logical(trace, allow_na = FALSE)

    if (!missing(power_nsim_local)) {
        lifecycle::deprecate_stop(
            when = "0.2.0",
            what = I("`power_nsim_local`"),
            with = I("`control_nsim_local()`")
        )
    }

    if (!missing(power_nsim_global)) {
        lifecycle::deprecate_stop(
            when = "0.2.0",
            what = I("`power_nsim_global`"),
            with = I("`control_nsim_global()`")
        )
    }

    # Create objective functions
    if (trial_success$m != ncol(pvals)) {
        stop(
            "`multigrain_trial_success` dimensions (m) must match the number of columns in `pvals`.", # nolint
            call. = FALSE
        )
    }

    .validate_start_graphs(start_graph, ncol(pvals), where = "graph_optimise")

    control <- control_prepare(
        control,
        pvals = pvals,
        trace = trace
    )

    ga_result <- NULL
    x0_for_local <- NULL

    if (global_search) {
        ga_result <- .graph_optimise_ga(
            pvals = pvals,
            graph_constraint = graph_constraint,
            trial_success = trial_success,
            alpha = alpha,
            nsim = control$nsim_global,
            start_graph = start_graph,
            global_opts = control$global_opt,
            num_threads = num_threads,
            verbose = verbose
        )

        x0_for_local <- pmin(
            pmax(
                ga_result$ga_output@solution[1, ],
                0
            ),
            1
        )
    } else if (!.is_default_start_graph(start_graph)) {
        x0_for_local <- .build_start_matrix(
            graph_constraint,
            start_graph
        )[1, , drop = TRUE]
    }

    local_result <- .graph_optimise_local(
        pvals = pvals,
        graph_constraint = graph_constraint,
        trial_success = trial_success,
        alpha = alpha,
        nsim = control$nsim_local,
        local_opts = control$local_opt,
        num_threads = num_threads,
        x0 = x0_for_local,
        verbose = verbose
    )

    best_graph_result <- choose_graph(ga_result, local_result)

    clean_graph <- prune_graph(
        pvals = pvals,
        hyp_weight = best_graph_result$hyp_weight,
        trans_matrix = best_graph_result$trans_matrix,
        alpha = alpha,
        graph_constraint = graph_constraint,
        trial_success = trial_success,
        gamma = 1,
        verbose = verbose
    )

    if (verbose) {
        cli::cli_progress_step(
            "Evaluating trial success of pruned graph"
        )
    }

    final_power <- calc_power_pvals(
        pvals = pvals,
        alpha = alpha,
        hyp_weight = clean_graph$hyp_weight,
        trans_matrix = clean_graph$trans_matrix,
        custom_power = list(trial_success = trial_success)
    )

    # Apply hypothesis names
    hyp_names <- names(graph_constraint$hyp_constraint)
    names(clean_graph$hyp_weight) <- hyp_names
    dimnames(clean_graph$trans_matrix) <- list(hyp_names, hyp_names)

    graph_optimal(
        hyp_weight = clean_graph$hyp_weight,
        trans_matrix = clean_graph$trans_matrix,
        constraints = graph_constraint,
        trial_success = trial_success,
        power = final_power,
        solution = list(
            opt_source = best_graph_result$source,
            graph_valid = c(
                local = isTRUE(local_result$is_graph_valid),
                global = if (!is.null(ga_result)) {
                    isTRUE(ga_result$is_graph_valid)
                } else {
                    NA
                }
            )
        ),
        global_search = global_search,
        control = control,
        global_output = ga_result$ga_output,
        local_output = local_result$local_output,
        start_graph = start_graph
    )
}

#' @export
#' @rdname graph_optimise
#' @usage NULL
graph_optimize <- graph_optimise


#' Optimise graph-based multiple testing procedures using [GA::ga()]
#'
#' This internal function optimises the hypothesis weights and transition matrix
#' of a graph-based multiple testing procedure  using a genetic algorithm using
#' [GA::ga()] for global optimisation. The function produces a graph where a
#' specified objective function - a *trial success measure* - is maximised under
#' given constraints on the graph structure, conditional on a p-value
#' distribution supplied.
#'
#' @inheritParams graph_optimise
#'
#' @noRd
.graph_optimise_ga <- function(
    pvals,
    graph_constraint,
    trial_success,
    alpha,
    nsim, # defaults to control$nsim_global
    global_opts, # defaults to control$global_opt
    num_threads = 1L,
    start_graph = list(
        list(
            hyp_weight = NULL,
            trans_matrix = NULL
        )
    ),
    verbose = TRUE
) {
    if (verbose) {
        cli::cli_progress_step("Running global optimization")
    }

    x0 <- .build_start_matrix(graph_constraint, start_graph)

    pvals_sampled <- pvals[sample(nsim), ]

    immutable_global_args <- list(
        type = "real-valued",
        fitness = create_obj_func(
            m = trial_success$m,
            power_criterion = trial_success$func,
            hyp_constraint = graph_constraint$hyp_constraint,
            trans_constraint = graph_constraint$trans_constraint,
            alpha = alpha,
            pvals = pvals_sampled,
            num_threads = num_threads
        ),
        lower = rep(0, ncol(x0)),
        upper = rep(1, ncol(x0)),
        population = .cauchy_population,
        mutation = .make_cauchy_mutation_multi(
            p_param_mutate = 0.1,
            scale = 1.0
        ),
        optim = TRUE,
        suggestions = x0
    )

    ga_args <- utils::modifyList(immutable_global_args, global_opts)

    ga_res <- do.call(GA::ga, ga_args)
    best_raw <- ga_res@solution[1, ]
    sol <- param_to_solution(best_raw, graph_constraint, process = TRUE)

    is_valid <- is_graph_valid(sol$hyp_weight, sol$trans_matrix)

    ga_subset_power <- NULL
    ga_trial_success <- NULL

    if (is_valid) {
        if (verbose) {
            cli::cli_progress_step(
                "Evaluating trial success of globally optimised graph"
            )
        }
        ga_subset_power <- calc_power_pvals(
            pvals_sampled, # calculated using nsim_global samples
            hyp_weight = sol$hyp_weight,
            alpha = alpha,
            trans_matrix = sol$trans_matrix,
            custom_power = list(trial_success = trial_success)
        )

        rej_mat <- graph_shortcut(
            pvals = pvals, # calculated on full sample
            alpha = alpha,
            w = sol$hyp_weight,
            G = sol$trans_matrix
        )

        ga_trial_success <- trial_success$func(rej_mat)
    }

    list(
        ga_hyp_weight = sol$hyp_weight,
        ga_trans_matrix = sol$trans_matrix,
        ga_trial_success = ga_trial_success,
        ga_subset_power = ga_subset_power,
        is_graph_valid = is_valid,
        ga_output = ga_res
    )
}


# Helper functions for modified hybrid-FES optimisation
# Custom mutation function
# DEPRECATED from v0.2.0 (only included for benchmarking)
.cauchy_mutation <- function(object, parent) {
    esch_mutation_rcpp(parent, object@lower, object@upper)
}

# Custom population function
.cauchy_population <- function(object) {
    esch_population_rcpp(
        popSize = object@popSize,
        nBits = length(object@lower),
        lower = object@lower,
        upper = object@upper
    )
}


#' Optimise graph-based multiple testing procedures
#'
#' Optimises the hypothesis weights and transition matrix of a graph-based
#' multiple testing procedure. It uses [nloptr::nloptr()] local optimisation. It
#' produces a graph where a specified objective function
#' - a *trial success measure* - is maximised under given constraints on the
#' graph structure, conditional on a p-value distribution supplied.
#'
#' @inheritParams graph_optimise
#'
#' @noRd
.graph_optimise_local <- function(
    pvals,
    graph_constraint,
    trial_success,
    alpha,
    local_opts,
    num_threads = 1L,
    nsim = nrow(pvals),
    x0 = NULL,
    verbose = TRUE
) {
    if (verbose) {
        cli::cli_progress_step("Running local optimization")
    }

    if (is.null(x0)) {
        x0 <- create_start_params(graph_constraint)
    }

    pvals_sampled <- pvals[sample(nsim), ]

    obj_fun <- create_obj_func(
        trial_success$m,
        trial_success$func,
        graph_constraint$hyp_constraint,
        graph_constraint$trans_constraint,
        alpha,
        pvals_sampled,
        num_threads
    )
    nlopt_obj_func <- function(x) {
        -obj_fun(x)
    }

    nlopt_result <- nloptr::nloptr(
        x0 = x0,
        eval_f = nlopt_obj_func,
        lb = rep(0, length(x0)),
        ub = rep(1, length(x0)),
        opts = local_opts
    )

    sol <- param_to_solution(
        nlopt_result$solution,
        graph_constraint,
        process = TRUE
    )
    sol <- repair_graph(sol$hyp_weight, sol$trans_matrix, graph_constraint)

    is_valid <- is_graph_valid(sol$hyp_weight, sol$trans_matrix)

    local_subset_power <- NULL
    local_trial_success <- NULL

    if (is_valid) {
        if (verbose) {
            cli::cli_progress_step(
                "Evaluating trial success of locally optimised graph"
            )
        }
        local_subset_power <- calc_power_pvals(
            pvals_sampled, # calculated using nsim_local samples
            hyp_weight = sol$hyp_weight,
            alpha = alpha,
            trans_matrix = sol$trans_matrix,
            custom_power = list(trial_success = trial_success)
        )

        rej_mat <- graph_shortcut(
            pvals = pvals, # calculated on full sample
            alpha = alpha,
            w = sol$hyp_weight,
            G = sol$trans_matrix
        )

        local_trial_success <- trial_success$func(rej_mat)
    }

    list(
        local_hyp_weight = sol$hyp_weight,
        local_trans_matrix = sol$trans_matrix,
        local_trial_success = local_trial_success,
        local_subset_power = local_subset_power,
        is_graph_valid = is_valid,
        local_output = nlopt_result
    )
}
