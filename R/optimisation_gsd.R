#' Optimise graph-based multiple testing procedures for a group sequential
#' design
#'
#' Group sequential counterpart of [graph_optimise()]. Optimise the hypothesis
#' weights and transition matrix of a graph-based multiple testing procedure
#' for a design with several analyses, using [nloptr::nloptr()] local
#' optimisation and, optionally, a genetic algorithm using [GA::ga()] for
#' global optimisation.
#'
#' The output is a graph where a specified objective function - the
#' *trial success measure* - is maximised under given constraints on the graph
#' structure, conditional on the repeated p-values supplied. Because the gain
#' sees the analysis at which each hypothesis was rejected, an early rejection
#' can be worth more than a late one.
#'
#' @details The search, the encoding and the pruning are those of
#'   [graph_optimise()]. Only the evaluation changes: the group sequential
#'   kernel runs the fixed-sample cascade once per analysis, carrying the
#'   updated graph forward, and reports the analysis at which each hypothesis
#'   was declared rejected. Those decision times are what `trial_success`
#'   scores.
#'
#' @param pvals A `multigrain_pvals_gsd` object of repeated (or, per
#'   hypothesis, sequential) p-values, created with [transform_pvalues_gsd()].
#' @param graph_constraint A `multigrain_graph_constraint` object containing
#'   constraints on the graph's weights and transition matrix. Created with
#'   [graph_constraint()].
#' @param trial_success A `multigrain_trial_success_gsd` object defining the
#'   trial success measure (gain function) to maximise. Created with
#'   [trial_success_gsd()]. A fixed-sample object from [trial_success()] is
#'   not accepted; write a rejection-only gain such as `r1 + r2` with
#'   [trial_success_gsd()] instead.
#' @inheritParams rlang::args_dots_empty
#' @param alpha A single numeric value representing the overall one-sided
#'   significance level, or `NULL` (the default) to use the level of `pvals`.
#'   A supplied value must be greater than `0` and at most the level the
#'   repeated p-values were built up to: above that they are capped at 1, so
#'   no allocation beyond it could ever reject.
#' @param start_graph Optional. Initial list of graphs suggested. Each graph is
#'   defined as a list containing starting weight vector `hyp_weight` and
#'   transition matrix `trans_matrix`. If `NULL`, default starting values are
#'   generated (currently with a Bonferroni-Holm graph).
#' @param global_search A logical indicating whether to perform a global
#'   optimisation before the local optimisation. Defaults to `TRUE`.
#' @param num_threads Number of threads to use for parallel execution of the
#'   shortcut algorithm. On shared systems (HPC clusters, login nodes), always
#'   explicitly set `num_threads` based on your resource allocation. Default is
#'   `1L` (serial execution).
#' @param control An optional `multigrain_control` object can be used to set
#'   various graph optimisation parameters. Created with [multigrain_control()].
#' @param verbose An optional string controlling verbosity (`"detail"` >
#'   `"info"` > `"silent"`). Verbosity can also be set at package level with the
#'   `multigrain_verbosity` option (see [multigrain_verbosity()]):
#'     * `"info"` (default): will only show milestones / informational messages
#'     highlighting the progress of the optimisation at coarse-grained level.
#'     * `"detail"`: will show milestones and information about fine-grained
#'     optimisation events.
#'     * `"silent"`: no information about the progress of the optimisation is
#'     printed to the console. Errors and warnings are still thrown normally.
#'
#' @returns A `multigrain_graph_optimal` object containing:
#'
#' * `hyp_weight`: Optimised hypothesis weights (numeric vector).
#' * `trans_matrix`: Optimised transition matrix (numeric matrix).
#' * `constraints`: List of constraints used in the optimisation for weights
#'   and transition matrix.
#' * `trial_success`: The trial success function used in the optimisation.
#' * `power`: The [calc_power_pvals_gsd()] list for the optimised graph, so it
#'   carries power by analysis, the decision-time distribution and the value
#'   of the trial success measure as well as the fixed-sample power metrics.
#' * `solution`: A list containing:
#'   * `opt_source`: Source of the optimal solution (`local` or `global`).
#'   * `graph_valid`: Named logical vector indicating validity of the local
#'     and global solutions (`c("local" = TRUE/FALSE,
#'     "global" = TRUE/FALSE)`).
#' * `global_search`: `TRUE` or `FALSE` indicating whether a global
#'   optimisation was performed.
#' * `control`: A modified [multigrain_control()] object used. The values
#'   passed on by the user are complemented with contextual defaults.
#' * `global_output`: Output from the genetic algorithm if global
#'   optimisation was performed.
#' * `local_output`: Output from the NLOPT optimisation.
#' * `start_graph`: Initial starting values used in the optimisation.
#'
#' @seealso [graph_optimise()] for the fixed-sample version,
#'   [transform_pvalues_gsd()] for the p-value object and
#'   [trial_success_gsd()] for the gain.
#'
#' @export
#' @examples
#' # two hypotheses, an interim at 50% information and a final analysis
#' corr <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#'
#' set.seed(1)
#' raw <- simulate_pvalues_gsd(
#'     power_nominal = c(0.9, 0.8),
#'     corr_matrix = corr,
#'     info_frac = c(0.5, 1),
#'     nsim = 2000
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
#'     r1 * d(t1) + r2 * d(t2),
#'     d = c(1, 0.75)
#' )
#'
#' \donttest{
#' result <- graph_optimise_gsd(
#'     pvals = pvals,
#'     graph_constraint = graph_constraint_free(2),
#'     trial_success = gain,
#'     num_threads = 2
#' )
#' }
graph_optimise_gsd <- function(
    pvals,
    graph_constraint,
    trial_success,
    ...,
    alpha = NULL,
    start_graph = list(
        list(
            hyp_weight = NULL,
            trans_matrix = NULL
        )
    ),
    global_search = TRUE,
    num_threads = 1L,
    control = multigrain_control(),
    verbose = multigrain_verbosity()
) {
    check_pvals_gsd(pvals)
    check_graph_constraint(graph_constraint)
    check_trial_success_gsd(trial_success)
    rlang::check_dots_empty()
    rlang::check_number_whole(num_threads, min = 1)
    check_control(control)
    check_logical(global_search, allow_na = FALSE)

    if (isTRUE(verbose)) {
        verbose <- "info"
    } else if (isFALSE(verbose)) {
        verbose <- "silent"
    } else {
        rlang::check_string(verbose)
        verbose <- rlang::arg_match(verbose, values = verbosity_levels)
    }

    alpha <- .gsd_check_alpha(alpha, pvals)

    # Both dimensions of the gain are settled before the first evaluation: a
    # gain carrying discount tables indexes a C array with the decision time
    # (design record, section 10 item 6).
    .gsd_check_gain_dims(trial_success, pvals)

    .validate_start_graphs(
        start_graph,
        m = pvals$m
    )

    control <- control_prepare_dims(
        control,
        nsim = pvals$nsim,
        m = pvals$m,
        verbose = verbose
    )

    ga_result <- NULL
    x0_for_local <- NULL

    if (global_search) {
        ga_result <- .graph_optimise_ga_gsd(
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

    local_result <- .graph_optimise_local_gsd(
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

    clean_graph <- prune_graph_gsd(
        pvals = pvals,
        hyp_weight = best_graph_result$hyp_weight,
        trans_matrix = best_graph_result$trans_matrix,
        alpha = alpha,
        graph_constraint = graph_constraint,
        trial_success = trial_success,
        gamma = 1,
        verbose = verbose
    )

    if (verbose != "silent") {
        cli::cli_progress_step(
            "Evaluating trial success of pruned graph"
        )
    }

    final_power <- calc_power_pvals_gsd(
        pvals = pvals,
        hyp_weight = clean_graph$hyp_weight,
        trans_matrix = clean_graph$trans_matrix,
        alpha = alpha,
        custom_power = list(trial_success = trial_success)
    )

    # Apply hypothesis names
    gc_names <- graph_constraint_get_names(graph_constraint)
    names(clean_graph$hyp_weight) <- gc_names
    dimnames(clean_graph$trans_matrix) <- list(gc_names, gc_names)

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
#' @rdname graph_optimise_gsd
#' @usage NULL
graph_optimize_gsd <- graph_optimise_gsd


#' Prepare the control object from dimensions rather than a matrix
#'
#' `control_prepare()` calibrates the control object from the number of rows
#' and columns of a p-value matrix. A `multigrain_pvals_gsd` object is a
#' three-dimensional array wrapped in an S3 class, so the group sequential
#' optimiser passes the two dimensions directly. Everything else, including
#' the defaults injected and the `verbose` handling, is that of
#' `control_prepare()`.
#'
#' @inheritParams control_nsim_local
#' @param nsim (whole number) Number of simulated trials.
#' @param m (whole number) Number of hypotheses.
#' @inheritParams graph_optimise_gsd
#' @inheritParams rlang::args_error_context
#'
#' @returns A modified [multigrain_control()].
#'
#' @noRd
control_prepare_dims <- function(
    ctrl,
    nsim,
    m,
    verbose = "silent",
    call = rlang::caller_env()
) {
    check_control(ctrl, call = call)
    rlang::check_number_whole(nsim, min = 1, call = call)
    rlang::check_number_whole(m, min = 1, call = call)
    rlang::check_string(verbose, call = call)

    # `control_prepare()` takes both dimensions from `dim()`, which is always
    # integer; coerce so that the two paths agree bit for bit.
    nsim <- as.integer(nsim)
    m <- as.integer(m)

    default_ctrl <- default_control()

    # calibrate the default number of simulations with the dimensions
    default_ctrl$nsim_local <- nsim
    default_ctrl$nsim_global <- min(5e4L, nsim)

    default_ctrl$global_opt$popSize <- min(max(40L * m, 200L), 500L)

    # calibrate user-supplied nsim values
    ctrl <- adjust_nsim_local(ctrl, nsim, call = call)
    ctrl <- adjust_nsim_global(ctrl, nsim, call = call)

    if (verbose == "detail") {
        # adjusting the default verbosity allows the user-set values to pass
        default_ctrl$local_opt$print_level <- 1L
        default_ctrl$global_opt$monitor <- TRUE
    }

    ctrl$nsim_local <- ctrl$nsim_local %||% default_ctrl$nsim_local
    ctrl$nsim_global <- ctrl$nsim_global %||% default_ctrl$nsim_global

    ctrl$global_opt <- purrr::list_modify(
        default_ctrl$global_opt,
        !!!ctrl$global_opt
    )

    ctrl$local_opt <- purrr::list_modify(
        default_ctrl$local_opt,
        !!!ctrl$local_opt
    )

    ctrl
}


# Subsample simulated trials, keeping the `multigrain_pvals_gsd` object and
# its three-dimensional array (a `drop = TRUE` slip would silently turn an
# nsim by m by 1 array into a matrix; design record, section 8 item 12).
.sample_pvals_gsd <- function(pvals, nsim) {
    idx <- sample.int(pvals$nsim, size = nsim, replace = FALSE)

    pvals$pvals <- pvals$pvals[idx, , , drop = FALSE]
    pvals$nsim <- nsim

    pvals
}


#' Optimise group sequential graphs using [GA::ga()]
#'
#' Group sequential twin of `.graph_optimise_ga()`. The fitness is the
#' compiled gain applied to the decision times returned by the group
#' sequential kernel; the chosen graph is then re-evaluated on the full sample.
#'
#' @inheritParams graph_optimise_gsd
#'
#' @noRd
.graph_optimise_ga_gsd <- function(
    pvals,
    graph_constraint,
    trial_success,
    nsim, # defaults to control$nsim_global
    global_opts, # defaults to control$global_opt
    alpha = 0.025,
    num_threads = 1L,
    start_graph = list(
        list(
            hyp_weight = NULL,
            trans_matrix = NULL
        )
    ),
    verbose = c("info", "detail", "silent")
) {
    verbose <- rlang::arg_match(verbose)

    if (verbose != "silent") {
        cli::cli_progress_step("Running global optimization")
    }

    x0 <- .build_start_matrix(graph_constraint, start_graph)

    pvals_sampled <- .sample_pvals_gsd(pvals, nsim)

    immutable_global_args <- list(
        type = "real-valued",
        fitness = create_obj_func_gsd(
            m = trial_success$m,
            power_criterion = trial_success$func,
            hyp_constraint = graph_constraint$hyp_constraint,
            trans_constraint = graph_constraint$trans_constraint,
            alpha = alpha,
            pvals = .gsd_kernel_matrix(pvals_sampled),
            K = pvals$K,
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
        if (verbose != "silent") {
            cli::cli_progress_step(
                "Evaluating trial success of globally optimised graph"
            )
        }
        ga_subset_power <- calc_power_pvals_gsd(
            pvals_sampled, # calculated using nsim_global samples
            hyp_weight = sol$hyp_weight,
            trans_matrix = sol$trans_matrix,
            alpha = alpha,
            custom_power = list(trial_success = trial_success)
        )

        time_mat <- graph_shortcut_gsd(
            pvals = .gsd_kernel_matrix(pvals), # calculated on full sample
            alpha = alpha,
            w = sol$hyp_weight,
            G = sol$trans_matrix,
            K = pvals$K
        )$time

        ga_trial_success <- trial_success$func(time_mat)
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


#' Optimise group sequential graphs using [nloptr::nloptr()]
#'
#' Group sequential twin of `.graph_optimise_local()`.
#'
#' @inheritParams graph_optimise_gsd
#'
#' @noRd
.graph_optimise_local_gsd <- function(
    pvals,
    graph_constraint,
    trial_success,
    local_opts,
    alpha = 0.025,
    num_threads = 1L,
    nsim = pvals$nsim,
    x0 = NULL,
    verbose = c("info", "detail", "silent")
) {
    verbose <- rlang::arg_match(verbose)

    if (verbose != "silent") {
        cli::cli_progress_step("Running local optimization")
    }

    if (is.null(x0)) {
        x0 <- create_start_params(graph_constraint)
    }

    pvals_sampled <- .sample_pvals_gsd(pvals, nsim)

    obj_fun <- create_obj_func_gsd(
        m = trial_success$m,
        power_criterion = trial_success$func,
        hyp_constraint = graph_constraint$hyp_constraint,
        trans_constraint = graph_constraint$trans_constraint,
        alpha = alpha,
        pvals = .gsd_kernel_matrix(pvals_sampled),
        K = pvals$K,
        num_threads = num_threads
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
        if (verbose != "silent") {
            cli::cli_progress_step(
                "Evaluating trial success of locally optimised graph"
            )
        }
        local_subset_power <- calc_power_pvals_gsd(
            pvals_sampled, # calculated using nsim_local samples
            hyp_weight = sol$hyp_weight,
            trans_matrix = sol$trans_matrix,
            alpha = alpha,
            custom_power = list(trial_success = trial_success)
        )

        time_mat <- graph_shortcut_gsd(
            pvals = .gsd_kernel_matrix(pvals), # calculated on full sample
            alpha = alpha,
            w = sol$hyp_weight,
            G = sol$trans_matrix,
            K = pvals$K
        )$time

        local_trial_success <- trial_success$func(time_mat)
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
