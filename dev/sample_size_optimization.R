#' Create objective function used for sample size optimisation when a target
#' trial success measure is known and specified
#'
#' Power objective is embedded via `power_criterion`.
#' Returns fitness = (N_max - N + ts) unless constraints violated,
#' in which case returns a penalized negative value.
#'
#' @param power_criterion function; trial_success()$func
#' @param N_max scalar; upper bound on N
#' @param trial_success_target scalar or NULL; required trial success (e.g. 0.7)
#' @param power_constraint optional numeric vector of length m with per-hyp
#'        marginal power targets (NA for "no target" on that hypothesis)
#' @return function(x, alpha, z_matrix, mu_standard, hyp_constraint, trans_constraint)
#'
#' @noRd
create_fitness_min_N <- function(
    power_criterion,
    N_max,
    trial_success_target = NULL,
    power_constraint = NULL
) {
    GRAPH_PENALTY_HARD <- 100
    # GRAPH_PENALTY_SOFT <- 10
    POWER_PENALTY <- 10
    SOFT_FLOOR <- 0.75 #  nolint
    PENALTY_SCALE <- 2 * N_max

    stopifnot(length(N_max) == 1L, is.finite(N_max), N_max > 0)
    if (!is.null(trial_success_target)) {
        stopifnot(
            length(trial_success_target) == 1L,
            is.finite(trial_success_target)
        )
    }

    force(power_criterion)
    require_ts <- !is.null(trial_success_target)

    if (is.null(power_constraint)) {
        power_idx <- integer()
        target_power <- numeric()
    } else {
        power_idx <- which(!is.na(power_constraint))
        target_power <- power_constraint[power_idx]
    }
    n_target <- length(power_idx)

    function(
        x,
        alpha,
        z_matrix,
        mu_standard,
        hyp_constraint,
        trans_constraint
    ) {
        theta <- split_theta_minN(x, hyp_constraint, trans_constraint)
        N <- theta$N
        w <- recover_full_weights(theta$w_pars, hyp_constraint)
        G <- recover_full_trans_matrix(theta$g_pars, trans_constraint)

        graph_violation <- graph_violation_score_cpp(w, G)
        if (graph_violation > 0) {
            return(-graph_violation * GRAPH_PENALTY_HARD)
        }

        if (any(w < 1e-4)) {
            w[w < 1e-4] <- 0
        }
        if (any(G < 1e-5)) {
            G[G < 1e-5] <- 0
        }

        # vt <- graph_violation_tiered_cpp(w, G, sum_floor = SOFT_FLOOR)
        # if (vt$hard > 0) return(-as.numeric(vt$hard) - GRAPH_PENALTY_HARD)
        # # if (vt$soft > 0) return(-as.numeric(vt$soft) - GRAPH_PENALTY_SOFT)

        if (any(w < 1e-4)) {
            w[w < 1e-4] <- 0
        }
        if (any(G < 1e-5)) {
            G[G < 1e-5] <- 0
        }
        w_norm <- normalise_sum(w)
        G_norm <- t(apply(G, 1, function(row) {
            s <- sum(row)
            if (s <= 0) row else normalise_sum(row / s)
        }))

        weighting_strategy <- calc_local_weights(w_norm, G_norm)

        ncp <- mu_standard * sqrt(N)
        obs_z <- sweep(z_matrix, 2, ncp, "+")
        rejMatrix <- apply_ctp_z(obs_z, alpha, weighting_strategy)

        # Power-side penalties (only if graph passes tiered checks)
        penalty <- 0

        ts_val <- power_criterion(rejMatrix)
        if (require_ts && ts_val < trial_success_target) {
            penalty <- penalty + (trial_success_target - ts_val)
        }

        if (n_target > 0L) {
            marg <- colMeans(rejMatrix[, power_idx, drop = FALSE])
            short <- pmax(target_power - marg, 0)
            if (any(short > 0)) penalty <- penalty + sum(short)
        }

        if (penalty > 0) {
            return(N - PENALTY_SCALE * penalty - POWER_PENALTY)
        }

        N_max - N + ts_val
    }
}


#' Split encoded parameter vector into N, hypothesis weights, and G params
#' to be used with `recover_full_weights()` and `recover_full_trans_matrix()`
#'
#' @param theta (numeric) Stacked sample-size-optimisation vector.
#' @param hyp_constraint (numeric) Vector; NA marks free weights, otherwise fixed.
#' @param trans_constraint (numeric) Matrix; NA marks free edges, otherwise fixed.
#'
#' @return list(N = numeric(1), w_pars = numeric(), g_pars = numeric()).
#' @noRd
split_theta_minN <- function(theta, hyp_constraint, trans_constraint) {
    stopifnot(is.numeric(theta), length(theta) >= 1L)

    N <- theta[1]
    rest <- theta[-1L]

    free_w <- sum(is.na(hyp_constraint))
    w_len <- max(free_w - 1L, 0L)

    if (w_len) {
        w_pars <- rest[seq_len(w_len)]
        g_pars <- rest[-seq_len(w_len)]
    } else {
        w_pars <- numeric(0)
        g_pars <- rest
    }

    list(N = N, w_pars = w_pars, g_pars = g_pars)
}


#' Global optimisation of graphical test to *minimize N* subject to trial success
#'
#' Uses GA global search only (no local refinement) to find the smallest N (up to N_max)
#' that meets a required trial-success target and optional marginal power constraints.
#'
#' @inheritParams optimise_N
#'
#' @noRd
.optimise_N_ga <- function(
    ncp_vector,
    sigma,
    graph_constraint,
    trial_success,
    N_max,
    alpha = 0.025,
    trial_success_target = NULL,
    power_constraint = NULL,
    nsim_global = 5e4,
    ga_options = ga_N_options(),
    start_graph = list(
        list(
            hyp_weight = NULL,
            trans_matrix = NULL
        )
    ),
    verbose = TRUE
) {
    if (!verbose) {
        ga_options$monitor <- FALSE
    }

    m <- length(ncp_vector)
    stopifnot(
        is.matrix(sigma),
        nrow(sigma) == m,
        ncol(sigma) == m,
        trial_success$m <= m
    )

    fitness_func <- create_fitness_min_N(
        power_criterion = trial_success$func,
        trial_success_target = trial_success_target,
        N_max = N_max,
        power_constraint = power_constraint
    )

    # Quasi-random standard normal MVN samples under null
    z_null <- rqmvnorm_qr(nsim_global, rep(0, m), sigma = sigma)

    x0_graph <- .build_start_matrix(graph_constraint, start_graph)

    N_suggest <- pmax(2, floor(N_max / 2))
    suggestions <- cbind(N_suggest, x0_graph)

    default_ga <- list(
        type = "real-valued",
        fitness = fitness_func,
        lower = c(2, rep(0, ncol(x0_graph))),
        upper = c(N_max, rep(1, ncol(x0_graph))),
        population = GA::gaControl("real-valued")$population,
        mutation = GA::gaControl("real-valued")$mutation,
        pcrossover = 0.7,
        pmutation = 0.3,
        maxiter = 1e5,
        run = 200,
        monitor = verbose,
        optim = TRUE,
        optimArgs = list(method = "Nelder-Mead", poptim = 0.2, pressel = 0.6),
        hyp_constraint = graph_constraint$hyp_constraint,
        trans_constraint = graph_constraint$trans_constraint,
        alpha = alpha,
        z_matrix = z_null,
        mu_standard = ncp_vector,
        suggestions = suggestions
    )

    ga_options_renamed <- ga_options |>
        rlang::set_names(rename_ga_option) |>
        unclass()

    ga_args <- utils::modifyList(default_ga, ga_options_renamed)

    ga_res <- do.call(GA::ga, ga_args)
    best_raw <- ga_res@solution[1, ]

    sol <- param_to_solution_minN(best_raw, graph_constraint, process = TRUE)
    is_valid <- FALSE
    is_valid <- is_graph_valid(sol$hyp_weight, sol$trans_matrix)

    global_power <- NULL
    if (is_valid) {
        ncp <- ncp_vector * sqrt(sol$N)
        z <- sweep(z_null, 2, ncp, "+")
        pvals <- 1 - stats::pnorm(z)
        global_power <- calc_power_pvals(
            pvals = pvals,
            alpha = alpha,
            hyp_weight = sol$hyp_weight,
            trans_matrix = sol$trans_matrix,
            custom_power = list(trial_success = trial_success)
        )
    }

    obj <- list(
        N = sol$N,
        hyp_weight = sol$hyp_weight,
        trans_matrix = sol$trans_matrix,
        constraints = graph_constraint,
        trial_success = trial_success,
        power = global_power,
        opt_settings = list(
            global_search = TRUE,
            popSize = ga_args$popSize,
            max_generations = ga_args$maxiter,
            max_stagnation = ga_args$run,
            nsim_global = nsim_global,
            x0 = suggestions,
            elitism = ga_res@elitism,
            pcrossover = ga_res@pcrossover,
            pmutation = ga_res@pmutation,
            N_max = N_max
        ),
        global_output = list(
            iter = ga_res@iter,
            convergence = ga_res@summary,
            raw_solution = ga_res@solution,
            fitness = summary(ga_res)$fitness
        ),
        solution = list(
            opt_source = "GA_minN",
            graph_valid = c(global = is_valid)
        ),
        x0 = suggestions
    )

    class(obj) <- "graph_optimal"
    obj
}


param_to_solution_minN <- function(
    optimised_params,
    graph_constraint,
    process = FALSE
) {
    theta <- split_theta_minN(
        optimised_params,
        graph_constraint$hyp_constraint,
        graph_constraint$trans_constraint
    )

    w_sol <- recover_full_weights(theta$w_pars, graph_constraint$hyp_constraint)
    G_sol <- recover_full_trans_matrix(
        theta$g_pars,
        graph_constraint$trans_constraint
    )

    if (process) {
        w_sol[w_sol < 1e-4] <- 0
        w_sol[w_sol > 1e-4 & w_sol < 1e-3] <- 1e-3
        w_sol <- normalise_sum(w_sol)

        G_sol[G_sol < 1e-5] <- 0
        G_sol[G_sol > 1e-5 & G_sol < 1e-3] <- 1e-3
        G_sol <- t(apply(G_sol, 1, function(x) normalise_sum(x / sum(x))))
    }

    list(N = theta$N, hyp_weight = w_sol, trans_matrix = G_sol)
}


#' Optimise sample size \eqn{N} for graphical multiple testing under a
#' trial-success or marginal power target
#'
#' @description
#' `optimise_N()` searches for the **smallest** total sample size \eqn{N}
#' (the sum of each participants in each treatment allocation, up to a
#' user-defined `N_max`) and a valid graphical test - defined by hypothesis
#' weights \eqn{\mathbf w} and transition matrix \eqn{\mathbf G} - that achieves
#' a user-specified *trial success* target and (optionally) per-hypothesis
#' marginal power constraints, while preserving strong FWER control via the
#' graphical test.
#'
#' The routine performs a two-stage optimisation:
#'
#' 1. a **global** [GA::ga()] optimisation of
#'    \eqn{N}, \eqn{\mathbf w}, \eqn{\mathbf G}
#'    to meet the target with the smallest possible \eqn{N}, followed by
#'
#' 2. a **local** derivative-free improvement of the graph to yield the
#'    maximum trial success measure at the chosen \eqn{N}.
#'
#' @param ncp_vector Numeric vector of length \eqn{m} where \eqn{m} is the
#'   number of hypotheses under FWER-control. The per-hypothesis
#'   **standardised treatment effect per \eqn{\sqrt{N}}** is used to form the
#'   non-centrality parameters (NCPs) at sample size \eqn{N}:
#'   \deqn{\Delta_i(N) = \texttt{ncp\_vector}_i \sqrt{N}.}
#' @param sigma Numeric \eqn{m \times m} correlation matrix of the test
#'   statistics under the design scenario (positive semi-definite) where \eqn{m}
#'   is the number of hypotheses under FWER-control.
#' @param graph_constraint A `graph_constraint` object containing constraints on
#'  the graph's weights and transition matrix. See [graph_constraint()].
#' @param trial_success A `trial_success` object defining the trial success
#'   measure (power objective) function to maximize. See [trial_success()].
#' @param N_max Positive scalar. Upper bound for sample size in the global GA
#'   search.
#' @param alpha Numeric scalar. One-sided familywise significance level;
#'   default `0.025`.
#' @param trial_success_target Optional numeric scalar. If
#'   supplied, the GA penalises solutions that do not achieve at least this
#'   trial success measure value.
#' @param power_constraint Optional numeric vector of length \eqn{m} with
#'   **per-hypothesis marginal power targets** in \eqn{[0,1]}; use `NA` to
#'   indicate no constraint for a hypothesis. Shortfalls are penalised in the
#'   GA fitness.
#' @param nsim_global Integer. Number of simulations used to evaluate the
#'   objective during the global (GA) phase. Default `5e4`.
#' @param nsim_local Integer. Number of simulations used to evaluate the
#'   objective during the local refinement phase at \eqn{N}. Default `1e5`.
#' @param nlopt_options (`nlopt_options`) a list of options to be passed to
#'   [nloptr::nloptr()]. Must be the result of a call to `nlopt_options()`.
#' @param ga_options (`ga_options`) a list of options to be passed to
#'   [GA::ga()]. This must be the result of a call to [ga_N_options()].
#' @param start_graph Optional. Initial list of graphs suggested. Each graph is
#'   defined as a list containing starting weight vector `w` and transition
#'   matrix `G`. If `NULL`, default starting values are generated (currently
#'   with a Bonferroni-Holm graph).
#' @param verbose Logical. If `TRUE` (default), prints progress messages from
#'   the GA and local optimisation.
#'
#' @details
#' Let \eqn{m} be the number of elementary hypotheses. We assume the vector of
#' test statistics is approximately multivariate normal with unit variances and
#' correlation `sigma`, and mean vector
#' \eqn{\boldsymbol\Delta(N) = \texttt{ncp\_vector} \cdot \sqrt{N}}.
#' The graphical test is evaluated via closed testing with weighted Bonferroni
#' intersection tests derived from \eqn{\mathbf w} and \eqn{\mathbf G}; thus,
#' any valid graph produced maintains **strong FWER control** at level
#' \eqn{\alpha}.
#'
#' The global GA maximises a fitness of the form
#' \deqn{\text{fitness} = N_{\max} - N + \text{trial-success} - \text{penalties},}
#' hence favouring smaller \eqn{N} subject to meeting `trial_success_target`
#' and any `power_constraint`. After selecting \eqn{N}, a local
#' derivative-free optimiser refines \eqn{\mathbf w,\mathbf G} at that
#' \eqn{N} to maximise `trial_success`.
#'
#' @return
#' An object of class `graph_optimal` with the following elements:
#' \itemize{
#'   \item `N` - the selected sample size \eqn{N}.
#'   \item `hyp_weight` - optimised hypothesis weights \eqn{\mathbf w}.
#'   \item `trans_matrix` - optimised transition matrix \eqn{\mathbf G}.
#'   \item `constraints` - the `graph_constraint` used.
#'   \item `trial_success` - the trial-success definition supplied.
#'   \item `power` - summary of objective evaluated at the chosen solution
#'         (same structure as in [optimise_graph()]).
#'   \item `global_opt_power` - power/trial-success metrics for the GA solution
#'         evaluated on the full Monte Carlo set.
#'   \item `local_opt_power` - power/trial-success metrics for the locally
#'         refined solution.
#'   \item `solution` - list with `opt_source` (`"GA_minN"` or `"nlopt"`) and
#'         `graph_valid` flags for global/local solutions.
#'   \item `global_output`, `local_output` - raw outputs from GA / NLopt.
#'   \item `opt_settings` - a collection of settings actually used (including
#'         `nsim_global`, `nsim_local`, `N_max`, GA and NLopt
#'         parameters).
#'   \item `x0` - starting parameter vector used for the local optimisation.
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' alpha <- 0.025
#' m <- 4
#' # Example: per-sqrt(N) NCPs (from design-standardised effects)
#' mu_per_sqrtN <- c(0.10, 0.08, 0.07, 0.06)
#' Sigma <- matrix(c(
#'   1,   0.5, 0.45, 0.225,
#'   0.5, 1,   0.225,0.45,
#'   0.45,0.225,1,   0.5,
#'   0.225,0.45,0.5, 1
#' ), nrow = m, byrow = TRUE)
#'
#' # Composite trial-success: (H1 & H3) OR (H2 & H4)
#' ts <- trial_success((r1 && r3) || (r2 && r4))
#'
#' # Constrain weights/edges (example)
#' gc <- graph_constraint(
#'   hyp_constraint = c(1, 0, 0, 0),
#'   trans_constraint = rbind(
#'     c(0, NA, NA, 0),
#'     c(0, 0,  NA, NA),
#'     c(0, NA, 0,  NA),
#'     c(0, 0,  1,  0)
#'   )
#' )
#'
#' res <- optimise_N(
#'   ncp_vector           = mu_per_sqrtN,
#'   sigma                = Sigma,
#'   graph_constraint     = gc,
#'   trial_success        = ts,
#'   alpha                = alpha,
#'   N_max                = 1500,
#'   trial_success_target = 0.30,
#'   power_constraint     = c(NA, 0.70, NA, NA),
#'   nsim_global          = 1e4,
#'   nsim_local           = 1e4,
#'   verbose              = TRUE
#' )
#'
#' res$N
#' res$hyp_weight
#' res$trans_matrix
#' res$power$trial_success
#' }
#'
#' @export
#' @aliases optimize_N
optimise_N <- function(
    ncp_vector,
    sigma,
    graph_constraint,
    trial_success,
    N_max,
    alpha = 0.025,
    trial_success_target = NULL,
    power_constraint = NULL,
    nsim_global = 5e4,
    nsim_local = 1e5,
    nlopt_options = nlopt_options(),
    ga_options = ga_N_options(),
    start_graph = list(
        list(
            hyp_weight = NULL,
            trans_matrix = NULL
        )
    ),
    verbose = FALSE
) {
    # this pattern is needed to avoid a recursive default argument reference error
    # basically we can't have the name of the argument `nlopt_options` and the
    # name of the function giving its default value `nlopt_options()`
    # "promise already under evaluation: recursive default argument reference or
    # earlier problems?"
    # luckily rlang quosures allow us to capture the promise and force its
    # evaluation (we can't do this in base R)
    nlopt_options <- rlang::enquo(nlopt_options)
    nlopt_options <- rlang::eval_tidy(nlopt_options)

    assert_ga_options(ga_options)
    assert_nlopt_options(nlopt_options)

    warn_ga_options_type(
        ga_options,
        expected_type = "N"
    )

    if (verbose) {
        nlopt_options$print_level <- 1
        ga_options$monitor <- TRUE
    }

    m <- length(ncp_vector)
    .validate_start_graphs(
        start_graph = start_graph,
        m = m,
        where = "optimise_N"
    )

    cli::cli_progress_step("Running global optimisation")

    ga_obj <- .optimise_N_ga(
        ncp_vector = ncp_vector,
        sigma = sigma,
        graph_constraint = graph_constraint,
        trial_success = trial_success,
        alpha = alpha,
        trial_success_target = trial_success_target,
        N_max = N_max,
        power_constraint = power_constraint,
        nsim_global = nsim_global,
        ga_options = ga_options,
        start_graph = start_graph,
        verbose = verbose
    )

    N_star <- ga_obj$N

    if (N_star > N_max) {
        # nolint
        warning(
            paste("Could not find graph with sample size less than N =", N_max),
            call. = FALSE
        )
    }

    z_local <- rqmvnorm_qr(
        nsim_local,
        rep(0, length(ncp_vector)),
        sigma = sigma
    )
    z_local <- sweep(z_local, 2, ncp_vector * sqrt(N_star), "+")
    pvals_local <- 1 - stats::pnorm(z_local)

    x0_for_local <- create_start_params(
        graph_constraint,
        w0 = as.vector(ga_obj$hyp_weight),
        g0 = ga_obj$trans_matrix
    )

    cli::cli_progress_step("Running local optimisation")

    local_graph <- .optimise_graph_local(
        pvals = pvals_local,
        graph_constraint = graph_constraint,
        trial_success = trial_success,
        alpha = alpha,
        nsim = nsim_local,
        nlopt_options = nlopt_options,
        x0 = x0_for_local
    )

    cli::cli_progress_step(
        "Evaluating trial success of globally optimised graph"
    )

    all_pval_power_ga <- calc_power_pvals(
        pvals = pvals_local,
        alpha = alpha,
        hyp_weight = ga_obj$hyp_weight,
        trans_matrix = ga_obj$trans_matrix,
        custom_power = list(trial_success = trial_success)
    )

    cli::cli_progress_step(
        "Evaluating trial success of locally optimised graph"
    )

    all_pval_power_local <- calc_power_pvals(
        pvals = pvals_local,
        alpha = alpha,
        hyp_weight = local_graph$hyp_weight,
        trans_matrix = local_graph$trans_matrix,
        custom_power = list(trial_success = trial_success)
    )

    ga_valid <- isTRUE(ga_obj$solution$graph_valid["global"])
    local_valid <- isTRUE(local_graph$solution$graph_valid["local"])

    pick_ga <- ga_valid &&
        (!local_valid ||
            (local_valid &&
                all_pval_power_ga$trial_success >
                    all_pval_power_local$trial_success))

    chosen <- if (pick_ga) ga_obj else local_graph
    other <- if (pick_ga) local_graph else ga_obj

    cli::cli_progress_step("Pruning redundant weights")

    pruned <- prune_graph(
        pvals = pvals_local,
        hyp_weight = chosen$hyp_weight,
        trans_matrix = chosen$trans_matrix,
        alpha = alpha,
        graph_constraint = graph_constraint,
        trial_success = trial_success,
        gamma = 1
    )

    chosen$hyp_weight <- pruned$hyp_weight
    names(chosen$hyp_weight) <- names(graph_constraint$hyp_constraint)

    chosen$trans_matrix <- pruned$trans_matrix
    dimnames(chosen$trans_matrix) <- list(
        names(graph_constraint$hyp_constraint),
        names(graph_constraint$hyp_constraint)
    )

    chosen$solution$graph_valid <- c(
        chosen$solution$graph_valid,
        other$solution$graph_valid
    )

    chosen$global_opt_power <- all_pval_power_ga
    chosen$local_opt_power <- all_pval_power_local
    chosen$power <- if (pick_ga) {
        chosen$global_opt_power
    } else {
        chosen$local_opt_power
    }

    if (is.null(chosen$global_output) && !is.null(ga_obj$global_output)) {
        chosen$global_output <- ga_obj$global_output
    }

    if (is.null(chosen$local_output) && !is.null(local_graph$local_output)) {
        chosen$local_output <- local_graph$local_output
    }

    chosen$opt_settings$global_search <- TRUE
    chosen$opt_settings$nsim_local <- nsim_local
    chosen$opt_settings$nsim_global <- nsim_global
    chosen$opt_settings$N_max <- N_max

    if (is.null(chosen$x0)) {
        chosen$x0 <- x0_for_local
    }
    chosen$N <- N_star

    class(chosen) <- "graph_optimal"

    chosen
}

#' @export
#' @rdname optimise_N
#' @usage NULL
optimize_N <- optimise_N
