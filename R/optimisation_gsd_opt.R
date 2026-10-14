# The driver of the joint search over the graph and the optimised
# Hwang-Shih-DeCani spending parameters (dev/hsd_opt_design_revised.md,
# sections 1.3 to 1.5 and 5.2; build plan, phase P3). `graph_optimise_gsd()`
# hands a `multigrain_pvals_gsd_spec` to `.graph_optimise_gsd_opt()` as its
# first statement; nothing here is reached on the fixed-spending path.
#
# Stage G is a genetic algorithm on the slice bank with one gene per
# parameter appended to the graph genes; Stage L refines the graph by COBYLA
# at the GA's parameter values and then sweeps every lattice node of each
# parameter, re-optimising the graph at each node from three candidate starts
# (the warm-started profile). All gains of Stage L are on one pinned row
# order, `idx_local`, so they are paired by simulated trial.

# Driver ------------------------------------------------------------------

#' Optimise a graph and its spending parameters
#'
#' The body of [graph_optimise_gsd()] when `pvals` is a
#' `multigrain_pvals_gsd_spec`. Validation mirrors the fixed path on the spec
#' (whose `nsim`, `m`, `K` and `alpha` fields are those of the p-values it
#' resolves to); the default population size counts each spending parameter
#' as one more hypothesis (build plan, reconciliation 5).
#'
#' `choose_graph()` is reused unchanged: the GA winner is re-scored on the
#' full sample at its own parameter values and Stage L's winner at its own,
#' so the comparison is on gains. When the GA wins, the reported spending
#' parameters are the GA's, and the reported profile is still Stage L's,
#' which started there (reconciliation 11).
#'
#' @inheritParams graph_optimise_gsd
#' @param spec A `multigrain_pvals_gsd_spec` object.
#' @param call The error context, the frame of [graph_optimise_gsd()].
#'
#' @returns A `multigrain_graph_optimal` object with the fields of the fixed
#'   path plus `spending_param`, `hyp_spending_param`, `spending`,
#'   `spending_bounds`, `spending_profile`, `spending_search` and
#'   `profile_converged`.
#' @noRd
.graph_optimise_gsd_opt <- function(
    spec,
    graph_constraint,
    trial_success,
    alpha = NULL,
    stopping_rule = NULL,
    start_graph = list(
        list(
            hyp_weight = NULL,
            trans_matrix = NULL
        )
    ),
    global_search = TRUE,
    num_threads = 1L,
    control = multigrain_control(),
    verbose = multigrain_verbosity(),
    call = rlang::caller_env()
) {
    check_graph_constraint(graph_constraint, call = call)
    check_trial_success_gsd(trial_success, call = call)
    check_stopping_rule_gsd(stopping_rule, allow_null = TRUE, call = call)
    rlang::check_number_whole(num_threads, min = 1, call = call)
    check_control(control, call = call)
    check_logical(global_search, allow_na = FALSE, call = call)

    if (isTRUE(verbose)) {
        verbose <- "info"
    } else if (isFALSE(verbose)) {
        verbose <- "silent"
    } else {
        rlang::check_string(verbose, call = call)
        verbose <- rlang::arg_match(
            verbose,
            values = verbosity_levels,
            error_call = call
        )
    }

    alpha <- .gsd_check_alpha(alpha, spec, call = call)

    .gsd_check_graph_dims(
        graph_constraint_get_m(graph_constraint),
        spec,
        "{.arg graph_constraint}",
        call = call
    )
    .gsd_check_gain_dims(trial_success, spec, call = call)

    .validate_start_graphs(start_graph, m = spec$m, call = call)

    stop_rule <- if (is.null(stopping_rule)) {
        NULL
    } else {
        .gsd_stop_table(
            .gsd_check_stopping_rule_dims(stopping_rule, spec, call = call),
            spec$m,
            spec$K,
            call = call
        )
    }

    d <- length(spec$param$names)
    p_graph <- ncol(.build_start_matrix(graph_constraint, start_graph))

    user_pop <- control$global_opt$popSize
    control <- control_prepare_dims(
        control,
        nsim = spec$nsim,
        m = spec$m,
        verbose = verbose,
        call = call
    )
    if (is.null(user_pop)) {
        control$global_opt$popSize <- min(max(40L * (spec$m + d), 200L), 500L)
    }

    ga_result <- NULL

    if (global_search) {
        ga_result <- .graph_optimise_ga_gsd_opt(
            spec = spec,
            graph_constraint = graph_constraint,
            trial_success = trial_success,
            alpha = alpha,
            nsim = control$nsim_global,
            global_opts = control$global_opt,
            num_threads = num_threads,
            start_graph = start_graph,
            verbose = verbose,
            stop_rule = stop_rule,
            p_graph = p_graph
        )
        x0 <- pmin(
            pmax(ga_result$ga_output@solution[1L, seq_len(p_graph)], 0),
            1
        )
        theta0 <- ga_result$ga_spending_param
    } else {
        x0 <- .hsd_start_genes(graph_constraint, start_graph)
        theta0 <- spec$param$start
    }

    local_result <- .graph_optimise_local_gsd_opt(
        spec = spec,
        graph_constraint = graph_constraint,
        trial_success = trial_success,
        alpha = alpha,
        nsim = control$nsim_local,
        local_opts = control$local_opt,
        num_threads = num_threads,
        x0 = x0,
        theta0 = theta0,
        archive = ga_result$archive,
        verbose = verbose,
        stop_rule = stop_rule,
        p_graph = p_graph
    )

    best_graph_result <- choose_graph(ga_result, local_result)
    theta_star <- if (best_graph_result$source == "global") {
        ga_result$ga_spending_param
    } else {
        local_result$local_spending_param
    }
    theta_star <- stats::setNames(
        .hsd_resolve_nodes(spec$param, theta_star),
        spec$param$names
    )

    pv_star <- resolve_pvalues_gsd(spec, theta_star)

    clean_graph <- prune_graph_gsd(
        pvals = pv_star,
        hyp_weight = best_graph_result$hyp_weight,
        trans_matrix = best_graph_result$trans_matrix,
        alpha = alpha,
        graph_constraint = graph_constraint,
        trial_success = trial_success,
        gamma = 1,
        verbose = verbose,
        stop_rule = stop_rule
    )

    if (verbose != "silent") {
        cli::cli_progress_step(
            "Evaluating trial success of pruned graph"
        )
    }

    final_power <- calc_power_pvals_gsd(
        pvals = pv_star,
        hyp_weight = clean_graph$hyp_weight,
        trans_matrix = clean_graph$trans_matrix,
        alpha = alpha,
        custom_power = list(trial_success = trial_success),
        stopping_rule = stopping_rule
    )

    gc_names <- graph_constraint_get_names(graph_constraint)
    names(clean_graph$hyp_weight) <- gc_names
    dimnames(clean_graph$trans_matrix) <- list(gc_names, gc_names)

    hyp_spending_param <- stats::setNames(
        unname(theta_star)[spec$group],
        gc_names
    )

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
        start_graph = start_graph,
        stopping_rule = stopping_rule,
        spending = list(
            spending_param = theta_star,
            hyp_spending_param = hyp_spending_param,
            spending = pv_star$spending,
            spending_bounds = spec$param$bounds,
            spending_step = spec$param$step,
            spending_group = spec$group,
            spending_profile = local_result$spending_profile,
            spending_search = .hsd_search_table(spec, ga_result$archive),
            profile_converged = local_result$profile_converged
        )
    )
}

# Graph genes of the local stage's start when there is no global search: the
# user's first start graph, encoded as `.build_start_matrix()` encodes user
# rows (a missing weight vector or matrix is filled as it fills them), else
# the default start. Unlike the fixed path, a user start graph is honoured
# (build plan, reconciliation 9).
.hsd_start_genes <- function(graph_constraint, start_graph) {
    first <- if (length(start_graph)) start_graph[[1L]] else NULL
    if (is.null(first$hyp_weight) && is.null(first$trans_matrix)) {
        return(as.numeric(create_start_params(graph_constraint)))
    }

    m <- graph_constraint_get_m(graph_constraint)
    w0 <- if (is.null(first$hyp_weight)) {
        rep(1, m)
    } else {
        as.numeric(first$hyp_weight)
    }
    G0 <- if (is.null(first$trans_matrix)) {
        start_mat <- matrix(0, m, m)
        if (m > 1L) {
            start_mat[row(start_mat) != col(start_mat)] <- 1 / (m - 1)
        }
        start_mat
    } else {
        as.matrix(first$trans_matrix)
    }

    .hsd_encode_graph(graph_constraint, w0, G0)
}

# The `spending_search` table: per parameter and node, the GA's visits and
# best fitness on its subsample (never comparable with Stage L's gains, hence
# the column name). All `NA` without a global search.
.hsd_search_table <- function(spec, archive) {
    param <- spec$param
    pieces <- lapply(seq_along(param$names), function(r) {
        lattice <- param$lattice[[r]]
        n_nodes <- length(lattice)
        if (is.null(archive)) {
            visits <- rep(NA_integer_, n_nodes)
            best <- rep(NA_real_, n_nodes)
        } else {
            visits <- archive$visits[[r]]
            best <- archive$best[[r]]
            best[!is.finite(best)] <- NA_real_
        }
        data.frame(
            param = rep(param$names[[r]], n_nodes),
            value = lattice,
            visits = visits,
            best_fitness_subsample = best
        )
    })

    out <- do.call(rbind, pieces)
    rownames(out) <- NULL
    out
}


# Stage G -----------------------------------------------------------------

#' Stage G: the genetic algorithm on the slice bank
#'
#' Twin of `.graph_optimise_ga_gsd()` with one gene per spending parameter
#' appended to the graph genes. The random number stream is consumed in the
#' same order: the start matrix (no draws), the subsample `idx_ga` (the call
#' `.sample_pvals_gsd()` makes), then `GA::ga()`. The GA arguments are the
#' fixed driver's. The winner's parameter genes map to lattice nodes; its
#' graph genes are decoded as the fixed driver decodes them (no
#' `repair_graph()`), and a valid winner is re-scored on the full sample at
#' its own parameter values.
#'
#' @inheritParams .graph_optimise_gsd_opt
#' @param nsim Number of simulated trials of the GA subsample.
#' @param global_opts The `global_opt` list of the prepared control.
#' @param stop_rule The stopping-rule table, or `NULL`.
#' @param p_graph Number of graph genes.
#'
#' @returns The fixed driver's list plus `ga_spending_param` (named lattice
#'   values), `archive` (a plain list of `best`, `visits` and `genes` per
#'   parameter) and `idx_ga`.
#' @noRd
.graph_optimise_ga_gsd_opt <- function(
    spec,
    graph_constraint,
    trial_success,
    alpha,
    nsim,
    global_opts,
    num_threads,
    start_graph,
    verbose,
    stop_rule,
    p_graph
) {
    if (verbose != "silent") {
        cli::cli_progress_step("Running global optimization")
    }

    param <- spec$param
    d <- length(param$names)
    n_nodes <- lengths(param$lattice, use.names = FALSE)

    x0 <- .build_start_matrix(graph_constraint, start_graph)

    idx_ga <- sample.int(spec$nsim, size = nsim, replace = FALSE)

    fitness <- create_obj_func_gsd_opt(
        spec,
        bank = .hsd_bank(spec, idx_ga),
        power_criterion = trial_success$func,
        hyp_constraint = graph_constraint$hyp_constraint,
        trans_constraint = graph_constraint$trans_constraint,
        alpha = alpha,
        p_graph = p_graph,
        num_threads = num_threads,
        stop_rule = stop_rule
    )

    start_q <- vapply(
        seq_len(d),
        function(r) match(param$start[[r]], param$lattice[[r]]),
        integer(1L)
    )
    start_genes <- (start_q - 0.5) / n_nodes
    suggestions <- cbind(
        x0,
        matrix(start_genes, nrow = nrow(x0), ncol = d, byrow = TRUE)
    )

    immutable_global_args <- list(
        type = "real-valued",
        fitness = fitness,
        lower = rep(0, p_graph + d),
        upper = rep(1, p_graph + d),
        population = .cauchy_population,
        mutation = .make_cauchy_mutation_multi(
            p_param_mutate = 0.1,
            scale = 1.0
        ),
        optim = TRUE,
        suggestions = suggestions
    )

    ga_args <- utils::modifyList(immutable_global_args, global_opts)

    ga_res <- do.call(GA::ga, ga_args)

    # the archive leaves as a plain list. The closure is not called again:
    # its working matrix is dropped, and so is the `@call` slot, which holds
    # the closure and with it the bank (`graph_optimal()` drops that slot in
    # any case, so the result is unchanged)
    state <- attr(fitness, "archive")
    archive <- list(
        best = state$best,
        visits = state$visits,
        genes = state$genes
    )
    state$mat <- NULL
    attr(ga_res, "call") <- NULL
    rm(fitness, state)

    best_raw <- ga_res@solution[1L, ]
    node_q <- .hsd_node_index(best_raw[p_graph + seq_len(d)], n_nodes)
    theta_ga <- stats::setNames(
        vapply(
            seq_len(d),
            function(r) param$lattice[[r]][[node_q[[r]]]],
            double(1L)
        ),
        param$names
    )

    sol <- param_to_solution(
        best_raw[seq_len(p_graph)],
        graph_constraint,
        process = TRUE
    )

    is_valid <- is_graph_valid(sol$hyp_weight, sol$trans_matrix)

    ga_subset_power <- NULL
    ga_trial_success <- NULL

    if (is_valid) {
        if (verbose != "silent") {
            cli::cli_progress_step(
                "Evaluating trial success of globally optimised graph"
            )
        }
        pv_ga <- resolve_pvalues_gsd(spec, theta_ga)

        ga_subset_power <- .calc_power_pvals_gsd_core(
            .hsd_subset_pvals(pv_ga, idx_ga), # the GA's subsample
            hyp_weight = sol$hyp_weight,
            trans_matrix = sol$trans_matrix,
            alpha = alpha,
            custom_power = list(trial_success = trial_success),
            stop_rule = stop_rule
        )

        res <- graph_shortcut_gsd(
            pvals = .gsd_kernel_matrix(pv_ga), # full sample
            alpha = alpha,
            w = sol$hyp_weight,
            G = sol$trans_matrix,
            K = spec$K,
            stop_rule = stop_rule
        )

        ga_trial_success <- trial_success$func(res$time)
    }

    list(
        ga_hyp_weight = sol$hyp_weight,
        ga_trans_matrix = sol$trans_matrix,
        ga_trial_success = ga_trial_success,
        ga_subset_power = ga_subset_power,
        is_graph_valid = is_valid,
        ga_output = ga_res,
        ga_spending_param = theta_ga,
        archive = archive,
        idx_ga = idx_ga
    )
}


# Stage L -----------------------------------------------------------------

#' Stage L: graph refinement, then the warm-started profile
#'
#' Draws `idx_local` once (the call `.sample_pvals_gsd()` makes) and keeps
#' that row order for every evaluation below (build plan, reconciliation 6).
#' Step 1 is `.hsd_refine_node()` from `x0` at `theta0`: the fixed local
#' driver's `nloptr` call on the pinned rows. Step 2 sweeps each parameter's
#' lattice outward from the current node (`.hsd_profile_axis()`) and moves
#' to the axis's best node when it beats the current gain; for `d >= 2` the
#' sweep is repeated while any axis moved, at most three rounds, and an axis
#' whose profile no longer passes through the final optimum is recomputed
#' once at the end. Differences and standard errors are against the final
#' optimum; a recompute can expose a negative difference, which is kept and
#' reported as `profile_converged = FALSE`.
#'
#' @inheritParams .graph_optimise_ga_gsd_opt
#' @param local_opts The `local_opt` list of the prepared control.
#' @param x0 Graph genes of the start.
#' @param theta0 Lattice values of the start, one per parameter.
#' @param archive The GA archive (a plain list) or `NULL`.
#' @param density `1` profiles every node; `2` every second node from the
#'   current one, plus both bounds.
#'
#' @returns The fixed local driver's list plus `local_spending_param`,
#'   `spending_profile`, `profile_converged` and `idx_local`.
#' @noRd
.graph_optimise_local_gsd_opt <- function(
    spec,
    graph_constraint,
    trial_success,
    alpha,
    nsim,
    local_opts,
    num_threads,
    x0,
    theta0,
    archive,
    verbose,
    stop_rule,
    p_graph,
    density = 1L
) {
    if (verbose != "silent") {
        cli::cli_progress_step("Running local optimization")
    }

    stopifnot(length(x0) == p_graph)

    idx_local <- sample.int(spec$nsim, size = nsim, replace = FALSE)

    ctx <- .hsd_local_context(
        spec,
        graph_constraint = graph_constraint,
        trial_success = trial_success,
        alpha = alpha,
        num_threads = num_threads,
        stop_rule = stop_rule,
        local_opts = local_opts,
        idx = idx_local,
        archive = archive,
        density = density
    )

    theta0 <- stats::setNames(
        .hsd_resolve_nodes(spec$param, theta0),
        spec$param$names
    )

    # step 1: COBYLA on the graph genes at theta0
    kmat0 <- .gsd_kernel_matrix(.hsd_node_array(ctx, theta0))
    step1 <- .hsd_refine_node(
        x0,
        obj = .hsd_objective(ctx, kmat0),
        kmat = kmat0,
        ctx = ctx
    )
    rm(kmat0)

    current <- list(
        theta = theta0,
        hyp_weight = step1$hyp_weight,
        trans_matrix = step1$trans_matrix,
        gain = step1$gain,
        blocks = step1$blocks,
        version = 0L
    )

    # step 2: the profile sweep, in rounds for d >= 2
    if (verbose != "silent") {
        cli::cli_progress_step(paste(
            "Profiling spending parameters:",
            toString(spec$param$names)
        ))
    }
    profiled <- .hsd_profile(current, ctx)
    current <- profiled$current

    theta_star <- current$theta
    is_valid <- is_graph_valid(current$hyp_weight, current$trans_matrix)

    local_subset_power <- NULL
    local_trial_success <- NULL

    if (is_valid) {
        if (verbose != "silent") {
            cli::cli_progress_step(
                "Evaluating trial success of locally optimised graph"
            )
        }
        local_subset_power <- .calc_power_pvals_gsd_core(
            .hsd_node_array(ctx, theta_star), # the pinned rows
            hyp_weight = current$hyp_weight,
            trans_matrix = current$trans_matrix,
            alpha = alpha,
            custom_power = list(trial_success = trial_success),
            stop_rule = stop_rule
        )

        res <- graph_shortcut_gsd(
            pvals = .gsd_kernel_matrix(resolve_pvalues_gsd(spec, theta_star)),
            alpha = alpha,
            w = current$hyp_weight,
            G = current$trans_matrix,
            K = spec$K,
            stop_rule = stop_rule
        )

        local_trial_success <- trial_success$func(res$time)
    }

    list(
        local_hyp_weight = current$hyp_weight,
        local_trans_matrix = current$trans_matrix,
        local_trial_success = local_trial_success,
        local_subset_power = local_subset_power,
        is_graph_valid = is_valid,
        local_output = step1$nloptr_output,
        local_spending_param = theta_star,
        spending_profile = profiled$profile,
        profile_converged = profiled$converged,
        idx_local = idx_local
    )
}

# The axis rounds of step 2 and the final recompute of stale axes. An axis's
# profile is stale when the optimum changed (node or graph) after it was
# computed, so that its row at the final node may not carry the final graph.
.hsd_profile <- function(current, ctx) {
    param <- ctx$spec$param
    d <- length(param$names)
    max_rounds <- if (d == 1L) 1L else 3L

    rows <- vector("list", d)
    axis_version <- integer(d)
    moved_any <- FALSE

    for (round in seq_len(max_rounds)) {
        moved_any <- FALSE
        for (r in seq_len(d)) {
            # the version the axis was computed against: if its own sweep
            # improves the optimum (a better graph at the current node, or a
            # move), the axis is stale and is recomputed below against the
            # final optimum, so that S_r and the profile agree at theta*
            version_before <- current$version
            rows[[r]] <- .hsd_profile_axis(r, current, ctx)
            gains <- vapply(rows[[r]], `[[`, double(1L), "gain")
            best <- which.max(gains)
            if (gains[[best]] > current$gain) {
                top <- rows[[r]][[best]]
                moved_any <- moved_any || top$value != current$theta[[r]]
                current$theta[[r]] <- top$value
                current$hyp_weight <- top$hyp_weight
                current$trans_matrix <- top$trans_matrix
                current$gain <- top$gain
                current$blocks <- top$blocks
                current$version <- current$version + 1L
            }
            axis_version[[r]] <- version_before
        }
        if (!moved_any) {
            break
        }
    }

    converged <- d == 1L || !moved_any

    for (r in which(axis_version != current$version)) {
        rows[[r]] <- .hsd_profile_axis(r, current, ctx, final = TRUE)
    }

    profile_df <- .hsd_profile_table(rows, param$names, current$blocks)
    if (any(profile_df$diff < 0)) {
        converged <- FALSE
    }

    list(current = current, profile = profile_df, converged = converged)
}

# The profile rows of one axis, `r`: every node (or every second node and
# the bounds, `density = 2`) of its lattice with the other parameters held
# at `current$theta`, walked outward from the current node in both
# directions. Each row is `.hsd_profile_node()` warm-started from the
# previous row in its direction. Returned in lattice order.
.hsd_profile_axis <- function(r, current, ctx, final = FALSE) {
    lattice <- ctx$spec$param$lattice[[r]]
    n_nodes <- length(lattice)
    at <- match(current$theta[[r]], lattice)

    keep <- seq_len(n_nodes)
    if (identical(as.integer(ctx$density), 2L)) {
        keep <- keep[(keep - at) %% 2L == 0L | keep %in% c(1L, n_nodes)]
    }

    rows <- vector("list", n_nodes)
    # in the final pass the row at the current node is the optimum by
    # definition: no refinement there, so S_r, the profile and the
    # reported gain agree exactly at theta* (design 2.5)
    rows[[at]] <- .hsd_profile_node(
        r,
        at,
        current,
        prev = NULL,
        ctx = ctx,
        refine = !final
    )
    walks <- list(keep[keep > at], rev(keep[keep < at]))
    for (walk in walks) {
        prev <- rows[[at]]
        for (q in walk) {
            rows[[q]] <- .hsd_profile_node(r, q, current, prev = prev, ctx)
            prev <- rows[[q]]
        }
    }

    rows[keep]
}

# One profile row: node `q` of axis `r`. The candidates are the current
# optimum's graph (first, so it wins ties and the row at the current node
# carries it unless something beats it), the previous row's graph and the
# GA archive's best graph at this node; each is evaluated on this node's
# array, COBYLA starts from the best, and the row records the best of the
# candidates and the refined graph.
.hsd_profile_node <- function(r, q, current, prev, ctx, refine = TRUE) {
    value <- ctx$spec$param$lattice[[r]][[q]]
    theta_q <- current$theta
    theta_q[[r]] <- value

    kmat <- .gsd_kernel_matrix(.hsd_node_array(ctx, theta_q))

    candidates <- list(
        star = list(
            hyp_weight = current$hyp_weight,
            trans_matrix = current$trans_matrix
        )
    )
    if (!is.null(prev)) {
        candidates$prev <- prev[c("hyp_weight", "trans_matrix")]
    }
    archived <- .hsd_archive_graph(ctx, r, q)
    if (!is.null(archived)) {
        candidates$archive <- archived
    }

    evaluated <- lapply(candidates, function(g) {
        .hsd_eval_graph(kmat, g$hyp_weight, g$trans_matrix, ctx)
    })
    cand_gain <- vapply(evaluated, `[[`, double(1L), "gain")

    if (!refine) {
        return(list(
            value = value,
            gain = cand_gain[["star"]],
            gain_fixed_graph = cand_gain[["star"]],
            n_eval = 0L,
            hyp_weight = current$hyp_weight,
            trans_matrix = current$trans_matrix,
            blocks = .hsd_block_gains(evaluated$star$time, ctx)
        ))
    }

    from <- candidates[[which.max(cand_gain)]]

    refined <- .hsd_refine_node(
        .hsd_encode_graph(
            ctx$graph_constraint,
            from$hyp_weight,
            from$trans_matrix
        ),
        obj = .hsd_objective(ctx, kmat),
        kmat = kmat,
        ctx = ctx
    )

    graphs <- c(candidates, list(refined = refined))
    gains <- c(cand_gain, refined = refined$gain)
    pick <- which.max(gains)
    blocks <- if (names(gains)[[pick]] == "refined") {
        refined$blocks
    } else {
        .hsd_block_gains(evaluated[[pick]]$time, ctx)
    }

    list(
        value = value,
        gain = gains[[pick]],
        gain_fixed_graph = cand_gain[["star"]],
        n_eval = refined$n_eval,
        hyp_weight = graphs[[pick]]$hyp_weight,
        trans_matrix = graphs[[pick]]$trans_matrix,
        blocks = blocks
    )
}

# The GA archive's best graph at node `q` of axis `r`, decoded as the GA
# winner is decoded, or NULL for an unvisited node, an absent archive or an
# invalid decode.
.hsd_archive_graph <- function(ctx, r, q) {
    if (is.null(ctx$archive)) {
        return(NULL)
    }
    genes <- ctx$archive$genes[[r]][q, ]
    if (anyNA(genes)) {
        return(NULL)
    }
    sol <- param_to_solution(genes, ctx$graph_constraint, process = TRUE)
    if (!is_graph_valid(sol$hyp_weight, sol$trans_matrix)) {
        return(NULL)
    }

    sol
}

# The `spending_profile` data frame from the rows of every axis, with the
# paired differences and standard errors against the optimum's block gains
# and the per-row graphs as the attribute `graphs`.
.hsd_profile_table <- function(rows, param_names, star_blocks) {
    pieces <- lapply(seq_along(rows), function(r) {
        axis_rows <- rows[[r]]
        paired <- lapply(axis_rows, function(row) {
            .hsd_block_diff(star_blocks, row$blocks)
        })
        data.frame(
            param = rep(param_names[[r]], length(axis_rows)),
            value = vapply(axis_rows, `[[`, double(1L), "value"),
            gain = vapply(axis_rows, `[[`, double(1L), "gain"),
            gain_fixed_graph = vapply(
                axis_rows,
                `[[`,
                double(1L),
                "gain_fixed_graph"
            ),
            diff = vapply(paired, `[[`, double(1L), "diff"),
            se = vapply(paired, `[[`, double(1L), "se"),
            n_eval = vapply(axis_rows, `[[`, integer(1L), "n_eval")
        )
    })

    out <- do.call(rbind, pieces)
    rownames(out) <- NULL
    attr(out, "graphs") <- unlist(
        lapply(rows, function(axis_rows) {
            lapply(axis_rows, `[`, c("hyp_weight", "trans_matrix"))
        }),
        recursive = FALSE
    )

    out
}


# Stage L helpers ---------------------------------------------------------

# Everything a Stage L evaluation needs, in one list: the spec, the graph
# constraint, the gain, the level, the kernel options, the local options, the
# pinned rows `idx`, the GA archive (or NULL), the profile density and the
# number of batch-means blocks.
.hsd_local_context <- function(
    spec,
    graph_constraint,
    trial_success,
    alpha,
    local_opts,
    idx,
    num_threads = 1L,
    stop_rule = NULL,
    archive = NULL,
    density = 1L,
    n_block = 100L
) {
    list(
        spec = spec,
        graph_constraint = graph_constraint,
        trial_success = trial_success,
        alpha = alpha,
        num_threads = num_threads,
        stop_rule = stop_rule,
        local_opts = local_opts,
        idx = idx,
        archive = archive,
        density = density,
        n_block = n_block
    )
}

# `pvals` restricted to rows `idx`, in that order, keeping the class and the
# three-dimensional array (as `.sample_pvals_gsd()` does, without drawing).
.hsd_subset_pvals <- function(pvals, idx) {
    pvals$pvals <- pvals$pvals[idx, , , drop = FALSE]
    pvals$nsim <- length(idx)

    pvals
}

# The resolved `multigrain_pvals_gsd` at `theta` on the pinned rows.
.hsd_node_array <- function(ctx, theta) {
    .hsd_subset_pvals(resolve_pvalues_gsd(ctx$spec, theta), ctx$idx)
}

# The unchanged fixed-path closure on a kernel matrix, so that step 1 is the
# fixed local driver's optimisation (P3 oracle (a)).
.hsd_objective <- function(ctx, kmat) {
    create_obj_func_gsd(
        m = ctx$trial_success$m,
        power_criterion = ctx$trial_success$func,
        hyp_constraint = ctx$graph_constraint$hyp_constraint,
        trans_constraint = ctx$graph_constraint$trans_constraint,
        pvals = kmat,
        K = ctx$spec$K,
        alpha = ctx$alpha,
        num_threads = ctx$num_threads,
        stop_rule = ctx$stop_rule
    )
}

# Gain and decision times of the graph `(w, G)` on the kernel matrix `kmat`;
# an invalid graph has gain `-Inf`, so it is never chosen.
.hsd_eval_graph <- function(kmat, w, G, ctx) {
    res <- graph_shortcut_gsd(
        pvals = kmat,
        alpha = ctx$alpha,
        w = w,
        G = G,
        K = ctx$spec$K,
        stop_rule = ctx$stop_rule
    )
    gain <- if (is_graph_valid(w, G)) ctx$trial_success$func(res$time) else -Inf

    list(gain = gain, time = res$time)
}

#' COBYLA on the graph genes at one node
#'
#' The `nloptr::nloptr()` call of `.graph_optimise_local_gsd()` on the
#' closure `obj`, from `x0`, with the local options; the solution is decoded
#' with `param_to_solution(process = TRUE)` and `repair_graph()` and
#' evaluated by the kernel on `kmat` (the matrix `obj` scores). Evaluations
#' are counted by a wrapper around the objective.
#'
#' @param x0 Graph genes of the start.
#' @param obj The closure from `.hsd_objective(ctx, kmat)`.
#' @param kmat The kernel matrix of the node's array on the pinned rows.
#' @param ctx The Stage L context.
#'
#' @returns A list with `x` (the solution genes), `hyp_weight`,
#'   `trans_matrix`, `gain`, `time`, `blocks`, `n_eval` and
#'   `nloptr_output`.
#' @noRd
.hsd_refine_node <- function(x0, obj, kmat, ctx) {
    counter <- new.env(parent = emptyenv())
    counter$n <- 0L
    nlopt_obj_func <- function(x) {
        counter$n <- counter$n + 1L
        -obj(x)
    }

    nlopt_result <- nloptr::nloptr(
        x0 = x0,
        eval_f = nlopt_obj_func,
        lb = rep(0, length(x0)),
        ub = rep(1, length(x0)),
        opts = ctx$local_opts
    )

    constraint <- ctx$graph_constraint
    sol <- param_to_solution(nlopt_result$solution, constraint, process = TRUE)
    sol <- repair_graph(sol$hyp_weight, sol$trans_matrix, constraint)

    ev <- .hsd_eval_graph(kmat, sol$hyp_weight, sol$trans_matrix, ctx)

    list(
        x = nlopt_result$solution,
        hyp_weight = sol$hyp_weight,
        trans_matrix = sol$trans_matrix,
        gain = ev$gain,
        time = ev$time,
        blocks = .hsd_block_gains(ev$time, ctx),
        n_eval = counter$n,
        nloptr_output = nlopt_result
    )
}

# Graph genes of `(w, G)` for a warm start: the user-row branch of
# `.build_start_matrix()`, projection onto the constraint and then
# `create_start_params(sum_to_one_constraint = FALSE)`. Never
# `.build_start_matrix(...)[1, ]`, whose first row is the Holm seed.
.hsd_encode_graph <- function(gc, w, G) {
    projected <- closest_graph_to_constraints(gc, w, G)

    as.numeric(create_start_params(
        gc,
        w0 = projected$hyp_weight,
        G0 = projected$trans_matrix,
        sum_to_one_constraint = FALSE
    ))
}

# Gains on `min(B, nrow(time))` contiguous blocks of rows of a decision-time
# matrix (batch means, design section 1.4).
# nolint start: object_name_linter. B is the design's symbol
.hsd_block_gains <- function(
    time,
    ctx,
    func = ctx$trial_success$func,
    B = ctx$n_block
) {
    n <- nrow(time)
    n_block <- min(B, n)
    cuts <- floor(seq(0, n, length.out = n_block + 1L))

    vapply(
        seq_len(n_block),
        function(b) {
            func(time[seq.int(cuts[[b]] + 1, cuts[[b + 1L]]), , drop = FALSE])
        },
        double(1L)
    )
}
# nolint end

# Paired difference `a - b` and its batch-means standard error from the block
# gains of two evaluations on the same rows.
.hsd_block_diff <- function(blocks_a, blocks_b) {
    diffs <- blocks_a - blocks_b

    list(diff = mean(diffs), se = stats::sd(diffs) / sqrt(length(diffs)))
}

# The paired difference and standard error of two decision-time matrices on
# the same rows (design section 1.4): `diff` is the mean of the block
# differences of `func`, `se` their standard deviation over `sqrt(B)`.
# nolint start: object_name_linter. B is the design's symbol
.hsd_block_se <- function(time_a, time_b, func, B = 100) {
    .hsd_block_diff(
        .hsd_block_gains(time_a, func = func, B = B),
        .hsd_block_gains(time_b, func = func, B = B)
    )
}
# nolint end

# The connected run of nodes around the value `at` whose `diff` is at most
# `tol`, as `c(lower, upper)` values (the plateau convention of design
# section 1.5).
# the argument order is the brief's contract
# nolint next: function_argument_linter.
.hsd_plateau <- function(value, diff, tol = 5e-4, at) {
    ord <- order(value)
    value <- value[ord]
    diff <- diff[ord]

    lo <- match(at, value)
    hi <- lo
    while (lo > 1L && diff[[lo - 1L]] <= tol) {
        lo <- lo - 1L
    }
    while (hi < length(value) && diff[[hi + 1L]] <= tol) {
        hi <- hi + 1L
    }

    value[c(lo, hi)]
}


# Print and summary -------------------------------------------------------

# The spending block of `print.multigrain_graph_optimal()` for a result of
# the optimised-spending path (design section 1.5; build plan,
# reconciliation 8): one line per hypothesis, then one gain-profile block per
# parameter with the re-optimised and the fixed-graph plateaus.
.hsd_cat_print <- function(x) {
    hyp_labels <- x$spending
    m <- length(hyp_labels)
    hyp_names <- names(x$hyp_weight)
    named <- length(hyp_names) == m && all(nzchar(hyp_names))
    if (!named) {
        hyp_names <- .gsd_hyp_labels(m)
    }

    notes <- .hsd_spending_notes(x)
    hyp_lines <- paste0("  ", hyp_names, ": ", hyp_labels)
    has_note <- nzchar(notes$hyp)
    hyp_lines[has_note] <- paste0(
        hyp_lines[has_note],
        "   ",
        notes$hyp[has_note]
    )

    cli::cat_line("\nSpending functions:")
    cli::cat_line(c(hyp_lines, notes$extra))

    prof <- x$spending_profile
    param_names <- names(x$spending_param)
    for (name in param_names) {
        rows <- prof[prof$param == name, , drop = FALSE]
        if (nrow(rows) == 0L) {
            next
        }
        at <- x$spending_param[[name]]
        reopt_at <- if (at %in% rows$value) {
            at
        } else {
            rows$value[[which.min(rows$diff)]]
        }
        reopt <- .hsd_plateau(rows$value, rows$diff, at = reopt_at)

        # the fixed-graph plateau is measured from the largest fixed-graph
        # gain on the axis, max(gain_fixed_graph) - gain_fixed_graph, not
        # from the optimum's gain as the `diff` column is
        fixed_diff <- max(rows$gain_fixed_graph) - rows$gain_fixed_graph
        fixed_at <- if (at %in% rows$value) {
            at
        } else {
            rows$value[[which.max(rows$gain_fixed_graph)]]
        }
        fixed <- .hsd_plateau(rows$value, fixed_diff, at = fixed_at)

        cli::cat_line(
            "Gain profile in ",
            name,
            " (graph re-optimised at each value):"
        )
        cli::cat_line(
            "  within 5e-4 of the maximum on ",
            .hsd_fmt_interval(reopt),
            "; fixed-graph plateau ",
            .hsd_fmt_interval(fixed)
        )
    }

    if (isFALSE(x$profile_converged)) {
        rounds <- if (length(param_names) > 1L) " in 3 rounds" else ""
        cli::cat_line(
            "  profile did not converge",
            rounds,
            "; the plateau is indicative"
        )
    }

    invisible(x)
}

# The spending part of `summary.multigrain_graph_optimal()`: per parameter,
# its profile rows and the GA's coverage of its lattice.
.hsd_cat_summary <- function(object) {
    prof <- object$spending_profile
    coverage <- object$spending_search

    for (name in names(object$spending_param)) {
        cli::cat_line(
            cli::style_underline(
                paste0("\nSpending-parameter profile (", name, ")")
            ),
            ":"
        )
        rows <- prof[prof$param == name, , drop = FALSE]
        prof_table <- data.frame(
            value = rows$value,
            gain = round(rows$gain, 4),
            gain_fixed_graph = round(rows$gain_fixed_graph, 4),
            diff = round(rows$diff, 4),
            se = signif(rows$se, 2),
            n_eval = rows$n_eval
        )
        print(prof_table, row.names = FALSE)

        visits <- coverage$visits[coverage$param == name]
        if (length(visits) == 0L || all(is.na(visits))) {
            cli::cat_line("no global search")
        } else {
            cli::cat_line(sprintf(
                "GA visited %d of %d nodes",
                sum(visits > 0L, na.rm = TRUE),
                length(visits)
            ))
        }
    }

    invisible(object)
}

# The bracketed note after each hypothesis's spending label, and a line per
# parameter whose bounds could not go on a first member. The result does not
# store which parameter each hypothesis uses, only its value
# (`hyp_spending_param`), so the map is inferred by `.hsd_infer_groups()`
# unless a `spending_group` field is present. A short-circuited tagged
# hypothesis is one with no parameter whose label is `sfHSD(param = v)`.
.hsd_spending_notes <- function(x) {
    hyp_labels <- x$spending
    m <- length(hyp_labels)
    theta <- x$spending_param
    param_names <- names(theta)
    hyp_theta <- unname(as.numeric(x$hyp_spending_param))

    groups <- .hsd_infer_groups(hyp_theta, unname(theta), x$spending_group)

    notes <- rep("", m)
    described <- rep(FALSE, length(theta))
    for (i in seq_len(m)) {
        candidates <- groups[[i]]
        if (length(candidates) == 1L) {
            r <- candidates[[1L]]
            if (described[[r]]) {
                notes[[i]] <- sprintf("[optimised: %s]", param_names[[r]])
            } else {
                notes[[i]] <- sprintf(
                    "[optimised: %s, %s]",
                    param_names[[r]],
                    .hsd_fmt_lattice(x, r)
                )
                described[[r]] <- TRUE
            }
        } else if (length(candidates) > 1L) {
            notes[[i]] <- sprintf(
                "[optimised: %s]",
                paste(param_names[candidates], collapse = " or ")
            )
        } else if (is.na(hyp_theta[[i]])) {
            # nolint next: nonportable_path_linter. a regex, not a path
            pattern <- "^sfHSD\\(param = ([^ ()]+)\\)$"
            if (grepl(pattern, hyp_labels[[i]])) {
                notes[[i]] <- sprintf(
                    "[fixed at %s: one look at full information]",
                    sub(pattern, "\\1", hyp_labels[[i]])
                )
            }
        }
    }

    extra <- character()
    for (r in which(!described)) {
        extra <- c(
            extra,
            sprintf("  %s: %s", param_names[[r]], .hsd_fmt_lattice(x, r))
        )
    }

    list(hyp = notes, extra = extra)
}

# Candidate parameters of each hypothesis: a list of integer vectors, empty
# for a hypothesis without a parameter. With `given` (an integer group map)
# that map is used. Otherwise every map consistent with the values is
# enumerated, under the spec's numbering (parameters are numbered by their
# first member, every parameter has a member); a hypothesis gets every
# parameter some consistent map gives it, so it has one candidate whenever
# the values settle its parameter.
.hsd_infer_groups <- function(hyp_theta, theta, given = NULL) {
    m <- length(hyp_theta)
    out <- rep(list(integer()), m)
    if (!is.null(given)) {
        for (i in which(!is.na(given))) {
            out[[i]] <- as.integer(given[[i]])
        }
        return(out)
    }

    active <- which(!is.na(hyp_theta))
    d <- length(theta)
    state <- new.env(parent = emptyenv())
    state$maps <- list()
    walk <- function(pos, opened, acc) {
        if (pos > length(active)) {
            if (opened == d) {
                state$maps[[length(state$maps) + 1L]] <- acc
            }
            return(invisible())
        }
        value <- hyp_theta[[active[[pos]]]]
        for (r in seq_len(min(opened + 1L, d))) {
            if (isTRUE(abs(theta[[r]] - value) <= 1e-12)) {
                walk(pos + 1L, max(opened, r), c(acc, r))
            }
        }
    }
    walk(1L, 0L, integer())
    maps <- state$maps

    for (j in seq_along(active)) {
        i <- active[[j]]
        out[[i]] <- if (length(maps)) {
            sort(unique(vapply(maps, `[[`, integer(1L), j)))
        } else {
            which(abs(theta - hyp_theta[[i]]) <= 1e-12)
        }
    }

    out
}

# "bounds a to b, step s" for parameter `r`; the step is the spacing of its
# lattice in `spending_search` (the full lattice), else in the profile.
.hsd_fmt_lattice <- function(x, r) {
    name <- names(x$spending_param)[[r]]
    bounds <- x$spending_bounds
    bounds_r <- if (name %in% rownames(bounds)) {
        bounds[name, ]
    } else {
        bounds[r, ]
    }
    out <- sprintf(
        "bounds %s to %s",
        format(bounds_r[[1L]], digits = 15),
        format(bounds_r[[2L]], digits = 15)
    )

    nodes <- x$spending_search$value[x$spending_search$param == name]
    if (length(nodes) < 2L) {
        nodes <- x$spending_profile$value[x$spending_profile$param == name]
    }
    nodes <- sort(unique(nodes))
    if (length(nodes) >= 2L) {
        out <- paste0(out, ", step ", format(min(diff(nodes)), digits = 15))
    }

    out
}

# "[lower, upper]" with each end formatted on its own.
.hsd_fmt_interval <- function(range) {
    ends <- vapply(range, format, character(1L), digits = 15)
    paste0("[", ends[[1L]], ", ", ends[[2L]], "]")
}
