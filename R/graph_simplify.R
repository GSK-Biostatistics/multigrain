# Probability that a stage-2 mutation performs a zeroing move rather than a
# Cauchy perturbation, split equally between the two zeroing variants (see
# `.make_cauchy_mutation_multi()`).
#
# OPEN, design record item O1. The value 0.2 and the equal split are the
# record's proposal, not a measured optimum, and should be experimented on.
# The evidence needed to close it: on m = 4 examples, the fraction of stage-2
# generations in which the elite's edge count decreases, for p_zero in
# {0.1, 0.2, 0.4}. Nothing else depends on the value, so re-tuning it is a
# one-line change.
.simplify_p_zero <- 0.2


#' Simplify an optimised graph within a trial success budget
#'
#' Takes a graph from [graph_optimise()] and searches for one with fewer edges,
#' giving up at most a stated fraction of its trial success measure. Optimised
#' graphs often carry several small transition weights that contribute almost
#' nothing but that a reviewer still has to reason about; `graph_simplify()`
#' trades a little expected gain for a graph that is easier to present.
#'
#' @details
#' The graph you pass in is the *reference*. Its trial success measure is
#' recomputed on the `pvals` you supply, and the returned graph is guaranteed
#' to reach at least `(1 - gain_tolerance)` times that value on the same
#' simulated trials. Among the graphs that clear that floor, the search prefers
#' the one with fewest edges, and among graphs with the same number of edges,
#' the one with the highest trial success.
#'
#' Family-wise error control is unaffected. A Bonferroni-based graphical
#' procedure controls the family-wise error rate strongly for any weights
#' summing to at most one and any transition matrix with non-negative entries,
#' zero diagonal and row sums at most one. Removing an edge and renormalising
#' its row preserves all three, so only the objective changes, never the class
#' of procedures searched.
#'
#' ## Edge removal during the global search
#'
#' The mutation the genetic algorithm normally uses selects each parameter with
#' probability 0.1 and replaces it with a draw from a Cauchy distribution
#' centred on its current value. That moves weight around well, but it is a
#' poor way to *remove* an edge: "removed" means the decoded entry falls below
#' `1e-5`, and a Cauchy draw lands there about once in 100,000 tries. Crossover
#' does not help, because a convex combination of two parents is zero only
#' where both parents are already zero, and the local search is a continuous
#' method that never aims for the threshold. Left alone, the search would
#' *prefer* sparser graphs but almost never *produce* one to compare.
#'
#' So each mutation call now first tosses a coin. If it comes up, the call
#' performs a zeroing move instead, of one of two kinds. It either picks one
#' free transition parameter that is currently above the threshold and sets it
#' to exactly zero, removing that edge outright; or it picks a row and rescales
#' its parameters so that the row's *derived* entry -- the one computed as one
#' minus the others, which no parameter can reach -- falls just below the
#' threshold. Without the second kind, one edge per row would be unreachable by
#' the global search and removable only by the pruning step at the end. If the
#' coin does not come up, or there is nothing left to zero, the call falls
#' through to the ordinary Cauchy perturbation.
#'
#' With `global_search = FALSE` no genetic algorithm runs, the local search
#' starts from the reference, and edges are removed only by pruning. That is
#' cheaper, but it cannot find removals that only become affordable once the
#' remaining weights are re-tuned, nor move to a different topology.
#'
#' ## What to expect
#'
#' Edge costs are lumpy, so the fraction of the budget actually spent varies:
#' the cap is the guarantee, the spend is whatever the graph allows. On a
#' four-hypothesis example with seven edges, a budget of `1e-3` bought one edge
#' for 9% of the budget; at `1e-2` the same, because the next cheapest edge
#' cost more than the budget that remained; at `2e-2`, two edges for 61% of it.
#' `gain_tolerance = 0` is well defined and useful: it asks for the fewest
#' edges among graphs whose trial success does not fall at all.
#'
#' Called after [graph_optimise()], the whole two-stage run costs roughly 1.3
#' to 1.6 times a single optimisation. Trying several values of
#' `gain_tolerance` does not repeat the first optimisation.
#'
#' Simplifying an already-simplified graph compounds the cap: the second call
#' measures its budget against the first call's result, not against the
#' original graph.
#'
#' @param graph_optimal A `multigrain_graph_optimal` object, as returned by
#'   [graph_optimise()]. This is the reference graph.
#'
#'   The trial success measure is compiled, and a compiled function does not
#'   survive being saved and reloaded. `graph_simplify()` rebuilds a dead
#'   function automatically from the objective stored on the object. A live
#'   function is checked against that objective before it is used; inconsistent
#'   or invalid objects produce an error.
#' @param pvals A numeric matrix of p-values, one row per simulated trial and
#'   one column per hypothesis. It is supplied again because the optimised
#'   object does not store it. It need not be the matrix used in stage 1: the
#'   cap is always measured on the matrix given here, so `gain_reference` is
#'   recomputed on it and may differ from the object's stored
#'   `$power$trial_success`.
#' @inheritParams rlang::args_dots_empty
#' @param gain_tolerance A single number in `[0, 1]`. The fraction of the
#'   reference graph's trial success measure, on the supplied `pvals`, that may
#'   be given up in total. Defaults to `1e-3`.
#' @param alpha A single numeric value representing the overall one-sided
#'   significance level. Defaults to the value stored on `graph_optimal`; an
#'   explicit value overrides it. Objects saved by a version that did not store
#'   `alpha` must supply it here.
#' @param global_search A logical indicating whether to run the genetic
#'   algorithm before the local optimisation. Defaults to `TRUE`.
#' @inheritParams graph_optimise num_threads
#' @param control An optional `multigrain_control` object. Defaults to the
#'   settings stored on `graph_optimal` with the genetic algorithm's `run`
#'   halved, since the search starts in the right basin and does not need
#'   stage 1's floor on cost. When the global search can run, the `GA::ga()`
#'   options `mutation` and `suggestions` are reserved: simplification requires
#'   its support-changing mutation and reference-first seed population. Setting
#'   either option in an explicit `control` object is an error. If they are
#'   inherited from `graph_optimal`, they are ignored with a warning. They are
#'   irrelevant and left alone when no global search runs.
#' @inheritParams graph_optimise verbose
#'
#' @returns A `multigrain_graph_optimal` object, so that [print()], [summary()],
#'   [plot()] and [calc_power_pvals()] work on it unchanged. Its `sparsity`
#'   element is a list containing:
#'
#'   * `gain_tolerance`: the cap, as supplied.
#'   * `reference`: a list with the reference graph's `hyp_weight` and
#'     `trans_matrix`, so the two graphs can be plotted side by side.
#'   * `gain_reference`: the reference graph's trial success on `pvals`.
#'   * `gain`: the returned graph's trial success, equal to
#'     `$power$trial_success`.
#'   * `gain_loss`: `gain_reference - gain`.
#'   * `gain_loss_fraction`: `gain_loss / gain_reference`.
#'   * `budget`: `gain_tolerance * gain_reference`.
#'   * `n_edges_reference`, `n_edges`: total non-zero transition entries,
#'     including any pinned by the constraint, before and after.
#'   * `n_edges_free_reference`, `n_edges_free`: non-zero *free* entries,
#'     before and after.
#'   * `prune_loss`: the exact loss of trial success across the removals that
#'     pruning accepted. Negative if pruning raised the measure.
#'   * `source`: `"global"`, `"local"` or `"reference"`, saying which stage
#'     produced the returned graph.
#'
#' @seealso [graph_optimise()] for the first-stage optimisation.
#'
#' @export
#' @examples
#'
#' pvals <- simulate_pvalues(
#'   power_nominal = c(0.9, 0.85, 0.8, 0.75),
#'   corr_matrix = diag(4),
#'   nsim = 5000
#' )
#'
#' ts <- trial_success(r1 + r2 + r3 + r4)
#'
#' \donttest{
#' optimal <- graph_optimise(
#'   pvals = pvals,
#'   graph_constraint = graph_constraint_free(4),
#'   trial_success = ts
#' )
#'
#' # Give up at most 0.5% of the trial success measure for a smaller graph
#' simpler <- graph_simplify(optimal, pvals, gain_tolerance = 5e-3)
#'
#' simpler$sparsity$n_edges_reference
#' simpler$sparsity$n_edges
#' }
graph_simplify <- function(
    graph_optimal,
    pvals,
    ...,
    gain_tolerance = 1e-3,
    alpha = NULL,
    global_search = TRUE,
    num_threads = 1L,
    control = NULL,
    verbose = multigrain_verbosity()
) {
    check_graph_optimal(graph_optimal)
    check_double_matrix(pvals)
    rlang::check_dots_empty()
    rlang::check_number_decimal(gain_tolerance, min = 0, max = 1)
    rlang::check_number_decimal(alpha, min = 0, max = 1, allow_null = TRUE)
    check_logical(global_search, allow_na = FALSE)
    rlang::check_number_whole(num_threads, min = 1)
    check_control(control, allow_null = TRUE)
    verbose <- .resolve_verbose(verbose)

    constraints <- graph_optimal$constraints
    trial_success <- graph_optimal$trial_success

    if (is.null(constraints) || is.null(trial_success)) {
        cli::cli_abort(
            "{.arg graph_optimal} must carry its constraint and trial-success \\
            objects."
        )
    }
    trial_success <- .restore_trial_success(trial_success)
    graph_optimal$trial_success <- trial_success
    if (trial_success$m != ncol(pvals)) {
        cli::cli_abort(
            "{.arg pvals} must have one column per hypothesis of \\
            {.arg graph_optimal}."
        )
    }

    alpha <- alpha %||% graph_optimal$alpha
    if (is.null(alpha)) {
        cli::cli_abort(
            "{.arg graph_optimal} does not store {.arg alpha}; supply it \\
            explicitly."
        )
    }

    ref <- list(
        hyp_weight = unname(graph_optimal$hyp_weight),
        trans_matrix = unname(graph_optimal$trans_matrix)
    )

    ref_power <- calc_power_pvals(
        pvals,
        hyp_weight = ref$hyp_weight,
        trans_matrix = ref$trans_matrix,
        alpha = alpha,
        custom_power = list(trial_success = trial_success)
    )
    u_ref <- ref_power$trial_success

    if (u_ref <= 0) {
        cli::cli_abort(
            "The reference graph has non-positive trial success on \\
            {.arg pvals}; there is nothing to trade."
        )
    }

    u_range <- .trial_success_range(trial_success)
    free_mask <- is.na(constraints$trans_constraint)
    n_removable <- sum(pmax(rowSums(free_mask) - 1L, 0L))

    control_supplied <- !is.null(control)
    if (!control_supplied) {
        control <- graph_optimal_get_control(graph_optimal) %||%
            multigrain_control()
        if (!is.null(control$global_opt$run)) {
            control$global_opt$run <- max(1L, control$global_opt$run %/% 2L)
        }
    }

    if (global_search && n_removable > 0L) {
        reserved <- intersect(
            names(control$global_opt),
            c("mutation", "suggestions")
        )
        if (length(reserved) > 0L) {
            reserved_text <- paste0("`", reserved, "`", collapse = ", ")
            if (control_supplied) {
                cli::cli_abort(c(
                    "{.arg control} sets reserved global optimisation \\
                    options: {reserved_text}.",
                    i = "{.fn graph_simplify} supplies its own support-changing \\
                    mutation and reference-first suggestions."
                ))
            }
            cli::cli_warn(c(
                "The stored control contains global optimisation options that \\
                {.fn graph_simplify} cannot reuse: {reserved_text}.",
                i = "Using the simplification mutation and reference-first \\
                suggestions instead."
            ))
            control$global_opt[reserved] <- NULL
        }
    }
    control <- control_prepare(control, pvals = pvals, verbose = verbose)

    if (n_removable == 0L) {
        cli::cli_warn(
            "No removable edges under this constraint; returning the \\
            reference graph."
        )
        return(.simplify_result(
            graph_optimal = graph_optimal,
            ref = ref,
            final = ref,
            final_power = ref_power,
            u_ref = u_ref,
            gain_tolerance = gain_tolerance,
            prune_loss = 0,
            result_source = "reference",
            control = control,
            alpha = alpha,
            global_search = global_search,
            ga_result = NULL,
            local_result = NULL
        ))
    }

    objective_args <- list(
        gain_tolerance = gain_tolerance,
        ref_graph = ref,
        u_range = u_range
    )

    seeds <- .build_simplify_seeds(
        graph_constraint = constraints,
        ref_graph = ref,
        pop_size = control$global_opt$popSize,
        start_graph = graph_optimal$start_graph,
        population = if (is_ga(graph_optimal$global_output)) {
            graph_optimal$global_output@population
        }
    )

    ga_result <- NULL
    x0 <- .encode_graph(constraints, ref$hyp_weight, ref$trans_matrix)

    if (global_search) {
        ga_result <- .graph_optimise_ga(
            pvals = pvals,
            graph_constraint = constraints,
            trial_success = trial_success,
            nsim = control$nsim_global,
            global_opts = control$global_opt,
            alpha = alpha,
            num_threads = num_threads,
            verbose = verbose,
            suggestions = seeds,
            p_zero = .simplify_p_zero,
            objective_args = objective_args
        )
        x0 <- pmin(pmax(ga_result$ga_output@solution[1, ], 0), 1)
    }

    local_result <- .graph_optimise_local(
        pvals = pvals,
        graph_constraint = constraints,
        trial_success = trial_success,
        nsim = control$nsim_local,
        local_opts = control$local_opt,
        alpha = alpha,
        num_threads = num_threads,
        x0 = x0,
        verbose = verbose,
        objective_args = objective_args
    )

    best <- suppressWarnings(choose_graph(ga_result, local_result))

    threshold_full <- (1 - gain_tolerance) * u_ref
    edge_price <- (u_range[["max"]] - threshold_full) + 1

    score <- function(hyp_weight, trans_matrix) {
        u <- calc_power_pvals(
            pvals,
            hyp_weight = hyp_weight,
            trans_matrix = trans_matrix,
            alpha = alpha,
            custom_power = trial_success
        )$custom_power
        .lexico(
            u = u,
            n_edges = sum(trans_matrix[free_mask] != 0),
            threshold = threshold_full,
            edge_price = edge_price,
            n_free = sum(free_mask)
        )
    }

    # The searches can return a graph that is not valid at all (both stages
    # failing is warned about by choose_graph()). Pruning it would abort, so
    # skip straight to the reference, which is feasible on this sample by
    # construction.
    pruned <- NULL
    if (is_graph_valid(best$hyp_weight, best$trans_matrix)) {
        pruned <- prune_graph(
            pvals = pvals,
            hyp_weight = best$hyp_weight,
            trans_matrix = best$trans_matrix,
            trial_success = trial_success,
            graph_constraint = constraints,
            alpha = alpha,
            gamma = 1,
            threshold = threshold_full,
            verbose = verbose
        )
    }

    use_ref <- is.null(pruned) ||
        score(ref$hyp_weight, ref$trans_matrix) >=
            score(pruned$hyp_weight, pruned$trans_matrix)

    final <- if (use_ref) ref else pruned
    result_source <- if (use_ref) "reference" else best$source

    if (verbose != "silent") {
        cli::cli_progress_step("Evaluating trial success of simplified graph")
    }

    final_power <- calc_power_pvals(
        pvals,
        hyp_weight = final$hyp_weight,
        trans_matrix = final$trans_matrix,
        alpha = alpha,
        custom_power = list(trial_success = trial_success)
    )

    .simplify_result(
        graph_optimal = graph_optimal,
        ref = ref,
        final = final,
        final_power = final_power,
        u_ref = u_ref,
        gain_tolerance = gain_tolerance,
        prune_loss = if (use_ref) 0 else pruned$prune_loss,
        result_source = result_source,
        control = control,
        alpha = alpha,
        global_search = global_search,
        ga_result = ga_result,
        local_result = local_result
    )
}


#' Assemble the `multigrain_graph_optimal` returned by `graph_simplify()`
#'
#' Applies the hypothesis names, builds the `sparsity` element of the design
#' record's section 3.9, and calls the constructor. `power` is computed by the
#' caller from the graph actually returned and passed in here, so the reported
#' trial success always describes the graph the user sees.
#'
#' @param graph_optimal The reference object.
#' @param ref (list) The reference `hyp_weight` and `trans_matrix`, unnamed.
#' @param final (list) The graph being returned.
#' @param final_power (list) `calc_power_pvals()` output for `final`.
#' @param u_ref (numeric) Reference trial success on the supplied `pvals`.
#' @param gain_tolerance (numeric) The cap, as supplied.
#' @param prune_loss (numeric) Loss across the removals pruning accepted.
#' @param result_source (character) `"global"`, `"local"` or
#'   `"reference"`.
#' @param control The prepared control object.
#' @param alpha (numeric) The significance level used.
#' @param global_search (logical) Whether the GA was run.
#' @param ga_result,local_result The internal optimiser results, or `NULL`.
#'
#' @returns A `multigrain_graph_optimal`.
#' @noRd
.simplify_result <- function(
    graph_optimal,
    ref,
    final,
    final_power,
    u_ref,
    gain_tolerance,
    prune_loss,
    result_source,
    control,
    alpha,
    global_search,
    ga_result,
    local_result
) {
    constraints <- graph_optimal$constraints
    free_mask <- is.na(constraints$trans_constraint)
    gc_names <- graph_constraint_get_names(constraints)

    hyp_weight <- final$hyp_weight
    trans_matrix <- final$trans_matrix
    names(hyp_weight) <- gc_names
    dimnames(trans_matrix) <- list(gc_names, gc_names)

    ref_weight <- ref$hyp_weight
    ref_matrix <- ref$trans_matrix
    names(ref_weight) <- gc_names
    dimnames(ref_matrix) <- list(gc_names, gc_names)

    gain <- final_power$trial_success

    sparsity <- list(
        gain_tolerance = gain_tolerance,
        reference = list(
            hyp_weight = ref_weight,
            trans_matrix = ref_matrix
        ),
        gain_reference = u_ref,
        gain = gain,
        gain_loss = u_ref - gain,
        gain_loss_fraction = (u_ref - gain) / u_ref,
        budget = gain_tolerance * u_ref,
        n_edges_reference = sum(ref$trans_matrix != 0),
        n_edges = sum(final$trans_matrix != 0),
        n_edges_free_reference = sum(ref$trans_matrix[free_mask] != 0),
        n_edges_free = sum(final$trans_matrix[free_mask] != 0),
        prune_loss = prune_loss,
        source = result_source
    )

    opt_source <- if (identical(result_source, "reference")) {
        "reference"
    } else {
        paste0("simplify:", result_source)
    }

    graph_optimal(
        hyp_weight = hyp_weight,
        trans_matrix = trans_matrix,
        constraints = constraints,
        trial_success = graph_optimal$trial_success,
        power = final_power,
        solution = list(
            opt_source = opt_source,
            graph_valid = c(
                local = if (is.null(local_result)) {
                    NA
                } else {
                    isTRUE(local_result$is_graph_valid)
                },
                global = if (is.null(ga_result)) {
                    NA
                } else {
                    isTRUE(ga_result$is_graph_valid)
                }
            )
        ),
        global_search = global_search,
        control = control,
        global_output = ga_result$ga_output,
        local_output = local_result$local_output,
        start_graph = graph_optimal$start_graph,
        alpha = alpha,
        sparsity = sparsity
    )
}
