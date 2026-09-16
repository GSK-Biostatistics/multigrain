#' Choose the best graph from GA and local optimisation results
#'
#' Compares the global (GA) result and local (COBYLA) result based on graph
#' validity and the objective evaluated on the full p-value matrix. Returns
#' the chosen `hyp_weight`, `trans_matrix`, and a source label.
#'
#' The objective is `ga_objective` / `local_objective`: the lexicographic
#' score when `graph_simplify()` set a threshold, and the raw trial success
#' otherwise, so this comparison is unchanged for `graph_optimise()`. Results
#' that predate those elements fall back to the trial-success values.
#'
#' @param ga_result List returned by `.graph_optimise_ga()`, or `NULL` when
#'   `global_search == FALSE`.
#' @param local_result List returned by `.graph_optimise_local()`.
#'
#' @returns A list with:
#'   - `hyp_weight`: numeric vector from the chosen result
#'   - `trans_matrix`: numeric matrix from the chosen result
#'   - `source`: character scalar, one of `"local"`, `"global"`
#'
#' @noRd
choose_graph <- function(ga_result, local_result) {
    if (is.null(ga_result)) {
        return(list(
            hyp_weight = local_result$local_hyp_weight,
            trans_matrix = local_result$local_trans_matrix,
            source = "local"
        ))
    }

    ga_valid <- isTRUE(ga_result$is_graph_valid)
    local_valid <- isTRUE(local_result$is_graph_valid)

    ga_score <- ga_result$ga_objective %||% ga_result$ga_trial_success
    local_score <- local_result$local_objective %||%
        local_result$local_trial_success

    # GA only overrides local if
    #  - GA is valid and beats local on the objective
    #  - Local is invalid but GA produced a valid graph
    use_ga <- ga_valid && (!local_valid || ga_score > local_score)

    if (use_ga) {
        return(list(
            hyp_weight = ga_result$ga_hyp_weight,
            trans_matrix = ga_result$ga_trans_matrix,
            source = "global"
        ))
    }

    if (!local_valid && !ga_valid) {
        cli::cli_warn("Neither GA nor local produced a valid graph.")
    }

    # Default: return local
    list(
        hyp_weight = local_result$local_hyp_weight,
        trans_matrix = local_result$local_trans_matrix,
        source = "local"
    )
}
