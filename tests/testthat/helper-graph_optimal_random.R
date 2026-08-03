#' Generate a random graph optimal
#'
#' A helper for testing. It is not exported.
#'
#' @param ... Arguments passed down to [graph_random()].
#'
#' @returns A minimal `multigrain_graph_optimal`. To be used only for tests.
#'
#' @noRd
graph_optimal_random <- function(...) {
    random_graph <- graph_random(...)

    graph_optimal(
        hyp_weight = random_graph$hyp_weight,
        trans_matrix = random_graph$trans_matrix
    )
}
