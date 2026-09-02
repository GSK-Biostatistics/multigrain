# deprecate_warn should be replaced with deprecate_stop after v0.3.0

#' Optimise graph (deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' Deprecated to adopt a more consistent naming convention. Please use
#' [graph_optimise()] instead.
#'
#' @inheritParams graph_optimise
#'
#' @keywords internal
optimise_graph <- function(...) {
    lifecycle::deprecate_warn(
        when = "0.3.0",
        what = "optimise_graph()",
        with = "graph_optimise()"
    )

    graph_optimise(...)
}

#' @keywords internal
#' @rdname optimise_graph
#' @usage NULL
optimize_graph <- function(...) {
    lifecycle::deprecate_warn(
        when = "0.3.0",
        what = "optimize_graph()",
        with = "graph_optimise()"
    )

    graph_optimise(...)
}
