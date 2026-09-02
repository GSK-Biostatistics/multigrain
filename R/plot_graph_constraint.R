#' Autoplot method for `multigrain_graph_constraint` objects
#'
#' @param object A `multigrain_graph_constraint` object.
#' @inheritParams autoplot.multigrain_graph_optimal root digits title
#' @param ... Additional arguments.
#'
#' @returns a `ggraph` / `ggplot` object.
#'
#' @export
#' @examples
#' library(ggplot2)
#' gc1 <- graph_constraint(
#'     c(NA, NA, 0, 0, 0),
#'     matrix(
#'        c(
#'            0, 0.8, 0.2, 0, 0,
#'            NA, 0, 0, NA, 0,
#'            NA, 0.8, 0, 0.1, NA,
#'            0.9, NA, NA, 0, NA,
#'            0, 0, 0, 1, 0
#'        ),
#'        nrow = 5,
#'        byrow = TRUE
#'    )
#' )
#'
#' autoplot(gc1)
#'
#' # If you want to control which nodes are treated as root, you can pass them
#' # as a numeric vector via the `root` argument. For example, we want nodes 1
#' # and 3 to be plotted as root:
#' ggplot2::autoplot(gc1, root = c(1, 3))
autoplot.multigrain_graph_constraint <- function(
    object,
    ...,
    root = NULL,
    digits = NULL,
    title = NULL
) {
    check_integerish(root, allow_null = TRUE)
    rlang::check_number_whole(digits, min = 0, max = 3, allow_null = TRUE)
    rlang::check_string(title, allow_null = TRUE)

    hyp_constraint <- object$hyp_constraint
    trans_constraint <- object$trans_constraint

    if (rlang::is_null(digits)) {
        digits <- dplyr::case_when(
            length(hyp_constraint) > 8 ~ 1L,
            length(hyp_constraint) > 4 ~ 2L,
            .default = 3L
        )
    }

    fan_strength <- dplyr::case_when(
        length(hyp_constraint) > 8 ~ 1.2,
        length(hyp_constraint) > 4 ~ 1.3,
        .default = 1.4
    )

    edges <- create_edges(trans_constraint)
    nodes <- create_nodes(hyp_constraint, edges = edges)

    graph_for_plot <- igraph::graph_from_data_frame(
        edges,
        nodes,
        directed = TRUE
    )

    layout_df <- create_layout(
        graph_for_plot,
        nodes = nodes,
        edges = edges,
        root = root
    )

    output <- plot_graph(
        graph = graph_for_plot,
        layout = layout_df,
        digits = digits,
        fan_strength = fan_strength
    )

    if (!rlang::is_null(title)) {
        output <- output +
            ggplot2::ggtitle(title) +
            # position the title centrally
            ggplot2::theme(
                plot.title = ggplot2::element_text(
                    hjust = 0.5
                )
            )
    }

    output
}

#' @rdname autoplot.multigrain_graph_constraint
#' @param x A `multigrain_graph_constraint` object.
#'
#' @export
plot.multigrain_graph_constraint <- function(
    x,
    ...,
    root = NULL,
    digits = NULL,
    title = NULL
) {
    autoplot(
        x,
        ...,
        root = root,
        digits = digits,
        title = title
    )
}
