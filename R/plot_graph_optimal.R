#' Autoplot method for `multigrain_graph_optimal` objects
#'
#' @param object A `multigrain_graph_optimal` object.
#' @param root (integer-like) An optional argument containing a vector of nodes
#'   to plot as root.
#' @param digits Number of decimal places to round to (between 0 and 3).
#'   Defaults to `NULL` which will choose the number of digits based on how
#'   crowded the graph plot will be:
#'     * a single digit if more than 8 hypotheses
#'     * 2 digits if 5 to 8 hypotheses
#'     * 3 digits if 4 or fewer hypotheses
#' @param title (string) Optional plot title.
#' @param ... Additional arguments.
#'
#' @returns a `ggraph` / `ggplot` object.
#'
#' @export
#' @examples
#' library(ggplot2)
#'
#' random_graph <- graph_optimal_random(5)
#' autoplot(random_graph)
#'
#' # If you want to control which nodes are treated as root, you can pass them
#' # as a numeric vector via the `root` argument. For example, we want nodes 1
#' # and 3 to be plotted as root:
#' autoplot(random_graph, root = c(1, 3))
#'
#' # control the rounding with the `digits` argument.
#' autoplot(random_graph, digits = 2)
#'
autoplot.multigrain_graph_optimal <- function(
    object,
    ...,
    root = NULL,
    digits = NULL,
    title = NULL
) {
    check_integerish(root, allow_null = TRUE)
    rlang::check_number_whole(digits, min = 0, max = 3, allow_null = TRUE)
    rlang::check_string(title, allow_null = TRUE)

    hyp_weight <- object$hyp_weight
    trans_matrix <- object$trans_matrix

    if (rlang::is_null(digits)) {
        digits <- dplyr::case_when(
            length(hyp_weight) > 8 ~ 1L,
            length(hyp_weight) > 4 ~ 2L,
            .default = 3L
        )
    }

    fan_strength <- dplyr::case_when(
        length(hyp_weight) > 8 ~ 1.2,
        length(hyp_weight) > 4 ~ 1.3,
        .default = 1.4
    )

    # create the nodes and edges data frames
    edges <- create_edges(trans_matrix)
    nodes <- create_nodes(hyp_weight, edges = edges)

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

#' Plot graph
#'
#' A wrapper around [ggraph::ggraph()].
#'
#' @param graph (igraph) an `igraph` object to plot.
#' @param layout (data.frame) a data.frame containing the nodes layout. The
#'   output of `create_layout()`.
#' @param digits (integerish) number of digits to use for rounding the values
#'   before plotting. Defaults to 2.
#' @param fan_strength (numeric) fan strength. Passed to
#'   [ggraph::geom_edge_fan()].
#'
#' @returns a `ggraph`, `ggplot` object.
#' @noRd
plot_graph <- function(graph, layout, digits = 2, fan_strength = 1.4) {
    output <- ggraph::ggraph(
        graph,
        # custom layout data.frame
        layout = layout
    ) +
        # layer for the edges
        ggraph::geom_edge_fan(
            ggplot2::aes(
                # end and start points for the arrows
                end_cap = ggraph::circle(12, "mm"),
                start_cap = ggraph::circle(12, "mm"),
                label = round(.data$value, digits)
            ),
            # arrow shape
            arrow = ggplot2::arrow(
                type = "closed",
                length = ggplot2::unit(4, "mm")
            ),
            # controls the width of the edge fan
            strength = fan_strength,
            # vertical shift
            label_dodge = ggplot2::unit(3, "mm"),
            angle_calc = "along",
            label_size = 4,
            # controls where along the edge the label is drawn (between 0 and 1)
            # if 0.5 the diagonal labels will overlap
            label_pos = 0.45
        ) +
        # layer for the nodes themselves
        ggraph::geom_node_point(
            ggplot2::aes(
                colour = as.character(.data$optimised)
            ),
            size = 30,
            alpha = 0.4,
            show.legend = FALSE
        ) +
        ggraph::geom_node_text(
            ggplot2::aes(
                #`name` is the name of the graph attribute
                label = glue::glue(
                    "{name}
                    w: {round(weight, digits)}"
                )
            ),
            size = 5
        ) +
        ggplot2::coord_cartesian(
            clip = "off"
        ) +
        ggraph::theme_graph(
            base_family = "sans"
        )

    output
}

# nolint start: line_length_linter

#' Create layout
#'
#' Create the layout for the graph nodes. In the future this could be a
#' `create_layout()` method for `multigrain_graph_optimal` or
#' `multigrain_graph_constraint`:
#' https://www.data-imaginist.com/posts/2017-02-06-ggraph-introduction-layouts/index.html#adding-support-for-new-data-sources
#'
#' @param graph An `igraph` object.
#' @param nodes A `data.frame` of nodes, the output of `create_nodes()`.
#' @param edges A `data.frame` of edges, the output of `create_edges()`.
#' @param root A numeric vector indicating which nodes to be regarded as root in
#'   the tree layout. Passed down to [igraph::layout_as_tree()].
#'
#' @returns a `data.frame` with 2 columns (`x` and `y`) representing the
#'   coordinates for the nodes.
#' @noRd
create_layout <- function(graph, nodes, edges, root = NULL) {
    if (is.null(root)) {
        root <- estimate_root(nodes, edges)
    }

    tree_layout <- igraph::layout_as_tree(graph, root = root)
    colnames(tree_layout) <- c("x", "y")

    output <- tibble::as_tibble(tree_layout) |>
        dplyr::mutate(
            node = nodes$hypothesis
        ) |>
        tibble::column_to_rownames(
            var = "node"
        )

    output
}
# nolint end

#' Create graph edges
#'
#' Prepares a `data.frame` of edges for use with
#' [igraph::graph_from_data_frame()].
#'
#' @param trans_matrix transition matrix
#'
#' @returns A `data.frame` with 3 columns:
#'   * `from`: (character) the starting node of the edge
#'   * `to`: (character) the ending node of the edge
#'   * `value`: (numeric) the value of the edge
#' @noRd
#'
#' @examples
#' random_graph <- graph_random(3)
#' trans_matrix <- random_graph$trans_matrix
#'
#' edges_df <- multigrain:::create_edges(trans_matrix)
create_edges <- function(trans_matrix) {
    trans_matrix |>
        tibble::as_tibble(
            rownames = "from"
        ) |>
        tidyr::pivot_longer(
            -"from",
            names_to = "to"
        ) |>
        # we do not plot edges with value 0
        dplyr::filter(
            .data$value != 0 | is.na(.data$value)
        )
}

#' Create graph nodes
#'
#' Prepares a `data.frame` of nodes for use with
#' [igraph::graph_from_data_frame()].
#'
#' @param hyp_weight (numeric) A numeric vector of hypotheses weights.
#' @param edges (`data.frame`) A data.frame containing edges, the output of
#' `create_edges()`. It should have 3 columns: `from`, `to`, and `value`.
#'
#' @returns a `data.frame` with 4 columns:
#'   * `hypothesis`: (character) hypothesis / node name
#'   * `weight`: (numeric) hypothesis / node weight
#'   * `root`: (logical) indicates whether a node is root
#'   * `level`: (numeric) indicates the level of a node
#' @noRd
#'
#' @examples
#' random_graph <- graph_random(3)
#' trans_matrix <- random_graph$trans_matrix
#' hyp_weight <- random_graph$hyp_weight
#'
#' edges_df <- multigrain:::create_edges(trans_matrix)
#' nodes_df <- multigrain:::create_nodes(hyp_weight, edges_df)
create_nodes <- function(hyp_weight, edges) {
    optimised_nodes <- names(which(hyp_weight != 0 | is.na(hyp_weight)))

    output <- hyp_weight |>
        tibble::enframe(
            name = "hypothesis",
            value = "weight"
        ) |>
        dplyr::mutate(
            optimised = .data$hypothesis %in% optimised_nodes
        )

    output
}

#' Estimate root
#'
#' Root are generally considered the nodes with a weight different from 0. If
#' there are more than 2 nodes considered as root, `estimate_root()` tries to
#' figure out which ones are closer to the root. If we imagine a tree, we can
#' assume that nodes with fewer outgoing edges are closer to the root.
#'
#' @param nodes (data.frame) a data.frame of nodes. the output of
#'   `create_nodes()`.
#' @param edges (data.frame) a data.frame of edges. the output of
#'   `create_edges()`.
#'
#' @returns a numeric vector representing the indices of the root nodes
#' @noRd
estimate_root <- function(nodes, edges) {
    node_names <- nodes$hypothesis

    root <- nodes |>
        dplyr::filter(.data$optimised) |>
        dplyr::pull(.data$hypothesis)

    if (length(root) > 2) {
        # we assume the nodes with fewer outgoing edges are more important
        # (closer to the root)
        root <- edges |>
            dplyr::group_by(
                .data$from
            ) |>
            dplyr::summarise(
                edges = dplyr::n()
            ) |>
            dplyr::ungroup() |>
            dplyr::filter(.data$from %in% root) |>
            dplyr::filter(
                .data$edges == min(.data$edges)
            ) |>
            dplyr::pull(.data$from)
    }

    which(node_names %in% root)
}

#' @rdname autoplot.multigrain_graph_optimal
#' @param x A `multigrain_graph_optimal` object.
#' @export
plot.multigrain_graph_optimal <- function(
    x,
    ...,
    root = NULL,
    digits = NULL,
    title = NULL
) {
    autoplot(
        object = x,
        ...,
        root = root,
        digits = digits,
        title = title
    )
}
