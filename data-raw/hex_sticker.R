## code to prepare `hexsticker` dataset goes here

library(hexSticker)
library(multigrain)
library(ggplot2)
library(ggraph)
library(dplyr)
library(tidygraph)
library(magick)


tm <- matrix(
    c(
        0, 0, 0, 1, 0,
        0.4, 0, 0.4, 0.1, 0.1,
        0, 0, 0, 0, 1,
        0, 0, 0, 0, 1,
        0, 0, 0, 1, 0
    ),
    nrow = 5,
    byrow = TRUE
)

dimnames(tm) <- list(
    c("H1", "H2", "H3", "H4", "H5"),
    c("H1", "H2", "H3", "H4", "H5")
)

hw <- c(H1 = 0.3, H2 = 0.4, H3 = 0.3, H4 = 0, H5 = 0)

edges <- multigrain:::create_edges(tm)
nodes <- multigrain:::create_nodes(hw, edges)


graph_for_plot <- igraph::graph_from_data_frame(
    edges,
    nodes,
    directed = TRUE
)

layout_df <- multigrain:::create_layout(
    graph_for_plot,
    nodes = nodes,
    edges = edges #,
    # root = root
)

layout_df$x[1] <- -0.5
layout_df$x[3] <- 0.5
layout_df$x[4] <- -0.25
layout_df$x[5] <- 0.25
layout_df$y[2] <- 1.6

hex_graph <- ggraph::ggraph(
    graph_for_plot,
    # custom layout data.frame
    layout = layout_df
) +
    # layer for the edges
    ggraph::geom_edge_fan(
        ggplot2::aes(
            # end and start points for the arrows
            end_cap = ggraph::circle(9, "mm"),
            start_cap = ggraph::circle(9, "mm")
        ),
        # arrow shape
        arrow = ggplot2::arrow(
            type = "closed",
            length = ggplot2::unit(4, "mm")
        ),
        # controls the width of the edge fan
        strength = 1.5
    ) +
    # layer for the nodes themselves
    ggraph::geom_node_point(
        ggplot2::aes(
            colour = as.character(.data$optimised)
        ),
        size = 20,
        alpha = 0.4,
        show.legend = FALSE
    ) +
    ggraph::geom_node_text(
        ggplot2::aes(
            #`name` is the name of the graph attribute
            label = name
        ),
        size = 18
    ) +
    ggplot2::coord_cartesian(
        clip = "off"
    ) +
    theme_transparent()


sticker <- sticker(
    hex_graph,
    package = "multigrain",
    # filename = "man/figures/logo.png",
    s_x = 1,
    s_y = 1.15,
    s_width = 1.4,
    s_height = 1.1,
    h_size = 5,
    h_color = "#F36633",
    h_fill = "lightyellow",
    # p_family = "FiraCode-Retina",
    p_color = "purple",
    # p_color = "#F36633",
    p_size = 40,
    p_y = 0.45,
    dpi = 300,
    url = "https://insert-the-new-URL-here/multigrain/",
    u_size = 9,
    u_color = "#F36633",
    u_y = 0.06
)

sticker

ggplot2::ggsave(
    filename = "man/figures/logo.png",
    plot = sticker,
    units = "mm"
)
