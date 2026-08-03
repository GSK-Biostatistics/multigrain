# vdiffr::expect_doppelganger() tests are a bit flaky and always fail on CI.
test_that("graph plotting with graph_custom_power", {
    graph_custom_power <- readRDS(test_path("data", "graph_custom_power.rds"))

    graph_custom_power_plot <- autoplot(graph_custom_power)

    expect_s3_class(
        graph_custom_power_plot,
        c(
            "ggraph",
            "ggplot2::ggplot",
            "ggplot",
            "ggplot2::gg",
            "S7_object",
            "gg"
        )
    )

    skip_on_ci()
    vdiffr::expect_doppelganger(
        title = "Optimal Graph for Custom Power Metric",
        autoplot(graph_custom_power)
    )
})

test_that("graph plotting with graph_average_power", {
    graph_average_power <- readRDS(test_path("data", "graph_average_power.rds"))

    graph_average_power_plot <- autoplot(graph_average_power)
    expect_s3_class(
        graph_average_power_plot,
        c(
            "ggraph",
            "ggplot2::ggplot",
            "ggplot",
            "ggplot2::gg",
            "S7_object",
            "gg"
        )
    )

    skip_on_ci()
    vdiffr::expect_doppelganger(
        title = "Optimal Graph for Average Power Metric",
        autoplot(graph_average_power)
    )
})

test_that("plot ferber_et_al_2011 graph", {
    names <- glue::glue("H{1:9}")

    hyp_weight <- c(0.2, 0.2, 0.3, 0.3, 0, 0, 0, 0, 0)
    names(hyp_weight) <- names

    trans_matrix <- matrix(
        c(
            0, 1 / 3, 1 / 3, 1 / 3,     0,     0,     0,     0,     0,
            1 / 3,     0, 1 / 3, 1 / 3,     0,     0,     0,     0,     0,
            1 / 6, 1 / 6,     0, 1 / 6, 1 / 4, 1 / 4,     0,     0,     0,
            1 / 6, 1 / 6, 1 / 6,     0,     0,     0, 1 / 6, 1 / 6, 1 / 6,
            0,     0,   0.2,     0,     0,   0.8,     0,     0,     0,
            0,     0,   0.2,     0,     0.8,   0,     0,     0,     0,
            0,     0,     0,   0.2,     0,     0,     0,   0.4,   0.4,
            0,     0,    0,    0.2,     0,     0,   0.4,     0,   0.4,
            0,     0,    0,    0.2,     0,     0,   0.4,   0.4,     0
        ),
        nrow = 9,
        byrow = TRUE
    )

    dimnames(trans_matrix) <- list(names, names)

    ferber_et_al_graph <- graph_optimal(
        hyp_weight = hyp_weight,
        trans_matrix = trans_matrix
    )

    expect_no_error(
        autoplot(ferber_et_al_graph)
    )

    skip_on_ci()
    vdiffr::expect_doppelganger(
        title = "Graph from Ferber et al. 2011",
        autoplot(ferber_et_al_graph)
    )
})

test_that("plot bretz_et_al_2009 graph", {
    names <- c("H11", "H12", "H13", "H21", "H22", "H23")

    hyp_weight <- c(1, 0, 0, 0, 0, 0)
    names(hyp_weight) <- names

    trans_matrix <- matrix(
        c(
            0, 0.5, 0,   0.5,   0,   0,
            0,   0, 0.2, 0.8,   0,   0,
            0,   0,   0,   0, 0.8, 0.2,
            0, 0.8,   0,   0, 0.2,   0,
            0,   0, 0.8,   0,   0, 0.2,
            0,   0,   0,   0,   0,   0
        ),
        nrow = 6,
        byrow = TRUE
    )

    dimnames(trans_matrix) <- list(names, names)

    bretz_et_al_2009_graph <- graph_optimal(
        hyp_weight = hyp_weight,
        trans_matrix = trans_matrix
    )

    expect_no_error(
        autoplot(bretz_et_al_2009_graph)
    )

    skip_on_ci()
    vdiffr::expect_doppelganger(
        title = "Graph from Bretz et al. 2009",
        autoplot(bretz_et_al_2009_graph)
    )
})

test_that("plot random graph", {
    set.seed(1)
    random_graph <- graph_optimal_random(5)

    expect_no_error(
        autoplot(random_graph)
    )

    skip_on_ci()
    vdiffr::expect_doppelganger(
        title = "Random Graph",
        autoplot(random_graph)
    )
})

test_that("users can control the number of digits and the root", {
    set.seed(1)
    random_graph <- graph_optimal_random(9)

    expect_no_error(
        autoplot(random_graph, digits = 2)
    )

    skip_on_ci()
    vdiffr::expect_doppelganger(
        title = "Random Graph with digits",
        autoplot(random_graph, digits = 2)
    )

    vdiffr::expect_doppelganger(
        title = "Random Graph with root",
        autoplot(random_graph, root = c(1, 2))
    )
})

test_that("create_edges() works", {
    set.seed(1L)
    trans_matrix <- graph_random(4)$trans_matrix

    edges_df <- create_edges(trans_matrix)

    expect_snapshot(
        create_edges(trans_matrix)
    )

    expect_s3_class(
        edges_df,
        "tbl_df"
    )

    expect_named(
        edges_df,
        c("from", "to", "value")
    )

    expect_identical(
        unique(edges_df$from),
        c("H1", "H2", "H3", "H4")
    )

    h1_h3_edge_value <- edges_df |>
        dplyr::filter(
            from == "H1",
            to == "H3"
        ) |>
        dplyr::pull(value)

    expect_identical(
        h1_h3_edge_value,
        trans_matrix["H1", "H3"]
    )
})

test_that("create_nodes() works", {
    set.seed(1L)
    random_graph <- graph_random(4)
    trans_matrix <- random_graph$trans_matrix
    hyp_weight <- random_graph$hyp_weight

    edges_df <- create_edges(trans_matrix)

    nodes <- create_nodes(hyp_weight, edges = edges_df)

    expect_snapshot(
        create_nodes(
            hyp_weight,
            edges = edges_df
        )
    )

    expect_s3_class(
        nodes,
        "tbl_df"
    )

    expect_named(
        nodes,
        c("hypothesis", "weight", "optimised")
    )

    expect_identical(
        unique(nodes$hypothesis),
        c("H1", "H2", "H3", "H4")
    )

    expect_identical(
        nodes$weight,
        hyp_weight,
        ignore_attr = "names"
    )
})

test_that("modify node position manually", {
    names <- glue::glue("H{1:9}")

    hyp_weight <- c(0.2, 0.2, 0.3, 0.3, 0, 0, 0, 0, 0)
    names(hyp_weight) <- names

    trans_matrix <- matrix(
        c(
            0, 1 / 3, 1 / 3, 1 / 3,     0,     0,     0,     0,     0,
            1 / 3,     0, 1 / 3, 1 / 3,     0,     0,     0,     0,     0,
            1 / 6, 1 / 6,     0, 1 / 6, 1 / 4, 1 / 4,     0,     0,     0,
            1 / 6, 1 / 6, 1 / 6,     0,     0,     0, 1 / 6, 1 / 6, 1 / 6,
            0,     0,   0.2,     0,     0,   0.8,     0,     0,     0,
            0,     0,   0.2,     0,     0.8,   0,     0,     0,     0,
            0,     0,     0,   0.2,     0,     0,     0,   0.4,   0.4,
            0,     0,    0,    0.2,     0,     0,   0.4,     0,   0.4,
            0,     0,    0,    0.2,     0,     0,   0.4,   0.4,     0
        ),
        nrow = 9,
        byrow = TRUE
    )

    dimnames(trans_matrix) <- list(names, names)

    ferber_et_al_graph <- graph_optimal(
        hyp_weight = hyp_weight,
        trans_matrix = trans_matrix
    )

    plot <- autoplot(ferber_et_al_graph)

    plot$data$x[1] <- -1.75
    plot$data$x[2] <- 0.75
    plot$data$y[8] <- -1

    skip_on_ci()
    vdiffr::expect_doppelganger(
        title = "Ferber et al. manual adjustment",
        plot
    )
})

test_that("another custom graph", {
    names <- glue::glue("H{1:8}")

    hyp_weight <- c(0.883, 0.117, 0, 0, 0, 0, 0, 0)
    names(hyp_weight) <- names

    trans_matrix <- matrix(
        c(
            0,     0.751, 0.093, 0,     0.156, 0,     0,     0,
            0.555, 0,     0,     0.193, 0,     0.252, 0,     0,
            0,     0.746, 0,     0,     0.252, 0,     0.002, 0,
            0.795, 0,     0,     0,     0,     0.204, 0,     0.001,
            0,     0.702, 0.298, 0,     0,     0,     0,     0,
            0.966, 0,     0,     0.033, 0,     0,     0,     0.001,
            0,     0.061, 0.919, 0,     0.021, 0,     0,     0,
            0.578, 0,     0,     0.029, 0,     0.393, 0,     0
        ),
        nrow = 8,
        byrow = TRUE
    )

    dimnames(trans_matrix) <- list(names, names)

    custom_graph <- graph_optimal(
        hyp_weight = hyp_weight,
        trans_matrix = trans_matrix
    )

    expect_no_error(
        autoplot(custom_graph)
    )

    skip_on_ci()
    vdiffr::expect_doppelganger(
        title = "Custom graph with 8 hypotheses",
        autoplot(custom_graph)
    )
})

test_that("autoplot and plot produce the same output: graph_average_power", {
    graph_average_power <- readRDS(test_path("data", "graph_average_power.rds"))

    expect_no_error(
        plot(graph_average_power)
    )

    # there are some non-deterministic elements in the plot objects themselves,
    # so we cannot compare them (mainly due to {igraph}). `equivalent_ggplot2()`
    # writes them to an .svg file and then compares the md5sums
    expect_true(
        equivalent_ggplot2(
            plot(graph_average_power),
            autoplot(graph_average_power)
        )
    )

    skip_on_ci()
    vdiffr::expect_doppelganger(
        title = "Graph average power with plot()",
        plot(graph_average_power)
    )
})

test_that("autoplot and plot produce the same output: graph_custom_power", {
    graph_custom_power <- readRDS(test_path("data", "graph_custom_power.rds"))

    expect_no_error(
        plot(graph_custom_power)
    )

    # there are some non-deterministic elements in the plot objects themselves,
    # so we cannot compare them (mainly due to {igraph}). `equivalent_ggplot2()`
    # writes them to an .svg file and then compares the md5sums
    expect_true(
        equivalent_ggplot2(
            plot(graph_custom_power),
            autoplot(graph_custom_power)
        )
    )

    skip_on_ci()
    vdiffr::expect_doppelganger(
        title = "Graph custom power with plot()",
        plot(graph_custom_power)
    )
})

test_that("plot and autoplot with random graph", {
    set.seed(1)
    random_graph <- graph_optimal_random(5)

    expect_no_error(
        plot(random_graph)
    )

    expect_true(
        equivalent_ggplot2(
            plot(random_graph),
            autoplot(random_graph)
        )
    )

    skip_on_ci()
    vdiffr::expect_doppelganger(
        title = "Random graph with plot",
        plot(random_graph)
    )
})

test_that("plot and autoplot with user-supplied title", {
    set.seed(1)
    random_graph <- graph_optimal_random(6)

    expect_no_error(
        plot(
            random_graph,
            title = "Random 6-hypotheses graph"
        )
    )

    expect_true(
        equivalent_ggplot2(
            plot(
                random_graph,
                title = "Random 6-hypotheses graph"
            ),
            autoplot(
                random_graph,
                title = "Random 6-hypotheses graph"
            )
        )
    )

    skip_on_ci()
    vdiffr::expect_doppelganger(
        title = "Random graph and title with plot",
        plot(
            random_graph,
            title = "Random 6-hypotheses graph"
        )
    )

    vdiffr::expect_doppelganger(
        title = "Random graph and title with autoplot",
        autoplot(
            random_graph,
            title = "Random 6-hypotheses graph"
        )
    )
})

test_that("plot with input checks", {
    set.seed(1)
    rg <- graph_optimal_random(6)

    # root must be an integerish vector
    expect_error(
        plot(rg, root = c(1.2, 3.1)),
        "`root` must be integer-like or `NULL`, not a double vector."
    )

    expect_error(
        plot(rg, root = TRUE),
        "`root` must be integer-like or `NULL`, not `TRUE`."
    )

    expect_error(
        plot(rg, digits = TRUE),
        "`digits` must be a whole number or `NULL`, not `TRUE`."
    )

    expect_error(
        plot(rg, digits = 4),
        "`digits` must be a whole number between 0 and 3 or `NULL`"
    )

    expect_error(
        plot(rg, title = TRUE),
        "`title` must be a single string or `NULL`, not `TRUE`."
    )
})
