test_that("graph_constraint plotting with autoplot", {
    gc1 <- graph_constraint(
        c(NA, NA, 0, 0, 0),
        matrix(
            c(
                0, 0.8, 0.2, 0, 0,
                NA, 0, 0, NA, 0,
                NA, 0.8, 0, 0.1, NA,
                0.9, NA, NA, 0, NA,
                0, 0, 0, 1, 0
            ),
            nrow = 5,
            byrow = TRUE
        )
    )

    expect_no_error(
        autoplot(gc1)
    )

    expect_s3_class(
        autoplot(gc1),
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
        title = "Graph constraint autoplot",
        autoplot(gc1)
    )
})

test_that("graph_constraint plotting with plot", {
    gc1 <- graph_constraint(
        c(NA, NA, 0, 0, 0),
        matrix(
            c(
                0, 0.8, 0.2, 0, 0,
                NA, 0, 0, NA, 0,
                NA, 0.8, 0, 0.1, NA,
                0.9, NA, NA, 0, NA,
                0, 0, 0, 1, 0
            ),
            nrow = 5,
            byrow = TRUE
        )
    )

    expect_no_error(
        plot(gc1)
    )

    expect_true(
        equivalent_ggplot2(
            plot(gc1),
            autoplot(gc1)
        )
    )

    skip_on_ci()
    vdiffr::expect_doppelganger(
        title = "Graph constraint plot",
        plot(gc1)
    )
})

test_that("plot and autoplot with user-supplied title", {
    gc1 <- graph_constraint(
        c(NA, NA, 0, 0, 0),
        matrix(
            c(
                0, 0.8, 0.2, 0, 0,
                NA, 0, 0, NA, 0,
                NA, 0.8, 0, 0.1, NA,
                0.9, NA, NA, 0, NA,
                0, 0, 0, 1, 0
            ),
            nrow = 5,
            byrow = TRUE
        )
    )

    expect_no_error(
        plot(
            gc1,
            title = "Graph constraint"
        )
    )

    expect_true(
        equivalent_ggplot2(
            plot(
                gc1,
                title = "Graph constraint plot with title"
            ),
            autoplot(
                gc1,
                title = "Graph constraint plot with title"
            )
        )
    )

    skip_on_ci()
    vdiffr::expect_doppelganger(
        title = "Graph constraint title with plot",
        plot(
            gc1,
            title = "GC + title and plot"
        )
    )

    vdiffr::expect_doppelganger(
        title = "Graph constraint title with autoplot",
        autoplot(
            gc1,
            title = "GC + title and autoplot"
        )
    )
})


test_that("plot with input checks", {
    gc1 <- graph_constraint(
        c(NA, NA, 0, 0, 0),
        matrix(
            c(
                0, 0.8, 0.2, 0, 0,
                NA, 0, 0, NA, 0,
                NA, 0.8, 0, 0.1, NA,
                0.9, NA, NA, 0, NA,
                0, 0, 0, 1, 0
            ),
            nrow = 5,
            byrow = TRUE
        )
    )

    # root must be an integerish vector
    expect_error(
        plot(gc1, root = c(1.2, 3.1)),
        "`root` must be integer-like or `NULL`, not a double vector."
    )

    expect_error(
        plot(gc1, root = TRUE),
        "`root` must be integer-like or `NULL`, not `TRUE`."
    )

    expect_error(
        plot(gc1, digits = TRUE),
        "`digits` must be a whole number or `NULL`, not `TRUE`."
    )

    expect_error(
        plot(gc1, digits = 4),
        "`digits` must be a whole number between 0 and 3 or `NULL`"
    )

    expect_error(
        plot(gc1, title = TRUE),
        "`title` must be a single string or `NULL`, not `TRUE`."
    )
})
