test_that("control_global can add and remove options", {
    ctrl <- multigrain_control()
    ctrl <- ctrl |> control_global(x = 1)

    expect_identical(ctrl$global_opt, list(x = 1))

    ctrl <- ctrl |> control_global(x = NULL)

    expect_identical(ctrl$global_opt, list())
})

test_that("control_global can add and remove multiple options", {
    ctrl <- multigrain_control()
    ctrl <- ctrl |>
        control_global(
            pcrossover = 0.2,
            pmutation = 0.8,
            maxiter = 1e5,
            popSize = 200,
            run = 200,
            monitor = FALSE
        )

    expect_identical(
        ctrl$global_opt,
        list(
            pcrossover = 0.2,
            pmutation = 0.8,
            maxiter = 1e+05,
            popSize = 200,
            run = 200,
            monitor = FALSE
        )
    )

    ctrl <- ctrl |>
        control_global(
            pcrossover = 0.2,
            pmutation = 0.8,
            maxiter = 1e5,
            popSize = 200,
            run = 200,
            monitor = FALSE
        )

    expect_identical(
        ctrl$global_opt,
        list(
            pcrossover = 0.2,
            pmutation = 0.8,
            maxiter = 1e5,
            popSize = 200,
            run = 200,
            monitor = FALSE
        )
    )

    ctrl <- ctrl |>
        control_global(
            run = NULL,
            monitor = TRUE
        )

    expect_identical(
        ctrl$global_opt,
        list(
            pcrossover = 0.2,
            pmutation = 0.8,
            maxiter = 1e5,
            popSize = 200,
            monitor = TRUE
        )
    )
})
