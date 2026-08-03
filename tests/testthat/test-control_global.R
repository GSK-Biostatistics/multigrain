# to ensure consistent cli output (prevents local vs GHA whitespace differences)
cli::start_app()
on.exit(cli::stop_app(), add = TRUE)

test_that("control_global can add and remove options", {
    ctrl <- multigrain_control()

    ctrl <- control_global(ctrl, x = 1)
    expect_identical(ctrl$global_opt, list(x = 1))

    ctrl <- control_global(ctrl, x = NULL)
    expect_identical(ctrl$global_opt, list())
})

test_that("control_global can add and remove multiple options", {
    ctrl <- multigrain_control()

    expect_snapshot(
        control_global(
            ctrl,
            pcrossover = 0.2,
            pmutation = 0.8,
            maxiter = 1e5,
            popSize = 200,
            run = 200,
            monitor = FALSE
        )
    )

    ctrl <- control_global(
        ctrl,
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

    ctrl <- control_global(
        ctrl,
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
