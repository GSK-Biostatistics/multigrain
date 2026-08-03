# to ensure consistent cli output (prevents local vs GHA whitespace differences)
cli::start_app()
on.exit(cli::stop_app(), add = TRUE)

test_that("control_local can add and remove options", {
    ctrl <- multigrain_control()

    ctrl <- control_local(ctrl, x = 1)
    expect_identical(ctrl$local_opt, list(x = 1))

    ctrl <- control_local(ctrl, x = NULL)
    expect_identical(ctrl$local_opt, list())
})

test_that("control_local can add and remove multiple options", {
    ctrl <- multigrain_control()

    expect_snapshot(
        control_local(
            ctrl,
            algorithm = "NLOPT_LN_COBYLA",
            xtol_rel = 5e-8,
            xtol_abs = 5e-9,
            maxeval = 5000
        )
    )

    ctrl <- control_local(
        ctrl,
        algorithm = "NLOPT_LN_COBYLA",
        xtol_rel = 5e-8,
        xtol_abs = 5e-9,
        maxeval = 5000
    )

    expect_identical(
        ctrl$local_opt,
        list(
            algorithm = "NLOPT_LN_COBYLA",
            xtol_rel = 5e-8,
            xtol_abs = 5e-9,
            maxeval = 5000
        )
    )

    ctrl <- control_local(
        ctrl,
        maxeval = NULL,
        xtol_abs = NULL
    )
    expect_identical(
        ctrl$local_opt,
        list(
            algorithm = "NLOPT_LN_COBYLA",
            xtol_rel = 5e-8
        )
    )
})
