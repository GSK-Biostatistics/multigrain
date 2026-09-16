# to ensure consistent cli output (prevents local vs GHA whitespace differences)
cli::start_app()
on.exit(cli::stop_app(), add = TRUE)

test_that("control_nsim_local can override default nsim_local value", {
    ctrl <- multigrain_control()

    expect_null(ctrl$nsim_local)

    ctrl <- control_nsim_local(ctrl, 1000)
    expect_identical(ctrl$nsim_local, 1000)
})

test_that("control_nsim_global can override default nsim_global value", {
    ctrl <- multigrain_control()

    expect_null(ctrl$nsim_global)

    ctrl <- control_nsim_global(ctrl, 1000)
    expect_identical(ctrl$nsim_global, 1000)
})

test_that("adjust_nsim_local", {
    ctrl_local <- multigrain_control() |>
        control_nsim_local(12000)

    expect_warning(
        adjust_nsim_local(ctrl_local, 10000),
        "(`nsim_local`) is greater than the number of rows in `pvals`",
        fixed = TRUE
    )

    expect_snapshot_warning(
        ctrl_local <- adjust_nsim_local(ctrl_local, 10000)
    )

    expect_identical(
        ctrl_local$nsim_local,
        10000
    )
})

test_that("adjust_nsim_global", {
    ctrl_global <- multigrain_control() |>
        control_nsim_global(13000)

    expect_warning(
        adjust_nsim_global(ctrl_global, 5000),
        "(`nsim_global`) is greater than the number of rows in `pvals`",
        fixed = TRUE
    )

    expect_snapshot_warning(
        ctrl_global <- adjust_nsim_global(ctrl_global, 5000)
    )

    expect_identical(
        ctrl_global$nsim_global,
        5000
    )
})

test_that("adjust_nsim_local does not adjust an empty nsim_local", {
    empty_ctrl <- multigrain_control()

    expect_identical(
        adjust_nsim_local(empty_ctrl, 10000),
        empty_ctrl
    )
})

test_that("adjust_nsim_global does not adjust an empty nsim_global", {
    empty_ctrl <- multigrain_control()

    expect_identical(
        adjust_nsim_global(empty_ctrl, 10000),
        empty_ctrl
    )
})
