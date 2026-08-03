# to ensure consistent cli output (prevents local vs GHA whitespace differences)
cli::start_app()
on.exit(cli::stop_app(), add = TRUE)

test_that("control_nsim_local can override default nsim_local value", {
    ctrl <- multigrain_control()

    expect_null(ctrl$nsim_local)

    ctrl <- control_nsim_local(ctrl, 1000)
    expect_identical(ctrl$nsim_local, 1000)
})

test_that("control_nsim_local complains with incorrect nsim_local", {
    ctrl <- multigrain_control()

    # nsim_local must be a positive integer scalar
    expect_error(
        control_nsim_local(ctrl, "foo"),
        '`nsim_local` must be a whole number, not the string "foo".'
    )

    expect_snapshot(error = TRUE, {
        control_nsim_local(ctrl, "foo")
    })

    expect_error(
        control_nsim_local(ctrl, c(100, 200)),
        "`nsim_local` must be a whole number, not a double vector."
    )

    expect_snapshot(error = TRUE, {
        control_nsim_local(ctrl, c(100, 200))
    })
})

test_that("control_nsim_global can override default nsim_global value", {
    ctrl <- multigrain_control()

    expect_null(ctrl$nsim_global)

    ctrl <- control_nsim_global(ctrl, 1000)
    expect_identical(ctrl$nsim_global, 1000)
})

test_that("control_nsim_global complains with incorrect nsim_global", {
    ctrl <- multigrain_control()

    # nsim_local must be a positive integer scalar
    expect_error(
        control_nsim_global(ctrl, "foo"),
        '`nsim_global` must be a whole number, not the string "foo".'
    )

    expect_snapshot(error = TRUE, {
        control_nsim_global(ctrl, "foo")
    })

    expect_error(
        control_nsim_global(ctrl, c(100, 200)),
        "`nsim_global` must be a whole number, not a double vector."
    )

    expect_snapshot(error = TRUE, {
        control_nsim_global(ctrl, c(100, 200))
    })
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
