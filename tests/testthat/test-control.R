# to ensure consistent cli output (prevents local vs GHA whitespace differences)
cli::start_app()
on.exit(cli::stop_app(), add = TRUE)

test_that("multigrain_control", {
    control <- multigrain_control()

    expect_s3_class(
        control,
        "multigrain_control"
    )

    expect_snapshot(
        multigrain_control()
    )

    expect_named(
        control,
        c(
            "nsim_local",
            "nsim_global",
            "local_opt",
            "global_opt"
        )
    )
})

test_that("new_multigrain_control() (low-level constructor)", {
    control <- new_multigrain_control()

    expect_s3_class(
        control,
        "multigrain_control"
    )

    expect_snapshot(
        new_multigrain_control()
    )

    expect_named(
        control,
        c(
            "nsim_local",
            "nsim_global",
            "local_opt",
            "global_opt"
        )
    )
})

test_that("is_control()", {
    multigrain_ctrl <- new_multigrain_control()

    expect_true(
        is_control(multigrain_ctrl)
    )
})

test_that("check_control() gives useful error", {
    expect_snapshot(error = TRUE, {
        check_control(1)
    })
})

test_that("check_control() with correct input", {
    multigrain_ctrl <- new_multigrain_control()

    expect_snapshot(
        check_control(
            multigrain_ctrl
        )
    )
})

test_that("check_control() with allow_null", {
    expect_snapshot(
        check_control(
            NULL,
            allow_null = TRUE
        )
    )

    expect_snapshot(error = TRUE, {
        check_control(NULL, allow_null = FALSE)
    })
})

test_that("multigrain_control print method", {
    test_pvals <- sample.int(100, size = 100 * 3, replace = TRUE) / 100
    test_pvals <- matrix(test_pvals, ncol = 3)

    expect_snapshot({
        multigrain_control() |>
            control_prepare(
                pvals = test_pvals
            )
    })

    expect_null(print.multigrain_control(NULL))
})
