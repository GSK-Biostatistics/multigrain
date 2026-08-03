# to ensure consistent cli output (prevents local vs GHA whitespace differences)
cli::start_app()
on.exit(cli::stop_app(), add = TRUE)

test_that("default_control", {
    expect_snapshot(
        print.default(
            default_control()
        )
    )

    expect_snapshot(
        print.default(
            default_control()
        )
    )

    ctrl <- default_control()

    expect_false(ctrl$global_opt$monitor)
    expect_identical(ctrl$local_opt$print_level, 0L)

    expect_s3_class(
        default_control(),
        "multigrain_control"
    )
})

test_that("control_prepare injects the expected defaults", {
    empty_ctrl <- multigrain_control()

    test_pvals <- sample.int(100, size = 60 * 3, replace = TRUE) / 100
    test_pvals <- matrix(test_pvals, ncol = 3)

    expect_snapshot(
        control_prepare(
            empty_ctrl,
            pvals = test_pvals
        )
    )

    prep_ctrl <- control_prepare(
        empty_ctrl,
        pvals = test_pvals
    )
    expect_identical(prep_ctrl$nsim_local, 60L)
    expect_identical(prep_ctrl$nsim_global, 60L)

    expect_identical(
        prep_ctrl$local_opt,
        list(
            algorithm = "NLOPT_LN_COBYLA",
            xtol_rel = 5e-08,
            xtol_abs = 5e-09,
            maxeval = 5000,
            print_level = 0L
        )
    )

    expect_identical(
        prep_ctrl$global_opt,
        list(
            pcrossover = 0.2,
            pmutation = 0.8,
            maxiter = 1e+05,
            popSize = 200L,
            run = 200,
            monitor = FALSE,
            optimArgs = list(
                method = "Nelder-Mead",
                poptim = 0.2,
                pressel = 0.6
            )
        )
    )
})

test_that("control_prepare uses pvals for calibration", {
    empty_ctrl <- multigrain_control()
    test_pvals <- sample.int(100, size = 10000 * 3, replace = TRUE) / 100
    test_pvals <- matrix(test_pvals, ncol = 3)

    prep_ctrl <- control_prepare(empty_ctrl, pvals = test_pvals)

    expect_identical(prep_ctrl$nsim_local, 10000L)
    expect_identical(prep_ctrl$nsim_global, 10000L)
    expect_identical(prep_ctrl$global_opt$popSize, 200L)

    test_pvals <- sample.int(100, size = 60000 * 3, replace = TRUE) / 100
    test_pvals <- matrix(test_pvals, ncol = 3)

    prep_ctrl <- control_prepare(empty_ctrl, pvals = test_pvals)

    expect_identical(prep_ctrl$nsim_local, 60000L)
    expect_identical(prep_ctrl$nsim_global, 50000L)

    test_pvals <- sample.int(100, size = 60000 * 3, replace = TRUE) / 100
    test_pvals <- matrix(test_pvals, ncol = 100)

    prep_ctrl <- control_prepare(empty_ctrl, pvals = test_pvals)
    expect_identical(prep_ctrl$global_opt$popSize, 500L)
})

test_that("control_prepare complains with too large nsim values", {
    test_pvals <- sample.int(100, size = 10000 * 3, replace = TRUE) / 100
    test_pvals <- matrix(test_pvals, ncol = 3)

    # nsim_local
    ctrl_local <- multigrain_control() |>
        control_nsim_local(12000)

    suppressWarnings(
        prep_ctrl_local <- control_prepare(ctrl_local, pvals = test_pvals)
    )

    expect_identical(
        prep_ctrl_local$nsim_local,
        10000L
    )

    # nsim_global
    ctrl_global <- multigrain_control() |>
        control_nsim_global(13000)

    suppressWarnings(
        prep_ctrl_global <- control_prepare(ctrl_global, pvals = test_pvals)
    )

    expect_identical(
        prep_ctrl_global$nsim_global,
        10000L
    )
})

test_that("control_prepare with verbose", {
    empty_ctrl <- multigrain_control()

    test_pvals <- sample.int(100, size = 60 * 3, replace = TRUE) / 100
    test_pvals <- matrix(test_pvals, ncol = 3)

    prep_ctrl_quiet <- control_prepare(
        empty_ctrl,
        pvals = test_pvals
    )

    expect_identical(prep_ctrl_quiet$local_opt$print_level, 0L)
    expect_false(prep_ctrl_quiet$global_opt$monitor)

    prep_ctrl_detail <- control_prepare(
        empty_ctrl,
        pvals = test_pvals,
        verbose = "detail"
    )

    expect_identical(prep_ctrl_detail$local_opt$print_level, 1L)

    expect_true(
        prep_ctrl_detail$global_opt$monitor
    )
})

test_that("control_prepare does not overwrite user-set verbosity", {
    user_ctrl <- multigrain_control() |>
        control_local(print_level = 3) |>
        control_global(monitor = FALSE)

    test_pvals <- sample.int(100, size = 60 * 3, replace = TRUE) / 100
    test_pvals <- matrix(test_pvals, ncol = 3)

    prep_ctrl <- control_prepare(
        user_ctrl,
        pvals = test_pvals
    )

    expect_identical(
        user_ctrl$global_opt$monitor,
        prep_ctrl$global_opt$monitor
    )

    expect_identical(
        user_ctrl$local_opt$print_level,
        prep_ctrl$local_opt$print_level
    )

    # verbose should basically have no impact on user-supplied args
    prep_ctrl2 <- control_prepare(
        user_ctrl,
        pvals = test_pvals,
        verbose = "detail"
    )

    expect_identical(
        user_ctrl$global_opt$monitor,
        prep_ctrl2$global_opt$monitor
    )

    expect_identical(
        user_ctrl$local_opt$print_level,
        prep_ctrl2$local_opt$print_level
    )
})

test_that("control_prepare with control_global and optimArgs", {
    # optimArgs are a bit trickier to update since we do not want to overwrite
    # all values with the defaults

    ctrl <- multigrain_control() |>
        control_global(
            optimArgs = list(
                method = "foo",
                poptim = 0.7
            )
        )

    expect_identical(
        ctrl$global_opt$optimArgs,
        list(
            method = "foo",
            poptim = 0.7
        )
    )

    test_pvals <- sample.int(100, size = 60 * 3, replace = TRUE) / 100
    test_pvals <- matrix(test_pvals, ncol = 3)

    ctrl_prep <- control_prepare(
        ctrl,
        pvals = test_pvals
    )

    expect_identical(
        ctrl_prep$global_opt$optimArgs,
        list(
            method = "foo",
            poptim = 0.7,
            pressel = 0.6
        )
    )

    ctrl <- multigrain_control() |>
        control_global(
            optimArgs = list(
                method = "foo",
                poptim = 0.7,
                bar = "baz"
            )
        ) |>
        control_prepare(pvals = test_pvals)

    expect_identical(
        ctrl$global_opt$optimArgs,
        list(
            method = "foo",
            poptim = 0.7,
            pressel = 0.6,
            bar = "baz"
        )
    )
})
