test_that("cat_control_opt", {
    expect_snapshot(
        cat_control_opt("foo", "bar")
    )

    expect_snapshot(
        cat_control_opt("foo", NULL)
    )
})

test_that("cat_optim_args", {
    ctrl <- multigrain_control() |>
        control_prepare(pvals = pvals_fixture)

    expect_snapshot({
        cat_optim_args(ctrl$global_opt$optimArgs)
    })

    expect_null(cat_optim_args(NULL))
})
