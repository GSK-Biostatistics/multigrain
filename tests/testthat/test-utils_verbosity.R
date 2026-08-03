test_that("multigrain_verbosity() defaults to 'info'", {
    withr::local_options(list(multigrain_verbosity = NULL))
    expect_identical(multigrain_verbosity(), "info")
})

test_that("multigrain_verbosity() validates the value it finds", {
    withr::local_options(list(multigrain_verbosity = TRUE))
    expect_snapshot(error = TRUE, multigrain_verbosity())

    withr::local_options(list(multigrain_verbosity = "foo"))
    expect_snapshot(error = TRUE, multigrain_verbosity())

    withr::local_options(list(multigrain_verbosity = 1))
    expect_snapshot(error = TRUE, multigrain_verbosity())
})

test_that("multigrain_verbosity() with other values", {
    withr::local_options(list(multigrain_verbosity = "detail"))
    expect_identical(multigrain_verbosity(), "detail")

    withr::local_options(list(multigrain_verbosity = "silent"))
    expect_identical(multigrain_verbosity(), "silent")
})
