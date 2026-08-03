test_that("ga check functions", {
    f <- function(x) abs(x) + cos(x)

    fitness <- function(x) -f(x)
    res <- GA::ga(
        type = "real-valued",
        fitness = fitness,
        lower = -20,
        upper = 20,
        monitor = NULL
    )

    expect_true(is_ga(res))
    expect_no_error(check_ga(res))

    expect_no_error(check_ga(NULL, allow_null = TRUE))
    expect_error(
        check_ga("foo"),
        '`"foo"` must be a GA object, not the string "foo"'
    )
})
