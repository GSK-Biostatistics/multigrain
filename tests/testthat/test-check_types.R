test_that("check_double", {
    expect_no_error(
        check_double(c(0.1, 0.2, 0.3, NA, NA))
    )

    expect_no_error(
        check_double(
            NULL,
            allow_null = TRUE
        )
    )

    expect_error(
        check_double(
            NULL,
            allow_null = FALSE
        )
    )

    expect_error(
        check_double(1:4),
        "`1:4` must be a double, not an integer vector."
    )

    expect_error(
        check_double("foo"),
        '`"foo"` must be a double, not the string "foo".'
    )

    # allow a vector of NAs
    expect_no_error(
        check_double(
            c(NA, NA),
            allow_na = TRUE
        )
    )

    expect_no_error(
        check_double(
            c(0.5, NA, NA),
            allow_na = TRUE
        )
    )

    expect_error(
        check_double(c(NA, NA)),
        "`c(NA, NA)` must be a double, not a logical vector.",
        fixed = TRUE
    )
})

test_that("check_double_matrix", {
    expect_no_error(
        check_double_matrix(
            matrix(
                c(
                    0.1, 0.2, 0.3, 0.4
                )
            )
        )
    )

    expect_no_error(
        check_double_matrix(
            NULL,
            allow_null = TRUE
        )
    )

    expect_error(
        check_double_matrix(
            NULL,
            allow_null = FALSE
        )
    )

    expect_error(
        check_double_matrix(1:4),
        "`1:4` must be a double matrix, not an integer vector."
    )

    expect_error(
        check_double_matrix("foo"),
        '`"foo"` must be a double matrix, not the string "foo".'
    )

    expect_error(
        check_double_matrix(
            matrix(
                letters[1:9],
                nrow = 3
            )
        ),
        "`matrix(letters[1:9], nrow = 3)` must be a double matrix",
        fixed = TRUE
    )
})

test_that("check_integerish", {
    # integer-ish is less strict
    expect_no_error(check_integerish(100))
    expect_no_error(check_integerish(100L))
    expect_error(
        check_integerish("foo"),
        '`"foo"` must be integer-like, not the string "foo".'
    )
    expect_error(
        check_integerish(TRUE),
        "`TRUE` must be integer-like, not `TRUE`."
    )

    # the max integer value supported in R is 2^31 - 1
    # integer-like, but cannot actually be cast or coerced to integer
    expect_no_error(check_integerish(2^31))

    expect_no_error(check_integerish(NULL, allow_null = TRUE))
})
