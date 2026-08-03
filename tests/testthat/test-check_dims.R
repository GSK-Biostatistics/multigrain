test_that("check_length", {
    # does not check vector type, only length
    expect_snapshot(error = TRUE, {
        check_length("foo", 2)
    })

    expect_snapshot(error = TRUE, {
        check_length(1:3, 2)
    })

    expect_no_error(
        check_length(NULL, 2, allow_null = TRUE)
    )

    expect_no_error(
        check_length(1:3, 3)
    )
})

test_that("check_dim", {
    # does not check vector type, only dim
    expect_snapshot(error = TRUE, {
        check_dim(
            matrix(c("foo", "bar")),
            2
        )
    })

    expect_snapshot(error = TRUE, {
        check_dim(
            matrix(1:6, ncol = 2),
            2
        )
    })

    expect_no_error(
        check_dim(NULL, 2, allow_null = TRUE)
    )

    expect_no_error(
        check_dim(
            matrix(1:9, ncol = 3),
            3
        )
    )
})
