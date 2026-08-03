test_that("str_trunc_light", {
    expect_snapshot({
        str_trunc_light(
            c(
                0.182978723405255,
                0,
                0.348936170213766,
                0.251063829788234,
                0.217021276596745
            )
        )
    })

    expect_identical(
        str_trunc_light(
            c(
                0.182978723405255,
                0,
                0.348936170213766,
                0.251063829788234,
                0.217021276596745
            )
        ),
        c(
            "0.18...",
            "0",
            "0.34...",
            "0.25...",
            "0.21..."
        )
    )

    # transformation to scientific notation should not happen
    expect_snapshot({
        str_trunc_light(
            c(
                0.182978723405255,
                0.000000000013766
            )
        )
    })

    expect_identical(
        str_trunc_light(
            c(
                0.182978723405255,
                0.000000000013766
            )
        ),
        c(
            "0.18...",
            "0.00..."
        )
    )

    # NAs are contagious
    expect_snapshot({
        str_trunc_light(
            c(
                0.182978723405255,
                0,
                NA
            )
        )
    })

    expect_identical(
        str_trunc_light(
            c(
                0.182978723405255,
                0,
                NA
            )
        ),
        c(
            "0.18...",
            "0",
            NA
        )
    )
})

test_that("pull_row & pull_rows", {
    expect_identical(
        pull_row(
            row = 2,
            x = matrix(1:9, nrow = 3)
        ),
        c(2L, 5L, 8L)
    )

    expect_identical(
        pull_rows(
            matrix(1:9, nrow = 3)
        ),
        list(
            c(1L, 4L, 7L),
            c(2L, 5L, 8L),
            c(3L, 6L, 9L)
        )
    )
})

test_that("print_offending_row", {
    expect_snapshot({
        print_offending_row(
            1,
            matrix(1:9, nrow = 3),
            type = "sum"
        )
    })

    expect_snapshot({
        print_offending_row(
            1,
            matrix(1:9, nrow = 3),
            type = "value",
            ref = 0
        )
    })

    expect_snapshot({
        print_offending_row(
            1,
            matrix(1:9, nrow = 3),
            type = "value",
            ref = 1
        )
    })

    expect_snapshot({
        print_offending_row(
            1,
            matrix(1:9, nrow = 3),
            type = "na"
        )
    })
})

test_that("offending_rows_bullets", {
    expect_snapshot({
        offending_rows_bullets(
            c(1, 2),
            matrix(1:9, nrow = 3),
            type = "sum"
        )
    })
})
