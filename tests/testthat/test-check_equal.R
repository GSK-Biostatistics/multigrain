test_that("check_equal() and check_not_equal() work", {
    expect_false(
        check_equal(c(1 - 10e-13, NA), 1)
    )

    expect_true(
        check_equal(1 - 10e-13, 1)
    )

    expect_false(
        check_equal(1 - 10e-13, 1, tolerance = 10e-14)
    )

    test_trans_constraint <- matrix_tolerance(5, tolerance = 10e-13)

    expect_true(
        check_equal(rowSums(test_trans_constraint), 1)
    )

    expect_false(
        check_equal(rowSums(test_trans_constraint), 1, tolerance = 0)
    )

    expect_true(
        check_not_equal(rowSums(test_trans_constraint), 1, tolerance = 0)
    )

    expect_true(
        check_equal(rowSums(test_trans_constraint), 1)
    )
})

test_that("vec_check_equal() and vec_check_not_equal()", {
    expect_identical(
        vec_check_equal(c(1 - 10e-13, NA), 1),
        c(TRUE, FALSE)
    )

    expect_identical(
        vec_check_not_equal(c(1 - 10e-13, NA), 1),
        c(FALSE, TRUE)
    )
})
