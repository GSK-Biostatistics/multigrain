# to ensure consistent cli output (prevents local vs GHA whitespace differences)
cli::start_app()
on.exit(cli::stop_app(), add = TRUE)

test_that("assert_trans_constr() works", {
    # error when complete rows do not add up to 1 - row 1
    expect_snapshot(error = TRUE, {
        assert_trans_constr(
            new_graph_constraint(
                c(NA, NA, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.1, 0,
                        NA, 0, 0, NA,
                        0, 1, 0, 0,
                        1, 0, 0, 0
                    ),
                    nrow = 4,
                    byrow = TRUE
                )
            )
        )
    })

    # complains with non-zero values on the diagonal
    expect_snapshot(error = TRUE, {
        assert_trans_constr(
            new_graph_constraint(
                c(NA, NA, 0, 0),
                matrix(
                    c(
                        0, 0.8, NA, NA,
                        NA, 0.1, 0, NA,
                        0, 1, 0, 0,
                        1, 0, 0, 0
                    ),
                    nrow = 4,
                    byrow = TRUE
                )
            )
        )
    })

    # error: trans_constraint is not square
    expect_snapshot(error = TRUE, {
        assert_trans_constr(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.2, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.9, 0, 0.1, 0,
                        1, 0, 0, 0, 0
                    ),
                    nrow = 4,
                    byrow = TRUE
                )
            )
        )
    })

    # error: row 4 has a value greater than 1
    expect_snapshot(error = TRUE, {
        assert_trans_constr(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.2, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.9, 0, 0.1, 0,
                        0, 0, 0, 0, 1.7,
                        1, 0, 0, 0, 0
                    ),
                    nrow = 5,
                    byrow = TRUE
                )
            )
        )
    })

    # error: row 4 has a single NA value
    expect_snapshot(error = TRUE, {
        assert_trans_constr(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.2, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.9, 0, 0.1, 0,
                        0, 0, NA, 0, 0.9,
                        1, 0, 0, 0, 0
                    ),
                    nrow = 5,
                    byrow = TRUE
                )
            )
        )
    })

    # error when complete rows do not add up to 1 - row 1
    expect_snapshot(error = TRUE, {
        assert_trans_constr(
            new_graph_constraint(
                c(NA, NA, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.1, 0,
                        NA, 0, 0, NA,
                        0, 1, 0, 0,
                        1, 0, 0, 0
                    ),
                    nrow = 4,
                    byrow = TRUE
                )
            )
        )
    })

    # it does nothing when valid
    expect_no_error(
        assert_trans_constr(
            new_graph_constraint(
                c(NA, NA, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.1, 0.1,
                        NA, 0, 0, NA,
                        0, 1, 0, 0,
                        1, 0, 0, 0
                    ),
                    nrow = 4,
                    byrow = TRUE
                )
            )
        )
    )
})

test_that("assert_trans_constr_diagonal() works", {
    # complains with non-zero values on the diagonal
    expect_snapshot(error = TRUE, {
        assert_trans_constr_diagonal(
            new_graph_constraint(
                c(NA, NA, 0, 0),
                matrix(
                    c(
                        0, 0.8, NA, NA,
                        NA, 0.1, 0, NA,
                        0, 1, 0, 0,
                        1, 0, 0, 0
                    ),
                    nrow = 4,
                    byrow = TRUE
                )
            )
        )
    })

    # complains with NAs on the diagonal
    expect_snapshot(error = TRUE, {
        assert_trans_constr_diagonal(
            new_graph_constraint(
                c(NA, NA, 0, 0),
                matrix(
                    c(
                        0, 0.8, NA, NA,
                        NA, NA, 0, NA,
                        0, 1, 0, 0,
                        1, 0, 0, 0
                    ),
                    nrow = 4,
                    byrow = TRUE
                )
            )
        )
    })

    # it does nothing when all the values on the diagonal are 0
    expect_no_error(
        assert_trans_constr_diagonal(
            new_graph_constraint(
                c(NA, NA, 0, 0),
                matrix(
                    c(
                        0, 0.8, NA, NA,
                        NA, 0, 0, NA,
                        0, 1, 0, 0,
                        1, 0, 0, 0
                    ),
                    nrow = 4,
                    byrow = TRUE
                )
            )
        )
    )
})

test_that("assert_trans_constr_square() works", {
    # no error: trans_constraint is square
    expect_no_error(
        assert_trans_constr_square(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.2, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.9, 0, 0.1, 0,
                        0.9, NA, NA, 0, 0,
                        1, 0, 0, 0, 0
                    ),
                    nrow = 5,
                    byrow = TRUE
                )
            )
        )
    )

    # error: tc is not square
    expect_snapshot(error = TRUE, {
        assert_trans_constr_square(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.2, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.9, 0, 0.1, 0,
                        1, 0, 0, 0, 0
                    ),
                    nrow = 4,
                    byrow = TRUE
                )
            )
        )
    })
})

test_that("assert_trans_constr_values() works", {
    # error: row 4 has a value greater than 1
    expect_snapshot(error = TRUE, {
        assert_trans_constr_values(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.2, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.9, 0, 0.1, 0,
                        0, 0, 0, 0, 1.7,
                        1, 0, 0, 0, 0
                    ),
                    nrow = 5,
                    byrow = TRUE
                )
            )
        )
    })

    # error: row 4 has a value less than 0
    expect_snapshot(error = TRUE, {
        assert_trans_constr_values(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.2, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.9, 0, 0.1, 0,
                        0, NA, NA, 0, -0.1,
                        1, 0, 0, 0, 0
                    ),
                    nrow = 5,
                    byrow = TRUE
                )
            )
        )
    })

    # all good - nothing flagged
    expect_no_error(
        assert_trans_constr_values(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.2, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.9, 0, 0.1, 0,
                        0.3, NA, NA, 0, 0,
                        1, 0, 0, 0, 0
                    ),
                    nrow = 5,
                    byrow = TRUE
                )
            )
        )
    )
})

test_that("assert_trans_constr_row_na() works", {
    # error: row 4 has a single NA value
    expect_snapshot(error = TRUE, {
        assert_trans_constr_row_na(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.2, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.9, 0, 0.1, 0,
                        0, 0, NA, 0, 0.9,
                        1, 0, 0, 0, 0
                    ),
                    nrow = 5,
                    byrow = TRUE
                )
            )
        )
    })

    # all good - nothing flagged
    expect_no_error(
        assert_trans_constr_row_na(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.2, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.9, 0, 0.1, 0,
                        0.3, NA, NA, 0, 0,
                        1, 0, 0, 0, 0
                    ),
                    nrow = 5,
                    byrow = TRUE
                )
            )
        )
    )
})

test_that("assert_trans_constr_row_sums() works", {
    # error when complete rows do not add up to 1 - row 1
    expect_snapshot(error = TRUE, {
        assert_trans_constr_row_sums(
            new_graph_constraint(
                c(NA, NA, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.1, 0,
                        NA, 0, 0, NA,
                        0, 1, 0, 0,
                        1, 0, 0, 0
                    ),
                    nrow = 4,
                    byrow = TRUE
                )
            )
        )
    })

    # error when complete rows have a sum greater than 1
    expect_snapshot(error = TRUE, {
        assert_trans_constr_row_sums(
            new_graph_constraint(
                c(NA, NA, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.2, 0,
                        0.2, 0, 0, 0.9,
                        0, 1, 0, 0,
                        1, 0, 0, 0
                    ),
                    nrow = 4,
                    byrow = TRUE
                )
            )
        )
    })

    # error when incomplete rows have a sum greater than 1
    expect_snapshot(error = TRUE, {
        assert_trans_constr_row_sums(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.2, 0, 0,
                        0.2, 0, NA, 0.9, NA,
                        0, 1, 0, 0, 0,
                        1, 0, 0, 0, 0,
                        1, 0, 0, 0, 0
                    ),
                    nrow = 5,
                    byrow = TRUE
                )
            )
        )
    })

    # error when incomplete rows add up to 1
    expect_snapshot(error = TRUE, {
        assert_trans_constr_row_sums(
            new_graph_constraint(
                c(NA, NA, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.1, 0.1,
                        NA, 0, NA, 1,
                        0, 1, 0, 0,
                        1, 0, 0, 0
                    ),
                    nrow = 4,
                    byrow = TRUE
                )
            )
        )
    })

    # it does nothing when complete rows add up to 1
    expect_no_error(
        assert_trans_constr_row_sums(
            new_graph_constraint(
                c(NA, NA, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.1, 0.1,
                        NA, 0, 0, NA,
                        0, 1, 0, 0,
                        1, 0, 0, 0
                    ),
                    nrow = 4,
                    byrow = TRUE
                )
            )
        )
    )
})

test_that("diagnose_trans_constr_diagonal() works", {
    # success: all elements on diagonal are zero
    expect_snapshot(
        diagnose_trans_constr_diagonal(
            graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.2, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.9, 0, 0.1, 0,
                        0.8, NA, NA, 0, 0,
                        1, 0, 0, 0, 0
                    ),
                    nrow = 5,
                    byrow = TRUE
                )
            )
        )
    )

    # danger: some elements on the diagonal are not zero
    expect_snapshot(
        diagnose_trans_constr_diagonal(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.2, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.9, 0, 0.1, 0,
                        1, 0, 0, 0, 0,
                        0.9, 0, 0, 0, 0.1
                    ),
                    nrow = 5,
                    byrow = TRUE
                )
            )
        )
    )

    # danger: some elements on the diagonal are NA
    expect_snapshot(
        diagnose_trans_constr_diagonal(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.2, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.9, 0, 0.1, 0,
                        0.8, 0, NA, NA, 0,
                        1, 0, 0, 0, 0
                    ),
                    nrow = 5,
                    byrow = TRUE
                )
            )
        )
    )
})

test_that("diagnose_trans_constr_square() works", {
    # success: trans_constraint is square
    expect_snapshot(
        diagnose_trans_constr_square(
            graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.2, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.9, 0, 0.1, 0,
                        0.8, NA, NA, 0, 0,
                        1, 0, 0, 0, 0
                    ),
                    nrow = 5,
                    byrow = TRUE
                )
            )
        )
    )

    # danger: trans_constraint is not square
    expect_snapshot(
        diagnose_trans_constr_square(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.2, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.9, 0, 0.1, 0,
                        1, 0, 0, 0, 0
                    ),
                    nrow = 4,
                    byrow = TRUE
                )
            )
        )
    )
})

test_that("diagnose_trans_constr_row_sums() works", {
    # danger: row 1 is complete and sum is 0.9 (less than 1)
    expect_snapshot(
        diagnose_trans_constr_row_sums(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.1, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.8, 0, 0.2, 0,
                        0.9, 0, 0.1, 0, 0,
                        1, 0, 0, 0, 0
                    ),
                    nrow = 5,
                    byrow = TRUE
                )
            )
        )
    )

    # danger: row 3 is complete and sum is 1.1 (greater than 1)
    expect_snapshot(
        diagnose_trans_constr_row_sums(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.2, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.9, 0, 0.2, 0,
                        0.9, NA, NA, 0, 0,
                        1, 0, 0, 0, 0
                    ),
                    nrow = 5,
                    byrow = TRUE
                )
            )
        )
    )

    # danger: row 4 is incomplete and sum is 1.1 (greater than 1)
    expect_snapshot(
        diagnose_trans_constr_row_sums(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.2, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.9, 0, 0.1, 0,
                        0.9, 0.2, NA, 0, NA,
                        1, 0, 0, 0, 0
                    ),
                    nrow = 5,
                    byrow = TRUE
                )
            )
        )
    )

    # danger: row 4 is incomplete and sum is 1
    expect_snapshot(
        diagnose_trans_constr_row_sums(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.2, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.9, 0, 0.1, 0,
                        0.3, NA, NA, 0, 0.7,
                        1, 0, 0, 0, 0
                    ),
                    nrow = 5,
                    byrow = TRUE
                )
            )
        )
    )

    # all good - nothing flagged
    expect_snapshot(
        diagnose_trans_constr_row_sums(
            graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.2, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.9, 0, 0.1, 0,
                        0.3, NA, NA, 0, 0.6,
                        1, 0, 0, 0, 0
                    ),
                    nrow = 5,
                    byrow = TRUE
                )
            )
        )
    )
})

test_that("diagnose_trans_constr_gt1() works as expected", {
    # danger: row 4 has a value greater than 1 and its sum is also greater
    # than 1
    expect_snapshot(
        diagnose_trans_constr_gt1(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.2, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.9, 0, 0.1, 0,
                        1.1, NA, NA, 0, 0,
                        1, 0, 0, 0, 0
                    ),
                    nrow = 5,
                    byrow = TRUE
                )
            )
        )
    )

    # no values greater than 1 in the transition matrix
    expect_snapshot(
        diagnose_trans_constr_gt1(
            graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.2, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.9, 0, 0.1, 0,
                        0.3, 0, NA, 0, NA,
                        1, 0, 0, 0, 0
                    ),
                    nrow = 5,
                    byrow = TRUE
                )
            )
        )
    )
})

test_that("diagnose_trans_constr_lt0() works", {
    # danger: row 4 has a value less than 0
    expect_snapshot(
        diagnose_trans_constr_lt0(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.2, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.9, 0, 0.1, 0,
                        -0.1, NA, NA, 0, 0,
                        1, 0, 0, 0, 0
                    ),
                    nrow = 5,
                    byrow = TRUE
                )
            )
        )
    )

    # all good - nothing flagged
    expect_snapshot(
        diagnose_trans_constr_lt0(
            graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.2, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.9, 0, 0.1, 0,
                        0.3, NA, NA, 0, 0,
                        1, 0, 0, 0, 0
                    ),
                    nrow = 5,
                    byrow = TRUE
                )
            )
        )
    )
})

test_that("diagnose_trans_constr_row_na() works", {
    # danger: rows 1 and 4 are incomplete with a single value to optimise and
    # should be flagged
    expect_snapshot(
        diagnose_trans_constr_row_na(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, NA, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.8, 0, 0.2, 0,
                        0.9, 0, NA, 0, 0,
                        1, 0, 0, 0, 0
                    ),
                    nrow = 5,
                    byrow = TRUE
                )
            )
        )
    )

    # all good - nothing to flag. the trans_constraint rows are either complete
    # or with at least 2 values to optimise
    expect_snapshot(
        diagnose_trans_constr_row_na(
            graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.2, 0, 0,
                        NA, 0, 0, NA, NA,
                        0, 0.9, 0, 0.1, 0,
                        0.9, 0, NA, 0, NA,
                        1, 0, 0, 0, 0
                    ),
                    nrow = 5,
                    byrow = TRUE
                )
            )
        )
    )
})

test_that("diagnose_trans_constr() works", {
    # reports rows not adding up to 1 - row 1
    expect_snapshot(
        diagnose_trans_constr(
            new_graph_constraint(
                c(NA, NA, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.1, 0,
                        NA, 0, 0, NA,
                        0, 1, 0, 0,
                        1, 0, 0, 0
                    ),
                    nrow = 4,
                    byrow = TRUE
                )
            )
        )
    )

    # reports non-zero values on the diagonal
    expect_snapshot(
        diagnose_trans_constr(
            new_graph_constraint(
                c(NA, NA, 0, 0),
                matrix(
                    c(
                        0, 0.8, NA, NA,
                        NA, 0.1, 0, NA,
                        0, 1, 0, 0,
                        1, 0, 0, 0
                    ),
                    nrow = 4,
                    byrow = TRUE
                )
            )
        )
    )

    # reports trans_constraint not being square
    expect_snapshot(
        diagnose_trans_constr(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.2, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.9, 0, 0.1, 0,
                        1, 0, 0, 0, 0
                    ),
                    nrow = 4,
                    byrow = TRUE
                )
            )
        )
    )

    # reports row 4 having a value greater than 1 and its sum being greater
    # than 1
    expect_snapshot(
        diagnose_trans_constr(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.2, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.9, 0, 0.1, 0,
                        0, 0, 0, 0, 1.7,
                        1, 0, 0, 0, 0
                    ),
                    nrow = 5,
                    byrow = TRUE
                )
            )
        )
    )

    # reports row 4 having a single NA value
    expect_snapshot(
        diagnose_trans_constr(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.2, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.9, 0, 0.1, 0,
                        0, 0, NA, 0, 0.9,
                        1, 0, 0, 0, 0
                    ),
                    nrow = 5,
                    byrow = TRUE
                )
            )
        )
    )

    # reports complete rows not adding up to 1 - row 1
    expect_snapshot(
        diagnose_trans_constr(
            new_graph_constraint(
                c(NA, NA, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.1, 0,
                        NA, 0, 0, NA,
                        0, 1, 0, 0,
                        1, 0, 0, 0
                    ),
                    nrow = 4,
                    byrow = TRUE
                )
            )
        )
    )

    # it does nothing when valid
    expect_snapshot(
        diagnose_trans_constr(
            new_graph_constraint(
                c(NA, NA, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.1, 0.1,
                        NA, 0, 0, NA,
                        0, 1, 0, 0,
                        1, 0, 0, 0
                    ),
                    nrow = 4,
                    byrow = TRUE
                )
            )
        )
    )

    # reports danger an all trans_constraint checks
    expect_snapshot(
        diagnose_trans_constr(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.3, 0, 0,
                        NA, 0, -0.1, NA, 1.2,
                        0, 0.9, NA, 0.1, 0,
                        1, 0, 0, NA, NA
                    ),
                    nrow = 4,
                    byrow = TRUE
                )
            )
        )
    )
})

test_that("validate_graph_constraint() picks up transition matrix issues", {
    # fails when complete rows do not add up to 1 - row 1
    expect_snapshot(error = TRUE, {
        validate_graph_constraint(
            new_graph_constraint(
                c(NA, NA, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.1, 0,
                        NA, 0, 0, NA,
                        0, 1, 0, 0,
                        1, 0, 0, 0
                    ),
                    nrow = 4,
                    byrow = TRUE
                )
            ),
            diagnose = TRUE
        )
    })

    # it does nothing when valid
    expect_snapshot(
        validate_graph_constraint(
            new_graph_constraint(
                c(NA, NA, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.1, 0.1,
                        NA, 0, 0, NA,
                        0, 1, 0, 0,
                        1, 0, 0, 0
                    ),
                    nrow = 4,
                    byrow = TRUE
                )
            )
        )
    )
})

test_that("assert_trans_constr_values() with tolerance", {
    # row 4 has a value greater than 1, but it is within tolerance
    expect_snapshot(
        assert_trans_constr_values(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.2, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.9, 0, 0.1, 0,
                        0, 0, 0, 0, 1 + 10e-13,
                        1, 0, 0, 0, 0
                    ),
                    nrow = 5,
                    byrow = TRUE
                )
            )
        )
    )

    # row 4 has a value greater than 1, but it is within tolerance
    # which errors when reducing the tolerance
    expect_snapshot(error = TRUE, {
        assert_trans_constr_values(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.2, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.9, 0, 0.1, 0,
                        0, 0, 0, 0, 1 + 10e-13,
                        1, 0, 0, 0, 0
                    ),
                    nrow = 5,
                    byrow = TRUE
                ),
                tolerance = 0
            )
        )
    })
})

test_that("trans_constraint row sum validation with small differences", {
    # build a trans_constraint with small floating point differences
    test_trans_constraint <- matrix_tolerance(5, tolerance = 10e-13)

    # row sums are not identical to 1
    expect_false(
        identical(
            rowSums(test_trans_constraint),
            1
        )
    )

    # but they are equal to 1 (with the allowed tolerance)
    expect_equal(
        rowSums(test_trans_constraint),
        rep(1, 5),
        tolerance = 10e-12
    )

    # graph_constraint() should work (be tolerant of this small difference)
    # complete rows and sum is 1 (with tolerance)
    expect_no_error(
        graph_constraint(
            trans_constraint = test_trans_constraint
        )
    )

    expect_no_error(
        assert_trans_constr_row_sums(
            new_graph_constraint(
                trans_constraint = test_trans_constraint
            )
        )
    )

    expect_snapshot(
        diagnose_trans_constr_row_sums(
            new_graph_constraint(
                trans_constraint = test_trans_constraint
            )
        )
    )

    # test with one row equal to 1 (even with tolerance 0)
    test_trans_constraint[1, ] <- c(0, 0, 0, 0, 1)

    # diagnose will pick up on the other rows
    expect_snapshot(
        diagnose_trans_constr_row_sums(
            new_graph_constraint(
                trans_constraint = test_trans_constraint,
                tolerance = 0
            )
        )
    )

    # graph_constraint() fails with tolerance = 0
    expect_snapshot(error = TRUE, {
        graph_constraint(
            trans_constraint = test_trans_constraint,
            tolerance = 0
        )
    })

    expect_snapshot(error = TRUE, {
        assert_trans_constr_row_sums(
            new_graph_constraint(
                trans_constraint = test_trans_constraint,
                tolerance = 0
            )
        )
    })

    expect_snapshot(
        diagnose_trans_constr_row_sums(
            new_graph_constraint(
                trans_constraint = test_trans_constraint,
                tolerance = 0
            )
        )
    )
})

test_that("incomplete trans_constraint rows with tolerance", {
    # tc is a matrix where the rows have 2 NAs, but they add up to 1
    tc <- matrix_tolerance_na(5)

    # confirm row sums are not identical to 1
    expect_false(
        identical(
            rowSums(tc, na.rm = TRUE),
            1
        )
    )

    # but they are equal to 1 (with the allowed tolerance)
    expect_equal(
        rowSums(tc, na.rm = TRUE),
        rep(1, 5),
        tolerance = 10e-12
    )

    # we expect an error indicating an incomplete transition row with sum 1
    # but no error indicating the sum to be greater than 1
    expect_snapshot(error = TRUE, {
        assert_trans_constr_row_sums(
            new_graph_constraint(
                trans_constraint = tc
            )
        )
    })

    # the above is clearer with the diagnosis
    expect_snapshot(
        diagnose_trans_constr_row_sums(
            new_graph_constraint(
                trans_constraint = tc
            )
        )
    )

    # with tolerance = 0 the validation fails with sum greater than 1
    expect_snapshot(error = TRUE, {
        assert_trans_constr_row_sums(
            new_graph_constraint(
                trans_constraint = tc,
                tolerance = 0
            )
        )
    })
})

test_that("diagnose_trans_constr_gt1() with tolerance", {
    # row 4 has a value greater than 1, but is is within tolerance so it should
    # not be flagged
    expect_snapshot(
        diagnose_trans_constr_gt1(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.2, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.9, 0, 0.1, 0,
                        1 + 10e-13, NA, NA, 0, 0,
                        1, 0, 0, 0, 0
                    ),
                    nrow = 5,
                    byrow = TRUE
                )
            )
        )
    )

    # lowering the tolerance results in the value being flagged
    expect_snapshot(
        diagnose_trans_constr_gt1(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.2, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.9, 0, 0.1, 0,
                        1 + 10e-13, NA, NA, 0, 0,
                        1, 0, 0, 0, 0
                    ),
                    nrow = 5,
                    byrow = TRUE
                ),
                tolerance = 0
            )
        )
    )
})
