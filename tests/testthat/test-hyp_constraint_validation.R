# to ensure consistent cli output (prevents local vs GHA whitespace differences)
cli::start_app()
on.exit(cli::stop_app(), add = TRUE)

test_that("assert_hyp_constraint_sum() works", {
    # no error when hw is incomplete, the sum of the weights is less than 1 and
    # there are enough weights to optimise
    expect_no_error(
        assert_hyp_constraint_sum(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0)
            )
        )
    )

    # it does nothing when the hw is complete and sum of the weights adds to 1
    expect_no_error(
        assert_hyp_constraint_sum(
            new_graph_constraint(
                c(1, 0, 0, 0, 0)
            )
        )
    )

    expect_no_error(
        assert_hyp_constraint_sum(
            new_graph_constraint(
                c(0.8, 0.2, 0, 0, 0)
            )
        )
    )

    # it complains when hc is complete and the weights add up to more than 1
    expect_snapshot(error = TRUE, {
        assert_hyp_constraint_sum(
            new_graph_constraint(
                c(0.8, 0.3, 0, 0, 0)
            )
        )
    })

    # it complains when hc is complete and the weights do not add up to 1
    expect_snapshot(error = TRUE, {
        assert_hyp_constraint_sum(
            new_graph_constraint(
                c(0.7, 0.2, 0, 0, 0)
            )
        )
    })

    # it complains when hc is incomplete and the weights add up to 1
    expect_snapshot(error = TRUE, {
        assert_hyp_constraint_sum(
            new_graph_constraint(
                c(0.7, 0.3, NA, NA, 0)
            )
        )
    })

    # the check function returns the input invisibly
    a <- new_graph_constraint(c(NA, NA, 0, 0, 0))

    expect_identical(
        suppressMessages(assert_hyp_constraint_sum(a)),
        a
    )
})

test_that("sum(hc) == 1 is tolerant of small floating point differences", {
    # make the sum(hc) == 1 condition less harsh
    a <- c(4.0003, 6.511, 2.3333, 0.003, 2.3)
    tolerance <- 10e-13
    b <- (a / sum(a)) + tolerance

    # sum of b is not identical to 1
    expect_false(sum(b) == 1)

    # but it is equal to it (with the allowed tolerance)
    # nolint start: expect_identical_linter
    expect_equal(
        sum(b),
        1
    )
    # nolint end

    # hyp_constraint() should work (be tolerant of this small difference)
    expect_no_error(
        graph_constraint(b)
    )

    expect_no_error(
        assert_hyp_constraint_sum(
            new_graph_constraint(b)
        )
    )

    expect_snapshot(
        diagnose_hyp_constraint_sum(
            new_graph_constraint(b)
        )
    )

    # graph_constraint fails with tolerance = 0
    expect_snapshot(error = TRUE, {
        graph_constraint(
            b,
            tolerance = 0
        )
    })

    expect_snapshot(error = TRUE, {
        assert_hyp_constraint_sum(
            new_graph_constraint(
                b,
                tolerance = 0
            )
        )
    })

    expect_snapshot(
        diagnose_hyp_constraint_sum(
            new_graph_constraint(b, tolerance = 0)
        )
    )
})

test_that("assert_hyp_constraint_values() works", {
    # no error - no values greater than 1
    expect_no_error(
        assert_hyp_constraint_values(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0)
            )
        )
    )

    # error - hc value greater than 1
    expect_snapshot(error = TRUE, {
        assert_hyp_constraint_values(
            new_graph_constraint(
                c(NA, NA, 0, 0.2, 1.9)
            )
        )
    })

    # error - hc value less than 0
    expect_snapshot(error = TRUE, {
        assert_hyp_constraint_values(
            new_graph_constraint(
                c(NA, NA, 0, 0.2, -0.1)
            )
        )
    })

    # the check function returns the input invisibly
    a <- new_graph_constraint(c(NA, NA, 0, 0, 0))

    expect_identical(
        suppressMessages(assert_hyp_constraint_values(a)),
        a
    )
})

test_that("assert_hyp_constraint_na() works", {
    expect_no_error(
        assert_hyp_constraint_na(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0)
            )
        )
    )

    # error picking up on a single NA value in hyp_constraint
    expect_snapshot(error = TRUE, {
        assert_hyp_constraint_na(
            new_graph_constraint(
                c(NA, 0, 0, 0.2, 0)
            )
        )
    })

    # the check function returns the input invisibly
    a <- new_graph_constraint(c(NA, NA, 0, 0, 0))

    expect_identical(
        suppressMessages(assert_hyp_constraint_na(a)),
        a
    )
})

test_that("assert_hyp_constraint() works", {
    # no error when hw is incomplete, the sum of the weights is less than 1 and
    # there are enough weights to optimise
    expect_no_error(
        assert_hyp_constraint(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0)
            )
        )
    )

    # it does nothing when the hw is complete and sum of the weights adds to 1
    expect_no_error(
        assert_hyp_constraint(
            new_graph_constraint(
                c(1, 0, 0, 0, 0)
            )
        )
    )

    # error when hw is complete and the weights add up to more than 1
    expect_snapshot(error = TRUE, {
        assert_hyp_constraint(
            new_graph_constraint(
                c(0.8, 0.3, 0, 0, 0)
            )
        )
    })

    # error - hw value greater than 1
    expect_snapshot(error = TRUE, {
        assert_hyp_constraint(
            new_graph_constraint(
                c(NA, NA, 0, 0.2, 1.9)
            )
        )
    })

    # error picking up on a single NA value in hyp_constraints
    expect_snapshot(error = TRUE, {
        assert_hyp_constraint(
            new_graph_constraint(
                c(NA, 0, 0, 0.2, 0)
            )
        )
    })

    # the check function returns the input invisibly
    a <- new_graph_constraint(c(NA, NA, 0, 0, 0))

    expect_identical(
        suppressMessages(assert_hyp_constraint(a)),
        a
    )
})

test_that("assert_hyp_constraint() works with empty vector", {
    expect_no_error(
        assert_hyp_constraint(
            new_graph_constraint()
        )
    )
})

test_that("diagnose_hyp_constraint_sum() works", {
    # incomplete hyp_constraints - nothing flagged - sum(hyp_constraints) <= 1
    expect_snapshot(
        diagnose_hyp_constraint_sum(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0)
            )
        )
    )

    # incomplete hyp_constraints - sum(hyp_constraints) > 1 is flagged
    expect_snapshot(
        diagnose_hyp_constraint_sum(
            new_graph_constraint(
                c(NA, NA, 0, 0.2, 0.9)
            )
        )
    )

    # incomplete hyp_constraints & sum(hyp_constraints) == 1 is flagged
    expect_snapshot(
        diagnose_hyp_constraint_sum(
            new_graph_constraint(
                c(NA, NA, 0, 0.2, 0.8)
            )
        )
    )

    # complete hyp_constraints & sum(hyp_constraints) != 1 is flagged
    expect_snapshot(
        diagnose_hyp_constraint_sum(
            new_graph_constraint(
                c(0, 0, 0, 0.2, 0.5)
            )
        )
    )

    # complete hyp_constraints & sum(hyp_constraints) == 1 nothing flagged
    expect_snapshot(
        diagnose_hyp_constraint_sum(
            new_graph_constraint(
                c(0, 0.2, 0.1, 0.2, 0.5)
            )
        )
    )

    # the check function returns the input invisibly
    a <- new_graph_constraint(c(0, 0.2, 0.1, 0.2, 0.5))

    expect_identical(
        suppressMessages(diagnose_hyp_constraint_sum(a)),
        a
    )
})

test_that("diagnose_hyp_constraint_gt1() works", {
    # nothing flagged - all values less than 1
    expect_snapshot(
        diagnose_hyp_constraint_gt1(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0)
            )
        )
    )

    # values greater than 1 are flagged
    expect_snapshot(
        diagnose_hyp_constraint_gt1(
            new_graph_constraint(
                c(NA, NA, 0, 0.2, 1.9)
            )
        )
    )

    # the check function returns the input invisibly
    a <- new_graph_constraint(c(0, 0.2, 0.1, 0.2, 0.5))

    expect_identical(
        suppressMessages(diagnose_hyp_constraint_gt1(a)),
        a
    )
})

test_that("diagnose_hyp_constraint_lt0() works", {
    # success - no values less than 0
    expect_snapshot(
        diagnose_hyp_constraint_lt0(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0)
            )
        )
    )

    # hyp_constraints value less than 0 is flagged
    expect_snapshot(
        diagnose_hyp_constraint_lt0(
            new_graph_constraint(
                c(NA, NA, 0, 0.2, -0.1)
            )
        )
    )

    # the check function returns the input invisibly
    a <- new_graph_constraint(c(0, 0.2, 0.1, 0.2, 0.5))

    expect_identical(
        suppressMessages(diagnose_hyp_constraint_lt0(a)),
        a
    )
})

test_that("diagnose_hyp_constraint_na() works", {
    # success - at least 2 NAs in hyp_constraints
    expect_snapshot(
        diagnose_hyp_constraint_na(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0)
            )
        )
    )

    # a single NA hyp_constraints value is flagged
    expect_snapshot(
        diagnose_hyp_constraint_na(
            new_graph_constraint(
                c(NA, 0, 0, 0.2, 0.7)
            )
        )
    )

    # no NAs in hyp_constraints is flagged as inform
    expect_snapshot(
        diagnose_hyp_constraint_na(
            new_graph_constraint(
                c(0.1, 0, 0, 0.2, 0.7)
            )
        )
    )

    # the check function returns the input invisibly
    a <- new_graph_constraint(c(0, 0.2, 0.1, 0.2, 0.5))

    expect_identical(
        suppressMessages(diagnose_hyp_constraint_na(a)),
        a
    )
})

test_that("diagnose_hyp_constraint() works", {
    # success: nothing flagged
    expect_snapshot(
        diagnose_hyp_constraint(
            new_graph_constraint(
                c(NA, NA, 0, 0, 0)
            )
        )
    )

    # danger: all flagged
    expect_snapshot(
        diagnose_hyp_constraint(
            new_graph_constraint(
                c(NA, 0, 0.9, 1.1, -0.1)
            )
        )
    )

    # the check function returns the input invisibly
    a <- new_graph_constraint(c(0, 0.2, 0.1, 0.2, 0.5))

    expect_identical(
        suppressMessages(diagnose_hyp_constraint(a)),
        a
    )
})

test_that("diagnose_hyp_constraint() works with empty vector", {
    # hc is regarded as complete and is flagged
    expect_snapshot(
        diagnose_hyp_constraint(
            new_graph_constraint()
        )
    )
})

test_that("validate_graph_constraint() picks up hypothesis weight issues", {
    tc <- matrix(
        c(
            0, 0, NA, NA, 0.5,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1,
            1, 0, 0, 0, 0
        ),
        nrow = 5,
        byrow = TRUE
    )

    # no error when hc is incomplete, the sum of the weights is less than 1 and
    # there are enough weights to optimise
    expect_no_error(
        validate_graph_constraint(
            graph_constraint(
                c(NA, NA, 0, 0, 0),
                tc
            )
        )
    )

    # it does nothing when the hc is complete and sum of the weights adds to 1
    expect_no_error(
        validate_graph_constraint(
            graph_constraint(
                c(1, 0, 0, 0, 0),
                tc
            )
        )
    )

    # error when hc is complete and the weights add up to more than 1
    expect_snapshot(error = TRUE, {
        validate_graph_constraint(
            new_graph_constraint(
                c(0.8, 0.3, 0, 0, 0),
                tc
            )
        )
    })

    # error - hc value greater than 1
    expect_snapshot(error = TRUE, {
        validate_graph_constraint(
            new_graph_constraint(
                c(NA, NA, 0, 0.2, 1.9),
                tc
            )
        )
    })

    # error picking up on a single NA value in hyp_constraints
    expect_snapshot(error = TRUE, {
        validate_graph_constraint(
            new_graph_constraint(
                c(NA, 0, 0, 0.2, 0),
                tc
            )
        )
    })

    # success: nothing flagged
    expect_snapshot(
        validate_graph_constraint(
            graph_constraint(
                c(NA, NA, 0, 0, 0),
                tc
            ),
            diagnose = TRUE
        )
    )

    # danger: all flagged + error as validate_hyp_constraint() calls
    # assert_hyp_constraint() # nolint
    expect_snapshot(error = TRUE, {
        validate_graph_constraint(
            new_graph_constraint(
                c(NA, 0, 0.9, 1.1, -0.1),
                tc
            ),
            diagnose = TRUE
        )
    })
})

test_that("assert_hyp_constraint_values() with tolerance", {
    # no error - one value greater than 1, but within tolerance
    expect_snapshot(
        assert_hyp_constraint_values(
            new_graph_constraint(
                c(NA, NA, 0, 0, 1 + 10e-13)
            )
        )
    )

    # no error - one value greater than 1, but with tolerance = 0
    expect_snapshot(error = TRUE, {
        assert_hyp_constraint_values(
            new_graph_constraint(
                c(NA, NA, 0, 0, 1 + 10e-13),
                tolerance = 0
            )
        )
    })
})

test_that("diagnose_hyp_constraint_gt1() with tolerance", {
    # nothing flagged: some values are greater than 1 but within tolerance
    expect_snapshot(
        diagnose_hyp_constraint_gt1(
            new_graph_constraint(
                c(NA, NA, 0, 0.2, 1 + 10e-13)
            )
        )
    )

    # some values are greater than 1 + value flagged as not within tolerance
    expect_snapshot(
        diagnose_hyp_constraint_gt1(
            new_graph_constraint(
                c(NA, NA, 0, 0.2, 1 + 10e-13),
                tolerance = 0
            )
        )
    )
})
