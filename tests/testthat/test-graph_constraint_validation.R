test_that("validate_graph_constraint() works", {
    expect_snapshot(
        validate_graph_constraint(
            graph_constraint(
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
            ),
            diagnose = TRUE
        )
    )

    # the hc has 4 elements while the tc is 5 * 5
    # the consistency assertion should error: the tc is inconsistent with the hc
    expect_snapshot(error = TRUE, {
        validate_graph_constraint(
            new_graph_constraint(
                c(NA, NA, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.1, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.9, 0, 0.2, 0,
                        1, 0, NA, 0, 0,
                        1.1, 0, 0, 0, 0
                    ),
                    nrow = 5,
                    byrow = TRUE
                )
            )
        )
    })
})

test_that("assert_hc_tc_consistency() works", {
    expect_no_error(
        assert_hc_tc_consistency(
            graph_constraint(
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

    # the hc has 4 elements while the tc is 5 * 5
    # the consistency assertion should error: the tc is inconsistent with the hc
    expect_snapshot(error = TRUE, {
        assert_hc_tc_consistency(
            new_graph_constraint(
                c(NA, NA, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.1, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.9, 0, 0.2, 0,
                        1, 0, NA, 0, 0,
                        1.1, 0, 0, 0, 0
                    ),
                    nrow = 5,
                    byrow = TRUE
                )
            )
        )
    })
})

test_that("diagnose_hc_tc_consistency() works", {
    # hyp_constraint and trans_constraint have similar dimensions
    # the consistency diagnosis should not flag anything
    expect_snapshot(
        diagnose_hc_tc_consistency(
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

    # hyp_constraint has 4 elements while the trans_constraint is 5 * 5
    # the consistency diagnosis should flag the trans_constraint as being
    # inconsistent with hyp_constraint
    expect_snapshot(
        diagnose_hc_tc_consistency(
            new_graph_constraint(
                c(NA, NA, 0, 0),
                matrix(
                    c(
                        0, 0.8, 0.2, 0, 0,
                        NA, 0, 0, NA, 0,
                        0, 0.9, 0, 0.1, 0,
                        1, 0, 0, 0, 0,
                        1, 0, 0, 0, 0
                    ),
                    nrow = 5,
                    byrow = TRUE
                )
            )
        )
    )

    # hyp_constraint has 4 elements while the trans_constraint is 5 columns * 4
    # rows the consistency diagnosis should flag
    expect_snapshot(
        diagnose_hc_tc_consistency(
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

test_that("assert_graph_constraint_not_full() fails when no optimisable vals", {
    expect_snapshot(error = TRUE, {
        assert_graph_constraint_not_full(
            new_graph_constraint(
                c(1, 0, 0),
                rbind(
                    c(0, 1, 0),
                    c(0.3, 0, 0.7),
                    c(1, 0, 0)
                )
            )
        )
    })
})

test_that("graph_constraint() fails when no optimisable values", {
    expect_snapshot(error = TRUE, {
        graph_constraint(
            c(1, 0, 0),
            rbind(
                c(0, 1, 0),
                c(0.3, 0, 0.7),
                c(1, 0, 0)
            ),
            diagnose = TRUE
        )
    })
})


test_that("diagnose_graph_constraint_not_full() reports no optimisable vals", {
    expect_snapshot(
        diagnose_graph_constraint_not_full(
            new_graph_constraint(
                c(1, 0, 0),
                rbind(
                    c(0, 1, 0),
                    c(0.3, 0, 0.7),
                    c(1, 0, 0)
                )
            )
        )
    )
})
