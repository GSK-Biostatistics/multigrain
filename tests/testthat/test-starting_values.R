test_that("default_start_weights works correctly", {
    # All weights are NA
    hyp_weight1 <- rep(NA, 5)
    expect_identical(
        default_start_weights(hyp_weight1),
        rep(0.2, 5)
    )

    # Some weights defined
    hyp_weight2 <- c(NA, NA, NA, NA, 0.3, NA, NA)
    expect_equal(
        default_start_weights(hyp_weight2),
        c(rep(0.7 / 6, 4), 0.3, rep(0.7 / 6, 2))
    )

    # Edge case - weights are NA
    hyp_weight3 <- c(0.5, 0.3, 0.2)
    expect_equal(default_start_weights(hyp_weight3), c(0.5, 0.3, 0.2))

    ## Test edge case - only one NA
    hyp_weight4 <- c(0.1, 0.3, NA, 0.2)
    expect_equal(default_start_weights(hyp_weight4), c(0.1, 0.3, 0.4, 0.2))
})

test_that("default_start_edges distributes transition weights correctly", {
    # All entries are NA
    trans_matrix1 <- matrix(NA, nrow = 3, ncol = 3, byrow = TRUE)
    diag(trans_matrix1) <- 0
    expected_G1 <- matrix(1 / 2, nrow = 3, ncol = 3, byrow = TRUE)
    diag(expected_G1) <- 0
    expect_identical(default_start_edges(trans_matrix1), expected_G1)

    # Some entries are defined
    trans_matrix2 <- matrix(
        c(
            0, NA, 0.8,
            NA, 0, NA,
            NA, 0, 0
        ),
        nrow = 3,
        byrow = TRUE
    )
    expected_G2 <- matrix(
        c(
            0, 0.2, 0.8,
            0.5, 0, 0.5,
            1.0, 0, 0
        ),
        nrow = 3,
        byrow = TRUE
    )

    # nolint start: expect_identical_linter
    expect_equal(
        default_start_edges(trans_matrix2),
        expected_G2
    )
    # nolint end

    # No entries are NA
    trans_matrix3 <- matrix(c(0.5, 0.5, 0, 0.5, 0.5, 0, 0.5, 0.5, 0), nrow = 3)
    expect_identical(default_start_edges(trans_matrix3), trans_matrix3)
})


test_that("create_start_params integrates inputs correctly", {
    # Test no constraints for 8 hypotheses
    graph_constraint8 <- graph_constraint_free(8)
    start8 <- create_start_params(graph_constraint8)
    expect_identical(
        as.vector(start8),
        c(rep(1 / 8, 7), rep(1 / 7, 6 * 8))
    )

    expect_identical(
        attr(start8, "w0_opt"),
        rep(1 / 8, 7)
    )

    expect_identical(
        unlist(attr(start8, "G0_opt")),
        rep(1 / 7, 6 * 8)
    )

    # Test creating default parameters with user defined constraints
    hyp_constraint1 <- c(NA, 0.3, NA)
    trans_constraint1 <- matrix(
        c(
            0, NA, NA,
            NA, 0, NA,
            NA, NA, 0
        ),
        nrow = 3,
        byrow = TRUE
    )

    graph_constraint1 <- new_graph_constraint(
        hyp_constraint = hyp_constraint1,
        trans_constraint = trans_constraint1
    )

    result1 <- create_start_params(graph_constraint1)
    expect_equal(attr(result1, "w0_opt"), 0.35)
    expect_equal(unlist(attr(result1, "G0_opt")), c(0.5, 0.5, 0.5))
    expect_equal(as.vector(result1), c(0.35, 0.5, 0.5, 0.5))

    # Test creating default parameters 2
    hyp_constraint2 <- c(NA, 0.3, NA, NA)
    trans_constraint2 <- matrix(
        c(
            0, NA, 0.1, NA,
            NA, 0, NA, NA,
            NA, NA, 0, 0,
            NA, 0.6, 0.2, 0
        ),
        nrow = 4,
        byrow = TRUE
    )

    graph_constraint2 <- new_graph_constraint(
        hyp_constraint = hyp_constraint2,
        trans_constraint = trans_constraint2
    )

    result2 <- create_start_params(graph_constraint2)

    expect_identical(
        attr(result2, "w0_opt"),
        c(0.7 / 3, 0.7 / 3)
    )

    expect_identical(
        unlist(attr(result2, "G0_opt")),
        c(0.9 / 2, 1 / 3, 1 / 3, 1 / 2)
    )

    expect_identical(
        as.vector(result2),
        c(
            attr(result2, "w0_opt"),
            unlist(
                attr(
                    result2,
                    "G0_opt"
                )
            )
        )
    )

    # Test with custom weights and transition matrix 1 (graph_constraint 1)
    custom_w0 <- c(0.7, 0.3, 0)
    custom_G0 <- matrix(
        c(
            0, 0.7, 0.3,
            0, 0, 1,
            0.7, 0.3, 0
        ),
        nrow = 3,
        byrow = TRUE
    )

    graph_constraint3 <- new_graph_constraint(
        hyp_constraint = hyp_constraint1,
        trans_constraint = trans_constraint1
    )

    result3 <- create_start_params(
        graph_constraint3,
        w0 = custom_w0,
        G0 = custom_G0
    )
    expect_equal(attr(result3, "w0_opt"), 0.7)
    expect_equal(unlist(attr(result3, "G0_opt")), c(0.7, 0, 0.7))

    # Check starting values yielded when there are no parameters to optimise

    # Check error return when elements in w0 do not match hyp_constraint
    wrong_w0 <- c(0.4, 0.4, 0.2)
    expect_error(
        create_start_params(graph_constraint1, w0 = wrong_w0, G0 = custom_G0),
        "Fixed weights in weight constraint vector must match elements in w0"
    )

    # Check error return when elements in G0 do not match trans_matrix
    hyp_constraint2 <- c(NA, 0.3, NA, NA)
    trans_constraint2 <- matrix(
        c(
            0, NA, 0.1, NA,
            NA, 0, NA, NA,
            NA, NA, 0, 0,
            NA, 0.6, 0.2, 0
        ),
        nrow = 4,
        byrow = TRUE
    )

    graph_constraint <- new_graph_constraint(
        hyp_constraint = hyp_constraint2,
        trans_constraint = trans_constraint2
    )

    custom_G2 <- matrix(
        c(
            0, 0.2, 0.1, 0.7,
            0.1, 0, 0.5, 0.4,
            1, 0, 0, 0,
            0.4, 0.4, 0.2, 0
        ),
        nrow = 4,
        byrow = TRUE
    )

    expect_error(
        create_start_params(graph_constraint, G0 = custom_G2),
        "Fixed elements in transition constraint matrix must match elements in G0" # nolint
    )

    # Test error handling
    expect_error(
        create_start_params(
            graph_constraint1,
            w0 = c(0.3, 0.3, 0.6)
        ),
        "Weights must sum to 1."
    )
    expect_error(
        create_start_params(
            graph_constraint1,
            G0 = matrix(
                c(
                    0, 0.5, 0.6,
                    1, 0, 0,
                    1, 0, 0
                ),
                nrow = 3,
                byrow = TRUE
            )
        ),
        "Rows of transition matrix must sum to 1."
    )
})

test_that("create_start_params with fully fixed transition rows; free_g == 0", {
    # Row 1 and row 3 are fully fixed (0 free params); row 2 has 2 free params
    hyp_constraint <- rep(NA_real_, 3)
    trans_constraint <- matrix(
        c(
            0, 0.5, 0.5,
            NA, 0, NA,
            0.5, 0.5, 0
        ),
        nrow = 3,
        byrow = TRUE
    )
    gc <- new_graph_constraint(
        hyp_constraint = hyp_constraint,
        trans_constraint = trans_constraint
    )

    # Should not error — previously crashed with mixed 0/negative indices
    result <- create_start_params(gc)

    # Only row 2 contributes free params (2 free, encoded as 2-1 = 1 param)
    expect_length(attr(result, "G0_opt"), 1)

    # Weight params: 3 free, encoded as 3-1 = 2
    expect_length(attr(result, "w0_opt"), 2)
})

test_that("roundtrip create_start_params->recover_full_weights/trans_matrix", {
    constr <- graph_constraint(
        hyp_constraint = rep(NA_real_, 3),
        trans_constraint = matrix(
            c(
                0, 1, 0,
                NA, 0, NA,
                NA, NA, 0
            ),
            nrow = 3,
            byrow = TRUE
        )
    )

    x <- create_start_params(constr)

    free_w <- sum(is.na(constr$hyp_constraint))
    w <- recover_full_weights(x[1:free_w - 1], constr$hyp_constraint)
    G <- recover_full_trans_matrix(x[free_w:length(x)], constr$trans_constraint)

    expect_identical(sum(w), 1)
    expect_identical(rowSums(G), rep(1, 3))
    expect_true(all(G >= 0))
})
