test_that(".build_start_matrix produces a valid output - no constraints", {
    gc <- graph_constraint_free(4)
    start_graph <- list(
        list(
            hyp_weight = c(0.1, 0.5, 0.4, 0),
            trans_matrix = gMCPLite::parallelGatekeeping()@m
        ),
        list(
            hyp_weight = c(0.6, 0.2, 0.2, 0.2),
            trans_matrix = gMCPLite::parallelGatekeeping()@m
        ),
        list(
            hyp_weight = c(0.6, 0.2, 0.2, 0.2),
            trans_matrix = gMCPLite::generalSuccessive(
                gamma = 0.4,
                delta = 0.2
            )@m
        )
    )

    s <- .build_start_matrix(gc, start_graph)
    expect_true(is.matrix(s))
    expect_gt(nrow(s), 1)
    expect_identical(ncol(s), length(create_start_params(gc)))
    expect_true(all(is.finite(s)))
    # no duplicate rows
    expect_identical(nrow(unique(s)), nrow(s))
    expect_identical(nrow(s), 5L)
    expect_identical(ncol(s), 11L)
})

test_that(".build_start_matrix produces a valid output - WITH constraints", {
    gc <- graph_constraint(
        hyp_constraint = c(NA, 0.2, 0.2, NA),
        trans_constraint = rbind(
            c(0, 0.8, NA, NA),
            c(0.8, 0, NA, NA),
            c(0.8, NA, 0, NA),
            c(0.8, NA, NA, 0)
        )
    )
    start_graph <- list(
        list(
            hyp_weight = c(0.1, 0.5, 0.4, 0),
            trans_matrix = gMCPLite::parallelGatekeeping()@m
        ),
        list(
            hyp_weight = c(0.6, 0.2, 0.2, 0.2),
            trans_matrix = gMCPLite::parallelGatekeeping()@m
        ),
        list(
            hyp_weight = c(0.6, 0.2, 0.2, 0.2),
            trans_matrix = gMCPLite::generalSuccessive(
                gamma = 0.4,
                delta = 0.2
            )@m
        )
    )

    s <- .build_start_matrix(gc, start_graph)
    expect_true(is.matrix(s))
    expect_gt(nrow(s), 1)
    expect_identical(ncol(s), length(create_start_params(gc)))
    expect_true(all(is.finite(s)))
    # no duplicate rows
    expect_identical(nrow(unique(s)), nrow(s))
    expect_identical(nrow(s), 5L)
    expect_identical(ncol(s), 5L)
})

test_that(".build_start_matrix with default start graph", {
    gc <- graph_constraint_free(5)

    expect_snapshot({
        .build_start_matrix(
            gc,
            start_graph = NULL
        )
    })

    expect_snapshot({
        .build_start_matrix(
            gc,
            start_graph = list(
                list(
                    hyp_weight = NULL,
                    trans_matrix = NULL
                )
            )
        )
    })
})

test_that(".is_default_start_graph", {
    expect_true(
        .is_default_start_graph(NULL)
    )

    expect_true(
        .is_default_start_graph(
            list(
                list(
                    hyp_weight = NULL,
                    trans_matrix = NULL
                )
            )
        )
    )

    expect_false(
        .is_default_start_graph("foo")
    )
})

test_that(".validate_start_graphs with default start graphs", {
    expect_no_error(
        .validate_start_graphs(NULL)
    )

    expect_no_error(
        .validate_start_graphs(
            list(
                list(
                    hyp_weight = NULL,
                    trans_matrix = NULL
                )
            )
        )
    )
})

test_that(".validate_start_graphs with other start graphs", {
    hyp_w <- c(0.1, 0.2, NA, NA, NA)
    trans_m <- matrix(
        rep_len(
            c(0.1, NA),
            length.out = 25
        ),
        nrow = 5,
        ncol = 5
    )
    diag(trans_m) <- 0

    expect_no_error(
        .validate_start_graphs(
            list(
                list(
                    hyp_weight = hyp_w,
                    trans_matrix = trans_m
                )
            ),
            m = 5
        )
    )

    expect_snapshot(error = TRUE, {
        .validate_start_graphs(
            list(
                list(
                    hyp_weight = c(0.1, 0.2, NA, NA),
                    trans_matrix = trans_m
                )
            ),
            m = 5
        )
    })

    expect_snapshot(error = TRUE, {
        .validate_start_graphs(
            list(
                list(
                    hyp_weight = hyp_w,
                    trans_matrix = hyp_w
                )
            ),
            m = 5
        )
    })

    trans_m_4 <- matrix(
        rep_len(
            c(0.1, NA),
            length.out = 16
        ),
        nrow = 4,
        ncol = 4
    )
    diag(trans_m) <- 0

    expect_snapshot(error = TRUE, {
        .validate_start_graphs(
            list(
                list(
                    hyp_weight = hyp_w,
                    trans_matrix = trans_m_4
                )
            ),
            m = 5
        )
    })
})
