test_that("new_graph_constraint() can create an empty object", {
    empty_gc <- new_graph_constraint()

    expect_s3_class(empty_gc, "multigrain_graph_constraint")

    expect_named(
        empty_gc,
        c("hyp_constraint", "trans_constraint")
    )

    expect_true(
        rlang::is_empty(empty_gc$hyp_constraint)
    )
    expect_true(
        rlang::is_empty(empty_gc$trans_constraint)
    )
})

test_that("new_graph_constraint() works", {
    # can create a named graph_constraint
    expect_snapshot(
        new_graph_constraint(
            hyp_constraint = c(NA, NA, 0, 0),
            trans_constraint = matrix(
                c(
                    NA, NA, 0, 0,
                    NA, NA, 0, 0,
                    NA, NA, 0, 0,
                    NA, NA, 0, 0
                ),
                nrow = 4,
                byrow = TRUE
            ),
            names = letters[1:4]
        )
    )

    expect_s3_class(
        new_graph_constraint(
            hyp_constraint = c(NA, NA, 0, 0),
            trans_constraint = matrix(
                c(
                    NA, NA, 0, 0,
                    NA, NA, 0, 0,
                    NA, NA, 0, 0,
                    NA, NA, 0, 0
                ),
                nrow = 4,
                byrow = TRUE
            ),
            names = letters[1:4]
        ),
        "multigrain_graph_constraint"
    )

    gc <- new_graph_constraint(
        hyp_constraint = c(NA, NA, 0, 0),
        trans_constraint = matrix(
            c(
                NA, NA, 0, 0,
                NA, NA, 0, 0,
                NA, NA, 0, 0,
                NA, NA, 0, 0
            ),
            nrow = 4,
            byrow = TRUE
        ),
        names = letters[1:4]
    )

    expect_named(
        gc$hyp_constraint,
        letters[1:4]
    )

    expect_identical(
        dimnames(gc$trans_constraint),
        list(
            letters[1:4],
            letters[1:4]
        )
    )

    # can create an unnamed `multigrain_graph_constraint`
    expect_snapshot(
        new_graph_constraint(
            hyp_constraint = c(NA, NA, 0, 0),
            trans_constraint = matrix(
                c(
                    NA, NA, 0, 0,
                    NA, NA, 0, 0,
                    NA, NA, 0, 0,
                    NA, NA, 0, 0
                ),
                nrow = 4,
                byrow = TRUE
            )
        )
    )

    expect_s3_class(
        new_graph_constraint(
            hyp_constraint = c(NA, NA, 0, 0),
            trans_constraint = matrix(
                c(
                    NA, NA, 0, 0,
                    NA, NA, 0, 0,
                    NA, NA, 0, 0,
                    NA, NA, 0, 0
                ),
                nrow = 4,
                byrow = TRUE
            )
        ),
        "multigrain_graph_constraint"
    )

    gc_no_names <- new_graph_constraint(
        hyp_constraint = c(NA, NA, 0, 0),
        trans_constraint = matrix(
            c(
                NA, NA, 0, 0,
                NA, NA, 0, 0,
                NA, NA, 0, 0,
                NA, NA, 0, 0
            ),
            nrow = 4,
            byrow = TRUE
        )
    )

    expect_named(
        gc_no_names$hyp_constraint,
        NULL
    )

    expect_null(
        dimnames(gc_no_names$trans_constraint)
    )
})

test_that("new_graph_constraint() complains with undesired inputs", {
    expect_snapshot(error = TRUE, {
        new_graph_constraint(
            hyp_constraint = "a"
        )
    })

    expect_snapshot(error = TRUE, {
        new_graph_constraint(
            hyp_constraint = c(NA, NA, 0, 0, 0),
            trans_constraint = "b"
        )
    })

    expect_snapshot(error = TRUE, {
        new_graph_constraint(
            hyp_constraint = c(NA, NA, 0, 0, 0),
            trans_constraint = c(NA, NA, 0, 0, 0)
        )
    })

    expect_snapshot(error = TRUE, {
        new_graph_constraint(
            hyp_constraint = c(NA, NA, 0, 0, 0),
            trans_constraint = matrix(
                c(
                    NA, NA,
                    0, 0
                ),
                nrow = 2,
                byrow = TRUE
            ),
            names = c(NA, NA, 0, 0, 0)
        )
    })

    expect_snapshot(error = TRUE, {
        new_graph_constraint(
            hyp_constraint = c(NA, NA, 0, 0),
            trans_constraint = matrix(
                c(
                    0, NA, NA, 0,
                    NA, 0, NA, 0,
                    NA, NA, 0, 0,
                    NA, NA, 0, 0
                ),
                nrow = 4,
                byrow = TRUE
            ),
            tolerance = "a"
        )
    })
})

test_that("graph_constraint() when one of the inputs is NULL", {
    # null `trans_constraint`
    expect_snapshot(
        graph_constraint(
            hyp_constraint = c(NA, 0.4, NA)
        )
    )

    # null `hyp_constraint`
    expect_snapshot(
        graph_constraint(
            trans_constraint = matrix(
                c(
                    0, NA, NA, 0,
                    NA, 0, NA, 0,
                    NA, NA, 0, 0,
                    NA, NA, 0, 0
                ),
                nrow = 4,
                byrow = TRUE
            )
        )
    )
})

test_that("graph_constraint() errors with both inputs NULL", {
    # i.e. when both `hyp_constraint` & `trans_constraint` are `NULL`
    expect_snapshot(error = TRUE, {
        graph_constraint()
    })
})

test_that("graph_constraint() errors trans_constraint not numeric", {
    expect_snapshot(error = TRUE, {
        graph_constraint(
            trans_constraint = matrix(
                letters[1:16],
                nrow = 4,
                byrow = TRUE
            )
        )
    })
})

test_that("graph_constraint complains when anything is passed via `...`", {
    expect_snapshot(error = TRUE, {
        graph_constraint(
            hyp_constraint = c(NA, NA, 0, 0),
            trans_constraint = matrix(
                c(
                    0, NA, NA, 0,
                    NA, 0, NA, 0,
                    NA, NA, 0, 0,
                    NA, NA, 0, 0
                ),
                nrow = 4,
                byrow = TRUE
            ),
            c("a", "b", "c", "d")
        )
    })
})

test_that("graph_constraint() errors when tolerance not positive numeric", {
    expect_snapshot(error = TRUE, {
        graph_constraint(
            hyp_constraint = c(NA, NA, 0, 0),
            trans_constraint = matrix(
                c(
                    0, NA, NA, 0,
                    NA, 0, NA, 0,
                    NA, NA, 0, 0,
                    NA, NA, 0, 0
                ),
                nrow = 4,
                byrow = TRUE
            ),
            tolerance = "a"
        )
    })

    expect_snapshot(error = TRUE, {
        graph_constraint(
            hyp_constraint = c(NA, NA, 0, 0),
            trans_constraint = matrix(
                c(
                    0, NA, NA, 0,
                    NA, 0, NA, 0,
                    NA, NA, 0, 0,
                    NA, NA, 0, 0
                ),
                nrow = 4,
                byrow = TRUE
            ),
            tolerance = -1e-6
        )
    })
})

test_that("graph_constraint with names length mismatch", {
    expect_error(
        graph_constraint(
            c(NA, NA, 0, 0),
            names = letters[1:5]
        ),
        "`names` must have 4 elements. It has 5."
    )
})


test_that("graph_constraint_free() works", {
    expect_s3_class(
        graph_constraint_free(4),
        "multigrain_graph_constraint"
    )

    expect_snapshot(
        graph_constraint_free(4)
    )
})

test_that("graph_constraint_free() works with 2 hypotheses", {
    # 2-hypotheses is an edge case -> test it
    expect_snapshot(
        graph_constraint_free(2)
    )
})

test_that("graph_constraint_free() complains", {
    # with non-numeric input
    expect_error(
        graph_constraint_free("a"),
        '`num_hyp` must be a whole number, not the string "a".',
        fixed = TRUE
    )

    # or when the input is a double
    expect_error(
        graph_constraint_free(2.2),
        "`num_hyp` must be a whole number, not the number 2.2.",
        fixed = TRUE
    )

    # or when the input is not scalar
    expect_error(
        graph_constraint_free(c(2, 3)),
        "`num_hyp` must be a whole number, not a double vector.",
        fixed = TRUE
    )

    # input is less than 2
    expect_error(
        graph_constraint_free(1),
        "`num_hyp` must be a whole number larger than or equal to 2",
        fixed = TRUE
    )
})

test_that("graph_constraint_free complains with names length mismatch", {
    expect_error(
        graph_constraint_free(
            5,
            names = letters[1:4]
        ),
        "`names` must have 5 elements. It has 4."
    )
})

test_that("is_graph_constraint()", {
    expect_true(
        is_graph_constraint(
            graph_constraint_free(
                4
            )
        )
    )

    expect_true(
        is_graph_constraint(
            graph_constraint(
                c(0, 0, 0, NA, NA)
            )
        )
    )

    expect_false(
        is_graph_constraint(
            "foo"
        )
    )

    expect_true(
        is_graph_constraint(
            new_graph_constraint()
        )
    )
})

test_that("check_graph_constraint", {
    expect_no_error(
        check_graph_constraint(
            graph_constraint(
                hyp_constraint = c(NA, 0.4, NA)
            )
        )
    )

    expect_no_error(
        check_graph_constraint(
            NULL,
            allow_null = TRUE
        )
    )

    expect_error(
        check_graph_constraint(
            NULL
        ),
        "`NULL` must be a multigrain graph constraint object, not `NULL`."
    )

    expect_error(
        check_graph_constraint("foo"),
        '`"foo"` must be a multigrain graph constraint object'
    )

    expect_error(
        check_graph_constraint(2),
        "`2` must be a multigrain graph constraint object, not the number 2."
    )
})

test_that("hyp_constraint_free() works", {
    expect_identical(
        hyp_constraint_free(2),
        c(NA_real_, NA_real_)
    )

    expect_identical(
        hyp_constraint_free(3),
        c(NA_real_, NA_real_, NA_real_)
    )
})

test_that("trans_constraint_free() works", {
    expect_identical(
        trans_constraint_free(2),
        matrix(
            c(
                0, 1,
                1, 0
            ),
            byrow = TRUE,
            ncol = 2
        )
    )

    expect_identical(
        trans_constraint_free(3),
        matrix(
            c(
                0, NA, NA,
                NA, 0, NA,
                NA, NA, 0
            ),
            byrow = TRUE,
            ncol = 3
        )
    )
})

test_that("graph_constraint: users can update hyp_constraint", {
    gc <- graph_constraint(
        hyp_constraint = c(NA, NA, 0, 0),
        trans_constraint = matrix(
            c(
                0, NA, NA, 0,
                NA, 0, NA, 0,
                NA, NA, 0, 0,
                NA, NA, 0, 0
            ),
            nrow = 4,
            byrow = TRUE
        )
    )

    # disable the linter as assignment is what we actually want to test
    expect_snapshot(error = TRUE, {
        gc$hyp_constraint <- "A"
    })

    expect_snapshot(error = TRUE, {
        gc[["hyp_constraint"]] <- c(1, NA, 0, 0)
    })

    # tolerance works
    expect_snapshot(error = TRUE, {
        gc["hyp_constraint", tolerance = 0] <- c(1 + 10e-12, 0, 0, 0)
    })

    expect_snapshot(
        gc["hyp_constraint"] <- c(1 + 10e-13, 0, 0, 0)
    )

    # diagnose works
    expect_snapshot(error = TRUE, {
        # fmt: skip
        gc["hyp_constraint", tolerance = 0, diagnose = TRUE] <-
            c(1 + 10e-12, 0, 0, 0)
    })

    expect_snapshot(
        # fmt: skip
        gc["hyp_constraint", tolerance = 10e-12, diagnose = TRUE] <-
            c(1 + 10e-13, 0, 0, 0)
    )

    expect_equal(
        gc$hyp_constraint,
        c(1 + 10e-13, 0, 0, 0),
        # don't check names
        ignore_attr = TRUE
    )

    expect_s3_class(gc, "multigrain_graph_constraint")
})

test_that("graph_constraint: users can update trans_constraint", {
    gc <- graph_constraint(
        hyp_constraint = c(NA, NA, 0, 0)
    )

    new_tc <- matrix(
        c(
            0, NA, NA, 0,
            0.7, 0, 0.3, 0,
            NA, NA, 0, 0.1,
            NA, NA, 0.2, 0
        ),
        nrow = 4,
        byrow = TRUE
    )

    expect_snapshot(error = TRUE, {
        gc$trans_constraint <- "A"
    })

    expect_snapshot(
        gc[["trans_constraint"]] <- new_tc
    )

    expect_equal(
        gc$trans_constraint,
        new_tc,
        # don't check names
        ignore_attr = TRUE
    )

    tolerance_tc <- matrix(
        c(
            0, 0, 0, 1 + 10e-13,
            0.7, 0, 0.3, 0,
            NA, NA, 0, 0.1,
            NA, NA, 0.2, 0
        ),
        nrow = 4,
        byrow = TRUE
    )

    # tolerance works
    expect_snapshot(error = TRUE, {
        gc[["trans_constraint", tolerance = 0]] <- tolerance_tc
    })

    expect_snapshot(
        gc[["trans_constraint", tolerance = 10e-12]] <- tolerance_tc
    )

    # diagnosis works
    expect_snapshot(error = TRUE, {
        gc[["trans_constraint", tolerance = 0, diagnose = TRUE]] <- tolerance_tc
    })

    expect_snapshot(
        gc[[
            "trans_constraint",
            tolerance = 10e-12,
            diagnose = TRUE
        ]] <- tolerance_tc
    )

    expect_equal(
        gc$trans_constraint,
        tolerance_tc,
        # don't check names
        ignore_attr = TRUE
    )

    expect_s3_class(gc, "multigrain_graph_constraint")
})

test_that("graph_constraint: update incoming names are preferred", {
    gc <- graph_constraint(
        hyp_constraint = c(NA, NA, 0, 0)
    )

    # a named matrix changes the names
    new_tc <- matrix(
        c(
            0, NA, NA, 0,
            0.7, 0, 0.3, 0,
            NA, NA, 0, 0.1,
            NA, NA, 0.2, 0
        ),
        nrow = 4,
        byrow = TRUE
    )

    dimnames(new_tc) <- list(
        letters[10:13],
        letters[10:13]
    )

    expect_snapshot(
        gc[["trans_constraint"]] <- new_tc
    )

    expect_identical(
        gc$trans_constraint,
        new_tc
    )

    expect_named(
        gc$hyp_constraint,
        letters[10:13]
    )

    expect_s3_class(gc, "multigrain_graph_constraint")

    # a named hyp_constraint changes the names
    gc <- graph_constraint(
        hyp_constraint = c(NA, NA, 0, 0)
    )

    expect_no_error(
        gc$hyp_constraint <- c(A1 = 0, A2 = 0, A3 = 0.3, A4 = 0.7)
    )

    expect_identical(
        dimnames(gc$trans_constraint),
        list(
            c("A1", "A2", "A3", "A4"),
            c("A1", "A2", "A3", "A4")
        )
    )
})

test_that("graph_constraint print and summary methods", {
    gc <- graph_constraint(
        hyp_constraint = c(NA, NA, 0, 0),
        trans_constraint = matrix(
            c(
                0, NA, NA, 0,
                NA, 0, NA, 0,
                NA, NA, 0, 0,
                NA, NA, 0, 0
            ),
            nrow = 4,
            byrow = TRUE
        )
    )

    expect_snapshot(print(gc))
    expect_snapshot(summary(gc))

    expect_null(print.multigrain_graph_constraint(NULL))
    expect_null(summary.multigrain_graph_constraint(NULL))
})

test_that("set methods inherit the original tolerance if unspecified", {
    expect_snapshot(error = TRUE, {
        graph_constraint(
            hyp_constraint = c(NA, NA, 0, 0),
            tolerance = "a"
        )
    })

    gc <- graph_constraint(
        hyp_constraint = c(NA, NA, 0, 0),
        tolerance = 1e-2
    )

    expect_identical(
        attr(gc, "tolerance"),
        1e-2
    )

    # no error since 1e-3 < tolerance (1e-2)
    expect_no_error(
        gc[["hyp_constraint"]] <- c(1 + 1e-3, 0, 0, 0)
    )

    # tolerance is unchanged (persists from the first gc definition)
    expect_identical(
        attr(gc, "tolerance"),
        1e-2
    )

    # this errors since we have modified the tolerance
    expect_snapshot(error = TRUE, {
        gc[["hyp_constraint", tolerance = 1e-5]] <- c(1 + 1e-3, 0, 0, 0)
    })

    # this works since the value is within tolerance, but we have also updated
    # the tolerance of the graph_constraint object
    expect_no_error(
        gc[["hyp_constraint", tolerance = 1e-5]] <- c(1 + 1e-6, 0, 0, 0)
    )

    expect_identical(
        attr(gc, "tolerance"),
        1e-5
    )

    gc <- graph_constraint(
        hyp_constraint = c(NA, NA, 0, 0),
        trans_constraint = matrix(
            c(
                0, NA, NA, 0,
                0.7, 0, 0.3 + 1e-13, 0,
                NA, NA, 0, 0.1,
                NA, NA, 0.2, 0
            ),
            nrow = 4,
            byrow = TRUE
        ),
        tolerance = 1e-13
    )

    # if we modify the tolerance during the update of hyp_constraint, the
    # validation fails for the transition matrix
    expect_snapshot(error = TRUE, {
        gc[["hyp_constraint", tolerance = 0]] <- c(1, 0, 0, 0)
    })
})


##---- Unit tests for closest_graph_to_constraints() ---##
test_that("weights: pin fixed entries and preserve sum=1", {
    set.seed(1)
    m <- 4
    hc <- c(NA, 0.2, NA, NA) # fix w2 = 0.2
    tc <- matrix(NA_real_, m, m)
    diag(tc) <- 0

    gc <- graph_constraint(
        hyp_constraint = hc,
        trans_constraint = tc
    )

    w_in <- runif(m)
    w_in <- w_in / sum(w_in)
    G_in <- matrix(runif(m * m), m, m)
    diag(G_in) <- 0
    G_in <- sweep(G_in, 1, rowSums(G_in), "/")

    out <- closest_graph_to_constraints(gc, w = w_in, G = G_in)

    expect_equal(out$hyp_weight[2], 0.2, tolerance = 1e-12)
    expect_equal(sum(out$hyp_weight), 1, tolerance = 1e-12)
    expect_true(all(out$hyp_weight >= 0))
    expect_true(all(abs(rowSums(out$trans_matrix) - 1) < 1e-12))
    expect_identical(diag(out$trans_matrix), rep(0, m))
})

test_that("weights: split evenly when all free proposals are zero", {
    m <- 4
    hc <- c(0.3, NA, NA, NA) # fixed mass 0.3; surplus 0.7
    tc <- matrix(NA_real_, m, m)
    diag(tc) <- 0

    gc <- graph_constraint(hyp_constraint = hc, trans_constraint = tc)

    w_in <- c(0, 0, 0, 0) # free part has zero mass
    G_in <- diag(0, m) # any valid G (row sums will be handled)

    out <- closest_graph_to_constraints(gc, w = w_in, G = G_in)

    expect_equal(out$hyp_weight[1], 0.3, tolerance = 1e-12)
    expect_true(all(abs(out$hyp_weight[2:4] - (0.7 / 3)) < 1e-12))
    expect_equal(sum(out$hyp_weight), 1, tolerance = 1e-12)
})

test_that("weights: returns pinned values when no free entries and sum==1", {
    m <- 3
    hc <- c(0.2, 0.3, 0.5) # no free entries; already sums to 1
    tc <- matrix(NA_real_, m, m)
    diag(tc) <- 0

    gc <- graph_constraint(hyp_constraint = hc, trans_constraint = tc)

    w_in <- c(0.9, 0.05, 0.05) # ignored (no free)
    G_in <- diag(0, m)

    out <- closest_graph_to_constraints(gc, w = w_in, G = G_in)
    expect_equal(out$hyp_weight, hc, tolerance = 1e-12)
})


test_that("weights: tolerates tiny rounding of fixed entries", {
    m <- 3
    eps <- 1e-14
    hc <- c(0.4 + eps, NA, NA)
    tc <- matrix(NA_real_, m, m)
    diag(tc) <- 0

    gc <- graph_constraint(hyp_constraint = hc, trans_constraint = tc)

    w_in <- runif(m)
    w_in <- w_in / sum(w_in)
    G_in <- diag(0, m)

    out <- closest_graph_to_constraints(gc, w = w_in, G = G_in)
    expect_lt(abs(out$hyp_weight[1] - (0.4 + eps)), 1e-12)
    expect_equal(sum(out$hyp_weight), 1, tolerance = 1e-12)
})

test_that("G: pin fixed off-diagonal and preserve row sums = 1", {
    set.seed(2)
    m <- 4
    hc <- rep(NA_real_, m)
    tc <- matrix(NA_real_, m, m)
    diag(tc) <- 0
    tc[2, 3] <- 0.6 # fix G[2,3] = 0.6 in row 2

    gc <- graph_constraint(hyp_constraint = hc, trans_constraint = tc)

    w_in <- runif(m)
    w_in <- w_in / sum(w_in)
    G_in <- matrix(runif(m * m), m, m)
    diag(G_in) <- 0
    G_in <- sweep(G_in, 1, rowSums(G_in), "/")

    out <- closest_graph_to_constraints(gc, w = w_in, G = G_in)

    expect_equal(out$trans_matrix[2, 3], 0.6, tolerance = 1e-12)
    expect_true(all(abs(rowSums(out$trans_matrix) - 1) < 1e-12))
    expect_true(all(out$trans_matrix >= 0))
    # w unchanged (no constraints on w)
    expect_equal(out$hyp_weight, w_in, tolerance = 0, ignore_attr = TRUE)
})


test_that("G: handles rows with zero free mass by even split", {
    m <- 4
    hc <- rep(NA_real_, m)
    tc <- matrix(NA_real_, m, m)
    diag(tc) <- 0
    tc[1, 2] <- 0.9 # leaves 0.1 for the rest of row 1

    gc <- graph_constraint(hyp_constraint = hc, trans_constraint = tc)

    w_in <- c(0.4, 0.3, 0.2, 0.1)
    G_in <- matrix(0, m, m)
    diag(G_in) <- 0 # free proposals all zero in row 1

    out <- closest_graph_to_constraints(gc, w = w_in, G = G_in)

    expect_equal(out$trans_matrix[1, 2], 0.9, tolerance = 1e-12)
    # remaining two off-diagonals (3,4) split the 0.1 equally
    expect_equal(sum(out$trans_matrix[1, c(3, 4)]), 0.1, tolerance = 1e-12)
    expect_true(all(abs(rowSums(out$trans_matrix) - 1) < 1e-12))
})


test_that("projection is idempotent when inputs already satisfy constraints", {
    m <- 4
    hc <- c(NA, 0.5, NA, NA)
    tc <- matrix(NA_real_, m, m)
    diag(tc) <- 0
    tc[3, 1] <- 0.2

    gc <- graph_constraint(hyp_constraint = hc, trans_constraint = tc)

    # Construct w/g that already satisfy:
    w_in <- c(0.3, 0.5, 0.1, 0.1)
    G_in <- matrix(0, m, m)
    diag(G_in) <- 0
    G_in[3, 1] <- 0.2
    G_in[3, 2] <- 0.5
    G_in[3, 4] <- 0.3
    # Fill remaining rows with valid distributions
    set.seed(5)
    for (i in setdiff(seq_len(m), 3)) {
        tmp <- runif(m)
        tmp[i] <- 0
        tmp <- tmp / sum(tmp)
        G_in[i, ] <- tmp
        G_in[i, i] <- 0
    }

    out <- closest_graph_to_constraints(gc, w = w_in, G = G_in)
    expect_equal(out$hyp_weight, w_in, tolerance = 1e-12)
    expect_equal(out$trans_matrix, G_in, tolerance = 1e-12)
    expect_true(all(abs(rowSums(out$trans_matrix) - 1) < 1e-12))
})

test_that("build_hyp_names", {
    expect_identical(
        build_hyp_names(4),
        c("H1", "H2", "H3", "H4")
    )
})

test_that("get graph constraint names", {
    gc <- graph_constraint_free(6)

    expect_identical(
        graph_constraint_get_names(gc),
        c("H1", "H2", "H3", "H4", "H5", "H6")
    )
})

test_that("get graph constraint m", {
    gc <- graph_constraint_free(5)

    expect_identical(
        graph_constraint_get_m(gc),
        5L
    )
})

test_that("get graph constraint tolerance", {
    gc <- graph_constraint_free(5)

    expect_identical(
        graph_constraint_get_tolerance(gc),
        sqrt(.Machine$double.eps)
    )

    gc_custom <- graph_constraint(
        hyp_constraint = c(NA, NA, NA),
        tolerance = 1e-6
    )

    expect_identical(
        graph_constraint_get_tolerance(gc_custom),
        1e-6
    )

    expect_error(
        graph_constraint_get_tolerance(list()),
        class = "rlang_error"
    )
})

test_that("graph constraint attributes", {
    gc <- graph_constraint_free(5)

    expect_identical(
        attr(gc, "m"),
        5L
    )

    expect_identical(
        attr(gc, "tolerance"),
        sqrt(.Machine$double.eps)
    )
})
