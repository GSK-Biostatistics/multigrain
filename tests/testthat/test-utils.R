# to ensure consistent cli output (prevents local vs GHA whitespace differences)
cli::start_app()
on.exit(cli::stop_app(), add = TRUE)

### Tests for is_graph_valid()
test_that("returns TRUE for a valid 3x3 graph", {
    hyp_weight <- c(0.4, 0.3, 0.3)
    trans_matrix <- matrix(
        c(
            0,   0.5, 0.5,
            0.3, 0,   0.7,
            0.6, 0.4, 0
        ),
        nrow = 3,
        byrow = TRUE
    )

    expect_true(is_graph_valid(hyp_weight, trans_matrix))
})

test_that("trans_matrix must be a matrix", {
    hyp_weight <- c(10.5, 0.5)
    not_a_matrix <- c(1, 0, 0, 1)

    expect_warning(
        is_graph_valid(hyp_weight, not_a_matrix),
        regexp = "not a matrix"
    )
    expect_false(is_graph_valid(hyp_weight, not_a_matrix)) |>
        suppressWarnings()
})

test_that("trans_matrix must be square", {
    hyp_weight <- c(0.5, 0.5)
    trans_matrix <- matrix(0, nrow = 2, ncol = 3)

    expect_warning(
        is_graph_valid(hyp_weight, trans_matrix),
        regexp = "not a square"
    )
    expect_false(is_graph_valid(hyp_weight, trans_matrix)) |>
        suppressWarnings()
})

test_that("hyp_weight length must match matrix dimension", {
    hyp_weight <- c(0.2, 0.3, 0.5) # length 3
    trans_matrix <- matrix(
        c(0, 1, 1, 0),
        nrow = 2,
        byrow = TRUE
    )

    expect_warning(
        is_graph_valid(hyp_weight, trans_matrix),
        regexp = "Length of `hyp_weight`"
    )
    expect_false(is_graph_valid(hyp_weight, trans_matrix)) |>
        suppressWarnings()
})

test_that("diagonal entries of trans_matrix must be zero", {
    hyp_weight <- c(0.6, 0.4)
    trans_matrix <- matrix(
        c(
            0.1, 0.9,
            0.9, 0
        ),
        nrow = 2,
        byrow = TRUE
    )

    expect_warning(
        is_graph_valid(hyp_weight, trans_matrix),
        regexp = "Diagonals of `trans_matrix` are not all zero"
    )
    expect_false(is_graph_valid(hyp_weight, trans_matrix)) |>
        suppressWarnings()
})

test_that("hyp_weight elements must lie in [0, 1]", {
    hyp_weight <- c(-0.1, 1.1) # outside the [0, 1] interval
    trans_matrix <- matrix(
        c(
            0, 1,
            1, 0
        ),
        nrow = 2,
        byrow = TRUE
    )

    expect_warning(
        is_graph_valid(hyp_weight, trans_matrix),
        regexp = "hyp_weight.*outside the interval \\[0, 1\\]"
    )
    expect_false(is_graph_valid(hyp_weight, trans_matrix)) |>
        suppressWarnings()
})


test_that("trans_matrix elements must lie in [0, 1]", {
    hyp_weight <- c(0.5, 0.5)
    trans_matrix <- matrix(
        c(
            0, 1.1,  # 1.1 invalid
            0.9, 0
        ),
        nrow = 2,
        byrow = TRUE
    )

    expect_warning(
        is_graph_valid(hyp_weight, trans_matrix),
        regexp = "trans_matrix.*outside the interval \\[0, 1\\]"
    )
    expect_false(is_graph_valid(hyp_weight, trans_matrix)) |>
        suppressWarnings()
})


test_that("hyp_weight must sum to 1 within tolerance", {
    hyp_weight_bad <- c(0.5, 0.6) # sums to 1.1
    trans_matrix <- matrix(
        c(
            0, 1,
            1, 0
        ),
        nrow = 2,
        byrow = TRUE
    )

    expect_warning(
        is_graph_valid(hyp_weight_bad, trans_matrix),
        regexp = "`hyp_weight` does not sum to 1"
    )
    expect_false(is_graph_valid(hyp_weight_bad, trans_matrix)) |>
        suppressWarnings()
})

test_that("row sums must be 1 when sum_to_one_constraint = TRUE", {
    hyp_weight <- c(0.5, 0.5)
    trans_matrix <- matrix(
        c(
            0, 1,
            0.8, 0 # row 2 sums to 0.8
        ),
        nrow = 2,
        byrow = TRUE
    )

    expect_warning(
        is_graph_valid(hyp_weight, trans_matrix, sum_to_one_constraint = TRUE),
        regexp = "rows of `trans_matrix` do not sum to 1"
    )
    expect_false(
        is_graph_valid(
            hyp_weight,
            trans_matrix,
            sum_to_one_constraint = TRUE
        )
    ) |>
        suppressWarnings()
})


test_that("row sums can differ from 1 when sum_to_one_constraint = FALSE", {
    hyp_weight <- c(0.5, 0.5)
    trans_matrix <- matrix(
        c(
            0, 1,
            0, 0 # fixed-sequence style: last row has sum 0
        ),
        nrow = 2,
        byrow = TRUE
    )

    expect_true(
        is_graph_valid(hyp_weight, trans_matrix, sum_to_one_constraint = FALSE)
    )
})


test_that("works for a valid fixed-sequence style graph", {
    hyp_weight <- c(1, 0, 0)
    trans_matrix <- matrix(
        c(
            0, 1, 0,
            0, 0, 1,
            0, 0, 0
        ),
        nrow = 3,
        byrow = TRUE
    )

    expect_true(
        is_graph_valid(hyp_weight, trans_matrix, sum_to_one_constraint = FALSE)
    )

    expect_warning(
        is_graph_valid(hyp_weight, trans_matrix, sum_to_one_constraint = TRUE),
        regexp = "rows of `trans_matrix` do not sum to 1"
    )
    expect_false(
        is_graph_valid(
            hyp_weight,
            trans_matrix,
            sum_to_one_constraint = TRUE
        )
    ) |>
        suppressWarnings()
})


test_that("bullets_with_header", {
    expect_snapshot(
        bullets_with_header("foo", list(x = 1, y = 2))
    )

    expect_snapshot(
        bullets_with_header("foo", list())
    )
})

test_that("as_simple", {
    a <- list(
        x = 1,
        y = list(
            z = letters[1:5],
            u = mtcars[1:5, ]
        ),
        v = "foo"
    )

    expect_identical(
        as_simple(a),
        "<list>"
    )

    expect_identical(
        as_simple(a$x),
        "1"
    )

    expect_identical(
        as_simple(a$y$z),
        "<character>"
    )

    expect_identical(
        as_simple(a$v),
        '\"foo\"'
    )

    expect_identical(
        as_simple(list()),
        "<unset>"
    )
})

test_that("modify_list works", {
    my_list <- list(
        x = 1,
        y = "foo"
    )

    # modify existing elements
    expect_identical(
        modify_list(
            my_list,
            x = 2
        ),
        list(
            y = "foo",
            x = 2
        )
    )

    # add new elements
    expect_identical(
        modify_list(
            my_list,
            bar = 2
        ),
        list(
            x = 1,
            y = "foo",
            bar = 2
        )
    )

    # empty dots return the input object unchanged
    expect_identical(
        modify_list(my_list),
        my_list
    )

    # attempting to set an element to NULL removes it
    expect_identical(
        modify_list(
            my_list,
            x = NULL
        ),
        list(
            y = "foo"
        )
    )
})

test_that("modify_list complains", {
    my_list <- list(
        x = 1,
        y = "foo"
    )

    expect_error(
        modify_list(
            my_list,
            2
        ),
        "All components of `...` must be named."
    )
})

test_that("modify_list strips names from empty output", {
    my_list <- list(
        x = 1,
        y = "foo"
    )

    # names are stripped from an empty output
    expect_identical(
        modify_list(
            my_list,
            x = NULL,
            y = NULL
        ),
        list()
    )
})


# --- normalise_sum tests ---
test_that("normalise_sum errors when fixed elements exceed target", {
    x <- c(0.6, 0.5, 0.1)
    expect_error(
        normalise_sum(x, fixed_idx = c(1L, 2L), target = 1),
        "exceeds target"
    )
})

test_that("normalise_sum basic: result sums to 1", {
    tol <- sqrt(.Machine$double.eps)
    x <- c(0.2, 0.3, 0.500000001)
    result <- normalise_sum(x, tolerance = tol)
    expect_equal(sum(result), 1, tolerance = tol)
    expect_true(all(result >= 0))
})

test_that("normalise_sum is idempotent", {
    tol <- sqrt(.Machine$double.eps)
    x <- c(0.25, 0.25, 0.5)
    once <- normalise_sum(x, tolerance = tol)
    twice <- normalise_sum(once, tolerance = tol)
    expect_equal(once, twice, tolerance = tol)
})

test_that("normalise_sum with fixed_idx preserves fixed elements", {
    tol <- sqrt(.Machine$double.eps)
    x <- c(0.25, 0.45, 0.30)
    fixed <- c(1L, 3L)
    result <- normalise_sum(x, fixed_idx = fixed, tolerance = tol)
    expect_equal(result[1], 0.25, tolerance = tol)
    expect_equal(result[3], 0.30, tolerance = tol)
    expect_equal(sum(result), 1, tolerance = tol)
})

test_that("normalise_sum with all-zero input returns all zeros", {
    tol <- sqrt(.Machine$double.eps)
    x <- c(0, 0, 0)
    result <- normalise_sum(x, tolerance = tol)
    expect_equal(result, c(0, 0, 0), tolerance = tol)
})

test_that("normalise_sum with all elements fixed returns input unchanged", {
    tol <- sqrt(.Machine$double.eps)
    x <- c(0.3, 0.3, 0.4)
    result <- normalise_sum(x, fixed_idx = 1:3, tolerance = tol)
    expect_equal(result, x, tolerance = tol)
})

test_that("normalise_sum with non-default target", {
    tol <- sqrt(.Machine$double.eps)
    x <- c(0.2, 0.3, 0.5)
    result <- normalise_sum(x, target = 2, tolerance = tol)
    expect_equal(sum(result), 2, tolerance = tol)
})

test_that("normalise_sum diagonal protection: fixed_idx includes self-loop", {
    tol <- sqrt(.Machine$double.eps)
    row <- c(0, 0.6, 0.4)
    fixed <- 1L
    result <- normalise_sum(row, fixed_idx = fixed, tolerance = tol)
    expect_equal(result[1], 0, tolerance = tol)
    expect_equal(sum(result), 1, tolerance = tol)
})

test_that("normalise_sum additive fallback corrects tiny residual", {
    # Vector where proportional scaling + complement leaves a tiny residual
    x <- c(1 / 3, 1 / 3, 1 / 3 + 1e-16)
    result <- normalise_sum(x)
    expect_equal(sum(result), 1, tolerance = sqrt(.Machine$double.eps))
})

test_that("normalise_sum clamps negative anchor to zero", {
    # Fixed elements nearly exhaust the target, free elements are tiny
    # After scaling, the anchor could go negative due to rounding
    x <- c(0.99, 0.009, 0.001)
    result <- normalise_sum(x, fixed_idx = 1L, target = 1)
    expect_true(all(result >= 0))
    expect_equal(sum(result), 1, tolerance = sqrt(.Machine$double.eps))
})

test_that("normalise_sum stress test: 1000 random vectors", {
    tol <- sqrt(.Machine$double.eps)
    set.seed(42)
    expect_no_error(
        for (i in seq_len(1000)) {
            n <- sample(2:10, 1)
            x <- runif(n)
            result <- normalise_sum(x, tolerance = tol)
            stopifnot(
                all.equal(sum(result), 1, tolerance = tol),
                result >= 0
            )
        }
    )
})


test_that("normalise_sum stress test with fixed_idx: 1000 random vectors", {
    tol <- sqrt(.Machine$double.eps)
    set.seed(123)
    expect_no_error(
        for (i in seq_len(1000)) {
            n <- sample(3:10, 1)
            x <- runif(n)
            x <- x / sum(x) # start from a valid probability vector
            n_fixed <- sample(0:(n - 1), 1)
            fixed_idx <- if (n_fixed > 0) {
                sort(sample.int(n, n_fixed))
            } else {
                integer(0)
            }
            free_idx <- setdiff(seq_len(n), fixed_idx)
            # Only perturb free elements so fixed sum stays < 1
            x[free_idx] <- x[free_idx] * runif(length(free_idx), 0.8, 1.2)
            result <- normalise_sum(x, fixed_idx = fixed_idx, tolerance = tol)
            stopifnot(
                all.equal(sum(result), 1, tolerance = tol),
                result >= 0
            )
            if (length(fixed_idx) > 0) {
                stopifnot(all.equal(
                    result[fixed_idx],
                    x[fixed_idx],
                    tolerance = tol
                ))
            }
        }
    )
})


test_that("normalise_sum problematic vector test", {
    problematic_vector <- readRDS(
        test_path(
            "data",
            "normalise_exactly_vec.rds"
        )
    )

    # Problematic vector only doesn't sum strictly on macOS.
    skip_on_os(c("linux", "windows"))
    expect_false(sum(problematic_vector) == 1)

    expect_no_warning(
        normalised_vec <- normalise_sum(problematic_vector)
    )

    expect_equal(sum(normalised_vec), 1, tolerance = sqrt(.Machine$double.eps))
})

test_that("normalise_sum complains when anything is passed via `...`", {
    # After scaling, the anchor could go negative due to rounding
    x <- c(0.99, 0.009, 0.001)
    expect_snapshot(error = TRUE, {
        normalise_sum(x, 1L)
    })
})

test_that("calc_ncp", {
    set.seed(1)
    expect_snapshot({
        calc_ncp(power = c(0.8, 0.9))
    })

    set.seed(2)
    expect_snapshot({
        calc_ncp(power = c(0.8, 0.9), alpha = 0.01)
    })
})
