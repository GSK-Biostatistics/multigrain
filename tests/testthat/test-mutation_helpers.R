# ---- .cauchy_perturb tests ----

test_that(".cauchy_perturb returns numeric scalar", {
    set.seed(42)
    val <- .cauchy_perturb(0.5, 0, 1)
    expect_type(val, "double")
    expect_length(val, 1)
})

test_that(".cauchy_perturb returns value within bounds", {
    set.seed(42)
    vals <- replicate(500, .cauchy_perturb(0.5, 0, 1, scale = 1.0))
    expect_true(all(vals >= 0 & vals <= 1))
})

test_that(".cauchy_perturb respects asymmetric bounds", {
    set.seed(42)
    vals <- replicate(500, .cauchy_perturb(0.3, 0.2, 0.4, scale = 1.0))
    expect_true(all(vals >= 0.2 & vals <= 0.4))
})

test_that(".cauchy_perturb perturbs around centre, not uniformly", {
    set.seed(42)
    vals <- replicate(500, .cauchy_perturb(0.5, 0, 1, scale = 0.01))
    expect_lt(mean(abs(vals - 0.5)), 0.05)
})

test_that(".cauchy_perturb falls back to uniform under extreme rejection", {
    set.seed(42)
    vals <- replicate(50, .cauchy_perturb(0.0, 0, 1, scale = 100))
    expect_true(all(vals >= 0 & vals <= 1))
})

test_that(".cauchy_perturb handles centre at exact boundary", {
    set.seed(42)
    vals_low <- replicate(500, .cauchy_perturb(0.0, 0, 1, scale = 1.0))
    expect_true(all(vals_low >= 0 & vals_low <= 1))

    vals_high <- replicate(500, .cauchy_perturb(1.0, 0, 1, scale = 1.0))
    expect_true(all(vals_high >= 0 & vals_high <= 1))
})

test_that(".cauchy_perturb handles degenerate zero-width interval", {
    val <- .cauchy_perturb(0.5, 0.5, 0.5, scale = 1.0)
    expect_equal(val, 0.5)
})

# ---- Helper: minimal S4 mock for mutation tests ----
# The mutation closure only accesses @population, @lower, @upper.
# Defining a local class avoids depending on GA's S4 class registration
# (which requires GA to be *attached*, not just imported).

setClass(
    "mock_ga",
    representation(
        population = "matrix",
        lower = "numeric",
        upper = "numeric"
    )
)

make_mock_ga <- function(population_matrix, lower, upper) {
    new("mock_ga", population = population_matrix, lower = lower, upper = upper)
}

# ---- .make_cauchy_mutation_multi tests ----

test_that(".make_cauchy_mutation_multi returns a closure", {
    mut_fn <- .make_cauchy_mutation_multi(p_param_mutate = 0.1, scale = 1.0)
    expect_type(mut_fn, "closure")
})

test_that("mutation output has correct length and respects bounds", {
    set.seed(42)
    d <- 10
    pop <- matrix(rep(0.5, d), nrow = 1)
    obj <- make_mock_ga(pop, lower = rep(0, d), upper = rep(1, d))

    mut_fn <- .make_cauchy_mutation_multi(p_param_mutate = 0.5, scale = 1.0)
    result <- mut_fn(obj, 1L)

    expect_length(result, d)
    expect_true(all(result >= 0 & result <= 1))
})

test_that("at least one parameter is always mutated", {
    set.seed(42)
    d <- 20
    parent <- rep(0.5, d)
    pop <- matrix(parent, nrow = 1)
    obj <- make_mock_ga(pop, lower = rep(0, d), upper = rep(1, d))

    mut_fn <- .make_cauchy_mutation_multi(p_param_mutate = 1e-10, scale = 1.0)

    expect_no_error(
        for (i in 1:30) {
            result <- mut_fn(obj, 1L)
            n_changed <- sum(result != parent)
            stopifnot(n_changed >= 1)
        }
    )
})

test_that("p_param_mutate = 1.0 mutates all parameters", {
    set.seed(42)
    d <- 20
    parent <- rep(0.5, d)
    pop <- matrix(parent, nrow = 1)
    obj <- make_mock_ga(pop, lower = rep(0, d), upper = rep(1, d))

    mut_fn <- .make_cauchy_mutation_multi(p_param_mutate = 1.0, scale = 1.0)

    changed_counts <- replicate(10, {
        r <- mut_fn(obj, 1L)
        sum(r != parent)
    })
    expect_gt(mean(changed_counts), d * 0.75)
})

test_that("mutation is perturbation-based, not replacement-based", {
    set.seed(42)
    d <- 20
    parent <- rep(0.5, d)
    pop <- matrix(parent, nrow = 1)
    obj <- make_mock_ga(pop, lower = rep(0, d), upper = rep(1, d))

    mut_fn <- .make_cauchy_mutation_multi(p_param_mutate = 1.0, scale = 0.001)
    result <- mut_fn(obj, 1L)

    expect_true(all(abs(result - parent) < 0.1))
})

test_that("scale parameter affects perturbation magnitude", {
    set.seed(42)
    d <- 50
    parent <- rep(0.5, d)
    pop <- matrix(parent, nrow = 1)
    obj <- make_mock_ga(pop, lower = rep(0, d), upper = rep(1, d))

    small_fn <- .make_cauchy_mutation_multi(p_param_mutate = 1.0, scale = 0.001)
    large_fn <- .make_cauchy_mutation_multi(p_param_mutate = 1.0, scale = 1.0)

    small_diffs <- replicate(20, {
        r <- small_fn(obj, 1L)
        mean(abs(r - parent))
    })
    large_diffs <- replicate(20, {
        r <- large_fn(obj, 1L)
        mean(abs(r - parent))
    })

    expect_lt(mean(small_diffs), mean(large_diffs))
})

test_that("mutation selects correct parent from multi-row population", {
    set.seed(42)

    d <- 5
    pop <- matrix(
        c(
            rep(0.1, d),
            rep(0.9, d)
        ),
        nrow = 2,
        byrow = TRUE
    )
    obj <- make_mock_ga(pop, lower = rep(0, d), upper = rep(1, d))

    mut_fn <- .make_cauchy_mutation_multi(p_param_mutate = 1.0, scale = 0.001)

    result <- mut_fn(obj, 2L)
    expect_true(all(abs(result - 0.9) < 0.1))
})

test_that("mutation preserves unmutated parameters exactly", {
    set.seed(42)
    d <- 20
    parent <- seq(0.05, 0.95, length.out = d)
    pop <- matrix(parent, nrow = 1)
    obj <- make_mock_ga(pop, lower = rep(0, d), upper = rep(1, d))

    mut_fn <- .make_cauchy_mutation_multi(p_param_mutate = 1e-10, scale = 1.0)
    result <- mut_fn(obj, 1L)

    n_changed <- as.integer(sum(result != parent))
    expect_identical(n_changed, 1L)
    unchanged <- which(result == parent)
    expect_identical(result[unchanged], parent[unchanged])
})

test_that("mutation with boundary parent values stays in bounds", {
    set.seed(42)
    d <- 10
    parent <- c(0, 0, 0, 1, 1, 1, 0.5, 0.5, 0, 1)
    pop <- matrix(parent, nrow = 1)
    obj <- make_mock_ga(pop, lower = rep(0, d), upper = rep(1, d))

    mut_fn <- .make_cauchy_mutation_multi(p_param_mutate = 1.0, scale = 1.0)

    expect_no_error(
        for (i in 1:100) {
            result <- mut_fn(obj, 1L)
            stopifnot(result >= 0 & result <= 1)
        }
    )
})


# ---- .g_param_rows / .zeroing_row_targets ----

test_that(".g_param_rows matches recover_full_trans_matrix parameter order", {
    # A constrained example: row 1 has a pinned non-zero entry alongside free
    # ones, rows 3 and 4 are fully fixed and contribute no parameters.
    hc <- c(NA, NA, 0, 0)
    tc <- rbind(
        c(0, NA, NA, 0),
        c(NA, 0, 0, NA),
        c(0, 1, 0, 0),
        c(1, 0, 0, 0)
    )
    gc <- graph_constraint(hyp_constraint = hc, trans_constraint = tc)

    param_rows <- .g_param_rows(gc$hyp_constraint, gc$trans_constraint)

    # 1 weight parameter (2 free weights - 1), then 1 parameter each for rows
    # 1 and 2, and none for the fully fixed rows 3 and 4.
    expect_identical(param_rows, c(NA_integer_, 1L, 2L))
    expect_length(param_rows, length(create_start_params(gc)))

    # Each G parameter must land in the row .g_param_rows() claims: give each
    # one a unique marker and find it in the decoded matrix.
    g_idx <- which(!is.na(param_rows))
    for (p in seq_along(g_idx)) {
        x <- rep(0.2, length(param_rows))
        marker <- 0.371
        x[g_idx[p]] <- marker
        theta <- split_theta(x, gc$hyp_constraint)
        G_dec <- recover_full_trans_matrix(theta$g_pars, gc$trans_constraint)
        expect_identical(
            which(abs(G_dec - marker) < 1e-12, arr.ind = TRUE)[1, "row"],
            param_rows[g_idx[p]],
            ignore_attr = TRUE
        )
    }
})

test_that(".g_param_rows handles a fully fixed weight vector", {
    hc <- c(0.5, 0.5, 0, 0)
    gc <- graph_constraint(
        hyp_constraint = hc,
        trans_constraint = trans_constraint_free(4)
    )
    param_rows <- .g_param_rows(gc$hyp_constraint, gc$trans_constraint)
    # No free weights, so no leading NA entries.
    expect_false(anyNA(param_rows))
    expect_identical(param_rows, rep(1:4, each = 2L))
})

test_that(".zeroing_row_targets subtracts the row's pinned entries", {
    tc <- rbind(
        c(0, 0.5, NA, NA),
        c(NA, 0, NA, NA),
        c(NA, NA, 0, NA),
        c(NA, NA, NA, 0)
    )
    targets <- .zeroing_row_targets(tc)
    expect_equal(targets[1], 1 - 5e-6 - 0.5, tolerance = 1e-15)
    expect_equal(targets[2:4], rep(1 - 5e-6, 3), tolerance = 1e-15)

    # A row whose fixed entries already sum to one has no room for the move.
    tc_full <- rbind(c(0, 1, 0, 0), c(NA, 0, NA, NA), c(NA, NA, 0, NA),
        c(NA, NA, NA, 0))
    expect_lt(.zeroing_row_targets(tc_full)[1], 0)
})


# ---- the zeroing move ----

# Verbatim copy of the mutation closure as it stood before `p_zero` was added.
# Used to prove that p_zero = 0 consumes the same random numbers.
# nolint start: object_usage_linter. `.cauchy_perturb` is in the package
# namespace at test time, but codetools cannot see it from this file.
cauchy_mutation_pre_change <- function(p_param_mutate = 0.1, scale = 1.0) {
    force(p_param_mutate)
    force(scale)
    function(object, parent) {
        parent_vec <- as.numeric(object@population[parent, ])
        d <- length(parent_vec)
        lower <- object@lower
        upper <- object@upper

        mutate_mask <- stats::runif(d) < p_param_mutate

        if (!any(mutate_mask)) {
            mutate_mask[sample.int(d, 1L)] <- TRUE
        }

        for (j in which(mutate_mask)) {
            parent_vec[j] <- .cauchy_perturb(
                parent_vec[j],
                lower[j],
                upper[j],
                scale
            )
        }

        parent_vec
    }
}
# nolint end

test_that("p_zero = 0 draws the same random numbers as before the move", {
    gc <- graph_constraint_free(4)
    param_rows <- .g_param_rows(gc$hyp_constraint, gc$trans_constraint)
    d <- length(param_rows)
    pop <- matrix(seq(0.05, 0.95, length.out = d), nrow = 1)
    obj <- make_mock_ga(pop, lower = rep(0, d), upper = rep(1, d))

    new_fn <- .make_cauchy_mutation_multi(
        p_param_mutate = 0.1,
        scale = 1.0,
        p_zero = 0,
        param_rows = param_rows,
        row_target = .zeroing_row_targets(gc$trans_constraint)
    )
    old_fn <- cauchy_mutation_pre_change(p_param_mutate = 0.1, scale = 1.0)

    new_out <- withr::with_seed(99, lapply(1:50, function(i) new_fn(obj, 1L)))
    old_out <- withr::with_seed(99, lapply(1:50, function(i) old_fn(obj, 1L)))
    expect_identical(new_out, old_out)

    # ... and the RNG state afterwards is the same, so nothing downstream of a
    # mutation call shifts either.
    seed_new <- withr::with_seed(99, {
        for (i in 1:50) new_fn(obj, 1L)
        .Random.seed
    })
    seed_old <- withr::with_seed(99, {
        for (i in 1:50) old_fn(obj, 1L)
        .Random.seed
    })
    expect_identical(seed_new, seed_old)
})

test_that("p_zero = 0 returns the plain Cauchy closure unwrapped", {
    fn <- .make_cauchy_mutation_multi(p_zero = 0)
    # The body is the Cauchy closure itself: no coin is drawn.
    expect_false(grepl("p_zero", deparse1(body(fn)), fixed = TRUE))
})

test_that("p_zero = 1 always zeroes a parameter or a row's derived entry", {
    gc <- graph_constraint_free(4)
    param_rows <- .g_param_rows(gc$hyp_constraint, gc$trans_constraint)
    row_target <- .zeroing_row_targets(gc$trans_constraint)
    d <- length(param_rows)
    g_idx <- which(!is.na(param_rows))

    # Every parameter is above 1e-5 and every row sums positive, so the move
    # can never fall through to the Cauchy branch.
    parent <- rep(0.3, d)
    pop <- matrix(parent, nrow = 1)
    obj <- make_mock_ga(pop, lower = rep(0, d), upper = rep(1, d))

    fn <- .make_cauchy_mutation_multi(
        p_param_mutate = 0.1,
        scale = 1.0,
        p_zero = 1,
        param_rows = param_rows,
        row_target = row_target
    )

    variants <- withr::with_seed(11, vapply(1:200, function(i) {
        out <- fn(obj, 1L)
        changed <- which(out != parent)

        zeroed <- length(changed) == 1L &&
            changed %in% g_idx &&
            identical(out[changed], 0)
        if (zeroed) {
            return("zero_param")
        }

        rows_touched <- unique(param_rows[changed])
        rescaled <- length(rows_touched) == 1L &&
            !is.na(rows_touched) &&
            !anyNA(param_rows[changed]) &&
            abs(
                sum(out[g_idx[param_rows[g_idx] == rows_touched]]) -
                    row_target[rows_touched]
            ) < 1e-12
        if (rescaled) {
            return("zero_derived")
        }

        "other"
    }, character(1)))

    expect_setequal(unique(variants), c("zero_param", "zero_derived"))
    # Roughly an equal split between the two variants.
    expect_gt(min(table(variants)), 60L)
})

test_that("the derived-entry move respects a row's pinned non-zero entry", {
    # A11: row 1 pins G[1, 2] = 0.5 and leaves columns 3 and 4 free, so the
    # rescale target must be 1 - 5e-6 - 0.5, not 1 - 5e-6.
    tc <- rbind(
        c(0, 0.5, NA, NA),
        c(NA, 0, NA, NA),
        c(NA, NA, 0, NA),
        c(NA, NA, NA, 0)
    )
    gc <- graph_constraint(
        hyp_constraint = c(NA, NA, NA, NA),
        trans_constraint = tc
    )
    param_rows <- .g_param_rows(gc$hyp_constraint, gc$trans_constraint)
    row_target <- .zeroing_row_targets(gc$trans_constraint)
    d <- length(param_rows)
    g_idx <- which(!is.na(param_rows))
    row1 <- g_idx[param_rows[g_idx] == 1L]

    parent <- rep(0.3, d)
    pop <- matrix(parent, nrow = 1)
    obj <- make_mock_ga(pop, lower = rep(0, d), upper = rep(1, d))

    fn <- .make_cauchy_mutation_multi(
        p_param_mutate = 0.1,
        scale = 1.0,
        p_zero = 1,
        param_rows = param_rows,
        row_target = row_target
    )

    seen <- 0L
    withr::with_seed(23, {
        for (i in 1:400) {
            out <- fn(obj, 1L)
            touched <- which(out != parent)
            # Row 1 has a single parameter, so zeroing it and
            # rescaling it touch the same index; the target sum is
            # what tells the derived-entry move apart.
            if (!setequal(touched, row1) ||
                abs(sum(out[row1]) - row_target[1]) > 1e-12) {
                next
            }
            seen <- seen + 1L
            theta <- split_theta(out, gc$hyp_constraint)
            G_dec <- recover_full_trans_matrix(
                theta$g_pars,
                gc$trans_constraint
            )
            # Derived entry lands just under the 1e-5 zeroing threshold, and
            # stays positive so the objective does not penalise it.
            # Absolute, not relative: at 5e-6 a relative tolerance of
            # 1e-12 is tighter than double precision allows.
            expect_lt(abs(G_dec[1, 4] - 5e-6), 1e-12)
            expect_true(all(G_dec >= 0))
            expect_identical(G_dec[1, 2], 0.5)
        }
    })
    expect_gt(seen, 0L)
})

test_that("the zeroing move falls through when there is nothing to zero", {
    gc <- graph_constraint_free(4)
    param_rows <- .g_param_rows(gc$hyp_constraint, gc$trans_constraint)
    d <- length(param_rows)
    g_idx <- which(!is.na(param_rows))

    # All transition parameters already zero: no candidate to zero, and every
    # row sums to zero, so both branches must fall back to Cauchy.
    parent <- rep(0, d)
    parent[is.na(param_rows)] <- 0.3
    pop <- matrix(parent, nrow = 1)
    obj <- make_mock_ga(pop, lower = rep(0, d), upper = rep(1, d))

    fn <- .make_cauchy_mutation_multi(
        p_param_mutate = 1.0,
        scale = 1.0,
        p_zero = 1,
        param_rows = param_rows,
        row_target = .zeroing_row_targets(gc$trans_constraint)
    )

    withr::with_seed(5, {
        for (i in 1:50) {
            out <- fn(obj, 1L)
            expect_length(out, d)
            expect_true(all(out >= 0 & out <= 1))
            # A Cauchy perturbation moved something; a zeroing move could not.
            expect_gt(sum(out != parent), 0L)
        }
    })
})
