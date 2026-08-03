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
    pop <- matrix(c(
    rep(0.1, d),
    rep(0.9, d)
  ), nrow = 2, byrow = TRUE)
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
