# ---- Helper: minimal S4 mock for mutation tests ----
# The mutation closure only accesses @population.
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

make_mock_ga <- function(population_matrix) {
    d <- ncol(population_matrix)
    new(
        "mock_ga",
        population = population_matrix,
        lower = rep(0, d),
        upper = rep(1, d)
    )
}

# ---- Helper: constraints and parents ----

# Free 4-hypothesis graph: 3 weight parameters, 2 per transition row.
free_4 <- graph_constraint_free(4)
free_4_h <- unname(free_4$hyp_constraint)
free_4_t <- unname(free_4$trans_constraint)

# Constrained 4-hypothesis graph: a fixed weight, a fully fixed row and a row
# with one fixed entry. Groups: weights (latent 1:2, c = 0.8), row 1 (3:4),
# row 3 (5, c = 0.7) and row 4 (6:7).
pinned_h <- c(NA, NA, 0.2, NA)
pinned_t <- matrix(NA_real_, 4, 4)
diag(pinned_t) <- 0
pinned_t[2, ] <- c(1, 0, 0, 0)
pinned_t[3, 4] <- 0.3

decode_graph <- function(x, hyp_constraint, trans_constraint) {
    theta <- split_theta(x, hyp_constraint)
    list(
        hyp_weight = recover_full_weights(theta$w_pars, hyp_constraint),
        trans_matrix = recover_full_trans_matrix(theta$g_pars, trans_constraint)
    )
}

# A valid parent with some empty entries in every group.
random_parent <- function(groups, n_latent) {
    x <- numeric(n_latent)
    for (g in groups) {
        k <- length(g$latent) + 1L
        s <- stats::rexp(k)
        s[stats::runif(k) < 0.3] <- 0
        if (all(s == 0)) {
            s[sample.int(k, 1L)] <- 1
        }
        x <- .sm_from_shares(x, g, s / sum(s))
    }
    x
}

# ---- .reflect01 tests ----

test_that(".reflect01 folds any number into [0, 1] and is the identity inside", {
    v <- c(-1e6, -3.7, -1, -0.25, 0, 1, 1.25, 2, 2.5, 7.3, 1e6)
    out <- .reflect01(v)
    expect_true(all(out >= 0 & out <= 1))
    expect_equal(.reflect01(c(-0.25, 1.25, 2.5)), c(0.25, 0.75, 0.5))

    inside <- seq(0.001, 0.999, length.out = 50)
    expect_identical(.reflect01(inside), inside)

    expect_true(is.nan(.reflect01(Inf)))
})

# ---- .simplex_project tests ----

test_that(".simplex_project lands on the simplex and fixes points on it", {
    v <- c(0.7, -0.2, 0.6, 0.1)
    p <- .simplex_project(v)
    expect_true(all(p >= 0))
    expect_equal(sum(p), 1)

    on_simplex <- c(0.5, 0.3, 0.2, 0)
    expect_equal(.simplex_project(on_simplex), on_simplex)
})

# ---- .sm_groups tests ----

test_that(".sm_groups follows the decoder's layout", {
    groups <- .sm_groups(pinned_h, pinned_t)

    expect_length(groups, 4)
    expect_identical(lapply(groups, `[[`, "latent"), list(1:2, 3:4, 5L, 6:7))
    expect_equal(vapply(groups, `[[`, numeric(1), "c"), c(0.8, 1, 0.7, 1))
    expect_identical(
        vapply(groups, `[[`, numeric(1), "threshold"),
        c(1e-4, 1e-5, 1e-5, 1e-5)
    )
})

test_that(".sm_groups skips groups with no free mass but keeps their positions", {
    h <- c(1, 0, 0, 0)
    t <- matrix(0, 4, 4)
    t[1, ] <- c(0, 1, NA, NA) # two free entries, c = 0: one latent, no group
    t[4, ] <- c(NA, NA, NA, 0) # three free entries, c = 1

    groups <- .sm_groups(h, t)

    expect_length(groups, 1)
    expect_identical(groups[[1]]$latent, 2:3)
    expect_equal(groups[[1]]$c, 1)
})

# ---- .sm_merge tests ----

test_that("a merge empties exactly one active entry, the smaller of its pair", {
    set.seed(42)
    g <- list(latent = 1:5, c = 1, threshold = 1e-5)
    s <- c(0.4, 0.3, 0.2, 0.1, 0, 0)

    for (i in 1:200) {
        out <- .sm_merge(g, s)
        changed <- which(out != s)
        expect_length(changed, 2)
        emptied <- changed[out[changed] == 0]
        expect_length(emptied, 1)
        expect_true(s[emptied] > 0)
        expect_true(s[emptied] <= max(s[changed]))
        expect_equal(sum(out), sum(s))
    }
})

# ---- .sm_transfer tests ----

test_that("a transfer keeps the pair's total and can add only its recipient", {
    set.seed(42)
    g <- list(latent = 1:4, c = 1, threshold = 1e-5)
    s <- c(0.5, 0.3, 0.2, 0, 0)
    active_before <- which(s > 0)

    for (i in 1:500) {
        out <- .sm_transfer(g, s, scale = 0.2)
        changed <- which(out != s)
        expect_lte(length(changed), 2)
        expect_true(all(out[changed] > 0))
        expect_equal(sum(out), sum(s))
        # At most one entry becomes active, and none is lost
        expect_lte(length(setdiff(which(out > 0), active_before)), 1)
        expect_true(all(active_before %in% which(out > 0)))
    }
})

# ---- .make_split_merge_mutation tests ----

test_that(".make_split_merge_mutation returns a closure", {
    mut_fn <- .make_split_merge_mutation(free_4_h, free_4_t)
    expect_type(mut_fn, "closure")
})

test_that(".make_split_merge_mutation checks its settings", {
    expect_error(.make_split_merge_mutation(free_4_h, free_4_t, p_del = 1.5))
    expect_error(.make_split_merge_mutation(free_4_h, free_4_t, p_del = -0.1))
    expect_error(.make_split_merge_mutation(free_4_h, free_4_t, scale = 0))
    expect_error(.make_split_merge_mutation(free_4_h, free_4_t, scale = c(1, 2)))
})

test_that("SM keeps the graph valid, keeps fixed entries and changes one group", {
    set.seed(42)
    cases <- list(
        list(h = free_4_h, t = free_4_t, d = 11L),
        list(h = pinned_h, t = pinned_t, d = 7L)
    )

    for (cs in cases) {
        groups <- .sm_groups(cs$h, cs$t)
        mut_fn <- .make_split_merge_mutation(cs$h, cs$t)
        pop <- t(replicate(300, random_parent(groups, cs$d)))
        obj <- make_mock_ga(pop)

        for (r in seq_len(nrow(pop))) {
            child <- mut_fn(obj, r)
            graph <- decode_graph(child, cs$h, cs$t)

            expect_length(child, cs$d)
            expect_true(is_graph_valid(graph$hyp_weight, graph$trans_matrix))
            expect_identical(
                graph$hyp_weight[!is.na(cs$h)],
                cs$h[!is.na(cs$h)]
            )
            expect_identical(
                graph$trans_matrix[!is.na(cs$t)],
                cs$t[!is.na(cs$t)]
            )

            n_changed <- sum(vapply(
                groups,
                function(g) !identical(child[g$latent], pop[r, g$latent]),
                logical(1)
            ))
            expect_identical(n_changed, 1L)
        }
    }
})

test_that("SM never returns its parent, at any p_del", {
    set.seed(42)
    groups <- .sm_groups(free_4_h, free_4_t)
    parent <- random_parent(groups, 11L)
    obj <- make_mock_ga(matrix(parent, nrow = 1))

    for (p_del in c(0, 0.45, 1)) {
        mut_fn <- .make_split_merge_mutation(free_4_h, free_4_t, p_del = p_del)
        n_null <- sum(replicate(1000, identical(mut_fn(obj, 1L), parent)))
        expect_identical(n_null, 0L)
    }
})

test_that("same seed gives the same child", {
    groups <- .sm_groups(pinned_h, pinned_t)
    set.seed(1)
    obj <- make_mock_ga(matrix(random_parent(groups, 7L), nrow = 1))
    mut_fn <- .make_split_merge_mutation(pinned_h, pinned_t)

    set.seed(99)
    first <- mut_fn(obj, 1L)
    set.seed(99)
    second <- mut_fn(obj, 1L)

    expect_identical(first, second)
})

test_that("SM picks positions uniformly: the derived entry is not special", {
    set.seed(42)
    h <- rep(NA_real_, 6)
    t <- matrix(0, 6, 6) # transitions fixed, so only the weights mutate
    g <- .sm_groups(h, t)[[1]]
    parent <- .sm_from_shares(numeric(5), g, rep(1 / 6, 6))
    obj <- make_mock_ga(matrix(parent, nrow = 1))

    for (p_del in c(0, 1)) {
        mut_fn <- .make_split_merge_mutation(h, t, p_del = p_del)
        hits <- integer(6)
        for (i in 1:6000) {
            child <- mut_fn(obj, 1L)
            # Rounding can move the derived entry by an ulp; ignore that
            step <- .sm_free_entries(child, g) - .sm_free_entries(parent, g)
            hits <- hits + (abs(step) > 1e-12)
        }
        expect_true(all(abs(hits / mean(hits) - 1) < 0.15))
    }
})

test_that("scale changes only the size of a transfer, not which entries move", {
    groups <- .sm_groups(free_4_h, free_4_t)
    set.seed(1)
    parent <- random_parent(groups, 11L)
    obj <- make_mock_ga(matrix(parent, nrow = 1))

    small_fn <- .make_split_merge_mutation(
        free_4_h,
        free_4_t,
        p_del = 0,
        scale = 0.01
    )
    large_fn <- .make_split_merge_mutation(
        free_4_h,
        free_4_t,
        p_del = 0,
        scale = 1
    )

    small_steps <- large_steps <- numeric(200)
    for (i in 1:200) {
        set.seed(i)
        small <- small_fn(obj, 1L)
        set.seed(i)
        large <- large_fn(obj, 1L)
        expect_identical(
            which(abs(small - parent) > 1e-12),
            which(abs(large - parent) > 1e-12)
        )
        small_steps[i] <- sum(abs(small - parent))
        large_steps[i] <- sum(abs(large - parent))
    }

    expect_lt(mean(small_steps), mean(large_steps))
})

test_that("SM repairs an invalid parent into a valid child", {
    set.seed(42)
    # Every group's stored parameters sum to more than its free mass
    parent <- rep(0.9, 11)
    obj <- make_mock_ga(matrix(parent, nrow = 1))
    mut_fn <- .make_split_merge_mutation(free_4_h, free_4_t)

    for (i in 1:50) {
        graph <- decode_graph(mut_fn(obj, 1L), free_4_h, free_4_t)
        expect_true(is_graph_valid(graph$hyp_weight, graph$trans_matrix))
    }
})

test_that("SM returns the parent when no group can be mutated", {
    h <- c(1, 0, 0, 0)
    t <- matrix(0, 4, 4)
    t[1, ] <- c(0, 1, NA, NA)
    parent <- 0
    obj <- make_mock_ga(matrix(parent, nrow = 1))
    mut_fn <- .make_split_merge_mutation(h, t)

    expect_identical(mut_fn(obj, 1L), parent)
})

test_that("mutation selects correct parent from multi-row population", {
    groups <- .sm_groups(free_4_h, free_4_t)
    set.seed(1)
    pop <- rbind(random_parent(groups, 11L), random_parent(groups, 11L))
    obj <- make_mock_ga(pop)
    mut_fn <- .make_split_merge_mutation(free_4_h, free_4_t)

    child <- mut_fn(obj, 2L)

    # One group moves; the other four are those of parent 2
    n_kept <- sum(vapply(
        groups,
        function(g) identical(child[g$latent], pop[2, g$latent]),
        logical(1)
    ))
    expect_identical(n_kept, length(groups) - 1L)
})
