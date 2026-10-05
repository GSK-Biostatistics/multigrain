# The canonical joint model of the design record, section 4.7:
#   E[Z_ik]            = Delta_i sqrt(t_ik)
#   Corr(Z_ik, Z_jl)   = rho_ij sqrt(min(t_ik, t_jl) / max(t_ik, t_jl))

test_that("dimensions and the info_frac attribute follow a vector info_frac", {
    corr <- matrix(c(1, 0.5, 0.5, 1), nrow = 2L)

    set.seed(1)
    pvals <- simulate_pvalues_gsd(
        c(0.8, 0.9),
        corr_matrix = corr,
        info_frac = c(1 / 3, 2 / 3, 1),
        nsim = 500
    )

    expect_identical(dim(pvals), c(500L, 2L, 3L))
    expect_type(pvals, "double")
    expect_identical(
        attr(pvals, "info_frac"),
        matrix(c(1 / 3, 2 / 3, 1), nrow = 2L, ncol = 3L, byrow = TRUE)
    )
    expect_true(all(pvals > 0 & pvals < 1))
})

test_that("dimensions and the info_frac attribute follow a matrix info_frac", {
    corr <- matrix(c(1, 0.5, 0.5, 1), nrow = 2L)
    info <- rbind(c(1 / 3, 2 / 3, 1), c(NA, 0.5, 1))

    set.seed(1)
    pvals <- simulate_pvalues_gsd(
        c(0.8, 0.9),
        corr_matrix = corr,
        info_frac = info,
        nsim = 500
    )

    expect_identical(dim(pvals), c(500L, 2L, 3L))
    expect_identical(attr(pvals, "info_frac"), info)
})

test_that("empirical moments match the canonical joint model", {
    skip_on_cran()

    # the one run above nsim = 1e4 in this phase: the gate asks for 5e-3 on
    # correlations, which needs a Monte Carlo standard error near 1e-3
    nsim <- 1e6
    rho <- 0.5
    power <- c(0.8, 0.9)
    corr <- matrix(c(1, rho, rho, 1), nrow = 2L)
    # H1 is analysed three times; H2 matures at analysis 2
    info <- rbind(c(1 / 3, 2 / 3, 1), c(0.5, 1, 1))

    set.seed(20260916)
    pvals <- simulate_pvalues_gsd(
        power,
        corr_matrix = corr,
        info_frac = info,
        nsim = nsim
    )

    z <- stats::qnorm(pvals, lower.tail = FALSE)
    dim(z) <- c(nsim, 6L)
    hyp <- rep(1:2, times = 3L)
    look <- rep(1:3, each = 2L)
    info_vec <- info[cbind(hyp, look)]

    expected_corr <- corr[cbind(rep(hyp, 6L), rep(hyp, each = 6L))] *
        sqrt(
            outer(info_vec, info_vec, pmin) / outer(info_vec, info_vec, pmax)
        )
    dim(expected_corr) <- c(6L, 6L)

    corr_error <- max(abs(stats::cor(z) - expected_corr))
    expect_lt(corr_error, 5e-3)

    # within-hypothesis sqrt(t_k / t_l) and cross-hypothesis rho sqrt(min/max)
    # are the two halves of that matrix; assert the headline entries directly
    expect_equal(stats::cor(z[, 1L], z[, 3L]), sqrt(0.5), tolerance = 5e-3)
    expect_equal(stats::cor(z[, 1L], z[, 5L]), sqrt(1 / 3), tolerance = 5e-3)
    expect_equal(
        stats::cor(z[, 1L], z[, 2L]),
        rho * sqrt((1 / 3) / 0.5),
        tolerance = 5e-3
    )
    expect_equal(
        stats::cor(z[, 3L], z[, 6L]),
        rho * sqrt(2 / 3),
        tolerance = 5e-3
    )

    mean_error <- max(abs(
        colMeans(z) - calc_ncp(power, alpha = 0.025)[hyp] * sqrt(info_vec)
    ))
    expect_lt(mean_error, 5e-3)
})

test_that("a matured hypothesis has identical columns after maturity", {
    corr <- matrix(c(1, 0.5, 0.5, 1), nrow = 2L)

    set.seed(2)
    pvals <- simulate_pvalues_gsd(
        c(0.8, 0.9),
        corr_matrix = corr,
        info_frac = rbind(c(1 / 3, 2 / 3, 1), c(0.5, 1, 1)),
        nsim = 500
    )

    expect_identical(pvals[, 2L, 2L], pvals[, 2L, 3L])
    expect_false(identical(pvals[, 1L, 2L], pvals[, 1L, 3L]))
})

test_that("an analysis without data is NA, and the transform accepts it", {
    corr <- matrix(c(1, 0.5, 0.5, 1), nrow = 2L)

    set.seed(3)
    pvals <- simulate_pvalues_gsd(
        c(0.8, 0.9),
        corr_matrix = corr,
        info_frac = rbind(c(1 / 3, 2 / 3, 1), c(NA, 0.5, 1)),
        nsim = 200
    )

    expect_true(all(is.na(pvals[, 2L, 1L])))
    expect_false(anyNA(pvals[, 2L, 2:3]))
    expect_false(anyNA(pvals[, 1L, ]))

    # the info_frac attribute is what `transform_pvalues_gsd()` reads, so the
    # transform needs no information fractions of its own
    expect_no_error(
        transformed <- transform_pvalues_gsd(
            pvals,
            spending = gsDesign::sfLDOF,
            grid_size = 64L
        )
    )

    expect_s3_class(transformed, "multigrain_pvals_gsd")
    expect_identical(dim(transformed$pvals), c(200L, 2L, 3L))
    # "no data" becomes "cannot reject"
    expect_identical(transformed$pvals[, 2L, 1L], rep(1, 200L))
})

test_that("invalid inputs are rejected", {
    corr <- matrix(c(1, 0.5, 0.5, 1), nrow = 2L)

    expect_error(
        simulate_pvalues_gsd(
            c(0.8, 0.9),
            corr_matrix = corr,
            info_frac = c(0.7, 0.5),
            nsim = 10
        ),
        "non-decreasing"
    )
    expect_error(
        simulate_pvalues_gsd(
            c(0.8, 0.9),
            corr_matrix = corr,
            info_frac = rbind(c(0.5, 1), c(NA, NA)),
            nsim = 10
        ),
        "Hypothesis 2 has no analysis"
    )
    expect_error(
        simulate_pvalues_gsd(
            c(0.8, 0.9),
            corr_matrix = corr,
            info_frac = rbind(c(0.5, 1), c(-0.5, 1)),
            nsim = 10
        ),
        "must be positive"
    )
    expect_error(
        simulate_pvalues_gsd(
            c(0.8, 0.9),
            corr_matrix = diag(3),
            info_frac = c(0.5, 1),
            nsim = 10
        ),
        "must be a 2 by 2 matrix"
    )
    expect_error(
        simulate_pvalues_gsd(
            c(0.8, 0.9, 0.9),
            corr_matrix = corr,
            info_frac = c(0.5, 1),
            nsim = 10
        ),
        "must be a 3 by 3 matrix"
    )
    expect_error(
        simulate_pvalues_gsd(
            c(0.8, 0.9),
            corr_matrix = corr,
            info_frac = rbind(c(0.5, 1), c(0.5, 1), c(0.5, 1)),
            nsim = 10
        ),
        "one row per hypothesis"
    )
    expect_error(
        simulate_pvalues_gsd(c(0.8, 0.9), corr_matrix = corr, nsim = 10),
        "is absent but must be supplied"
    )
})

test_that("K = 1 at full information reproduces simulate_pvalues()", {
    power <- c(0.8, 0.85, 0.9)
    corr <- matrix(
        c(
            1, 0.5, 0.5,
            0.5, 1, 0.5,
            0.5, 0.5, 1
        ),
        nrow = 3L
    )

    gsd <- withr::with_seed(
        11,
        simulate_pvalues_gsd(
            power,
            corr_matrix = corr,
            info_frac = 1,
            nsim = 1000
        )
    )
    fixed <- withr::with_seed(
        11,
        simulate_pvalues(power, corr_matrix = corr, nsim = 1000)
    )

    # `identical()` holds, not merely `all.equal()`: with a single analysis at
    # full information the mean vector is `ncp * sqrt(1)` and the covariance is
    # `corr_matrix * 1`, both bit-for-bit the arguments `simulate_pvalues()`
    # passes, so the two `rmvnorm()` calls are the same call
    expect_identical(gsd[,, 1L], fixed)
})
