test_that("simulate_pvalues", {
    nominal_power <- c(0.8, 0.85, 0.9)
    corr <- matrix(
        c(
            1, 0.5, 0.5,
            0.5, 1, 0.5,
            0.5, 0.5, 1
        ),
        nrow = 3
    )

    num_simulations <- 100

    set.seed(1)

    pvals <- simulate_pvalues(
        nominal_power,
        corr_matrix = corr,
        nsim = num_simulations
    )

    expect_identical(
        dim(pvals),
        c(100L, 3L)
    )

    expect_true(is.matrix(pvals))
})

test_that("simulate_pvalues complains when users pass anything via `...`", {
    nominal_power <- c(0.8, 0.85, 0.9)
    alpha_level <- 0.025

    expect_snapshot(error = TRUE, {
        simulate_pvalues(
            nominal_power,
            alpha_level
        )
    })
})
