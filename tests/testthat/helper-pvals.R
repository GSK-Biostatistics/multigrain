pvals_fixture <- withr::with_seed(5, {
    ncp_vec <- rep(stats::qnorm(1 - 0.025), 6) -
        stats::qnorm(1 - c(0.9, 0.8, 0.6, 0.85, 0.85, 0.5))
    corr <- matrix(0.2, nrow = 6, ncol = 6)
    diag(corr) <- 1
    sims <- mvtnorm::rmvnorm(5e4, mean = ncp_vec, sigma = corr)
    stats::pnorm(sims, lower.tail = FALSE)
})

# High-power p-values (independent, mean z = 3), so most trials reject
# several hypotheses
pvals_high_power <- withr::with_seed(11, {
    z <- matrix(stats::rnorm(2000 * 6, mean = 3), ncol = 6)
    stats::pnorm(z, lower.tail = FALSE)
})
