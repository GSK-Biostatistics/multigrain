## code to prepare `example_graph_optimal` dataset goes here
ctrl <- multigrain_control() |>
    control_global(run = 7) |>
    control_nsim_local(2e4)

set.seed(10)

power_vec <- c(0.9, 0.8, 0.6, 0.85, 0.85, 0.5)
ncp_vec <- rep(stats::qnorm(1 - 0.025), 6) - stats::qnorm(1 - power_vec)

corr_test <- matrix(
        0.2,
        nrow = 6,
        ncol = 6,
        byrow = TRUE
    )
diag(corr_test) <- 1

pvals <- withr::with_seed(5, {
    sims <- mvtnorm::rmvnorm(2^20, mean = ncp_vec, sigma = corr_test)
    pvals <- stats::pnorm(sims, lower.tail = FALSE)
})

pvals_4m <- pvals[, 1:4]

m <- 4
no_constr <- graph_constraint_free(m)

conjunctive_4m_power <- trial_success(
    r1 && r2 && r3 && r4,
    verbose = FALSE
)

alpha <- 0.025

graph_optimal_example <- graph_optimise(
    pvals = pvals_4m,
    alpha = alpha,
    graph_constraint = no_constr,
    trial_success = conjunctive_4m_power,
    control = ctrl,
    num_threads = cran_cores(),
    global_search = TRUE
)

usethis::use_data(
    graph_optimal_example,
    overwrite = TRUE
)
