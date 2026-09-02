m <- 4
alpha <- 0.025
power_vector <- c(0.93, 0.91, 0.90, 0.85)

corr_mat <- matrix(
    0.2,
    nrow = m,
    ncol = m,
    byrow = TRUE
)
diag(corr_mat) <- 1

sims <- mvtnorm::rmvnorm(
    5e4,
    sigma = as.matrix(corr_mat),
    mean = calc_ncp(power_vector)
)
pvals_4m <- stats::pnorm(sims, lower.tail = FALSE)

no_constr <- graph_constraint_free(m)

conjunctive_4m_power <- trial_success(r1 && r2 && r3 && r4, verbose = FALSE)

test_that("deprecation message for optimise_graph and optimise_graph", {
    ctrl <- multigrain_control() |>
        control_global(run = 7) |>
        control_nsim_global(2e4) |>
        control_nsim_local(2e4)

    set.seed(20)

    expect_snapshot_warning({
        graph_optimise_4m_result <- optimise_graph(
            pvals = pvals_4m,
            alpha = alpha,
            graph_constraint = no_constr,
            trial_success = conjunctive_4m_power,
            control = ctrl,
            verbose = FALSE
        )
    })

    expect_s3_class(graph_optimise_4m_result, "multigrain_graph_optimal")

    expect_snapshot_warning({
        graph_optimise_4m_result2 <- optimize_graph(
            pvals = pvals_4m,
            alpha = alpha,
            graph_constraint = no_constr,
            trial_success = conjunctive_4m_power,
            control = ctrl,
            verbose = FALSE
        )
    })

    expect_s3_class(graph_optimise_4m_result2, "multigrain_graph_optimal")
})
