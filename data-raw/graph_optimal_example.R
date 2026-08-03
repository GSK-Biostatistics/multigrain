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

hyp_weight_constr <- c(NA, NA, 0, 0)
trans_matrix_constr <- matrix(
    c(
        0, NA, NA, 0,
        NA, 0, 0, NA,
        0, 1, 0, 0,
        1, 0, 0, 0
    ),
    nrow = 4,
    byrow = TRUE
)

my_constraint <- graph_constraint(
    hyp_constraint = hyp_weight_constr,
    trans_constraint = trans_matrix_constr
)

custom_power <- trial_success(0.25 * (2 * (r1 && r2) + r1 * r3 + r2 * r4))

graph_optimal_example <- graph_optimise(
    pvals = pvals_4m,
    graph_constraint = my_constraint,
    trial_success = custom_power,
    control = ctrl,
    num_threads = cran_cores(),
    global_search = TRUE,
    verbose = "detail"
)

usethis::use_data(
    graph_optimal_example,
    overwrite = TRUE
)
