# Graph average power ----------------------------------------------------

# Storing power vector and correlation of test statistics
power_vector <- c(0.97, 0.91, 0.86, 0.83)
corr_matrix <- matrix(
    c(
        1, 0.5, 0.8, 0.4,
        0.5, 1, 0.4, 0.8,
        0.8, 0.4, 1, 0.5,
        0.4, 0.8, 0.5, 1
    ),
    nrow = 4
)

# Creating power objective functions as targets for optimisation
average_power <- trial_success(0.25 * (r1 + r2 + r3 + r4))
custom_power <- trial_success(0.25 * (2 * (r1 && r2) + r1 * r3 + r2 * r4))

# Creating a `multigrain_graph_constraint` object to enforce logical constraints
# on order of rejections

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

print(my_constraint)

pvals <- simulate_pvalues(
    power_nominal = power_vector,
    alpha = 0.025
)

ctrl <- multigrain_control()

# Find optimal graph based on average power objective function
graph_average_power <- graph_optimise(
    pvals = pvals,
    graph_constraint = my_constraint,
    trial_success = average_power,
    alpha = 0.025, # one-sided
    num_threads = cran_cores(),
    control = ctrl,
    trace = TRUE
)

saveRDS(
    graph_average_power,
    file = file.path(
        "tests",
        "testthat",
        "data",
        "graph_average_power.rds"
    )
)


# Graph custom power -----------------------------------------------------

# Find optimal graph based on custom power objective function (defined earlier)
graph_custom_power <- graph_optimise(
    pvals = pvals,
    graph_constraint = my_constraint,
    trial_success = custom_power,
    alpha = 0.025, # one-sided
    num_threads = cran_cores(),
    control = ctrl,
    trace = TRUE
)

saveRDS(
    graph_custom_power,
    file = file.path(
        "tests",
        "testthat",
        "data",
        "graph_custom_power.rds"
    )
)
