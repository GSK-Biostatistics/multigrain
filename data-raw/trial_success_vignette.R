## Code to precompute optimised graphs for vignettes/articles/trial_success.qmd
## Run from the package root: Rscript data-raw/trial_success_vignette.R

library(multigrain)

alpha <- 0.025
cores <- parallel::detectCores()

ctrl <- multigrain_control() |>
    control_nsim_global(1e5) |>
    control_global(run = 1000)

# -------------------------------------------------------------------------
# Example 1: FIBRONEER-ILD (step-function for success on at least one dose)
# -------------------------------------------------------------------------

n_per_arm <- 360
delta <- c(0.250, 0.150, 0.230, 0.145)
ncp_dose <- delta * sqrt(n_per_arm)
power_dose <- stats::pnorm(ncp_dose - stats::qnorm(1 - alpha))

corr_dose <- matrix(
    c(
        1.0, 0.4, 0.5, 0.2,
        0.4, 1.0, 0.2, 0.5,
        0.5, 0.2, 1.0, 0.4,
        0.2, 0.5, 0.4, 1.0
    ),
    nrow = 4,
    byrow = TRUE
)

gc_dose <- graph_constraint(
    hyp_constraint = c(NA, 0, NA, 0),
    trans_constraint = matrix(
        c(
            0,  NA, NA, NA,
            0,  0,  NA, NA,
            NA, NA, 0,  NA,
            NA, NA, 0,  0
        ),
        nrow = 4,
        byrow = TRUE
    )
)

win <- trial_success((r1 && r2) || (r3 && r4))

set.seed(2026)
pvals_dose <- simulate_pvalues(
    power_nominal = power_dose,
    alpha = alpha,
    corr_matrix = corr_dose,
    nsim = 1e6
)

graph_dose <- graph_optimise(
    pvals = pvals_dose,
    graph_constraint = gc_dose,
    trial_success = win,
    alpha = alpha,
    control = ctrl,
    num_threads = cores
)

saveRDS(
    graph_dose,
    file = "vignettes/articles/data/trial-success-dose-graph.rds"
)

# -------------------------------------------------------------------------
# Example 2: ECZTRA 3 (hurdle with relative secondary values)
# -------------------------------------------------------------------------

scores <- c(H3 = 85, H4 = 50, H5 = 15)
s_all <- 180

kappa <- (s_all - 100) / 100
v_base <- 1 / (1 + kappa)
v <- kappa / (1 + kappa) * scores / sum(scores)

psi_ecztra <- sprintf(
    "r1 * r2 * (%.3f + %.3f * r3 + %.3f * r4 + %.3f * r5)",
    v_base, v[["H3"]], v[["H4"]], v[["H5"]]
)
gain_ecztra <- trial_success(!!psi_ecztra)

power_ecztra <- c(0.90, 0.90, 0.80, 0.70, 0.85)

corr_ecztra <- matrix(
    c(
        1.0, 0.8, 0.5, 0.5, 0.8,
        0.8, 1.0, 0.5, 0.5, 0.8,
        0.5, 0.5, 1.0, 0.5, 0.5,
        0.5, 0.5, 0.5, 1.0, 0.5,
        0.8, 0.8, 0.5, 0.5, 1.0
    ),
    nrow = 5,
    byrow = TRUE
)

gc_ecztra <- graph_constraint(
    hyp_constraint = c(1, 0, 0, 0, 0),
    trans_constraint = matrix(
        c(
            0, 1, 0,  0,  0,
            0, 0, NA, NA, NA,
            0, 0, 0,  NA, NA,
            0, 0, NA, 0,  NA,
            0, 0, NA, NA, 0
        ),
        nrow = 5,
        byrow = TRUE
    )
)

set.seed(2027)
pvals_ecztra <- simulate_pvalues(
    power_nominal = power_ecztra,
    alpha = alpha,
    corr_matrix = corr_ecztra,
    nsim = 1e6
)

graph_ecztra <- graph_optimise(
    pvals = pvals_ecztra,
    graph_constraint = gc_ecztra,
    trial_success = gain_ecztra,
    alpha = alpha,
    control = ctrl,
    num_threads = cores
)

saveRDS(
    graph_ecztra,
    file = "vignettes/articles/data/trial-success-ecztra-graph.rds"
)
