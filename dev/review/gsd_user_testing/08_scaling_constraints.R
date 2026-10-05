# 08_scaling_constraints.R
# Beyond two analyses and three hypotheses, graph constraints, and
# reproducibility of the search.
#
#   D. Three analyses, two hypotheses, a three-level discount.
#   E. A constrained graph on the three-hypothesis design of 02.
#   F. Four hypotheses, two analyses: how run time grows.
#   G. Same seed twice, and one thread against two.
#
# Run from the repository root:
#   Rscript dev/review/gsd_user_testing/08_scaling_constraints.R

suppressPackageStartupMessages({
    library(multigrain)
    library(gsDesign)
})

here <- "dev/review/gsd_user_testing"
if (!dir.exists(here)) {
    here <- "."
}
n_threads <- as.integer(Sys.getenv("GSD_TEST_THREADS", "2"))
alpha <- 0.025

kernel_times <- function(pv, w, G) {
    x <- pv$pvals
    dim(x) <- c(pv$nsim, pv$m * pv$K)
    multigrain:::graph_shortcut_gsd_parallel(
        x,
        pv$alpha,
        w,
        G,
        pv$K,
        n_threads,
        1000L
    )$time
}
ldof_plain <- function(a, t) gsDesign::sfLDOF(a, t)$spend
out <- list()

# =============================================================================
# D. Three analyses
# =============================================================================
cat(
    "\n==================== D. Three analyses, two hypotheses ====================\n"
)
corr2 <- matrix(c(1, 0.5, 0.5, 1), 2)
infoD <- c(1 / 3, 2 / 3, 1)
set.seed(20261011)
tD <- system.time({
    rawD <- simulate_pvalues_gsd(
        c(0.80, 0.90),
        alpha = alpha,
        corr_matrix = corr2,
        info_frac = infoD,
        nsim = 1e5
    )
    pvD <- transform_pvalues_gsd(rawD, spending = sfLDOF, alpha = alpha)
})[["elapsed"]]
summary(pvD)
gainD <- trial_success_gsd(0.6 * d(t1) + 0.4 * d(t2), d = c(1, 0.9, 0.8))
set.seed(1)
elD <- system.time(
    resD <- graph_optimise_gsd(
        pvD,
        graph_constraint_free(2),
        gainD,
        num_threads = n_threads,
        verbose = "silent"
    )
)[["elapsed"]]
G2 <- rbind(c(0, 1), c(1, 0))
w_grid <- seq(0, 1, by = 0.01)
profD <- vapply(
    w_grid,
    function(w1) gainD$func(kernel_times(pvD, c(w1, 1 - w1), G2)),
    0
)
cat(sprintf("Simulate and transform %.1f s; optimise %.0f s\n", tD, elD))
cat(sprintf(
    "Optimiser: w1 = %.4f, gain %.5f. Grid: w1 = %.2f, gain %.5f, plateau [%.2f, %.2f]\n",
    resD$hyp_weight[1],
    resD$power$trial_success,
    w_grid[which.max(profD)],
    max(profD),
    min(w_grid[profD > max(profD) - 5e-4]),
    max(w_grid[profD > max(profD) - 5e-4])
))
print(resD$power[c("local_power_by_analysis", "time_distribution")], digits = 4)

# graphicalMCP on a subsample (its K = 3 boundaries use a randomised integral)
nD <- 0L
dD <- 0L
rD <- 0L
if (requireNamespace("graphicalMCP", quietly = TRUE)) {
    wD <- unname(resD$hyp_weight)
    GD <- unname(resD$trans_matrix)
    tauD <- kernel_times(pvD, wD, GD)
    graph <- graphicalMCP::graph_create(wD, GD)
    set.seed(5)
    idx <- sample.int(1e5, 200)
    for (s in idx) {
        p_s <- matrix(rawD[s, , ], 2, 3, dimnames = list(names(graph$hypotheses), NULL))
        o <- suppressWarnings(graphicalMCP::graph_test_shortcut_gsd(
            graph = graph,
            p = p_s,
            alpha = alpha,
            info_frac = infoD,
            spending_fn = ldof_plain,
            look_back = FALSE
        ))
        rej <- unname(o$outputs$rejected)
        nD <- nD + sum(rej != (tauD[s, ] > 0))
        dD <- dD +
            sum(as.integer(unname(o$outputs$decision_at))[rej] != tauD[s, rej])
        rD <- rD + sum(rej)
    }
    cat(sprintf(
        "graphicalMCP on 200 trials: %d disagreements in 400 reject or retain decisions; %d in the decision times of %d rejections\n",
        nD,
        dD,
        rD
    ))
}
out$three_looks <- list(
    hyp_weight = resD$hyp_weight,
    gain = resD$power$trial_success,
    seconds = elD,
    grid_w1 = w_grid[which.max(profD)],
    grid_gain = max(profD),
    oracle_disagreements = c(nD, dD, rD)
)

# =============================================================================
# E. Constrained graph
# =============================================================================
cat("\n==================== E. Constrained graph ====================\n")
sav <- readRDS(file.path(here, "02_three_hyp_two_stage_results.rds"))
ds <- sav$design
set.seed(ds$seed)
raw3 <- simulate_pvalues_gsd(
    ds$power_nom,
    alpha = alpha,
    corr_matrix = ds$corr,
    info_frac = ds$info_frac,
    nsim = ds$nsim
)
pv3 <- transform_pvalues_gsd(raw3, spending = sfLDOF, alpha = alpha)
gain3 <- trial_success_gsd(
    !!ds$value[1] *
        d(t1) +
        r1 * (!!ds$value[2] * d(t2) + !!ds$value[3] * d(t3)),
    d = c(1, ds$delta),
    verbose = "silent"
)
# H1 is a gatekeeper with all the initial weight; H2 and H3 recycle to each
# other; only the split of H1's level between H2 and H3 is free.
con <- graph_constraint(
    hyp_constraint = c(1, 0, 0),
    trans_constraint = rbind(c(0, NA, NA), c(0, 0, 1), c(0, 1, 0))
)
print(con)
set.seed(1)
elE <- system.time(
    resE <- graph_optimise_gsd(
        pv3,
        con,
        gain3,
        num_threads = n_threads,
        verbose = "silent"
    )
)[["elapsed"]]
cat(sprintf(
    "Optimised in %.0f s: gain %.5f (free optimum in 02: %.5f)\n",
    elE,
    resE$power$trial_success,
    sav$optimised$discounted$power$trial_success
))
cat("w =", round(resE$hyp_weight, 4), "\n")
print(round(resE$trans_matrix, 4))
fixed_ok <- isTRUE(all.equal(unname(resE$hyp_weight), c(1, 0, 0))) &&
    isTRUE(all.equal(
        unname(resE$trans_matrix[2:3, ]),
        rbind(c(0, 0, 1), c(0, 1, 0))
    )) &&
    isTRUE(all.equal(sum(resE$trans_matrix[1, ]), 1))
cat("Fixed entries respected and H1's row sums to 1:", fixed_ok, "\n")
out$constrained <- list(
    hyp_weight = resE$hyp_weight,
    trans_matrix = resE$trans_matrix,
    gain = resE$power$trial_success,
    seconds = elE,
    fixed_ok = fixed_ok
)

# =============================================================================
# F. Four hypotheses
# =============================================================================
cat(
    "\n==================== F. Four hypotheses, two analyses ====================\n"
)
corr4 <- matrix(0.5, 4, 4)
diag(corr4) <- 1
set.seed(20261012)
raw4 <- simulate_pvalues_gsd(
    c(0.90, 0.85, 0.80, 0.75),
    alpha = alpha,
    corr_matrix = corr4,
    info_frac = c(0.5, 1),
    nsim = 1e5
)
pv4 <- transform_pvalues_gsd(raw4, spending = sfLDOF, alpha = alpha)
gain4 <- trial_success_gsd(
    0.4 * d(t1) + 0.3 * d(t2) + 0.2 * d(t3) + 0.1 * d(t4),
    d = c(1, 0.8),
    verbose = "silent"
)
set.seed(1)
elF <- system.time(
    resF <- graph_optimise_gsd(
        pv4,
        graph_constraint_free(4),
        gain4,
        num_threads = n_threads,
        verbose = "silent"
    )
)[["elapsed"]]
G4_holm <- matrix(1 / 3, 4, 4)
diag(G4_holm) <- 0
G4_seq <- rbind(c(0, 1, 0, 0), c(0, 0, 1, 0), c(0, 0, 0, 1), c(1, 0, 0, 0))
cmpF <- c(
    optimised = resF$power$trial_success,
    holm = gain4$func(kernel_times(pv4, rep(0.25, 4), G4_holm)),
    fixed_sequence = gain4$func(kernel_times(pv4, c(1, 0, 0, 0), G4_seq))
)
cat(sprintf("Optimised in %.0f s\n", elF))
print(round(cmpF, 5))
cat("w =", round(resF$hyp_weight, 3), "\n")
print(round(resF$trans_matrix, 3))
out$four_hyp <- list(
    hyp_weight = resF$hyp_weight,
    trans_matrix = resF$trans_matrix,
    gains = cmpF,
    seconds = elF
)

# =============================================================================
# G. Reproducibility
# =============================================================================
cat("\n==================== G. Reproducibility ====================\n")
set.seed(20261013)
rawG <- simulate_pvalues_gsd(
    c(0.80, 0.90),
    alpha = alpha,
    corr_matrix = corr2,
    info_frac = c(0.5, 1),
    nsim = 2e4
)
pvG <- transform_pvalues_gsd(rawG, spending = sfLDOF, alpha = alpha)
gainG <- trial_success_gsd(
    0.6 * d(t1) + 0.4 * d(t2),
    d = c(1, 0.8),
    verbose = "silent"
)
run <- function(seed, threads) {
    set.seed(seed)
    r <- graph_optimise_gsd(
        pvG,
        graph_constraint_free(2),
        gainG,
        num_threads = threads,
        verbose = "silent"
    )
    list(w = r$hyp_weight, G = r$trans_matrix, gain = r$power$trial_success)
}
a <- run(1, 2L)
b <- run(1, 2L)
c1 <- run(1, 1L)
cat(
    "Same seed, two threads, run twice: identical result:",
    identical(a, b),
    "\n"
)
cat(
    "Same seed, one thread against two: identical result:",
    identical(a, c1),
    "\n"
)
cat(sprintf(
    "w1 = %.6f (two threads), %.6f (repeat), %.6f (one thread)\n",
    a$w[1],
    b$w[1],
    c1$w[1]
))
out$reproducibility <- list(
    same_seed = identical(a, b),
    threads = identical(a, c1)
)

timing <- c(three_looks = elD, constrained = elE, four_hyp = elF)
cat("\nOptimisation seconds (", n_threads, " threads):\n", sep = "")
print(round(timing, 1))
out$timing <- timing
saveRDS(out, file.path(here, "08_scaling_constraints_results.rds"))
