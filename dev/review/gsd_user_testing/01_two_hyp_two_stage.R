# 01_two_hyp_two_stage.R
# User test of the group sequential extension (branch gsd-build).
# A 2-stage, 2-endpoint trial: simulate, transform, define a time-discounted
# gain with trial_success_gsd(), optimise the graph, and compare with simple
# graphs and with a grid over the single free weight.
#
# Run from the repository root with an *installed* build of the branch:
#   Rscript dev/review/gsd_user_testing/01_two_hyp_two_stage.R

suppressPackageStartupMessages({
    library(multigrain)
    library(gsDesign)
})

out_dir <- "dev/review/gsd_user_testing"
if (!dir.exists(out_dir)) out_dir <- "."
n_threads <- as.integer(Sys.getenv("GSD_TEST_THREADS", "2"))

# ---- design ---------------------------------------------------------------
# H1 is the more valuable endpoint but the less well powered one; H2 is worth
# less and is better powered. One interim at 50% information, then the final.
alpha     <- 0.025
power_nom <- c(H1 = 0.80, H2 = 0.90)      # final-analysis power at full alpha
corr      <- matrix(c(1, 0.5, 0.5, 1), 2)
info_frac <- c(0.5, 1)
nsim      <- 1e5
value     <- c(0.6, 0.4)                  # value of each claim
delta     <- 0.8                          # a final-analysis claim keeps 80%

# ---- simulate and transform ----------------------------------------------
set.seed(20261005)
t_sim <- system.time(
    raw <- simulate_pvalues_gsd(
        power_nominal = power_nom, alpha = alpha,
        corr_matrix = corr, info_frac = info_frac, nsim = nsim
    )
)
t_tr <- system.time(
    pvals <- transform_pvalues_gsd(raw, spending = sfLDOF, alpha = alpha)
)
summary(pvals)

# ---- gains ------------------------------------------------------------------
# The discounted gain is the object of the exercise; the undiscounted one is
# optimised as well so that the effect of the discount on the graph is visible.
gain_disc <- trial_success_gsd(
    !!value[1] * d(t1) + !!value[2] * d(t2),
    d = c(1, delta)
)
gain_flat <- trial_success_gsd(
    !!value[1] * r1 + !!value[2] * r2,
    K = 2
)
print(gain_disc)

# ---- optimise -------------------------------------------------------------
optimise <- function(gain, seed = 1) {
    set.seed(seed)
    elapsed <- system.time(
        res <- graph_optimise_gsd(
            pvals = pvals,
            graph_constraint = graph_constraint_free(2),
            trial_success = gain,
            num_threads = n_threads
        )
    )[["elapsed"]]
    list(res = res, elapsed = elapsed)
}
opt_disc <- optimise(gain_disc)
opt_flat <- optimise(gain_flat)
res <- opt_disc$res

print(res)
summary(res)
cat("\nPower fields that print() and summary() do not show:\n")
print(res$power[c("local_power_by_analysis", "mean_decision_look",
                  "time_distribution")])

# ---- comparators on the same draws -----------------------------------------
G_full <- rbind(c(0, 1), c(1, 0))
eval_graph <- function(w, G = G_full) {
    calc_power_pvals_gsd(
        pvals, hyp_weight = w, trans_matrix = G,
        custom_power = list(disc = gain_disc, flat = gain_flat)
    )
}
graphs <- list(
    `optimised (discounted gain)`   = list(w = res$hyp_weight, G = res$trans_matrix),
    `optimised (undiscounted gain)` = list(w = opt_flat$res$hyp_weight,
                                           G = opt_flat$res$trans_matrix),
    `Holm`                          = list(w = c(0.5, 0.5), G = G_full),
    `fixed sequence H1 -> H2`       = list(w = c(1, 0), G = G_full),
    `fixed sequence H2 -> H1`       = list(w = c(0, 1), G = G_full)
)
comp <- do.call(rbind, lapply(names(graphs), function(nm) {
    g <- graphs[[nm]]
    p <- eval_graph(g$w, g$G)
    data.frame(
        graph = nm, w1 = g$w[1], w2 = g$w[2],
        gain_discounted = p$disc, gain_undiscounted = p$flat,
        power_H1 = p$local_power[1], power_H2 = p$local_power[2],
        interim_H1 = p$local_power_by_analysis[1, 1],
        interim_H2 = p$local_power_by_analysis[2, 1]
    )
}))
cat("\nComparators (same simulated trials):\n")
print(comp, digits = 4, row.names = FALSE)

# ---- gain profile over w1 (two hypotheses: one free weight) ----------------
w_grid  <- seq(0, 1, by = 0.01)
prof    <- lapply(w_grid, function(w1) eval_graph(c(w1, 1 - w1)))
profile <- data.frame(
    w1 = w_grid,
    gain_discounted   = vapply(prof, function(p) p$disc, 0),
    gain_undiscounted = vapply(prof, function(p) p$flat, 0)
)
report_profile <- function(col, opt, label) {
    g <- profile[[col]]
    i <- which.max(g)
    plateau <- range(profile$w1[g > g[i] - 5e-4])
    cat(sprintf(
        "%s: grid max %.5f at w1 = %.2f (within 5e-4 for w1 in [%.2f, %.2f]); optimiser %.5f at w1 = %.4f\n",
        label, g[i], profile$w1[i], plateau[1], plateau[2],
        opt$res$power$trial_success, opt$res$hyp_weight[1]
    ))
}
cat("\nGrid over w1 with full recycling:\n")
report_profile("gain_discounted", opt_disc, "discounted  ")
report_profile("gain_undiscounted", opt_flat, "undiscounted")

timing <- c(simulate = t_sim[["elapsed"]], transform = t_tr[["elapsed"]],
            optimise_discounted = opt_disc$elapsed,
            optimise_undiscounted = opt_flat$elapsed)
cat("\nElapsed seconds (", n_threads, " threads):\n", sep = "")
print(round(timing, 1))

saveRDS(
    list(
        design = list(alpha = alpha, power_nom = power_nom, corr = corr,
                      info_frac = info_frac, nsim = nsim, value = value,
                      delta = delta, spending = "sfLDOF", seed = 20261005),
        optimised = list(
            discounted = list(hyp_weight = res$hyp_weight,
                              trans_matrix = res$trans_matrix,
                              power = res$power, solution = res$solution),
            undiscounted = list(hyp_weight = opt_flat$res$hyp_weight,
                                trans_matrix = opt_flat$res$trans_matrix,
                                power = opt_flat$res$power,
                                solution = opt_flat$res$solution)
        ),
        comparators = comp, profile = profile, timing = timing,
        threads = n_threads,
        session = c(R = R.version.string,
                    multigrain = as.character(utils::packageVersion("multigrain")),
                    gsDesign = as.character(utils::packageVersion("gsDesign")))
    ),
    file.path(out_dir, "01_two_hyp_two_stage_results.rds")
)
