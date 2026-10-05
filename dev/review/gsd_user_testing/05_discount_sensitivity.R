# 05_discount_sensitivity.R
# How the discount for a late claim and the spending function move the optimal
# split of alpha between two hypotheses. Grid over w1 with full recycling; one
# kernel run per grid point gives every discount, because an additive gain is
# linear in the decision-time distribution.
#
# Run from the repository root:
#   Rscript dev/review/gsd_user_testing/05_discount_sensitivity.R

suppressPackageStartupMessages({
    library(multigrain)
    library(gsDesign)
})

here <- "dev/review/gsd_user_testing"
if (!dir.exists(here)) here <- "."
n_threads <- as.integer(Sys.getenv("GSD_TEST_THREADS", "2"))

corr   <- matrix(c(1, 0.5, 0.5, 1), 2)
G_full <- rbind(c(0, 1), c(1, 0))
w_grid <- seq(0, 1, by = 0.02)
deltas <- c(1, 0.8, 0.5, 0.2)
pocock <- function(a, t) gsDesign::sfLDPocock(a, t)

profile_case <- function(power_nom, spending, value, nsim = 2e5, seed = 20261008) {
    set.seed(seed)
    raw <- simulate_pvalues_gsd(power_nom, corr_matrix = corr,
                                info_frac = c(0.5, 1), nsim = nsim)
    pv  <- transform_pvalues_gsd(raw, spending = spending)
    td  <- lapply(w_grid, function(w1) {
        calc_power_pvals_gsd(pv, c(w1, 1 - w1), G_full)$time_distribution
    })
    sapply(deltas, function(dl) {
        vapply(td, function(x) sum(value * (x[, 2] + dl * x[, 3])), 0)
    })
}

cases <- list(
    `power 0.80/0.90, value 0.6/0.4, LDOF/LDOF (design of 01)` =
        list(c(0.80, 0.90), sfLDOF, c(0.6, 0.4)),
    `power 0.80/0.90, value 0.6/0.4, Pocock/Pocock` =
        list(c(0.80, 0.90), pocock, c(0.6, 0.4)),
    `power 0.85/0.85, value 0.5/0.5, LDOF/LDOF (symmetric)` =
        list(c(0.85, 0.85), sfLDOF, c(0.5, 0.5)),
    `power 0.85/0.85, value 0.5/0.5, Pocock/Pocock (symmetric)` =
        list(c(0.85, 0.85), pocock, c(0.5, 0.5))
)

out <- list()
for (nm in names(cases)) {
    cs <- cases[[nm]]
    g  <- profile_case(cs[[1]], cs[[2]], cs[[3]])
    colnames(g) <- paste0("delta=", deltas)
    tab <- apply(g, 2, function(col) {
        i <- which.max(col)
        c(w1_opt = w_grid[i], gain_opt = col[i],
          gain_holm = col[w_grid == 0.5],
          gain_w1_0 = col[1], gain_w1_1 = col[length(col)],
          plateau_lo = min(w_grid[col > col[i] - 5e-4]),
          plateau_hi = max(w_grid[col > col[i] - 5e-4]))
    })
    cat("\n== ", nm, "\n", sep = "")
    print(round(tab, 4))
    out[[nm]] <- list(profile = cbind(w1 = w_grid, g), table = tab)
}

# The symmetric LDOF case in full: with a strong discount the profile has a
# dip at the Holm point w1 = 0.5 and two maxima near the ends.
sym <- out[["power 0.85/0.85, value 0.5/0.5, LDOF/LDOF (symmetric)"]]$profile
cat("\nSymmetric LDOF/LDOF profile at selected w1:\n")
print(round(sym[round(w_grid * 100) %in% c(0, 4, 10, 20, 30, 40, 50, 60, 70, 80, 90, 96, 100), ], 4))

# Is the dip at the Holm point more than Monte Carlo noise, and does the search
# leave it? Paired difference of per-trial gains between w1 = 0.10 and 0.50 on
# the same draws, then the optimiser started as usual (the default start graph
# is Bonferroni-Holm) with and without the global search.
set.seed(20261008)
raw_s <- simulate_pvalues_gsd(c(0.85, 0.85), corr_matrix = corr,
                              info_frac = c(0.5, 1), nsim = 2e5)
pv_s  <- transform_pvalues_gsd(raw_s, spending = sfLDOF)
x_s   <- pv_s$pvals
dim(x_s) <- c(pv_s$nsim, 4L)
per_trial <- function(w1, dl = 0.5) {
    tau <- multigrain:::graph_shortcut_gsd(x_s, 0.025, c(w1, 1 - w1), G_full, 2L)$time
    dd  <- c(0, 1, dl)
    0.5 * dd[tau[, 1] + 1L] + 0.5 * dd[tau[, 2] + 1L]
}
dif <- per_trial(0.10) - per_trial(0.50)
cat(sprintf("\nSymmetric LDOF, delta = 0.5: gain(w1 = 0.10) - gain(w1 = 0.50) = %.5f (paired s.e. %.5f, z = %.1f)\n",
            mean(dif), stats::sd(dif) / sqrt(length(dif)), mean(dif) / (stats::sd(dif) / sqrt(length(dif)))))
gain_s <- trial_success_gsd(0.5 * d(t1) + 0.5 * d(t2), d = c(1, 0.5), verbose = "silent")
set.seed(1)
loc <- graph_optimise_gsd(pv_s, graph_constraint_free(2), gain_s,
                          global_search = FALSE, num_threads = n_threads, verbose = "silent")
set.seed(1)
glo <- graph_optimise_gsd(pv_s, graph_constraint_free(2), gain_s,
                          num_threads = n_threads, verbose = "silent")
cat(sprintf("Optimiser from the default start: local only w1 = %.3f (gain %.5f); global then local w1 = %.3f (gain %.5f)\n",
            loc$hyp_weight[1], loc$power$trial_success,
            glo$hyp_weight[1], glo$power$trial_success))
out$symmetric_check <- list(
    dip = mean(dif), dip_se = stats::sd(dif) / sqrt(length(dif)),
    local_only = loc$hyp_weight, global = glo$hyp_weight
)

# Why: the interim boundary under O'Brien-Fleming-type spending is strongly
# convex in the level a hypothesis holds, so splitting alpha costs far more
# than proportionally at the interim. Under Pocock-type spending it is close
# to linear.
lv <- c(0.025, 0.0125, 0.00625)
bd <- function(sf) sapply(lv, function(a) {
    inc <- diff(c(0, sf(a, c(0.5, 1))$spend))
    stats::pnorm(gsBound1(theta = 0, I = c(0.5, 1), a = rep(-20, 2), probhi = inc)$b,
                 lower.tail = FALSE)[1]
})
cat("\nInterim nominal boundary at levels 0.025, 0.0125, 0.00625 (information fraction 0.5):\n")
print(rbind(LDOF = bd(sfLDOF), Pocock = bd(sfLDPocock)), digits = 4)

saveRDS(out, file.path(here, "05_discount_sensitivity_results.rds"))
