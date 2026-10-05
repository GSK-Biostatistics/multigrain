# 07_features.R
# Features of the GSD family that 01 and 02 do not touch.
#
#   A. Manuscript Example 5: PFS complete at the interim, OS at information
#      fractions 0.7 and 1. All 21 (r, delta) cells against the paper's
#      reference, and the optimiser at three cells.
#   B. Look-back on and off on the three-hypothesis design of 02.
#   C. Mixed design: different spending functions and information fractions
#      per hypothesis and one endpoint with no data at the interim.
#
# Run from the repository root:
#   Rscript dev/review/gsd_user_testing/07_features.R

suppressPackageStartupMessages({
    library(multigrain)
    library(gsDesign)
})

here <- "dev/review/gsd_user_testing"
if (!dir.exists(here)) here <- "."
source(file.path(here, "00_reference.R"))
n_threads <- as.integer(Sys.getenv("GSD_TEST_THREADS", "2"))
alpha <- 0.025

kernel_times <- function(pv, w, G) {
    x <- pv$pvals
    dim(x) <- c(pv$nsim, pv$m * pv$K)
    multigrain:::graph_shortcut_gsd_parallel(x, pv$alpha, w, G, pv$K, n_threads, 1000L)$time
}
out <- list()

# =============================================================================
# A. Example 5
# =============================================================================
cat("\n==================== A. Example 5 (PFS and OS) ====================\n")
ref5 <- readRDS("tests/testthat/data/gsd_example5_reference.rds")
inp  <- ref5$inputs
info5 <- rbind(PFS = c(1, 1), OS = c(inp$t_os, 1))
pow5  <- stats::pnorm(c(inp$Delta_P, inp$Delta_O) - stats::qnorm(1 - alpha))
corr5 <- matrix(c(1, inp$rho, inp$rho, 1), 2)

set.seed(20261009)
raw5 <- simulate_pvalues_gsd(pow5, alpha = alpha, corr_matrix = corr5,
                             info_frac = info5, nsim = 1e5)
pv5  <- transform_pvalues_gsd(raw5, spending = sfLDOF, alpha = alpha)
summary(pv5)
G2 <- rbind(c(0, 1), c(1, 0))

# decision-time distribution on a grid over the PFS weight
w_grid <- seq(0, 1, by = 0.005)
td <- lapply(w_grid, function(w1) {
    calc_power_pvals_gsd(pv5, c(w1, 1 - w1), G2)$time_distribution
})
egain <- function(x, r, delta) {
    v <- c(1, r) / (1 + r)
    sum(v * (x[, 2] + delta * x[, 3]))
}
cells <- ref5$w_star
cells$gain_at_w_star <- cells$grid_max <- cells$w1_grid <- NA_real_
for (j in seq_len(nrow(cells))) {
    g <- vapply(td, egain, 0, r = cells$r[j], delta = cells$delta[j])
    cells$w1_grid[j]        <- w_grid[which.max(g)]
    cells$grid_max[j]       <- max(g)
    cells$gain_at_w_star[j] <- g[which.min(abs(w_grid - cells$w1_star[j]))]
}
cells$gap_at_w_star <- cells$grid_max - cells$gain_at_w_star
cells$minus_paper   <- cells$gain_at_w_star - cells$Egain_star
cat("\nAll 21 cells: paper's optimal weight w1_star and expected gain (N = 5e6, classical\nO'Brien-Fleming boundary) against this run (N = 1e5, sfLDOF)\n")
print(cells[order(cells$delta, cells$r), ], digits = 4, row.names = FALSE)
cat(sprintf("largest gain lost by using the paper's weight in place of this run's grid optimum: %.1e\n",
            max(cells$gap_at_w_star)))
cat(sprintf("gain at the paper's weight minus the paper's gain: range %.4f to %.4f\n",
            min(cells$minus_paper), max(cells$minus_paper)))

# The paper's boundary. The offset above grows as delta falls, which points at
# the interim boundary: the paper used the classical O'Brien-Fleming boundary
# (0.0082 at the interim), not sfLDOF (0.0074). The design record notes that
# the classical boundary can be supplied as a spending function through
# gsDesign(). As written there it fails twice inside the transform: on the
# single information fraction of 1 that the spend-at-1 check passes in
# (gsDesign needs at least two analyses), and, once that case is guarded, with
# "f() values at end points not of opposite sign" at some allocated levels
# between about 4e-14 and 3e-11. Those levels fall back to sfLDOF here; no
# graph gives a hypothesis a level that small in this example.
of_spend <- function(a, t) {
    if (length(t) == 1L) return(a * (t >= 1))
    r <- try(
        cumsum(gsDesign::gsDesign(k = length(t), test.type = 1, alpha = a,
                                  sfu = "OF", timing = t)$upper$spend),
        silent = TRUE
    )
    if (inherits(r, "try-error")) gsDesign::sfLDOF(a, t)$spend else r
}
grid_levels <- exp(seq(log(1e-14), log(alpha), length.out = 1024))
of_fails <- vapply(grid_levels, function(a) {
    inherits(try(gsDesign::gsDesign(k = 2, test.type = 1, alpha = a, sfu = "OF",
                                    timing = c(inp$t_os, 1)), silent = TRUE), "try-error")
}, TRUE)
cat(sprintf("\ngsDesign(sfu = \"OF\") fails at %d of the 1024 table levels, from %.1e to %.1e\n",
            sum(of_fails), min(grid_levels[of_fails]), max(grid_levels[of_fails])))
t_of <- system.time(
    pv5_of <- transform_pvalues_gsd(raw5, spending = list(PFS = sfLDOF, OS = of_spend),
                                    alpha = alpha)
)[["elapsed"]]
cat(sprintf("Transform with the classical boundary: %.0f s\n", t_of))
summary(pv5_of)
td_of <- lapply(w_grid, function(w1) {
    calc_power_pvals_gsd(pv5_of, c(w1, 1 - w1), G2)$time_distribution
})
cells$of_gain_at_w_star <- NA_real_
for (j in seq_len(nrow(cells))) {
    g <- vapply(td_of, egain, 0, r = cells$r[j], delta = cells$delta[j])
    cells$of_gain_at_w_star[j] <- g[which.min(abs(w_grid - cells$w1_star[j]))]
}
cells$of_minus_paper <- cells$of_gain_at_w_star - cells$Egain_star
cat("With the paper's boundary, gain at the paper's weight minus the paper's gain:\n")
print(cells[order(cells$delta, cells$r), c("r", "delta", "Egain_star", "of_gain_at_w_star",
                                            "of_minus_paper")], digits = 4, row.names = FALSE)
cat(sprintf("range %.4f to %.4f (one standard error of a gain at N = 1e5 is roughly 0.0005 to 0.001)\n",
            min(cells$of_minus_paper), max(cells$of_minus_paper)))

# The remaining offset has one sign in every cell because the 21 cells share
# one set of draws. Two fresh sets of 1e6 trials with the paper's boundary
# settle whether it is noise.
big <- do.call(rbind, lapply(c(101, 102), function(seed) {
    set.seed(seed)
    raw_b <- simulate_pvalues_gsd(pow5, alpha = alpha, corr_matrix = corr5,
                                  info_frac = info5, nsim = 1e6)
    pv_b  <- transform_pvalues_gsd(raw_b, spending = list(PFS = sfLDOF, OS = of_spend),
                                   alpha = alpha)
    dif <- vapply(seq_len(nrow(cells)), function(j) {
        x <- calc_power_pvals_gsd(pv_b, c(cells$w1_star[j], 1 - cells$w1_star[j]),
                                  G2)$time_distribution
        egain(x, cells$r[j], cells$delta[j]) - cells$Egain_star[j]
    }, 0)
    data.frame(seed = seed, nsim = 1e6, min = min(dif), mean = mean(dif), max = max(dif))
}))
cat("\nWith the paper's boundary at N = 1e6, gain at the paper's weight minus the paper's gain over the 21 cells:\n")
print(big, digits = 3, row.names = FALSE)

# the optimiser at three cells
opt5 <- do.call(rbind, lapply(list(c(1, 1), c(4, 0.75), c(8, 0.5)), function(cd) {
    r <- cd[1]; delta <- cd[2]
    v <- c(1, r) / (1 + r)
    gain <- trial_success_gsd(!!v[1] * d(t1) + !!v[2] * d(t2), d = c(1, delta),
                              verbose = "silent")
    set.seed(1)
    el <- system.time(
        res <- graph_optimise_gsd(pv5, graph_constraint_free(2), gain,
                                  num_threads = n_threads, verbose = "silent")
    )[["elapsed"]]
    j <- which(cells$r == r & cells$delta == delta)
    data.frame(r = r, delta = delta, w1_optimiser = res$hyp_weight[1],
               gain_optimiser = res$power$trial_success,
               w1_grid = cells$w1_grid[j], gain_grid = cells$grid_max[j],
               w1_paper = cells$w1_star[j], gain_paper = cells$Egain_star[j],
               seconds = el)
}))
cat("\nOptimiser at three cells:\n")
print(opt5, digits = 4, row.names = FALSE)

# the "PFS rejected late" mechanism, and the reference on all trials
w5 <- c(0.2, 0.8)
tau5 <- kernel_times(pv5, w5, G2)
ref5t <- ref_test(unclass(raw5), w5, G2, t1 = c(1, inp$t_os), sf = ref_sf_ldof, alpha = alpha)
cat(sprintf("\nGraph w = (0.2, 0.8): reference against transform + kernel, %d of %d decision times differ\n",
            sum(tau5 != ref5t), length(tau5)))
cat(sprintf("PFS declared at the final analysis although its data were complete at the interim: %.4f of trials\n",
            mean(tau5[, 1] == 2)))
out$example5 <- list(cells = cells, optimiser = opt5, n_diff = sum(tau5 != ref5t),
                     paper_boundary_1e6 = big, of_fails = sum(of_fails),
                     pfs_late = mean(tau5[, 1] == 2))

# =============================================================================
# B. Look-back
# =============================================================================
cat("\n==================== B. Look-back ====================\n")
sav <- readRDS(file.path(here, "02_three_hyp_two_stage_results.rds"))
ds  <- sav$design
set.seed(ds$seed)
raw3 <- simulate_pvalues_gsd(ds$power_nom, alpha = alpha, corr_matrix = ds$corr,
                             info_frac = ds$info_frac, nsim = ds$nsim)
pv_off <- transform_pvalues_gsd(raw3, spending = sfLDOF, alpha = alpha)
pv_on  <- transform_pvalues_gsd(raw3, spending = sfLDOF, alpha = alpha, look_back = TRUE)
pv_mix <- transform_pvalues_gsd(raw3, spending = sfLDOF, alpha = alpha,
                                look_back = c(FALSE, TRUE, TRUE))
gain3 <- trial_success_gsd(
    !!ds$value[1] * d(t1) + r1 * (!!ds$value[2] * d(t2) + !!ds$value[3] * d(t3)),
    d = c(1, ds$delta), verbose = "silent"
)
G_holm <- matrix(0.5, 3, 3); diag(G_holm) <- 0
graphs <- list(
    `optimised graph of 02` = list(w = unname(sav$optimised$discounted$hyp_weight),
                                   G = unname(sav$optimised$discounted$trans_matrix)),
    `Holm` = list(w = rep(1, 3) / 3, G = G_holm)
)
lb <- do.call(rbind, lapply(names(graphs), function(nm) {
    g <- graphs[[nm]]
    t_off <- kernel_times(pv_off, g$w, g$G)
    t_on  <- kernel_times(pv_on, g$w, g$G)
    t_mix <- kernel_times(pv_mix, g$w, g$G)
    both  <- t_off > 0
    data.frame(
        graph = nm,
        gain_off = gain3$func(t_off), gain_mixed = gain3$func(t_mix),
        gain_on = gain3$func(t_on),
        trials_changed = mean(rowSums(t_off != t_on) > 0),
        rejections_lost = sum(t_off > 0 & t_on == 0),
        rejections_gained = sum(t_off == 0 & t_on > 0),
        decided_later = sum(t_on[both] > t_off[both]),
        decided_earlier = sum(t_on[both] < t_off[both])
    )
}))
print(lb, digits = 5, row.names = FALSE)
cat("Look-back should never lose a rejection or delay one: rejections_lost and decided_later must be 0.\n")

# FWER with look-back on, global null, Holm
set.seed(777)
raw0 <- simulate_pvalues_gsd(rep(alpha, 3), alpha = alpha, corr_matrix = ds$corr,
                             info_frac = ds$info_frac, nsim = 1e6)
fw <- vapply(list(off = FALSE, on = TRUE), function(z) {
    p0 <- transform_pvalues_gsd(raw0, spending = sfLDOF, alpha = alpha, look_back = z)
    mean(rowSums(kernel_times(p0, rep(1, 3) / 3, G_holm) > 0) > 0)
}, 0)
cat(sprintf("FWER under the global null, Holm, 1e6 trials (s.e. about 0.00016): look-back off %.5f, on %.5f\n",
            fw[["off"]], fw[["on"]]))
out$look_back <- list(table = lb, fwer = fw)

# =============================================================================
# C. Mixed design
# =============================================================================
cat("\n==================== C. Mixed spending, information fractions and a late endpoint ====================\n")
infoC <- rbind(H1 = c(0.5, 1), H2 = c(0.6, 1), H3 = c(NA, 1))
pocock <- function(a, t) gsDesign::sfLDPocock(a, t)
set.seed(20261010)
rawC <- simulate_pvalues_gsd(ds$power_nom, alpha = alpha, corr_matrix = ds$corr,
                             info_frac = infoC, nsim = 1e5)
pvC  <- transform_pvalues_gsd(rawC, spending = list(sfLDOF, pocock, sfLDOF), alpha = alpha)
summary(pvC)

# the simulator's cross-hypothesis correlation at unequal information fractions
z <- stats::qnorm(unclass(rawC), lower.tail = FALSE)
emp <- c(stats::cor(z[, 1, 1], z[, 2, 1]), stats::cor(z[, 1, 1], z[, 2, 2]),
         stats::cor(z[, 1, 2], z[, 3, 2]), stats::cor(z[, 2, 1], z[, 3, 2]))
thy <- 0.5 * c(sqrt(0.5 / 0.6), sqrt(0.5), 1, sqrt(0.6))
cat("\nCorrelation of test statistics, simulated against rho * sqrt(min t / max t):\n")
print(data.frame(pair = c("H1 look 1, H2 look 1", "H1 look 1, H2 look 2",
                          "H1 look 2, H3 look 2", "H2 look 1, H3 look 2"),
                 simulated = round(emp, 4), theory = round(thy, 4)), row.names = FALSE)
cat("H3 raw p-values at the interim are all NA:", all(is.na(rawC[, 3, 1])),
    "; H3 repeated p-values at the interim are all 1:", all(pvC$pvals[, 3, 1] == 1), "\n")

set.seed(1)
elC <- system.time(
    resC <- graph_optimise_gsd(pvC, graph_constraint_free(3), gain3,
                               num_threads = n_threads, verbose = "silent")
)[["elapsed"]]
cat(sprintf("\nOptimised in %.0f s: gain %.5f, w = (%s)\n", elC, resC$power$trial_success,
            paste(round(resC$hyp_weight, 3), collapse = ", ")))
print(round(resC$trans_matrix, 3))
print(resC$power[c("local_power_by_analysis", "time_distribution")], digits = 4)

wC <- unname(resC$hyp_weight); GC <- unname(resC$trans_matrix)
tauC <- kernel_times(pvC, wC, GC)
refC <- ref_test(unclass(rawC), wC, GC, t1 = c(0.5, 0.6, NA),
                 sf = list(ref_sf_ldof, ref_sf_pocock, ref_sf_ldof), alpha = alpha)
cat(sprintf("Reference against transform + kernel on the optimised graph: %d of %d decision times differ\n",
            sum(tauC != refC), length(tauC)))
tauH <- kernel_times(pvC, rep(1, 3) / 3, G_holm)
refH <- ref_test(unclass(rawC), rep(1, 3) / 3, G_holm, t1 = c(0.5, 0.6, NA),
                 sf = list(ref_sf_ldof, ref_sf_pocock, ref_sf_ldof), alpha = alpha)
cat(sprintf("The same for Holm (every hypothesis starts with weight): %d of %d decision times differ\n",
            sum(tauH != refH), length(tauH)))
cat("H3 ever rejected at the interim:", any(tauC[, 3] == 1) || any(tauH[, 3] == 1), "\n")
out$mixed <- list(hyp_weight = resC$hyp_weight, trans_matrix = resC$trans_matrix,
                  gain = resC$power$trial_success, seconds = elC,
                  n_diff_opt = sum(tauC != refC), n_diff_holm = sum(tauH != refH),
                  corr = data.frame(simulated = emp, theory = thy))

saveRDS(out, file.path(here, "07_features_results.rds"))
