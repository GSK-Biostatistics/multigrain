# 04_stopping_rules.R
# Prototype for the stopping-criteria gap. No package code is changed.
#
# The kernel always analyses every look, so a hypothesis can be "rejected" at
# the final analysis of a trial that would in practice have stopped at the
# interim. This script shows, on a two-hypothesis two-look design:
#
#   1. An efficacy stopping rule that depends only on which hypotheses have
#      been rejected so far is a deterministic edit of the decision-time
#      matrix: find the stopping look S, then void every rejection after S.
#      Three implementations agree exactly: (a) that edit applied in R to the
#      kernel output, (b) the same rule written in the existing
#      trial_success_gsd() grammar, (c) an independent reference that really
#      stops analysing the trial.
#   2. What the stopping rule does to the optimal graph and to the value of
#      the design, under two readings of "when is a claim worth its value".
#   3. That the spending function of the secondary endpoint, a fixed input
#      today, becomes a first-order design choice once the trial can stop.
#   4. That a non-binding futility rule can be applied as a mask on the
#      transformed p-values before the search.
#
# Run from the repository root:
#   Rscript dev/review/gsd_user_testing/04_stopping_rules.R

suppressPackageStartupMessages({
    library(multigrain)
    library(gsDesign)
})

here <- "dev/review/gsd_user_testing"
if (!dir.exists(here)) here <- "."
source(file.path(here, "00_reference.R"))
n_threads <- as.integer(Sys.getenv("GSD_TEST_THREADS", "2"))

# ---- design: H1 primary (drives stopping), H2 key secondary ----------------
alpha     <- 0.025
power_nom <- c(H1 = 0.90, H2 = 0.80)
corr      <- matrix(c(1, 0.5, 0.5, 1), 2)
info_frac <- c(0.5, 1)
nsim      <- 2e5
value     <- c(0.6, 0.4)
delta     <- 0.8
G_full    <- rbind(c(0, 1), c(1, 0))

set.seed(20261007)
raw <- simulate_pvalues_gsd(power_nom, alpha = alpha, corr_matrix = corr,
                            info_frac = info_frac, nsim = nsim)

pocock <- function(a, t) gsDesign::sfLDPocock(a, t)
spending <- list(
    `H1 LDOF, H2 LDOF`   = list(pkg = list(sfLDOF, sfLDOF),
                                ref = list(ref_sf_ldof, ref_sf_ldof)),
    `H1 LDOF, H2 Pocock` = list(pkg = list(sfLDOF, pocock),
                                ref = list(ref_sf_ldof, ref_sf_pocock))
)
pv <- lapply(spending, function(s) {
    transform_pvalues_gsd(raw, spending = s$pkg, alpha = alpha)
})
for (nm in names(pv)) {
    cat("\n", nm, "\n", sep = "")
    summary(pv[[nm]])
}

kernel_times <- function(p, w, G = G_full) {
    x <- p$pvals
    dim(x) <- c(p$nsim, p$m * p$K)
    multigrain:::graph_shortcut_gsd_parallel(x, alpha, w, G, p$K, n_threads, 1000L)$time
}

# ---- the stopping rule as an edit of the decision-time matrix --------------
# stop_when(tau, k) returns, per trial, whether the trial stops after look k
# given the rejections declared at looks <= k. Here: stop once the primary
# has been rejected.
stop_when_primary <- function(tau, k) tau[, 1] > 0 & tau[, 1] <= k

apply_stopping <- function(tau, stop_when, K = 2L) {
    S <- rep(K, nrow(tau))
    for (k in rev(seq_len(K - 1L))) S[stop_when(tau, k)] <- k
    tau[tau > S] <- 0L            # nothing can be rejected after the trial stops
    list(time = tau, stop_look = S)
}

# ---- 1. three implementations of "stop once H1 is rejected" ----------------
w_test <- c(0.8, 0.2)
d_tab  <- c(0, 1, delta)
gain_r <- function(tau) mean(value[1] * d_tab[tau[, 1] + 1L] + value[2] * d_tab[tau[, 2] + 1L])

# (b) existing grammar: H2 at the final analysis only counts if the trial is
# still running, that is if H1 was not rejected at the interim
gain_stop <- trial_success_gsd(
    !!value[1] * d(t1) + !!value[2] * ((t2 == 1) + !!delta * (t2 == 2 && t1 != 1)),
    d = c(1, delta)
)
gain_nostop <- trial_success_gsd(
    !!value[1] * d(t1) + !!value[2] * d(t2),
    d = c(1, delta), verbose = "silent"
)

cat("\n1. Three implementations of the stopping rule, graph w = (0.8, 0.2), full recycling\n")
n_ref <- 5e4
for (nm in names(pv)) {
    tau  <- kernel_times(pv[[nm]], w_test)
    edit <- apply_stopping(tau, stop_when_primary)
    g_a  <- gain_r(edit$time)
    g_b  <- gain_stop$func(tau)
    tau_ref <- ref_test(unclass(raw)[seq_len(n_ref), , , drop = FALSE], w_test, G_full,
                        t1 = info_frac[1], sf = spending[[nm]]$ref, alpha = alpha,
                        stop_fun = function(tr) tr[1] > 0)
    cat(sprintf(
        "   %-20s (a) edit in R %.8f   (b) existing grammar %.8f   difference %.1e\n",
        nm, g_a, g_b, g_a - g_b
    ))
    cat(sprintf(
        "   %-20s (c) reference that stops, first %d trials: %d of %d decision times differ from (a)\n",
        "", n_ref, sum(tau_ref != edit$time[seq_len(n_ref), ]), length(tau_ref)
    ))
}

# ---- 2 and 3. profiles over w1 under three value models ---------------------
# All three gains are linear in the joint distribution of (t1, t2), so one
# kernel run per graph gives every model.
#   M0  claim-level timing, no stopping (what 01 optimises)
#   M1  claim-level timing, stop once H1 is rejected
#   M2  trial-level timing, stop once H1 is rejected: every claim is realised
#       when the trial ends, so the discount applies to the stopping look
joint <- function(tau) {
    table(factor(tau[, 1], 0:2), factor(tau[, 2], 0:2)) / nrow(tau)
}
psi <- list(
    M0 = outer(0:2, 0:2, function(a, b) value[1] * d_tab[a + 1] + value[2] * d_tab[b + 1]),
    M1 = outer(0:2, 0:2, function(a, b) {
        value[1] * d_tab[a + 1] + value[2] * ifelse(a == 1 & b == 2, 0, d_tab[b + 1])
    }),
    M2 = outer(0:2, 0:2, function(a, b) {
        ifelse(a == 1, value[1] + value[2] * (b == 1),
               delta * (value[1] * (a > 0) + value[2] * (b > 0)))
    })
)
w_grid <- seq(0, 1, by = 0.02)
profiles <- list()
for (nm in names(pv)) {
    jt <- lapply(w_grid, function(w1) joint(kernel_times(pv[[nm]], c(w1, 1 - w1))))
    profiles[[nm]] <- data.frame(
        spending = nm, w1 = w_grid,
        M0 = vapply(jt, function(j) sum(j * psi$M0), 0),
        M1 = vapply(jt, function(j) sum(j * psi$M1), 0),
        M2 = vapply(jt, function(j) sum(j * psi$M2), 0),
        power_H1 = vapply(jt, function(j) sum(j[2:3, ]), 0),
        power_H2_nostop = vapply(jt, function(j) sum(j[, 2:3]), 0),
        power_H2_stop   = vapply(jt, function(j) sum(j[, 2:3]) - j[2, 3], 0),
        p_stop_interim  = vapply(jt, function(j) sum(j[2, ]), 0)
    )
}
prof <- do.call(rbind, profiles)

summarise <- function(nm, model) {
    p <- profiles[[nm]]
    g <- p[[model]]
    i <- which.max(g)
    near <- range(p$w1[g > g[i] - 5e-4])
    data.frame(
        spending = nm, model = model,
        w1_opt = p$w1[i], plateau = sprintf("[%.2f, %.2f]", near[1], near[2]),
        gain_opt = g[i], gain_holm = g[p$w1 == 0.5], gain_fixed_seq = g[p$w1 == 1],
        power_H1 = p$power_H1[i],
        power_H2 = if (model == "M0") p$power_H2_nostop[i] else p$power_H2_stop[i],
        p_stop_interim = if (model == "M0") 0 else p$p_stop_interim[i]
    )
}
tab <- do.call(rbind, lapply(names(pv), function(nm) {
    do.call(rbind, lapply(c("M0", "M1", "M2"), function(mo) summarise(nm, mo)))
}))
cat("\n2. Optimal w1 (grid, full recycling) and gain under each value model\n")
print(tab, digits = 4, row.names = FALSE)

cat("\n3. Cost of ignoring the stopping rule when choosing the graph\n")
for (nm in names(pv)) {
    p <- profiles[[nm]]
    for (mo in c("M1", "M2")) {
        i0 <- which.max(p$M0)
        cat(sprintf(
            "   %-20s %s: M0-optimal graph (w1 = %.2f) scores %.4f under %s; the %s optimum scores %.4f (regret %.4f). M0 reported %.4f for it.\n",
            nm, mo, p$w1[i0], p[[mo]][i0], mo, mo, max(p[[mo]]),
            max(p[[mo]]) - p[[mo]][i0], p$M0[i0]
        ))
    }
}

# ---- end to end: optimise under the stopping rule with today's package -----
cat("\n4. graph_optimise_gsd() with the stopping rule written into the gain (H1 LDOF, H2 Pocock)\n")
set.seed(1)
t_opt <- system.time(
    res <- graph_optimise_gsd(
        pvals = pv[["H1 LDOF, H2 Pocock"]],
        graph_constraint = graph_constraint_free(2),
        trial_success = gain_stop,
        num_threads = n_threads,
        verbose = "silent"
    )
)[["elapsed"]]
p_b <- profiles[["H1 LDOF, H2 Pocock"]]
cat(sprintf(
    "   optimiser: w1 = %.4f, gain %.5f (%.0f s); grid: w1 = %.2f, gain %.5f\n",
    res$hyp_weight[1], res$power$trial_success, t_opt,
    p_b$w1[which.max(p_b$M1)], max(p_b$M1)
))
cat("   What calc_power_pvals_gsd() reports for that graph ignores the stopping rule:\n")
cat(sprintf("   reported power H2 = %.4f; power H2 under the stopping rule = %.4f\n",
            res$power$local_power[2],
            mean(apply_stopping(kernel_times(pv[["H1 LDOF, H2 Pocock"]],
                                             res$hyp_weight, res$trans_matrix),
                                stop_when_primary)$time[, 2] > 0)))

# ---- 5. giving early stopping a value of its own ----------------------------
# M1 plus a saving c every time the trial stops at the interim, on the scale of
# the claim values (the total claim value is 1). In today's grammar this is
# `+ c * (t1 == 1)`. Without such a term the optimiser is paid only for claims,
# and under M1 it protects H2 by making an interim rejection of H1 less likely.
cat("\n5. Optimal w1 under M1 plus a saving c for stopping at the interim\n")
saving_tab <- do.call(rbind, lapply(names(pv), function(nm) {
    p <- profiles[[nm]]
    do.call(rbind, lapply(c(0, 0.02, 0.05, 0.10, 0.20), function(cc) {
        g <- p$M1 + cc * p$p_stop_interim
        i <- which.max(g)
        data.frame(spending = nm, saving = cc, w1_opt = p$w1[i], gain_opt = g[i],
                   p_stop_interim = p$p_stop_interim[i],
                   power_H1 = p$power_H1[i], power_H2 = p$power_H2_stop[i])
    }))
}))
print(saving_tab, digits = 4, row.names = FALSE)

# ---- 6. non-binding futility as a mask on the transformed p-values ---------
# Rule: stop after the interim when the primary looks unpromising, here when
# its raw interim p-value is above 0.3. The rule reads raw data only, so it
# does not depend on the graph, and being non-binding it leaves the efficacy
# boundaries alone. It can therefore be applied once, before the search, by
# setting every repeated p-value at the final analysis to 1 in the trials it
# stops. Checked against the reference stopping those trials for real, alone
# and together with the efficacy rule of section 1.
cat("\n6. Futility (stop if the raw interim p-value of H1 exceeds 0.3)\n")
futile <- unclass(raw)[, 1, 1] > 0.3
cat(sprintf("   trials stopped for futility: %.4f\n", mean(futile)))
mask_futility <- function(p) {
    p$pvals[futile, , 2] <- 1
    p
}
pv_fut <- lapply(pv, mask_futility)
for (nm in names(pv)) {
    tau_m <- kernel_times(pv_fut[[nm]], w_test)
    ref_f <- ref_test(unclass(raw)[seq_len(n_ref), , , drop = FALSE], w_test, G_full,
                      t1 = info_frac[1], sf = spending[[nm]]$ref, alpha = alpha,
                      stop_after_1 = futile[seq_len(n_ref)])
    both_m <- apply_stopping(tau_m, stop_when_primary)$time
    ref_b  <- ref_test(unclass(raw)[seq_len(n_ref), , , drop = FALSE], w_test, G_full,
                       t1 = info_frac[1], sf = spending[[nm]]$ref, alpha = alpha,
                       stop_fun = function(tr) tr[1] > 0,
                       stop_after_1 = futile[seq_len(n_ref)])
    cat(sprintf(
        "   %-20s futility only: %d of %d decision times differ from the reference; with the efficacy rule: %d of %d\n",
        nm, sum(ref_f != tau_m[seq_len(n_ref), ]), length(ref_f),
        sum(ref_b != both_m[seq_len(n_ref), ]), length(ref_b)
    ))
}
fut_tab <- do.call(rbind, lapply(names(pv), function(nm) {
    jt <- lapply(w_grid, function(w1) joint(kernel_times(pv_fut[[nm]], c(w1, 1 - w1))))
    do.call(rbind, lapply(c("M0", "M1"), function(mo) {
        g  <- vapply(jt, function(j) sum(j * psi[[mo]]), 0)
        g0 <- profiles[[nm]][[mo]]
        data.frame(spending = nm, model = mo,
                   w1_opt_no_futility = w_grid[which.max(g0)], gain_no_futility = max(g0),
                   w1_opt_futility = w_grid[which.max(g)], gain_futility = max(g),
                   regret_of_ignoring = max(g) - g[which.max(g0)])
    }))
}))
cat("   Optimal w1 and gain with and without the futility rule:\n")
print(fut_tab, digits = 4, row.names = FALSE)

saveRDS(
    list(design = list(alpha = alpha, power_nom = power_nom, corr = corr,
                       info_frac = info_frac, nsim = nsim, value = value,
                       delta = delta, seed = 20261007),
         table = tab, profiles = prof, saving = saving_tab, futility = fut_tab,
         optimiser = list(hyp_weight = res$hyp_weight, trans_matrix = res$trans_matrix,
                          gain = res$power$trial_success, elapsed = t_opt)),
    file.path(here, "04_stopping_rules_results.rds")
)
