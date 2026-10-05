# 02_three_hyp_two_stage.R
# User test of the group sequential extension (branch gsd-build).
# The same exercise as 01 with three hypotheses: a primary endpoint H1 and two
# secondaries H2 and H3 that only carry value if the primary is rejected.
#
# Run from the repository root with an *installed* build of the branch:
#   Rscript dev/review/gsd_user_testing/02_three_hyp_two_stage.R

suppressPackageStartupMessages({
    library(multigrain)
    library(gsDesign)
})

out_dir <- "dev/review/gsd_user_testing"
if (!dir.exists(out_dir)) {
    out_dir <- "."
}
n_threads <- as.integer(Sys.getenv("GSD_TEST_THREADS", "2"))

# ---- design ---------------------------------------------------------------
alpha <- 0.025
power_nom <- c(H1 = 0.90, H2 = 0.80, H3 = 0.70)
corr <- matrix(0.5, 3, 3)
diag(corr) <- 1
info_frac <- c(0.5, 1)
nsim <- 1e5
value <- c(0.5, 0.3, 0.2)
delta <- 0.8

# ---- simulate and transform ----------------------------------------------
set.seed(20261006)
t_sim <- system.time(
    raw <- simulate_pvalues_gsd(
        power_nominal = power_nom,
        alpha = alpha,
        corr_matrix = corr,
        info_frac = info_frac,
        nsim = nsim
    )
)
t_tr <- system.time(
    pvals <- transform_pvalues_gsd(raw, spending = sfLDOF, alpha = alpha)
)
summary(pvals)

# ---- gain -------------------------------------------------------------------
# Hurdle on the primary: a secondary claim has value only if H1 is rejected.
# Each claim is discounted by the analysis at which it was declared.
gain_disc <- trial_success_gsd(
    !!value[1] * d(t1) + r1 * (!!value[2] * d(t2) + !!value[3] * d(t3)),
    d = c(1, delta)
)
gain_flat <- trial_success_gsd(
    !!value[1] * r1 + r1 * (!!value[2] * r2 + !!value[3] * r3),
    K = 2
)
print(gain_disc)

# ---- optimise -------------------------------------------------------------
optimise <- function(gain, seed = 1) {
    set.seed(seed)
    elapsed <- system.time(
        res <- graph_optimise_gsd(
            pvals = pvals,
            graph_constraint = graph_constraint_free(3),
            trial_success = gain,
            num_threads = n_threads
        )
    )[["elapsed"]]
    list(res = res, elapsed = elapsed)
}
opt_disc <- optimise(gain_disc, seed = 1)
opt_disc2 <- optimise(gain_disc, seed = 2) # second seed: search stability
opt_flat <- optimise(gain_flat, seed = 1)
res <- opt_disc$res

print(res)
summary(res)
cat("\nPower fields that print() and summary() do not show:\n")
print(res$power[c(
    "local_power_by_analysis",
    "mean_decision_look",
    "time_distribution"
)])

# ---- comparators on the same draws -----------------------------------------
eval_graph <- function(w, G) {
    calc_power_pvals_gsd(
        pvals,
        hyp_weight = w,
        trans_matrix = G,
        custom_power = list(disc = gain_disc, flat = gain_flat)
    )
}
G_holm <- matrix(0.5, 3, 3)
diag(G_holm) <- 0
G_seq <- rbind(c(0, 1, 0), c(0, 0, 1), c(1, 0, 0))
G_gate <- rbind(c(0, 0.5, 0.5), c(0, 0, 1), c(0, 1, 0))
graphs <- list(
    `optimised (discounted, seed 1)` = list(
        w = res$hyp_weight,
        G = res$trans_matrix
    ),
    `optimised (discounted, seed 2)` = list(
        w = opt_disc2$res$hyp_weight,
        G = opt_disc2$res$trans_matrix
    ),
    `optimised (undiscounted)` = list(
        w = opt_flat$res$hyp_weight,
        G = opt_flat$res$trans_matrix
    ),
    `Holm` = list(w = rep(1, 3) / 3, G = G_holm),
    `fixed sequence H1 -> H2 -> H3` = list(w = c(1, 0, 0), G = G_seq),
    `H1 gate, then Holm on H2, H3` = list(w = c(1, 0, 0), G = G_gate)
)
comp <- do.call(
    rbind,
    lapply(names(graphs), function(nm) {
        g <- graphs[[nm]]
        p <- eval_graph(g$w, g$G)
        data.frame(
            graph = nm,
            gain_discounted = p$disc,
            gain_undiscounted = p$flat,
            power_H1 = p$local_power[1],
            power_H2 = p$local_power[2],
            power_H3 = p$local_power[3],
            interim_H1 = p$local_power_by_analysis[1, 1],
            interim_H2 = p$local_power_by_analysis[2, 1],
            interim_H3 = p$local_power_by_analysis[3, 1]
        )
    })
)
cat("\nComparators (same simulated trials):\n")
print(comp, digits = 4, row.names = FALSE)
cat("\nOptimised graphs:\n")
for (nm in names(graphs)[1:3]) {
    cat(nm, ": w =", paste(round(graphs[[nm]]$w, 3), collapse = ", "), "\n")
    print(round(graphs[[nm]]$G, 3))
}

# ---- independent optimality check: random search ----------------------------
# 20,000 random graphs (Dirichlet weights, uniform transition splits, no alpha
# discarded) evaluated on the same draws. The optimiser should not be beaten
# by more than Monte Carlo ties.
random_graph <- function() {
    w <- stats::rexp(3)
    w <- w / sum(w)
    u <- stats::runif(3)
    G <- rbind(c(0, u[1], 1 - u[1]), c(u[2], 0, 1 - u[2]), c(u[3], 1 - u[3], 0))
    list(w = w, G = G)
}
kernel_mat <- pvals$pvals
dim(kernel_mat) <- c(nsim, 3 * 2)
gain_of <- function(g) {
    out <- multigrain:::graph_shortcut_gsd_parallel(
        kernel_mat,
        alpha,
        g$w,
        g$G,
        2L,
        n_threads,
        1000L
    )
    gain_disc$func(out$time)
}
set.seed(99)
t_rand <- system.time({
    cand <- replicate(20000, random_graph(), simplify = FALSE)
    rand_gain <- vapply(cand, gain_of, 0)
})
best <- cand[[which.max(rand_gain)]]
cat(sprintf(
    "\nRandom search over 20,000 graphs: best gain %.5f; optimiser %.5f (seed 1), %.5f (seed 2)\n",
    max(rand_gain),
    opt_disc$res$power$trial_success,
    opt_disc2$res$power$trial_success
))
cat("Best random graph: w =", paste(round(best$w, 3), collapse = ", "), "\n")
print(round(best$G, 3))
cat("Quantiles of the random-search gain:\n")
print(round(stats::quantile(rand_gain, c(0, 0.5, 0.9, 0.99, 0.999, 1)), 5))

timing <- c(
    simulate = t_sim[["elapsed"]],
    transform = t_tr[["elapsed"]],
    optimise_discounted_seed1 = opt_disc$elapsed,
    optimise_discounted_seed2 = opt_disc2$elapsed,
    optimise_undiscounted = opt_flat$elapsed,
    random_search_20000 = t_rand[["elapsed"]]
)
cat("\nElapsed seconds (", n_threads, " threads):\n", sep = "")
print(round(timing, 1))

keep <- function(o) {
    list(
        hyp_weight = o$res$hyp_weight,
        trans_matrix = o$res$trans_matrix,
        power = o$res$power,
        solution = o$res$solution
    )
}
saveRDS(
    list(
        design = list(
            alpha = alpha,
            power_nom = power_nom,
            corr = corr,
            info_frac = info_frac,
            nsim = nsim,
            value = value,
            delta = delta,
            spending = "sfLDOF",
            seed = 20261006
        ),
        optimised = list(
            discounted = keep(opt_disc),
            discounted_seed2 = keep(opt_disc2),
            undiscounted = keep(opt_flat)
        ),
        comparators = comp,
        random_search = list(
            best_gain = max(rand_gain),
            best_graph = best,
            quantiles = stats::quantile(
                rand_gain,
                c(0, 0.5, 0.9, 0.99, 0.999, 1)
            )
        ),
        timing = timing,
        threads = n_threads,
        session = c(
            R = R.version.string,
            multigrain = as.character(utils::packageVersion("multigrain")),
            gsDesign = as.character(utils::packageVersion("gsDesign"))
        )
    ),
    file.path(out_dir, "02_three_hyp_two_stage_results.rds")
)
