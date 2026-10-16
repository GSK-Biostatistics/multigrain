# B2-lite: per-evaluation breakdown of the objective closure (informs O5).
#
# Usage (from the repo root), one R process per build:
#   Rscript dev/review/shortcut_table/bench_objective.R <lib> <out_rds> [N]
#
# Synthetic stand-ins for the curated problems: unconstrained graphs with
# m = 5, 6, 8, 10, a random dense candidate, N = 5e4 trials (the GA
# subsample) and two trial-success criteria. Times one call of the closure
# from create_obj_func() and, separately, the kernel alone and the
# power criterion alone; the remainder is decoding and constraint checks.

args <- commandArgs(TRUE)
lib <- args[1]
out_rds <- args[2]
N <- if (length(args) >= 3) as.numeric(args[3]) else 5e4

source("dev/review/shortcut_table/common.R")
build <- load_build(lib)
ns <- asNamespace("multigrain")
alpha <- 0.025

time_call <- function(f, iterations = 50L) {
    bm <- bench::mark(f(), iterations = iterations, check = FALSE, memory = FALSE, filter_gc = FALSE)
    as.numeric(bm$median)
}

criteria <- list(
    average = function(m) paste0("(", paste0("r", seq_len(m), collapse = " + "), ") / ", m),
    disjunctive = function(m) paste0("r", seq_len(m), collapse = " || ")
)

rows <- list()
for (m in c(5, 6, 8, 10)) {
    pvals <- sim_pvals(N, m, "high", seed = 6e4 + m)
    gc <- multigrain::graph_constraint_free(m)
    set.seed(7e4 + m)
    gr <- multigrain::graph_random(m)
    x0 <- ns$create_start_params(gc, w0 = gr$hyp_weight, G0 = gr$trans_matrix)
    theta <- ns$split_theta(x0, gc$hyp_constraint)
    w <- ns$recover_full_weights(theta$w_pars, gc$hyp_constraint)
    G <- ns$recover_full_trans_matrix(theta$g_pars, gc$trans_constraint)
    w[w < 1e-4] <- 0
    G[G < 1e-5] <- 0
    for (crit_name in names(criteria)) {
        crit_str <- criteria[[crit_name]](m)
        ts <- multigrain::trial_success(!!crit_str, verbose = "silent")
        obj <- ns$create_obj_func(
            m,
            power_criterion = ts$func,
            hyp_constraint = gc$hyp_constraint,
            trans_constraint = gc$trans_constraint,
            pvals = pvals,
            alpha = alpha,
            num_threads = 1L
        )
        rej <- ns$graph_shortcut(pvals, alpha, w, G)
        rows[[length(rows) + 1]] <- data.frame(
            build = build$build,
            m = m,
            N = N,
            criterion = crit_name,
            objective_value = obj(x0),
            eval_sec = time_call(function() obj(x0)),
            kernel_sec = time_call(function() ns$graph_shortcut(pvals, alpha, w, G)),
            criterion_sec = time_call(function() ts$func(rej)),
            mean_rej = sum(rej) / N
        )
        cat(sprintf("m=%2d %-11s eval=%.2fms kernel=%.2fms criterion=%.2fms\n", m, crit_name,
                    1e3 * rows[[length(rows)]]$eval_sec, 1e3 * rows[[length(rows)]]$kernel_sec,
                    1e3 * rows[[length(rows)]]$criterion_sec))
    }
}
res <- do.call(rbind, rows)
res$other_sec <- pmax(res$eval_sec - res$kernel_sec - res$criterion_sec, 0)
res$kernel_share <- res$kernel_sec / res$eval_sec
saveRDS(res, out_rds)
