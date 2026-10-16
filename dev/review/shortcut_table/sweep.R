# Equivalence sweep: forced "sequential" vs "table" kernels (design 6.2).
#
# Usage (from the repo root):
#   Rscript dev/review/shortcut_table/sweep.R <lib> <out_dir> [m_max]
#
# For each m in 2..m_max (default 16): 50 dense and 50 sparse random graphs,
# N = 5e4 synthetic p-values at three power levels. Counts indicator
# differences between the forced kernels (serial) and spot-checks the
# parallel table kernel (num_threads = 4) on the first 5 graphs of each cell.
# Any difference is replayed in R and classified as a floating-point tie
# (relative gap <= 1e-12 between some p_j and its critical value) or a defect.

args <- commandArgs(TRUE)
lib <- args[1]
out_dir <- args[2]
m_max <- if (length(args) >= 3) as.integer(args[3]) else 16L
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

source("dev/review/shortcut_table/common.R")
build <- load_build(lib)

N <- 5e4
alpha <- 0.025
n_graphs <- 50
tie_tol <- 1e-12

rows <- list()
diff_log <- list()
t_start <- Sys.time()

for (m in 2:m_max) {
    for (power in c("high", "low", "very_high")) {
        pvals <- sim_pvals(N, m, power, seed = 1000 * m + match(power, c("high", "low", "very_high")))
        for (type in c("dense", "sparse")) {
            n_diff <- 0
            n_tie <- 0
            n_defect <- 0
            n_par_diff <- 0
            rej_sum <- 0
            for (g in seq_len(n_graphs)) {
                seed <- 1e5 * m + 1e3 * g + (type == "sparse")
                gr <- if (type == "dense") dense_graph(m, seed) else sparse_graph(m, seed)

                out_seq <- shortcut_forced("sequential", pvals, alpha, gr$w, gr$G)
                out_tab <- shortcut_forced("table", pvals, alpha, gr$w, gr$G)
                rej_sum <- rej_sum + sum(out_seq)

                if (g <= 5) {
                    out_par <- shortcut_forced(
                        "table",
                        pvals,
                        alpha,
                        gr$w,
                        gr$G,
                        num_threads = 4
                    )
                    n_par_diff <- n_par_diff + sum(out_par != out_tab)
                }

                d <- which(out_seq != out_tab, arr.ind = TRUE)
                if (nrow(d) > 0) {
                    n_diff <- n_diff + nrow(d)
                    crit <- multigrain:::graph_shortcut_table(alpha, gr$w, gr$G)
                    for (t in unique(d[, 1])) {
                        rp <- replay_trial(pvals[t, ], gr$w, gr$G, alpha, crit)
                        is_tie <- rp$min_rel_gap <= tie_tol
                        n_j <- sum(d[, 1] == t)
                        if (is_tie) n_tie <- n_tie + n_j else n_defect <- n_defect + n_j
                        diff_log[[length(diff_log) + 1]] <- list(
                            m = m,
                            power = power,
                            type = type,
                            graph = g,
                            trial = t,
                            p = pvals[t, ],
                            w = gr$w,
                            G = gr$G,
                            rej_seq_kernel = out_seq[t, ],
                            rej_tab_kernel = out_tab[t, ],
                            replay = rp,
                            class = if (is_tie) "fp_tie" else "defect"
                        )
                    }
                }
            }
            rows[[length(rows) + 1]] <- data.frame(
                m = m,
                power = power,
                type = type,
                graphs = n_graphs,
                indicators = n_graphs * N * m,
                mean_rej_per_trial = rej_sum / (n_graphs * N),
                diffs = n_diff,
                fp_ties = n_tie,
                defects = n_defect,
                parallel_diffs = n_par_diff
            )
            cat(sprintf(
                "m=%2d %-4s %-6s rej/trial=%.2f diffs=%d (ties=%d, defects=%d) par_diffs=%d  [%.0fs]\n",
                m, power, type, rej_sum / (n_graphs * N), n_diff, n_tie, n_defect,
                n_par_diff, as.numeric(difftime(Sys.time(), t_start, units = "secs"))
            ))
        }
    }
}

summary_df <- do.call(rbind, rows)
saveRDS(
    list(build = build, summary = summary_df, diffs = diff_log, N = N, alpha = alpha),
    file.path(out_dir, "sweep.rds")
)

sink(file.path(out_dir, "sweep.txt"))
cat("Build:", build$build, "(", build$version, ")\n")
cat("N =", N, " alpha =", alpha, " graphs per cell =", n_graphs, "\n\n")
print(summary_df, row.names = FALSE)
cat("\nTotal indicators:", format(sum(summary_df$indicators), big.mark = ","), "\n")
cat("Total differences:", sum(summary_df$diffs), "\n")
cat("  floating-point ties:", sum(summary_df$fp_ties), "\n")
cat("  defects:", sum(summary_df$defects), "\n")
cat("Parallel (4 threads) vs serial table differences:", sum(summary_df$parallel_diffs), "\n")
if (length(diff_log) > 0) {
    cat("\nDifferences:\n")
    for (d in diff_log) {
        cat(sprintf(
            "m=%d %s %s graph=%d trial=%d class=%s min_rel_gap=%.3g\n",
            d$m, d$power, d$type, d$graph, d$trial, d$class, d$replay$min_rel_gap
        ))
    }
}
sink()
cat("Done in", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), "min\n")
