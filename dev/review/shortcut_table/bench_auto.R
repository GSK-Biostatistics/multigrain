# B5-local: "auto" against the forced kernels, timed directly, for the
# fast-shortcut build. Checks that the kernel choice and option reads add no
# measurable overhead, at small calls and either side of the switch boundary.
#
# Usage (from the repo root):
#   Rscript dev/review/shortcut_table/bench_auto.R <lib> <out_rds>

args <- commandArgs(TRUE)
lib <- args[1]
out_rds <- args[2]

source("dev/review/shortcut_table/common.R")
build <- load_build(lib)
alpha <- 0.025
cells <- data.frame(
    m = c(3, 8, 8, 16, 17, 18, 19, 20),
    N = c(1e4, 1e4, 5e4, 5e4, 5e4, 2e5, 2e5, 1e6)
)

rows <- list()
for (i in seq_len(nrow(cells))) {
    m <- cells$m[i]
    N <- cells$N[i]
    for (power in c("high", "low")) {
        pvals <- sim_pvals(N, m, power, seed = 8e4 + 100 * m + (power == "low"))
        gr <- dense_graph(m, seed = 9e5 + m)
        iters <- if (N >= 1e6) 10L else 20L
        bm <- bench::mark(
            auto = shortcut_forced("auto", pvals, alpha, gr$w, gr$G),
            sequential = shortcut_forced("sequential", pvals, alpha, gr$w, gr$G),
            table = shortcut_forced("table", pvals, alpha, gr$w, gr$G),
            iterations = iters,
            check = TRUE,
            memory = FALSE,
            filter_gc = FALSE
        )
        t <- setNames(as.numeric(bm$median), as.character(bm$expression))
        rows[[length(rows) + 1]] <- data.frame(
            m = m,
            N = N,
            power = power,
            auto_kernel = multigrain:::graph_shortcut_kernel(N, m),
            auto_ms = 1e3 * t[["auto"]],
            sequential_ms = 1e3 * t[["sequential"]],
            table_ms = 1e3 * t[["table"]],
            auto_vs_sequential = t[["auto"]] / t[["sequential"]],
            auto_vs_best = t[["auto"]] / min(t[["sequential"]], t[["table"]])
        )
        print(rows[[length(rows)]], row.names = FALSE, digits = 3)
    }
}
res <- do.call(rbind, rows)
saveRDS(res, out_rds)
