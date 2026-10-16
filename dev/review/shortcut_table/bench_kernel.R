# B1: serial kernel microbenchmark (decides O1, the switch rule).
#
# Usage (from the repo root), one R process per build:
#   Rscript dev/review/shortcut_table/bench_kernel.R <lib> <out_rds> \
#       [m_list] [N_list] [graphs] [iterations]
# e.g.
#   Rscript .../bench_kernel.R lib/v1.0.0 results/b1_v1.0.0.rds
#   Rscript .../bench_kernel.R lib/fast-shortcut results/b1_fast.rds
#
# A build with graph_shortcut_kernel() (the fast-shortcut build) is timed with
# the kernel forced to "sequential" and to "table", and the table build alone
# is timed as a forced-table call on a single trial. Other builds (v1.0.0) are
# timed as they are. Every cell records a hash of the rejection matrix so
# builds and kernels can be checked for identical output.

args <- commandArgs(TRUE)
lib <- args[1]
out_rds <- args[2]
parse_num <- function(x) as.numeric(strsplit(x, ",", fixed = TRUE)[[1]])
m_list <- if (length(args) >= 3) parse_num(args[3]) else c(3, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20)
N_list <- if (length(args) >= 4) parse_num(args[4]) else c(1e4, 5e4, 2e5, 1e6)
n_graphs <- if (length(args) >= 5) as.integer(args[5]) else 2L
iterations <- if (length(args) >= 6) as.integer(args[6]) else 20L

source("dev/review/shortcut_table/common.R")
build <- load_build(lib)
has_table <- exists("graph_shortcut_kernel", envir = asNamespace("multigrain"))
kernels <- if (has_table) c("sequential", "table") else "default"
gs <- multigrain:::graph_shortcut
alpha <- 0.025

time_call <- function(f) {
    bm <- bench::mark(f(), iterations = iterations, check = FALSE, memory = FALSE, filter_gc = FALSE)
    c(median = as.numeric(bm$median), min = as.numeric(bm$min))
}

cat(sprintf("Build %s (%s), kernels: %s, iterations %d, graphs %d\n",
            build$build, build$version, paste(kernels, collapse = "/"), iterations, n_graphs))
rows <- list()
t_start <- Sys.time()

for (m in m_list) {
    for (N in N_list) {
        for (power in c("high", "low")) {
            pvals <- sim_pvals(N, m, power, seed = 2e4 + 100 * m + (power == "low"))
            for (g in seq_len(n_graphs)) {
                gr <- dense_graph(m, seed = 3e5 + 100 * m + g)
                build_sec <- NA_real_
                auto_kernel <- NA_character_
                if (has_table) {
                    auto_kernel <- multigrain:::graph_shortcut_kernel(N, m)
                    p1 <- pvals[1, , drop = FALSE]
                    build_sec <- time_call(function() {
                        shortcut_forced("table", p1, alpha, gr$w, gr$G)
                    })[["median"]]
                }
                for (k in kernels) {
                    f <- if (k == "default") {
                        function() gs(pvals, alpha, gr$w, gr$G)
                    } else {
                        function() shortcut_forced(k, pvals, alpha, gr$w, gr$G)
                    }
                    out <- f()
                    tm <- time_call(f)
                    rows[[length(rows) + 1]] <- data.frame(
                        build = build$build,
                        version = build$version,
                        kernel = k,
                        m = m,
                        N = N,
                        power = power,
                        graph = g,
                        median_sec = tm[["median"]],
                        min_sec = tm[["min"]],
                        build_sec = if (k == "table") build_sec else NA_real_,
                        auto_kernel = auto_kernel,
                        mean_rej = sum(out) / N,
                        hash = rlang::hash(unname(out))
                    )
                    rm(out)
                }
            }
        }
        res <- do.call(rbind, rows)
        saveRDS(res, out_rds)
        last <- res[res$m == m & res$N == N, ]
        cat(sprintf("m=%2d N=%7.0f  %s  [%.1f min]\n", m, N,
                    paste(sprintf("%s=%.2fms", last$kernel, 1e3 * last$median_sec)[seq_len(min(4, nrow(last)))], collapse = " "),
                    as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
    }
}
cat("Done in", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), "min\n")
