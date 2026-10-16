# B4: parallel scaling and grain size (decides O2) for the fast-shortcut build.
#
# Usage (from the repo root):
#   Rscript dev/review/shortcut_table/bench_threads.R <lib> <out_rds> \
#       [threads] [m_list] [N_list] [grains]
#
# Times graph_shortcut_parallel() with the kernel forced to "sequential"
# (auto grain) and to "table" (auto grain and each fixed grain). Compare the
# kernels only at the same thread count.

args <- commandArgs(TRUE)
lib <- args[1]
out_rds <- args[2]
parse_num <- function(x) as.numeric(strsplit(x, ",", fixed = TRUE)[[1]])
threads <- if (length(args) >= 3) parse_num(args[3]) else c(1, 2, 4, 6, 8, 12)
m_list <- if (length(args) >= 4) parse_num(args[4]) else c(8, 12, 16)
N_list <- if (length(args) >= 5) parse_num(args[5]) else c(5e4, 1e6)
grains <- if (length(args) >= 6) parse_num(args[6]) else c(64, 256, 1024, 4096, 10000)

source("dev/review/shortcut_table/common.R")
build <- load_build(lib)
alpha <- 0.025
gsp <- multigrain:::graph_shortcut_parallel

cat(sprintf("Build %s; cores: %d logical, %d physical\n", build$build,
            parallel::detectCores(), parallel::detectCores(logical = FALSE)))

time_call <- function(f, iterations) {
    bm <- bench::mark(f(), iterations = iterations, check = FALSE, memory = FALSE, filter_gc = FALSE)
    as.numeric(bm$median)
}

rows <- list()
t_start <- Sys.time()
for (m in m_list) {
    gr <- dense_graph(m, seed = 4e5 + m)
    for (N in N_list) {
        pvals <- sim_pvals(N, m, "high", seed = 5e4 + m)
        iterations <- if (N >= 1e6) 10L else 20L
        ref <- NULL
        for (nt in threads) {
            configs <- rbind(
                data.frame(kernel = "sequential", grain = -1),
                data.frame(kernel = "table", grain = c(-1, grains))
            )
            for (i in seq_len(nrow(configs))) {
                k <- configs$kernel[i]
                gsz <- configs$grain[i]
                f <- function() {
                    old <- options(multigrain.shortcut_kernel = k)
                    on.exit(options(old))
                    gsp(pvals, alpha, gr$w, gr$G, num_threads = nt, grain_size = gsz)
                }
                out <- f()
                if (is.null(ref)) ref <- out
                rows[[length(rows) + 1]] <- data.frame(
                    m = m,
                    N = N,
                    threads = nt,
                    kernel = k,
                    grain = gsz,
                    median_sec = time_call(f, iterations),
                    identical_to_first = identical(out, ref)
                )
            }
            res <- do.call(rbind, rows)
            saveRDS(res, out_rds)
            last <- res[res$m == m & res$N == N & res$threads == nt, ]
            cat(sprintf("m=%2d N=%7.0f T=%2d seq=%.2fms table(auto)=%.2fms best_table=%.2fms@%s [%.1f min]\n",
                        m, N, nt, 1e3 * last$median_sec[1], 1e3 * last$median_sec[2],
                        1e3 * min(last$median_sec[-1]), last$grain[-1][which.min(last$median_sec[-1])],
                        as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
        }
    }
}
cat("Done in", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), "min\n")
