# Builder check (design 6.2): critical-value table from the fast builder
# against calc_local_weights() and graphicalMCP::graph_generate_weights().
#
# Usage (from the repo root):
#   Rscript dev/review/shortcut_table/builder_check.R <lib> <out_dir>
#
# For m = 2..10, 20 dense and 20 sparse random graphs. Rows are matched by
# their membership columns, not by assumed row order. Reports the maximum
# absolute difference on the weight scale (table / alpha) and on the
# critical-value scale (table vs alpha * weights).

args <- commandArgs(TRUE)
lib <- args[1]
out_dir <- args[2]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

source("dev/review/shortcut_table/common.R")
build <- load_build(lib)
has_gmcp <- requireNamespace("graphicalMCP", quietly = TRUE)

alpha <- 0.025

# Weights from a (2^m - 1) x 2m membership/weight matrix, placed in the
# table's row order (row mask + 1, bit j - 1 = H_j).
to_table_order <- function(lw, m) {
    member <- lw[, seq_len(m), drop = FALSE]
    mask <- as.vector(member %*% 2^(seq_len(m) - 1))
    out <- matrix(0, 2^m, m)
    out[mask + 1, ] <- unname(lw[, m + seq_len(m), drop = FALSE])
    stopifnot(!anyDuplicated(mask), length(mask) == 2^m - 1)
    out
}

rows <- list()
for (m in 2:10) {
    for (type in c("dense", "sparse")) {
        max_w_clw <- 0
        max_w_gmcp <- 0
        max_c_clw <- 0
        max_full_row <- 0
        for (g in 1:20) {
            seed <- 7e5 + 1e3 * m + g
            gr <- if (type == "dense") dense_graph(m, seed) else sparse_graph(m, seed)
            tab <- multigrain:::graph_shortcut_table(alpha, gr$w, gr$G)

            clw <- to_table_order(multigrain:::calc_local_weights(gr$w, gr$G), m)
            max_w_clw <- max(max_w_clw, abs(tab / alpha - clw))
            max_c_clw <- max(max_c_clw, abs(tab - alpha * clw))
            max_full_row <- max(max_full_row, abs(tab[2^m, ] - gr$w * alpha))

            if (has_gmcp) {
                graph <- graphicalMCP::graph_create(gr$w, gr$G)
                gw <- to_table_order(graphicalMCP::graph_generate_weights(graph), m)
                max_w_gmcp <- max(max_w_gmcp, abs(tab / alpha - gw))
            }
        }
        rows[[length(rows) + 1]] <- data.frame(
            m = m,
            type = type,
            graphs = 20,
            max_abs_w_vs_calc_local_weights = max_w_clw,
            max_abs_w_vs_graphicalMCP = if (has_gmcp) max_w_gmcp else NA_real_,
            max_abs_crit_vs_alpha_clw = max_c_clw,
            full_set_row_vs_w_alpha = max_full_row
        )
    }
}

summary_df <- do.call(rbind, rows)
saveRDS(list(build = build, summary = summary_df), file.path(out_dir, "builder_check.rds"))

sink(file.path(out_dir, "builder_check.txt"))
cat("Build:", build$build, "(", build$version, ")  alpha =", alpha, "\n")
cat("graphicalMCP available:", has_gmcp, "\n\n")
print(summary_df, row.names = FALSE, digits = 3)
cat("\nOverall max |table/alpha - calc_local_weights|:", format(max(summary_df$max_abs_w_vs_calc_local_weights), digits = 3), "\n")
cat("Overall max |table/alpha - graphicalMCP|:", format(max(summary_df$max_abs_w_vs_graphicalMCP), digits = 3), "\n")
cat("Full-set row exactly w * alpha:", all(summary_df$full_set_row_vs_w_alpha == 0), "\n")
sink()
print(summary_df, row.names = FALSE, digits = 3)
