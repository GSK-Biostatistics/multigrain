# Compare two clean benchmark result sets (e.g. results/run2 and results/run3)
# on every headline quantity in the report: B1 speedups, the same-code
# control, break-even, G4, B4 table/sequential ratios, B2 per-evaluation
# speedups and the B5 auto check.
#
# Usage (from the repo root):
#   Rscript dev/review/shortcut_table/compare_runs.R <dir_a> <dir_b>
# Run analyse_b1.R on both directories first.

args <- commandArgs(TRUE)
dir_a <- args[1]
dir_b <- args[2]
lab <- basename(c(dir_a, dir_b))

rel <- function(a, b) b / a - 1
summ <- function(x) {
    sprintf("median %+.1f%%, IQR [%+.1f%%, %+.1f%%], range [%+.1f%%, %+.1f%%]",
            100 * median(x), 100 * quantile(x, 0.25), 100 * quantile(x, 0.75),
            100 * min(x), 100 * max(x))
}

# ---- B1 -------------------------------------------------------------------
an <- lapply(c(dir_a, dir_b), function(d) readRDS(file.path(d, "b1_analysis.rds")))
key <- c("m", "N", "power")
b1 <- merge(an[[1]]$agg, an[[2]]$agg, by = key, suffixes = c(".a", ".b"))
cat("== B1 (", nrow(b1), "cells ) change from", lab[1], "to", lab[2], "==\n")
cat("v1.0.0 time:        ", summ(rel(b1$v1_sec.a, b1$v1_sec.b)), "\n")
cat("table time:         ", summ(rel(b1$table_sec.a, b1$table_sec.b)), "\n")
cat("speedup v1/table:   ", summ(rel(b1$speedup_vs_v1.a, b1$speedup_vs_v1.b)), "\n")
cat("same-code control (seq/v1):", lab[1], round(median(b1$seq_vs_v1.a), 3),
    "|", lab[2], round(median(b1$seq_vs_v1.b), 3), "\n")
big <- b1[abs(rel(b1$speedup_vs_v1.a, b1$speedup_vs_v1.b)) > 0.10,
          c(key, "speedup_vs_v1.a", "speedup_vs_v1.b")]
cat("cells whose speedup moved by more than 10%:", nrow(big), "\n")
if (nrow(big) > 0) print(transform(big, speedup_vs_v1.a = round(speedup_vs_v1.a, 2),
                                  speedup_vs_v1.b = round(speedup_vs_v1.b, 2)), row.names = FALSE)
cat("kappa range meeting no-regression gate:",
    lab[1], an[[1]]$kappa_noregress, "|", lab[2], an[[2]]$kappa_noregress, "\n")
same_be <- identical(an[[1]]$break_even$table_wins, an[[2]]$break_even$table_wins)
cat("break-even pattern (table wins per m, N) identical:", same_be, "\n")

# ---- B4 -------------------------------------------------------------------
b4 <- lapply(c(dir_a, dir_b), function(d) {
    x <- readRDS(file.path(d, "b4_fast-shortcut.rds"))
    s <- x[x$kernel == "sequential", c("m", "N", "threads", "median_sec")]
    t <- x[x$kernel == "table" & x$grain == -1, c("m", "N", "threads", "median_sec")]
    out <- merge(s, t, by = c("m", "N", "threads"), suffixes = c("_seq", "_tab"))
    out$ratio <- out$median_sec_seq / out$median_sec_tab
    out$identical <- all(x$identical_to_first)
    out
})
b4m <- merge(b4[[1]], b4[[2]], by = c("m", "N", "threads"), suffixes = c(".a", ".b"))
cat("\n== B4 (", nrow(b4m), "cells ) seq/table ratio at equal threads ==\n")
cat("ratio change:", summ(rel(b4m$ratio.a, b4m$ratio.b)), "\n")
cat("ratio range:", lab[1], paste(round(range(b4m$ratio.a), 2), collapse = "-"),
    "|", lab[2], paste(round(range(b4m$ratio.b), 2), collapse = "-"), "\n")
cat("all outputs identical:", lab[1], b4[[1]]$identical[1], "|", lab[2], b4[[2]]$identical[1], "\n")
print(transform(b4m[b4m$threads %in% c(1, 4, 12), c("m", "N", "threads", "ratio.a", "ratio.b")],
                ratio.a = round(ratio.a, 2), ratio.b = round(ratio.b, 2)), row.names = FALSE)

# ---- B2 -------------------------------------------------------------------
b2 <- lapply(c(dir_a, dir_b), function(d) {
    v <- readRDS(file.path(d, "b2_v1.0.0.rds"))
    f <- readRDS(file.path(d, "b2_fast-shortcut.rds"))
    out <- merge(v[, c("m", "criterion", "eval_sec", "objective_value")],
                 f[, c("m", "criterion", "eval_sec", "objective_value", "kernel_sec")],
                 by = c("m", "criterion"), suffixes = c("_v1", "_fast"))
    out$speedup <- out$eval_sec_v1 / out$eval_sec_fast
    out$kernel_share_fast <- out$kernel_sec / out$eval_sec_fast
    out$same_value <- out$objective_value_v1 == out$objective_value_fast
    out
})
b2m <- merge(b2[[1]], b2[[2]], by = c("m", "criterion"), suffixes = c(".a", ".b"))
cat("\n== B2 per-evaluation speedup ==\n")
print(data.frame(m = b2m$m, criterion = b2m$criterion,
                 v1_ms.a = round(1e3 * b2m$eval_sec_v1.a, 2), fast_ms.a = round(1e3 * b2m$eval_sec_fast.a, 2),
                 speedup.a = round(b2m$speedup.a, 2),
                 v1_ms.b = round(1e3 * b2m$eval_sec_v1.b, 2), fast_ms.b = round(1e3 * b2m$eval_sec_fast.b, 2),
                 speedup.b = round(b2m$speedup.b, 2),
                 share.b = round(b2m$kernel_share_fast.b, 2)), row.names = FALSE)
cat("objective values identical across builds:", all(b2m$same_value.a), all(b2m$same_value.b), "\n")

# ---- B5 -------------------------------------------------------------------
f5 <- c(file.path(dir_a, "b5_auto_kappa2.rds"), file.path(dir_a, "b5_auto.rds"))
f5a <- f5[file.exists(f5)][1]
b5 <- list(readRDS(f5a), readRDS(file.path(dir_b, "b5_auto.rds")))
b5m <- merge(b5[[1]], b5[[2]], by = c("m", "N", "power"), suffixes = c(".a", ".b"))
cat("\n== B5 auto check ==\n")
print(data.frame(m = b5m$m, N = b5m$N, power = b5m$power, picks = b5m$auto_kernel.b,
                 vs_seq.a = round(b5m$auto_vs_sequential.a, 3), vs_seq.b = round(b5m$auto_vs_sequential.b, 3),
                 vs_best.a = round(b5m$auto_vs_best.a, 2), vs_best.b = round(b5m$auto_vs_best.b, 2)),
      row.names = FALSE)
cat("max auto/sequential:", lab[1], round(max(b5m$auto_vs_sequential.a), 3),
    "|", lab[2], round(max(b5m$auto_vs_sequential.b), 3), "\n")
