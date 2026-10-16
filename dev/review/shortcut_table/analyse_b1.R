# Analysis of B1 (serial kernel benchmark): identity across builds and
# kernels, speedups, break-even, calibration of TABLE_KAPPA, and the local
# B5/G4 checks.
#
# Usage (from the repo root):
#   Rscript dev/review/shortcut_table/analyse_b1.R <results_dir> [kappa]
# Reads b1_v1.0.0.rds, b1_fast-shortcut.rds and any
# b1_calibration*.rds (extra m values, fast-shortcut build only). If kappa is
# given, the B5/G4 checks use it; otherwise they use the calibrated value.

args <- commandArgs(TRUE)
dir <- args[1]
kappa_in <- if (length(args) >= 2) as.numeric(args[2]) else NA_real_

v1 <- readRDS(file.path(dir, "b1_v1.0.0.rds"))
fast <- readRDS(file.path(dir, "b1_fast-shortcut.rds"))
cal_files <- list.files(dir, pattern = "^b1_calibration.*[.]rds$", full.names = TRUE)
fast_all <- do.call(rbind, c(list(fast), lapply(cal_files, readRDS)))

key <- c("m", "N", "power", "graph")
wide <- function(d) {
    seq <- d[d$kernel == "sequential", c(key, "median_sec", "hash", "mean_rej")]
    tab <- d[d$kernel == "table", c(key, "median_sec", "hash", "build_sec")]
    names(seq)[5:6] <- c("seq_sec", "seq_hash")
    names(tab)[5:6] <- c("table_sec", "table_hash")
    merge(seq, tab, by = key)
}
fw <- wide(fast_all)

# ---- Identity -------------------------------------------------------------
idt <- merge(
    v1[, c(key, "median_sec", "hash")],
    wide(fast),
    by = key
)
names(idt)[names(idt) == "median_sec"] <- "v1_sec"
names(idt)[names(idt) == "hash"] <- "v1_hash"
cat("== Identity (rejection-matrix hashes) ==\n")
cat("cells:", nrow(idt),
    " v1.0.0 == fast sequential:", sum(idt$v1_hash == idt$seq_hash),
    " v1.0.0 == fast table:", sum(idt$v1_hash == idt$table_hash), "\n")
cat("calibration cells (fast only): seq == table:",
    sum(fw$seq_hash == fw$table_hash), "of", nrow(fw), "\n\n")

# ---- Speed table (median over graphs) ---------------------------------------
agg <- aggregate(
    cbind(v1_sec, seq_sec, table_sec, build_sec, mean_rej) ~ m + N + power,
    data = idt,
    FUN = median
)
agg$speedup_vs_v1 <- agg$v1_sec / agg$table_sec
agg$seq_vs_v1 <- agg$seq_sec / agg$v1_sec
agg$build_share <- agg$build_sec / agg$table_sec
agg <- agg[order(agg$power, agg$N, agg$m), ]
fmt <- transform(
    agg,
    v1_ms = round(1e3 * v1_sec, 2),
    seq_ms = round(1e3 * seq_sec, 2),
    table_ms = round(1e3 * table_sec, 2),
    build_ms = round(1e3 * build_sec, 3),
    rej = round(mean_rej, 2),
    speedup = round(speedup_vs_v1, 2),
    seq_over_v1 = round(seq_vs_v1, 3),
    build_share = round(build_share, 3)
)[, c("power", "N", "m", "rej", "v1_ms", "seq_ms", "table_ms", "build_ms",
      "build_share", "speedup", "seq_over_v1")]
cat("== B1: serial medians (ms), table speedup over v1.0.0 ==\n")
print(fmt, row.names = FALSE)

# ---- Break-even and kappa ---------------------------------------------------
# Two criteria for the "auto" rule 2^m <= kappa * N:
#   vs_seq:  auto never more than 5% slower than the sequential kernel
#            (design gate G4: no regression against v1.0.0);
#   vs_best: auto within 5% of the faster forced kernel (design B5).
fw$ratio <- 2^fw$m / fw$N
fw$table_wins <- fw$table_sec < fw$seq_sec
chosen_sec <- function(kappa) ifelse(fw$ratio <= kappa, fw$table_sec, fw$seq_sec)
slowdown_best <- function(kappa) chosen_sec(kappa) / pmin(fw$table_sec, fw$seq_sec) - 1
slowdown_seq <- function(kappa) chosen_sec(kappa) / fw$seq_sec - 1
kappas <- 2^seq(-6, 6, by = 0.25)
ks <- data.frame(
    kappa = kappas,
    max_vs_seq = vapply(kappas, function(k) max(slowdown_seq(k)), 0),
    max_vs_best = vapply(kappas, function(k) max(slowdown_best(k)), 0),
    total_sec = vapply(kappas, function(k) sum(chosen_sec(k)), 0)
)
ok <- ks[ks$max_vs_seq <= 0.05, ]
kappa_noregress <- if (nrow(ok) > 0) max(ok$kappa) else NA_real_
kappa_minimax <- ks$kappa[which.min(ks$max_vs_best)]
cat("\n== Break-even ==\n")
be <- aggregate(table_wins ~ m + N, data = fw, FUN = mean)
be <- be[order(be$N, be$m), ]
be$ratio <- 2^be$m / be$N
print(be[be$m >= 12, ], row.names = FALSE)
cat("\nLargest ratio 2^m/N where table wins in every graph/power:",
    format(max(fw$ratio[ave(fw$table_wins, fw$m, fw$N, FUN = all) == 1]), digits = 3), "\n")
cat("Smallest ratio where sequential wins in some graph/power:",
    format(min(fw$ratio[!fw$table_wins]), digits = 3), "\n")
cat("Largest kappa (2^(k/4) grid) with auto <= 5% slower than sequential everywhere:",
    format(kappa_noregress, digits = 4), "\n")
cat("Kappa minimising the worst slowdown vs the faster kernel:",
    format(kappa_minimax, digits = 4), "\n")
ks_print <- ks[ks$kappa >= 0.25 & ks$kappa <= 16, ]
ks_print[2:3] <- round(100 * ks_print[2:3], 1)
names(ks_print)[2:3] <- c("max_vs_seq_pct", "max_vs_best_pct")
print(ks_print, row.names = FALSE)

# ---- Local B5 / G4 ------------------------------------------------------------
kappa <- if (is.na(kappa_in)) kappa_noregress else kappa_in
cols <- c(key, "seq_sec", "table_sec", "ratio")
cat("\n== Local checks with kappa =", format(kappa, digits = 4), "==\n")
s_seq <- slowdown_seq(kappa)
s_best <- slowdown_best(kappa)
cat("G4-local: max slowdown of auto vs sequential:", round(100 * max(s_seq), 1), "%\n")
print(fw[order(-s_seq), ][1:3, cols], row.names = FALSE)
cat("B5-local: max slowdown of auto vs faster kernel:", round(100 * max(s_best), 1), "%\n")
print(fw[order(-s_best), ][1:5, cols], row.names = FALSE)
cat("\n== G4-local: serial table speedup over sequential, N = 5e4 ==\n")
g4 <- aggregate(cbind(seq_sec, table_sec) ~ m + power, data = fw[fw$N == 5e4, ], FUN = median)
g4$speedup <- round(g4$seq_sec / g4$table_sec, 2)
print(g4[g4$m %in% c(5, 6, 8), c("m", "power", "speedup")], row.names = FALSE)

saveRDS(list(identity = idt, agg = agg, kappa_grid = ks, kappa_noregress = kappa_noregress,
             kappa_minimax = kappa_minimax, break_even = be, fw = fw),
        file.path(dir, "b1_analysis.rds"))
