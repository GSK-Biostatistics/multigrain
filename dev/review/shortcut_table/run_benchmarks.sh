#!/usr/bin/env bash
# Runs the local benchmarks one after another (never concurrently, so
# timings are not contended). Run from the repo root:
#   bash dev/review/shortcut_table/run_benchmarks.sh <lib_v1.0.0> <lib_fast> <out_dir> [Rscript]
set -u
LIB_V1="$1"; LIB_FAST="$2"; OUT="$3"; RS="${4:-Rscript}"
D=dev/review/shortcut_table
mkdir -p "$OUT"
"$RS" $D/bench_kernel.R    "$LIB_V1"   "$OUT/b1_v1.0.0.rds"        > "$OUT/b1_v1.0.0.log" 2>&1
"$RS" $D/bench_kernel.R    "$LIB_FAST" "$OUT/b1_fast-shortcut.rds" > "$OUT/b1_fast-shortcut.log" 2>&1
"$RS" $D/bench_threads.R   "$LIB_FAST" "$OUT/b4_fast-shortcut.rds" > "$OUT/b4_fast-shortcut.log" 2>&1
"$RS" $D/bench_objective.R "$LIB_V1"   "$OUT/b2_v1.0.0.rds"        > "$OUT/b2_v1.0.0.log" 2>&1
"$RS" $D/bench_objective.R "$LIB_FAST" "$OUT/b2_fast-shortcut.rds" > "$OUT/b2_fast-shortcut.log" 2>&1
"$RS" $D/bench_kernel.R    "$LIB_FAST" "$OUT/b1_calibration.rds" 13,15,17,19 > "$OUT/b1_calibration.log" 2>&1
echo "all benchmarks done"
