# Table Kernel for `graph_shortcut()`: Local Review Report

*16 October 2026. Branch `fast-shortcut`. The kernel and tests are in commit `0344474` (built with κ = 1). The calibrated κ = 2, the switch-rule test, DESCRIPTION, NEWS and this report are in `abcf908`. The scripts and results in `dev/review/shortcut_table/` are not yet committed. Design: `debug/improve-shortcut.md` (git-ignored).*

*Machine: Windows 11 laptop, 10 physical / 12 logical cores, on AC power with sleep held off. R 4.6.1, Rtools45 g++ `-O2 -mfpmath=sse -msse2` (no FMA). v1.0.0 (`7884824`) and this branch were installed from source tarballs into separate libraries and timed in separate R processes.*

Tags as in the design: **[measured]** produced by running code here; **[derived]** follows from the algorithm or the measurements.

---

## 0. Summary

* **Same answers.** **[measured]**
  * The table kernel and the original kernel disagreed on **0 of 2.0 × 10⁹** rejection indicators: m = 2–16, dense and sparse graphs, three power levels.
  * v1.0.0 and this build produced identical rejection matrices in every benchmark cell.
  * Objective values from `create_obj_func()` are identical between the builds.
* **Faster wherever it is used.** **[measured]** Serial, per call:
  * 1.5–2.7× faster than v1.0.0 for m = 3–6 and 3.3–4.4× for m = 8–10 at high power;
  * 1.3–2.3× for m = 3–10 at low power;
  * up to 6.7× for m = 14–18 at large N.

  One GA objective evaluation (N = 5 × 10⁴) went from 3.7 ms to 1.7 ms at m = 5 (2.2×) and from 16.8 ms to 3.6 ms at m = 10 (4.7×).
* **Switch rule (O1).** **[measured]**
  * `"auto"` uses the table when 2^m ≤ 2N. That means up to m = 16 in the GA phase and up to m = 20 at N = 10⁶; above that it falls back to the original kernel.
  * With κ = 2, `"auto"` was never slower than the original kernel in any cell.
  * The cost is that in a few high-power cells near the boundary the table would have been 2.2–2.5× faster than what `"auto"` picks (§5).
* **Threads and grain size (O2).** **[measured]**
  * The gain holds at equal thread counts: 2.5–8.9× for m = 8–16 with 1–12 threads.
  * The exception is m = 16 with N = 5 × 10⁴, where the serial table build caps the gain at about 1.1× on 12 threads.
  * The table walk is insensitive to grain size, so the existing heuristic is kept.
* **Not measured here:** the curated problems, full optimisation runs (B3), the HPC nodes, and more than 12 threads.

---

## 1. Measurement Conditions

Two benchmark passes were spoiled and are kept only for the record.

| Result set | What happened | Status |
|---|---|---|
| `results/run1_update_noise/` | The laptop entered Modern Standby at 19:43 (idle timeout) and rebooted for an update at 20:07, mid-run. | **Discarded.** Same-code timings disagreed by up to ±40%. |
| `results/run2/b1_*_battery.*` | On battery from 20:44 to 21:03. | **Discarded.** |
| `results/run2/` (all other files) | On AC, no standby. A keep-awake request was held from 21:07. | **Used below.** |
| `results/run3/` | Full repeat on the committed build `abcf908`, 21:54–22:28. On AC, keep-awake held, and no Kernel-Power events in the window. | **Reproduces run 2** (below). |

**Check that the clean sets are clean.** v1.0.0 and this build's sequential kernel run the same code. Their B1 timings had a median ratio of 0.983 in run 2 (IQR 0.972–0.995, range 0.82–1.09 over 88 cells) and 1.003 in run 3 (IQR 0.991–1.015) **[measured]**.

**Reproducibility** (`compare_runs.R`, `results/run3/compare_run2_run3.txt`) **[measured]**. Run 3 changes no conclusion:

* B1 speedups over v1.0.0 moved by a median of −1.8% (IQR −3.6% to 0%). Only 2 of 88 cells moved by more than 10%, both at m = 5; for example 2.20× became 1.94× at high power with N = 5e4.
* The break-even pattern is identical, so the same range of κ meets the gate.
* B4 table/sequential ratios moved by a median of −0.1%. The one large move was at m = 8, N = 5e4, 12 threads (2.5× → 3.6×), on calls under 1 ms.
* B2 per-evaluation speedups were 2.2–4.5×, against 2.0–4.7× in run 2.
* `"auto"` made the same choice in all 16 B5 cells; worst auto/sequential was 1.022.
* Identity held again in every cell.

Read the single-run figures below as accurate to within about ±10%.

---

## 2. Gates

| Gate | Result |
|---|---|
| **G1.** Package tests | **Pass.** `test-graph_shortcut_table.R`: 16 tests, 504 expectations. `test-RcppExports.R` (gMCPLite parity), `test-objective_function.R`, `test-calc_power.R`, `test-post_optim_processing.R`, `test-optimisation.R` and `test-graph_optimal.R` all pass. The only two skips were already there (`apply_ctp` was removed earlier). Snapshots are unchanged; the CRLF-only rewrites were reverted. |
| **G2.** Sweep | **Pass.** 0 differences, 0 floating-point ties (§3). |
| **G3.** B3 identity on the curated problems | **Not run here**; left for the HPC suite. Local identity: v1.0.0 equals both kernels in all 176 B1 cells; sequential equals table in 240 B1 and calibration cells, all 252 B4 configurations and all 16 B5 cells. |
| **G4.** Speed | **Met at high power, not at low power.** See below. |

**G4 in detail.** The gate asks for at least 2× serial speedup at m ∈ {5, 6, 8} with N = 5 × 10⁴, and for `"auto"` never more than 5% slower than sequential.

* Speedup over sequential:
  * high power: **2.22×, 2.52×, 3.39×** (meets the gate);
  * low power: **1.65×, 1.66×, 1.95×**, with 0.9–1.2 rejections per trial (misses it).
  * The design's prototype profile corresponds to "high" here.
* `"auto"` vs sequential: never slower. Worst ratio 1.003 in direct timing (§5).

---

## 3. Equivalence (Design §6.2)

**Sweep** (`sweep.R`, output in `results/sweep.txt`) **[measured]**

* **Graphs:** for each m = 2–16, 50 dense graphs and 50 sparse ones (about half the edges and a third of the weights zeroed, and some rows of G summing to 0.7).
* **Trials:** N = 5 × 10⁴, equicorrelated (ρ = 0.5) one-sided p-values at three power levels.
  * Mean rejections per trial at m = 16 were 1.6–1.8 at low power, 5.7–6.1 at high power and 11.2–11.5 at very high power.
* **Result:** 2.025 × 10⁹ indicators compared between forced `"sequential"` and `"table"`, with **0 differences**. The 4-thread table output matched the serial output in every spot check.
* **Not covered:** the delivered graphs of the curated problems, which were not available.

**Builder check** (`builder_check.R`, `results/builder_check.txt`) **[measured]**

* For m = 2–10, 20 dense and 20 sparse graphs each, the table divided by α matches `calc_local_weights()` and `graphicalMCP::graph_generate_weights()` to a maximum absolute difference of **3.7 × 10⁻¹⁵**. On the critical-value scale the difference is at most 9 × 10⁻¹⁷.
* The full-set row is exactly `w * alpha`.

**Why the agreement is exact so often** **[derived]**

* The builder works in α-space with the sequential kernel's own expressions. So a set reached by rejecting in increasing index order, the common case, gets bit-identical critical values.
* A disagreement needs both of the following:
  * a rejection path out of index order;
  * a p-value within rounding (about 10⁻¹⁶ relative) of a critical value.

**Non-zero diagonals** **[measured]** The sequential kernel uses `g_kk` on the first rejection and the table ignores it. Callers cannot pass a non-zero diagonal:

* `calc_power_pvals()` rejects one through `is_graph_valid()`;
* the optimiser fixes the diagonal at 0;
* pruning forces it to 0.

---

## 4. B1: Serial Kernel Benchmark

`bench_kernel.R`, median of 20 iterations, 2 random dense graphs per m, two power levels **[measured]**.

**Speedup of the table kernel over v1.0.0, per call**

| m | high, N=1e4 | high, 5e4 | high, 2e5 | high, 1e6 | low, 1e4 | low, 5e4 | low, 2e5 | low, 1e6 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 1.5 | 1.6 | 1.6 | 1.5 | 1.8 | 1.4 | 1.4 | 1.3 |
| 4 | 1.8 | 1.9 | 1.8 | 1.8 | 1.4 | 1.5 | 1.4 | 1.4 |
| 5 | 2.2 | 2.2 | 2.3 | 2.2 | 1.6 | 1.7 | 1.6 | 1.6 |
| 6 | 2.7 | 2.7 | 2.5 | 2.6 | 1.7 | 1.7 | 1.6 | 1.6 |
| 8 | 3.3 | 3.4 | 3.5 | 3.4 | 2.0 | 2.0 | 2.1 | 2.0 |
| 10 | 3.9 | 4.2 | 4.4 | 4.3 | 2.1 | 2.1 | 2.2 | 2.3 |
| 12 | 4.1 | 5.2 | 5.3 | 5.2 | 1.9 | 2.6 | 2.5 | 2.6 |
| 14 | 3.2 | 5.6 | 6.3 | 6.7 | 1.2 | 2.5 | 3.2 | 3.3 |
| 16 | *1.0* | 3.4 | 5.5 | 6.6 | *0.3* | 1.2 | 2.4 | 3.2 |
| 18 | *0.3* | *1.5* | 3.8 | 6.6 | *0.1* | *0.5* | 1.4 | 3.0 |
| 20 | *0.1* | *0.4* | *1.5* | 4.2 | *0.0* | *0.1* | *0.5* | 1.6 |

Italics: cells where `"auto"` (κ = 2) uses the original kernel, so the user sees a speedup of about 1.0.

**Rejections per trial and build time**

* Rejections per trial: 1.2 → 7.2 at high power and 0.7 → 2.0 at low power, as m goes from 3 to 20.
* Table build alone: 0.02 ms at m = 8, 0.24 ms at m = 12, 1.3 ms at m = 14, 8.5 ms at m = 16, 41 ms at m = 18, 196 ms at m = 20. Below m = 8 the build is lost in the call overhead (about 6 µs).
* At N = 5 × 10⁴ the build is at most 9% of table time for m ≤ 12. It is 19–25% at m = 14, 53–65% at m = 16, and at least 80% from m = 18.
* As the design predicted, the speedup grows with rejections per trial: each rejection costs the original kernel an O(m²) graph update.

**Example times per call at N = 5 × 10⁴, high power (v1.0.0 → table):** m = 5: 3.61 → 1.64 ms; m = 8: 9.92 → 2.88 ms; m = 10: 15.85 → 3.82 ms; m = 12: 25.0 → 4.8 ms.

---

## 5. O1: The Switch Rule

**Rule:** use the table when m ≤ `table_max_m` (default 20) and **2^m ≤ κN, with κ = 2**. Kernel choice does not depend on `num_threads`, so output never depends on thread count.

**Break-even** (B1 plus a calibration pass at m = 13, 15, 17 and 19; `results/run2/b1_analysis.txt`) **[measured]**

* The table won in every graph at both power levels up to 2^m/N = **1.64**.
* Sequential first won at 2^m/N = **2.62**, at low power only. At high power the table still won up to about 5.
* So the break-even ratio depends on rejections per trial, which is not known before the call.
* Any κ from 1.64 up to 2.62 keeps `"auto"` no slower than sequential on this grid. κ = 2 is the round value in that range.

**Direct check at the boundary** (`bench_auto.R`, `results/run2/b5_auto_kappa2.*`) **[measured]**

| m | N | power | auto picks | auto (ms) | sequential (ms) | table (ms) |
|---:|---:|---|---|---:|---:|---:|
| 3 | 1e4 | high | table | 0.19 | 0.25 | 0.19 |
| 8 | 5e4 | high | table | 2.85 | 9.54 | 2.72 |
| 16 | 5e4 | low | table | 13.4 | 17.1 | 13.6 |
| 17 | 5e4 | high | sequential | 60.2 | 60.4 | **26.8** |
| 17 | 5e4 | low | sequential | 17.5 | 17.8 | 23.3 |
| 18 | 2e5 | low | table | 64.4 | 84.4 | 64.1 |
| 19 | 2e5 | high | sequential | 339 | 338 | **133** |
| 19 | 2e5 | low | sequential | 110 | 110 | 116 |
| 20 | 1e6 | low | table | 349 | 515 | 338 |

* `"auto"` times the same as the kernel it picks, so the option reads cost nothing measurable, and it is never slower than sequential.
* `bench::mark(check = TRUE)` confirmed identical output in every cell.
* The cost of a power-blind rule is in bold: the table would have been 2.2–2.5× faster in those two high-power cells.

**Recommendation for O1**

* Keep κ = 2. It meets the design's no-regression gate, and the curated problems (m ≤ 10) lie far inside the table region, where 2^m/N ≤ 0.02 in the GA phase.
* If m ≥ 16 becomes important, consider a **pilot rule**:
  * run the sequential kernel on the first few hundred trials, whose output is needed anyway;
  * estimate rejections per trial from them;
  * choose the kernel for the remaining trials with the cost model build ≈ c_b · 2^m m² against sequential ≈ c_s · N m² (1 + r). Here c_b ≈ 0.4–0.5 ns and c_s ≈ 0.6 ns (median, range 0.45–0.71); the design's prototype machine gave about 0.9 and 1 ns.

  The choice then depends on the data but still not on threads, and both kernels give the same output.
* **Recheck κ on the HPC nodes.** The break-even is set by the machine-specific ratio of build cost to per-trial cost **[derived]**.

---

## 6. B4: Threads and Grain Size (O2)

`bench_threads.R`, one dense graph, high power. Times are table / sequential in ms at equal thread counts **[measured]**.

| m | N | 1 thread | 4 threads | 12 threads |
|---:|---:|---|---|---|
| 8 | 5e4 | 2.55 / 10.9 (4.3×) | 0.87 / 4.15 (4.8×) | 0.86 / 2.14 (2.5×) |
| 8 | 1e6 | 55.2 / 224 (4.1×) | 20.2 / 87.7 (4.3×) | 13.1 / 50.3 (3.8×) |
| 12 | 5e4 | 4.96 / 29.3 (5.9×) | 1.98 / 11.4 (5.7×) | 1.26 / 6.28 (5.0×) |
| 12 | 1e6 | 92.1 / 585 (6.4×) | 35.5 / 240 (6.8×) | 21.1 / 135 (6.4×) |
| 16 | 5e4 | 16.3 / 59.3 (3.6×) | 12.0 / 22.7 (1.9×) | 11.5 / 13.0 (1.1×) |
| 16 | 1e6 | 157 / 1165 (7.4×) | 58.6 / 519 (8.9×) | 39.0 / 264 (6.8×) |

* **Scaling.** Both kernels gain about 4× from 1 to 12 threads on this laptop.
* **Where the gain shrinks.** At m = 16 and N = 5 × 10⁴, the serial table build (about 8–10 ms) does not shrink with threads, so the gain falls to 1.1× at 12 threads.
  * This is the cost of a thread-blind rule.
  * At higher thread counts on HPC nodes, m = 16 in the GA phase could be slightly slower with the table **[derived]**.
  * A parallel build would remove this (the subtrees under each first-removed hypothesis are independent). Not done here.
* **Grain size.** For the table walk, fixed grains of 64, 256, 1024, 4096 and 10,000 were all within noise of each other. The auto grain was within 3% of the best fixed grain at the median (worst 1.26×, a single outlier).

**Recommendation for O2:** keep the heuristic, with m² operations per trial and a target of 10⁵ per chunk (`TABLE_TARGET_OPS`).

---

## 7. B2-lite: One Objective Evaluation (O5)

`bench_objective.R`: synthetic, unconstrained, with a random dense candidate, N = 5 × 10⁴ and serial execution **[measured]**.

| m | rejections/trial | v1.0.0 eval (ms) | fast eval (ms) | speedup | kernel share, v1.0.0 → fast |
|---:|---:|---:|---:|---:|---|
| 5 | 2.0 | 3.72 | 1.68 | 2.2× | 90% → 82% |
| 6 | 2.6 | 5.28 | 2.05 | 2.6× | 91% → 85% |
| 8 | 3.5 | 9.96 | 2.84 | 3.5× | 99% → 85% |
| 10 | 4.0 | 16.76 | 3.57 | 4.7× | 97% → 91% |

* Figures are for the average-power criterion; the disjunctive criterion is within 10%. Objective values are identical between builds.
* After the change, the kernel is still most of an evaluation. The rest is decoding and constraint checks (about 0.2 ms) and `power_criterion()` (0.1–0.2 ms).

**O5:** the next bottleneck is still the kernel itself. The fixed R overhead (about 0.3–0.4 ms) would matter only for small m.

---

## 8. Recommendations for O1–O5

* **O1.** κ = 2, with the rule 2^m ≤ 2N (§5). Recheck on HPC; consider the pilot rule if m ≥ 16 matters.
* **O2.** Keep the grain heuristic (§6). Consider a parallel table build only if m ≈ 16 with many threads matters.
* **O3.** Keep `multigrain.shortcut_kernel` and `multigrain.shortcut_table_max_m`. The HPC A/B (B3) needs them, they cost nothing measurable, and they are undocumented.
* **O4.** Defer to the parametric/Simes design.
  * `multigrain:::graph_shortcut_table()` already exposes the table, but in α-space and in bitmask row order (bit j−1 is H_j).
  * A `calc_local_weights()`-format export would divide by α, or build with α = 1, and reorder the rows.
* **O5.** The kernel is still about 85% of an evaluation (§7). No new bottleneck has emerged.

## 9. Not Measured, and Caveats

* The curated problems, their delivered graphs, and full 5,000-generation runs (B3).
* The HPC node types, more than 12 threads, and other compilers or flags.
  * Builds with FMA contraction (for example `-march=native`) may change the last bits of the two kernels differently. Any non-identical B3 result should first be checked with `common.R::replay_trial()`.
* Memory: the table is 2^m · m · 8 bytes per call. That is 1.75 MB at m = 14, 8.4 MB at m = 16 and 168 MB at m = 20, allocated and freed on each call.

## 10. Files

* **Scripts**, in `dev/review/shortcut_table/`:
  * `common.R`, `sweep.R`, `builder_check.R`;
  * `bench_kernel.R`, `bench_threads.R`, `bench_objective.R`, `bench_auto.R`;
  * `analyse_b1.R`, `compare_runs.R`, `run_benchmarks.sh`.
* **Outputs**, in `results/`:
  * `sweep.*` and `builder_check.*`;
  * `run2/`, the clean benchmark set used in the tables;
  * `run3/`, the repeat on `abcf908`, with `compare_run2_run3.txt`;
  * `run1_update_noise/` and `run2/*_battery.*`, discarded and kept for the record.
