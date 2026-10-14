# Phase D discovery report — joint optimisation of HSD spending parameters

Branch `gsd-build-spending` (HEAD `9b09385` when the batteries ran), 2026-10-14. R 4.3.3,
gsDesign 3.11.0, GA 3.2.5, nloptr 2.2.1, installed `multigrain` 0.3.0 (fixed path only; P1's
uncommitted files were not installed and were not used). `nsim = 1e5` throughout,
`num_threads = 2` in every optimiser and kernel call, at most two heavy R processes at once on
4 cores. Seeds: M0/M1 raw draws `set.seed(20261007)`; §2 design raw draws `set.seed(1)`; outer
optimiser runs `set.seed(s)`, s ∈ {1, 2, 3} (§2 design {1, 2}, plus s = 3 at ϑ = −4); d04 staged
runs `set.seed(1)`, `set.seed(2)`; GA bank subsample `set.seed(11)`; GA runs
`set.seed(100 + s)`, s = 1..5.

Phase D replaces the design's derived cost and accuracy numbers (`dev/hsd_opt_design_revised.md`
§4.1–4.3) with measurements, and settles by the build plan's §3 rules the defaults that P3 needs:
the lattice step, whether the warm-started profile is dense by default, the GA population
constant, whether the warm-started profile is good enough to stand in for an outer loop of full
optimisations (the gate the local stage rests on), whether building the eager tables needs a
progress step, and the COBYLA budget at profile nodes. Every number below is tagged
**[measured]**; the design's **[derived]** values are kept beside them where they exist.

Designs (brief D-b): **M0/M1** — m = 2, powers (0.90, 0.80), corr 0.5, t = (0.5, 1),
gain `trial_success_gsd(0.6 * d(t1) + 0.4 * d(t2), d = c(1, 0.8))`, M1 adds
`stopping_rule_gsd(r1)`, `graph_constraint_free(2)`, H1 `sfLDOF`, H2 `sfHSD(γ)`.
**§2 design (S2)** — m = 3, powers (0.9, 0.85, 0.8), independent, t = (0.5, 1), gain
`trial_success_gsd(r1 + r2 + r3)`, `graph_constraint_free(3)`, H1 `sfLDOF`, H2 and H3 one shared
`sfHSD(γ)`. Lattice `seq(-4, 1, 0.25)` (21 nodes).

---

## Section 1 — Where the transform's time goes

Source: `d01_transform_split_output.txt` (D1, D2) and `d03_memory_and_copy_output.txt` (D5).

**Table D1a — `.gsd_boundary_table()` per hypothesis, HSD(−2), median of 5** [measured]

| G | K = 2: s per table | ms per `gsBound1()` | K = 3: s per table | ms per `gsBound1()` |
|---:|---:|---:|---:|---:|
| 128 | 0.109 | 0.852 | 0.145 | 1.133 |
| 256 | 0.200 | 0.781 | 0.309 | 1.207 |
| 512 | 0.428 | 0.836 | 0.613 | 1.197 |
| 1024 | 0.887 | 0.866 | 1.222 | 1.193 |

Linear fit 0.877 ms × G (K = 2) and 1.198 ms × G (K = 3), intercept ≈ 0 [measured]; design
§4.1: 1.06 ms per call, 1.1 s / 1.5 s per G = 1024 table [derived].

**Table D1b — `.gsd_invert()` on G = 1024 tables, median of 11** [measured]

| N | K | ms per hypothesis | ms per look | m = 2 | m = 3 |
|---:|---:|---:|---:|---:|---:|
| 1e4 | 2 | 3.40 | 1.70 | 6.8 ms | 10.2 ms |
| 1e5 | 2 | 7.73 | 3.87 | 15.5 ms | 23.2 ms |
| 1e4 | 3 | 3.08 | 1.03 | 6.1 ms | 9.2 ms |
| 1e5 | 3 | 11.64 | 3.88 | 23.3 ms | 34.9 ms |

Full `transform_pvalues_gsd()`, m = 3, K = 2, N = 1e5, G = 1024: **2.85 s** [measured] (3 tables
+ 3 inversions predict 2.68 s; job prompt 3.2 s). The transform is ≈ 99 % tables; inverting one
hypothesis on 1e5 rows costs 8–12 ms [measured] against the 50 ms per profile node assumed in
§4.2 [derived]. In d02/d04/d05 a full per-node transform cost a median 1.9 s (m = 2) and 2.8 s
(m = 3) [measured].

**Table D2 — interpolation error vs G = 1024 (and vs the exact boundary), K = 2, t = (0.5, 1),
200 log-uniform p in [1e-6, 0.025], worst case over ϑ ∈ {−4, −2, 1} and look** [measured]

| G | max rel. error vs 1024 | max rel. error vs exact |
|---:|---:|---:|
| 128 | 1.38e-4 | 1.39e-4 |
| 256 | 3.44e-5 | 3.51e-5 |
| 512 | 8.18e-6 | 8.33e-6 |
| 1024 | 0 | 1.94e-6 |

The error sits at the final look (interim ≈ 1e-14) and is largest at ϑ = 1. This is the number
roxygen quotes for `grid_size = 512` (8.2e-6).

**Table D5 — column copy vs kernel, N_GA = 5e4, m = 3, K = 2, 2 threads** (d03) [measured]

| Operation | median ms | p90 ms |
|---|---:|---:|
| assign K = 2 columns (top level) | 0.145 | 0.183 |
| kernel call | 1.767 | 2.527 |
| closure eval (`<<-` assign + kernel) | 1.946 | 2.397 |

Copy / kernel = **8.2 %** (threshold 15 %): no C++ column copier [measured]; design §4.2 ≈ 0.3 ms
against 3 ms [derived]. Projected 16.4 % at m_act = 2 and 32.8 % at m_act = 4. `tracemem()`: the
kernel call never copies; `<<-` on a matrix built inside the factory never copies (one copy in
total if the matrix came in as a factory argument); a plain `<-` inside the closure copies on
every evaluation. Bank at defaults (M = 21, K = 2, N_GA = 5e4, m_act = 2) = 33.6 MB.

**Eager-table cost at the plan's sizes** [derived from D1a]: one schedule, K = 2, M = 21,
G = 1024: 21 × 0.887 s = **18.6 s** (+ 0.9 s for the fixed `sfLDOF` table = 19.5 s); K = 3:
25.7 s; two schedules K = 2: 37.3 s; two schedules K = 3: 51.3 s. Design §4.1: 23 s / 32 s
[derived].

## Section 2 — How many evaluations the optimiser makes

**Table D6 — §2 design, `sfLDOF` × 3, drivers in `graph_optimise_gsd()` order with a counting
gain** [measured] (`d04_eval_counts_and_budget_output.txt`)

| Run | Stage | GA generations / COBYLA iterations | evaluations | wall (s) |
|---|---|---:|---:|---:|
| seed 1 | GA | 366 | 59 871 | 184.5 |
| | local (COBYLA) | 112 | 116 | 0.73 |
| | prune | – | 10 | 0.49 |
| | final | – | 1 | 0.04 |
| seed 2 | GA | 499 | 81 864 | 227.0 |
| | local | 94 | 98 | 0.54 |
| | prune | – | 10 | 0.31 |
| | final | – | 1 | 0.03 |
| `global_search = FALSE` (Holm start) | local | 96 | 100 | 0.58 |
| | prune | – | 10 | 0.24 |
| | final | – | 1 | 0.02 |

Staged totals 185.7 s (seed 1), 227.8 s (seed 2), 0.8 s (no GA). The one-call
`graph_optimise_gsd()` at seed 1 took **172.7 s** and is `identical()` to the staged run on
`hyp_weight` and `trans_matrix` (gain 2.448070); the same holds for `global_search = FALSE`
(0.9 s, gain 2.448120). The seed-2 one-call check was cut for wall time. GA counts include the
GA's full-sample re-score (one call) and every Nelder–Mead evaluation. The GA is **99–100 %** of
the run [measured] (design §4.1: ≈ 90 % [derived]); 2.8–3.1 ms per GA evaluation on 5e4 rows
including GA overhead [measured] (design: ≈ 3 ms [derived]). Each GA evaluation count is
≥ 200 × popSize = 4e4 as the design reasoned; generations are 366–499, not ≈ 200, because the
`run = 200` window keeps restarting on small improvements. The full outer runs of d02 on this
design took 126–276 s (GA generations 262–547) [measured], against README 4.5's 97–154 s for
its (different-gain) §2 design.

**Table D7 — warm-started COBYLA `E_C` on the §2 design (H1 `sfLDOF`, H2/H3 HSD(γ)), full
sample, from the γ = −2 optimum of the d02 profile** [measured]

| node γ | evaluations | COBYLA wall (s) | ms / eval | gain at start | gain after |
|---:|---:|---:|---:|---:|---:|
| −2 (own optimum) | 95 | 0.65 | 6.9 | 2.415870 | 2.415870 |
| −2.25 | 107 | 0.61 | 5.7 | 2.419680 | 2.419810 |
| −1.75 | 124 | 0.74 | 5.9 | 2.410530 | 2.411170 |
| −3 | 103 | 0.65 | 6.3 | 2.430290 | 2.431170 |
| −1 | 104 | 0.62 | 6.0 | 2.391980 | 2.392550 |
| −4 | 93 | 0.44 | 4.8 | 2.439440 | 2.439510 |
| 1 | 127 | 0.79 | 6.2 | 2.311530 | 2.314170 |

`E_C` = **93–127, median 105** at 5 graph genes [measured] against 200–800 [derived, §4.2]; at
1 graph gene (M0, M1) the d02 profile needed 23–31 evaluations per node. One profile node then
costs ≈ 15 ms inversion (2 active hypotheses) + 3 × ≈ 6 ms candidate kernel calls + 0.45–0.8 s
COBYLA ≈ **0.5–0.8 s** [measured parts] against 1.5–5.5 s [derived]; a 21-node d = 1 profile
≈ **11–17 s** on the §2 design and ≈ 2.5 s on M0/M1 [measured: d02 refinement totals 12.6 s
and 2.0 s, transform excluded] against 0.5–2 min [derived].

## Section 3 — Well-ordering across [−10, 2]

**Table D3** (d01 (b), G = 1024) [measured]: every node of `seq(-10, 2, 0.25)` (49 nodes) on
the four schedules (0.5, 1), (1/3, 2/3, 1), (0.25, 0.5, 0.75, 1), (0.6, 1) **passes** both
well-ordering and invertibility; the minimum number of distinct positive boundaries in any column
is 1024 (the `.hsd_check_invertible()` threshold is 2). Corners ϑ = −10, t = (0.1, 1) and ϑ = 2,
t = (0.9, 1) pass, with smallest boundaries 7.8e-19 and 7.4e-16; no boundary ≤ 0 anywhere. 198
tables in 118.8 s on 2 cores. Identities: `sfHSD(0.025, c(0.5, 1), 0)$spend` = (0.0125, 0.025)
(`param = 0` accepted); `max |sfHSD(1, t, g) / closed form − 1| = 0`; proportionality in `a`
between 1e-14 and 0.025 holds to 3.3e-16. **No node in [−10, 2] fails**; the plan's d01 pause
condition does not fire.

| schedule | min distinct | min boundary | min positive boundary |
|---|---:|---:|---:|
| (0.1, 1) corner | 1024 | 7.80e-19 | 7.80e-19 |
| (0.25, 0.5, 0.75, 1) | 1024 | 5.08e-18 | 5.08e-18 |
| (0.5, 1) | 1024 | 6.69e-17 | 6.69e-17 |
| (0.6, 1) | 1024 | 1.83e-16 | 1.83e-16 |
| (0.9, 1) corner | 1024 | 7.43e-16 | 7.43e-16 |
| (1/3, 2/3, 1) | 1024 | 1.23e-17 | 1.23e-17 |

## Section 4 — Warm-started profile versus the outer loop

Prototype (d02): one `transform_pvalues_gsd()` per node with `hsd_fn(γ)`; outer reference = a
full default `graph_optimise_gsd()` at γ ∈ {−4, −3, −2, −1, 0, 1}; profile = COBYLA
(`multigrain_control()` local options) on the full sample at each node, started from the better
of the previous node's optimum and the overall best graph (the refined optimum at the start
node), encoded by `closest_graph_to_constraints()` + `create_start_params(sum_to_one_constraint
= FALSE)`, decoded by `param_to_solution(process = TRUE)` + `repair_graph()`; P̃ = best of the
candidates and the refined graph (so P̃ ≥ S). Start node = node of the best seed-1 outer result.
No GA-archive start (the design's third candidate) was used, so this is a conservative test.

**Table D4 — M0** (start node −3) [measured]

| node | P̃ | S (fixed graph) | Û seed 1 | Û seed 2 | Û seed 3 | max Û | gap = max Û − P̃ | in plateau |
|---:|---:|---:|---:|---:|---:|---:|---:|:---:|
| −4 | 0.710446 | 0.710442 | 0.710438 | 0.710446 | 0.710446 | 0.710446 | +0.0 | no |
| −3 | 0.711060 | 0.711060 | 0.711060 | 0.711081 | 0.711060 | 0.711081 | +2.08e-5 | yes |
| −2 | 0.710876 | 0.710781 | 0.710895 | 0.710883 | 0.710868 | 0.710895 | +1.88e-5 | yes |
| −1 | 0.709550 | 0.709211 | 0.709570 | 0.709567 | 0.709570 | 0.709570 | +1.92e-5 | no |
| 0 | 0.706727 | 0.706069 | 0.706752 | 0.706737 | 0.706718 | 0.706752 | +2.56e-5 | no |
| 1 | 0.701950 | 0.700924 | 0.701950 | 0.701935 | 0.701953 | 0.701953 | +2.4e-6 | no |

argmax P̃ = −2.75 (0.711126); plateau(P̃, 5e-4) = [−3.75, −1.75]; plateau(S) = [−3.75, −1.75].
Max gap inside the plateau **2.08e-5**, outside 2.56e-5. **PASS.**

**Table D4 — M1** (start node 0) [measured]

| node | P̃ | S | Û seed 1 | Û seed 2 | Û seed 3 | max Û | gap | in plateau |
|---:|---:|---:|---:|---:|---:|---:|---:|:---:|
| −4 | 0.680199 | 0.680199 | 0.680144 | 0.680167 | 0.680203 | 0.680203 | +4.0e-6 | no |
| −3 | 0.685203 | 0.685183 | 0.685224 | 0.685146 | 0.685268 | 0.685268 | +6.48e-5 | no |
| −2 | 0.689636 | 0.689625 | 0.689614 | 0.689636 | 0.689636 | 0.689636 | 0 | no |
| −1 | 0.692275 | 0.692223 | 0.692275 | 0.692265 | 0.692269 | 0.692275 | 0 | no |
| 0 | 0.692769 | 0.692769 | 0.692769 | 0.692731 | 0.692792 | 0.692792 | +2.32e-5 | yes |
| 1 | 0.690615 | 0.690490 | 0.690487 | 0.690609 | 0.690610 | 0.690610 | −4.4e-6 | no |

argmax P̃ = −0.25 (0.693038); plateau(P̃) = [−0.75, 0]; plateau(S) = [−0.75, 0]. Max gap inside
**2.32e-5**, outside 6.48e-5. **PASS.**

**Table D4 — §2 design** (start node −4; outer seeds 1, 2, plus seed 3 at −4 only) [measured]

| node | P̃ | S | Û seed 1 | Û seed 2 | Û seed 3 | max Û | gap | in plateau |
|---:|---:|---:|---:|---:|---:|---:|---:|:---:|
| −4 | 2.439860 | 2.439860 | 2.439860 | 2.440100 | 2.440230 | 2.440230 | +3.70e-4 | yes |
| −3 | 2.430660 | 2.430540 | 2.431300 | 2.431330 | – | 2.431330 | +6.70e-4 | no |
| −2 | 2.415870 | 2.415720 | 2.416620 | 2.416510 | – | 2.416620 | +7.50e-4 | no |
| −1 | 2.392920 | 2.391920 | 2.392650 | 2.392850 | – | 2.392850 | −7.0e-5 | no |
| 0 | 2.359420 | 2.358500 | 2.359540 | 2.359720 | – | 2.359720 | +3.00e-4 | no |
| 1 | 2.314060 | 2.311610 | 2.314450 | 2.314420 | – | 2.314450 | +3.90e-4 | no |

argmax P̃ = −4 (the lower bound; 2.439860); plateau(P̃) = plateau(S) = [−4, −4]. Max gap inside
**3.70e-4** (one node), outside up to **7.5e-4**. **PASS, with a thin margin** — see below.

Reading of the §2 result. The gain here is a sum of three rejection indicators (scale 0–3), and
the full optimiser's own seed-to-seed spread at one node is 2.4e-4–3.7e-4 (−4: seeds 1, 2, 3 =
2.43986, 2.44010, 2.44023) [measured] — ten times the 4e-5 of README 4.2. The gap at −4 is
entirely that spread (P̃(−4) is the seed-1 outer value, which COBYLA did not improve), not a
failure of warm-starting. Outside the plateau the profile is 3–8e-4 below the best outer run at
−3, −2, 0, 1, so the warm-started profile is a lower bound tighter than S but not tight to 5e-4
on this 5-gene design away from its argmax; the GA-archive start that P3 adds is the design's
remedy for exactly this, and was not part of the prototype.

**README 7.5 reproduction (M1)** [measured]: P̃ at γ = −4 (LDOF-like H2) 0.680199 with
w1 = 0.439; profile argmax γ = −0.25 (between LDOF-like and Pocock-like; HSD(1) ≈ Pocock) with
P̃ = 0.693038 and w1 = 0.453; gain +0.0128 (+1.9 %). README 7.5 (nsim 2e5, `sfLDOF` vs
`sfLDPocock` on H2, w1 grid): 0.6707 at w1 0.32 → 0.6903 at w1 0.46 (+0.0196). Under M1 the
profile does move H2 toward early spending, but its maximum is interior (−0.25), not at the
Pocock end (P̃(1) = 0.690615, 2.4e-3 below); w1 moves only 0.44 → 0.45 here because the HSD(−4)
optimum already sits at w1 ≈ 0.44 (the w1 plateau under M1 is wide: outer runs at γ = −4 give
w1 0.398–0.444 within 6e-5). M0 for contrast: argmax −2.75 (0.711126, w1 0.956), gain falls
steadily toward γ = 1 (0.701950) — early spending on H2 is worth nothing without the stopping rule,
as README 7.5 found.

**Cost per node** [measured]

| design | 21 full transforms (prototype) | refinement, 21 nodes | evaluations per node (median) | ms per COBYLA evaluation | outer loop (6 nodes × seeds) |
|---|---:|---:|---:|---:|---:|
| M0 | 39.9 s (1.91 s / node) | 2.0 s | 26 | 2.7 | 18 runs, 1183 s (median 64 s) |
| M1 | 37.0 s (1.76 s / node) | 2.0 s | 27 | 2.9 | 18 runs, 1177 s (median 63 s) |
| S2 | 59.0 s (2.81 s / node) | 12.6 s | 102 | 5.5 | 13 runs, 2430 s (median 186 s) |

In the package path the transform column becomes one inversion per node (15–23 ms, D1b), since
the spec's tables are built once; the profile is then ≈ 2–3 s (m = 2) and ≈ 13 s (m = 3).

## Section 5 — GA with a parameter gene

Prototype bank closure (d04 `bank_*`): 21 node arrays on `idx_ga = sample.int(1e5, 5e4)`
(seed 11); genes = 1 graph gene (`graph_constraint_free(2)`) + `u ∈ [0, 1]`,
`q = min(21, 1 + floor(u · 21))`; H2's two columns copied from slice q into a 5e4 × 4 working
matrix built inside the factory; decode/penalty block copied from `create_obj_func_gsd()`.
`GA::ga()` with the fixed driver's arguments (`.cauchy_population`,
`.make_cauchy_mutation_multi(0.1, 1)`, `optim = TRUE`, default `optimArgs`, `popSize 200`,
`run 200`, `maxiter 1e5`, `pcrossover 0.2`, `pmutation 0.8`), suggestions = the two
`.build_start_matrix()` rows with `u = (9 − 0.5)/21` (start node −2) appended. Without the gene:
`.graph_optimise_ga_gsd()` at γ = −2. Seeds `set.seed(100 + s)`. "full gain" = the winner
re-scored on the full 1e5 rows at its own node. Building the bank: 21 transforms ≈ 39 s
(prototype; tables would be prebuilt), slicing 0.13 s; H1's columns identical across nodes.

**Table D8 — GA with vs without the parameter gene, 5 seeds (popSize 300: 3 seeds)** [measured]

| design | configuration | generations (median; range) | evaluations (median) | wall s (median) | best node found (all seeds) | best fitness, 5e4 subsample | full-sample gain at own node | d02 profile argmax (P̃) |
|---|---|---:|---:|---:|---|---:|---:|---|
| M0 | gene, popSize 200 | 205 (203–251) | 36 324 | 68.6 | −2.75 × 5 | 0.711028 | 0.711128 | **−2.75** (0.711126) |
| M0 | gene, popSize 300 | 204 (203–224) | 52 971 | 95.9 | −2.75 × 3 | 0.711028 | 0.711129 | |
| M0 | no gene (γ = −2), popSize 200 | 227 (216–252) | 39 799 | 61.7 | (−2 fixed) | 0.710638–0.711586 | 0.710846–0.710891 | |
| M1 | gene, popSize 200 | 215 (205–223) | 37 575 | 70.9 | −0.75 × 5 | 0.692647 | 0.692738 | **−0.25** (0.693038), plateau [−0.75, 0] |
| M1 | gene, popSize 300 | 207 (201–210) | 53 931 | 100.8 | −0.75 × 3 | 0.692647 | 0.692738 | |
| M1 | no gene (γ = −2), popSize 200 | 220 (210–380) | 38 414 | 60.4 | (−2 fixed) | 0.689040–0.690414 | 0.689474–0.689619 | |

- **Generations with / without the gene at popSize 200 (medians): 0.90 (M0), 0.98 (M1)**
  [measured]; rule threshold 1.5. Evaluations likewise ≈ 0.9–1.0×; wall time ≈ 1.1–1.2× (the
  column copy and the decode of one more gene: 68.6 vs 61.7 s on M0, 70.9 vs 60.4 s on M1).
  Design §3.10/§4.2: 1.2–1.5× the GA time [derived].
- **Node found.** On M0 every gene run found the profile argmax −2.75. On M1 every gene run
  found −0.75, the best node on the 5e4 GA subsample; on the full sample it is the plateau's
  lower end (P̃(−0.75) = 0.692730, 3.0e-4 below the argmax −0.25) — inside the plateau, and the
  Stage L sweep that starts there reaches −0.25. Every gene run visited all 21 nodes, with the
  winner visited ≈ 2.5–3× as often as its neighbours (archive coverage complete).
- **Population.** popSize 300 found the same nodes, the same fitness to 1e-6, in the same number
  of generations and ≈ 1.45× the evaluations and wall time: more population does not help the
  gene at d = 1. The popSize rule `min(max(40(m + d), 200), 500)` gives 200 for m + d = 3
  (inert here, as the plan noted).
- Against the plain GA at the start node, the gene GA improves the full-sample gain by
  +2.4e-4 (M0) and +3.1e-3 (M1) — the joint search moves off a fixed γ = −2 as intended.

## Section 6 — Lattice step

d05: at the d02 argmax g*, fine grid `seq(g* − 0.5, g* + 0.5, 0.025)` (41 nodes), one transform
and a COBYLA from the g* optimum at each, full sample. Coarse nodes = fine nodes at multiples of
Δ; the window is intersected with the default bounds [−4, 1] because no lattice step can reach
nodes outside them (the unbounded window is reported beside it).

**Table D9 — loss of a coarse step = F − best coarse node in the window** [measured]

| design | g* | F (fine max in bounds) | Δ | coarse nodes | best coarse node | loss | loss, unbounded window |
|---|---:|---:|---:|---:|---:|---:|---:|
| M1 | −0.25 | 0.693060 at −0.200 | 0.5 | 2 | −0.50 | **1.31e-4** | 1.31e-4 |
| M1 | | | 0.25 | 5 | −0.25 | **2.12e-5** | 2.12e-5 |
| M1 | | | 0.125 | 9 | −0.25 | **2.12e-5** | 2.12e-5 |
| S2 | −4 | 2.439860 at −4.000 | 0.5 | 2 | −4 | 0 | 0 |
| S2 | | | 0.25 | 3 | −4 | 0 | 0 |
| S2 | | | 0.125 | 5 | −4 | 0 | 0 |

The §2 design is uninformative for Δ: its argmax is the lower bound, which is a node of every
lattice. Its fine grid shows the gain still rising below the bound (2.442500 at γ = −4.5,
+2.6e-3 over −4) — the default `bounds = c(-4, 1)` cuts off this design's optimum (see
section 10). On M1 the fine curve is ragged at the 1e-5–5e-5 scale between neighbouring fine
nodes (Monte Carlo and COBYLA noise), so the 2.1e-5 loss at Δ ≤ 0.25 is at the noise floor.
Cost: M1 74 s transforms + 4.6 s refinement; S2 109 s + 24 s.

---

## Decisions

1. **Default `step` = 0.25.** Rule: coarsest Δ ∈ {0.5, 0.25, 0.125} with loss ≤ 4e-5 on both
   d05 designs, else 0.25. Measured (D9): M1 loss 1.31e-4 at Δ = 0.5 (fails), **2.12e-5 at
   Δ = 0.25** and 2.12e-5 at 0.125; §2 design 0 at every Δ (argmax on the bound, uninformative).
   Coarsest passing step 0.25 → `sf_hsd_opt(step = 0.25)` stands; lattice M = 21 on [−4, 1].
2. **Profile density: full 21-node sweep by default** (no every-second-node thinning; no
   `control` opt-in needed). Rule: 21 × (inversion + 3 kernel calls + E_C × per-eval) ≤ GA stage
   time on the same design. Measured on the §2 design: 21 × (15–23 ms inversion [D1b] + 3 × ≈ 6
   ms + E_C 93–127 [D7] × 4.8–6.9 ms) ≈ **11–17 s** (prototype refinement total 12.6 s, d02)
   against a GA stage of **184.5–227.0 s** (D6) — ≈ 6–8 % of it. On M0/M1 ≈ 2–3 s against a GA
   of ≈ 60–70 s.
3. **popSize constant: 40 stands** (`min(max(40 (m + d), 200), 500)`). Rule: change to 60 only if
   the GA with the gene needs > 1.5× the generations without it at equal popSize. Measured (D8):
   0.90 (M0) and 0.98 (M1). popSize 300 gave identical answers at 1.45× the cost.
4. **Warm-start gate: PASS → Stage L as designed.** Rule: max over d02 nodes inside the plateau
   of (max-seed full gain − profile gain) ≤ 5e-4 on all three designs. Measured (D4): **M0
   2.08e-5, M1 2.32e-5, §2 design 3.70e-4** (one outer node inside each of the M1 and §2
   plateaus, two for M0). The §2 margin is thin and equals the full optimiser's own seed spread
   at that node; outside the plateau the §2 gap reaches 7.5e-4 (M0 2.6e-5, M1 6.5e-5). The
   §4.3-2 fallback is not triggered; the GA-archive third start (absent from the prototype)
   should narrow the off-plateau gap.
5. **Table-building progress step: yes, a `cli::cli_progress_step()` per node gated on
   verbosity** (default `grid_size` stays 1024; roxygen quotes the D2 error for
   `grid_size = 512`: max relative error 8.2e-6). Rule: eager-table time at m = 3, K = 2, M = 21,
   G = 1024 > 25 % of the d04 default run (also checking K = 3 and a two-schedule spec). Measured
   inputs: tables 19.5 s at the primary size (18.6 s tagged + 0.9 s fixed; D1a) = **11 %** of the
   172.7 s default run — below; K = 3 one schedule 25.7 s = 15 % — below; two schedules K = 2
   37.3 s = 22 % — below; **two schedules K = 3 51.3 s = 30 % — above.** The primary case does
   not fire the rule; the rule's own extra check does, so the progress step is recommended (it
   costs nothing on the fixed path). If the orchestrator reads the rule as primary-case-only,
   the answer is "no"; the numbers are here either way. (The K = 3 comparison is against the
   K = 2 run; a K = 3 default run, not measured, would be longer and the share lower.)
6. **COBYLA budget: same `local_opt` at profile nodes** (`NLOPT_LN_COBYLA`, `xtol_rel 5e-8`,
   `xtol_abs 5e-9`, `maxeval 5000`); no profile-specific `maxeval`. Measured: 23–35 evaluations
   per node at 1 graph gene, 83–128 at 5 graph genes (d02, D7, d05); the maximum is 2.6 % of
   `maxeval`, so every node stopped on the tolerances, not the budget.

The five defaults for P3: **step 0.25; full 21-node profile; popSize constant 40 (rule
unchanged); Stage L as designed (gate passed); progress step per node on table building
(verbosity-gated)**; plus the unchanged COBYLA options.

## Pre-existing issues (reported, not fixed)

1. **`start_graph` ignored when `global_search = FALSE`** (plan reconciliation 9).
   `.build_start_matrix()` always puts the Holm seed in row 1 (`R/optimisation_start.R:9`,
   `seeds[[1]] <- create_start_params(graph_constraint)`), the fixed-sequence seed in row 2 and
   the user's graphs after them; with `global_search = FALSE` the drivers take
   `.build_start_matrix(graph_constraint, start_graph)[1, , drop = TRUE]` as COBYLA's start
   (`R/optimisation_gsd.R:229-233`, `R/optimisation.R:188-192`), i.e. always Holm. A user-supplied
   `start_graph` therefore has no effect on the local-only search, in both the fixed-sample and the
   group sequential optimiser. Not fixed; the opt path encodes `start_graph[[1]]` itself
   (`.hsd_encode_graph()`).
2. `pvals[sample(nsim), ]` in `R/optimisation.R` permutes the first `nsim` rows rather than
   sampling `nsim` of all rows (from `dev/implement_gsd.md`; not replicated in `_gsd` code, whose
   `.sample_pvals_gsd()` uses `sample.int(pvals$nsim, nsim)`).
3. Unary minus in `parse_and_transform()` (`R/trial_success.R`) indexes a second argument that
   does not exist (from `dev/implement_gsd.md`).

## What the design did not anticipate

- **Optimiser spread on the §2 design is ≈ 4e-4, not 4e-5.** With the additive gain
  `r1 + r2 + r3` (scale 0–3) the full default optimiser's seed-to-seed spread at one γ is
  2.4e-4–3.7e-4 [measured], so the absolute 5e-4 plateau tolerance and gate are only just above
  the optimiser's own noise on this design; the 4e-5 of README 4.2 was measured on an m = 2
  gain on a 0–1 scale. The d02 gate passes there by 1.3e-4.
- **The §2 optimum lies outside the default bounds.** P̃ is monotone decreasing in γ over
  [−4, 1] and the fine grid rises to 2.4425 at γ = −4.5. The default `bounds = c(-4, 1)` makes
  the reported optimum a boundary value for this undiscounted design; the print method should
  probably say when ϑ* sits on a bound (P4 decision, not taken here).
- **The default run is slower than README 4.5's numbers.** One-call default run on the §2
  design 172.7 s; staged 185.7 s and 227.8 s; d02 outer runs 126–276 s; GA generations 366–547,
  well above the ≈ 200 the design's evaluation estimate assumed. The design's "today 2–2.5 min"
  is about 3–3.8 min here. Projected new-path total (d = 1, m = 3): GA 1.0–1.5 × 3–3.8 min
  (D8: 0.90–0.98× the generations, ≈ 1.1–1.2× the wall time with the gene at m = 2) + 20 s tables + ≈ 15 s profile ≈ 3.5–6
  min — inside the §4.2 projection of 3.5–6.5 min, so the plan's "> 2×" pause does not fire.
- **`E_C` is 4–8× smaller than derived** (≈ 105 at 5 genes, ≈ 27 at 1 gene, vs 200–800), so
  the profile is ≈ 5–8 % of the GA stage, not comparable to it.
- **HSD(−4) is not `sfLDOF`.** README 7.4/7.5 used `sfLDOF` on H2; at γ = −4 the M1 optimum
  is 0.6802 at w1 ≈ 0.44 against README's 0.6707 at w1 0.32 with `sfLDOF` (and nsim 2e5, a
  different draw), so the 7.5 comparison is qualitative.
- **Batteries must not be edited while they run.** `Rscript` reads its file lazily; an edit
  made to `d02_warm_start_profile.R` while the M0/M1 job ran made that job fail with a parse
  error *after* it had saved both designs (no data lost). Scripts were not edited afterwards
  while running.
- The plan's popSize rule is inert for every design here (m + d ≤ 4 → 200), as the plan noted.

## Scripts and run times

All under `dev/review/hsd_opt/`, `Rscript`-runnable from the repo root.

- `d_helpers.R` — shared fixtures and prototype helpers (designs, counting gain, warm-started
  refinement, plateau); sourced by d02/d04/d05.
- `d01_transform_split.R` — 157.5 s (other agent; output `d01_transform_split_output.txt`).
- `d02_warm_start_profile.R` — parts `M0` 1229 s, `M1` 1216 s, `S2` 1791 s (resumed run with
  `D02_SEEDS_S2=1,2`, after a first run that did 3 outer runs in 772 s); `report` merges into
  `d02_warm_start_profile_results.rds` / `_output.txt`. Logs `d02_log_*.txt`.
- `d03_memory_and_copy.R` — 7.7 s (other agent).
- `d04_eval_counts_and_budget.R` — parts `counts` 626 s, `bank_M0` 1023 s, `bank_M1`
  1090 s; `report` writes `d04_eval_counts_and_budget_results.rds` / `_output.txt`. Logs
  `d04_log_*.txt`.
- `d05_grid_step_loss.R` — 218 s.

Wall time for D-b (13:58 to 15:14) ≈ 76 min on two lanes (≈ 2 h 10 min of single-process compute).
Cuts made for wall time: §2 outer loop at seeds {1, 2} (seed 3 only at γ = −4) instead of
{1, 2, 3} (5 of 18 runs dropped); one-call identity check of D6 only at seed 1 (and
`global_search = FALSE`), not seed 2; the popSize-300 configuration at 3 seeds, not 5.
