<!--
Provenance: build plan written by a Claude Code session on 2026-10-14 from
dev/hsd_opt_design_revised.md (the design) and dev/implement_hsd_opt.md (the job prompt),
after reading the code on gsd-build at 3e55fb5. User decisions of the same day: Phase D runs
straight through and resolves its open defaults by the design's own rules; implementer agents
run on Opus, one per phase, reporting rather than committing; the orchestrating session gates,
reviews, commits and pushes. This file is the live record of the build: it is updated when a
Phase D decision rule resolves or a gate forces a change.
-->

# Build plan: joint optimisation of HSD spending parameters — branch `gsd-build-spending`

## Context

`dev/hsd_opt_design_revised.md` (on `gsd-build` at `3e55fb5`) is the specification: optionally
optimise, jointly with the graph, Hwang–Shih–DeCani spending parameters for hypotheses the user
tags with `sf_hsd_opt()`, while the fixed-spending path stays byte-identical. Job prompt:
`dev/implement_hsd_opt.md`; house rules: `dev/implement_gsd.md` ("Read first", "Constraints") and
the lessons in `dev/gsd_progress.md` ("Stopping rules" entry).

User decisions (2026-10-14): **Phase D runs straight through** and resolves its open defaults by
the design's own rules (§4.2, §4.3 item 5, §3.10); pause only if a gate fails or a measurement
contradicts the design. **Implementers are Opus**, one agent per phase, each reports and never
commits; **Fable orchestrates**: briefs, re-runs every gate from a clean build, reviews the diff,
commits explicit paths, pushes. No pull request.

Deliverables beyond code: `dev/review/hsd_opt_discovery.md` + scripts, record §4.12, NEWS,
pkgdown, article subsection, **`dev/callflows/graph_optimise_gsd_opt.md`** (call flow),
**`dev/hsd_opt_explained.md`** (layperson explanation), `dev/gsd_progress.md` entry.

Environment: R 4.3.3, gsDesign 3.11.0, GA 3.2.5, nloptr 2.2.1, testthat 3.3.2, roxygen2 8.1.1,
lintr 3.4.0; `quarto` R package and `quarto`/`air` binaries missing (article chunks
`#| eval: false` with pasted output, as before); 4 cores, 15 GB. Clean build:
`Sys.setenv(NOT_CRAN="true"); devtools::clean_dll(); devtools::install(quick=TRUE); devtools::load_all()`.

Facts verified during planning that change the design's text: `gsDesign::sfHSD(0.025, c(0.5,1), 0)$spend`
is accepted and linear (`0.0125 0.0250`); the corners −10 and 2 return ordinary doubles. The
"fewer than two distinct positive boundaries" abort is in `.gsd_invert()`
(`R/transform_pvalues_gsd.R:365-376`), **not** in `.gsd_boundary_table()`, so an eager table build
does not by itself prove invertibility. The popSize rule `min(max(40(m+d),200),500)` gives 200
for every benchmark design (m+d ≤ 5). The README §7.4 script uses `nsim = 2e5`; Phase D uses
1e5 (job prompt) and re-derives the M0/M1 reference values on its own draws.

---

## 1. Reconciliations — design vs. code (binding for the agents)

| # | Code fact | Decision |
|---|---|---|
| 1 | `.gsd_spend()` calls `spending(alpha, t)` with two args; `sf_hsd_opt()` returns `function(alpha, t, param)`. | Spec builder wraps per node: `.hsd_node_fn(f, node) = function(a, t) f(a, t, param = node)` before every `.gsd_check_spending()` / `.gsd_boundary_table()` call. `.gsd_spending_list()` refuses a `multigrain_sf_hsd_opt` function on the fixed path unless `allow_opt = TRUE` (set only by the spec builder). Table builds are wrapped in `rlang::try_fetch()` so a well-ordering abort is re-raised naming hypothesis, node value and the "narrow `bounds`" hint. **Plus** `.hsd_check_invertible(tab)`: `sum(!duplicated(col) & col > 0) >= 2` per column at build time, so `.gsd_invert()`'s abort cannot fire inside the search. |
| 2 | `create_obj_func_gsd()` decodes with `split_theta()`/`recover_*()`, `recover_full_trans_matrix()` `stop()`s on a length mismatch; GA winner decoded with `param_to_solution(process=TRUE)` **without** `repair_graph()`. | Parameter genes **appended at the end**: `x = c(graph genes [p_graph], u_1..u_d)`, `p_graph = ncol(.build_start_matrix(gc, start_graph))`. Closure slices `xg <- x[seq_len(p_graph)]`, `u <- pmin(pmax(x[p_graph + seq_len(d)], 0), 1)`; decode/penalty block **copied verbatim** from `create_obj_func_gsd()` (house rule: duplicate, don't share; the P2 oracle makes it safe). GA winner sliced the same way, then decoded exactly as `.graph_optimise_ga_gsd()` does. `lower`/`upper` = `rep(0/1, p_graph + d)`; `suggestions` = `.build_start_matrix()` rows with the start-node gene `(q_start − 0.5)/M_r` appended. |
| 3 | `.gsd_kernel_matrix()` reads `$pvals`; dimension helpers read `$nsim,$m,$K,$alpha`. | Spec stores raw draws under **`$raw`**, never `$pvals`; class `multigrain_pvals_gsd_spec` only (no inheritance). `check_pvals_gsd()` refuses it with a hint; a kernel-reaching call fails on `NULL` dims. |
| 4 | `check_pvals_gsd(pvals)` is the first statement of `graph_optimise_gsd()` (:150). | New first statement: `if (inherits(pvals, "multigrain_pvals_gsd_spec")) { rlang::check_dots_empty(); return(.graph_optimise_gsd_opt(spec = pvals, <all other args>)) }` — dots check inside the branch so the fixed path's check order is untouched. Handoff in `transform_pvalues_gsd()` between `.gsd_spending_list()` (:126) and the `.gsd_check_spending()` loop (:128), passing `spending_label`; the spec builder runs the spending check itself on fixed entries. |
| 5 | popSize rule `min(max(40L*m,200L),500L)` in `control_prepare_dims()` (R/optimisation_gsd.R:356); user `popSize` wins. | `control_prepare_dims()` is **not edited**. In `.graph_optimise_gsd_opt()`: `user_pop <- control$global_opt$popSize`; after `control_prepare_dims(control, nsim, m = spec$m, verbose)`, `if (is.null(user_pop)) control$global_opt$popSize <- min(max(40L*(spec$m + d), 200L), 500L)`. |
| 6 | `.sample_pvals_gsd()` always draws a permutation, even at `nsim_local == N`; `.graph_optimise_local_gsd()` draws its own subsample internally. | Stage L draws `idx_local <- sample.int(N, nsim_local)` **once in the opt driver** (same RNG order as today) and subsets every resolved node array by it (`drop = FALSE`). Step 1 of Stage L is therefore `.hsd_refine_node()` (an `nloptr` call identical to the local driver's, on the pinned rows), not the local driver itself; the P3 oracle (a) compares against a by-hand `nloptr` call with the same closure, `x0`, options and seed. |
| 7 | `.gsd_spending_labels()` is deparse-based; `print.multigrain_pvals_gsd` never shows spending (only `summary()`). | New `.hsd_labels()`: tagged hypotheses → `sfHSD(param = <value>)` on a resolved object, `sfHSD(param = <name> in [lo, hi])` on the spec; fixed ones → list name if all named, else the i-th deparsed argument of a captured `list(...)` call, else the existing fallback. Fixed `print()` untouched; the spec gets its own `print()`. |
| 8 | `print.multigrain_graph_optimal` order: class, weights, matrix, trial-success fn, [stopping rule], [value]; `_snaps/graph_optimal.md` GSD block must stay byte-identical. | New blocks **after the "Trial success function" block and before the stopping-rule block** (design §1.5 "after the trial success function"; `R/graph_optimal.R:235-241`), via `.hsd_cat_print(x)` guarded by `!is.null(x$spending_param)`; in `summary()` via `.hsd_cat_summary()` after `summary(object$trial_success)` and before the stopping-rule summary (:272-276). Both helpers live in `R/optimisation_gsd_opt.R`; `R/graph_optimal.R` gets only the argument and two guarded calls. |
| 9 | `.build_start_matrix(...)[1, ]` is always the Holm seed → `global_search=FALSE` ignores a user `start_graph` (R/optimisation_gsd.R:229-233; R/optimisation.R:188-192). | Pre-existing: **report** (discovery report "Pre-existing issues", progress entry, record §10 item 16); **do not fix, do not replicate**. The opt path honours `start_graph[[1]]` by encoding it with `.hsd_encode_graph()` = `closest_graph_to_constraints()` + `create_start_params(gc, w0, G0, sum_to_one_constraint = FALSE)` (the user-row branch of `.build_start_matrix()`). |
| 10 | `graph_optimal()` attaches `stopping_rule` only when non-NULL (:89-91). | `graph_optimal()` gains `spending = NULL` after `stopping_rule`; when non-NULL, `out[names(spending)] <- spending` attaches `spending_param`, `hyp_spending_param`, `spending`, `spending_bounds`, `spending_profile`, `spending_search`, `profile_converged`. Fixed path never passes it. |
| 11 | `choose_graph()` compares `ga_trial_success` vs `local_trial_success` on full-sample gains. | Reused unchanged; both winners are scored on the full sample at their own ϑ, so the comparison is on gains. If the GA wins, ϑ* = ϑ_GA and the reported profile is still Stage L's (it started at ϑ_GA). Stated in roxygen. |

Fixed by design: lattice `step ∈ {0.5, 0.25, 0.125}`; default `bounds = c(-4, 1)`, user bounds within `[-10, 2]`, lower < upper; `start = NULL` → node nearest the midpoint; one recycled `sf_hsd_opt()` = one shared parameter; clamping, not penalties; `B = min(100, nrow)` batch-means blocks; plateau tolerance `5e-4`; ≤ 3 profile rounds for `d ≥ 2`; `max_bank_gb = 4` (warn > 2 GB).

---

## 2. Branch and orchestration

**Step 0 — branch and plan on the branch (first action after approval).**
`git checkout -b gsd-build-spending 3e55fb5`; save **this plan** verbatim as
`dev/hsd_opt_build_plan.md` (with a two-line provenance header: written by this session on
2026-10-14 from `dev/hsd_opt_design_revised.md`; user decisions of the same day); commit
"Add the build plan for the HSD spending-parameter optimisation"; `git push -u origin
gsd-build-spending`. Every later phase works and commits on this branch. (Design text says
`gsd-build-hsc-opt`; the user's name wins. Report paths stay `dev/review/hsd_opt_discovery.md`,
`dev/review/hsd_opt/`.) Neither branch exists yet. The plan file on the branch is updated
whenever a Phase D decision rule resolves or a gate forces a change, so it stays the live record.

**Prerequisite gate (orchestrator).** Clean build at `3e55fb5`; `test_file()` on
`test-transform_pvalues_gsd.R`, `test-optimisation_gsd.R`, `test-graph_optimal.R`,
`test-check_gsd.R`, `test-objective_function_gsd.R`; record counts (the "unchanged" baseline);
`Rscript dev/gsd_identity_check.R` green.

**Per-phase cycle.**
1. Orchestrator writes the brief to `dev/review/hsd_opt/briefs/<phase>.md` (scope, files,
   functions, gate checklist, applicable reconciliations, "report, don't commit", "complete files,
   not diffs") and launches **one Opus agent** (`model: "opus"`, background).
2. Agent implements, runs its own gate, returns a report (files, gate output, divergences, what it
   could not verify). Report saved beside the brief.
3. Orchestrator: read the diff; clean build; every gate `test_file()`; `lintr::lint()` on touched
   files (only `%||%` backport warnings accepted); `LC_ALL=C.UTF-8 Rscript -e 'devtools::document()'`
   twice (second a no-op); `git status` for CRLF snapshot churn (`git checkout -- <snap>`);
   `Rscript dev/gsd_identity_check.R` after P1–P4 (never concurrently with a GA battery: it installs
   twice).
4. `git add <explicit paths>`; commit; push.

**Pause for the user only if:** a gate still fails after one fix round with the same agent; the
identity script fails; a d01 node in `[-10, 2]` fails well-ordering or invertibility; the d02 gate
fails and its stated fallback also fails; d04 shows the default run > 2× the §4.2 projection
(3.5–6.5 min); an implementer proposes a design change.

**Parallelism.** d01 first (minutes). P1 starts once d01 is in (it needs: bounds valid across the
range, `param = 0` accepted). d02/d04/d05 run as two background jobs at `num_threads = 2` while
P1 and P2 are implemented and gated (their tests are `nsim ≤ 2000`, seconds). P3's brief needs
d02 (warm-start gate), d04 (profile density, popSize) and d05 (Δ), so P3 waits for the discovery
report. P4, P5 serial.

**Commits:** one per phase (D, P1, P2, P3, P4, P5) + a final progress-entry commit. Because the
worker has restarted twice today, every phase commits as soon as its gate passes.

---

## 3. Phase D — discovery (Opus; `nsim = 1e5`, `num_threads = 2`, under `dev/review/hsd_opt/`)

Fixtures: README §7.4/7.5 shared design (m=2, powers 0.90/0.80, corr 0.5, `info_frac c(0.5,1)`,
alpha 0.025, values (0.6,0.4), `d = c(1,0.8)`, full-recycling G, seed 20261007; M0 no stop, M1
`stopping_rule_gsd(r1)`) at nsim 1e5, and the README §2 design (m=3, K=2, sfLDOF). Each script
`Rscript`-runnable from the repo root, writing `<name>_results.rds` and `<name>_output.txt`.

| Script | Measures | Tables | Est. |
|---|---|---|---|
| `d01_transform_split.R` | (a) `.gsd_boundary_table()` time at G∈{128,256,512,1024}×K∈{2,3} and `.gsd_invert()` time at N∈{1e4,1e5}×K×m∈{2,3}; max relative interpolation error vs G=1024 at 200 random p-values, ϑ∈{−4,−2,1}. (b) Well-ordering **and invertibility** at every node of `seq(-10, 2, 0.25)` on schedules (0.5,1), (1/3,2/3,1), (0.25,0.5,0.75,1), (0.6,1); corners ϑ=−10, t=(0.1,1) and ϑ=2, t=(0.9,1); `sfHSD(1,t,ϑ)$spend` identity and linearity in `a` at 1e-14 and 0.025; `param = 0`. | D1 ms per table / per inversion; D2 error by G and ϑ; D3 pass/fail per node×schedule, min distinct boundaries | 3–4 min |
| `d02_warm_start_profile.R` | M0, M1, §2 design: outward warm-started profile over 21 nodes (prototype: resolve per node with `function(a,t) sfHSD(a,t,g)$spend`, `create_obj_func_gsd()` + `nloptr` from the better of neighbour's optimum and current best) vs a full default `graph_optimise_gsd()` at nodes {−4,−3,−2,−1,0,1} × 3 seeds. | D4 per design: node, profile gain, full gain per seed, gap, plateau, gate verdict | ~50 min (54 runs, two at a time) |
| `d03_memory_and_copy.R` | `tracemem()` on an N_GA×mK matrix across a kernel call and a column assignment; copy vs kernel time at `nsim_global = 5e4`. | D5 copy ms, kernel ms, ratio | 1 min |
| `d04_eval_counts_and_budget.R` | §2 design: counter-wrapped `trial_success$func`, drivers in `graph_optimise_gsd()` order (GA, local, `choose_graph`, `prune_graph_gsd`, `calc_power_pvals_gsd`) for seeds 1, 2 and `global_search = FALSE`; GA `@iter`, counts and times per stage. Warm-started COBYLA counts `E_C` at ϑ∈{−2.25,−1.75,−3,−1,−4,1} from the −2 optimum. Prototype bank closure (gene appended, clamp, node map, column copy); GA with/without the gene on M0/M1, 5 seeds, popSize 200 (note the rule bites only at m+d ≥ 6). | D6 counts/times by stage; D7 `E_C` by node; D8 GA ± gene: generations, evaluations, wall, best node vs d02 argmax | ~30 min |
| `d05_grid_step_loss.R` | At the d02 optimum on M1 and §2: warm-started COBYLA on Δ=0.025 within ±0.5; loss of Δ∈{0.5,0.25,0.125} = fine max − best coarse node. | D9 loss by Δ and design | ~8 min |

**Decision rules (written into the report's "Decisions" section with the number that decided each):**
- **Default `step`:** coarsest Δ ∈ {0.5, 0.25, 0.125} with loss ≤ 4e-5 (README 4.2 optimiser spread) on both d05 designs; else 0.25.
- **Profile density:** full 21-node d=1 profile is the default if its cost (21 × (inversion + 3 kernel calls + `E_C` × per-eval)) ≤ the GA stage time on the same design; else every second node by default (`control` opt-in for the full sweep). Record which.
- **popSize:** `min(max(40(m+d),200),500)` stands unless D8 shows the GA with the gene needs > 1.5× the generations without it at equal popSize → constant becomes 60.
- **Warm-start gate:** max over d02 nodes inside the reported plateau of (full gain − profile gain) ≤ 5e-4 on all three designs → Stage L as designed; else the §4.3-2 fallback (second start per node from a coarser GA archive); if that fails too → pause.
- **Table building:** if eager-table time at m=3, K=2, M=21, G=1024 > 25 % of the d04 default run (also check K=3 and a two-schedule spec), the spec builder reports a `cli::cli_progress_step()` per node gated on verbosity; default `grid_size` unchanged; roxygen quotes the D2 error for `grid_size = 512`.
- **COBYLA budget:** same `local_opt` for profile nodes (oracle holds); no profile-specific `maxeval`.

**Resolved 2026-10-14 by Phase D (`dev/review/hsd_opt_discovery.md` "Decisions"):** default
`step` 0.25 (d05 loss 2.12e-5 at 0.25, 1.31e-4 at 0.5); full 21-node profile by default (6–8 % of
the GA stage); popSize constant 40 stands (gene cost 0.90–0.98× the generations); warm-start gate
PASS on M0, M1 and the §2 design (in-plateau gaps 2.1e-5, 2.3e-5, 3.7e-4) so Stage L is as
designed; progress step per node table, verbosity-gated (two schedules at K = 3: 30 % of a run);
same `local_opt` at profile nodes (E_C 23–128, 2.6 % of `maxeval`). No pause condition fired.
Noted for P3/P5: optimiser noise on the §2 design is ≈ 4e-4, close to the 5e-4 plateau tolerance;
the default bounds put the §2 optimum on the lower bound (gain still rising below −4).

Report also has "Pre-existing issues" (reconciliation 9 with both file:line refs; the two older
ones from `implement_gsd.md` restated). **Gate D:** every script runs end to end from a clean
build; one table per item of design §4.1–4.3; the Decisions section names the five defaults.

---

## 4. Phases P1–P5

**Shared test fixture** `tests/testthat/helper-hsd_opt.R`: `withr::with_seed(7, simulate_pvalues_gsd(c(0.9,0.85,0.8), corr_matrix = diag(3), info_frac = c(0.5,1), nsim = 2000))`;
d=1 spec: `spending = list(gsDesign::sfLDOF, sf_hsd_opt(id="sec", bounds=c(-4,-2), step=0.5), sf_hsd_opt(id="sec", bounds=c(-4,-2), step=0.5))`, `grid_size = 64L` (5 nodes); `hsd_fixed_fn(g) = function(a,t) gsDesign::sfHSD(a,t,g)$spend`;
tiny control as `test-optimisation_gsd.R:23-26` (`maxiter 20, run 5, popSize 30, maxeval 200, nsim_global 1000`);
d=2 fixture: m=4, K=2, nsim 1000, two groups with 3-node lattices (`bounds=c(-3,-2), step=0.5`), one short-circuited hypothesis (`info_frac` row `c(NA, 1)`).

### P1 — spec class, lattice, eager tables, resolve (Opus; overlaps Phase D after d01)

**New** `R/transform_pvalues_gsd_opt.R`: `sf_hsd_opt(id = NULL, bounds = c(-4, 1), step = 0.25, start = NULL)` [export], `is_sf_hsd_opt()`, `.has_hsd_opt()`, `.hsd_node_fn()`, `.hsd_groups()`, `.hsd_lattice()`, `.hsd_node_tables()` (dedup by (schedule, node)), `.hsd_check_invertible()`, `.gsd_invert_hyp()` (steps 4–7 of `.gsd_transform_hyp()` on a prebuilt table, duplicated), `.transform_pvalues_gsd_spec()`, `new_pvals_gsd_spec()` (fields `raw, nsim, m, K, alpha, info_frac, look_back, grid_size, fixed, fixed_tables, param (names, bounds d×2, step, start, lattice), group (length m, NA fixed/short-circuited), node_tables, spending, short_circuit, spending_label`), `print.multigrain_pvals_gsd_spec()` [export], `resolve_pvalues_gsd(spec, spending_param)` [export; off-lattice abort `abs(v − node) > 1e-12`; returns `new_pvals_gsd()`], `.hsd_labels()`.
**Touch** `R/transform_pvalues_gsd.R` (one guarded `if` between :126 and :128; `allow_opt` on `.gsd_spending_list()`), `R/check_gsd.R` (spec hint before `stop_input_type()`). `d = 0` after short-circuits → plain object + `cli::cli_inform()`. No RNG anywhere.
**Tests** `test-transform_pvalues_gsd_opt.R`:
- [ ] `sf_hsd_opt()` validation messages name the limits (snapshot: lower ≥ −10, upper ≤ 2, lower < upper, `step` set, group consistency)
- [ ] single recycled tag = one shared parameter (`d == 1`, `group == rep(1, m)`)
- [ ] untagged call untouched: class `multigrain_pvals_gsd`, `identical()` to a by-hand `.gsd_transform_hyp()` assembly
- [ ] all-short-circuited tagged spec → `multigrain_pvals_gsd` `identical()` to the fixed path, message emitted
- [ ] node tables shared across hypotheses on one schedule (`identical()` to per-hypothesis `.gsd_boundary_table()`)
- [ ] `resolve_pvalues_gsd(spec, node)` `identical()` on `$pvals,$tables,$info_frac,$look_back,$alpha` to `transform_pvalues_gsd(raw, spending = list(sfLDOF, hsd_fixed_fn(node), hsd_fixed_fn(node)), grid_size = 64L)` at all 5 nodes; labels `"sfHSD(param = -3.5)"` etc.
- [ ] off-lattice value refused; non-well-ordered node aborts naming hypothesis, node, hint (test double of class `multigrain_sf_hsd_opt`); `check_pvals_gsd()` refuses the spec with the hint (snapshot); spec `print()` succeeds
- [ ] unchanged: `test-transform_pvalues_gsd.R`, `test-check_gsd.R`, `test-optimisation_gsd.R` baselines; identity check; lintr; `document()` no-op

### P2 — objective closure and bank (Opus; after P1)

**New** `R/objective_function_gsd_opt.R`: `create_obj_func_gsd_opt(spec, idx, power_criterion, hyp_constraint, trans_constraint, K, alpha, num_threads, stop_rule, p_graph)` (working matrix = `.gsd_kernel_matrix()` of `resolve(spec, start)[idx]`; fixed columns set once; archive environment `best/visits/genes` per axis as `attr(fn, "archive")`, visits counted only for evaluations reaching the kernel), `.hsd_bank(spec, idx)` (per active hyp × node, `N_idx × K` via `.gsd_invert_hyp()`), `.hsd_node_index()` (clamp + `min(M, 1 + floor(u*M))`), `.check_bank_memory()`.
**Tests** `test-objective_function_gsd_opt.R`:
- [ ] node oracle: 20 random graph gene vectors × every node (d=1 fixture, and the full node product of the d=2 fixture): `identical(fn(c(x, u_q)), create_obj_func_gsd(..., pvals = .gsd_kernel_matrix(resolve(spec, node_q)[idx]))(x))`
- [ ] genes outside [0,1] → end nodes (`−0.3` → 1, `1.7` → M, `1` → M); archive matches a by-hand replay of 30 evaluations
- [ ] working matrix not copied per evaluation (`tracemem()`; generous bound, exact numbers in D5); memory guard warns/aborts at the stated sizes (synthetic byte counts)
- [ ] unchanged: `test-objective_function_gsd.R` (21); identity check; lintr

### P3 — driver (Opus; after P2 and Phase D) — the critical phase

**New** `R/optimisation_gsd_opt.R`: `.graph_optimise_gsd_opt()` (spec validation via the unchanged `.gsd_check_alpha`, `.gsd_check_graph_dims`, `.gsd_check_gain_dims`, `.validate_start_graphs`, stop table; control per reconciliation 5), `.graph_optimise_ga_gsd_opt()` (RNG: `idx_ga` → `GA::ga()`; winner sliced, decoded as the fixed GA driver, re-scored full-sample at ϑ_GA), `.graph_optimise_local_gsd_opt()` (`idx_local` once; step 1 `.hsd_refine_node()` at ϑ_GA from the GA genes, or from `.hsd_encode_graph(start_graph[[1]])`/Holm when `global_search = FALSE`; outward sweep per axis with three candidate starts; `gain` = max of the four values; fixed-graph column; `diff`, `se` by `.hsd_block_se()`; rounds for d ≥ 2 (max 3, `profile_converged`); stale axes recomputed once at the end — a recompute may expose `diff < 0`: keep the sign honest, set `profile_converged = FALSE`), `.hsd_profile_axis()`, `.hsd_refine_node()`, `.hsd_encode_graph()`, `.hsd_block_se()`, `.hsd_plateau()`; then `choose_graph()`, `prune_graph_gsd()` + `calc_power_pvals_gsd()` on `resolve(spec, ϑ*)`, `graph_optimal(..., spending = list(...))`. `spending_profile` columns `param, value, gain, gain_fixed_graph, diff, se, n_eval` (+ attr `graphs`); `spending_search` columns `param, value, visits, best_fitness_subsample`. Defaults for `step`, profile density, popSize constant from Phase D.
**Touch** `R/optimisation_gsd.R` (dispatch, reconciliation 4), `R/graph_optimal.R` (`spending = NULL` attachment only).
**Tests** `test-optimisation_gsd_opt.R` (tiny GA control):
- [ ] (a) warm-start oracle: `.hsd_refine_node(x0, sub)` `identical()` to a by-hand `nloptr::nloptr()` with the same closure, `x0`, `local_opts`, seed; `param_to_solution(.hsd_encode_graph(gc, w, G))` reproduces `(w, G)` for three graphs
- [ ] (b) `identical(res$power, calc_power_pvals_gsd(resolve_pvalues_gsd(spec, res$spending_param), res$hyp_weight, res$trans_matrix, custom_power = list(trial_success = gain), stopping_rule = rule))`
- [ ] (c) same seed twice → `identical()` on every component except `trial_success`
- [ ] (d) m=4, d=2 fields: `length(spending_param) == 2`, names = ids, lattice membership, `NA`s for fixed and short-circuited, `nrow(spending_profile) == M_1 + M_2`, `diff == 0 && se == 0` at ϑ*, `gain >= gain_fixed_graph`, `is.logical(profile_converged)`
- [ ] (e) `global_search = FALSE` with and without a user `start_graph` (start honoured via `.hsd_encode_graph()`); (f) stopping rule end to end: `stop_distribution` present, `stopping_rule` attached
- [ ] popSize rule applied; user `control_global(popSize = 31)` survives
- [ ] unchanged: `test-optimisation_gsd.R` (79, Figure 3b ran), `test-graph_optimal.R` (snapshot byte-identical), `test-choose_graph.R`; identity check; lintr
**Dev benchmark** `dev/review/hsd_opt/gate_p3_benchmarks.R` (orchestrator, nsim 1e5): joint gain ≥ max over the 21 d02 nodes of the outer-loop gain − 5e-4 on M1 and §2; d=2 coarse (nsim 2e4, 5×5) joint gain ≥ product-grid max − 5e-4; joint gain ≥ fixed-spending gain at the start node − 5e-4; wall time at m=3, K=2, nsim 1e5, d=1 vs §4.2 (pause if > 2×).

### P4 — reporting (Opus; after P3)

**Touch** `R/graph_optimal.R` (two guarded calls, reconciliation 8), `R/optimisation_gsd_opt.R` (`.hsd_cat_print()`: "Spending functions:" per-hypothesis label with `[optimised: <name>, bounds a to b, step s]` on a group's first member, `[optimised: <name>]` on later members; "Gain profile in <name> (graph re-optimised at each value): within 5e-4 of the maximum on [l, u]; fixed-graph plateau [l', u']"; `.hsd_cat_summary()`: profile table per parameter (`round(, 4)`) and "GA visited k of M nodes").
**Tests** **new file** `test-graph_optimal_hsd_opt.R` with its **own** snapshot `_snaps/graph_optimal_hsd_opt.md`: hand-built `graph_optimal()` objects (no optimiser run) — d=1 with a 5-row profile; d=2 with two axes, a short-circuited hypothesis, `profile_converged = FALSE`; `.hsd_plateau()` on a profile with a gap (connected run only).
**Gate:** new snapshots; `test-graph_optimal.R` baseline and `git diff --quiet -- tests/testthat/_snaps/graph_optimal.md`; `test-plot_graph_optimal.R` baseline; `test-optimisation_gsd_opt.R` gains "print() and summary() of a real result succeed"; identity check; lintr.

### P5 — documentation and record (Opus; after P4)

- roxygen: `sf_hsd_opt()`, `resolve_pvalues_gsd()`, `transform_pvalues_gsd()` (`@param spending` tagged case, `@returns` spec case), `graph_optimise_gsd()` (`@param pvals` spec case, `@returns` new fields, the `choose_graph` note of reconciliation 11); `LC_ALL=C.UTF-8 devtools::document()`; `_pkgdown.yml` "Group sequential designs" adds `sf_hsd_opt`, `resolve_pvalues_gsd` after `transform_pvalues_gsd`.
- `vignettes/articles/trial_success.qmd`: `### Optimising the spending parameter` after "Stopping for efficacy" (`#| label: hsd-opt`, `#| eval: false`), output pasted from `dev/review/hsd_opt/article_example.R` (README §7.5 design under M1, LDOF primary, `sf_hsd_opt()` on H2, nsim 1e5).
- `NEWS.md` `## New functionality`: `sf_hsd_opt()`, `resolve_pvalues_gsd()`, the optimised path and result fields. No bug-fix entry for reconciliation 9 (reported, not fixed).
- `dev/gsd_design_record.md`: **§4.12 Spending-parameter optimisation `[Rev 2026-10-14]`** between line 259 and `## 5.` (pattern of 4.11: **Decision.** LOCKED …; italic bullets *Tagging*, *Spec object and eager tables*, *Lattice and snapping*, *Search (Stage G, Stage L)*, *Reporting*, *Dispatch and touchpoints*, *Error control*, *Defaults (from Phase D)*; **Rationale.** with `[measured]` numbers; **Alternatives.**: blending, penalties, fixed-graph sweep, outer product grid, transform inside the objective, inheriting the class); `[Rev 2026-10-14]` paragraph in §4.10 listing the touchpoints; §10 item 16 (reconciliation 9 + anything d04 leaves open); 4.5 paragraph gets a pointer to 4.12.
- **`dev/callflows/graph_optimise_gsd_opt.md`** (house format): "# Call flow: `graph_optimise_gsd()` on a spending spec and its helpers"; intro naming the three new files and four touchpoints and the fixed helpers reused unchanged; "## The call chain at a glance" (1 handoff → spec builder; 2 dispatch; 3 validation and control; 4 Stage G; 5 Stage L step 1; 6 profile sweep; 7 `choose_graph()`; 8 resolve, prune, power, `graph_optimal()`; 9 display); `---`; one `## fn() -- desc` per function (§5.1 list + `.hsd_check_invertible`, `.hsd_encode_graph`, `.hsd_node_index`, `.hsd_cat_*`) with `### Purpose / Inputs and output / How it works / Why it is done this way / Edge cases`; `[Rev 2026-10-14]` paragraph in `dev/callflows/graph_optimise_gsd.md` pointing to the dispatch.
- **`dev/hsd_opt_explained.md`** (plain English, ~150 lines, one formula at most): what an alpha-spending function is (a schedule for how much of the 2.5 % may be used at each look); what gamma does (−4 saves nearly everything for the final analysis, 1 spends early; `h(0.5)` ≈ 0.12 vs 0.62); why it matters more once a trial can stop (README 7.5: 0.6707 → 0.6903 when H2 goes Pocock-like, 14× the Holm gap) and why jointly (w1 moved 0.32 → 0.46); how the search works (21 pre-checked candidate values; the GA tries graph and gamma together; the graph is re-tuned at every gamma to draw the profile); how to read the profile table and the two plateau intervals; what the result means and its limits (chosen on the same simulated trials that score it — check on fresh draws with `resolve_pvalues_gsd()`); what it does not change (error control, kernel, fixed-spending results).
- **Gate P5:** clean build; `document()` twice (no-op); `pkgdown::check_pkgdown()` on a scratch copy with `.qmd` → `.Rmd`; every GSD `test_file()` at baseline/new counts; identity check; lintr.

**Final commit:** `dev/gsd_progress.md` entry `## Spending-parameter optimisation (2026-10-14)` in the Stopping-rules format: Status table (Prerequisite, D, P1–P5: Phase | Commit | Gate | Review verdict), What was built (file by file), Decisions applied (Phase D rules and outcomes), Where diverged (reconciliations 1–11 vs the design), What was not verified (quarto render; `max_bank_gb` at real sizes; d ≥ 3), What the brief did not anticipate (invertibility at build time; popSize rule inert at m+d ≤ 5; `nsim = 2e5` in the 7.4 script; the local driver's internal subsample), Environment.

---

## 5. Risks and fallbacks

- **Eager tables 20–30 s per spec at G=1024** (more at K=3 / two schedules): progress step by the §3 rule; roxygen quotes the `grid_size = 512` error from D2; spec reuse across gains/rules is the structural mitigation.
- **d02 gate fails:** fallback second start per node from a coarser GA archive; if that fails, pause with the gap table.
- **GA gene adds no value at d=1:** expected; report; keep (needed for d ≥ 2; its archive seeds the profile).
- **COBYLA at profile nodes too slow:** density rule; `n_eval` column makes it visible.
- **Memory:** bank at defaults ≈ 34 MB; Stage L ≈ 5 MB per node array; guard exercised with synthetic sizes only.
- **Snapshot churn:** new snapshots in new files; never re-accept `graph_optimal.md`; `git diff --stat tests/testthat/_snaps` before every commit.
- **Core contention:** batteries at 2 threads, two at a time; identity script never concurrent with a battery.
- **Worker restarts:** briefs/reports on disk under `dev/review/hsd_opt/briefs/`; commit per phase.

## 6. End-to-end verification (orchestrator, after P5)

1. Clean build; `document()` twice; `git status` clean after the second.
2. `test_file()` on: `test-transform_pvalues_gsd.R`, `test-transform_pvalues_gsd_opt.R`, `test-check_gsd.R`, `test-objective_function_gsd.R`, `test-objective_function_gsd_opt.R`, `test-optimisation_gsd.R`, `test-optimisation_gsd_opt.R`, `test-graph_optimal.R`, `test-graph_optimal_hsd_opt.R`, `test-calc_power_gsd.R`, `test-post_optim_processing_gsd.R`, `test-stopping_rule_gsd.R`, `test-optimisation.R`, `test-calc_power.R`, `test-gsd_guards.R`, `test-plot_graph_optimal.R` — 0 fail, 0 skip (NOT_CRAN set); existing files at baseline counts, new files at their phase counts.
3. `Rscript dev/gsd_identity_check.R`: `.Random.seed` identical; every component identical except `trial_success$func`.
4. GSD fixed-path identity `dev/review/hsd_opt/gsd_fixed_path_identity.R` (two private libraries: `3e55fb5` via `git worktree` vs the working tree, as the identity script does): untagged `graph_optimise_gsd()` on the `test-optimisation_gsd.R` fixture with `stopping_rule_gsd(r1)`, seed 202 → `identical()` on every component except `trial_success`; untagged `transform_pvalues_gsd()` → `identical()` on the whole object.
5. Acceptance `dev/review/hsd_opt/acceptance_m3.R`: §2 design at nsim 1e4, `spending = list(gsDesign::sfLDOF, sf_hsd_opt(id = "sec"), sf_hsd_opt(id = "sec"))`, `stopping_rule_gsd(r1)`, `num_threads = 2`, default control → completes; `print()`/`summary()` show the spending block and profile; `identical(res$power, calc_power_pvals_gsd(resolve_pvalues_gsd(spec, res$spending_param), ...))`; `spending_param` on the lattice; wall time recorded.
6. `lintr::lint()` on every touched R/test file; `pkgdown::check_pkgdown()` on the scratch copy; `git log --stat` shows no `man/` hand edits and no change to `_snaps/graph_optimal.md`; branch pushed.
