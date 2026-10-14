# Brief P3: the driver — Stage G, Stage L, profile, dispatch

Branch `gsd-build-spending`, repo `/home/user/multigrain-spending`. You implement phase P3 of
`dev/hsd_opt_build_plan.md` — the critical phase. P1 and P2 are committed; Phase D is committed
and its decisions are in `dev/review/hsd_opt_discovery.md` "Decisions". **Report, do not
commit.** Complete files, not diffs.

## Defaults fixed by Phase D (read the Decisions section; these are binding)

<!-- ORCHESTRATOR FILLS FROM dev/review/hsd_opt_discovery.md -->
- Default `step` for `sf_hsd_opt()`: **{STEP}** (already in P1 as 0.25 — change the default
  there only if Phase D chose otherwise, and say so).
- Profile density: **{DENSITY}** (full sweep of every lattice node, or every second node by
  default with `control` opt-in for the full sweep).
- popSize rule: `min(max({POP_CONST} * (m + d), 200), 500)`.
- Warm-start gate verdict: **{GATE}**; if it failed, the fallback (second start per node from
  the GA archive at a coarser lattice) is part of this phase.
- Spec table-building progress step: **{PROGRESS}** (P1's builder is silent; if Phase D says a
  progress step is needed, add `cli::cli_progress_step()` per node in `.hsd_node_tables()`
  gated on `multigrain_verbosity() != "silent"` — the one allowed edit to P1's file).
- Profile COBYLA budget: same `local_opt` as the main pass (no profile-specific `maxeval`).

## Read first (~60 minutes; do not skip)

1. `dev/hsd_opt_build_plan.md` §1 (rows 2, 4, 5, 6, 9, 10, 11 bind you), §4 "P3" in full.
2. `dev/hsd_opt_design_revised.md` §1.3, §1.4, §1.5 (the result fields), §2.3, §2.5, §3
   items 1, 6, 7, 10, 11, §5.2 in full, §5.3 "P3".
3. `dev/review/hsd_opt_discovery.md` — Decisions, Section 4 (what the prototype profile did
   and its cost), Section 5 (GA with the gene).
4. `R/transform_pvalues_gsd_opt.R` (P1) and `R/objective_function_gsd_opt.R` (P2) in full,
   with `tests/testthat/helper-hsd_opt.R`.
5. `R/optimisation_gsd.R` in full — `graph_optimise_gsd()` (you add the dispatch as its first
   statement), `.graph_optimise_ga_gsd()` and `.graph_optimise_local_gsd()` (your `_opt`
   drivers mirror their RNG order, GA wiring, `nloptr` call, decode, re-score and return
   shapes exactly), `.sample_pvals_gsd()`, `control_prepare_dims()`.
6. `R/optimisation_start.R` (`.build_start_matrix()` — row 1 is the Holm seed, which is why
   warm starts never use it), `R/starting_values.R` (`create_start_params()`,
   `param_to_solution()`), `R/graph_constraint.R` `closest_graph_to_constraints()` (:432),
   `R/post_optim_processing_gsd.R` (`repair_graph()`, `prune_graph_gsd()`), `R/choose_graph.R`,
   `R/calc_power_gsd.R` (`calc_power_pvals_gsd()`), `R/graph_optimal.R` (`graph_optimal()`,
   the `stopping_rule` attachment at :89-91 — you add `spending = NULL` the same way),
   `R/control_prepare.R`, `R/check_gsd.R`.
7. `tests/testthat/test-optimisation_gsd.R` (fixtures, tiny GA control at :23-26, the oracle
   style) and `tests/testthat/test-graph_optimal.R` (the snapshot you must not disturb).
8. `dev/callflows/graph_optimise_gsd.md` "The call chain at a glance".

## What to build — new file `R/optimisation_gsd_opt.R`

### `.graph_optimise_gsd_opt(spec, graph_constraint, trial_success, alpha, stopping_rule, start_graph, global_search, num_threads, control, verbose)`
Mirror `graph_optimise_gsd()`'s validation **on the spec** in the same order: the checks on
`graph_constraint`, `trial_success`, `stopping_rule`, `num_threads`, `control`, `global_search`,
verbose normalisation, `.gsd_check_alpha(alpha, spec)`, `.gsd_check_graph_dims()`,
`.gsd_check_gain_dims()`, `.validate_start_graphs()`, `stop_rule <- .gsd_stop_table(…)`. Then
`user_pop <- control$global_opt$popSize`; `control <- control_prepare_dims(control, nsim =
spec$nsim, m = spec$m, verbose)`; `if (is.null(user_pop)) control$global_opt$popSize <-
min(max({POP_CONST}L * (spec$m + d), 200L), 500L)` (plan row 5). `d <- length(spec$param$names)`,
`p_graph <- ncol(.build_start_matrix(graph_constraint, start_graph))`.
Then Stage G (if `global_search`), Stage L, `choose_graph()`, final evaluation, result.

### `.graph_optimise_ga_gsd_opt(spec, graph_constraint, trial_success, alpha, nsim, global_opts, num_threads, start_graph, verbose, stop_rule, p_graph)`
RNG order: `x0 <- .build_start_matrix(…)` (no draws), then `idx_ga <- sample.int(spec$nsim,
nsim)` (the same call `.sample_pvals_gsd()` makes), then the bank (`.hsd_bank(spec, idx_ga)`),
the closure (`create_obj_func_gsd_opt(…, start_param = spec$param$start)`), suggestions = `x0`
with each start node's gene `(q_start − 0.5) / M_r` appended as extra columns, `lower = rep(0,
p_graph + d)`, `upper = rep(1, …)`, and `GA::ga()` with **exactly** the fixed driver's
`immutable_global_args` and `modifyList(…, global_opts)`. Winner: `sol <- ga_res@solution[1, ]`;
`u <- sol[p_graph + seq_len(d)]` → nodes via `.hsd_node_index()` → `theta_ga` (the node
values); graph genes `sol[seq_len(p_graph)]` decoded with `param_to_solution(process = TRUE)`
and **no** `repair_graph()` (as the fixed GA driver). Validity check as the fixed driver; if
valid, re-score: `.calc_power_pvals_gsd_core()` on the subsample at `theta_ga` (resolve, subset
`idx_ga`), then `graph_shortcut_gsd()` on the full `.gsd_kernel_matrix(resolve(spec,
theta_ga))` and `trial_success$func(res$time)`. Return the fixed driver's list plus
`ga_spending_param = theta_ga`, `archive = attr(fn, "archive")`, `idx_ga`.

### `.graph_optimise_local_gsd_opt(spec, graph_constraint, trial_success, alpha, nsim, local_opts, num_threads, x0, theta0, archive, verbose, stop_rule, p_graph, density)`
- `idx_local <- sample.int(spec$nsim, nsim)` **once** (same call as `.sample_pvals_gsd()`),
  pinned for everything below (plan row 6).
- Helper `node_array(theta)`: `resolve_pvalues_gsd(spec, theta)`, subset the array to
  `idx_local` (`drop = FALSE`), `$nsim <- length(idx_local)`, return the `multigrain_pvals_gsd`.
- Helper `objective(pv)`: `create_obj_func_gsd(spec$m, trial_success$func, hc, tc,
  .gsd_kernel_matrix(pv), spec$K, alpha, num_threads, stop_rule)` — the **unchanged** fixed
  closure, so the P3 oracle (a) is literal.
- **Step 1.** `.hsd_refine_node(x0, pv = node_array(theta0), …)`: `nloptr::nloptr(x0, eval_f
  = function(x) -obj(x), lb = 0, ub = 1, opts = local_opts)` exactly as
  `.graph_optimise_local_gsd()` does, decode with `param_to_solution(process = TRUE)` then
  `repair_graph()`, evaluate the decoded graph on the same array (`obj(encoded)`? no — evaluate
  via the kernel on `(w, G)` as the fixed driver's re-score does, on the `idx_local` rows), count
  evaluations with a wrapper counter. Returns `list(x, w, G, gain, n_eval, nloptr_output)`.
  `x0` when `global_search = FALSE`: `.hsd_encode_graph(gc, start_graph[[1]])` if the user gave
  one (honour it — plan row 9), else `create_start_params(gc)`; `theta0 = spec$param$start`.
- **Step 2, per axis r, outward from the current node in both directions** (`.hsd_profile_axis()`):
  at each node `q` (every node, or every second node plus the bounds if `density == 2`): build
  `pv_q`, `obj_q`; candidates = the previous node's refined `(w, G)`, the current overall best
  `(w*, G*)`, and the archive's best graph at `q` on this axis if it has one (`genes[q, ]`
  decoded; NA → skip); evaluate each on `pv_q` (kernel on `idx_local` rows); start
  `.hsd_refine_node()` from the best candidate's encoding (`.hsd_encode_graph()`); record
  `gain = max(candidate gains, refined gain)` with the corresponding graph, `gain_fixed_graph` =
  the `(w*, G*)` candidate's value, `n_eval`; the paired difference and SE against the current
  optimum via `.hsd_block_se()`. After the axis: move the current node to the argmax if it beats
  the current gain by > 0 (ties keep).
- **Rounds** for `d ≥ 2`: repeat the axis sweep while any axis moved, at most 3 rounds;
  `profile_converged = (no axis moved in the last round)`. At the end, any axis whose stored
  profile was computed through a stale `theta_{-r}` is recomputed once so every reported
  profile passes through `theta*`. For `d = 1` this is one sweep.
- `.hsd_block_se(time_a, time_b, func, B = 100)`: split rows into `min(B, nrow)` contiguous
  blocks; `func()` on each block's time matrix for both; `diff = mean(block diffs)`,
  `se = sd(block diffs) / sqrt(B)`. (Both time matrices come from kernel calls on the same
  `idx_local` rows, so they are paired.)
- `.hsd_encode_graph(gc, w, G)`: `proj <- closest_graph_to_constraints(gc, w, G)`;
  `create_start_params(gc, w0 = proj$hyp_weight, G0 = proj$trans_matrix, sum_to_one_constraint
  = FALSE)`; return the numeric vector. **Never** `.build_start_matrix(…)[1, ]`.
- `.hsd_plateau(value, diff, tol = 5e-4, at)`: the connected run of nodes around `at` with
  `diff <= tol`; returns `c(lo, hi)`.
- Returns the fixed local driver's list plus `local_spending_param = theta*`,
  `spending_profile` (data frame `param, value, gain, gain_fixed_graph, diff, se, n_eval`, with
  attribute `graphs` = list of `(w, G)` per row), `profile_converged`, `idx_local`.

### Final assembly (in `.graph_optimise_gsd_opt()`)
`best <- choose_graph(ga_result, local_result)` (unchanged); `theta_star <- if (best$source
== "global") ga_result$ga_spending_param else local_result$local_spending_param` — check the
exact `source` labels `choose_graph()` uses. `pv_star <- resolve_pvalues_gsd(spec, theta_star)`;
`prune_graph_gsd(pv_star, …)`; `calc_power_pvals_gsd(pv_star, …, custom_power =
list(trial_success = trial_success), stopping_rule = stopping_rule)`; names from the
constraint; `graph_optimal(…, stopping_rule = stopping_rule, spending = list(spending_param =
named theta_star, hyp_spending_param = length-m with NA for inactive, spending = pv_star$spending
(labels), spending_bounds = spec$param$bounds, spending_profile, spending_search (data frame
`param, value, visits, best_fitness_subsample` from the archive; all-NA visits if
`global_search = FALSE`), profile_converged))`.

### Touchpoints
- `R/optimisation_gsd.R`: first statement of `graph_optimise_gsd()`:
  ```r
  if (inherits(pvals, "multigrain_pvals_gsd_spec")) {
      rlang::check_dots_empty()
      return(.graph_optimise_gsd_opt(
          spec = pvals, graph_constraint = graph_constraint, trial_success = trial_success,
          alpha = alpha, stopping_rule = stopping_rule, start_graph = start_graph,
          global_search = global_search, num_threads = num_threads, control = control,
          verbose = verbose
      ))
  }
  ```
- `R/graph_optimal.R`: argument `spending = NULL` after `stopping_rule`; after the
  `stopping_rule` attachment: `if (!is.null(spending)) out[names(spending)] <- spending`.
  Validate it is a named list. **No display changes** (P4).
- Nothing else existing changes. `document()` twice.

## Tests — `tests/testthat/test-optimisation_gsd_opt.R`
The plan §4 P3 checklist (a)–(f) plus the popSize items, every item a `test_that()`, with the
P1/P2 fixtures and the tiny GA control (`popSize 20–30, maxiter 5–20, run 3–5, maxeval ≤ 200,
nsim_global ≤ 1000`). Snapshot nothing in `graph_optimal.md`. (a) is the gate's oracle: the
by-hand `nloptr::nloptr()` call must use the identical closure object, `x0`, `opts` and
`set.seed()` and be `identical()` in `$solution` and `$objective`. (b) is `identical()`, not
`all.equal()`. (c) runs the whole pipeline twice under `withr::with_seed()`.

## Dev benchmark — `dev/review/hsd_opt/gate_p3_benchmarks.R`
As the plan §4 P3 "Dev benchmark": run by the orchestrator, not a test; write it so it reads
`d02_warm_start_profile_results.rds` for the outer-loop reference gains, runs the real
`graph_optimise_gsd()` on the spec for M1 and the §2 design at `nsim = 1e5` (d = 1, 21 nodes),
and the d = 2 coarse check at `nsim = 2e4` (two tagged groups, 5 × 5 nodes; the reference is a
by-hand product-grid loop of fixed transforms + default `graph_optimise_gsd()` at each of the
25 pairs — expensive, so make the 25-pair reference optional behind a flag and default to the
9 pairs nearest the start node), reports the three inequalities and the wall time against the
design §4.2 projection. Do not run it yourself at `nsim = 1e5`; run it once at `nsim = 5e3`
with a `SMOKE = TRUE` flag to prove it executes.

## Gate (run before reporting; the orchestrator re-runs from a clean build)
```r
Sys.setenv(NOT_CRAN = "true"); devtools::load_all()
testthat::test_file("tests/testthat/test-optimisation_gsd_opt.R", package = "multigrain")
testthat::test_file("tests/testthat/test-optimisation_gsd.R", package = "multigrain")       # 79
testthat::test_file("tests/testthat/test-graph_optimal.R", package = "multigrain")          # 33, snapshot untouched
testthat::test_file("tests/testthat/test-choose_graph.R", package = "multigrain")
testthat::test_file("tests/testthat/test-objective_function_gsd_opt.R", package = "multigrain")
testthat::test_file("tests/testthat/test-transform_pvalues_gsd_opt.R", package = "multigrain") # 168
lintr::lint("R/optimisation_gsd_opt.R"); lintr::lint("R/optimisation_gsd.R"); lintr::lint("R/graph_optimal.R")
lintr::lint("tests/testthat/test-optimisation_gsd_opt.R")
```
`git diff --quiet -- tests/testthat/_snaps/graph_optimal.md` must hold. `load_all()` fine;
no `install()`/`clean_dll()`.

## Report
Files; gate counts; the smoke run of the benchmark script; divergences and why; what you could
not verify; what the design did not anticipate. Do not commit.
