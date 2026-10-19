# Families of hypotheses: gate results

*Branch `family-opt` (based on `main` `16f66d0`), 19 October 2026. Plan: `~/.claude/plans/zazzy-kindling-nebula.md`; design: `families_design.md`. All changes are uncommitted.*

## Decisions taken with Alex

1. **Raw p-values everywhere (Option A).** Every function takes raw `pvals` plus `families`. `.prepare_pvals()` runs right before each kernel call, and once, outside the closure, in `create_obj_func()`. `pvals` therefore always means raw p-values.
2. **q inside the kernel.** The family p-value is the largest within-family adjusted p-value, computed per trial in the kernel. This removes the design's separate `q` matrix and its row-index plumbing.
3. Within-family tests `"hommel"` (default), `"hochberg"` and `"holm"`. Serial and parallel kernels. `family_power` in `calc_power_pvals()`. D4 is an issue draft only (`debug/d4-issue.md`).

## F1: kernel and adjustment

`dev/review/families/verify.R` (output in `results/verify.txt`):

| Check | Result |
| --- | --- |
| Kernel vs literal Algorithm 5.7, 60 random cases × 400 trials per test and rejection order (m 2–12, M 1–6, non-contiguous members, row sums 1 or below 1) | 60/60 for each of Hommel, Hochberg and Holm, with lowest-index and random order |
| All-singleton families vs `graph_shortcut()`, 40 random graphs, m 2–10 | 40/40 `identical()` |
| `adjust_pvals_families()` vs `stats::p.adjust()`, k 1–7, 10,000 rows | `identical()` for all 21 method × k cells |
| Serial vs parallel (1–4 threads), 20 cases × 1e5 trials | 20/20 `identical()` |

The design's worked example (§1.4) is a unit test. It reproduces the rejected set H1, H6–H12, H17, H18 and the Hommel values of the table.

**Option A cost** (`timing.R`, m = 20, M = 6, the drawn allergy graph, serial, laptop on AC):

| N | adjust (s) | kernel (s) | ratio |
| --- | --- | --- | --- |
| 5e4 | 0.011 | 0.005 | 2.10 |
| 2e5 | 0.047 | 0.022 | 2.14 |
| 1e6 | 0.235 | 0.105 | 2.24 |

An earlier ad-hoc run on the same code gave 0.58 s and 0.26 s at N = 1e6 (ratio 2.2–2.4), so absolute times vary with machine state while the ratio is stable. One adjustment costs about two kernel calls, which is within the gate.

The re-adjustments only happen in scoring and pruning, about 30 calls per run for M = 6. At N = 2e5 that adds about 1.5–4 s to a run that took 262 s in the design's end-to-end test.

Two implementation choices in `adjust_pvals_families.cpp` halved the cost:
- an allocation-free insertion sort instead of `std::stable_sort`; it is still stable, so ties keep column order as `order()` does;
- a `no_init` output matrix.

## F2: integration

* **Full test suite.** Every `tests/testthat/test-*.R` file was run with `testthat::test_file()` against a private install: 444 `test_that()` blocks in 35 files, 0 failures, 0 errors, 0 warnings, and the 3 skips that were already there.
* **Existing snapshots unchanged.** New snapshot entries were added only for the families tests. One existing test was updated: `new_graph_optimal()` now lists the new `families` field.
* **Byte identity with `families = NULL`.** `main` and this branch were installed into separate libraries and the same seeded workloads were run in each. The outputs are `identical()` for:
  * `calc_power_pvals()`;
  * `create_obj_func()`, serial and 2 threads;
  * `prune_graph()`;
  * `graph_optimise()` with global and local search, serial and 2 threads: weights, matrix, power, solution, control, GA solution/fitness/summary, and nloptr solution/objective/iterations.
* **Families end to end** (m = 4, families {1}, {2,3}, {4}, nsim 1e4):
  * runs;
  * is reproducible with the same seed;
  * the reported power is `identical()` to `calc_power_pvals(families = )` on the returned graph;
  * `popSize` is calibrated on 3 nodes.
* **Lint and format.** `air format --check` passes on all R files. `lintr` is clean on the package and test files; the remaining lints are only script idioms (`library()`, `sink()`) in the dev scripts.

## F3: documentation

* `?family_spec` and the `families` parameter in `?graph_optimise` and `?calc_power_pvals` (`calc_power_pvals()` inherits it).
* `NEWS.md`, `_pkgdown.yml` (reference and the "Details" article list), `inst/WORDLIST`.
* Article `vignettes/articles/families.qmd`: a simplified ETHOS two-dose design, m = 12, M = 6. The optimised graph comes from `data-raw/families-vignette.R` (default control, 4 threads, about 14 minutes on the laptop; pruning took 1 s including the Option A re-adjustments) and is saved to `vignettes/articles/data/families-graph.rds`.
  * The article renders with Quarto against the dev build. `pkgdown::check_pkgdown()` reports no problems.
  * On an independent sample of 1e5 trials, the gain is 6.077 for the optimised family graph and 5.874 for the hierarchy as drawn.

## Deviations from the plan

* `check_families_dims()` has no leading dot, following the package's `check_*()` naming.
* `calc_power_pvals()` takes `families` before `call`. `graph_optimise()` takes it after `alpha`.
* The test helper's literal oracle is split into `alg57_oracle()`, `alg57_trial()` and `bretz_update()` to pass the cyclomatic-complexity lint.
