# Brief P4: `print()` and `summary()` of a spending-optimised result

Branch `gsd-build-spending`, repo `/home/user/multigrain-spending`. You implement phase P4 of
`dev/hsd_opt_build_plan.md`. P1–P3 are committed. **Report, do not commit.** Complete files.

## Read first (~30 minutes)

1. `dev/hsd_opt_build_plan.md` §1 row 8 (binding: insertion points), §4 "P4".
2. `dev/hsd_opt_design_revised.md` §1.5 (the exact `print()` block and the `summary()`
   contents), §3 item 8.
3. `R/graph_optimal.R` in full — `print.multigrain_graph_optimal()` (:217-249) and
   `summary.multigrain_graph_optimal()` (:254-285); `summarise_power_object()`.
4. `R/optimisation_gsd_opt.R` (P3): the `spending` list the driver attaches (`spending_param`,
   `hyp_spending_param`, `spending`, `spending_bounds`, `spending_profile`, `spending_search`,
   `profile_converged`) and `.hsd_plateau()`.
5. `R/stopping_rule_gsd.R` `.gsd_rule_cat_detail()` and `R/trial_success_gsd.R`'s
   `summary` method — the house style for `cli::cat_line` blocks.
6. `tests/testthat/test-graph_optimal.R` (:147-231: how the snapshot objects are hand-built with
   `graph_optimal()`; the GSD block) and `tests/testthat/_snaps/graph_optimal.md`.

## What to build

- In `R/optimisation_gsd_opt.R` (append; P3's file): `.hsd_cat_print(x)` and
  `.hsd_cat_summary(object)`.
  - `.hsd_cat_print(x)` writes, in the fixed `cli::cat_line` style:
    ```
    Spending functions:
      H1: sfLDOF
      H2: sfHSD(param = -2.25)   [optimised: gamma, bounds -4 to 1, step 0.25]
      H3: sfHSD(param = -2.25)   [optimised: gamma]
    Gain profile in gamma (graph re-optimised at each value):
      within 5e-4 of the maximum on [-3.25, -1.25]; fixed-graph plateau [-2.5, -2]
    ```
    Hypothesis names: use `names(x$hyp_weight)` when set (the driver sets them from the
    constraint), else `H<i>`. The `[optimised: …, bounds …, step …]` suffix goes on the first
    member of each group, `[optimised: <name>]` on later members; short-circuited tagged
    hypotheses get `[fixed at <value>: one look at full information]`. One "Gain profile"
    block per parameter; both intervals from `.hsd_plateau()` on `spending_profile` (the
    re-optimised `diff` column, and a fixed-graph `diff` computed as `max(gain_fixed_graph) -
    gain_fixed_graph` — state that in a comment). If `profile_converged` is `FALSE`, add a line
    "profile did not converge in 3 rounds; the plateau is indicative".
  - `.hsd_cat_summary(object)`: a header "Spending-parameter profile (<name>):" then the
    `spending_profile` rows of that parameter as a data frame with columns `value, gain,
    gain_fixed_graph, diff, se, n_eval`, `round(, 4)` on the gains and `signif(, 2)` on `se`,
    printed with `print(df, row.names = FALSE)`; then "GA visited k of M nodes" from
    `spending_search` (`sum(visits > 0)`), or "no global search" when visits are all `NA`.
- `R/graph_optimal.R`: in `print()`, immediately **after** the "Trial success function" block
  and **before** the stopping-rule block: `if (!is.null(x$spending_param)) .hsd_cat_print(x)`;
  in `summary()`, **after** `summary(object$trial_success)` and **before** the stopping-rule
  summary: `if (!is.null(object$spending_param)) .hsd_cat_summary(object)`. No other change.
  The existing snapshot `tests/testthat/_snaps/graph_optimal.md` must stay byte-identical.

## Tests — new file `tests/testthat/test-graph_optimal_hsd_opt.R`
With its own snapshot `tests/testthat/_snaps/graph_optimal_hsd_opt.md`. Objects hand-built
with `graph_optimal()` (no optimiser run — construct `hyp_weight`, `trans_matrix`, a
`graph_constraint_free()`, a `trial_success_gsd()`, a `power` list with `trial_success`, and
the `spending` list by hand):
- [ ] d = 1: m = 3, H1 fixed `sfLDOF`, H2/H3 group `sec`, bounds `c(-4, 1)`, step 0.25, a
  21-row profile with a plateau from −3.25 to −1.25 and a narrower fixed-graph plateau;
  `print()` and `summary()` snapshots.
- [ ] d = 2: m = 4, two axes (`gamma1` 5 nodes, `b` 3 nodes), H4 short-circuited,
  `profile_converged = FALSE`, no global search (`visits` NA); both snapshots.
- [ ] `.hsd_plateau()` on a profile with a gap (`diff` dips above `tol` then returns below):
  only the connected run around the argmax.
- [ ] a result with a stopping rule prints the spending block **before** the stopping-rule
  block (assert on captured output order).
- [ ] `test-graph_optimal.R` still 33 and `git diff --quiet -- tests/testthat/_snaps/graph_optimal.md`.
- [ ] in `test-optimisation_gsd_opt.R`, add one test: `print()` and `summary()` of a real
  (tiny) spec result run without error (`expect_output`).

## Gate
```r
Sys.setenv(NOT_CRAN = "true"); devtools::load_all()
testthat::test_file("tests/testthat/test-graph_optimal_hsd_opt.R", package = "multigrain")
testthat::test_file("tests/testthat/test-graph_optimal.R", package = "multigrain")      # 33
testthat::test_file("tests/testthat/test-optimisation_gsd_opt.R", package = "multigrain")
testthat::test_file("tests/testthat/test-optimisation_gsd.R", package = "multigrain")   # 79
testthat::test_file("tests/testthat/test-plot_graph_optimal.R", package = "multigrain") # 44
lintr::lint("R/optimisation_gsd_opt.R"); lintr::lint("R/graph_optimal.R")
lintr::lint("tests/testthat/test-graph_optimal_hsd_opt.R")
```
Snapshot files must have LF endings. No `install()`/`clean_dll()`.

## Report
Files; counts; divergences; not verified; unanticipated. Do not commit.
