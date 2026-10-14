# Brief P1: spec class, lattice, eager tables, `resolve_pvalues_gsd()`

Branch `gsd-build-spending`, repo `/home/user/multigrain-spending`. You implement phase P1 of
`dev/hsd_opt_build_plan.md`. **Report, do not commit.** Write complete files, not diffs.

## Read first (~45 minutes; do not skip)

1. `dev/hsd_opt_build_plan.md` §1 (reconciliations — rows 1, 3, 4, 7 bind you), §4 "Shared test
   fixture" and "P1".
2. `dev/hsd_opt_design_revised.md` §1.1, §1.2, §2.1, §3 items 3–5 and 9, §5.1, §5.2, Appendix A.1.
3. `dev/implement_gsd.md` "Read first" and "Constraints" (house style: 4-space indent, rlang
   standalone checks in `R/check_types.R`, `rlang::check_dots_empty()`, optional args after `...`
   named, `cli::cli_abort(call = call)`, `new_*` constructors, roxygen markdown, `#' @export`).
4. `R/transform_pvalues_gsd.R` in full: `transform_pvalues_gsd()`, `.gsd_transform_hyp()`,
   `new_pvals_gsd()`, `print`/`summary` methods, `.gsd_boundary_table()`, `.gsd_invert()`,
   `.gsd_check_well_ordered()`, `.gsd_spend()`, `.gsd_check_spending()`, `.gsd_looks()`,
   `.gsd_warn_matured()`, `.gsd_info_frac()`, `.gsd_look_back()`, `.gsd_spending_list()`,
   `.gsd_spending_labels()`.
5. `R/check_gsd.R` (`check_pvals_gsd()` :10-30, `.gsd_kernel_matrix()` :180).
6. `tests/testthat/test-transform_pvalues_gsd.R` (fixture style, `withr::with_seed`,
   `expect_snapshot(error = TRUE)` usage) and `tests/testthat/helper-gsd_reference.R`.
7. `NEWS.md` top entries (naming and argument conventions).

## What to build

### New file `R/transform_pvalues_gsd_opt.R`

- `sf_hsd_opt(id = NULL, bounds = c(-4, 1), step = 0.25, start = NULL)` **exported**. Returns
  `function(alpha, t, param) gsDesign::sfHSD(alpha, t, param)$spend` with class
  `c("multigrain_sf_hsd_opt", "function")` and attributes `id` (character or NULL — normalise
  with `as.character()`; a NULL `id` stays NULL), `bounds`, `step`, `start`. Validation with
  `cli::cli_abort` messages that **name the limits**: `bounds` finite numeric length 2,
  `bounds[1] >= -10`, `bounds[2] <= 2`, `bounds[1] < bounds[2]`; `step` one of
  `0.5, 0.25, 0.125` (exact match after `as.numeric`); `start` NULL or a finite number inside
  `bounds` (snapped later); `id` NULL or a length-1 non-NA string/number.
- `is_sf_hsd_opt(x)`; `.has_hsd_opt(spending)` (TRUE if `spending` is a tagged function or a
  list containing one).
- `.hsd_node_fn(f, node)`: `function(a, t) f(a, t, param = node)` — the two-argument wrapper
  every fixed helper receives (plan row 1).
- `.hsd_groups(spending, m, call)`: from the length-m list, the group map: integer vector
  `group` (length m, `NA` for untagged entries, 1..d for tagged ones; entries sharing a non-NULL
  `id` share a group; `id = NULL` entries each get their own group — **but a single tagged
  function passed alone and recycled by `.gsd_spending_list()` is one group**, design change 6:
  detect it by `identical()` of the list elements when the original `spending` was a single
  function). Group members must have `identical()` `bounds`, `step`, `start`, else abort naming
  the hypotheses. Returns `list(group, names, bounds (d×2), step (length d), start (length d))`
  where `names` are the ids, else `"gamma"` for d=1, `"gamma1", "gamma2", …` for d≥2 (ids take
  precedence; mixing named and unnamed groups keeps each group's own name).
- `.hsd_lattice(bounds, step)`: the multiples of `step` inside `[lower, upper]` plus both bounds,
  sorted unique. Because `step` is a binary fraction and bounds are user doubles, compute nodes as
  `step * seq(ceiling(lower / step), floor(upper / step))` and `union` the bounds; `unique()`
  after `round(, 12)` is NOT allowed — nodes must be exact doubles (`k * step` is exact).
- `.hsd_snap(value, lattice)`: nearest node (ties to the lower node).
- `.hsd_check_invertible(tab, hyp, node, call)`: for each column of `tab$bounds`,
  `sum(!duplicated(col) & col > 0) >= 2`, else `cli_abort` naming the hypothesis, the analysis,
  the node value, with the hint to narrow `bounds` (plan row 1).
- `.hsd_node_tables(spec_inputs, …)`: for every distinct `(schedule key, node)` among active
  tagged hypotheses, one table: `tab <- .gsd_boundary_table(t_row[looks], spending =
  .hsd_node_fn(f, node), alpha, grid_size, hyp = i)` inside `rlang::try_fetch(..., error =
  function(cnd) cli::cli_abort(c("Hypothesis {i} is not well ordered at {name} = {node}.",
  i = "Narrow {.arg bounds} of {.fn sf_hsd_opt}."), parent = cnd, call = call))`, then
  `.hsd_check_invertible()`. The schedule key is `paste(format(t_row[looks], digits = 15),
  collapse = ",")` together with `alpha` and `grid_size` (both constant within a spec, so the
  key is the schedule). Tables are stored once in a named list and referenced by key from each
  hypothesis. Also run `.gsd_check_spending(.hsd_node_fn(f, node), alpha, hyp = i)` at every
  node (cheap).
- `.gsd_invert_hyp(raw, t_row, tab, look_back)`: steps 4–7 of `.gsd_transform_hyp()` on a
  prebuilt table — the `.gsd_invert()` loop per look, matured copy-forward, `NA` looks = 1,
  `cummin` under look-back. **Duplicate the code** from `.gsd_transform_hyp()`; do not refactor
  that function. Returns `list(values, table = c(list(looks, maturity), tab))` with the same
  shape as `.gsd_transform_hyp()`.
- `.transform_pvalues_gsd_spec(pvals, info_frac, look_back, spending, alpha, grid_size,
  spending_label, call)`: called from `transform_pvalues_gsd()` after `.gsd_spending_list()`.
  It (1) runs `.gsd_check_spending()` on every **fixed** entry exactly as the fixed loop does;
  (2) transforms every fixed hypothesis with `.gsd_transform_hyp()` exactly as today and stores
  `values` and `table`; (3) applies the short-circuit rule to tagged hypotheses
  (`length(looks) == 1L && t_row[[looks]] >= 1`, from `.gsd_looks()`): such a hypothesis is
  transformed with `.gsd_transform_hyp()` using `.hsd_node_fn(f, start node)` (the table is NULL
  there anyway) and dropped from the active set with `group[i] <- NA`; a group left empty
  disappears (renumber groups 1..d); (4) if `d == 0` after this, returns the plain
  `multigrain_pvals_gsd` that the fixed path would build (`new_pvals_gsd()` with
  `.gsd_spending_labels()`-style labels, see `.hsd_labels()` below) and emits
  `cli::cli_inform()` saying every tagged hypothesis was short-circuited; (5) builds lattices,
  tables; (6) returns `new_pvals_gsd_spec(...)`. **Draws no random numbers.**
- `new_pvals_gsd_spec(...)`: fields `raw` (the input array — **never** a field called
  `pvals`), `nsim`, `m`, `K`, `alpha`, `info_frac`, `look_back`, `grid_size`, `fixed` (list of
  m: an `nsim × K` matrix for fixed/short-circuited hypotheses, NULL for active ones),
  `fixed_tables` (list of m: the `table` element for fixed/short-circuited, NULL otherwise),
  `param` (list: `names`, `bounds` d×2 with rownames = names, `step` length d, `start` length
  d (snapped nodes), `lattice` list of d numeric vectors), `group` (length m integer, NA for
  inactive), `schedule_key` (length m character, NA for inactive), `node_tables` (named list
  keyed `"<schedule_key>|<node>"`), `spending` (length-m labels from `.hsd_labels(…, spec =
  TRUE)`), `short_circuit` (logical length m), `spending_label`. Class
  `"multigrain_pvals_gsd_spec"` **only**.
- `resolve_pvalues_gsd(spec, spending_param)` **exported**: `spending_param` numeric of length
  d (named or not; if named, names must match `spec$param$names` in any order); each value must
  be on its lattice (`abs(v - node) <= 1e-12` for some node, else abort naming the parameter,
  the value and the nearest nodes); for each active hypothesis take the table for
  `(schedule_key, node)` and run `.gsd_invert_hyp()`; assemble the `nsim × m × K` array with the
  fixed matrices; return `new_pvals_gsd(pvals = out, nsim, m, n_look = K, alpha, info_frac,
  look_back, spending = .hsd_labels(…, values = spending_param), tables = …)` where `tables[[i]]`
  for an active hypothesis is `c(list(looks, maturity), tab)` exactly as the fixed path builds it.
  Keep `dimnames(pvals)` from `raw`. The result must be `identical()` to the fixed path's
  (gate below) — compare field by field while developing.
- `.hsd_labels(fixed_labels, group, names, bounds, values = NULL)`: length-m character. Fixed
  (and short-circuited) hypotheses: the label the fixed path gives — compute
  `.gsd_spending_labels(spending, label = spending_label)` **on the original list** once in the
  spec builder and keep its entries for those hypotheses; a short-circuited tagged hypothesis
  reads `sfHSD(param = <start node>)`. Active hypotheses: `sprintf("sfHSD(param = %s)",
  format(value, digits = 15))` on a resolved object; `sprintf("sfHSD(param = %s in [%s, %s])",
  name, lower, upper)` on the spec.
- `print.multigrain_pvals_gsd_spec(x, ...)` **exported** (`@export` + `@rdname` of
  `resolve_pvalues_gsd` or its own): class line; dimensions as the fixed `print` shows them;
  one line per hypothesis `H<i>: <label>` with `[short-circuited: one look at full
  information]` where applicable; one line per parameter `<name>: bounds [lo, hi], step s,
  M nodes, start v`; `tables built: n`. Match the fixed `print`'s `cli::cat_line` style.

### Touchpoints (each a few lines, guarded)

- `R/transform_pvalues_gsd.R`: between the `.gsd_spending_list()` call and the
  `.gsd_check_spending()` loop insert
  ```r
  if (.has_hsd_opt(spending)) {
      return(.transform_pvalues_gsd_spec(
          pvals, info_frac = info_frac, look_back = look_back, spending = spending,
          alpha = alpha, grid_size = grid_size, spending_label = spending_label
      ))
  }
  ```
  and give `.gsd_spending_list()` an argument `allow_opt = FALSE` after `call`: when FALSE and
  any element `is_sf_hsd_opt()`, that is fine **only** because the handoff above catches it
  first — so instead make the check in `.gsd_spend()`: if `inherits(spending,
  "multigrain_sf_hsd_opt")` abort with "an optimised spending function reached the fixed
  transform; this is a bug in multigrain" (`.internal = TRUE`). Keep the edit to
  `.gsd_spending_list()` out entirely if it is not needed; the goal is that a tagged function can
  never be called with two arguments silently.
- `R/check_gsd.R` `check_pvals_gsd()`: inside the `!missing` branch, before
  `stop_input_type()`: if `inherits(pvals, "multigrain_pvals_gsd_spec")`, abort with
  `"{.arg {arg}} is a spending spec, not a transformed p-value object."` and the hint
  `"Pass it to {.fn graph_optimise_gsd} to optimise the spending parameters, or call {.fn
  resolve_pvalues_gsd} with chosen values."`, `call = call`.

Nothing else under `R/` changes. `NAMESPACE` and `man/` are regenerated by
`LC_ALL=C.UTF-8 Rscript -e 'devtools::document()'` — run it, and run it twice (second a
no-op). Never hand-edit `man/`.

### Tests

- `tests/testthat/helper-hsd_opt.R`: the shared fixture of the plan (§4): seed-7 raw draws
  (m=3, K=2, nsim 2000, `corr_matrix = diag(3)`, `info_frac = c(0.5, 1)`); `hsd_opt_spec()`
  building the d=1 spec (`spending = list(gsDesign::sfLDOF, sf_hsd_opt(id = "sec", bounds =
  c(-4, -2), step = 0.5), sf_hsd_opt(id = "sec", bounds = c(-4, -2), step = 0.5))`,
  `grid_size = 64L`); `hsd_fixed_fn(g)`; the d=2 fixture (m=4, K=2, nsim 1000, two groups
  with 3-node lattices `bounds = c(-3, -2), step = 0.5`, one short-circuited hypothesis via an
  `info_frac` matrix with a `c(NA, 1)` row); the tiny GA control from
  `test-optimisation_gsd.R:23-26` for later phases. Build fixtures lazily in functions so
  `test_file()` on other files does not pay for them.
- `tests/testthat/test-transform_pvalues_gsd_opt.R`, the gate checklist of the plan §4 P1,
  every item as a `test_that()`. Snapshots for the abort messages
  (`expect_snapshot(error = TRUE)`), in a new `_snaps/transform_pvalues_gsd_opt.md`. m ≤ 4,
  nsim ≤ 2000, grid_size ≤ 128 everywhere.

## Gate (run yourself before reporting; the orchestrator re-runs it from a clean build)

```r
Sys.setenv(NOT_CRAN = "true"); devtools::load_all()
testthat::test_file("tests/testthat/test-transform_pvalues_gsd_opt.R", package = "multigrain")
testthat::test_file("tests/testthat/test-transform_pvalues_gsd.R", package = "multigrain")   # unchanged count
testthat::test_file("tests/testthat/test-check_gsd.R", package = "multigrain")               # unchanged count
testthat::test_file("tests/testthat/test-optimisation_gsd.R", package = "multigrain")        # unchanged count
lintr::lint("R/transform_pvalues_gsd_opt.R"); lintr::lint("R/transform_pvalues_gsd.R")
lintr::lint("R/check_gsd.R"); lintr::lint("tests/testthat/test-transform_pvalues_gsd_opt.R")
lintr::lint("tests/testthat/helper-hsd_opt.R")
```
`devtools::load_all()` is allowed for you (it compiles in place); **do not** run
`devtools::install()` or `clean_dll()`: the installed library is shared with a measurement
agent. Record the pass/fail/skip counts. If the fixed-path counts change, stop and report.

## Report (your final message)

Files created/modified (full paths); the gate output with counts; any place you diverged from
this brief or the plan and why; anything you could not verify; anything the design did not
anticipate (for the progress entry). Do not commit. Do not touch `dev/` except nothing — the
orchestrator writes the review files.
