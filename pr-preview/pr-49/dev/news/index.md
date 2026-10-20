# Changelog

## multigrain (development version)

### New functionality

- [`stopping_rule_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/stopping_rule_gsd.md)
  creates an efficacy stopping rule for a group sequential design: one
  Boolean condition over the rejection indicators `r1, r2, ...` per
  interim analysis (or a single condition applied at every interim).
  Pass it to
  [`graph_optimise_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/graph_optimise_gsd.md)
  or
  [`calc_power_pvals_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/calc_power_pvals_gsd.md)
  through their new `stopping_rule` argument and each simulated trial
  stops at the first interim analysis whose condition holds, so the gain
  and every reported power count only the rejections made before the
  trial stopped.
  [`calc_power_pvals_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/calc_power_pvals_gsd.md)
  also reports `stop_distribution`, the proportion of trials stopping at
  each analysis, and `mean_stop_look`. A trial that rejects every
  hypothesis stops at that analysis with or without a rule. Stopping
  only removes rejections, so error control is unaffected.
- [`sf_hsd_opt()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/sf_hsd_opt.md)
  tags a hypothesis whose Hwang-Shih-DeCani spending parameter is to be
  chosen by the optimiser rather than fixed in advance. Hypotheses with
  the same `id` share a parameter; `bounds` (default `c(-4, 1)`, within
  -10 to 2) and `step` (0.5, 0.25 or 0.125) set the lattice of values
  searched. Given one or more such entries,
  [`transform_pvalues_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/transform_pvalues_gsd.md)
  returns a `multigrain_pvals_gsd_spec` object holding the raw p-values
  and a boundary table for every lattice node, each checked for well
  ordering when the spec is built, so every value the search can try is
  a valid transform.
- [`resolve_pvalues_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/resolve_pvalues_gsd.md)
  turns a spec and a value per spending parameter into the
  `multigrain_pvals_gsd` object
  [`transform_pvalues_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/transform_pvalues_gsd.md)
  would have returned with those values fixed, from the stored tables.
  Use it to evaluate a fixed choice, to reproduce an optimised result’s
  final evaluation, or to check a chosen design on fresh simulated
  p-values.
- [`graph_optimise_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/graph_optimise_gsd.md)
  accepts a spec and optimises the spending parameters jointly with the
  graph: a genetic algorithm with one gene per parameter, then a COBYLA
  refinement of the graph at every lattice node of each parameter,
  warm-started from its neighbour. The result adds `spending_param`,
  `hyp_spending_param`, `spending`, `spending_bounds`, `spending_step`,
  `spending_group`, `spending_profile` (the gain at each node with the
  graph re-optimised and with it fixed, and the paired difference from
  the optimum with its standard error), `spending_search` and
  `profile_converged`; [`print()`](https://rdrr.io/r/base/print.html)
  shows the spending functions and the range of values within 5e-4 of
  the best gain, and [`summary()`](https://rdrr.io/r/base/summary.html)
  the profile. A result without optimised spending is unchanged.

### Bug fixes

- [`graph_optimise_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/graph_optimise_gsd.md)
  and
  [`calc_power_pvals_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/calc_power_pvals_gsd.md)
  now check that the graph has as many hypotheses as the p-values and
  abort with a message naming both numbers; previously the mismatch
  surfaced only as the internal error “G must be an m x m matrix”.
- [`graph_optimise()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/graph_optimisation.md)
  and
  [`calc_power_pvals()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/calc_power_pvals.md)
  now refuse a gain function created with
  [`trial_success_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/trial_success_gsd.md),
  which needs decision times; previously the object was accepted and
  every rejection was silently scored as an analysis-1 rejection. Use
  [`graph_optimise_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/graph_optimise_gsd.md)
  and
  [`calc_power_pvals_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/calc_power_pvals_gsd.md).

## multigrain 1.0.0

CRAN release: 2026-10-18

- Initial release.
- Added a NEWS.md file to track package changes.
