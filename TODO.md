# TODO

This text is not a task.

## Post-OSS-release

clean-up TODO.md

hexsticker

close(?) all remaining issues on internal site and capture them as TODOs

contributing guidelines (CONTRIBUTING.md)

issues templates (revisit what we currently have)

Add GitHub repo topics

Address [community
standards](https://github.com/GSK-Biostatistics/multigrain/community)

Code of conduct

Security policy

PR template

Repository admins accept content reports

## Pre-CRAN-release

## Bugs

Work on the website ~3d \#feat @john 2020-03-20

Fix the homepage ~1d \#bug @jane

## Improvements

Add option to generate random hypothesis weights or transition matrix
rows using a Dirichlet approach to simplex sampling.

Further testing for `R\random_graph.R` once Dirichlet simplex approach
explored.

plot.multigrain_graph_optimal is very complex. We could / should look
into having a [`format()`](https://rdrr.io/r/base/format.html) function
that handles the preprocessing

[`calc_power_pvals()`](https://gsk-biostatistics.github.io/multigrain/reference/calc_power_pvals.md)
: `custom_power` is somewhat similar to rlang’s list of quosures
(e.g. ensuring every element of the list is named and of the expected
type). Maybe we could follow a similar approach.

continue using param documentation inheritance

constructors for `multigrain_graph` and `multigrain_graph_optimal`, with
`multigrain_graph` as a parent and `multigrain_graph_optimal` to inherit
(e.g. plot method)

`graph_optimal` keeps the core elements and everything else gets moved
to attributes.

revisit `create_start_params()` when documenting `start_graph` design.

do not inform about the diagnosis mode when running with
`diagnosis = TRUE`.

- e.g. in tests
  `gc[["trans_constraint", tolerance = 0, diagnose = TRUE]] <- tolerance_tc`

revisit the idea of redesigning the optimisation interface, e.g., the
{httr2} approach:

- 2 states: unoptimised and optimised.
- create an unoptimised graph object. this object would collect all the
  characteristics of the unoptimised graph
- have a series of functions to operate on it (to modify its various
  components)
- a single optimisation function responsible for taking the object from
  the unoptimised state into the optimised one under the hood we could
  still keep the various functions, but the user will not need to deal
  with optimise_graph(), optimise_N(), etc.
- this should reduce some duplication in code between the several
  optimisation functions.
- should also make it easier for users since there is only 1
  optimisation function.

revisit the idea of redesigning `graph_constraint`

- [`graph_constraint()`](https://gsk-biostatistics.github.io/multigrain/reference/graph_constraint.md)
  needs at least `m` (which means `graph_constraint(4)` will be the
  equivalent of `graph_constraint_free(4)`)

  ``` r
  graph_constraint(4) |>
  graph_constraint_hyp(c(NA, NA, NA,0)) |>
  graph_constraint_trans(...) |>
  graph_constraint_prepare() |>
  ...
  ```

revisit cyclomatic complexity (the linter is set at 22 for now, reduce
it to 15)

## Performance

Port benchmarking suite

## Plotting

graph_constraint plotting: 3 colours (0, 1 and unconstrained)

design a better interface for custom plot layouts

- Some graphs (e.g. Ferber et al. 2011 in the test suite) need node
  positions nudged beyond what
  [`igraph::layout_as_tree()`](https://r.igraph.org/reference/layout_as_tree.html)
  gives. Today the only workaround is mutating the ggplot data post-hoc
  (`plot$data$x[1] <- …`).
- Let users supply a `layout` (`data.frame` / `list`) instead. Open
  design question: make `layout` an arg on
  [`graph_constraint()`](https://gsk-biostatistics.github.io/multigrain/reference/graph_constraint.md)
  so it persists through to `graph_optimal` and keeps the two plots
  consistent?

show epsilon edges in graph_optimal plotting

- `param_to_solution(process = TRUE)` snaps small values to `0.001` but
  doesn’t flag which entries were epsilon. Re-normalisation shifts the
  value away from exactly `0.001`, so the flag must be carried
  explicitly (e.g. a logical matrix in the `solution` slot).
- Update the plot method to render epsilon edges as dotted arrows via a
  `linetype` aesthetic on `geom_edge_fan()`. (`arrow` is not mappable
  per-edge — varying the arrowhead would need two filtered layers.)

graph plotting (inc. graph constraints) in grid format, circle, linear

- Offer alternative node arrangements beyond the current tree layout:
  grid, circle and linear.
- Investigate whether we can let the user control the arc / curvature of
  individual arrows.

`create_layout()` masks
[`ggraph::create_layout()`](https://ggraph.data-imaginist.com/reference/ggraph.html)

## Quality of life (aka spring cleaning)

switch back to `expect_equal()`. Follow Hadley’s advice: “Numerical
precision can also vary across platforms, so use expect_equal() unless
you have a specific reason for using expect_identical().”

investigate
[{bundle}](https://rstudio.github.io/bundle/articles/bundle.html).

explore use of GitHub copilot

printing `graph_optimal_random()` should only print `hyp_w` and
`trans_m`.

## Documentation

Write documentation for start_graph; e.g. a flowchart / excalidraw and
user workflows.

MCP server for {multigrain} alongside vignettes?

function documentation:

check all args are documented as sentences. They should start with a
capital letter and end in a period / full stop.

arguments:

`call` in
[`calc_power_pvals()`](https://gsk-biostatistics.github.io/multigrain/reference/calc_power_pvals.md):
this is an internal rlang plumbing argument and should not be exposed in
the user-facing `@param` documentation.

`sum_to_one_constraint` in
[`is_graph_valid()`](https://gsk-biostatistics.github.io/multigrain/reference/is_graph_valid.md):
the description uses a double negative (“both the transition matrix rows
are not constrained to sum-to-one”). Reword, e.g. “If `FALSE`, rows of
the transition matrix are not required to sum to 1 (e.g. for
fixed-sequence graphs).”

`@return` sections:

[`graph_random()`](https://gsk-biostatistics.github.io/multigrain/reference/graph_random.md):
returns a plain `list`, not an S3 object. Consider whether it should
return a named S3 class consistent with the rest of the package, or
document the intentional design choice explicitly. Should we have a
`multigrain_graph_random` class?

`@seealso` cross-references: add links to connect the natural user
workflow:

[`simulate_pvalues()`](https://gsk-biostatistics.github.io/multigrain/reference/simulate_pvalues.md)
→
[`calc_ncp()`](https://gsk-biostatistics.github.io/multigrain/reference/calc_ncp.md),
[`calc_power_pvals()`](https://gsk-biostatistics.github.io/multigrain/reference/calc_power_pvals.md),
[`graph_optimise()`](https://gsk-biostatistics.github.io/multigrain/reference/graph_optimise.md)

[`trial_success()`](https://gsk-biostatistics.github.io/multigrain/reference/trial_success.md)
→
[`graph_optimise()`](https://gsk-biostatistics.github.io/multigrain/reference/graph_optimise.md),
[`calc_power_pvals()`](https://gsk-biostatistics.github.io/multigrain/reference/calc_power_pvals.md)

[`graph_constraint()`](https://gsk-biostatistics.github.io/multigrain/reference/graph_constraint.md)
→
[`graph_constraint_free()`](https://gsk-biostatistics.github.io/multigrain/reference/graph_constraint_free.md),
[`graph_optimise()`](https://gsk-biostatistics.github.io/multigrain/reference/graph_optimise.md)

[`multigrain_control()`](https://gsk-biostatistics.github.io/multigrain/reference/multigrain_control.md)
→
[`control_nsim_local()`](https://gsk-biostatistics.github.io/multigrain/reference/control_nsim_local.md),
[`control_nsim_global()`](https://gsk-biostatistics.github.io/multigrain/reference/control_nsim_local.md),
[`control_local()`](https://gsk-biostatistics.github.io/multigrain/reference/control_local.md),
[`control_global()`](https://gsk-biostatistics.github.io/multigrain/reference/control_global.md)

check data object documentation, might have to switch to Rd (markdown
might not be supported)

use the `@family` tag to group function docs in the same manner as the
pkgdown website groups.

remove deprecated params? given this is going to be the first CRAN
release?

[`trial_success()`](https://gsk-biostatistics.github.io/multigrain/reference/trial_success.md)
vignette

pkgdown website preview for PRs

Add a non-technical introduction to the logic / business need / approach
behind {multigrain}

# Backlog (postponed tasks)

add a GHA workflow for formatting with Air

# Done ✓

add package-level verbosity control

make the `verbose` argument to
[`graph_optimise()`](https://gsk-biostatistics.github.io/multigrain/reference/graph_optimise.md),
[`trial_success()`](https://gsk-biostatistics.github.io/multigrain/reference/trial_success.md)
*et al.* an enumeration - `c("info", "detail", "silent")`

introduce `...` to separate required from optional arguments to relevant
functions:

[`calc_power_pvals()`](https://gsk-biostatistics.github.io/multigrain/reference/calc_power_pvals.md)

[`graph_constraint()`](https://gsk-biostatistics.github.io/multigrain/reference/graph_constraint.md)

[`graph_optimise()`](https://gsk-biostatistics.github.io/multigrain/reference/graph_optimise.md)

[`simulate_pvalues()`](https://gsk-biostatistics.github.io/multigrain/reference/simulate_pvalues.md)

[`normalise_sum()`](https://gsk-biostatistics.github.io/multigrain/reference/normalise_sum.md)

move the following functions and their tests to `R/graph_optimal` and
`tests/testthat/test-graph_optimal.R`, respectively.

`print.multigrain_graph_optimal()`

`summary.multigrain_graph_optimal()`

`is_graph_optimal()`

`check_graph_optimal()`

add a `multigrain_graph_optimal` constructor.

~~This task has been declined~~ (declined)

## Documentation for OSS release

functions documentation:

titles and descriptions:

added a description for
[`graph_optimal_get_control()`](https://gsk-biostatistics.github.io/multigrain/reference/graph_optimal_get_control.md)

merged the documentation for
[`control_nsim_global()`](https://gsk-biostatistics.github.io/multigrain/reference/control_nsim_local.md)
and
[`control_nsim_local()`](https://gsk-biostatistics.github.io/multigrain/reference/control_nsim_local.md)
using the same title and adding descriptions

added a description for
[`graph_constraint_free()`](https://gsk-biostatistics.github.io/multigrain/reference/graph_constraint_free.md)

added a description for
[`autoplot.multigrain_graph_constraint()`](https://gsk-biostatistics.github.io/multigrain/reference/autoplot.multigrain_graph_constraint.md)

added a description for
[`autoplot.multigrain_graph_optimal()`](https://gsk-biostatistics.github.io/multigrain/reference/autoplot.multigrain_graph_optimal.md)

arguments:

`global_search` in
[`graph_optimise()`](https://gsk-biostatistics.github.io/multigrain/reference/graph_optimise.md):
“A logical indicate whether” → “A logical indicating whether”.

`.ctrl` in
[`control_local()`](https://gsk-biostatistics.github.io/multigrain/reference/control_local.md)
/
[`control_global()`](https://gsk-biostatistics.github.io/multigrain/reference/control_global.md)
and `graph_optimal` in
[`graph_optimal_get_control()`](https://gsk-biostatistics.github.io/multigrain/reference/graph_optimal_get_control.md):
missing terminal period.

argument style: `root` and `digits` in plot methods use parenthetical
type annotations (`(integer-like)`, `Number`) instead of prose
sentences. Align with the rest of the package style.

use @inheritParams more (document important arguments only once).

`alpha`: documented in
[`graph_optimise()`](https://gsk-biostatistics.github.io/multigrain/reference/graph_optimise.md)
and inherited in
[`calc_ncp()`](https://gsk-biostatistics.github.io/multigrain/reference/calc_ncp.md),
[`calc_power_pvals()`](https://gsk-biostatistics.github.io/multigrain/reference/calc_power_pvals.md),
and
[`simulate_pvalues()`](https://gsk-biostatistics.github.io/multigrain/reference/simulate_pvalues.md).

aligned `verbose` in
[`graph_optimise()`](https://gsk-biostatistics.github.io/multigrain/reference/graph_optimise.md)
with `verbose` in
[`trial_success()`](https://gsk-biostatistics.github.io/multigrain/reference/trial_success.md).

[`calc_power_pvals()`](https://gsk-biostatistics.github.io/multigrain/reference/calc_power_pvals.md)
inherits the documentation for `trans_matrix` and `hyp_weight` from
[`is_graph_valid()`](https://gsk-biostatistics.github.io/multigrain/reference/is_graph_valid.md).

`pvals`: the
[`calc_power_pvals()`](https://gsk-biostatistics.github.io/multigrain/reference/calc_power_pvals.md)
arg inherits the documentation from
[`graph_optimise()`](https://gsk-biostatistics.github.io/multigrain/reference/graph_optimise.md).

~~`graph_constraint`: descriptions have different focus across
[`graph_optimise()`](https://gsk-biostatistics.github.io/multigrain/reference/graph_optimise.md)
and
[`graph_random()`](https://gsk-biostatistics.github.io/multigrain/reference/graph_random.md).
Consolidate with `@inheritParams`.~~

`tolerance`: partial consolidation of the arg description across
[`graph_constraint()`](https://gsk-biostatistics.github.io/multigrain/reference/graph_constraint.md),
[`is_graph_valid()`](https://gsk-biostatistics.github.io/multigrain/reference/is_graph_valid.md),
and
[`normalise_sum()`](https://gsk-biostatistics.github.io/multigrain/reference/normalise_sum.md).

`sum_to_one_constraint`: in
[`calc_power_pvals()`](https://gsk-biostatistics.github.io/multigrain/reference/calc_power_pvals.md)
the documented default did not match the actual default.

`@return` sections:

[`multigrain_control()`](https://gsk-biostatistics.github.io/multigrain/reference/multigrain_control.md):
“an S3 object (a list)” — the parenthetical is redundant. Simplify to “A
`multigrain_control` object.”

[`control_global()`](https://gsk-biostatistics.github.io/multigrain/reference/control_global.md)
/
[`control_local()`](https://gsk-biostatistics.github.io/multigrain/reference/control_local.md):
`@return` is lowercased (“a modified…”);
[`control_nsim_global()`](https://gsk-biostatistics.github.io/multigrain/reference/control_nsim_local.md)
uses a capital. Align.

`autoplot` methods: “`ggraph` / `ggplot` object” — ggraph inherits from
ggplot, so “`ggraph` object” is sufficient and more precise.

examples:

[`graph_random()`](https://gsk-biostatistics.github.io/multigrain/reference/graph_random.md):
uses explicit [`print()`](https://rdrr.io/r/base/print.html) calls —
idiomatic R examples just evaluate the expression directly.
