# Optimise graph-based multiple testing procedures for a group sequential design

Group sequential counterpart of
[`graph_optimise()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/graph_optimisation.md).
Optimise the hypothesis weights and transition matrix of a graph-based
multiple testing procedure for a design with several analyses, using
[`nloptr::nloptr()`](https://astamm.github.io/nloptr/reference/nloptr.html)
local optimisation and, optionally, a genetic algorithm using
[`GA::ga()`](https://github.com/luca-scr/GA/reference/ga.html) for
global optimisation.

## Usage

``` r
graph_optimise_gsd(
  pvals,
  graph_constraint,
  trial_success,
  ...,
  alpha = NULL,
  stopping_rule = NULL,
  start_graph = list(list(hyp_weight = NULL, trans_matrix = NULL)),
  global_search = TRUE,
  num_threads = 1L,
  control = multigrain_control(),
  verbose = multigrain_verbosity()
)
```

## Arguments

- pvals:

  A `multigrain_pvals_gsd` object of repeated (or, per hypothesis,
  sequential) p-values, created with
  [`transform_pvalues_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/transform_pvalues_gsd.md),
  or a `multigrain_pvals_gsd_spec` object, created by
  [`transform_pvalues_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/transform_pvalues_gsd.md)
  when one or more spending functions are
  [`sf_hsd_opt()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/sf_hsd_opt.md),
  to optimise those spending parameters with the graph (see Details).

- graph_constraint:

  A `multigrain_graph_constraint` object containing constraints on the
  graph's weights and transition matrix. Created with
  [`graph_constraint()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/graph_constraint.md).

- trial_success:

  A `multigrain_trial_success_gsd` object defining the trial success
  measure (gain function) to maximise. Created with
  [`trial_success_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/trial_success_gsd.md).
  A fixed-sample object from
  [`trial_success()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/trial_success.md)
  is not accepted; write a rejection-only gain such as `r1 + r2` with
  [`trial_success_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/trial_success_gsd.md)
  instead.

- ...:

  These dots are for future extensions and must be empty.

- alpha:

  A single numeric value representing the overall one-sided significance
  level, or `NULL` (the default) to use the level of `pvals`. A supplied
  value must be greater than `0` and at most the level the repeated
  p-values were built up to: above that they are capped at 1, so no
  allocation beyond it could ever reject.

- stopping_rule:

  A `multigrain_stopping_rule_gsd` object from
  [`stopping_rule_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/stopping_rule_gsd.md),
  or `NULL` (the default) for a trial that runs to its final analysis
  unless it rejects every hypothesis earlier.

- start_graph:

  Optional. Initial list of graphs suggested. Each graph is defined as a
  list containing starting weight vector `hyp_weight` and transition
  matrix `trans_matrix`. If `NULL`, default starting values are
  generated (currently with a Bonferroni-Holm graph).

- global_search:

  A logical indicating whether to perform a global optimisation before
  the local optimisation. Defaults to `TRUE`.

- num_threads:

  Number of threads to use for parallel execution of the shortcut
  algorithm. On shared systems (HPC clusters, login nodes), always
  explicitly set `num_threads` based on your resource allocation.
  Default is `1L` (serial execution).

- control:

  An optional `multigrain_control` object can be used to set various
  graph optimisation parameters. Created with
  [`multigrain_control()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/multigrain_control.md).

- verbose:

  An optional string controlling verbosity (`"detail"` \> `"info"` \>
  `"silent"`). Verbosity can also be set at package level with the
  `multigrain_verbosity` option (see
  [`multigrain_verbosity()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/multigrain_verbosity.md)):

  - `"info"` (default): will only show milestones / informational
    messages highlighting the progress of the optimisation at
    coarse-grained level.

  - `"detail"`: will show milestones and information about fine-grained
    optimisation events.

  - `"silent"`: no information about the progress of the optimisation is
    printed to the console. Errors and warnings are still thrown
    normally.

## Value

A `multigrain_graph_optimal` object containing:

- `hyp_weight`: Optimised hypothesis weights (numeric vector).

- `trans_matrix`: Optimised transition matrix (numeric matrix).

- `constraints`: List of constraints used in the optimisation for
  weights and transition matrix.

- `trial_success`: The trial success function used in the optimisation.

- `power`: The
  [`calc_power_pvals_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/calc_power_pvals_gsd.md)
  list for the optimised graph, so it carries power by analysis, the
  decision-time distribution and the value of the trial success measure
  as well as the fixed-sample power metrics.

- `solution`: A list containing:

  - `opt_source`: Source of the optimal solution (`local` or `global`).

  - `graph_valid`: Named logical vector indicating validity of the local
    and global solutions
    (`c("local" = TRUE/FALSE, "global" = TRUE/FALSE)`).

- `global_search`: `TRUE` or `FALSE` indicating whether a global
  optimisation was performed.

- `control`: A modified
  [`multigrain_control()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/multigrain_control.md)
  object used. The values passed on by the user are complemented with
  contextual defaults.

- `global_output`: Output from the genetic algorithm if global
  optimisation was performed.

- `local_output`: Output from the NLOPT optimisation.

- `start_graph`: Initial starting values used in the optimisation.

- `stopping_rule`: the rule, when one was supplied.

When `pvals` is a spending spec, also:

- `spending_param`: the chosen value of each spending parameter (named
  numeric vector; the names are the `id`s of
  [`sf_hsd_opt()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/sf_hsd_opt.md),
  else `gamma`, or `gamma1`, `gamma2`, ... when there are several).

- `hyp_spending_param`: the value per hypothesis, `NA` for a hypothesis
  with fixed spending or one whose single analysis at full information
  makes spending irrelevant.

- `spending`: a spending label per hypothesis, such as
  `sfHSD(param = -0.25)` for a tagged one.

- `spending_bounds`, `spending_step`: the bounds (a matrix with one row
  per parameter) and the lattice step of each parameter.

- `spending_group`: the parameter of each hypothesis, an integer index
  into `spending_param` (`NA` for none).

- `spending_profile`: the warm-started profile (see Details), a data
  frame with columns `param`, `value`, `gain`, `gain_fixed_graph`,
  `diff`, `se` and `n_eval` (the number of COBYLA evaluations at that
  node), and the graph of each row as the attribute `graphs`.

- `spending_search`: the genetic algorithm's coverage of each lattice, a
  data frame with columns `param`, `value`, `visits` and
  `best_fitness_subsample`, the best gain seen at that node on the
  algorithm's subsample, which is not comparable with the profile's
  gains. All `NA` without a global search.

- `profile_converged`: `FALSE` when, with two or more parameters, the
  rounds of the profile did not settle within three, or when a node of a
  final recomputed profile beats the chosen value (a negative `diff`);
  the plateaus are then indicative.

## Details

The output is a graph where a specified objective function - the *trial
success measure* - is maximised under given constraints on the graph
structure, conditional on the repeated p-values supplied. Because the
gain sees the analysis at which each hypothesis was rejected, an early
rejection can be worth more than a late one.

The search, the encoding and the pruning are those of
[`graph_optimise()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/graph_optimisation.md).
Only the evaluation changes: the group sequential kernel runs the
fixed-sample cascade once per analysis, carrying the updated graph
forward, and reports the analysis at which each hypothesis was declared
rejected. Those decision times are what `trial_success` scores.

With a stopping rule the kernel stops each simulated trial at the first
interim whose condition holds, so the gain and every reported power see
only the rejections made before the trial stopped; a trial that rejects
every hypothesis stops at that analysis in any case.

**Optimised spending parameters.** When `pvals` is a spending spec, from
[`transform_pvalues_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/transform_pvalues_gsd.md)
with one or more
[`sf_hsd_opt()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/sf_hsd_opt.md)
entries, the Hwang-Shih-DeCani parameter of each tagged group is chosen
jointly with the graph. Every value the search can try is a node of the
parameter's lattice (the multiples of `step` inside `bounds`), and every
node's boundary table was built and checked for well ordering when the
spec was created, so the chosen design controls the error rate as any
fixed choice would. The search has two stages:

- the global stage runs the genetic algorithm on the graph genes with
  one gene per spending parameter appended, each mapped to a lattice
  node, on the `nsim_global` subsample. The repeated p-values of that
  subsample are computed once at every node. The default population size
  counts each parameter as one more hypothesis;

- the local stage refines the graph by COBYLA at the genetic algorithm's
  parameter values, then walks each parameter's lattice outward from
  there and re-optimises the graph at every node, starting COBYLA from
  the best of the current optimum's graph, the previous node's graph and
  the best graph the genetic algorithm found at that node (a
  *warm-started profile*). It moves to the best node; with two or more
  parameters the sweep is repeated until no parameter moves, at most
  three rounds. Every gain of this stage is computed on one ordering of
  the simulated trials, so differences between nodes are paired.

The profile is returned as `spending_profile`: for every node of each
parameter, with the other parameters held at their chosen values, the
gain with the graph re-optimised (`gain`), the gain of the chosen graph
kept fixed (`gain_fixed_graph`), the paired difference from the optimum
(`diff`) and its batch-means standard error over 100 blocks of trials
(`se`). At the chosen value `gain` equals `gain_fixed_graph` and `diff`
and `se` are 0. [`print()`](https://rdrr.io/r/base/print.html) condenses
the profile to two *plateaus*: the connected run of nodes around the
chosen value whose `diff` is at most `5e-4`, and the same run for the
fixed-graph gain, measured from its largest value. A wide plateau means
the parameter hardly matters for this design within that range; the
fixed-graph plateau is narrower when the best graph depends on the
parameter. [`summary()`](https://rdrr.io/r/base/summary.html) prints the
table and the number of nodes the genetic algorithm visited.

As on the fixed path, the better of the two stages' graphs is kept, each
scored on the full sample at its own parameter values. If the genetic
algorithm's graph is kept, its parameter values are reported and the
profile is still the local stage's, which started from them. Unlike the
fixed path, a `start_graph` is honoured when `global_search = FALSE`.

The parameters and the graph are chosen on the same simulated trials
that score them, so the reported gain is optimistic (the winner's
curse), a little more so than for the graph alone. To check the design,
simulate fresh p-values, build a spec from them with the same
`spending`, and evaluate the chosen graph with
[`calc_power_pvals_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/calc_power_pvals_gsd.md)
on `resolve_pvalues_gsd(new_spec, result$spending_param)`.

The boundary tables are the main cost of building a spec: about 0.9 s
per lattice node at the default `grid_size = 1024` with two analyses,
about 20 s for the default 21-node lattice. `grid_size = 512` halves
this; in the package's measurements its largest relative interpolation
error against `grid_size = 1024` was `8.2e-6`. A spec can be reused
across gains and stopping rules on the same draws.

## See also

[`graph_optimise()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/graph_optimisation.md)
for the fixed-sample version,
[`transform_pvalues_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/transform_pvalues_gsd.md)
for the p-value object and
[`trial_success_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/trial_success_gsd.md)
for the gain.

## Examples

``` r
# two hypotheses, an interim at 50% information and a final analysis
corr <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)

set.seed(1)
raw <- simulate_pvalues_gsd(
    power_nominal = c(0.9, 0.8),
    corr_matrix = corr,
    info_frac = c(0.5, 1),
    nsim = 2000
)

pvals <- transform_pvalues_gsd(
    raw,
    spending = gsDesign::sfLDOF,
    grid_size = 128L
)

# a rejection at the final analysis is worth 75% of one at the interim
gain <- trial_success_gsd(
    r1 * d(t1) + r2 * d(t2),
    d = c(1, 0.75)
)
#> ✔ Trial success function compiled and sourced successfully.

# \donttest{
result <- graph_optimise_gsd(
    pvals = pvals,
    graph_constraint = graph_constraint_free(2),
    trial_success = gain,
    stopping_rule = stopping_rule_gsd(r1),
    num_threads = 2
)
#> ℹ Running global optimization
#> ✔ Running global optimization [7.2s]
#> 
#> ℹ Evaluating trial success of globally optimised graph
#> ✔ Evaluating trial success of globally optimised graph [14ms]
#> 
#> ℹ Running local optimization
#> ✔ Running local optimization [11ms]
#> 
#> ℹ Evaluating trial success of locally optimised graph
#> ✔ Evaluating trial success of locally optimised graph [13ms]
#> 
#> ℹ Pruning redundant weights and edges
#> ✔ Pruning redundant weights and edges [8ms]
#> 
#> ℹ Evaluating trial success of pruned graph
#> ✔ Evaluating trial success of pruned graph [6ms]
#> 

# choose the spending parameter of H2 jointly with the graph
spec <- transform_pvalues_gsd(
    raw,
    spending = list(gsDesign::sfLDOF, sf_hsd_opt(bounds = c(-4, 0))),
    grid_size = 128L
)
result_spending <- graph_optimise_gsd(
    pvals = spec,
    graph_constraint = graph_constraint_free(2),
    trial_success = gain,
    stopping_rule = stopping_rule_gsd(r1),
    num_threads = 2
)
#> ℹ Running global optimization
#> ✔ Running global optimization [9.5s]
#> 
#> ℹ Evaluating trial success of globally optimised graph
#> ✔ Evaluating trial success of globally optimised graph [15ms]
#> 
#> ℹ Running local optimization
#> ✔ Running local optimization [13ms]
#> 
#> ℹ Profiling spending parameters: gamma
#> ✔ Profiling spending parameters: gamma [157ms]
#> 
#> ℹ Evaluating trial success of locally optimised graph
#> ✔ Evaluating trial success of locally optimised graph [16ms]
#> 
#> ℹ Pruning redundant weights and edges
#> ✔ Pruning redundant weights and edges [8ms]
#> 
#> ℹ Evaluating trial success of pruned graph
#> ✔ Evaluating trial success of pruned graph [7ms]
#> 
result_spending$spending_param
#> gamma 
#>  -0.5 
result_spending
#> <multigrain_graph_optimal>
#> Optimal graph found (given user-defined constraints on graph and computational resources):
#> 
#> Hypothesis weights:
#>     H1     H2 
#> 0.2513 0.7487 
#> 
#> Transition matrix:
#>    H1 H2
#> H1  0  1
#> H2  1  0
#> 
#> Trial success function:
#> r1 * d(t1) + r2 * d(t2)
#> 
#> Spending functions:
#>   H1: gsDesign::sfLDOF
#>   H2: sfHSD(param = -0.5)   [optimised: gamma, bounds -4 to 0, step 0.25]
#> Gain profile in gamma (graph re-optimised at each value):
#>   within 5e-4 of the maximum on [-0.75, -0.5]; fixed-graph plateau [-0.5, -0.5]
#> 
#> Stopping rule:
#> Every interim: r1
#> Hypotheses referenced: up to r1
#> 
#> Value of trial success measure:
#> 1.3134
# }
```
