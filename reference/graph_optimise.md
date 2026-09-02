# Optimise graph-based multiple testing procedures

Optimise the hypothesis weights and transition matrix of a graph-based
multiple testing procedure using
[`nloptr::nloptr()`](https://astamm.github.io/nloptr/reference/nloptr.html)
local optimisation and, optionally, a genetic algorithm using
[`GA::ga()`](https://github.com/luca-scr/GA/reference/ga.html) for
global optimisation.

## Usage

``` r
graph_optimise(
  pvals,
  graph_constraint,
  trial_success,
  ...,
  alpha = 0.025,
  start_graph = list(list(hyp_weight = NULL, trans_matrix = NULL)),
  global_search = TRUE,
  num_threads = 1L,
  control = multigrain_control(),
  verbose = multigrain_verbosity()
)
```

## Arguments

- pvals:

  A numeric matrix of p-values. Each row represents a simulated trial,
  and each column corresponds to a hypothesis.

- graph_constraint:

  A `multigrain_graph_constraint` object containing constraints on the
  graph's weights and transition matrix. Created with
  [`graph_constraint()`](https://gsk-biostatistics.github.io/multigrain/reference/graph_constraint.md).

- trial_success:

  A `multigrain_trial_success` object defining the trial success measure
  (power objective) function to maximize. Created with
  [`trial_success()`](https://gsk-biostatistics.github.io/multigrain/reference/trial_success.md).

- ...:

  These dots are for future extensions and must be empty.

- alpha:

  A single numeric value representing the overall one-sided significance
  level. Default is 0.025.

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
  [`multigrain_control()`](https://gsk-biostatistics.github.io/multigrain/reference/multigrain_control.md).

- verbose:

  An optional string controlling verbosity (`"detail"` \> `"info"` \>
  `"silent"`). Verbosity can also be set at package level with the
  `multigrain_verbosity` option (see
  [`multigrain_verbosity()`](https://gsk-biostatistics.github.io/multigrain/reference/multigrain_verbosity.md)):

  - `"info"` (default): will only show milestones / informational
    messages highlighting the progress of the optimisation at
    coarse-grained level.

  - `"detail"`: will show milestones and information about fine-grained
    optimisation events.

  - `"silent"`: no information about the progress of the optimisation is
    printed to the console. Errors and warnings are still thrown
    normally.

## Value

A `multigrain_graph_optimal` object containing: \* `hyp_weight`:
Optimised hypothesis weights (numeric vector). \* `trans_matrix`:
Optimised transition matrix (numeric matrix). \* `constraints`: List of
constraints used in the optimisation for weights and transition matrix.
\* `trial_success`: The trial success function used in the optimisation.
\* `power`: Power metrics for the optimised graph. \* `solution`: A list
containing: \* `opt_source`: Source of the optimal solution (`local` or
`global`). \* `graph_valid`: Named logical vector indicating validity of
the local and global solutions
(`c("local" = TRUE/FALSE, "global" = TRUE/FALSE)`). \* `global_search`:
`TRUE` or `FALSE` indicating whether a global optimisation was
performed. \* `control`: A modified
[`multigrain_control()`](https://gsk-biostatistics.github.io/multigrain/reference/multigrain_control.md)
object used. The values passed on by the user are complemented with
contextual defaults. \* `global_output`: Output from the genetic
algorithm if global optimisation was performed. \* `local_output`:
Output from the NLOPT optimisation. \* `start_graph`: Initial starting
values used in the optimisation.

## Details

The output is a graph where a specified objective function - the *trial
success measure* - is maximised under given constraints on the graph
structure, conditional on a p-value distribution supplied.

## Examples

``` r

# Generate test data
pvals <- simulate_pvalues(
  power_nominal = c(0.9, 0.85, 0.8, 0.75),
  corr_matrix = diag(4),
  nsim = 5000
)

# Create trial success function
ts <- trial_success(r1 + r2 + r3 + r4)
#> ✔ Trial success function compiled and sourced successfully.

# \donttest{
# Optimise graph
result <- graph_optimise(
  pvals = pvals,
  graph_constraint = graph_constraint_free(4),
  trial_success = ts,
  num_threads = 2
)
#> ℹ Running global optimization
#> ✔ Running global optimization [26.8s]
#> 
#> ℹ Evaluating trial success of globally optimised graph
#> ✔ Evaluating trial success of globally optimised graph [15ms]
#> 
#> ℹ Running local optimization
#> ✔ Running local optimization [87ms]
#> 
#> ℹ Evaluating trial success of locally optimised graph
#> ✔ Evaluating trial success of locally optimised graph [14ms]
#> 
#> ℹ Pruning redundant weights and edges
#> ✔ Pruning redundant weights and edges [20ms]
#> 
#> ℹ Evaluating trial success of pruned graph
#> ✔ Evaluating trial success of pruned graph [7ms]
#> 
# }
```
