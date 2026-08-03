# Package index

## Inputs

Create the optimisation inputs

- [`simulate_pvalues()`](https://gsk-biostatistics.github.io/multigrain/reference/simulate_pvalues.md)
  : Simulate raw p-values

- [`trial_success()`](https://gsk-biostatistics.github.io/multigrain/reference/trial_success.md)
  : Create a trial success function

- [`graph_constraint()`](https://gsk-biostatistics.github.io/multigrain/reference/graph_constraint.md)
  :

  Create a *graph constraint* object for optimisation procedures

- [`graph_constraint_free()`](https://gsk-biostatistics.github.io/multigrain/reference/graph_constraint_free.md)
  :

  Create an unconstrained *graph constraint*

## Advanced control

Control more advanced aspects of the global and local optimisation

- [`multigrain_control()`](https://gsk-biostatistics.github.io/multigrain/reference/multigrain_control.md)
  : Set parameters for graph optimisation
- [`control_nsim_local()`](https://gsk-biostatistics.github.io/multigrain/reference/control_nsim_local.md)
  [`control_nsim_global()`](https://gsk-biostatistics.github.io/multigrain/reference/control_nsim_local.md)
  : Modify the number of simulations
- [`control_local()`](https://gsk-biostatistics.github.io/multigrain/reference/control_local.md)
  : Modify local optimisation options
- [`control_global()`](https://gsk-biostatistics.github.io/multigrain/reference/control_global.md)
  : Modify global optimisation options

## Optimisation

Optimisation functions

- [`graph_optimise()`](https://gsk-biostatistics.github.io/multigrain/reference/graph_optimise.md)
  : Optimise graph-based multiple testing procedures

## Post-processing

Work with an optimised graph

- [`calc_power_pvals()`](https://gsk-biostatistics.github.io/multigrain/reference/calc_power_pvals.md)
  : Calculate power for a graph-based multiple test procedure using
  p-values
- [`is_graph_valid()`](https://gsk-biostatistics.github.io/multigrain/reference/is_graph_valid.md)
  : Check the validity of a graph-based MTP
- [`graph_optimal_get_control()`](https://gsk-biostatistics.github.io/multigrain/reference/graph_optimal_get_control.md)
  : Get the optimisation control

## Plotting

Plotting functions

- [`autoplot(`*`<multigrain_graph_constraint>`*`)`](https://gsk-biostatistics.github.io/multigrain/reference/autoplot.multigrain_graph_constraint.md)
  [`plot(`*`<multigrain_graph_constraint>`*`)`](https://gsk-biostatistics.github.io/multigrain/reference/autoplot.multigrain_graph_constraint.md)
  :

  Autoplot method for `multigrain_graph_constraint` objects

- [`autoplot(`*`<multigrain_graph_optimal>`*`)`](https://gsk-biostatistics.github.io/multigrain/reference/autoplot.multigrain_graph_optimal.md)
  [`plot(`*`<multigrain_graph_optimal>`*`)`](https://gsk-biostatistics.github.io/multigrain/reference/autoplot.multigrain_graph_optimal.md)
  :

  Autoplot method for `multigrain_graph_optimal` objects

## Helper functions

Other useful functions

- [`calc_ncp()`](https://gsk-biostatistics.github.io/multigrain/reference/calc_ncp.md)
  : Calculate non-centrality parameter
- [`normalise_sum()`](https://gsk-biostatistics.github.io/multigrain/reference/normalise_sum.md)
  : Normalise graph weights to sum to a target value
- [`graph_random()`](https://gsk-biostatistics.github.io/multigrain/reference/graph_random.md)
  : Generate random graph

## Options

Options consulted by multigrain

- [`multigrain_verbosity()`](https://gsk-biostatistics.github.io/multigrain/reference/multigrain_verbosity.md)
  : Multigrain verbosity

## Data

multigrain objects used for examples and tests

- [`graph_optimal_example`](https://gsk-biostatistics.github.io/multigrain/reference/graph_optimal_example.md)
  : Optimised graph example
