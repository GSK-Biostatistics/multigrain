# Modify the number of simulations

`control_nsim_local()` allows you to set the number of local
simulations.

`control_nsim_global()` allows you to set the number of global
simulations.

## Usage

``` r
control_nsim_local(ctrl, nsim_local)

control_nsim_global(ctrl, nsim_global)
```

## Arguments

- ctrl:

  A
  [multigrain_control](https://gsk-biostatistics.github.io/multigrain/dev/reference/multigrain_control.md)
  object.

- nsim_local:

  The number of simulations to use when evaluating the trial success
  function in local optimisation.

  - If set and lower than the number of rows in the `pvals` matrix,
    [`graph_optimise()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/graph_optimise.md)
    will use a random sample of `nsim_local` rows from `pvals` for local
    optimisation.

  - If unset or greater than or equal to the number of rows in the
    `pvals` matrix, *all* rows from `pvals` will be used.

- nsim_global:

  The number of simulations to use when evaluating the trial success
  function in global optimisation. If unset, the minimum between 50000
  and the number of sets of p-values will be used.

## Value

A modified
[multigrain_control](https://gsk-biostatistics.github.io/multigrain/dev/reference/multigrain_control.md).

## Examples

``` r
multigrain_control() |>
    control_nsim_local(10000)
#> <multigrain_control>
#> local simulations: 10000

multigrain_control() |>
    control_nsim_global(10000)
#> <multigrain_control>
#> global simulations: 10000
```
