# Set parameters for graph optimisation

There are two steps needed to fine tune the graph optimisation with
multigrain:

1.  Create a control object with `multigrain_control()`.

2.  Define its behaviour with `control_` functions:

    - [`control_nsim_local()`](https://gsk-biostatistics.github.io/multigrain/reference/control_nsim_local.md)
      to set the number of local simulations.

    - [`control_nsim_global()`](https://gsk-biostatistics.github.io/multigrain/reference/control_nsim_local.md)
      to set the number of global simulations.

    - [`control_local()`](https://gsk-biostatistics.github.io/multigrain/reference/control_local.md)
      to modify the local optimisation options.

    - [`control_global()`](https://gsk-biostatistics.github.io/multigrain/reference/control_global.md)
      to modify the global optimisation options.

Any unset parameters will automatically be set before running the
optimisation. There are predefined defaults, but they are calibrated
based on the `pvals` dimensions.

## Usage

``` r
multigrain_control()
```

## Value

A multigrain *control* object.

## Examples

``` r
multigrain_control()
#> <multigrain_control>
```
