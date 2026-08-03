# Modify local optimisation options

`control_local()` is for expert use only; it allows you to directly set
[`nloptr::nloptr()`](https://astamm.github.io/nloptr/reference/nloptr.html)
options to access features that are otherwise not available in
multigrain.

## Usage

``` r
control_local(.ctrl, ...)
```

## Arguments

- .ctrl:

  A
  [multigrain_control](https://gsk-biostatistics.github.io/multigrain/reference/multigrain_control.md)
  object.

- ...:

  \<[`dynamic-dots`](https://rlang.r-lib.org/reference/dyn-dots.html)\>
  Name-value pairs of
  [`nloptr::nloptr()`](https://astamm.github.io/nloptr/reference/nloptr.html)
  options and their values.

## Value

A modified
[multigrain_control](https://gsk-biostatistics.github.io/multigrain/reference/multigrain_control.md).

## Examples

``` r
# `control_local()` allows you to access `nloptr::nloptr()` options that are
# not otherwise exposed by multigrain. For example, in special cases you
# may want to try a different local optimisation algorithm. multigrain uses
# "NLOPT_LN_COBYLA", but if you're convinced you want to try other
# algorithms, you can access this (and other) `nloptr::nloptr()` options:
multigrain_control() |>
    control_local(algorithm = "NLOPT_LN_NEWUOA")
#> <multigrain_control>
#> local optimisation:
#> • algorithm: "NLOPT_LN_NEWUOA"
```
