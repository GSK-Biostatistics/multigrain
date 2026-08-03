# Modify global optimisation options

`control_global()` is for expert use only; it allows you to directly set
[`GA::ga()`](https://github.com/luca-scr/GA/reference/ga.html) options
to access features that are otherwise not available in multigrain.

## Usage

``` r
control_global(.ctrl, ...)
```

## Arguments

- .ctrl:

  A
  [multigrain_control](https://gsk-biostatistics.github.io/multigrain/dev/reference/multigrain_control.md)
  object.

- ...:

  \<[`dynamic-dots`](https://rlang.r-lib.org/reference/dyn-dots.html)\>
  Name-value pairs of
  [`GA::ga()`](https://github.com/luca-scr/GA/reference/ga.html) options
  and their values.

## Value

A modified
[multigrain_control](https://gsk-biostatistics.github.io/multigrain/dev/reference/multigrain_control.md).

## Examples

``` r
# `control_global()` allows you to access `GA::ga()` options that are not
# otherwise exposed by multigrain. For example, in special cases you may
# need to control the population size (the number of candidate graphs
# evaluated at each generation).
# multigrain makes some informed choices, but if you're convinced you want to
# try other values, you can access this (and other) [GA::ga()] options:
multigrain_control() |>
    control_global(pcrossover = 0.2)
#> <multigrain_control>
#> global optimisation:
#> • pcrossover: 0.2
```
