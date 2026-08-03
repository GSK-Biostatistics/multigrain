# Multigrain verbosity

`multigrain_verbosity()` returns the option named
`multigrain_verbosity`, which controls multigrain's verbosity. There are
3 possible levels ("detail" \> "info" \> "silent"):

- "info" (default): will only show milestones / informational messages
  highlighting the progress of the optimisation at coarse-grained level.

- "detail": will show milestones and information about fine-grained
  optimisation events.

- "silent": no information about the progress of the optimisation is
  printed to the console. Errors and warnings are still thrown normally.

## Usage

``` r
multigrain_verbosity()
```

## Value

A string indicating the verbosity level as set with options or "info" if
`multigrain_verbosity` is unset.

## Details

If the `multigrain_verbosity` option is unset, then the "info" level
will be used.

## Examples

``` r
multigrain_verbosity()
#> [1] "info"
```
