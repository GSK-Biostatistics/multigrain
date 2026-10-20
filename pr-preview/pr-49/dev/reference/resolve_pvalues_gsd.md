# Resolve a spending spec at chosen parameter values

When any spending function passed to
[`transform_pvalues_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/transform_pvalues_gsd.md)
is tagged with
[`sf_hsd_opt()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/sf_hsd_opt.md),
the transform returns a `multigrain_pvals_gsd_spec`: the raw p-values,
the hypotheses with fixed spending already transformed, and the boundary
table of every lattice node of every tagged parameter.
`resolve_pvalues_gsd()` picks one node per parameter and returns the
`multigrain_pvals_gsd` object that
[`transform_pvalues_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/transform_pvalues_gsd.md)
would have returned had the tagged hypotheses used
`function(a, t) gsDesign::sfHSD(a, t, value)$spend`, using the stored
tables, so it does not rebuild any boundary.

## Usage

``` r
resolve_pvalues_gsd(spec, spending_param)

# S3 method for class 'multigrain_pvals_gsd_spec'
print(x, ...)
```

## Arguments

- spec:

  A `multigrain_pvals_gsd_spec` object from
  [`transform_pvalues_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/transform_pvalues_gsd.md).

- spending_param:

  A numeric vector with one value per spending parameter of `spec`. If
  named, the names must be the parameter names (in any order). Every
  value must be a node of its parameter's lattice.

- x:

  A `multigrain_pvals_gsd_spec` object.

- ...:

  Not used.

## Value

A `multigrain_pvals_gsd` object; the `spending` label of each tagged
hypothesis reads, for example, `sfHSD(param = -2.5)`.
[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly.

## Details

Only lattice nodes can be resolved, since only their tables were built
and checked; a value off the lattice is refused, naming the nearest
nodes. Resolving is quick (an interpolation per hypothesis and
analysis), so it is the way to evaluate a fixed choice of parameters, to
reproduce the final evaluation of
[`graph_optimise_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/graph_optimise_gsd.md)
with `resolve_pvalues_gsd(spec, result$spending_param)`, or to check a
chosen design on fresh simulated p-values: build a new spec from them
with the same spending functions and resolve it at the chosen values.

[`print()`](https://rdrr.io/r/base/print.html) of a spec shows its
dimensions, the spending of each hypothesis, each parameter's bounds,
step, number of nodes and start, and the number of boundary tables
built.

## See also

[`sf_hsd_opt()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/sf_hsd_opt.md),
[`transform_pvalues_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/transform_pvalues_gsd.md),
[`graph_optimise_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/graph_optimise_gsd.md),
[`calc_power_pvals_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/calc_power_pvals_gsd.md)

## Examples

``` r
spending <- list(
    gsDesign::sfLDOF,
    sf_hsd_opt(id = "secondary", bounds = c(-4, -2), step = 0.5)
)

set.seed(1)
raw <- array(stats::runif(100 * 2 * 2), dim = c(100L, 2L, 2L))
spec <- transform_pvalues_gsd(
    raw,
    info_frac = c(0.5, 1),
    spending = spending,
    grid_size = 64L
)

spec
#> <multigrain_pvals_gsd_spec>
#> 100 simulation(s), 2 hypothes(es), 2 analys(es); alpha = 0.025
#> Information fractions (NA = no data at that analysis):
#>    analysis 1 analysis 2
#> H1        0.5          1
#> H2        0.5          1
#> Look-back: H1 = FALSE, H2 = FALSE
#> Spending:
#> H1: spending[[1]]
#> H2: sfHSD(param = secondary in [-4, -2])
#> Spending parameters:
#> secondary: bounds [-4, -2], step 0.5, 5 nodes, start -3
#> Boundary tables built: 5

resolve_pvalues_gsd(spec, c(secondary = -3.5))
#> <multigrain_pvals_gsd>
#> 100 simulation(s), 2 hypothes(es), 2 analys(es); alpha = 0.025
#> Information fractions (NA = no data at that analysis):
#>    analysis 1 analysis 2
#> H1        0.5          1
#> H2        0.5          1
#> Look-back: H1 = FALSE, H2 = FALSE
```
