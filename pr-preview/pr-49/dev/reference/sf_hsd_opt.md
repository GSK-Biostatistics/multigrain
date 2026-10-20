# Hwang-Shih-DeCani spending with an optimised parameter

`sf_hsd_opt()` tags a hypothesis whose Hwang-Shih-DeCani spending
parameter should be chosen by
[`graph_optimise_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/graph_optimise_gsd.md)
jointly with the graph, rather than fixed in advance. Pass it, alone or
inside the list of spending functions, to
[`transform_pvalues_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/transform_pvalues_gsd.md),
which then returns a `multigrain_pvals_gsd_spec` object instead of
transformed p-values.

## Usage

``` r
sf_hsd_opt(id = NULL, bounds = c(-4, 1), step = 0.25, start = NULL)
```

## Arguments

- id:

  `NULL`, or a single string or number naming the parameter. Hypotheses
  with the same `id` share one parameter, and the `id` is the name under
  which its value is reported.

- bounds:

  A numeric vector of length 2: the lowest and highest value of the
  parameter to search. The lower bound must be at least -10, the upper
  at most 2, and the lower below the upper. Default is `c(-4, 1)`.

- step:

  The spacing of the lattice of values searched: one of `0.5`, `0.25` or
  `0.125`. Default is `0.25`.

- start:

  `NULL` or a number inside `bounds`: the value the search starts from,
  snapped to the nearest lattice node. `NULL` starts from the node
  nearest the midpoint of `bounds`.

## Value

A spending function `function(alpha, t, param)` of class
`multigrain_sf_hsd_opt`, returning
`gsDesign::sfHSD(alpha, t, param)$spend`, with attributes `id`,
`bounds`, `step` and `start`.

## Details

The Hwang-Shih-DeCani parameter sets how quickly the level is spent: at
an interim analysis at half the information, a parameter of -4 has spent
about 12% of the level, a boundary like O'Brien and Fleming's, and a
parameter of 1 about 62%, close to Pocock-type spending; 0 spends in
proportion to information.

The parameter is searched over a lattice: the multiples of `step` inside
`bounds`, together with both bounds (21 nodes for the defaults). Every
node's boundary table is built, and checked for well ordering and
invertibility, when the spec is created, so every value the search can
try is a valid group sequential transform. Building the tables takes
about 0.9 s per node at the default `grid_size` of
[`transform_pvalues_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/transform_pvalues_gsd.md)
with two analyses; a progress bar is shown when it takes more than a
couple of seconds.

Tagged hypotheses with the same `id` share one parameter. Tagged
hypotheses without an `id` each get their own, except that a single
`sf_hsd_opt()` passed as the whole `spending` argument (and so recycled
across hypotheses) is one parameter shared by all of them. Hypotheses
sharing a parameter must use identical `bounds`, `step` and `start`.

A tagged hypothesis with a single analysis at full information is not
affected by its spending function and takes no part in the search; if
that leaves no parameter to optimise,
[`transform_pvalues_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/transform_pvalues_gsd.md)
returns an ordinary `multigrain_pvals_gsd` object with a message.

## See also

[`gsDesign::sfHSD()`](https://keaven.github.io/gsDesign//reference/sfHSD.html),
[`transform_pvalues_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/transform_pvalues_gsd.md),
[`graph_optimise_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/graph_optimise_gsd.md)
for the search and the result,
[`resolve_pvalues_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/resolve_pvalues_gsd.md)

## Examples

``` r
# one parameter shared by hypotheses 2 and 3, searched over -4, -3.5, ..., -2
spending <- list(
    gsDesign::sfLDOF,
    sf_hsd_opt(id = "secondary", bounds = c(-4, -2), step = 0.5),
    sf_hsd_opt(id = "secondary", bounds = c(-4, -2), step = 0.5)
)

set.seed(1)
raw <- array(stats::runif(100 * 3 * 2), dim = c(100L, 3L, 2L))
spec <- transform_pvalues_gsd(
    raw,
    info_frac = c(0.5, 1),
    spending = spending,
    grid_size = 64L
)
spec
#> <multigrain_pvals_gsd_spec>
#> 100 simulation(s), 3 hypothes(es), 2 analys(es); alpha = 0.025
#> Information fractions (NA = no data at that analysis):
#>    analysis 1 analysis 2
#> H1        0.5          1
#> H2        0.5          1
#> H3        0.5          1
#> Look-back: H1 = FALSE, H2 = FALSE, H3 = FALSE
#> Spending:
#> H1: spending[[1]]
#> H2: sfHSD(param = secondary in [-4, -2])
#> H3: sfHSD(param = secondary in [-4, -2])
#> Spending parameters:
#> secondary: bounds [-4, -2], step 0.5, 5 nodes, start -3
#> Boundary tables built: 5
```
