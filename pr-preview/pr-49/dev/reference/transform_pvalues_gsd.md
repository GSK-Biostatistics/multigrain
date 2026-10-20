# Transform group sequential p-values into repeated p-values

`transform_pvalues_gsd()` converts an array of nominal (raw) p-values
from a group sequential design into *repeated* p-values: for each
hypothesis and each analysis, the smallest level at which that p-value
would cross its alpha-spending boundary. Because the boundary is
monotone in the allocated level, the group sequential rejection rule
"\\p\_{i,k}\\ below the boundary for level \\w_i \alpha\\" becomes the
fixed-sample rule "\\p^r\_{i,k} \le w_i \alpha\\", so a graphical
procedure can be run on the transformed values without ever evaluating a
boundary.

## Usage

``` r
transform_pvalues_gsd(
  pvals,
  ...,
  info_frac = NULL,
  spending,
  alpha = 0.025,
  look_back = FALSE,
  grid_size = 1024L
)
```

## Arguments

- pvals:

  A numeric array of raw p-values with dimensions `c(nsim, m, K)`:
  simulations by hypotheses by analyses.

- ...:

  These dots are for future extensions and must be empty.

- info_frac:

  Information fractions. Either a numeric vector of length `K`, applied
  to every hypothesis, or an `m` by `K` numeric matrix with `NA` at
  analyses where a hypothesis has no data. Defaults to the `info_frac`
  attribute of `pvals`.

- spending:

  An alpha-spending function, or a list of `m` of them (one per
  hypothesis). A spending function takes a level and a vector of
  information fractions and returns the cumulative level spent, either
  as a plain numeric vector or as an object with a `spend` element (as
  the `gsDesign` spending functions do). It must spend the whole level
  by full information. Any entry may instead be
  [`sf_hsd_opt()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/sf_hsd_opt.md),
  which tags a hypothesis whose Hwang-Shih-DeCani parameter is to be
  optimised by
  [`graph_optimise_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/graph_optimise_gsd.md);
  a single
  [`sf_hsd_opt()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/sf_hsd_opt.md)
  passed on its own is recycled and is one parameter shared by every
  hypothesis. With any tagged entry the transform returns a spending
  spec instead (see Value).

- alpha:

  A number giving the overall one-sided significance level. Default is
  `0.025`.

- look_back:

  A logical of length 1 or `m`. `TRUE` switches a hypothesis to
  sequential p-values (see Details). Default is `FALSE`.

- grid_size:

  A whole number giving the number of allocated levels in each boundary
  table. Default is `1024`.

## Value

A `multigrain_pvals_gsd` object made up of:

- `pvals`: the transformed array, of the same dimensions as the input.

- `nsim`, `m`, `K`: simulations, hypotheses and analyses.

- `alpha`: the level the boundary tables were built up to.

- `info_frac`: the `m` by `K` matrix of information fractions.

- `look_back`: the per-hypothesis semantics actually used.

- `spending`: a label per spending function.

- `tables`: per hypothesis, the analyses that carry distinct
  information, the analysis at which it matured, and the grid of
  allocated levels with the nominal boundaries computed at each of them.

When any entry of `spending` is
[`sf_hsd_opt()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/sf_hsd_opt.md),
a `multigrain_pvals_gsd_spec` object instead. It holds the raw p-values
(as `raw`, not `pvals`), the hypotheses with fixed spending already
transformed, the lattice of each spending parameter and one boundary
table per lattice node, each checked for well ordering when the spec is
built. It is not a `multigrain_pvals_gsd` object and cannot be passed to
[`calc_power_pvals_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/calc_power_pvals_gsd.md):
pass it to
[`graph_optimise_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/graph_optimise_gsd.md)
to optimise the parameters with the graph, or to
[`resolve_pvalues_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/resolve_pvalues_gsd.md)
to obtain the transformed p-values at chosen values. If every tagged
hypothesis has a single analysis at full information, where spending has
no effect, an ordinary `multigrain_pvals_gsd` object is returned with a
message.

## Details

Boundaries are computed with
[`gsDesign::gsBound1()`](https://keaven.github.io/gsDesign//reference/gsBound.html)
on a log-spaced grid of allocated levels running from `1e-14` up to
`alpha`, and the resulting monotone table is inverted by log-log linear
interpolation.

Per hypothesis, writing \\D\\ for the analyses with a non-missing
information fraction and \\k^{mat}\\ for the first of those that reaches
full information (or the last analysis in \\D\\ if none does):

- an analysis without data (`NA` information fraction) before
  \\k^{mat}\\ is given a repeated p-value of 1, i.e. "cannot reject";

- analyses after \\k^{mat}\\ copy the repeated p-value of \\k^{mat}\\
  forward, because the statistic of a matured endpoint does not change
  (it may still be rejected later if recycling gives it more alpha). A
  non-missing raw p-value supplied after maturity is ignored, with a
  warning if it differs from the maturity value;

- a hypothesis with a single analysis at full information is returned
  unchanged, since there the boundary *is* the level;

- if `look_back` is `TRUE` for a hypothesis its repeated p-values are
  replaced by their running minimum over analyses (the *sequential*
  p-value), which allows a rejection on the strength of the evidence at
  any earlier analysis.

If a boundary table is not non-decreasing in the allocated level the
spending function is not "well ordered", the group sequential graphical
procedure is not valid for it, and the transform aborts.

## Examples

``` r
# two hypotheses analysed at 60% and 100% of the planned information
set.seed(1)
raw <- array(stats::runif(100 * 2 * 2), dim = c(100L, 2L, 2L))

pvals_gsd <- transform_pvalues_gsd(
    raw,
    info_frac = c(0.6, 1),
    spending = gsDesign::sfLDOF,
    grid_size = 128L
)

pvals_gsd
#> <multigrain_pvals_gsd>
#> 100 simulation(s), 2 hypothes(es), 2 analys(es); alpha = 0.025
#> Information fractions (NA = no data at that analysis):
#>    analysis 1 analysis 2
#> H1        0.6          1
#> H2        0.6          1
#> Look-back: H1 = FALSE, H2 = FALSE
```
