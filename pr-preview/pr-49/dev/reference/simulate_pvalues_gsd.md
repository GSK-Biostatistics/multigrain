# Simulate raw p-values for a group sequential design

`simulate_pvalues_gsd()` simulates raw p-values at several analyses
("looks") of a group sequential trial, under the canonical joint normal
model for the test statistics. It is the group sequential counterpart of
[`simulate_pvalues()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/simulate_pvalues.md):
`power_nominal` and `corr_matrix` describe the *final-look* statistics,
and the information fractions say how much of that information is
available at each analysis.

## Usage

``` r
simulate_pvalues_gsd(
  power_nominal,
  ...,
  alpha = 0.025,
  corr_matrix = diag(length(power_nominal)),
  info_frac,
  nsim = 1e+05
)
```

## Arguments

- power_nominal:

  A numeric vector of nominal power values for each hypothesis, at its
  final analysis.

- ...:

  These dots are for future extensions and must be empty.

- alpha:

  A number giving the overall one-sided significance level used to
  convert `power_nominal` into non-centrality parameters. Default is
  `0.025`.

- corr_matrix:

  A numeric matrix representing the correlation matrix \\\Sigma\\ of the
  test statistics at the final analysis.

- info_frac:

  Information fractions. Either a numeric vector of length `K`, applied
  to every hypothesis, or an `m` by `K` numeric matrix with `NA` at
  analyses where a hypothesis has no data. They must be positive and
  non-decreasing over the analyses where a hypothesis has data, and
  every hypothesis must have data at some analysis.

- nsim:

  An integer indicating the number of simulations to run. Default is
  `1e5`.

## Value

A numeric array of raw p-values with dimensions `c(nsim, m, K)`:
simulations by hypotheses by analyses. The `m` by `K` matrix of
information fractions is attached as an `info_frac` attribute, where
[`transform_pvalues_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/transform_pvalues_gsd.md)
picks it up.

## Details

The non-centrality parameter \\\Delta_i\\ of each hypothesis is obtained
from its nominal power with
[`calc_ncp()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/calc_ncp.md),
exactly as in
[`simulate_pvalues()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/simulate_pvalues.md).
Writing \\t\_{i,k}\\ for the information fraction of hypothesis \\i\\ at
analysis \\k\\, the test statistics follow the canonical joint
distribution

\$\$E\[Z\_{i,k}\] = \Delta_i \sqrt{t\_{i,k}}, \qquad
\mathrm{Corr}(Z\_{i,k}, Z\_{j,l}) = \rho\_{ij} \sqrt{\min(t\_{i,k},
t\_{j,l}) / \max(t\_{i,k}, t\_{j,l})},\$\$

where \\\rho\_{ij}\\ is the entry of `corr_matrix` for hypotheses \\i\\
and \\j\\. Two analyses of the same hypothesis at the same information
fraction are therefore perfectly correlated, and two different
hypotheses observed at the same information fraction have correlation
\\\rho\_{ij}\\. P-values are \\1 - \Phi(Z\_{i,k})\\.

Only the analyses that carry distinct information are simulated. Once a
hypothesis reaches full information its statistic no longer changes, so
the p-value at the analysis where it matured is copied forward to every
later analysis; drawing those columns instead would make the covariance
singular. An analysis where a hypothesis has no data (`NA` information
fraction) before it matures is `NA` in the returned array;
[`transform_pvalues_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/transform_pvalues_gsd.md)
turns that into a repeated p-value of 1, which is "cannot reject".

## Examples

``` r
# two hypotheses analysed at 60% and 100% of the planned information
corr <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)

set.seed(1)
pvals <- simulate_pvalues_gsd(
    c(0.8, 0.9),
    corr_matrix = corr,
    info_frac = c(0.6, 1),
    nsim = 100
)

dim(pvals)
#> [1] 100   2   2
attr(pvals, "info_frac")
#>      [,1] [,2]
#> [1,]  0.6    1
#> [2,]  0.6    1

# the second hypothesis has no data at the interim analysis, and the first
# matures there: its p-value is carried forward to the final analysis
set.seed(1)
pvals_na <- simulate_pvalues_gsd(
    c(0.8, 0.9),
    corr_matrix = corr,
    info_frac = rbind(c(1, 1), c(NA, 1)),
    nsim = 100
)

identical(pvals_na[, 1L, 1L], pvals_na[, 1L, 2L])
#> [1] TRUE
all(is.na(pvals_na[, 2L, 1L]))
#> [1] TRUE
```
