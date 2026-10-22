# Calculate non-centrality parameter

Calculate the non-centrality parameter of the test statistic that
corresponds to the statistical power for a given number of hypotheses in
a trial.

## Usage

``` r
calc_ncp(power, alpha = 0.025)
```

## Arguments

- power:

  A numeric vector indicating each hypothesis' statistical power
  (probability of rejecting the null hypothesis if it false).

- alpha:

  A single numeric value representing the overall one-sided significance
  level. Default is 0.025.

## Value

A numeric vector of non-centrality parameters corresponding to each
element in `power`.

## Details

`calc_ncp()` uses the quantiles of the standard normal distribution to
compute the non-centrality parameter. It is given by the difference
between the critical value for a standard normal distribution and the
quantile corresponding to the specified power.

## Examples

``` r
# Basic usage
calc_ncp(power = c(0.8, 0.9))
#> [1] 2.801585 3.241516

# Custom significance level
calc_ncp(power = c(0.8, 0.9), alpha = 0.01)
#> [1] 3.167969 3.607899
```
