# Simulate raw p-values

`simulate_pvalues()` simulates raw p-values from under the alternative
hypotheses and the assumption that the distribution of test statistics
is a multivariate normal distribution. The most important input
parameters are the nominal power of each hypothesis, the significance
level, and a correlation matrix.

## Usage

``` r
simulate_pvalues(
  power_nominal,
  ...,
  alpha = 0.025,
  corr_matrix = diag(length(power_nominal)),
  nsim = 1e+05
)
```

## Arguments

- power_nominal:

  A numeric vector of nominal power values for each hypothesis.

- ...:

  These dots are for future extensions and must be empty.

- alpha:

  A single numeric value representing the overall one-sided significance
  level. Default is 0.025.

- corr_matrix:

  A numeric matrix representing the correlation matrix \\\Sigma\\ of the
  test statistics.

- nsim:

  An integer indicating the number of simulations to run. Default is
  `1e5`.

## Value

A matrix where each row represents a set of simulated raw p-values for
each hypothesis.

## Details

It starts by calculating the non-centrality parameter \\\Delta\\, of the
a test statistic z - where \\z \sim N(\Delta, 1)\\ - based on the
nominal power (probability of rejection under the alternative,
unadjusted for multiplicity) and significance level. It then generates
raw p-values by simulating the test statistics using a multivariate
normal distribution with the given correlation matrix.

`simulate_pvalues()` assumes a point global alternative, using a vector
of transformed p-values \\(\Phi^{-1}(1-p_1), \ldots, \Phi^{-1}(1-p_m))\\
which follows a multivariate normal distribution with a correlation
matrix \\\Sigma\\. Here, \\\Phi^{-1}\\ signifies the inverse function of
the standard normal distribution.

This assumption holds, for example, when \\p_1, \ldots, p_m\\ are the
raw p-values derived from one-sided z-tests for distinct hypotheses.

Note: applying the transformation \\\Phi^{-1}(1-p_i)\\ to p-values from
two-sided tests does not generally result in a multivariate normal
distribution.

## Examples

``` r
# Define parameters for simulation
nominal_power <- c(0.8, 0.85, 0.9)
corr <- matrix(
  c(1, 0.5, 0.5,
    0.5, 1, 0.5,
    0.5, 0.5, 1),
   nrow = 3
 )
num_simulations <- 1000

# Simulate raw p-values
pvals <- simulate_pvalues(
    nominal_power,
    corr_matrix = corr,
    nsim = num_simulations
)
```
