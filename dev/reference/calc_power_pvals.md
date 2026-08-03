# Calculate power for a graph-based multiple test procedure using p-values

Calculate multiplicity-adjusted (marginal) power values for a
graph-based multiple test procedure applied to a matrix of p-values.
Provide both the `hyp_weight` vector and transition matrix
`trans_matrix` for the graph you wish to evaluate.

## Usage

``` r
calc_power_pvals(
  pvals,
  hyp_weight,
  trans_matrix,
  ...,
  alpha = 0.025,
  custom_power = NULL,
  sum_to_one_constraint = TRUE,
  call = rlang::caller_env()
)
```

## Arguments

- pvals:

  A numeric matrix of p-values. Each row represents a simulated trial,
  and each column corresponds to a hypothesis.

- hyp_weight:

  A numeric vector representing the weights of the nodes (hypotheses) in
  the graph. All elements must be in the interval \[0, 1\] and the
  vector must sum to 1.

- trans_matrix:

  A numeric matrix representing the transition matrix between nodes. All
  elements must be in the interval \[0, 1\], with diagonal elements
  equal to 0, and each row must sum to 1.

- ...:

  These dots are for future extensions and must be empty.

- alpha:

  A single numeric value representing the overall one-sided significance
  level. Default is 0.025.

- custom_power:

  A list of user-defined power functions. Alternatively, a single
  function or `multigrain_trial_success` object can be provided. There
  are two ways to specify this:

  - Anonymous functions can be provided to specify the success criteria.
    Functions must take one simulation's logical vector of results as an
    input, and return a scalar. For example, if one is interested in the
    power to reject hypotheses 1 and 3 one could specify:
    `f = function(x) {x[1] && x[3]}`. If the power of rejecting
    hypotheses 1 and 2 is also of interest one would use an (optionally
    named) list:
    `f = list( power1and3 = function(x) {x[1] && x[3]}, power1and2 = function(x) {x[1] && x[2]} )`.
    If the list has no names, the functions will be referenced as
    `"func1"`, `"func2"`, etc. in the output. The user can also provide
    a `multigrain_trial_success` object instead (resulting in a faster
    calculation).

  - Instead of anonymous functions, one can pass a
    `multigrain_trial_success` object (or within a list as with the
    anonymous functions). See
    [`trial_success()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/trial_success.md)
    and the examples.

- sum_to_one_constraint:

  A logical value controlling whether to allow graphs where transition
  matrix rows are not constrained to sum to one, for example in a fixed
  sequence. Defaults to `TRUE`.

- call:

  The execution environment of a currently running function, e.g.
  `caller_env()`. The function will be mentioned in error messages as
  the source of the error. See the `call` argument of
  [`abort()`](https://rlang.r-lib.org/reference/abort.html) for more
  information.

## Value

A list containing:

- `local_power`: The local power for each hypothesis: the proportion of
  simulations in which each hypothesis is rejected.

- `exp_rejections`: The expected number of rejections across all
  simulations.

- `disj_power`: Disjunctive power: the probability of rejecting at least
  one hypothesis.

- `conj_power`: Conjunctive power: the probability of rejecting all
  hypotheses.

- `"..."`: Additional results as specified by user-defined success
  functions in `custom_power`.

## Details

`calc_power_pvals()` calculates several power metrics:

- *local power* (probability to reject each individual hypothesis),

- *disjunctive power* (the probability to reject at least one
  hypothesis),

- *conjunctive power* (the probability to reject all hypotheses), and

- the expected number of rejections. Optionally, you can specify
  user-defined trial success criteria via `custom_power`.

## Examples

``` r

# First we simulate our p-value distribution

power_nominal <- c(0.90, 0.87, 0.73)
corr_matrix <- matrix(
  c(
    1, 0.2, 0.2,
    0.2, 1, 0.2,
    0.2, 0.2, 1
  ),
  nrow = 3,
  byrow = TRUE
)

pvals <- simulate_pvalues(
  power_nominal = power_nominal,
  corr_matrix = corr_matrix
)

# Second we construct our graph (in this case a fixed sequence)

hyp_weights <- c(1, 0, 0)
trans_matrix <- matrix(
  c(
    0, 1, 0,
    0, 0, 1,
    0, 0, 0
  ),
  nrow = 3,
  byrow = TRUE
)

# Third, we construct a list of metrics we wish to evaluate our graph with

trial_success_measure <- trial_success((r1 && r2) || r3)
#> ✔ Trial success function compiled and sourced successfully.
power_metrics <- list(
  trial_success = trial_success_measure,
  average_power = function(x) {
    x[1] + x[2] + x[3]
  }
)

# Finally we calculate the power of graph conditional on p-value distribution

result <- calc_power_pvals(
  pvals = pvals,
  hyp_weight = hyp_weights,
  trans_matrix = trans_matrix,
  custom_power = power_metrics,
  sum_to_one_constraint = FALSE # As third row of graph does not sum to 1
)

result
#> $local_power
#> [1] 0.90018 0.79120 0.60123
#> 
#> $exp_rejections
#> [1] 2.29261
#> 
#> $disj_power
#> [1] 0.90018
#> 
#> $conj_power
#> [1] 0.60123
#> 
#> $trial_success
#> [1] 0.7912
#> 
#> $average_power
#> [1] 2.29261
#> 
```
