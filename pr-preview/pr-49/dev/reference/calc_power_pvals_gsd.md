# Calculate power for a group sequential graph-based multiple test procedure

Group sequential counterpart of
[`calc_power_pvals()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/calc_power_pvals.md).
Applies a graph to an object of repeated p-values built by
[`transform_pvalues_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/transform_pvalues_gsd.md)
and reports power by hypothesis and by analysis, together with the
distribution of the analysis at which each hypothesis was declared
rejected.

## Usage

``` r
calc_power_pvals_gsd(
  pvals,
  hyp_weight,
  trans_matrix,
  ...,
  alpha = NULL,
  custom_power = NULL,
  stopping_rule = NULL,
  sum_to_one_constraint = TRUE,
  call = rlang::caller_env()
)
```

## Arguments

- pvals:

  A `multigrain_pvals_gsd` object of repeated (or, per hypothesis,
  sequential) p-values, created with
  [`transform_pvalues_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/transform_pvalues_gsd.md).

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
  level to test at, or `NULL` (the default) to use the level of `pvals`.
  A supplied value must be greater than `0` and at most the level the
  repeated p-values were built up to.

- custom_power:

  A list of user-defined power measures. Alternatively a single measure
  may be given on its own. Three kinds are accepted, and each receives a
  different matrix:

  - a `multigrain_trial_success_gsd` object from
    [`trial_success_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/trial_success_gsd.md)
    is evaluated by its compiled function on the integer matrix of
    decision times (simulations by hypotheses, `0` for never rejected);

  - a fixed-sample `multigrain_trial_success` object from
    [`trial_success()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/trial_success.md)
    is evaluated on the logical matrix of rejections, so a
    rejection-only gain reports the same number here as
    [`calc_power_pvals()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/calc_power_pvals.md)
    would;

  - a plain function is called once per simulated trial with that
    trial's integer vector of decision times, and the results are
    averaged.

  If the list has no names, the measures are reported as `"func1"`,
  `"func2"`, and so on.

- stopping_rule:

  A `multigrain_stopping_rule_gsd` object from
  [`stopping_rule_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/stopping_rule_gsd.md),
  or `NULL` (the default) for a trial that runs to its final analysis
  unless it rejects every hypothesis earlier.

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

- `local_power`: the proportion of simulations in which each hypothesis
  is rejected at any analysis.

- `local_power_by_analysis`: an `m` by `K` matrix whose entry \\(i, k)\\
  is the proportion of simulations in which hypothesis \\i\\ was
  rejected at analysis \\k\\ or earlier. Its last column is
  `local_power`.

- `exp_rejections`: the expected number of rejections.

- `disj_power`: the probability of rejecting at least one hypothesis.

- `conj_power`: the probability of rejecting all hypotheses.

- `mean_decision_look`: for each hypothesis, the mean analysis index
  over the simulations in which it was rejected, and `NA` when it never
  was.

- `time_distribution`: an `m` by `K + 1` matrix of proportions, with a
  `"never"` column followed by one column per analysis. Its rows sum to
  one, and one minus its `"never"` column equals `local_power` up to
  floating-point rounding.

- `stop_distribution`: a named numeric vector of length `K`, the
  proportion of simulations whose stopping look is each analysis; sums
  to one; without a rule the mass sits on analysis `K` except for trials
  that rejected every hypothesis earlier.

- `mean_stop_look`: the mean stopping look over the simulations.

- `"..."`: one entry per measure in `custom_power`.

Every field is computed from the stopped decision times, so all of them
honour the stopping rule: a hypothesis not rejected by the time a trial
stops counts as never rejected in that trial.

## Details

`calc_power_pvals_gsd()` runs the group sequential kernel once and
summarises both of its outputs: the rejection indicators and the
decision times (the analysis at which the procedure declared each
hypothesis rejected, `0` for never). Every quantity is on the
analysis-index scale; calendar times are not stored by the package.

## See also

[`calc_power_pvals()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/calc_power_pvals.md)
for the fixed-sample version and
[`graph_optimise_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/graph_optimise_gsd.md)
for the optimiser.

## Examples

``` r
# two hypotheses, an interim at 60% information and a final analysis
corr <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)

set.seed(1)
raw <- simulate_pvalues_gsd(
    power_nominal = c(0.9, 0.8),
    corr_matrix = corr,
    info_frac = c(0.6, 1),
    nsim = 500
)

pvals <- transform_pvalues_gsd(
    raw,
    spending = gsDesign::sfLDOF,
    grid_size = 128L
)

# a rejection at the final analysis is worth 75% of one at the interim
gain <- trial_success_gsd(
    d(t1) + d(t2),
    d = c(1, 0.75),
    verbose = "silent"
)

# stop the trial at the interim once H1 is rejected
result <- calc_power_pvals_gsd(
    pvals,
    hyp_weight = c(0.5, 0.5),
    trans_matrix = matrix(c(0, 1, 1, 0), nrow = 2),
    custom_power = list(gain = gain),
    stopping_rule = stopping_rule_gsd(r1)
)

result$local_power_by_analysis
#>    analysis 1 analysis 2
#> H1      0.338      0.852
#> H2      0.236      0.650
result$time_distribution
#>    never analysis 1 analysis 2
#> H1 0.148      0.338      0.514
#> H2 0.350      0.236      0.414
result$stop_distribution
#> analysis 1 analysis 2 
#>      0.338      0.662 
result$gain
#> [1] 1.27
```
