# Create a trial success function

Create a user-defined **trial-success utility** \\\psi\\ that assigns
value to each rejection pattern from a graphical multiple testing
procedure. The function compiles \\\psi\\ to fast C++ for
simulation/optimisation.

## Usage

``` r
trial_success(objective, verbose = multigrain_verbosity())
```

## Arguments

- objective:

  An expression or string encoding the trial-success utility \\\psi\\.
  The symbols `r1, r2, ...` refer to rejection indicators for the
  corresponding hypotheses. To inject values from your R session, use
  rlang's unquote operator `!!` (see Examples). Arithmetic and logical
  operators are allowed.

- verbose:

  An optional string controlling verbosity ("detail" \> "info" \>
  "silent"). Verbosity can also be set at package level with the
  `multigrain_verbosity` option (see
  [`multigrain_verbosity()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/multigrain_verbosity.md)):

  - `"info"` (default): will inform about the successful compilation of
    the trial success function.

  - `"detail"`: no additional information available, will have the same
    effect as `"info"`.

  - `"silent"`: silent, no information about the trial success function
    compilation. Errors and warnings are still thrown.

## Value

A `multigrain_trial_success` object made up of:

- `func`: compiled function that evaluates \\\psi\\ row-wise on a matrix
  of rejection indicators and returns the mean utility (i.e., expected
  trial success under the simulated scenario).

- `m`: number of hypotheses implied by `r1, ..., rm`.

- `objective`: the original utility expression (as a string).

## Details

In code, \\\psi\\ is written using symbols `r1, r2, ..., rm`, where each
`ri` is the binary indicator that hypothesis \\H_i\\ is rejected by the
chosen graph. You can combine these indicators with arithmetic
(`+ - * /`) and logical operators (`&&`, `||`, or the words `and`, `or`)
to reflect your design priorities (e.g., “any success”, “all successes”,
weighted composites).

The utility \\\psi\\ is evaluated on each simulated rejection pattern
and averaged, yielding the **expected trial success** for a given graph
and data-generating scenario. This lets you optimise graphs against the
utility that captures your clinical/regulatory goals, rather than a
single power summary.

## Examples

``` r

# Expected number of rejections (equals m × average power)
exp_rejs <- trial_success(r1 + r2 + r3 + r4)
#> ✔ Trial success function compiled and sourced successfully.

# \donttest{
# Disjunctive success: any rejection among four
disj_power <- trial_success(r1 || r2 || r3 || r4)
#> ✔ Trial success function compiled and sourced successfully.

# Conjunctive success: all four must be rejected
conj_power <- trial_success(r1 && r2 && r3 && r4)
#> ✔ Trial success function compiled and sourced successfully.

# Composite utility: require H1 AND H2 rejection to get H3 and H4 value
composite <- trial_success((r1 && r2) * (r3 + r4))
#> ✔ Trial success function compiled and sourced successfully.

# Weighted priorities (e.g., H1 gets weight 2 times H2 or H3)
weighted <- trial_success(2*r1 + r2 + r3)
#> ✔ Trial success function compiled and sourced successfully.

# Inject values from your environment with !!
w <- 2
trial_success(!!w * r1 + r2 + r3)
#> ✔ Trial success function compiled and sourced successfully.
#> <multigrain_trial_success>
#> 2 * r1 + r2 + r3

# Programmatic string input
expr_str <- sprintf("%s * r1 + r2", w)
trial_success(!!expr_str)
#> ✔ Trial success function compiled and sourced successfully.
#> <multigrain_trial_success>
#> 2 * r1 + r2

# }
```
