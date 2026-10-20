# Create an efficacy stopping rule for a group sequential design

Create a rule that stops a simulated group sequential trial for efficacy
at an interim analysis, given the hypotheses rejected so far. The rule
is one Boolean condition per interim analysis: when the condition for
interim analysis `k` is true, the trial stops after analysis `k` and no
hypothesis can be rejected at a later analysis. The final analysis
always ends the trial, so it has no condition.

## Usage

``` r
stopping_rule_gsd(condition)
```

## Arguments

- condition:

  The stopping condition, or one condition per interim analysis. To
  inject values from your R session, use rlang's injection operator
  [`!!`](https://rlang.r-lib.org/reference/injection-operator.html). One
  of:

  - an expression, such as `r1 && (r2 || r3)`, or a single string, such
    as `"r1 and (r2 or r3)"`: one condition, applied at every interim
    analysis;

  - a call to [`list()`](https://rdrr.io/r/base/list.html) with one
    expression or string per interim analysis, such as
    `list(r1, r1 && (r2 || r3))`, or a list of them injected with `!!`;

  - a character vector with one string per interim analysis, such as
    `c("r1", "r1 and (r2 or r3)")`, or one injected with `!!`.

## Value

A `multigrain_stopping_rule_gsd` object made up of:

- `conditions`: a list of the conditions, as language objects.

- `labels`: the conditions deparsed to strings, for display.

- `m_min`: the largest hypothesis index referred to by an `r<i>` symbol,
  the fewest hypotheses the rule can be used with.

- `n_interim`: `NULL` for a single condition applied at every interim
  analysis, otherwise the number of conditions supplied.

## Details

A condition is written using:

- `r1, r2, ..., rm`: the indicator that hypothesis \\H_i\\ has been
  rejected at the current interim analysis or at an earlier one;

- the logical operators `&&`, `||` and `!` (or the words `and`, `or` and
  `not` in a string);

- comparisons (`== != < <= > >=`) and arithmetic (`+ - * /`), so that a
  count of rejections such as `r1 + r2 + r3 >= 2` can be tested;

- numeric or logical scalars, and parentheses.

Every condition must refer to at least one `r<i>`. Decision times (`t1`,
...) and discount tables, which a group sequential gain can use (see
[`trial_success_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/trial_success_gsd.md)),
are not available to a stopping rule.

### When the rule is applied

The rule does not store the number of hypotheses `m` or of analyses `K`;
both are settled when it meets a `multigrain_pvals_gsd` object (see
[`transform_pvalues_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/transform_pvalues_gsd.md)).
The rule may then refer to no hypothesis beyond `m`. A single condition
is recycled to each of the `K - 1` interim analyses; otherwise one
condition per interim analysis is needed. A list of `K` conditions is
accepted with a warning: the final analysis always ends the trial, so
its condition is dropped. With `K = 1` there is no interim analysis, and
a single condition is accepted and has no effect. The rule is evaluated
once on every possible set of rejected hypotheses, `2^m` sets per
interim analysis, so it can be used with at most 16 hypotheses.

Stopping only removes rejections, so error control is unaffected: a
stopped trial keeps the rejections it has already declared and declares
none at a later analysis.

A trial that has rejected every hypothesis stops at that analysis
whether or not a rule is supplied, as nothing is left to test.

### Operator precedence

Operator precedence follows R's own rules. Parenthesise to make intent
explicit: write `r1 && (r2 || r3)` rather than relying on `&&` binding
tighter than `||`.

## See also

[`trial_success_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/trial_success_gsd.md)
for the group sequential gain,
[`transform_pvalues_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/transform_pvalues_gsd.md)
for the p-values the rule is settled against.

## Examples

``` r
# Stop at any interim analysis once H1 and at least one of H2 and H3 are
# rejected
stopping_rule_gsd(r1 && (r2 || r3))
#> <multigrain_stopping_rule_gsd>
#> Every interim: r1 && (r2 || r3)
#> Hypotheses referenced: up to r3

# One condition per interim analysis (K = 3): stop at the first interim
# once H1 is rejected, at the second only once H2 or H3 is rejected too
stopping_rule_gsd(list(r1, r1 && (r2 || r3)))
#> <multigrain_stopping_rule_gsd>
#> Interim 1: r1
#> Interim 2: r1 && (r2 || r3)
#> Hypotheses referenced: up to r3

# The first rule, written as a string
stopping_rule_gsd("r1 and (r2 or r3)")
#> <multigrain_stopping_rule_gsd>
#> Every interim: r1 && (r2 || r3)
#> Hypotheses referenced: up to r3

# Stop once at least two of three hypotheses are rejected
stopping_rule_gsd(r1 + r2 + r3 >= 2)
#> <multigrain_stopping_rule_gsd>
#> Every interim: r1 + r2 + r3 >= 2
#> Hypotheses referenced: up to r3
```
