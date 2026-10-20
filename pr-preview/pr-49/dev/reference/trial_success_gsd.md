# Create a trial success function for a group sequential design

Create a user-defined **trial-success utility** \\\psi\\ for a group
sequential design. It assigns value to each simulated trial from the
hypotheses that were rejected *and the analysis at which each rejection
was declared*, so that an early rejection can be worth more than a late
one. The function compiles \\\psi\\ to C++ for simulation and
optimisation, exactly as
[`trial_success()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/trial_success.md)
does for a fixed-sample design.

## Usage

``` r
trial_success_gsd(objective, ..., K = NULL, verbose = multigrain_verbosity())
```

## Arguments

- objective:

  An expression or string encoding the trial-success utility \\\psi\\
  over `r1, r2, ...` and `t1, t2, ...`. To inject values from your R
  session, use rlang's injection operator -
  [`!!`](https://rlang.r-lib.org/reference/injection-operator.html) -
  (see Examples). Use `!!` only inside `objective`; it has a different
  meaning in ordinary arguments (see **Common mistakes**).

- ...:

  Discount tables, each a named numeric vector of length `K` (one value
  per analysis). Supply only the values for analyses 1 through `K`; the
  value for "never rejected" (`ti = 0`) is always `0` and is supplied
  automatically — do not include it yourself (see **Common mistakes**).
  Every table must be named; the names are used as functions in
  `objective`.

- K:

  An optional whole number giving the number of analyses. When discount
  tables are supplied `K` is inferred from their common length and, if
  also given, must agree with it. When neither is given `K` is left
  `NULL` and checked against the p-values at first use.

- verbose:

  An optional string controlling verbosity (`"detail"` \> `"info"` \>
  `"silent"`). Verbosity can also be set at package level with the
  `multigrain_verbosity` option (see
  [`multigrain_verbosity()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/multigrain_verbosity.md)):

  - `"info"` (default): informs about successful compilation of the
    trial success function.

  - `"detail"`: same as `"info"` (no additional detail available).

  - `"silent"`: no compilation messages. Errors and warnings are still
    thrown.

## Value

A `multigrain_trial_success_gsd` object (which also inherits from
`multigrain_trial_success`) made up of:

- `func`: compiled function that evaluates \\\psi\\ row-wise on an
  integer matrix of decision times (simulations by hypotheses, `0` for
  never rejected) and returns the mean utility.

- `m`: number of hypotheses implied by the largest index among the
  `r<i>` and `t<i>` symbols.

- `K`: number of analyses, or `NULL` if unknown.

- `objective`: the original utility expression, as a display string
  (numbers shown to 15 significant digits; the compiled function uses
  the exact values).

- `cpp_code`: the generated C++ source.

- `tables`: the discount tables, as supplied.

## Details

In code, \\\psi\\ is written using:

- `r1, r2, ..., rm`: the binary indicator that hypothesis \\H_i\\ is
  rejected at any analysis;

- `t1, t2, ..., tm`: the analysis ("look") at which \\H_i\\ was
  rejected, an integer from `1` to `K`, or `0` if it was never rejected.
  `ri` is the same as `ti > 0`;

- arithmetic (`+ - * /`), logical operators (`&&`, `||`, or the words
  `and`, `or` in a string) between logical terms, comparisons
  (`== != < <= > >=`), and parentheses;

- discount tables: any numeric vector of length `K` passed by name
  through `...` can be applied to a decision time as `name(ti)`. The
  result is the `ti`-th element of the table, or `0` when `ti` is `0`
  (never rejected).

  ### Decision times

  Each `ti` records the analysis at which the graphical procedure
  *declared* hypothesis \\H_i\\ rejected. This may be later than the
  analysis at which the data for \\H_i\\ matured or first crossed a
  boundary. For example, suppose PFS data are final at analysis 1 but
  PFS carries zero initial weight. PFS cannot reject until OS rejects at
  analysis 2 and recycles its alpha to PFS. The result is `t_PFS = 2`,
  even though the PFS p-value was available at analysis 1.

  The rejection indicator `ri` is derived from `ti > 0` inside the
  compiled function, so `ri` and `ti` can never disagree: if `ti = 0`
  then `ri = FALSE`, and if `ti >= 1` then `ri = TRUE`.

  ### Discount tables

  A discount table maps each decision time to a value. You supply one
  value per analysis, for analyses 1 through `K`. The constructor
  automatically prepends `0.0` for `ti = 0` (never rejected), so you
  should **not** include a leading zero yourself. For a two-analysis
  design where a rejection at the final analysis is worth 75% of one at
  the interim, write `d = c(1, 0.75)`. Element names on the vector are
  ignored; only the argument name matters (the name you use to apply it
  in `objective`, e.g. `d(t1)`).

  ### Operator precedence

  Operator precedence follows R's own rules (unlike
  [`trial_success()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/trial_success.md),
  whose expression-input parser treats `&&` and `||` at a different
  precedence level). Always parenthesise to make intent explicit: write
  `(t1 == 1 && t2 == 1)` rather than relying on `&&` binding tighter
  than `+`.

## Common mistakes

**Including a leading zero in a discount table.** Writing
`d = c(0, 1, 0.75)` does not give a two-analysis table — it gives a
*three*-analysis table with `K = 3`, where analysis 1 has value 0 and
analysis 2 has value 1, so every claim is mispriced. The increase from
`0` to `1` triggers a warning, but the table still compiles. For
`K = 2`, write `d = c(1, 0.75)`.

**Using `!!` inside a discount table argument.** The injection operator
[`!!`](https://rlang.r-lib.org/reference/injection-operator.html) is
only meaningful inside `objective`, where
[`rlang::enexpr()`](https://rlang.r-lib.org/reference/defusing-advanced.html)
captures the expression. In an ordinary argument like
`d = c(1, !!delta)`, `!!` is R's double negation: `!!0.75` evaluates to
`TRUE` (i.e. `1`), so the table silently becomes `c(1, 1)` and the
discount is lost. Write `d = c(1, delta)` instead.

**Passing a GSD object to the fixed-sample optimiser.** This is an
error.
[`graph_optimise()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/graph_optimisation.md)
and
[`calc_power_pvals()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/calc_power_pvals.md)
supply a logical rejection matrix, not a decision-time matrix, so they
refuse a `multigrain_trial_success_gsd` object and point you at the
group sequential consumers,
[`graph_optimise_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/graph_optimise_gsd.md)
and
[`calc_power_pvals_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/calc_power_pvals_gsd.md).

## Interpreting the gain

The gain returned by `$func()` is the mean of \\\psi\\ over simulated
trials and is **not** a probability. Its scale depends entirely on how
you define `objective`. For instance, \\\psi = r_1 + r_2\\ at 60% power
per hypothesis returns roughly 1.2, not 0.6. Only relative values matter
for optimisation (the argmax is invariant to positive affine
transformations of \\\psi\\), but normalise for reporting if you want
interpretable numbers.

A discount table is expected to take values in `[0, 1]` and to assign no
more value to a later decision than to an earlier one. A table that
breaks either expectation still compiles — the package does not forbid
other value tables — but it warns once, naming the table and the
offending values or analyses. Both the range and the ordering are judged
with a tolerance of `sqrt(.Machine$double.eps)`, so a table computed
rather than typed is not flagged for rounding: a value of `1 + 1e-12` is
not reported as outside `[0, 1]`. The first value need not be `1`.

## See also

[`trial_success()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/trial_success.md)
for the fixed-sample version,
[`graph_optimise_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/graph_optimise_gsd.md)
and
[`calc_power_pvals_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/calc_power_pvals_gsd.md)
for the consumers of the resulting object.

## Examples

``` r
# Manuscript Example 5: PFS (H1) and OS (H2), two analyses; a rejection at
# the final analysis is worth 75% of one at the interim.
# d = c(1, 0.75): full value at analysis 1, 75% at analysis 2.
# The zero for "never rejected" is supplied automatically.
v_pfs <- 0.4
v_os  <- 0.6
gain <- trial_success_gsd(
    !!v_pfs * d(t1) + !!v_os * d(t2),
    d = c(1, 0.75)
)
#> ✔ Trial success function compiled and sourced successfully.
gain
#> <multigrain_trial_success_gsd/multigrain_trial_success>
#> 0.4 * d(t1) + 0.6 * d(t2)
#> Analyses (K): 2
#> Discount table d(t): 1.00, 0.75

# The gain on every possible (t1, t2) pair for K = 2:
#
#   t1 \ t2 |    0      1      2
#   --------+---------------------
#      0    | 0.000  0.600  0.450
#      1    | 0.400  1.000  0.850
#      2    | 0.300  0.900  0.750

# \donttest{
# Mixed r and t: OS value plus a bonus only if PFS was declared at
# the interim
trial_success_gsd(r2 + 0.5 * (r1 && t1 == 1), K = 2)
#> ✔ Trial success function compiled and sourced successfully.
#> <multigrain_trial_success_gsd/multigrain_trial_success>
#> r2 + 0.5 * (r1 && t1 == 1)
#> Analyses (K): 2

# Hurdle gain with a time-discounted base value: both claims required,
# and the package is worth less if the confirmatory claim is late
trial_success_gsd((r1 && r2) * d(t2), d = c(1, 0.75))
#> ✔ Trial success function compiled and sourced successfully.
#> <multigrain_trial_success_gsd/multigrain_trial_success>
#> (r1 && r2) * d(t2)
#> Analyses (K): 2
#> Discount table d(t): 1.00, 0.75

# Co-primary endpoints rejected at the same analysis
trial_success_gsd((t1 == 1 && t2 == 1) + 0.5 * (t1 == 2 && t2 == 2))
#> ✔ Trial success function compiled and sourced successfully.
#> <multigrain_trial_success_gsd/multigrain_trial_success>
#> (t1 == 1 && t2 == 1) + 0.5 * (t1 == 2 && t2 == 2)
#> Analyses (K): not set

# Rejection indicators work as in trial_success()
trial_success_gsd(r1 + r2, K = 2)
#> ✔ Trial success function compiled and sourced successfully.
#> <multigrain_trial_success_gsd/multigrain_trial_success>
#> r1 + r2
#> Analyses (K): 2
# }
```
