# Create a *trial success* function

Create a user-defined **trial-success utility** (gain function) \\\psi\\
that assigns a numerical value to each rejection pattern from a
graphical multiple testing procedure. The function compiles \\\psi\\ to
C++ for fast evaluation during simulation and optimisation.

## Usage

``` r
trial_success(objective, verbose = multigrain_verbosity())
```

## Arguments

- objective:

  An expression or string encoding the trial-success utility \\\psi\\.
  The symbols `r1, r2, ...` refer to rejection indicators for the
  corresponding hypotheses. The number of hypotheses \\m\\ is inferred
  as the largest index found (e.g., `r1 + r3` implies \\m = 3\\ and
  warns about the gap at `r2`). To inject values from your R session,
  use rlang's injection operator
  [`!!`](https://rlang.r-lib.org/reference/injection-operator.html) (see
  Examples).

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

A `multigrain_trial_success` object made up of:

- `func`: compiled function that evaluates \\\psi\\ row-wise on a
  logical matrix of rejection indicators (\\N \times m\\) and returns
  the mean utility (a single scalar).

- `m`: number of hypotheses implied by the `r1, ..., rm` symbols.

- `objective`: the original utility expression (as a string).

- `cpp_code`: the generated C++ source that was compiled.

## Details

In code, \\\psi\\ is written using symbols `r1, r2, ..., rm`, where each
`ri` is the binary indicator that hypothesis \\H_i\\ is rejected by the
chosen graph. You can combine these indicators with arithmetic
(`+ - * /`) and logical operators (`&&`, `||`, or the words `and`, `or`)
to reflect your design priorities (e.g., “any claim”, “all claims”,
weighted composites, gated secondaries).

### How it works

The expression is compiled to a C++ function that loops over the rows of
an \\N \times m\\ `LogicalMatrix` (one row per simulated trial, one
column per hypothesis). For each row, `r1` is `1` if the hypothesis was
rejected and `0` otherwise. The compiled function evaluates \\\psi\\ on
every row, sums the results, and divides by \\N\\: `$func()` returns one
scalar, never a vector.

### Injection with `!!`

Use rlang's
[`!!`](https://rlang.r-lib.org/reference/injection-operator.html) to
splice a value from your environment into the expression at build time.
The injected value becomes a compile-time literal in the C++ code:

    w <- 0.6
    trial_success(!!w * r1 + (1 - !!w) * r2)

Without `!!`, a bare symbol like `w` is treated as an unknown rejection
indicator and produces an informative error.

## Note

**Known limitations:**

- Unary minus (e.g., `trial_success(-r1 + r2)`) errors with
  `subscript out of bounds`. Use `trial_success(0 - r1 + r2)` or
  `trial_success(r2 - r1)` as a workaround.

- An index gap (e.g., `r1 + r3` with no `r2`) produces a **warning** and
  sets \\m = 3\\, not an error. The compiled function expects a matrix
  with 3 columns; column 2 is present but never referenced.

## Operators and precedence

**Expression input** accepts: `+`, `-`, `*`, `/`, `&&`, `||`, and
parentheses. The words `and` and `or` (case-insensitive) are accepted in
string input as aliases for `&&` and `||`.

**Precedence warning:** the internal parser rewrites `&&` and `||` to
placeholder operators before R parses the expression. These placeholders
bind tighter than `+`, `-`, `*`, and `/` — which differs from standard R
precedence. For example, `trial_success(r1 + r2 && r3)` compiles to
`r1 + (r2 && r3)`, not `(r1 + r2) && r3`. **Always use parentheses**
around logical sub-expressions to make the grouping explicit.

## String input

When `objective` is a character string, it is passed to the C++ code
generator **without** the operator whitelist check that expression input
receives. This means unsupported operators or functions (e.g., `^`,
`min`, `abs`) may compile but produce unexpected behaviour, or fail with
a raw g++ diagnostic rather than an R-level error message. Prefer
expression input unless you need to build the objective
programmatically.

## Interpreting the gain

The value returned by `$func()` is the **mean of \\\psi\\** across
simulated trials. It is **not** a probability and is not bounded by
\\\[0, 1\]\\. For example, `trial_success(r1 + r2 + r3 + r4)` with \\m =
4\\ hypotheses each at 60\\ 0.6. This is \\m\\ times average power.

Only **relative** values of the gain matter for optimisation: the graph
that maximises \\\psi\\ also maximises \\100\psi\\ and \\\psi + 5\\.
Normalise the scale when comparing or reporting results across different
gain functions.

## References

Spiers, A. D. V., Grayling, M. J., Wheeler, G. M., and Mander, A. P.
(2026). Gain-function optimisation of graphical multiple testing
procedures for confirmatory clinical trials. *arXiv:2609.19994v1*.
<https://arxiv.org/abs/2609.19994>

## See also

[`trial_success_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/trial_success_gsd.md)
for the group sequential design version, which adds decision-time
variables and discount tables.

## Examples

``` r
# Expected number of rejections (equals m * average power)
exp_rejs <- trial_success(r1 + r2 + r3 + r4)
#> ✔ Trial success function compiled and sourced successfully.

# \donttest{
# Disjunctive success: any rejection among four
disj_power <- trial_success(r1 || r2 || r3 || r4)
#> ✔ Trial success function compiled and sourced successfully.

# Conjunctive success: all four must be rejected
conj_power <- trial_success(r1 && r2 && r3 && r4)
#> ✔ Trial success function compiled and sourced successfully.

# FIBRONEER-style hurdle: dual-endpoint success for at least one dose.
# H1,H2 = high dose (primary, secondary); H3,H4 = low dose.
hurdle <- trial_success((r1 && r2) || (r3 && r4))
#> ✔ Trial success function compiled and sourced successfully.

# Gated secondary: H2 only earns value if H1 is also rejected
gated <- trial_success(r1 * (1 + 0.5 * r2))
#> ✔ Trial success function compiled and sourced successfully.

# Hurdle plus weighted incremental claims
composite <- trial_success((r1 && r2) * (1 + 0.3 * r3 + 0.2 * r4))
#> ✔ Trial success function compiled and sourced successfully.

# Composite utility: require H1 AND H2 rejection to get H3 and H4 value
gate_then_sum <- trial_success((r1 && r2) * (r3 + r4))
#> ✔ Trial success function compiled and sourced successfully.

# Weighted priorities (e.g., H1 gets weight 2 times H2 or H3)
weighted <- trial_success(2 * r1 + r2 + r3)
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
