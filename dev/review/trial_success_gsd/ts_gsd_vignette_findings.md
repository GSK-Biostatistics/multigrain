# Review of open `trial_success_gsd()` vignette questions

Date reviewed: 2026-09-18  
Branch: `gsd-build`  
Reviewed HEAD: `45bc997`

## Scope

This review addresses items 2-9 in
`debug/plan/ts_gsd_vignette_brief.md`. It checks the current implementation
against:

- `R/trial_success.R`
- `R/trial_success_gsd.R`
- `R/transform_pvalues_gsd.R`
- `R/optimisation.R`
- `R/calc_power.R`
- `R/graph_optimal.R`
- `src/graph_shortcut_gsd.cpp`
- the associated `testthat` files
- `dev/gsd_design_record.md`
- `dev/gsd_progress.md`
- `docs/spiers_gain_function_arxiv.tex`

Focused examples were also rerun against the current branch. P4 has not been
implemented: `graph_optimise_gsd()`, `create_obj_func_gsd()`,
`calc_power_pvals_gsd()`, and `prune_graph_gsd()` do not exist yet.

During the review, uncommitted edits appeared in `R/trial_success.R`,
`R/trial_success_gsd.R`, `_pkgdown.yml`, and the vignette files. The changes in
the two R files are roxygen-only; they do not alter the executable behavior
described below.

## Recommended decisions

| Item | Finding | Recommendation |
|---|---|---|
| 2. Logical precedence | Confirmed fixed-sample parser bug. The expression is silently changed after R has parsed it. | File the issue below and fix `trial_success()` to use R precedence. Mixed real/logical operands should then be rejected, as in `trial_success_gsd()`. |
| 3. Unary operators | Confirmed fixed-sample parser bug affecting unary `-` and unary `+`. | File the issue below and handle one-argument arithmetic calls explicitly. |
| 4. String validation | Confirmed validation gap. Strings are a route around the public grammar. | Parse and validate strings with the same whitelist used for expression input before generating C++. |
| 5. GSD gain in fixed consumers | Confirmed silent wrong-result path in both `graph_optimise()` and `calc_power_pvals()`. | Fixed-sample consumers should reject GSD gain objects now. Preserve explicit type-aware validation so future GSD consumers and neutral result containers can accept them. |
| 6. `look_back` | Fully implemented per hypothesis in the transform. | Document the implemented behavior. Do not describe v1 as `FALSE`-only. |
| 7. Calendar time | The manuscript uses calendar time, while the package stores only look indices. The design record says supplied times are stored, but there is no way to supply them. | Add optional analysis-time metadata to the GSD p-value object and propagate it into P4 reporting. Keep the kernel and discount lookup indexed by look. |
| 8. `K = NULL` | Not inherently unsafe. It occurs only for table-free gains. The real memory-safety risk is a table-backed gain used with a different `pvals$K`, or direct calls with invalid time values. | Keep `K = NULL` for table-free, K-agnostic gains, clarify the documentation, and make P4 enforce all dimension/range contracts before evaluation. |
| 9. Discount-table shape | Current code accepts increasing and out-of-range values without warning, contrary to the manuscript's definition of a discount function. | Warn for increases or values outside `[0, 1]`; do not require the first value to equal 1. |

### Implications for the current uncommitted documentation draft

The new roxygen and `vignettes/articles/trial_success.qmd` changes are broadly
consistent with the observed behavior, but they should be revised in four
places before being treated as the intended API:

1. The precedence and string-validation sections should label those behaviors
   as bugs or temporary known limitations, not stable parser conventions.
2. The GSD roxygen currently says non-monotone discount tables compile "by
   design". That conflicts with the decision in item 9 to add a warning.
3. The vignette should distinguish analysis indices from calendar times and
   state that the latter are not currently stored.
4. The `K = NULL` wording should describe a table-free gain as K-agnostic and
   reserve the exact K-match requirement for a non-NULL `gain$K`.

---

## 2. `&&` and `||` precedence in `trial_success()`

### Current behavior

`trial_success()` does not preserve R's parse tree for logical expressions.
`replace_r_indices()` first replaces `&&` and `||` with custom infix operators
`%AND%` and `%OR%`, then reparses the resulting string:

```r
fixed_expr <- gsub("&&", "%AND%", fixed_expr, fixed = TRUE)
fixed_expr <- gsub("||", "%OR%", fixed_expr, fixed = TRUE)
ast <- str2lang(fixed_expr)
```

Custom `%...%` operators have different precedence from `&&` and `||`.
Consequently, the second parse can have a different tree from the expression R
originally captured.

The current branch produces:

```r
trial_success(r1 + r2 && r3)
```

```cpp
total += double(x(i, 0)) +
         double(x(i, 1)) * double(x(i, 2));
```

This is `r1 + (r2 && r3)`. R parses the user's expression as
`(r1 + r2) && r3`. Under the parser's own type rules, that R interpretation
should fail because `r1 + r2` is real-valued and `&&` requires two Boolean
operands.

Similarly:

```r
trial_success(2 * r1 || r2)
```

currently compiles as:

```cpp
2.0 * std_min(double(1), r1 + r2)
```

instead of rejecting the real-valued left operand of `||`.

`trial_success_gsd()` parses native `&&` and `||`, follows R precedence, and
rejects both examples with:

```text
`&&` (AND) only allowed between booleans.
```

or:

```text
`||` (OR) only allowed between booleans.
```

There is a second silent difference for expressions containing both logical
operators. The placeholders have equal custom-infix precedence and associate
left to right, while R gives `&&` higher precedence than `||`. Therefore:

```r
r1 || r2 && r3
```

is interpreted by R and `trial_success_gsd()` as:

```r
r1 || (r2 && r3)
```

but by `trial_success()` as:

```r
(r1 || r2) && r3
```

### Assessment

This is a correctness bug, not a documentation convention. It can silently
change the utility being optimised, so the optimiser can select a graph for a
different objective from the one the user wrote.

The existing fixed-sample tests currently assert some of the incorrect
regroupings in `tests/testthat/test-trial_success.R`. Those expectations should
be changed as part of the fix rather than treated as compatibility requirements.

### Correct behavior

`trial_success()` should follow R precedence, matching
`trial_success_gsd()`. In particular:

- `(r1 + r2) && r3` should be rejected because the left operand is real.
- `(2 * r1) || r2` should be rejected for the same reason.
- `r1 + (r2 && r3)` should remain valid because the parentheses explicitly
  make the logical subexpression an arithmetic term.
- `r1 || r2 && r3` should mean `r1 || (r2 && r3)`.
- String aliases `and`, `or`, `AND`, and `OR` should still be supported, but
  should be normalised to native logical operators before parsing rather than
  to custom infix placeholders.

The full issue draft is in [GitHub issue draft 1](#github-issue-draft-1).

---

## 3. Unary `+` and `-` in `trial_success()`

### Current behavior

Both of these fail:

```r
trial_success(-r1 + r2)
trial_success(+r1 + r2)
```

with:

```text
subscript out of bounds
```

The error is produced in `parse_and_transform()`. Arithmetic calls are treated
as though they always have a left and right operand:

```r
left_type <- transformed_args[[1]]$type
right_type <- transformed_args[[2]]$type
```

Unary `+` and unary `-` have only one transformed argument, so indexing
`[[2]]` fails.

`trial_success_gsd()` already has the required arity check:

```r
if (length(transformed_args) == 1L) {
    return(list(expr = new_call, type = "real"))
}
```

and successfully compiles unary minus.

### Assessment

This is a parser bug. The raw `subscript out of bounds` error is also an
implementation leak: it gives no useful information about what the user did.

### Correct behavior

The fixed-sample parser should explicitly support unary `+` and `-`, classify
their result as real, and retain the existing binary arithmetic rules. This
should work for expression and string input and for nested forms such as:

```r
-(r1 + r2)
r1 + (-r2)
```

A unary arithmetic result used directly as an operand of `&&` or `||` should
still fail the Boolean-operand check.

The full issue draft is in [GitHub issue draft 2](#github-issue-draft-2).

---

## 4. Validation of string input to `trial_success()`

### Current behavior

`trial_success()` has two materially different input paths.

For language input:

```r
trial_success(r1 + r2)
```

`resolve_expr()` calls `validate_expr_symbols()` before deparsing. The
validator allows only:

- `r<i>` symbols
- numeric or logical literals
- `+`, `-`, `*`, `/`
- `&&`, `||`
- parentheses

For character input:

```r
trial_success("r1 + r2")
```

`resolve_expr()` returns the string unchanged. `replace_r_indices()` parses it,
but `parse_and_transform()` treats unrecognised calls as real-valued calls
instead of rejecting them.

As a result, these currently compile:

```r
trial_success("r1 == 1")
trial_success("sqrt(r1 + r2)")
trial_success("abs(r1 - r2)")
```

Other unsupported expressions reach C++ generation and fail with compiler
diagnostics rather than package-level validation errors, for example:

```r
trial_success("r1 ^ 2")
trial_success("min(r1, r2)")
trial_success("r1 & r2")
```

The exact result depends on whether the generated text happens to be valid C++.
That means string input is not merely less well validated; it exposes a larger
and accidental language.

`trial_success_gsd()` does not have this gap. It parses a string and then calls
`validate_expr_symbols_gsd()` on the parsed tree before transformation.

### Why this matters

String input is the documented programmatic interface. It should not have
different semantics or error behavior from an equivalent unquoted expression.
The current difference creates three problems:

1. Unsupported utilities may compile and run.
2. Typographical errors may appear as long Rcpp/g++ failures.
3. A utility generated programmatically can behave differently from the same
   utility written directly in R code.

### Recommended behavior

Parse and validate character input against the same grammar as language input
before C++ generation.

A sound implementation sequence is:

1. Normalise documented string-only word operators (`and`/`or`) to native
   `&&`/`||`.
2. Parse with `str2lang()`.
3. Run the same symbol/operator validator used for captured language input.
4. Transform that validated tree to C++.

This work fits naturally with the precedence fix because both favor one native
parse path. After the change:

- `"r1 == 1"` should receive the normal unsupported-operator message unless
  comparison support is deliberately added to the fixed-sample grammar.
- `"sqrt(r1 + r2)"`, `"r1 ^ 2"`, and similar inputs should receive concise R
  validation errors.
- valid programmatic strings and the documented word forms should continue to
  work.

---

## 5. GSD trial-success objects in fixed-sample consumers

### Current behavior

`trial_success_gsd()` returns:

```r
class = c(
    "multigrain_trial_success_gsd",
    "multigrain_trial_success"
)
```

The shared parent class makes:

```r
is_trial_success(gsd_gain)
```

true. `check_trial_success()` consequently accepts the GSD object.

This is unsafe because the compiled functions have different contracts:

| Object | Compiled input | Meaning |
|---|---|---|
| `multigrain_trial_success` | logical `N x m` matrix | whether each hypothesis was rejected |
| `multigrain_trial_success_gsd` | integer `N x m` matrix | look of rejection, with 0 meaning never |

`graph_optimise()` validates with `check_trial_success()` and then evaluates its
objective on the logical matrix returned by `graph_shortcut()`.
`calc_power_pvals()` has the same problem through `.auto_name_custom_power()`
and `.eval_custom_power()`, which use `is_trial_success()` and call
`item$func(rej_mat)`.

Rcpp silently converts the logical matrix to an integer matrix. Therefore every
`TRUE` becomes decision time 1 and every `FALSE` becomes 0.

A rerun on the current branch gave:

```text
gain on the true decision-time matrix: 0.5833333
gain on the corresponding logical matrix: 0.6666667
```

`calc_power_pvals()` accepted the GSD object and returned the latter value
without an error or warning.

### Assessment

This is a high-impact silent correctness problem. It can affect both
optimisation and post-hoc power calculation. A user can obtain a plausible
number that systematically values every rejection as an analysis-1 rejection.

### Correct behavior now

Fixed-sample consumers should reject a
`multigrain_trial_success_gsd` object until the GSD-specific P4 path exists.
The error should explain that:

- the supplied gain requires decision times;
- the current function supplies only fixed-sample rejection indicators; and
- the user must use the future GSD-specific consumer.

This protection is required at least in:

- `graph_optimise()`
- `calc_power_pvals()` and its `custom_power` normalisation/evaluation path

### Class-design recommendation

Changing `is_trial_success()` globally to return `FALSE` for the subclass would
prevent the current silent path, but it is probably too blunt. The design
intends both gain types to be storable in the common
`multigrain_graph_optimal` result class, and P4 needs to recognise the GSD
subclass explicitly.

Prefer an explicit expected-type contract, for example:

- a fixed-sample predicate/check that accepts the parent class only when the
  object is not GSD;
- a GSD predicate/check using `is_trial_success_gsd()`; and
- an "either gain type" check only in neutral containers.

P4 must then dispatch deliberately:

- fixed gain -> `rejected`
- GSD gain -> `time`

Inheritance can remain for shared presentation and conceptual grouping, but it
must not determine the matrix passed to `$func`.

---

## 6. `look_back` in v1

### Current behavior

`look_back` is implemented in `transform_pvalues_gsd()`. It:

- accepts a logical scalar or a logical vector of length `m`;
- defaults to `FALSE`;
- is stored on the `multigrain_pvals_gsd` object;
- applies a running minimum to the transformed p-values for each selected
  hypothesis.

The rerun confirmed that `look_back = TRUE` produces the running minimum of the
ordinary repeated p-values.

It is intentionally not an argument to `trial_success_gsd()`. It changes the
testing procedure and therefore the decision-time matrix produced by the
kernel; the gain function only evaluates that matrix.

### Correct documentation

The vignette should state that v1 supports per-hypothesis look-back. It should
also explain:

- `FALSE` follows the ordinary repeated-p-value rule and is the default.
- `TRUE` permits later rejection using stronger evidence from an earlier
  analysis after alpha has been recycled.
- This can change both whether a hypothesis is rejected and the analysis at
  which the rejection is declared.
- The protocol or analysis plan must choose the intended convention; the
  package does not choose it from the data.

Any prompt text saying v1 supports only `look_back = FALSE` is stale.

---

## 7. Calendar time

### Manuscript semantics

The manuscript defines analyses at calendar times:

```text
s1 < ... < sK
```

and defines each decision time `tau_i` as one of those calendar times or
infinity. The discount function is explicitly a function of calendar time:

```text
d(tau_i)
```

For Example 5, the two decision times are 37 and 47 months.

### Package semantics

The implementation records only integer look indices:

```text
0 = never, 1 = first analysis, ..., K = final analysis
```

This is the right representation for the kernel and generated C++:

- it is compact;
- it indexes the discount table directly; and
- it avoids putting calendar-time data in the optimisation hot path.

The user manually evaluates the intended calendar-time discount function at
each analysis and passes the resulting values positionally:

```r
d = c(d(37), d(47))
```

For Example 5 this is:

```r
d = c(1, 0.75)
```

This is enough to compute the correct expected gain.

### What is missing

No current argument accepts calendar times:

- `simulate_pvalues_gsd()` stores only `info_frac`.
- `transform_pvalues_gsd()` accepts and stores `info_frac`, spending,
  `look_back`, `alpha`, and `K`.
- `trial_success_gsd()` stores `K` and numeric discount values.
- element names on discount vectors are dropped by `as.double()`.

The current `multigrain_pvals_gsd` fields are:

```text
pvals, nsim, m, K, alpha, info_frac, look_back, spending, tables
```

Therefore the statement in design-record section 4.4 that the package stores
calendar times "if supplied, as labels" is not implemented: there is no supply
path and no storage field.

### Why this is a real problem

The missing metadata does not make the positional gain calculation wrong when
the user has correctly precomputed the table. It does create important
problems:

1. **The analysis schedule is not auditable.** The object cannot show that
   positions 1 and 2 meant 37 and 47 months.
2. **A discount table can be paired with the wrong schedule.** The values alone
   do not identify the calendar times at which they were elicited.
3. **Future P4 reporting is underspecified.** The design record proposes
   `mean_decision_time`. Averaging indices gives a mean analysis number, not a
   mean calendar time. With unequally spaced analyses, converting an average
   index afterwards is mathematically wrong.
4. **Decision-time distributions lack meaningful labels.** Columns can be
   labelled "look 1" and "look 2", but not "37 months" and "47 months".
5. **The package cannot reproduce the manuscript's stated output scale.**
   The manuscript's `tau` is calendar time, while the package's `t<i>` is an
   index. That difference is manageable only if it is explicit and preserved.

### Recommended design

Calendar time belongs to the GSD design/p-value object, not primarily to
`trial_success_gsd()`. The same analysis schedule governs simulation,
transformation, kernel output, all gain functions, and all reporting.

Add optional analysis-time metadata along this path:

1. `simulate_pvalues_gsd(..., analysis_times = NULL)` may attach it to the raw
   array, as it currently does for `info_frac`.
2. `transform_pvalues_gsd(..., analysis_times = NULL)` should accept it
   directly or inherit the simulator attribute, then store it on
   `multigrain_pvals_gsd`.
3. Validate that it has length `K`, contains finite values, and is strictly
   increasing. The API also needs a documented unit convention or a separate
   label such as `"months"`.
4. Keep the kernel output as integer look indices.
5. In P4, always report a clearly named index-scale quantity such as
   `mean_decision_look`. Report `mean_decision_time` only when calendar times
   are available, by mapping each nonzero decision index to its supplied time
   before averaging.
6. Label `time_distribution` with `"never"` and the supplied analysis times
   when available.

The discount table should remain positional. Its values are the precomputed
`d(s_k)` values used in the hot path. Optionally preserving vector names can
help presentation, but names alone are not a sufficient design contract and
should not replace time metadata on the p-value object.

Adding `times =` only to `trial_success_gsd()` is not recommended. It would tie
the schedule to one gain object and allow that gain to be used accidentally
with p-values generated for a different schedule.

---

## 8. Reachable `K = NULL`

### When `K` is `NULL`

`trial_success_gsd()` determines `K` in this order:

1. the common length of supplied discount tables;
2. an explicit `K` argument;
3. otherwise `NULL`.

Therefore a gain with a discount table cannot have `K = NULL`; table length
always determines it. Nullable `K` occurs only for table-free expressions:

```r
trial_success_gsd(r1 + r2)
trial_success_gsd(t1 == 1 && t2 == 1)
```

Both rejection-only and time-based table-free expressions can currently have
`K = NULL`.

### What `K` controls

`K` does not control the number of matrix columns. The matrix passed to the
compiled gain has one column per hypothesis, and each cell contains a decision
index. The generated function already checks that there are at least `m`
columns.

`K` instead describes the allowed range of each cell:

```text
0, 1, ..., K
```

It is also the length of every generated C++ discount array after the leading
zero for "never" is added.

### Why `K = NULL` is not itself unsafe

Without a discount table, generated code uses decision times only in arithmetic
and comparisons. There is no C-array lookup indexed by `t<i>`. A table-free
gain can therefore be genuinely K-agnostic:

```r
trial_success_gsd(r1 + r2)
```

has the same meaning for any number of analyses.

Even:

```r
trial_success_gsd(t1 == 1)
```

has a well-defined reusable meaning: reward rejection at the first analysis,
whatever the total number of analyses.

The current code also allows an impossible comparison such as `t1 == 3` to be
built with `K = NULL` and later evaluated on a two-look time matrix. It simply
returns false for all rows. That is a potential user mistake, but it is not a
memory-safety failure.

### The separate table-index risk

The dangerous case is a table-backed gain:

```r
trial_success_gsd(d(t1), d = c(1, 0.75))
```

which generates an unchecked C++ lookup:

```cpp
d_tab[t(i, 0)]
```

The valid indices are 0, 1, and 2. A value greater than 2, or `NA` represented
as `INT_MIN`, reads outside the array. The design record reports that review
experiments produced garbage values and process crashes.

This is why P4 must verify that:

- `gain$m == pvals$m`;
- `gain$K == pvals$K` whenever `gain$K` is not `NULL`;
- the reshaped transformed p-values contain no `NA`;
- only the GSD kernel's integer `time` matrix is passed to a GSD gain; and
- kernel times are in `0..pvals$K`.

Because a table always sets `gain$K`, the P4 equality check covers the
table-length mismatch.

### What "checked at first use" should mean

The current roxygen text says that when `K` is not supplied it is "checked
against the p-values at first use." That is imprecise. If `gain$K` is `NULL`,
there is no stored value to compare.

The intended contract should instead be:

- A table-free gain may be K-agnostic and leave `K` unset.
- At first use, a GSD consumer obtains the effective number of analyses from
  `pvals$K`.
- If the gain has a non-NULL `K`, the consumer requires exact equality.
- If the gain has NULL `K`, the consumer may use it with that p-value object,
  subject to the normal `m`, matrix, missing-value, and range checks.

### Recommendation

Do not require `K` merely because P4 is not implemented yet. Keep `K = NULL`
for table-free gains so rejection-only gains remain reusable across designs.
Update the roxygen wording to describe K-agnostic gains rather than promising a
comparison that cannot occur.

P4 should implement the checks above before the first optimisation evaluation
and before each public power calculation. It should not rely only on the
generated function's current column-count check.

There is a separate API-safety decision for direct calls to the documented
`$func` field. The safest design would expose a checked R wrapper and keep the
raw compiled function internal, or add a range/NA scan to the generated
function. The optimisation path can still use a trusted internal fast path
after validating the p-value object and relying on the kernel's `0..K`
contract.

If a stricter, fully self-describing gain object is preferred, requiring `K`
whenever any `t<i>` symbol is present is defensible. It is not necessary for
correctness, however, and would reject useful K-agnostic expressions such as
`t1 == 1`. The more important requirement is explicit validation at the
consumer boundary.

---

## 9. Monotonicity and range of discount tables

### Current behavior

`.gsd_gain_tables()` currently checks only that each table:

- is numeric;
- is non-empty;
- contains finite values;
- has the same length as every other table.

It does not check the interpretation promised by the documentation and
manuscript.

Both of these compile without a warning:

```r
trial_success_gsd(d(t1), d = c(0.5, 1))
trial_success_gsd(d(t1), d = c(1.2, -0.1))
```

The first rewards a later rejection more than an earlier rejection. The second
uses multipliers outside the manuscript's `[0, 1]` range.

### Expected semantics

The manuscript defines a discount function as:

- valued in `[0, 1]`;
- non-increasing in calendar time; and
- zero for never rejected.

The implementation supplies the final condition automatically by prepending
zero at C++ index 0. Monotonicity should therefore be checked only across the
user's analysis values, not between the implicit "never" value and analysis 1.

The first analysis value need not equal 1. A user may legitimately put the
entire gain on an absolute or already-discounted scale.

### Recommended warning

Warn, rather than error, when a table:

- contains a value below 0 or above 1; or
- increases from one analysis to a later analysis.

The warning should name the table and the offending positions or values. For
example:

```text
Discount table `d` increases from analysis 1 (0.5) to analysis 2 (1).
i A discount table normally assigns no more value to a later decision.
```

and:

```text
Discount table `d` contains values outside [0, 1]: 1.2, -0.1.
i The manuscript defines discount multipliers on [0, 1].
```

A single construction may report both properties. A small numerical tolerance
for monotonicity is reasonable if table values are computed rather than typed,
but it should be documented and should not hide material increases.

This warning also catches the common mistake:

```r
d = c(0, 1, 0.75)
```

where the user incorrectly supplies the "never" zero themselves. The increase
from 0 to 1 is a strong signal that the table is shifted by one position. A
more specific hint can say that `trial_success_gsd()` supplies `d(0) = 0`
automatically.

The warning policy should be documented as semantic guidance, not as a
mathematical restriction on all possible utility lookup tables. If the package
eventually wants to support arbitrary time-value tables, those should be named
and documented separately from discount tables.

---

## Source map

| Topic | Primary implementation or evidence |
|---|---|
| Fixed-sample expression/string split | `R/trial_success.R:121-131` |
| Fixed-sample grammar validation | `R/trial_success.R:146-193` |
| Placeholder substitution and reparse | `R/trial_success.R:287-312` |
| Fixed-sample transformer and unary failure | `R/trial_success.R:382-470` |
| Existing tests that encode regrouped precedence | `tests/testthat/test-trial_success.R:276-350` |
| GSD table and `K` validation | `R/trial_success_gsd.R:172-281` |
| GSD grammar and string revalidation | `R/trial_success_gsd.R:326-466`, `623-636` |
| GSD generated C++ and table lookup | `R/trial_success_gsd.R:472-560`, `690-768` |
| GSD class inheritance and predicate | `R/trial_success_gsd.R:162-164`, `550-561` |
| Fixed optimiser validation and objective contract | `R/optimisation.R:88-154`, `R/objective_function.R:24-89` |
| Fixed power custom-function dispatch | `R/calc_power.R:114-230` |
| Neutral result-container check | `R/graph_optimal.R:24-49` |
| `look_back` implementation | `R/transform_pvalues_gsd.R:58-59`, `94`, `205-208`, `595-613` |
| GSD p-value object fields | `R/transform_pvalues_gsd.R:224-247` |
| Kernel decision-time contract | `src/graph_shortcut_gsd.cpp:54-84`, `99-124`, `186-188` |
| Current GSD tests for nullable `K` and class acceptance | `tests/testthat/test-trial_success_gsd.R:143-181`, `496-505` |
| Design-record decision time and calendar-time claim | `dev/gsd_design_record.md:140-152` |
| Design-record P4 and `K` obligations | `dev/gsd_design_record.md:191-205`, `334-345` |
| Manuscript calendar-time discount semantics | `docs/spiers_gain_function_arxiv.tex:379-434`, `1002-1106` |

---

## GitHub issue draft 1

### Title

`trial_success()` changes `&&`/`||` precedence when generating C++

### Body

#### Summary

`trial_success()` replaces `&&` and `||` with custom `%AND%` and `%OR%`
operators before parsing the objective. Custom infix operators have different
precedence from R's native logical operators, so some objectives are silently
regrouped and compiled with different semantics from the expression the user
wrote.

This can change the trial-success utility supplied to graph optimisation
without an error or warning.

#### Minimal reproduction

```r
ts <- trial_success(r1 + r2 && r3, verbose = "silent")
cat(ts$cpp_code)
```

The generated expression is equivalent to:

```r
r1 + (r2 && r3)
```

but R parses the input as:

```r
(r1 + r2) && r3
```

The latter should fail the existing type rule because `r1 + r2` is real-valued
and `&&` accepts only Boolean operands.

A second example also compiles when it should fail:

```r
trial_success(2 * r1 || r2, verbose = "silent")
```

It is generated as:

```r
2 * (r1 || r2)
```

There is also a silent `&&`/`||` disagreement:

```r
trial_success(r1 || r2 && r3, verbose = "silent")
```

R means:

```r
r1 || (r2 && r3)
```

whereas the placeholder parser groups the equal-precedence custom operators
left to right:

```r
(r1 || r2) && r3
```

`trial_success_gsd()` already follows R precedence and can be used as a
reference for the intended behavior.

#### Root cause

`replace_r_indices()` in `R/trial_success.R` rewrites native logical operators
to `%AND%` and `%OR%` before calling `str2lang()`. The replacement changes the
parse tree because `%...%` operators do not have the precedence of `&&` and
`||`.

#### Expected behavior

`trial_success()` should preserve R's operator precedence.

- `r1 + r2 && r3` should be parsed as `(r1 + r2) && r3` and rejected by the
  Boolean-operand check.
- `2 * r1 || r2` should be parsed as `(2 * r1) || r2` and rejected.
- `r1 + (r2 && r3)` should remain valid.
- `r1 || r2 && r3` should mean `r1 || (r2 && r3)`.
- Documented string aliases `and`/`or` and case variants should keep working.

#### Suggested implementation direction

Normalise word aliases to native `&&`/`||`, parse once with R, and handle
native logical calls in `parse_and_transform()`, as
`trial_success_gsd()` does. Do not use custom infix placeholders for parsing.

#### Acceptance criteria

- [ ] Native and string forms use R's `&&`/`||` precedence.
- [ ] `trial_success(r1 + r2 && r3)` errors because `&&` receives a real
      operand.
- [ ] `trial_success(2 * r1 || r2)` errors because `||` receives a real
      operand.
- [ ] `trial_success(r1 + (r2 && r3))` compiles.
- [ ] `trial_success(r1 || r2 && r3)` matches
      `r1 || (r2 && r3)` on a truth-table matrix.
- [ ] `and`/`or` string aliases retain the same semantics as `&&`/`||`.
- [ ] Tests that currently lock in the regrouped expressions are replaced with
      tests for R precedence.
- [ ] Existing correctly parenthesised objectives retain their results.

#### Scope

This issue concerns the fixed-sample `trial_success()` parser. It should not
change `trial_success_gsd()`, whose native parsing already has the intended
precedence.

---

## GitHub issue draft 2

### Title

`trial_success()` fails on unary `+` and `-` with `subscript out of bounds`

### Body

#### Summary

The fixed-sample trial-success parser assumes every arithmetic call has two
operands. Unary `+` and unary `-` therefore fail with an internal
`subscript out of bounds` error instead of compiling.

#### Minimal reproduction

```r
trial_success(-r1 + r2, verbose = "silent")
# Error: subscript out of bounds

trial_success(+r1 + r2, verbose = "silent")
# Error: subscript out of bounds
```

`trial_success_gsd()` already handles unary arithmetic:

```r
trial_success_gsd(-r1 + r2, K = 2, verbose = "silent")
# compiles
```

#### Root cause

In `parse_and_transform()` in `R/trial_success.R`, the arithmetic branch
unconditionally reads:

```r
left_type <- transformed_args[[1]]$type
right_type <- transformed_args[[2]]$type
```

A unary call has only one transformed argument, so `[[2]]` is out of bounds.

#### Expected behavior

Unary `+` and unary `-` should compile and return a real-valued expression.
Binary `+` and `-` should retain their current behavior.

Examples that should work:

```r
trial_success(-r1 + r2)
trial_success(+r1 + r2)
trial_success(-(r1 + r2))
trial_success(r1 + (-r2))
trial_success("-r1 + r2")
```

If a unary arithmetic expression is used directly as an operand of `&&` or
`||`, the existing Boolean-operand validation should still reject it.

#### Suggested implementation direction

Mirror the arity handling in `.gsd_gain_transform_call()`:

```r
if (length(transformed_args) == 1L) {
    return(list(expr = new_call, type = "real"))
}
```

Then apply the existing binary arithmetic type combination only when two
operands are present.

#### Acceptance criteria

- [ ] Unary minus compiles for expression and string input.
- [ ] Unary plus compiles for expression and string input.
- [ ] Nested unary forms compile and evaluate correctly.
- [ ] Binary arithmetic output remains unchanged.
- [ ] Unary arithmetic used as a logical operand still fails the Boolean type
      check.
- [ ] Invalid arities do not produce raw `subscript out of bounds` errors.
- [ ] Tests cover both generated C++ and values on logical matrices.

#### Scope

This issue concerns the fixed-sample `trial_success()` parser.
`trial_success_gsd()` already implements the intended unary behavior.
