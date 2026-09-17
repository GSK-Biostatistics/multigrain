# Call flow: `trial_success_gsd()` and its helpers

This document walks through every function involved in creating a
compiled group sequential trial-success utility, in the order the code
executes them. All functions live in `R/trial_success_gsd.R` except
`combine_arithmetic()`, which is reused from `R/trial_success.R`.

Notation follows the design record: *m* hypotheses, *K* analyses
("looks"). The decision time of hypothesis *i* is an integer in 0..K,
where 0 means "never rejected". The rejection indicator of hypothesis
*i*, written `ri` in code, is the boolean `ti > 0`. A discount table is
a user-supplied numeric vector of length K, indexed by the decision
time; the table's value at index 0 is always 0. The gain function
(the utility) is written as an R expression over `r1, r2, ..., rm`,
`t1, t2, ..., tm`, arithmetic, comparisons, logical connectives, and
discount-table calls; it is compiled to a C++ function that takes an
integer matrix of decision times (simulations by hypotheses) and returns
the mean utility.


## The call chain at a glance

1. `trial_success_gsd()` normalises `verbose`, validates the discount
   tables with `.gsd_gain_tables()`, settles K with `.gsd_gain_K()`,
   captures the user's expression, and calls `resolve_expr_gsd()` to
   turn it into a display string.

2. `resolve_expr_gsd()` either deparses the captured language object
   (after `validate_expr_symbols_gsd()` checks it) or passes a string
   through.

3. `new_trial_success_gsd()` receives the display string, the captured
   language object, the validated tables and the settled K. It counts
   hypotheses with `count_unique_indices_gsd()`, translates the
   expression to a C++ body with `replace_indices_gsd()` (which calls
   `parse_and_transform_gsd()` recursively), wraps it in a row loop
   over the decision-time matrix, compiles it with `Rcpp::sourceCpp()`,
   and returns the S3 object.

4. `print.multigrain_trial_success_gsd()` and
   `summary.multigrain_trial_success_gsd()` display the object, both
   calling `.gsd_gain_cat_detail()` for the K and tables block.


---


## `trial_success_gsd()` -- exported entry point

### Purpose

This is the sole exported function in the file. It creates a compiled
trial-success utility for a group sequential design: the user writes an
R expression over rejection indicators, decision times, comparisons and
discount tables, and the function compiles it to a C++ function that
can be evaluated row-wise on the integer decision-time matrix produced
by the kernel (`graph_shortcut_gsd()`). It is called by the user (or,
in the pipeline, by the optimiser setup code) and calls
`.gsd_gain_tables()`, `.gsd_gain_K()`, `resolve_expr_gsd()`, and
`new_trial_success_gsd()`.

### Inputs and output

- `objective`: an R expression or a character string encoding the gain.
  When an expression, `rlang::enexpr()` captures it before evaluation,
  so the user can inject values from the calling environment with `!!`.
  When a string, the word forms `and`/`or` are accepted.

- `...`: named numeric vectors (discount tables). Each name becomes a
  function in the expression language (e.g. `d = c(1, 0.75)` makes
  `d(t1)` valid). Every table must be the same length; that length is K.

- `K`: an optional positive integer giving the number of analyses. When
  tables are supplied, K is inferred from their common length; if also
  given explicitly, it must agree. When neither is given, K is left
  `NULL` and checked at first use (section 10 item 6 of the design
  record).

- `verbose`: controls whether a success message is printed after
  compilation; accepts `TRUE`/`FALSE`, `"info"`, `"detail"`, or
  `"silent"`.

The function returns a `multigrain_trial_success_gsd` object (see
`new_trial_success_gsd()` below).

### How it works

1. The `verbose` argument is normalised: `TRUE` becomes `"info"`,
   `FALSE` becomes `"silent"`, and a string is matched against the
   package's `verbosity_levels`.

2. `.gsd_gain_tables(rlang::list2(...))` validates the discount tables
   from `...` and returns a named list of double vectors.

3. `.gsd_gain_K(tables, K)` settles K: it is the common table length
   when tables are present (checked against an explicit K), or the
   explicit K alone, or `NULL` when neither is available.

4. `rlang::enexpr(objective)` captures the unevaluated expression. If
   the user wrote bare code, this is a language object in which injected
   (`!!`) constants have already been spliced; if the user passed a
   string, this is a character scalar.

5. `resolve_expr_gsd(expr_lang, table_names)` converts the captured
   value to a display string. When the input is a language object, it is
   validated and deparsed; the deparsed string has at most 15 significant
   digits per constant (R's `deparse1` default). When the input is
   already a string, it is returned as-is.

6. `new_trial_success_gsd()` receives the display string and, when the
   original was a language object, the captured language object itself.
   The language object is the one used for compilation, so that injected
   numeric constants keep their full double precision rather than the
   15-digit truncation of the display string. The display string is what
   the object shows to the user.


---


## `.gsd_gain_tables()` -- validate discount tables

### Purpose

Validates the discount tables collected from `...` in
`trial_success_gsd()`. Called only from that entry point.

### Inputs and output

- `tables`: a named list of numeric vectors, as captured by
  `rlang::list2(...)`.

The function returns a named list of double vectors that all share one
length. An empty list is returned if no tables were supplied.

### How it works

If `tables` is empty, the function returns an empty named list
immediately. Otherwise it checks four conditions, each aborting with
a cli error on failure.

- Every element must be named and have a non-empty name.
- Names must be unique.
- Each name must be a valid C identifier (starts with a letter, contains
  only letters, digits and underscores), must not match the rejection
  indicator pattern `r<digits>` or the decision-time pattern
  `t<digits>`, and must not be one of the word operators `and` or `or`
  (case-insensitive). This matters because the name becomes
  `static const double <name>_tab[]` in C++, and the `r<i>`/`t<i>`
  patterns and word operators are recognised by pattern before table
  names during parsing.
- Each table must be a non-empty numeric vector of finite values.
- All tables must have the same length (one value per analysis).

Finally, every table is coerced to double with `as.double()`.

### Edge cases

An attempt to name a table `r1`, `t2`, `and`, or `or` is caught and
reported with a message explaining the restriction. A table containing
`NA`, `Inf` or `NaN` is rejected. Duplicate names are rejected.


---


## `.gsd_gain_K()` -- settle the number of analyses

### Purpose

Determines the value of K from the discount tables and the explicit K
argument, or leaves it `NULL` when neither is available. Called only from
`trial_success_gsd()`.

### Inputs and output

- `tables`: the validated list of tables (may be empty).
- `K`: the user's explicit K, or `NULL`.

Returns an integer K, or `NULL`.

### How it works

If the user supplied an explicit K, it is checked as a positive whole
number (`rlang::check_number_whole(K, min = 1)`) and coerced to integer.
If no tables were supplied, K (explicit or `NULL`) is returned
immediately. If tables are present, K is their common length (read from
`length(tables[[1L]])`); if an explicit K was also given and disagrees,
the function aborts. The returned value is always an integer or `NULL`.


---


## `resolve_expr_gsd()` -- expression capture

### Purpose

Converts the captured `objective` to a display string. It is the group
sequential twin of the fixed-sample `resolve_expr()`. Called by
`trial_success_gsd()`.

### Inputs and output

- `expr_lang`: the value captured by `rlang::enexpr(objective)` -- a
  language object if the user wrote bare code, or a character scalar if
  the user passed a string.
- `table_names`: names of the discount tables, used during validation.

Returns a length-1 character string.

### How it works

If `expr_lang` is a character vector, it must be of length 1; the
string is returned as-is (validation happens later, when
`replace_indices_gsd()` parses it). If it is a language object,
`validate_expr_symbols_gsd()` is called to check that every node uses
the allowed grammar, and then `deparse1(expr_lang, width.cutoff = 500)`
produces the display string. Any other type aborts.


---


## `validate_expr_symbols_gsd()` -- grammar check

### Purpose

After `!!` unquoting, the expression may contain only the symbols
`r<digit>` and `t<digit>`, numeric or logical scalars, the arithmetic
operators, logical connectives, comparisons, parentheses, and calls of
a discount table to a single `t<digit>` symbol. Any other symbol means
the user forgot to unquote a variable, and any other operator means the
expression is not in the gain language. Called by `resolve_expr_gsd()`
on the captured language object, and again by `replace_indices_gsd()` on
the parsed AST.

### Inputs and output

- `expr`: a language object (the post-unquoting expression).
- `table_names`: discount table names.

Returns invisible `NULL` on success; aborts otherwise.

### How it works

The function dispatches on the type of node.

- **Call**: dispatched to `.gsd_gain_validate_call()`.
- **Symbol**: dispatched to `.gsd_gain_validate_symbol()`.
- **Numeric or logical scalar**: accepted if length 1. A vector literal
  (for instance an injected `c(1, 2)`) aborts with a message suggesting
  a discount table instead.
- **Anything else**: aborts.


### `.gsd_gain_validate_call()`

If the call's function name is one of the table names, the call must
have exactly one argument and that argument must be a bare symbol
matching `^t\d+$`. A table call like `d(t1 + 1)` or `d(r1)` is
rejected. Otherwise, the function name must be in the allowed set:
`+`, `-`, `*`, `/`, `(`, `&&`, `||`, `==`, `!=`, `<`, `<=`, `>`, `>=`.
This is a superset of what the fixed-sample validator allows (the
fixed-sample version has no comparisons and no table calls). If the
operator is allowed, all sub-arguments are validated recursively.

### `.gsd_gain_validate_symbol()`

A bare symbol must match `^[rt]\d+$` (a rejection indicator or a
decision time). If the symbol is a table name used without parentheses
(e.g. `d` instead of `d(t1)`), the error message says the table must
be applied to a decision time. Any other symbol gets an error suggesting
`!!` injection.


---


## `new_trial_success_gsd()` -- constructor

### Purpose

The S3 constructor. Receives a validated expression string, an optional
captured language object, validated tables and a settled K. Counts the
hypotheses, translates the expression to C++, compiles it, and returns
the object. Called only from `trial_success_gsd()`.

### Inputs and output

- `expr_string`: the display string (from `resolve_expr_gsd()`).
- `expr_lang`: the captured language object, or `NULL` when the user
  passed a string. When present, this is used for compilation instead of
  the display string, so that injected constants keep their full double
  precision.
- `tables`: the validated named list of discount tables.
- `K`: the settled integer, or `NULL`.
- `verbose`: one of `"info"`, `"detail"`, `"silent"`.

Returns a list with class
`c("multigrain_trial_success_gsd", "multigrain_trial_success")`
containing `func`, `m`, `K`, `objective`, `cpp_code`, and `tables`.

### How it works

1. `count_unique_indices_gsd(expr_string)` finds the largest index among
   all `r<i>` and `t<i>` symbols in the display string; that is `m`.

2. `replace_indices_gsd(expr_lang %||% expr_string, table_names)` is
   called on the captured language object if one exists, otherwise on
   the display string. It returns a C++ expression string that can be
   dropped into the body of a `for` loop.

3. For each discount table, a `static const double` array is generated.
   Index 0 is always `0.0` (the "never rejected" case); indices 1
   through K hold the table values, written with `.gsd_gain_cpp_number()`
   to ensure exact double round-tripping. For example,
   `d = c(1, 0.75)` produces
   `static const double d_tab[] = {0.0, 1.0, 0.75};`.

4. The full C++ source is assembled by `sprintf()`. It includes Rcpp,
   defines `std_min` as `std::min` (so that the deparsed R call
   `std_min(...)` is valid C++), places the table arrays above the
   function, and generates `double powerFunc(IntegerMatrix t)`. The
   function body checks once that the matrix has at least `m` columns,
   then loops over rows, accumulating the gain expression into `total`,
   and returns `total / n`.

5. `Rcpp::sourceCpp(code = cpp_code, env = local_env)` compiles the
   source into a fresh environment. The `powerFunc` binding from that
   environment becomes the `func` field of the returned object.

6. Unless `verbose` is `"silent"`, a success message is printed.

7. The object is assembled with `structure()`.

### Why there are two representations of the expression

The `objective` field stores the deparsed display string, which shows
constants to 15 significant digits (R's `deparse1` default). The
compiled function uses the captured language object when one is
available, so injected constants keep their full double precision. For
example, `!!(1/3) * r1` produces a display string
`0.333333333333333 * r1` (15 digits), but the compiled C++ contains
`0.3333333333333333` (16 digits -- the fewest that round-trip for this
value). This distinction is noted in the design record section 4.6,
[Rev 2026-09-17] paragraph "Two protections".

### What the compiled function checks and does not check

The generated C++ checks, once per call, that the integer matrix `t` has
at least `m` columns. If not, it calls Rcpp's `stop()` with a message
naming the actual and expected column counts. Without this check, a
narrower matrix would read past the end of a row and crash R.

What is deliberately *not* checked is the value of each time entry. The
discount-table lookup `name_tab[t(i, idx)]` is an unchecked C array
index. If a time value exceeds K, the index reads past the end of the
static array and returns garbage. If a time is `NA` (which is
`INT_MIN` in C), the index is negative and the read crashes R. The
design record (section 10 item 6) explains why this is acceptable: only
the P2 kernel (`graph_shortcut_gsd()`) may feed the compiled function,
and that kernel emits only values in 0..K with no `NA`. The P4 code
(`calc_power_pvals_gsd()` and `create_obj_func_gsd()`) is required to
assert `gain$K == pvals$K` and `!anyNA()` before evaluating.


---


## `count_unique_indices_gsd()` -- count hypotheses

### Purpose

Determines the number of hypotheses implied by the expression, which is
the largest index among the `r<i>` and `t<i>` symbols. Called by
`new_trial_success_gsd()`.

### Inputs and output

- `expr_string`: the display string (the same one stored in
  `objective`).

Returns a positive integer.

### How it works

The function finds all word-boundary-delimited matches of `[rt][0-9]+`
in the string. If none are found, it aborts with a message saying the
expression must reference at least one rejection indicator or decision
time. The numeric suffix of each match is extracted with `sub("^[rt]",
"", ...)` and coerced to integer. The maximum of the unique indices is
`m`. If the sequence 1..m has gaps (e.g. the expression uses `r1` and
`r3` but not `r2`), a warning is issued listing the missing indices.

This is the group sequential twin of `count_unique_indices()` from
`R/trial_success.R`, extended to match `t<i>` as well as `r<i>`.


---


## `replace_indices_gsd()` -- expression to C++ body

### Purpose

Translates the gain expression into a C++ expression string over the row
index `i` and the `IntegerMatrix t`. Called by `new_trial_success_gsd()`.

### Inputs and output

- `expr`: either a character string (the display string, used when the
  user passed a string) or a language object (the captured expression,
  used when the user wrote bare code).
- `table_names`: names of the discount tables.

Returns a character string of C++ code.

### How it works

1. If `expr` is a string, the word operators `and`/`or` (case-insensitive)
   are replaced by `&&`/`||` with `gsub()`, and the result is parsed
   with `str2lang()` to obtain an AST. If `expr` is already a language
   object, it is used directly.

2. `validate_expr_symbols_gsd()` is called on the AST. This is the
   second validation: the first happened in `resolve_expr_gsd()` on the
   captured expression, but when the input was a string, that call was
   skipped, so the string is validated here after parsing.

3. `parse_and_transform_gsd(ast, table_names)` recursively rewrites each
   node of the AST into an R call object that deparses to valid C++.

4. `deparse1(transformed$expr)` produces the C++ string. Backticks
   inserted by R's deparser (around names that are not valid R
   identifiers, such as the literal `1e-05`) are stripped with `gsub()`.
   Spaces are then inserted around the `/` operator, as the fixed-sample
   parser does; the code gives no reason and the output is valid C++
   either way, so this is cosmetic.

### Why native parsing instead of placeholders

The fixed-sample parser in `replace_r_indices()` rewrites `&&`, `||`,
`and` and `or` to the R infix operators `%AND%` and `%OR%` before
calling `str2lang()`. That cannot be copied here because R's custom
infix operators (`%op%`) have higher precedence than the comparison
operators, and comparison operators are non-associative, so
`str2lang("t1 == 1 %AND% t2 == 1")` is a parse error. The GSD parser
therefore maps only the word forms to `&&`/`||` and lets R parse the
native operators; the `&&` and `||` calls are then handled directly in
the transformer.

The consequence is that the GSD grammar has R's own operator precedence,
which differs from the fixed-sample parser in two ways. First, because
`%op%` binds tighter than `+`, the fixed-sample parser reads
`r1 + r2 && r3` as `r1 + (r2 && r3)`, while the GSD parser reads it as
`(r1 + r2) && r3`, which the bool type rule rejects with an error asking
for parentheses. Second -- and this case is silent rather than an error
-- `||` and `&&` have equal precedence under the fixed-sample
placeholders (both are `%op%`, grouped left to right), so
`r1 || r2 && r3` is read as `(r1 || r2) && r3`, while
`trial_success_gsd()` reads it as R does: `r1 || (r2 && r3)`, because
`&&` binds tighter than `||` in standard R. Both compile without warning
and give different gains.

To verify, running `trial_success_gsd(r1 || r2 && r3, K = 2)` on the
time matrix `[1, 0, 0]` (only H1 rejected) returns 1, because
`r1 || (r2 && r3)` evaluates to `1 || (0 && 0) = 1`. The old
placeholder reading `(r1 || r2) && r3` would return 0 on the same row.
The design record (section 4.6, "Parsing and precedence", [Rev
2026-09-17]) documents this difference and notes that the test file
asserts the R reading.


---


## `parse_and_transform_gsd()` -- recursive AST rewriter

### Purpose

Recursively rewrites each node of the parsed expression into an R call
object that deparses to valid C++, while tracking the type of each
sub-expression as either `"bool"` (a 0/1 quantity) or `"real"` (a
general double). Called by `replace_indices_gsd()`.

### Inputs and output

- `node`: one node of the AST (a call, symbol, numeric scalar, or
  logical scalar).
- `table_names`: names of the discount tables.

Returns a list with two elements: `expr` (the rewritten R call object)
and `type` (`"bool"` or `"real"`).

### How it works

The function dispatches on the type of node.

- **Call**: dispatched to `.gsd_gain_transform_call()`.
- **Symbol**: dispatched to `.gsd_gain_transform_symbol()`.
- **Numeric scalar**: converted to a C++ double literal by
  `.gsd_gain_cpp_number()` and wrapped as a symbol so it deparses
  correctly. Type is `"real"`.
- **Logical scalar**: `TRUE` becomes `1.0`, `FALSE` becomes `0.0`. Type
  is `"bool"`.


---


## `.gsd_gain_transform_call()` -- transform a call node

### Purpose

Handles one call node in the AST. Dispatches to the appropriate C++
translation depending on the operator.

### How it works

The function inspects the call's function name (as text) and takes one
of five branches.

**Parentheses.** If the operator is `(`, the inside is transformed
recursively and returned directly (parentheses are elided from the
output; the deparsed R call already groups correctly).

**Discount table.** If the operator name is one of `table_names`, the
call is `name(t<i>)`. The hypothesis index is extracted from the symbol
name (e.g. `t1` gives `idx = 0`), and the call is rewritten as
`name_tab[t(i, idx)]` -- an array lookup into the static C array
declared by `new_trial_success_gsd()`. Index 0 of the array is `0.0`,
so a never-rejected hypothesis contributes zero from every table. Type
is `"real"`.

**Logical operators** (`&&` and `||`). Both operands are transformed
recursively. If either operand's type is not `"bool"`, the function
calls `stop()` with a message: `` `&&` (AND) only allowed between
booleans `` or `` `||` (OR) only allowed between booleans ``. This is
the bool/real type rule, which is the same as in the fixed-sample
parser. The emitted C++ is:

- `&&`: `A * B`, which is 1 when both are 1 and 0 otherwise. Type
  `"bool"`.
- `||`: `std_min(double(1), A + B)`, which clamps the sum at 1. Type
  `"bool"`.

Both forms work because the operands are guaranteed to be 0 or 1.

**Comparisons** (`==`, `!=`, `<`, `<=`, `>`, `>=`). Both operands are
transformed recursively; any type is accepted. The comparison is wrapped
as `double(A op B)`, producing a 0.0 or 1.0 in C++. Type is `"bool"`.

**Arithmetic** (`+`, `-`, `*`, `/`). If there is only one operand
(unary minus or plus, detected by `length(transformed_args) == 1`), the
result is `"real"`. If there are two operands, `combine_arithmetic()`
from `R/trial_success.R` determines the type: `bool * bool` is `"bool"`;
every other arithmetic combination is `"real"`. The call is rebuilt with
the transformed operands.


---


## `combine_arithmetic()` -- type rule for arithmetic (reused from `trial_success.R`)

### Purpose

Determines the result type when an arithmetic operator is applied to two
sub-expressions of known type. This function lives in `R/trial_success.R`
and is reused unchanged by the GSD parser. It is called by
`.gsd_gain_transform_call()` for binary arithmetic.

### How it works

If the operator is not one of `+`, `-`, `*`, `/`, the result is `"real"`
(this branch is never reached in the GSD parser because those are the
only arithmetic operators that pass validation). For the four arithmetic
operators, the rule is:

- `bool * bool` returns `"bool"` (this is the AND product).
- `bool + bool`, `bool - bool`, `bool / bool` return `"real"`.
- Any combination involving `"real"` returns `"real"`.

The function returns a list with `type` and an `error_msg` field (always
`NULL` in practice; the field is a vestige of an earlier design where
some combinations were errors).

### The bool/real type rule in full

Combining everything the parser enforces, the complete type rule is:

| Expression form         | Type    | Errors when                             |
|-------------------------|---------|-----------------------------------------|
| `r<i>`                  | bool    |                                         |
| `t<i>`                  | real    |                                         |
| numeric literal         | real    |                                         |
| `TRUE` / `FALSE`        | bool    |                                         |
| comparison (`==`, etc.) | bool    |                                         |
| `name(t<i>)`            | real    |                                         |
| `bool && bool`          | bool    | either operand is real                  |
| `bool \|\| bool`        | bool    | either operand is real                  |
| `bool * bool`           | bool    |                                         |
| any other arithmetic    | real    |                                         |
| unary `-` or `+`        | real    |                                         |

The critical constraint is that `&&` and `||` require both operands to
be `"bool"`. This means `r1 && t1` errors (t1 is real), `r1 || 2`
errors (a numeric literal is real), and `r1 || d(t1)` errors (a table
lookup is real). But `(t1 == 1) && (t2 == 2)` succeeds because both
comparisons are bool, and `0.5 * r1` succeeds as real.


---


## `.gsd_gain_transform_symbol()` -- transform a symbol node

### Purpose

Translates a bare symbol (`r<i>` or `t<i>`) into C++ matrix access.
Called by `parse_and_transform_gsd()`.

### How it works

The symbol text is matched against two patterns.

- `r<i>` (rejection indicator): becomes `double(t(i, idx) > 0)`, where
  `idx` is the 0-based column index (`i - 1`). The `double()` cast
  ensures the boolean comparison result is a C++ `double`, so it mixes
  correctly with other double terms in the expression. Type is `"bool"`.

- `t<i>` (decision time): becomes `double(t(i, idx))`, casting the
  integer matrix element to double. Type is `"real"`.

Any other symbol falls through to a default branch that returns the node
unchanged with type `"real"`. This branch cannot be reached after
validation; it is kept for symmetry with the fixed-sample parser.


---


## `.gsd_gain_cpp_number()` -- numeric constant as C++ literal

### Purpose

Converts an R numeric scalar to a C++ double literal string that
round-trips exactly: reading the string back with a C++ compiler
produces the same double. Called when `parse_and_transform_gsd()`
encounters a numeric literal, and when `new_trial_success_gsd()` writes
discount-table values.

### How it works

The function starts by formatting the value with `sprintf("%.15g", x)`
(15 significant digits). It then checks whether `as.numeric(txt)` is
equal to the original value. If not, it tries 16 digits, then 17
(enough to round-trip any IEEE 754 double). Once a round-tripping
representation is found, the function checks whether the result contains
a decimal point or an exponent marker (`e` or `E`). If it does not (the
value is an integer like `2`), `.0` is appended to make it a C++ double
literal (`2.0`). If it already contains a decimal point or exponent
(like `0.75` or `1e-05`), it is left as-is.

This avoids the bug that the fixed-sample parser has with exponent
notation: R deparses `1e-05` as the string `1e-05`, and the fixed-sample
parser's blanket "append `.0` if no decimal point" rule would produce
the invalid C++ literal `1e-05.0`. The GSD parser checks for both `.`
and `e`/`E` before appending.

For example: `0.75` stays `"0.75"` (15 digits suffice); `1/3` becomes
`"0.3333333333333333"` (16 digits needed to round-trip); `2` becomes
`"2.0"`; `1e-05` stays `"1e-05"`.


---


## `print.multigrain_trial_success_gsd()`

### Purpose

Prints a compact representation of the gain object. Dispatched by R's
S3 method dispatch for objects of class
`multigrain_trial_success_gsd`.

### How it works

If `x` is `NULL`, the function returns invisibly. Otherwise it prints
the class vector (e.g.
`<multigrain_trial_success_gsd/multigrain_trial_success>`), the
objective string on the next line, and then calls
`.gsd_gain_cat_detail(x)` for the K and tables block. It returns `x`
invisibly. Example output:

```
<multigrain_trial_success_gsd/multigrain_trial_success>
0.4 * d(t1) + 1 * d(t2)
Analyses (K): 2
Discount table d(t): 1.00, 0.75
```


---


## `summary.multigrain_trial_success_gsd()`

### Purpose

Prints a fuller summary of the gain object. Dispatched by R's S3 method
dispatch.

### How it works

If `object` is `NULL`, the function returns invisibly. Otherwise it
prints an underlined header "Trial success function (group sequential):",
the objective string, and `.gsd_gain_cat_detail(object)`. It returns
`object` invisibly.


---


## `.gsd_gain_cat_detail()` -- shared display tail

### Purpose

Prints the K line and the discount-table lines, shared by `print()` and
`summary()`.

### How it works

It prints `"Analyses (K): "` followed by the formatted K (or `"not set"`
if K is `NULL`). Then, for each table name, it prints
`"Discount table <name>(t): "` followed by the formatted values
separated by commas. Returns invisible `NULL`.


---


## `is_trial_success_gsd()`

### Purpose

Tests whether an object is a `multigrain_trial_success_gsd`.

### How it works

Returns `inherits(x, "multigrain_trial_success_gsd")`. Because the
class vector is
`c("multigrain_trial_success_gsd", "multigrain_trial_success")`, objects
of this class also pass `is_trial_success()` from the fixed-sample file,
which checks for `"multigrain_trial_success"`.


---


## How each symbol translates to C++

For a concrete reference, the table below shows how each element of the
gain language is translated. The translations were verified by running
`trial_success_gsd()` and inspecting the `cpp_code` field.

| Expression element | C++ translation                                     | Type |
|--------------------|-----------------------------------------------------|------|
| `r1`               | `double(t(i, 0) > 0)`                               | bool |
| `t1`               | `double(t(i, 0))`                                    | real |
| `t1 == 1`          | `double(double(t(i, 0)) == 1.0)`                     | bool |
| `A && B`           | `A * B`                                              | bool |
| `A \|\| B`         | `std_min(double(1), A + B)`                          | bool |
| `d(t1)`            | `d_tab[t(i, 0)]`                                     | real |
| numeric `0.4`      | `0.4`                                                | real |
| `TRUE`             | `1.0`                                                | bool |
| integer `2`        | `2.0`                                                | real |


## Complete C++ for Example 5

The gain `0.4 * d(t1) + 1 * d(t2)` with `d = c(1, 0.75)` compiles to:

```cpp
#include <Rcpp.h>
using namespace Rcpp;

#define std_min std::min

static const double d_tab[] = {0.0, 1.0, 0.75};

// [[Rcpp::export]]
double powerFunc(IntegerMatrix t) {
    if (t.ncol() < 2) {
        stop("the decision-time matrix has %d column(s) but the trial success function refers to 2 hypotheses", t.ncol());
    }
    int n = t.nrow();
    double total = 0.0;

    for (int i = 0; i < n; i++) {
        total += (0.4 * d_tab[t(i, 0)] + 1.0 * d_tab[t(i, 1)]);
    }

    return total / n;
}
```

The `d_tab` array has three entries: index 0 is `0.0` (never rejected),
index 1 is `1.0` (rejected at look 1, full value), index 2 is `0.75`
(rejected at look 2, discounted). The body reads the decision time of
each hypothesis from the integer matrix `t`, uses it as an index into
`d_tab`, multiplies by the hypothesis weight, and sums. The function
returns the mean over all simulated trials.

Evaluated on a three-row time matrix where row 1 has H1 rejected at
look 1 and H2 at look 2, row 2 has only H2 rejected at look 1, and
row 3 has no rejections, the function returns
(0.4 * 1.0 + 1 * 0.75 + 0.4 * 0.0 + 1 * 1.0 + 0 + 0) / 3 = 0.7167.
