# Call flow: `stopping_rule_gsd()` and its helpers

This document walks through the functions that create, display and settle an
efficacy stopping rule for a group sequential design (design record section
4.11), in the order the code executes them. All of them live in
`R/stopping_rule_gsd.R`. The truth table built from a settled rule
(`.gsd_stop_table()` in `R/check_gsd.R`) and the consumers that use it are
described in `dev/callflows/graph_optimise_gsd.md`.

Notation follows the design record: *m* hypotheses, *K* analyses ("looks"),
so *K - 1* interim analyses. A condition is a Boolean expression over the
rejection indicators `r1, r2, ..., rm`, evaluated after each interim analysis
on the hypotheses rejected so far (at that analysis or earlier). When it is
`TRUE` the trial stops there. The rule object does not know *m* or *K*; both
are settled when the rule meets a `multigrain_pvals_gsd` object.


## The call chain at a glance

1. `stopping_rule_gsd()` captures its argument with `rlang::enexpr()` and
   splits it into conditions with `.gsd_rule_capture()`, which uses
   `.gsd_rule_list_items()`, `.gsd_rule_list_item()`, `.gsd_rule_is_expr()`
   and `.gsd_rule_parse_string()`.

2. Each condition is deparsed for display, checked against the grammar by
   `.gsd_rule_validate()` (through `.gsd_rule_validate_call()` and
   `.gsd_rule_validate_symbol()`), and scanned for its hypothesis indices by
   `.gsd_rule_indices()`, which also requires at least one `r<i>`.

3. The object is returned: `conditions`, `labels`, `m_min`, `n_interim`.
   `print()` and `summary()` display it through `.gsd_rule_cat_detail()`.

4. In `graph_optimise_gsd()` and `calc_power_pvals_gsd()`,
   `check_stopping_rule_gsd()` checks the type, and
   `.gsd_check_stopping_rule_dims()` settles the rule against the p-values,
   returning exactly `K - 1` conditions for `.gsd_stop_table()`.


---


## `stopping_rule_gsd()` -- exported entry point

### Purpose

Creates a `multigrain_stopping_rule_gsd` object from one condition (applied at
every interim) or one condition per interim analysis. Called by the user.

### Inputs and output

- `condition`: an expression (`r1 && (r2 || r3)`), a single string
  (`"r1 and (r2 or r3)"`), a call to `list()` of expressions or strings, a
  literal `c("...")` of strings, or a list or character vector injected with
  `!!`.

Returns a list of class `multigrain_stopping_rule_gsd` with:

- `conditions`: a list of language objects (or scalar literals).
- `labels`: the conditions deparsed with `deparse1(width.cutoff = 500L)`.
- `m_min`: the largest hypothesis index referred to, the fewest hypotheses
  the rule can be used with.
- `n_interim`: `NULL` for a single condition to be recycled, otherwise the
  number of conditions supplied.

### How it works

1. `rlang::check_required(condition)`; then
   `.gsd_rule_capture(rlang::enexpr(condition))` returns the conditions and
   `n_interim`.

2. The labels are computed with `vapply(conditions, deparse1, ...)`.

3. For each condition, `.gsd_rule_validate()` checks the grammar and
   `.gsd_rule_indices()` returns its indices (aborting if there are none);
   `m_min` is the running maximum.

4. The object is assembled with `structure()`.


---


## `.gsd_rule_capture()` -- split the argument into conditions

### Purpose

Turns the captured argument into a list of conditions and `n_interim`.
Called by `stopping_rule_gsd()`.

### How it works

1. `.gsd_rule_list_items(cond_lang)` returns the elements if the argument is
   a list form. An empty `list()` aborts. Each element goes through
   `.gsd_rule_list_item()`, and `n_interim` is the number of elements, even
   when it is 1 (`list(r1)` is one condition for a design with one interim,
   not a recycled condition).

2. Otherwise, a character vector (a single string literal, or a vector
   injected with `!!`) is parsed element by element with
   `.gsd_rule_parse_string()`; an empty vector aborts. `n_interim` is `NULL`
   for a single string and the length otherwise.

3. Otherwise, anything `.gsd_rule_is_expr()` accepts is one condition with
   `n_interim = NULL`.

4. Anything else aborts, naming the type received.

### `.gsd_rule_list_items()`

Returns the elements when the argument is an injected list value, a call to
`list(...)`, or a call to `c(...)` whose arguments are all single string
literals; `NULL` otherwise. A `c()` call with anything else is left as an
expression and then fails the grammar as an unsupported function.

### `.gsd_rule_list_item()`

A single string is parsed with `.gsd_rule_parse_string()`; an expression is
returned as it is; anything else aborts naming the element's position.

### `.gsd_rule_is_expr()`

`TRUE` for a call, a symbol, or a numeric or logical value; a vector literal
passes here and is refused by the validator with a hint.

### `.gsd_rule_parse_string()`

Aborts on `NA`; maps the words `and`, `or` and `not` (case-insensitive,
word-bounded) to `&&`, `||` and `!`; parses with `str2lang()`, wrapping a
parse error in "Could not parse the stopping condition ...". Precedence is
therefore R's own, as in `trial_success_gsd()`.


---


## `.gsd_rule_validate()` -- grammar check

### Purpose

Stopping-rule twin of `validate_expr_symbols_gsd()`: after `!!` unquoting a
condition may contain only `r<i>` symbols, numeric or logical scalars, the
logical, comparison and arithmetic operators, and parentheses. Called by
`stopping_rule_gsd()` on each condition.

### How it works

It dispatches on the node:

- **Call**: `.gsd_rule_validate_call()`. The function must be a symbol in
  `gsd_rule_logic_ops` (`&&`, `||`, `!`), `gsd_rule_compare_ops`
  (`==`, `!=`, `<`, `<=`, `>`, `>=`), `gsd_rule_arith_ops` (`+`, `-`, `*`,
  `/`) or `(`; otherwise it aborts listing the allowed operators (interpolated
  as values, never literally, so a `<` cannot be read as cli markup) and
  noting that there are no discount tables in a stopping rule. Arguments are
  validated recursively.
- **Symbol**: `.gsd_rule_validate_symbol()`. A symbol matching
  `^r[1-9][0-9]*$` passes (so `r0` and `r01` do not). `t<i>` aborts with
  "decision times are not available to a stopping rule"; `s` aborts with
  "the stopping look is what the rule determines"; anything else aborts as an
  un-injected variable with a `!!` hint.
- **Numeric or logical value**: must have length 1; a vector aborts with a
  hint to give one condition per interim in a `list()` instead.
- **Anything else** aborts.


---


## `.gsd_rule_indices()` -- the hypotheses a condition refers to

### Purpose

Returns the integer indices of the `r<i>` symbols in a validated condition
and aborts when there are none. Called by `stopping_rule_gsd()`.

### How it works

`grep("^r[1-9][0-9]*$", all.names(expr), value = TRUE)`, then the numeric
suffixes. An empty result aborts: a condition that refers to no rejection
indicator cannot depend on the trial (so `TRUE` alone is refused).


---


## `print.multigrain_stopping_rule_gsd()`, `summary.multigrain_stopping_rule_gsd()` and `.gsd_rule_cat_detail()`

`print()` writes the class line, `summary()` an underlined "Stopping rule
(group sequential):" header; both return their argument invisibly and do
nothing for `NULL`. `.gsd_rule_cat_detail()` writes the shared body:
"Every interim: <label>" for a recycled condition, or one
"Interim k: <label>" line per condition, then
"Hypotheses referenced: up to r<m_min>". `print.multigrain_graph_optimal()`
calls `.gsd_rule_cat_detail()` directly for a result that carries a rule.

```
<multigrain_stopping_rule_gsd>
Interim 1: r1
Interim 2: r1 && (r2 || r3)
Hypotheses referenced: up to r3
```


---


## `is_stopping_rule_gsd()` and `check_stopping_rule_gsd()`

`is_stopping_rule_gsd()` is `inherits(x, "multigrain_stopping_rule_gsd")`.
`check_stopping_rule_gsd(stopping_rule, arg, call, allow_null = FALSE)`
returns invisibly for a rule (or for `NULL` when `allow_null` is `TRUE`) and
otherwise aborts with `rlang::stop_input_type()`, naming
`stopping_rule_gsd()`. Called with `allow_null = TRUE` by
`graph_optimise_gsd()`, `calc_power_pvals_gsd()` and `graph_optimal()`.


---


## `.gsd_check_stopping_rule_dims()` -- settle a rule against the p-values

### Purpose

Turns a rule into exactly one condition per interim analysis of a given
`multigrain_pvals_gsd` object, ready for `.gsd_stop_table(conditions,
pvals$m, pvals$K)`. Called by `graph_optimise_gsd()` and
`calc_power_pvals_gsd()` on the public path, once per call.

### Inputs and output

- `stopping_rule`: a `multigrain_stopping_rule_gsd` object.
- `pvals`: a `multigrain_pvals_gsd` object.
- `arg`, `call`: error context.

Returns a list of `K - 1` conditions (an empty list when `K = 1`).

### How it works

1. If `m_min > pvals$m` it aborts: the rule refers to `r<m_min>` but the
   p-values hold fewer hypotheses.

2. A single condition (`n_interim` is `NULL`) is recycled with
   `rep(conditions, K - 1)`; with `K = 1` this is the empty list, so the
   condition is accepted and has no effect.

3. `K - 1` conditions are returned as they are.

4. `K` conditions are accepted with a warning that the final analysis always
   ends the trial, naming the dropped condition, and the first `K - 1` are
   returned.

5. Any other number aborts, stating how many were given and how many are
   needed (or a single condition, not in a list).

### Why it is done this way

The rule is written without reference to a design, so the object stores no
`m` or `K` and the settling happens where the p-values are known. Accepting `K`
conditions with a warning rather than an error covers the natural mistake of
writing one condition per analysis.
