# Call flow: `transform_pvalues_gsd()` and its helpers

This document walks through every function involved in converting raw
group sequential p-values into repeated (or sequential) p-values, in the
order the code executes them. It is the readable companion to the
`transform` lane of the `gsd-workflow.excalidraw` diagram. All functions
live in `R/transform_pvalues_gsd.R`.

Notation follows the design record: *m* hypotheses, *K* analyses
("looks"), information fraction *t* in (0, 1] for each
hypothesis-by-look cell, spending function phi, nominal boundary
alpha-star, repeated p-value p^r, and sequential p-value p^s = running
minimum of p^r over looks. A "look" and an "analysis" are the same
thing. "Level" and "allocated level" both mean the alpha currently held
by a hypothesis under the graph.


## The call chain at a glance

1. `transform_pvalues_gsd()` receives the raw N-by-m-by-K array,
   validates it, and normalises every user-facing argument into its
   canonical form by calling the argument helpers `.gsd_info_frac()`,
   `.gsd_look_back()`, and `.gsd_spending_list()`.

2. For each of the *m* hypotheses, `.gsd_check_spending()` verifies that
   the spending function spends the whole level by full information,
   calling `.gsd_spend()` to evaluate it.

3. For each hypothesis, `.gsd_transform_hyp()` does the real work: it
   calls `.gsd_looks()` to find which analyses carry distinct
   information and when the hypothesis matures, calls
   `.gsd_warn_matured()` to warn about post-maturity data, then either
   short-circuits (single matured look) or builds a boundary table with
   `.gsd_boundary_table()`, checks the table with
   `.gsd_check_well_ordered()`, and inverts it per look with
   `.gsd_invert()`. The constant `gsd_grid_min` controls the floor of
   the grid and the floor of the repeated p-value.

4. `new_pvals_gsd()` assembles the transformed array and the per-hypothesis
   metadata into an S3 object of class `multigrain_pvals_gsd`.

5. `print.multigrain_pvals_gsd()` and `summary.multigrain_pvals_gsd()`
   display the result, with `summary` calling `.gsd_full_bounds()` to
   show the nominal boundaries at the full level. The display helpers
   `.gsd_hyp_labels()` and `.gsd_hyp_matrix()` format rows and columns.


---


## `transform_pvalues_gsd()` -- entry point

### Purpose

This is the sole exported function in the file and the only place in the
package where spending functions and `gsDesign` appear. It converts an
N-by-m-by-K array of raw (nominal) one-sided p-values into repeated
p-values, so that the downstream kernel can apply the familiar
fixed-sample rejection rule "p^r <= w * alpha" instead of comparing each
p-value against a spending boundary that depends on the current weight.
The caller is the user (or, in the pipeline, the optimiser setup code);
it calls the argument helpers, `.gsd_check_spending()`, and
`.gsd_transform_hyp()` for each hypothesis, then wraps the output with
`new_pvals_gsd()`.

### Inputs and output

- `pvals`: a three-dimensional numeric array with dimensions
  `c(nsim, m, K)` -- simulations by hypotheses by analyses. Values
  must lie in [0, 1]; `NA` is allowed (a look with no data). The
  function aborts if `pvals` is not three-dimensional or contains
  out-of-range values.

- `...`: enforced empty by `rlang::check_dots_empty()`.

- `info_frac`: either a numeric vector of length K (the same schedule
  for every hypothesis) or an m-by-K numeric matrix with `NA` where a
  hypothesis has no data at a look. If `NULL` (the default), the
  function falls back to `attr(pvals, "info_frac")`, which is what
  `simulate_pvalues_gsd()` attaches. If still `NULL` after the
  fallback, the function aborts.

- `spending`: a single spending function or a list of m spending
  functions. Each must accept `(alpha, t)` and return cumulative
  spend, either as a plain numeric vector or as an object with a
  `$spend` element (the convention used by `gsDesign` functions such as
  `sfLDOF`).

- `alpha`: overall one-sided significance level, a scalar in [0, 1].
  Default 0.025.

- `look_back`: a logical of length 1 or m. When `TRUE` for a
  hypothesis, the repeated p-values are replaced by their running
  minimum (the sequential p-value). Default `FALSE`.

- `grid_size`: a whole number (at least 2) giving the number of
  allocated levels in each boundary table. Default 1024.

The function returns a `multigrain_pvals_gsd` S3 object (see the
`new_pvals_gsd()` section below).

### How it works

1. The function captures the unevaluated `spending` argument as a label
   with `rlang::enquo()` and `rlang::as_label()` before any evaluation,
   so that the label stored in the output is the user's expression
   (for instance `"sfLDOF"`) rather than the function body.

2. It validates that `pvals` is a numeric (`check_double()`), that
   `alpha` is a number in [0, 1], and that `grid_size` is a whole
   number at least 2, then coerces `grid_size` to integer.

3. It reads `dim(pvals)` and aborts if the array is not
   three-dimensional. It also aborts if any non-`NA` value is outside
   [0, 1]. The three dimensions are unpacked as `n_sim`, `m`, and
   `n_look`.

4. If `info_frac` is `NULL`, it tries the `info_frac` attribute of
   `pvals`. It then calls `.gsd_info_frac()` to normalise the value
   into an m-by-K double matrix (broadcasting a vector, checking
   dimensions of a matrix). Similarly, `.gsd_look_back()` normalises
   `look_back` into a length-m logical, and `.gsd_spending_list()`
   normalises `spending` into a list of m functions.

5. A pre-check loop calls `.gsd_check_spending()` for each hypothesis,
   which evaluates the spending function at information fraction 1 and
   aborts if it does not spend the whole level (within a relative
   tolerance of 1e-6 * alpha).

6. The output array `out` is initialised to all 1s, with the same
   dimensions and dimnames as `pvals`. This is the "looks without data
   get p^r = 1" convention: any cell that the per-hypothesis transform
   does not overwrite already says "cannot reject". A list `tables` of
   length m is allocated for the per-hypothesis boundary metadata.

7. The main loop iterates over hypotheses. For hypothesis *i*, the raw
   data are extracted as an `n_sim`-by-`n_look` matrix (using
   `matrix()` to ensure the result is always a matrix even when
   `n_sim = 1`), and `.gsd_transform_hyp()` is called with that
   matrix, the i-th row of `info_frac`, the i-th spending function, and
   the scalar arguments. The returned `values` matrix is written into
   `out[, i, ]`, and the returned `table` is stored in `tables[[i]]`.

8. Finally, `new_pvals_gsd()` wraps the output array and all metadata
   into the S3 object, using `.gsd_spending_labels()` to produce a
   character label for each spending function.


---


## Argument helpers

These six small functions normalise user-facing arguments into the
canonical shapes the rest of the code expects. They are called only from
`transform_pvalues_gsd()` (or from the print/summary methods in the case
of the display helpers).


### `.gsd_info_frac(info_frac, m, n_look)`

Accepts either a numeric vector of length `n_look` (K) or an m-by-K
numeric matrix, and returns an m-by-K double matrix. If a vector is
given, it is broadcast to every hypothesis via
`matrix(info_frac, nrow = m, ncol = n_look, byrow = TRUE)`. The
function aborts if `info_frac` is `NULL` (with a message suggesting the
attribute fallback), non-numeric, the wrong length, or the wrong
matrix dimensions. It forces the storage mode to `"double"` before
returning.


### `.gsd_look_back(look_back, m)`

Accepts a logical (no `NA`) of length 1 or m, and returns a
length-m logical vector. A scalar is recycled to length m with `rep()`.
The function aborts if `look_back` is not logical (via `check_logical()`
from the rlang standalone), contains `NA`, or has the wrong length.


### `.gsd_spending_list(spending, m)`

Accepts a single function or a list of m functions, and returns a list
of m functions. A single function is wrapped with
`rep(list(spending), m)`. The function aborts if the input is neither a
function nor a list of functions, or if the list has the wrong length.


### `.gsd_spending_labels(spending, label)`

Produces a character vector of length m labelling each spending function
for display. Three cases, tried in order: if the list is named and every
name is non-empty, the names are returned; if there are at least two
distinct functions in the list (compared by identity), indexed labels
like `"sfLDOF[[1]]"` are returned; otherwise the captured label is
recycled for all m hypotheses.


### `.gsd_hyp_labels(m)` and `.gsd_hyp_matrix(x, m, n_look)`

`.gsd_hyp_labels()` returns `paste0("H", seq_len(m))`, a character
vector `c("H1", "H2", ...)`. `.gsd_hyp_matrix()` sets the dimnames of
an m-by-K matrix to hypothesis labels on rows and `"analysis 1"`,
`"analysis 2"`, etc. on columns, then returns it. Both are used by
`print()` and `summary()`.


---


## Pre-check: `.gsd_check_spending()` and `.gsd_spend()`

### `.gsd_check_spending(spending, alpha, hyp)`

**Purpose.** Validates that a spending function spends the full level by
full information -- the phi(a, 1) = a identity. If it does not, the
level that the boundary does not consume is silently wasted (the
hypothesis can never be rejected with that portion of its alpha). Called
once per hypothesis from the pre-check loop of
`transform_pvalues_gsd()`. It calls `.gsd_spend()`.

**How it works.** It evaluates `.gsd_spend(spending, alpha, t_look = 1)`
-- the spending function at the full level and a single information
fraction of 1 -- and checks that the returned value differs from `alpha`
by no more than `1e-6 * alpha` in absolute value. If it exceeds this
tolerance, the function aborts with a message naming the hypothesis, the
amount actually spent, and the intended level.


### `.gsd_spend(spending, alpha, t_look, hyp)`

**Purpose.** A thin wrapper that evaluates a spending function and
extracts the cumulative spend as a numeric vector, handling the two
return conventions. Called by `.gsd_check_spending()` and by
`.gsd_boundary_table()`.

**How it works.** It calls `spending(alpha, t_look)`. If the return
value is a list with a non-`NULL` element named `"spend"`, it extracts
that element (this handles `gsDesign` functions like `sfLDOF`, which
return a list with `$spend`). It then coerces to numeric with
`as.numeric()`. Finally it checks that the result has the same length as
`t_look` and contains no `NA`; if not, it aborts with a message about
the hypothesis and the expected length.


---


## Per-hypothesis core

This is the heart of the transform. For each hypothesis, the entry point
is `.gsd_transform_hyp()`, which orchestrates the sequence: determine
which looks matter, optionally short-circuit, build the boundary table,
check its monotonicity, invert it, copy matured values forward, and
optionally apply look-back.


### The constant `gsd_grid_min` (line 4)

```
gsd_grid_min <- 1e-14
```

This is the smallest allocated level represented in a boundary table. It
serves two purposes. First, it is the lower end of the log-spaced grid
of levels at which boundaries are computed. Second, it is the floor to
which a repeated p-value is clamped when the raw p-value falls below the
smallest positive entry in the boundary table.

The design record (section 4.1, "Rationale" paragraph on the floor, and
Appendix B, "Clamp safety") explains why this floor exists and why it is
1e-14 rather than zero. The kernel's weight-snapping rules can produce
allocations as small as 1e-4 * 1e-5 * 0.025 = 2.5e-11. If the floor
were zero, a row with a very small raw p-value (say 5e-11, which occurs
in about 0.9% of simulated rows at noncentrality 4) would receive a
repeated p-value of 0, which is below *any* positive allocation, causing
a spurious rejection. With the floor at 1e-14, the repeated p-value is
at least 1e-14, which is above zero but far below any realistic
allocation, so the spurious rejection is avoided. The residual
conservatism -- a non-rejection when the true repeated p-value lies
between 0 and 1e-14 and the allocation is in that same interval -- is
negligible because no realistic allocation is that small.


### `.gsd_transform_hyp(raw, t_row, spending, alpha, grid_size, look_back, hyp)`

**Purpose.** Implements steps 1 and 4 through 7 of design record
section 4.1 for a single hypothesis. Called from the main loop of
`transform_pvalues_gsd()` and calls `.gsd_looks()`,
`.gsd_warn_matured()`, `.gsd_boundary_table()`, `.gsd_invert()`. It
returns both the transformed values and the boundary table metadata.

**Inputs.**

- `raw`: an n_sim-by-n_look numeric matrix of raw p-values for this
  hypothesis.
- `t_row`: the information fraction row for this hypothesis (length
  n_look, may contain `NA`).
- `spending`: the spending function for this hypothesis.
- `alpha`, `grid_size`: passed through.
- `look_back`: a scalar logical for this hypothesis.
- `hyp`: hypothesis index, used in messages.

**Output.** A list with two elements: `values`, an n_sim-by-n_look
matrix of repeated (or sequential) p-values; and `table`, a list
containing `looks` (integer vector), `maturity` (integer scalar), and
optionally `grid` and `bounds` from the boundary table (or `NULL` if the
short-circuit path was taken).

**How it works.**

1. `.gsd_looks(t_row, hyp)` is called to determine which analyses carry
   distinct information. It returns `looks`, the integer indices of
   those analyses, and `maturity`, the analysis at which the hypothesis
   reaches full information (or the last available analysis if it never
   does). See the `.gsd_looks()` section below for details.

2. `.gsd_warn_matured(raw, hyp, maturity, n_look)` checks whether any
   raw p-values after the maturity analysis differ from the maturity
   value, and warns if so.

3. The output matrix `out` is initialised to all 1s. Analyses that are
   not in `looks` and are before maturity therefore receive p^r = 1
   ("cannot reject"), implementing the "looks without data" convention.

4. **Short-circuit branch.** If there is exactly one look in `looks`
   *and* its information fraction is at least 1, the boundary at level
   *a* is *a* itself (the spending function is exhausted in one look, so
   the boundary equals the level). No table is built (`tab` is set to
   `NULL`). The raw p-values at that single look are copied into every
   analysis from that look onward:
   `out[, seq.int(looks, n_look)] <- raw[, looks]`. This handles the
   K = 1 case and any hypothesis that matures at its first analysis,
   such as PFS in the manuscript's Example 5. The design record
   (section 4.1, step 1) notes that without this short-circuit, K = 1
   would not be bit-identical to the fixed-sample kernel, because
   `gsBound1` at t = 1 returns 0.012499999999999968 for level 0.0125,
   not 0.0125 exactly (verified by running `gsBound1` -- see
   Appendix B of the design record).

5. **Table branch.** Otherwise, `.gsd_boundary_table()` is called with
   the information fractions at the active looks (`t_row[looks]`) to
   build a monotone table of nominal boundaries against allocated
   levels. The returned `tab` has fields `grid` (the grid of allocated
   levels) and `bounds` (a grid_size-by-length(looks) matrix of
   boundaries).

6. A loop over the active looks inverts each column of the boundary
   table. For the l-th active look, `.gsd_invert(tab$bounds[, l],
   grid = tab$grid, p = raw[, looks[[l]]])` maps every raw p-value at
   that analysis to a repeated p-value. The result is written into
   `out[, looks[[l]]]`.

7. **Matured copy-forward.** If `maturity` is strictly less than
   `n_look`, every analysis after maturity receives the repeated
   p-value from the maturity analysis:
   `out[, seq.int(maturity + 1, n_look)] <- out[, maturity]`. The
   endpoint's data are complete, so its evidence does not change, but
   it may still be rejected later if recycling gives it more alpha.

8. **Look-back.** If `look_back` is `TRUE` and there are at least two
   analyses, a cumulative-minimum pass replaces each column from the
   second onward by `pmin(out[, k], out[, k - 1])`. This turns the
   repeated p-values into sequential p-values: the hypothesis can now
   be rejected on the strength of evidence from any earlier analysis.
   Because analyses without data already hold p^r = 1 (step 3), the
   running minimum correctly carries the last real sequential p-value
   through such gaps.

9. The function returns `list(values = out, table = c(list(looks, maturity), tab))`.


### `.gsd_looks(t_row, hyp)`

**Purpose.** Identifies which analyses carry distinct information for one
hypothesis and determines when it matures. Called by
`.gsd_transform_hyp()`.

**Inputs.** `t_row` is the information fraction row (length K, may
contain `NA`); `hyp` is the hypothesis index for error messages.

**Output.** A list with `looks` (integer vector of analysis indices
at or before maturity that have data) and `maturity` (integer, the
analysis index at which the hypothesis matures).

**How it works.**

1. `have <- which(!is.na(t_row))` finds the analyses with an
   information fraction. If none exist, the function aborts ("Hypothesis
   {hyp} has no analysis with an information fraction").

2. `t_have <- t_row[have]` extracts the non-`NA` fractions. It aborts
   if any are non-positive or if they are not non-decreasing (using
   `is.unsorted()`).

3. `full <- which(t_have >= 1)[1]` finds the first analysis (among
   those with data) that reaches full information. `maturity` is the
   original index of that analysis (in the K-length vector), or the
   last available analysis if no fraction reaches 1.

4. `looks` is the subset of `have` up to and including `maturity`.
   Analyses after maturity are excluded because their information
   fraction would be >= 1 again, and passing a repeated information
   fraction of 1 to `gsBound1` produces an unrejectable boundary (the
   spend increment is zero). The design record (Appendix B, first block
   and Appendix C item 1) documents this `gsBound1` edge case.

**Edge cases.** A hypothesis whose information fractions never reach 1
matures at the last analysis with data, and all its data-bearing
analyses are active. A hypothesis with data at only one analysis and
t >= 1 will yield `length(looks) == 1`, triggering the short-circuit in
`.gsd_transform_hyp()`. A hypothesis with `NA` at look 1 and data at
looks 2 and 3 simply has `looks = c(2, 3)` and look 1 keeps its
initialised p^r = 1.


### `.gsd_warn_matured(raw, hyp, maturity, n_look)`

**Purpose.** Warns the user if raw p-values supplied after maturity
differ from the value at the maturity analysis. The maturity value is
used regardless; the warning exists so the user knows their later data
were ignored. Called by `.gsd_transform_hyp()`.

**How it works.** If `maturity >= n_look`, there is nothing after
maturity, so the function returns `FALSE` invisibly. Otherwise, for each
analysis from `maturity + 1` to `n_look`, it compares the raw column
against `raw[, maturity]`. A value "differs" if it is not `NA` and
either the reference is `NA` or the two values are not equal. If any
simulation row differs at any post-maturity analysis, a single warning
is issued and the function returns `TRUE`. The loop exits on the first
differing analysis (it does not report *which* rows differ, only the
maturity index).


### `.gsd_boundary_table(t_look, spending, alpha, grid_size, hyp)`

**Purpose.** Steps 2 and 3 of design record section 4.1: builds a
monotone table of nominal boundaries against allocated levels, one row
per grid level and one column per active look. Called by
`.gsd_transform_hyp()` and calls `.gsd_spend()`, `gsDesign::gsBound1()`,
and `.gsd_check_well_ordered()`.

**Inputs.**

- `t_look`: the information fractions at the active looks (a numeric
  vector, no `NA`, positive, non-decreasing).
- `spending`: the spending function for this hypothesis.
- `alpha`: the overall level (the upper end of the grid).
- `grid_size`: the number of grid levels.
- `hyp`: hypothesis index for messages.

**Output.** A list with `grid` (a numeric vector of length `grid_size`
containing the log-spaced allocated levels) and `bounds` (a
grid_size-by-length(t_look) numeric matrix of nominal boundaries, where
`bounds[g, l]` is the largest p-value at look l that rejects when the
hypothesis holds level `grid[g]`).

**How it works.**

1. The grid of allocated levels is built as `exp(seq(log(gsd_grid_min),
   log(alpha), length.out = grid_size))`. This is a log-spaced grid
   from 1e-14 to alpha, meaning the grid points are equally spaced on
   the log scale. The design record (section 4.1, "Alternatives"
   paragraph) explains that a log grid gives uniform relative error
   under log-log interpolation, unlike a linear grid which wastes
   resolution where boundaries are large.

2. `lower` is set to `rep(-20, length(t_look))`. This is the lower
   (futility) bound on the z-scale, effectively negative infinity;
   `gsBound1` needs a lower bound to compute the upper boundary, and
   -20 on the standard-normal scale is far enough below zero that the
   probability of crossing it is negligible.

3. The `bounds` matrix is allocated as `grid_size` rows by
   `length(t_look)` columns, filled with `NA_real_`.

4. A loop over grid levels (g from 1 to `grid_size`) does two things at
   each level:

   a. `.gsd_spend()` evaluates the spending function at `grid_levels[g]`
      and the active information fractions, returning a vector of
      cumulative spend at each look. The spend increments
      `diff(c(0, cum_spend))` give the amount of level to spend at each
      look beyond what was spent at earlier looks.

   b. `gsDesign::gsBound1(theta = 0, I = t_look, a = lower, probhi =
      diff(c(0, cum_spend)))` computes the upper z-scale boundary at
      each look. `theta = 0` means the boundary is computed under the
      null hypothesis. `I` is the information-fraction vector. `a` is
      the lower boundary. `probhi` is the vector of incremental
      probabilities to spend. The call returns a list whose `$b`
      element is a vector of z-scale upper boundaries. These are
      converted to p-values on the probability scale by
      `pnorm(b, lower.tail = FALSE)` and stored as `bounds[g, ]`.

5. After the loop, `.gsd_check_well_ordered(bounds, hyp)` verifies
   that every column of `bounds` is non-decreasing in the grid index.
   If not, the spending function is not "well ordered" and the
   procedure is invalid.

**Why log-spaced.** At the low end of the grid, LDOF spending at three
looks produces boundaries that underflow to about 2.75e-89 at the first
look (verified by running `gsBound1` at grid levels below about 1e-6 --
see Appendix B of the design record). These duplicate, effectively
constant, entries are later removed by `.gsd_invert()`. A linear grid
would need orders of magnitude more points to achieve the same
resolution in the region where the boundaries are tiny.


### `.gsd_check_well_ordered(bounds, hyp)`

**Purpose.** Asserts that the boundary table is non-decreasing in each
column -- that is, as the allocated level increases, the nominal boundary
at every look is non-decreasing. This is condition (2) of Maurer and
Bretz (2013): if a spending function violates it, the sequentially
rejective graphical procedure is not valid, and the transform must stop.
Called by `.gsd_boundary_table()`.

**How it works.** For each column l of `bounds`, it computes a tolerance
as `sqrt(.Machine$double.eps) * max(column)` (approximately 1.49e-8
times the column maximum). It then checks whether any successive
difference `diff(column)` is less than the negative of this tolerance.
If so, the function aborts, naming the hypothesis and the analysis at
which the boundary decreases. Decreases smaller than the tolerance are
allowed, to accommodate floating-point noise from `gsBound1`.

The design record (section 10, item 9d) notes that the record says
"non-decreasing" while the code tolerates decreases below
`sqrt(.Machine$double.eps) * max(column)`, and flags this as a minor
review follow-up.


### `.gsd_invert(bounds, grid, p, floor, call)`

**Purpose.** Step 4 of design record section 4.1: inverts one column of
a boundary table to map raw p-values to repeated p-values. A repeated
p-value is the smallest allocated level at which the raw p-value would
cross its boundary. Because the boundary is monotone in the level, this
inversion is well defined. Called by `.gsd_transform_hyp()`, once per
active look.

**Inputs.**

- `bounds`: a numeric vector of length grid_size -- one column of the
  boundary table (the nominal boundaries at one look, indexed by grid
  level).
- `grid`: the vector of allocated levels (same length as `bounds`).
- `p`: a numeric vector of length n_sim -- the raw p-values to invert.
- `floor`: the value to clamp to below the table. Defaults to
  `gsd_grid_min` (1e-14).

**Output.** A numeric vector of length n_sim, the repeated p-values.

**How it works.**

1. Rows of the boundary table where the boundary is duplicated or
   non-positive are removed: `keep <- !duplicated(bounds) & bounds > 0`.
   This handles the underflow plateau at the low end of the LDOF table
   (the first-look boundary underflows to a constant 2.75e-89 across
   many grid levels, as described in section 4.1 step 4 of the design
   record). If fewer than two distinct positive boundaries survive, the
   function aborts, because interpolation needs at least two points.

2. The inversion is performed by log-log linear interpolation:
   `exp(stats::approx(log(bounds[keep]), log(grid[keep]), xout = log(p),
   rule = 1)$y)`. On the log scale, both the boundary values and the
   grid levels are monotone, and `approx` with `rule = 1` returns `NA`
   for any `xout` outside the range of the x-values (i.e. for p-values
   above or below the table).

3. Three clamp rules handle values outside the table and missing data:

   - If `p >= max(bounds)`, the raw p-value is above the largest
     boundary in the table (the boundary at alpha). No allocation up to
     alpha can reject it, so `out` is set to 1 ("never rejectable").

   - If `p < min(bounds[keep])`, the raw p-value is below the smallest
     positive boundary in the table. The repeated p-value is floored at
     `floor` (1e-14 by default) rather than set to 0, for the clamp
     safety reasons described above.

   - If `p` is `NA`, the repeated p-value is set to 1. This means a
     missing raw p-value is treated as "cannot reject", consistent with
     the all-1s initialisation in `.gsd_transform_hyp()`.

   All three clamp rules use `which()`, which silently skips `NA`
   indices, so they do not interfere with each other.

**Why log-log interpolation.** The design record (section 4.1,
"Alternatives") and Appendix B show that with 1024 grid points on
[1e-14, 0.025], the maximum relative error against `uniroot` at 200
random p-values is about 1.5e-6. Monte Carlo noise in the objective at
N = 1e5 is of order 1e-3, so the interpolation error is negligible.


---


## Result: `new_pvals_gsd()` and display methods

### `new_pvals_gsd(pvals, nsim, m, n_look, alpha, info_frac, look_back, spending, tables)`

**Purpose.** S3 constructor for the `multigrain_pvals_gsd` class. Called
once, at the end of `transform_pvalues_gsd()`.

**How it works.** It wraps its arguments in a named list and sets the
class to `"multigrain_pvals_gsd"` via `structure()`. Note that the
`n_look` parameter is stored under the name `K` in the returned object
(`K = n_look`); all other names match the parameter. The constructor has
default values (empty arrays, empty vectors) but
`transform_pvalues_gsd()` always passes explicit values.

The returned object contains:

- `pvals`: the N-by-m-by-K array of transformed p-values.
- `nsim`, `m`, `K`: the three dimension sizes.
- `alpha`: the overall level.
- `info_frac`: the m-by-K matrix of information fractions.
- `look_back`: the length-m logical of per-hypothesis semantics.
- `spending`: a character vector of length m, one label per hypothesis.
- `tables`: a list of length m, each element containing `looks`,
  `maturity`, and (unless the short-circuit applied) `grid` and
  `bounds`.


### `print.multigrain_pvals_gsd(x, ...)`

Prints a one-line class header via `cli`, a summary line with the
simulation count, hypothesis count, analysis count, and alpha, the
information-fraction matrix formatted by `.gsd_hyp_matrix()`, and the
look-back flag per hypothesis. If `x` is `NULL` it returns silently.


### `summary.multigrain_pvals_gsd(object, ...)`

Calls `print()` for the header, then adds a per-hypothesis detail table
(a data frame with columns for hypothesis label, active analyses,
maturity analysis, look-back flag, and spending label) and the nominal
boundaries at the full level (alpha) via `.gsd_full_bounds()`.


### `.gsd_full_bounds(object)`

**Purpose.** Extracts the nominal boundary of each hypothesis at each
analysis when the hypothesis holds the full level alpha. Used only by
`summary()`. It reads the boundary tables stored on the object rather
than recomputing them.

**How it works.** For each hypothesis, if the boundary table is `NULL`
(the short-circuit case), the boundary at every active look is simply
`alpha`. Otherwise, the last row of the `bounds` matrix
(`bounds[nrow(bounds), ]`) gives the boundaries at the largest grid
level, which is alpha. These values are placed into the active-look
columns. If the hypothesis matures before the last analysis, the
maturity boundary is copied forward to all later analyses (the same
logic as the matured copy-forward in `.gsd_transform_hyp()`). Analyses
before the first active look (where the hypothesis has no data) are left
as `NA`.
