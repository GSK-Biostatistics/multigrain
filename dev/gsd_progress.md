# GSD build progress: phases P0 to P3

Branch `gsd-build`, created 2026-09-16 from `8907be0` (head of `gsd-build-plan`). P0 to P2 were
an orchestrated session: one implementation agent per phase, the orchestrator re-running every
gate from a clean build before accepting it, and one adversarial review agent per phase. P3
(2026-09-17) was implemented and tested by the session itself, reviewed by the `gsd-reviewer`
agent before the gate commit, and documented by the `gsd-documenter` agent afterwards. Nothing
below is taken from an agent's report without the gate having been re-run.

## Status

| Phase | Commit | Gate (re-run before acceptance) | Review verdict |
|---|---|---|---|
| P0 Transform | `e4ec646` | `test-transform_pvalues_gsd.R`: 66 pass, 0 fail, 0 skip | clears gate; two record corrections, four minor follow-ups |
| P1 Simulator | `96dba9f` | `test-sim_pvals_gsd.R`: 29 pass, 0 fail, 0 skip (N = 1e6 block ran) | clears gate; one input-validation gap to decide |
| P2 Kernel | `d0b6556` | `test-graph_shortcut_gsd.R`: 64 pass, 0 fail, 0 skip; `test-RcppExports.R`: 65 pass, 2 pre-existing skips | clears gate; one record inaccuracy, two test-oracle gaps |
| P3 Gain | `0908406` | `test-trial_success_gsd.R`: 116 pass, 0 fail, 0 warn, 0 skip, stable across two runs (`load_all()`; no `src/` change); `test-trial_success.R`: 117 pass, 0 fail | clears gate; three findings, two fixed before the commit, one handed to P4 |

All three gates were re-run after `devtools::clean_dll()` and `devtools::install(quick = TRUE)`,
so no stale shared object was involved. `graphicalMCP` 0.3.0 was installed before P0, so every
oracle test executed; `skip_if_not_installed("graphicalMCP")` stays in the tests for CRAN.

The section 4.10 rule (new files only) holds. `git diff --stat 8907be0..d0b6556` touches
`DESCRIPTION` (two dependency lines), `NAMESPACE` (regenerated), `R/RcppExports.R` and
`src/RcppExports.cpp` (regenerated, the two new kernel wrappers only), and otherwise only new
files. `src/graph_shortcut.cpp` is byte-identical to `main`.

P3 keeps the rule: `git show --stat 0908406` touches `NAMESPACE` (one export, two S3 methods),
`man/trial_success_gsd.Rd` (roxygen) and otherwise only new files. `R/trial_success.R` is
byte-identical to the branch base; its only difference from `main` is the spelling commit
`324b7ca` that predates this work.

P4 (optimiser and post-processing) and P5 (documentation) have not been started.
The branch has not been pushed and no PR has been opened.

## What was built

### P0: `R/transform_pvalues_gsd.R`

Exported `transform_pvalues_gsd(pvals, ..., info_frac = NULL, spending, alpha = 0.025,
look_back = FALSE, grid_size = 1024L)` implementing all seven steps of record section 4.1.
Internals `.gsd_boundary_table()`, `.gsd_invert()`, `.gsd_check_spending()`, plus
`.gsd_looks()` (which P1 reuses so the simulator and the transform agree on which looks carry
distinct information), `.gsd_check_well_ordered()`, `.gsd_transform_hyp()`. Constructor
`new_pvals_gsd()` with `print()` and `summary()`. Floor constant `gsd_grid_min <- 1e-14`.
`DESCRIPTION` gained `Imports: gsDesign` and `Suggests: graphicalMCP (>= 0.3.0)`.
The two appendix scripts are kept as `dev/gsd_check1_tables.R` and `dev/gsd_check2_equivalence.R`.

Verified (implementer, then orchestrator, then reviewer independently):

- Both appendix scripts reproduce the record's numbers character for character on this machine
  (R 4.6.1, gsDesign 3.11.0), including 0/1200 mismatches in all four look-back configurations.
- Maurer and Bretz Table 1 boundaries to four significant figures at all four levels, cross-checked
  against `gsDesign::gsDesign()`.
- Grid inverse vs `uniroot()` at 200 random p-values: max relative error 6.7e-7 and 1.5e-6 at
  G = 1024, matching Appendix B to three significant figures.
- K = 1 short-circuit is bit-identical on 600 x 2 values including 100 rows with p in [1e-20, 1e-12].
- Values below the table map to exactly 1e-14, never 0; values at or above the top-of-grid boundary
  map to exactly 1.
- The defining property `p <= boundary(a)  <=>  p^r <= a` holds on 81,000 (row, look, allocation)
  comparisons with 0 disagreements, in both `<=` and strict `<` form.
- All seven named gsDesign spending families pass the well-ordering check; a family-switching
  function aborts naming the hypothesis at every switch level tried.
- Section 8 items 3 (transform half), 5, 6, 7, 8, 14, 15, 16 all pass.

### P1: `R/sim_pvals_gsd.R`

Exported `simulate_pvalues_gsd(power_nominal, ..., alpha = 0.025, corr_matrix =
diag(length(power_nominal)), info_frac, nsim = 1e5)`. One `mvtnorm::rmvnorm()` draw over the
union of (hypothesis, distinct look) pairs with the canonical joint model of section 4.7, placed
into an `N x m x K` array, matured columns copied forward, `NA` where a hypothesis has no data
before its first look, single attribute `info_frac` (m x K, NA-padded).

Verified:

- Reviewer's own covariance from the formula (m = 3, K = 3, one NA-padded row, one maturing early,
  rho = 0.4, N = 1e6): max correlation deviation 1.39e-3, max mean deviation 1.52e-3 (gate 5e-3).
- Final-look marginal rejection rate equals `power_nominal` within 1.2e-4 at N = 1e6, confirming
  the OPEN item 1 resolution applied (see below).
- K = 1, t = 1 output is `identical()` to `simulate_pvalues()` under the same seed, including the
  post-call `.Random.seed`, with default and non-identity correlation.
- Matured columns `identical()`; a hypothesis never reaching t = 1 has nothing copied.
- Regenerated `man/` from a clean `git archive` copy is line-for-line identical to the commit.

### P2: `src/graph_shortcut_gsd.cpp`

Internal `graph_shortcut_gsd(pvals, alpha, w, G, K)` and `graph_shortcut_gsd_parallel(pvals,
alpha, w, G, K, num_threads, grain_size)`, returning `list(rejected = N x m logical, time = N x m
integer, 0 = never)`. The cascade is the fixed-sample one verbatim (the reviewer's normalised
line-by-line comparison finds only the look loop, the rejected-flag guard, the time write, the
output types and `ops_per_row = K * m^3`). Test helper `tests/testthat/helper-gsd_reference.R`
carries the Appendix B direct implementation as `gsd_reference_direct()` and the `dim<-` reshape
`gsd_pvals_matrix()`.

Verified:

- Brief "Verification" claim 1: K = 1 `rejected` is `identical()` to `graph_shortcut()` on N = 1e4,
  m = 4, five random graphs, 2,200 rows below 1e-12 and 200 exact zeros; also through the full
  transform-reshape-kernel chain.
- Claim 2: Example 5 gives `time = c(2, 2)`; the oracle reproduces it only with
  `look_back = c(TRUE, FALSE)`, and the test asserts that `c(FALSE, FALSE)` does not.
- Claim 3: p = 5e-11 against a kernel-propagated allocation of 2.5e-11 does not reject; 1e-11 does.
- Claim 4: weight exactly 0 with p^r = 1e-14 or exactly 0 never rejects until alpha arrives, then
  rejects at exactly that look.
- Appendix B equivalence: 0/1200 and 0/1200 in FFF, TTT, TFT, FTF, with the record's local powers
  and mean decision times reproduced exactly; the reviewer's independent direct implementation
  agrees with both the kernel and the shipped helper on every cell.
- Beyond the reference: m = 4, K = 4, N = 2000, NA-padded and early-maturing rows, two zero
  initial weights, mixed look-back: 0/8000 mismatches. A hypothesis never reaching t = 1: 0/1500.
- Scan-order independence, decision-time semantics under look-back (time 3 not 1), graph-state
  carry-over across looks, the `sumrej == m - 1` update-skip path across looks: all pass.
- Parallel `identical()` to serial at 1, 2, 4, 8 threads and grain sizes 1, 7, 10000.
- Maurer and Bretz case study: H1, H2, H3 rejected at analysis 2, H4 retained, matching the oracle
  on the rejected set.
- Serial 1.5 ms vs 4 threads 1.0 ms per call at N = 1e4, m = 4, K = 3 (threading overhead
  dominates at this size, as expected).

### P3: `R/trial_success_gsd.R`

Exported `trial_success_gsd(objective, ..., K = NULL, verbose = multigrain_verbosity())`. The
expression may use `r<i>` (rejected at any analysis), `t<i>` (analysis of rejection, 0 = never),
`+ - * /`, `&& || and or`, `== != < <= > >=`, parentheses, numeric literals or `!!`-injected
values, and discount tables passed by name through `...` and applied as `name(t<i>)`. The
compiled function is `double powerFunc(IntegerMatrix t)` over the kernel's `time` matrix:
`r1` becomes `double(t(i, 0) > 0)`, `t1` becomes `double(t(i, 0))`, a comparison becomes
`double(A op B)`, `&&` and `||` emit the fixed-sample `A * B` and `std_min(double(1), A + B)`,
and `d(t1)` becomes `d_tab[t(i, 0)]` with `static const double d_tab[] = {0.0, ...}`. The
object has class `c("multigrain_trial_success_gsd", "multigrain_trial_success")` and fields
`func`, `m`, `K`, `objective`, `cpp_code`, `tables`, with `print()` and `summary()` methods and
the internal predicate `is_trial_success_gsd()`. Every function carries a call-flow comment;
the readable companion is `dev/callflows/trial_success_gsd.md`.

Internals, all new and `_gsd`-suffixed so `R/trial_success.R` is untouched: `.gsd_gain_tables()`
(names must be C identifiers other than `r<i>`, `t<i>`, `and`, `or`; finite numeric; one common
length), `.gsd_gain_K()`, `resolve_expr_gsd()`, `validate_expr_symbols_gsd()` with
`.gsd_gain_validate_call()` and `.gsd_gain_validate_symbol()`, `new_trial_success_gsd()`,
`count_unique_indices_gsd()`, `replace_indices_gsd()`, `parse_and_transform_gsd()` with
`.gsd_gain_transform_call()`, `.gsd_gain_transform_symbol()` and `.gsd_gain_cpp_number()`.
`combine_arithmetic()` from the fixed-sample file is reused unchanged.

Verified (session, then reviewer independently on its own time matrices):

- The four record 4.6 gains (Example 5; dual primaries with same-look bonus; co-primaries;
  PFS/OS) compile and return values `identical()` to R on a time matrix covering every
  `(t1, t2)` pair in 0..2, with R accumulating in double in row order as the C++ loop does.
- `d(0)` is exactly 0; `d = c(1/3, 2/3)` returns exactly `1/3` and `2/3`.
- `K`: inferred from tables; explicit `K` agreeing accepted, disagreeing errors; tables of
  different lengths error; no tables and no `K` gives `NULL`; `K = 0` and `K = 1.5` error.
- Section 8 item 10: `(t1 == 1) && r2` and `d(t1)` with `d = c(1, 0.75)` compile; `d(0)` is 0.
- Section 8 item 11 in its P3 form and brief claim 5: two compilations of the same expression
  give `identical()` `cpp_code` and `identical()` values on the hand matrix and on 1000 random
  rows.
- Twenty-five error paths give the intended message, including the `!!` hint, table misuse
  (`d(r1)`, `d(1)`, `d(t1 + 1)`, `d(t1, t2)`, bare `d`), forbidden table names, `&&`/`||` with
  a real operand, an injected vector literal, and expressions with no `r`/`t` symbol.
- Precedence: `t1 == 1 && t2 == 1` unquoted and as `"t1 == 1 and t2 == 1"`; all six comparison
  operators against literals and arithmetic; `r1 || r2 && r3` is `r1 || (r2 && r3)`.
- `t10`, `r12` index correctly with `m = 12`; tables named `t`, `n`, `i`, `total`, `std_min`,
  `d_1`, `d2` compile (the `_tab` suffix avoids the generated identifiers); `random` and `orbit`
  are not mangled by the `and`/`or` word mapping; `1e-05` and `1e+06` literals compile; unary
  minus works; `r1 + r2` equals `(t1 > 0) + (t2 > 0)` bitwise; `check_trial_success()` accepts
  the object.

## Decisions applied in P3 (need user confirmation)

5. **Native parsing instead of the placeholder copy.** Forced, see "Roadmap changes" below.
6. **`K` may be `NULL`** when neither a table nor the argument fixes it (record 10 item 6 as
   written). A stricter rule is a one-line change in `.gsd_gain_K()`.
7. **Table names are restricted** to C identifiers that are not `r<i>`, `t<i>`, `and`, `or`.
   The record did not specify; without the restriction the generated C++ would not compile or
   the name would be captured by the symbol grammar.
8. **Two extra validations** beyond the record: a table call must take a single bare `t<i>`
   (so `d(t1 + 1)` and `d(r1)` error rather than compiling to nonsense), and a literal must be a
   scalar (an injected vector such as `!!d * r1` is caught before compilation with a hint to
   pass it as a table).
9. **Two protections added after the review** (see below): the compiled function errors when
   the time matrix has fewer than `m` columns, and constants are compiled at full precision.

## Decisions applied by the orchestrator (need user confirmation)

1. **OPEN item 1, simulator inputs.** The record's own 4.7 signature (`power_nominal` and
   `corr_matrix` describing the final-look statistics, exactly as `simulate_pvalues()`) was
   implemented. The alternative (noncentrality parameters directly) was not. If the user prefers
   noncentrality, P1 needs a small change and its tests a re-derivation.
2. **Function names** (`transform_pvalues_gsd()`, `simulate_pvalues_gsd()`, `graph_shortcut_gsd()`)
   were taken from the user's phase list and are now in the API.
3. **`grid_size` default 1024**, per the record.
4. `corr_matrix` was given the default `diag(length(power_nominal))` to mirror `simulate_pvalues()`;
   the record's 4.7 shows no default.

## Roadmap changes

The plan's structure (P0 to P5, gates, LOCKED decisions) is unchanged. The record and the brief
have been corrected where the implementation showed a stated fact to be wrong or underspecified;
every edit is marked `[Rev 2026-09-16]` in `dev/gsd_design_record.md`.

Record corrections:

- **Table 2 is reproduced to about 3.5 significant figures, not four.** gsDesign gives 0.168214,
  0.131537 and 0.0116485 where the paper prints 0.1683, 0.1316 and 0.0117 (ADDPLAN vs gsDesign,
  deviations up to 8.5e-5). The record's own Appendix A output shows this. The test asserts the
  Appendix A gsDesign values, which is what the gate specifies. Sections Part I, 4.1 rationale and
  6 P0 corrected.
- **The graphicalMCP oracle agreement is absolute 1e-5, not relative, and only at K = 2.**
  `repeated_p()` roots to `tol = 1e-6` on the alpha scale, clamps to [1e-6, 1 - 1e-6], and at
  K = 3 is not reproducible run to run (randomised `pmvnorm`). Measured: package vs high-precision
  uniroot 4.0e-8 absolute; oracle vs the same 2.7e-5 absolute, 0.24 relative. The P0 test uses
  K = 2 and p in [1e-5, 1e-2] (max absolute difference 2.4e-7). Sections 6 P0 and 7 corrected.
- **`src/RcppExports.cpp` is also regenerated.** Section 4.10 listed only `R/RcppExports.R` and
  `DESCRIPTION`; the P4 no-change gate must allow both RcppExports files.
- **graphicalMCP conventions** the oracle tests must respect, added to 4.9: `spending_fn` must
  return a plain numeric vector (passing `gsDesign::sfLDOF` directly silently gives repeated
  p-values of 1 with only a coercion warning); `p` must carry no column names; `decision_at` is
  written for tested-but-not-rejected hypotheses, so compare on the rejected set only; the oracle
  uses `<=`; our repeated p-values are capped at 1 above `alpha` where the oracle reports the true
  value.
- **Section 6 P2** now names the helper `gsd_reference_direct()` and records its Appendix B
  restrictions (every hypothesis must reach t = 1; no NA in `tmat`).
- **Appendix B notes**: the `P <- lapply(...)` block simulates hypotheses independently (no
  cross-hypothesis rho), which is fine for the equivalence check but is not the joint model of 4.7;
  and `matrix(p[, kk], N, K - kk)` warns when `kk == K`.
- **Section 10**: items 1, 2, 3 marked closed; new items 9 to 11 added (below).

Brief corrections (`dev/implement_gsd.md`): no `CLAUDE.md` exists; `NOT_CRAN=true` is required
for `skip_on_cran()` gates run from `Rscript`; kernels are internal so `test_file()` needs
`package = "multigrain"` or `load_all()`; the graphicalMCP conventions above; both RcppExports
files are allowed changes; P0 to P2 marked done.

### P3 (2026-09-17), marked `[Rev 2026-09-17]`

- **The P3 parser cannot be "a copy with three additions".** The fixed-sample
  `replace_r_indices()` rewrites `&&`/`||`/`and`/`or` to `%AND%`/`%OR%` before `str2lang()`.
  R parses `%op%` above `==` and comparison operators are non-associative, so
  `str2lang("t1 == 1 %AND% t2 == 1")` is a parse error (verified, R 4.6.1); none of the
  supplement gains in 4.6 could compile. `trial_success_gsd()` maps only the word forms to
  `&&`/`||` and parses natively, handling `&&`/`||` calls in the transformer with the same bool
  rule and the same emitted C++. Consequence: R precedence. `r1 + r2 && r3` errors (real `&&`
  bool) where the fixed parser gives `r1 + (r2 && r3)`; `r1 || r2 && r3` is `r1 || (r2 && r3)`
  where the fixed parser gives `(r1 || r2) && r3` silently. Record 4.6 gained a "Parsing and
  precedence" paragraph; section 6 P3 and section 9 were corrected.
- **Numeric literals.** The fixed parser appends `.0` to any literal without a decimal point,
  so `1e-05` becomes the invalid C++ `1e-05.0` (verified). The GSD parser keeps literals that
  already carry a point or an exponent. Added to section 9 as a pre-existing issue.
- **Record 10 item 6 extended** with the P4 obligations that follow from the unchecked table
  index (below).

Brief corrections: P3 marked done; the parsing lesson; the P4 obligations from P3
(`gain$K == pvals$K`, `is_trial_success_gsd()` to tell the two gain classes apart, and that
cli messages must interpolate comparison operators as values because a literal `<` inside cli
markup is read as a delimiter and the message fails to format).

## P3 review findings and their disposition

The reviewer confirmed the parse-error claim and judged the native-parsing divergence forced
and minimal. Three findings:

1. **Unchecked table index.** `name_tab[t(i, idx)]` is a raw C array index. The reviewer
   demonstrated a segmentation fault on an `NA` time (`INT_MIN`), a segmentation fault on a time
   matrix with fewer columns than `m`, and garbage (0.56, then 1.1e272) for times above `K`.
   Fixed in the commit: the generated function checks `t.ncol() >= m` once per call and errors
   otherwise (tested). Not fixed, by decision: per-element checks of the time value would put a
   branch in the objective for a matrix only the P2 kernel is meant to produce; P4 must assert
   `gain$K == pvals$K` when `gain$K` is set and `!anyNA()` on the kernel input (record 10
   item 6 and 11; brief "Verification").
2. **Injected literals did not round-trip.** `!!(1/3) * r1` evaluated to 0.33333333333333298:
   the deparse that produces the display string carries 15 significant digits. Fixed in the
   commit: the captured language object is compiled directly (the deparsed string is display
   only), and every constant and table value is written with the fewest significant digits that
   round-trip (`.gsd_gain_cpp_number()`); `!!(1/3) * r1` now returns exactly `1/3` (tested).
3. **Silent coercion at the C++ boundary.** `func()` accepts a double or logical matrix and
   truncates, so a fixed-sample rejection matrix handed to a GSD gain returns a number rather
   than an error. Not changed in P3; the brief now tells P4 to dispatch on
   `is_trial_success_gsd()`.

The reviewer also noted that the record and brief were edited while it worked (the
`[Rev 2026-09-17]` amendments above); it reviewed against the original 4.6 text and found no
contradiction.

## Follow-ups surfaced by the reviews (not fixed; none violates a LOCKED decision)

Recorded as new section 10 items 9 to 11.

- `transform_pvalues_gsd()`: `alpha` bounds are inclusive (0 and 1 pass validation and fail inside
  gsDesign); a spending function returning NA gets a misleading length message; an unnamed list of
  spending functions prints deparsed code as the label in `summary()`; the well-ordering check
  tolerates decreases below `sqrt(.Machine$double.eps) * max(column)` (about 3.7e-10 absolute)
  where the record says strictly non-decreasing. All one-line fixes; the record should say whether
  the tolerance is wanted.
- `simulate_pvalues_gsd()`: a duplicated interior information fraction below 1 (e.g.
  `c(0.5, 0.5, 1)`) is accepted, draws a singular covariance, and then makes the transform abort
  with "fewer than two distinct positive boundaries". The record should either forbid it at input
  or state that it is tolerated. `alpha = 0` and `power_nominal > 1` are accepted as in
  `simulate_pvalues()`.
- Kernel: `NA` or `NaN` in the input matrix silently means "never reject" and no R caller
  validates the matrix yet. The transform never emits NA, so the documented pipeline is safe;
  `create_obj_func_gsd()` in P4 is the place for an `anyNA()` guard.
- Kernel index arithmetic is `int`, so `N * m * K` must stay below 2^31; ample at the record's
  largest intended run.
- The P1 commit carries `Co-Authored-By: Claude Opus 5` (the implementing model) rather than the
  session's Fable trailer; the `Claude-Session` line is correct on all three commits.

## Pre-existing issues noticed, not touched

As the brief instructed: `pvals[sample(nsim), ]` in `R/optimisation.R` permutes the first `nsim`
rows rather than sampling from all rows; unary minus in `parse_and_transform()` in
`R/trial_success.R` indexes a second argument that does not exist. Neither is replicated in the
`_gsd` code. Both should become GitHub issues.

Added in P3, both in `R/trial_success.R` and both verified by running `replace_r_indices()`:
the placeholder substitution gives `&&`/`||` a higher precedence than `+` and equal precedence
to each other, so an unquoted `r1 + r2 && r3` (parsed by R as `(r1 + r2) && r3`) is regrouped to
`r1 + (r2 && r3)` and `r1 || r2 && r3` to `(r1 || r2) && r3`, silently; and an exponent-form
literal such as `1e-05` is emitted as the invalid C++ `1e-05.0`. The GSD parser has neither.

## Environment

R 4.6.1 (Windows 11), gsDesign 3.11.0, graphicalMCP 0.3.0 (installed this session from CRAN),
mvtnorm 1.4.2, Rcpp 1.1.2, RcppParallel 6.2.1, testthat 3.3.2, devtools 2.5.2, lintr 3.4.0,
Rtools g++ 14.3.0. In the P0 to P2 session the `.claude/agents/*.md` definitions were not
registered as agent types, so each agent ran as a general-purpose Opus agent with the role text
embedded in its delegation. In the P3 session (2026-09-17) `gsd-reviewer` and `gsd-documenter`
were registered and ran as defined.
