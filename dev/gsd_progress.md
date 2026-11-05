# GSD build progress: phases P0 to P4

Branch `gsd-build`, created 2026-09-16 from `8907be0` (head of `gsd-build-plan`). P0 to P2 were
an orchestrated session: one implementation agent per phase, the orchestrator re-running every
gate from a clean build before accepting it, and one adversarial review agent per phase. P3
(2026-09-17) was implemented and tested by the session itself, reviewed by the `gsd-reviewer`
agent before the gate commit, and documented by the `gsd-documenter` agent afterwards. P4
(2026-09-18) returned to the orchestrated pattern (`gsd-implementer`, then `gsd-reviewer` and
`gsd-documenter`), with one difference: nothing was committed by the session or its agents; the
user reviews the working tree and commits. Nothing below is taken from an agent's report without
the gate having been re-run.

## Status

| Phase | Commit | Gate (re-run before acceptance) | Review verdict |
|---|---|---|---|
| P0 Transform | `e4ec646` | `test-transform_pvalues_gsd.R`: 66 pass, 0 fail, 0 skip | clears gate; two record corrections, four minor follow-ups |
| P1 Simulator | `96dba9f` | `test-sim_pvals_gsd.R`: 29 pass, 0 fail, 0 skip (N = 1e6 block ran) | clears gate; one input-validation gap to decide |
| P2 Kernel | `d0b6556` | `test-graph_shortcut_gsd.R`: 64 pass, 0 fail, 0 skip; `test-RcppExports.R`: 65 pass, 2 pre-existing skips | clears gate; one record inaccuracy, two test-oracle gaps |
| P3 Gain | `0908406` | `test-trial_success_gsd.R`: 116 pass, 0 fail, 0 warn, 0 skip, stable across two runs (`load_all()`; no `src/` change); `test-trial_success.R`: 117 pass, 0 fail | clears gate; three findings, two fixed before the commit, one handed to P4 |
| P4 Optimiser | uncommitted, awaiting user review (2026-09-18) | twelve files from a clean install (re-run after the review fixes), all 0 fail, 0 skip: `test-gsd_guards.R` 11, `test-trial_success_gsd.R` 133, `test-objective_function_gsd.R` 18, `test-calc_power_gsd.R` 37, `test-post_optim_processing_gsd.R` 15, `test-optimisation_gsd.R` 71 (Figure 3b block ran, max gain gap 2.3e-4), plus the untouched `test-optimisation.R` 91, `test-calc_power.R` 53, `test-post_optim_processing.R` 114, `test-trial_success.R` 117, `test-graph_shortcut_gsd.R` 64, `test-transform_pvalues_gsd.R` 66; identity script: `.Random.seed` identical, every component and every field of `trial_success` identical except the compiled `func` | clears gate after one fix; six findings, four fixed, one added, one recorded |

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

P4 (optimiser and post-processing) is built and gated in the working tree but not committed
(section "P4" below). P5 (documentation) has not been started. The branch has not been pushed
and no PR has been opened.

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

## P4 (2026-09-18): optimiser and post-processing

### Before the build: the record was amended

The P4 planning session read the Sol-5-6 review of `trial_success_gsd()`
(`dev/review/trial_success_gsd/ts_gsd_vignette_findings.md`) and put four decisions to the
user; the record was amended `[Rev 2026-09-18]` before the implementer started.

| Sol item | User decision | Landed in |
|---|---|---|
| 5. `graph_optimise()` / `calc_power_pvals()` silently accept a GSD gain and score every rejection as analysis 1 | Guard in the fixed consumers, separate commit; record 4.10 gains its one exception | change set 1 |
| 9. Discount tables: no warning for increases or values outside [0, 1] | Warn, not error; P3 amendment, separate commit | change set 2 |
| 7. Calendar times cannot be supplied | Index scale only: `mean_decision_look`, `time_distribution` labelled by analysis; calendar time is OPEN item 12 | change set 3 |
| 8. `K = NULL` for table-free gains | As recorded; consumers assert `m`, `K` (when set), `!anyNA()`, `max(time) <= K` | change set 3 |
| 2, 3, 4. Fixed-sample `trial_success()` parser bugs | Out of scope; OPEN item 13; issue drafts in the findings file | not in P4 |
| 6. `look_back` documentation | P5 | not in P4 |

The Figure 3b reference was recovered from the paper's code archive
(`origin/sample-size-pub:data-raw/Spiers_gain_function_code_for_outputs.zip`; copies in
`dev/review/fig3b/`, the 21-row table in `tests/testthat/data/gsd_example5_reference.rds`).
Running the paper's own procedure at N = 1e5 showed the record's gate ("optimal w_PFS within
0.02") cannot pass for anyone: the gain is flat near its maximum, the argmax wanders by up to
0.13 across three seeds while the gain at the reference optimum stays within 2.2e-4 of the grid
maximum (Appendix B, "Figure 3b argmax noise"). The gate was restated on the gain (record 6 P4).

### What was built

Three separable change sets, none committed. Files per set, for staging:

- **Change set 1 (guard).** `R/optimisation.R` (one hunk after `check_trial_success()`),
  `R/calc_power.R` (one hunk in `.auto_name_custom_power()`, top level and list elements),
  `NEWS.md` (bug-fix entry of record section 9), new `tests/testthat/test-gsd_guards.R`.
- **Change set 2 (discount-table warning).** `R/trial_success_gsd.R` (`.gsd_gain_warn_tables()`
  and `.gsd_gain_table_bullets()` called from `.gsd_gain_tables()`; roxygen: warning policy
  replaces "by design", the no-consumer `@note` removed, the fixed-optimiser paragraph now says
  it is an error), `man/trial_success_gsd.Rd`, `tests/testthat/test-trial_success_gsd.R` (seven
  tests appended; two pre-existing tests whose tables now legitimately warn are wrapped in
  `suppressWarnings()`).
- **Change set 3 (P4 proper).** New `R/check_gsd.R`, `R/objective_function_gsd.R`,
  `R/optimisation_gsd.R`, `R/calc_power_gsd.R`, `R/post_optim_processing_gsd.R`; new
  `man/graph_optimise_gsd.Rd`, `man/calc_power_pvals_gsd.Rd`; `NAMESPACE` (three exports);
  `man/trial_success.Rd` (was stale on the branch since `d9a3f8e`; regenerated as a side effect);
  new tests `test-objective_function_gsd.R`, `test-optimisation_gsd.R`, `test-calc_power_gsd.R`,
  `test-post_optim_processing_gsd.R`; `tests/testthat/data/gsd_example5_reference.rds`;
  `dev/gsd_identity_check.R`; `dev/review/fig3b/`; `dev/callflows/graph_optimise_gsd.md`;
  `dev/gsd_design_record.md` and `dev/implement_gsd.md` (planning-session and implementer edits,
  not separable). `DESCRIPTION` and `src/` unchanged. Running `devtools::document()` once per
  commit is simpler than splitting `man/` and `NAMESPACE` by hand.

Exported: `graph_optimise_gsd()` (alias `graph_optimize_gsd`), `calc_power_pvals_gsd()`.
Internal: `is_pvals_gsd()`, `check_pvals_gsd()`, `check_trial_success_gsd()`,
`.gsd_check_gain_dims()`, `.gsd_check_alpha()`, `.gsd_kernel_matrix()` (the `dim<-` reshape and
the single `anyNA()` abort), `create_obj_func_gsd()`, `control_prepare_dims()`,
`.sample_pvals_gsd()`, `.graph_optimise_ga_gsd()`, `.graph_optimise_local_gsd()`,
`.gsd_power_by_analysis()`, `.gsd_mean_decision_look()`, `.gsd_time_distribution()`,
`.auto_name_custom_power_gsd()`, `.eval_custom_power_gsd()`, `prune_graph_gsd()`,
`prune_hyp_weights_gsd()`, `prune_edges_gsd()`, `.try_prune_gsd()`. All are copies of their
fixed-sample originals with the kernel call, the `time` matrix, the `K` argument and the `_gsd`
validation as the only differences (plus an `as.integer()` coercion in `control_prepare_dims()`
so that it agrees bit for bit with `control_prepare()`). Reused unchanged: `split_theta()`,
`recover_full_*()`, `param_to_solution()`, `repair_graph()`, `.build_start_matrix()`,
`create_start_params()`, `choose_graph()`, `graph_optimal()`, `.redistribute_mass()`,
`.marginal_violated()`.

`calc_power_pvals_gsd()` returns `local_power`, `local_power_by_analysis` (m x K, cumulative),
`exp_rejections`, `disj_power`, `conj_power`, `mean_decision_look` (NA where never rejected),
`time_distribution` (m x (K+1), columns `never`, `analysis 1`, ...), then the custom entries.
`custom_power` dispatches on class: GSD gain on `time`, fixed-sample gain on `rejected`, plain
function on the time row. `alpha` defaults to `pvals$alpha` and may not exceed it.

### Verified by the orchestrator (clean `clean_dll()` + `install(quick = TRUE)`, `NOT_CRAN=true`)

Test counts are in the status table; every file 0 fail, 0 warn, 0 skip. Figure 3b block
(N = 1e5, seed 20260918, `sfLDOF`, grid step 0.005): gain gap at the paper's optimal w_PFS
against the package's own grid maximum, all 21 cells:

```
     r delta w1_star w1_hat          gap
  0.25  0.50  0.9970  1.000 1.900000e-05
  0.50  0.50  0.9850  0.990 2.666667e-05
  1.00  0.50  0.9555  0.960 2.750000e-05
  2.00  0.50  0.8800  0.920 7.666667e-05
  4.00  0.50  0.7250  0.740 9.000000e-05
  8.00  0.50  0.4685  0.350 4.611111e-05
 16.00  0.50  0.1220  0.145 1.358824e-04
  0.25  0.75  0.9920  0.990 0.000000e+00
  0.50  0.75  0.9805  0.960 3.500000e-05
  1.00  0.75  0.9125  0.920 6.875000e-05
  2.00  0.75  0.8250  0.800 9.083333e-05
  4.00  0.75  0.5860  0.580 9.500000e-06
  8.00  0.75  0.2620  0.185 1.450000e-04
 16.00  0.75  0.0710  0.130 1.391176e-04
  0.25  1.00  0.9850  0.990 2.600000e-05
  0.50  1.00  0.9510  0.960 8.333333e-05
  1.00  1.00  0.8690  0.800 1.250000e-04
  2.00  1.00  0.6725  0.740 2.300000e-04
  4.00  1.00  0.3775  0.345 1.300000e-04
  8.00  1.00  0.1220  0.185 1.177778e-04
 16.00  1.00  0.0275  0.040 1.552941e-04
delta = 1.00: argmax non-increasing in r: TRUE
delta = 0.75: argmax non-increasing in r: TRUE
delta = 0.50: argmax non-increasing in r: TRUE
r = 1, delta = 1.00: optimised w_PFS = 0.8432 (reference 0.8690)
r = 4, delta = 0.75: optimised w_PFS = 0.5804 (reference 0.5860)
r = 8, delta = 0.50: optimised w_PFS = 0.5711 (reference 0.4685)
```

Max gap 2.3e-4 against the 5e-4 tolerance, in the band Appendix B measured for the paper's own
code. The three `graph_optimise_gsd()` runs land within 5e-4 of the grid maximum on the gain;
the third is 0.10 from the reference in w_PFS, which is the flatness the revised gate exists
for. Compiled gains agree with the R formula on `time_distribution` to 3.2e-13.

`dev/gsd_identity_check.R` (installs `324b7ca` from a `git worktree` and the working tree into
two temporary libraries, runs the same seeded `graph_optimise()` in two `Rscript` processes),
orchestrator's run:

```
== graph_optimise() identity check ==
identical(result): FALSE
identical(.Random.seed): TRUE

Falling back to component comparison.
Top-level components that are not identical(): trial_success
  hyp_weight             identical: TRUE   all.equal: TRUE
  trans_matrix           identical: TRUE   all.equal: TRUE
  power                  identical: TRUE   all.equal: TRUE
  solution               identical: TRUE   all.equal: TRUE
  global_output@solution identical: TRUE   all.equal: TRUE
```

The `trial_success` element carries the function compiled by `Rcpp::sourceCpp()`, whose
environment holds an external pointer that differs between processes, so `identical()` on the
whole object cannot hold across processes whatever the code (record 4.10 amended).

Also re-run: `git diff --stat` touches, among pre-existing files, only `NAMESPACE`, `NEWS.md`,
`R/calc_power.R`, `R/optimisation.R` (change set 1), `R/trial_success_gsd.R` (set 2),
`man/trial_success.Rd`, `man/trial_success_gsd.Rd`, `tests/testthat/test-trial_success_gsd.R`
and the two `dev/` documents. The snapshot files rewritten with CRLF by the test run were
restored with `git checkout -- tests/testthat/_snaps`.

### Decisions applied in P4 (need user confirmation)

10. **`graph_optimise_gsd()` refuses a fixed-sample gain** rather than evaluating it on
    `rejected`. `calc_power_pvals_gsd()` accepts one as `custom_power`. Rationale: the optimiser
    has one objective and the `_gsd` gain language covers rejection-only gains; the power
    reporter is where mixing measures is useful.
11. **`alpha` may be below `pvals$alpha` but not above it** (record 4.7). A user who wants a
    larger level must re-transform.
12. **`mean_decision_look` is `NA`**, not `NaN`, for a hypothesis never rejected.
13. **The `max(time) <= K` assertion sits in `calc_power_pvals_gsd()` only**, not in the
    objective closure (record 10 item 6: the P2 kernel cannot emit more).

### Where the implementer diverged (all reported, none violating a LOCKED decision)

- `control_prepare_dims()` coerces `nsim` and `m` with `as.integer()` so its output is
  `identical()` to `control_prepare()` on a matrix (whose `dim()` is integer).
- Two pre-existing tests in `test-trial_success_gsd.R` (`d = c(1/3, 2/3)`, `e = c(2, 1)`) now
  legitimately warn and are wrapped in `suppressWarnings()`.
- cli plural markers with two quantities in one message abort with "Multiple quantities for
  pluralization"; the `m`/`K` mismatch messages spell the numbers out instead. Same trap as the
  literal `<`; added to the brief.
- `.gsd_check_alpha()` runs `rlang::check_number_decimal(max = pvals$alpha)` first, so a value
  above the bound errors with rlang's message (which names the bound); the explicit abort is
  reachable only at `alpha == 0`.
- `dev/gsd_identity_check.R` installs through a child `Rscript` with `R_LIBS` set:
  `devtools::install()` forwards `lib =` to `install_deps()`, not to the install, so a literal
  `lib =` would have installed `main` into the user's library.

### What the record did not anticipate

1. `identical()` on the whole `graph_optimise()` result cannot hold across processes (above);
   the record's fallback was written for `global_output`, which in fact compares identical.
   Record 4.10 and 6 P4 amended.
2. `summary()` of a GSD `multigrain_graph_optimal` does not show `local_power_by_analysis`,
   `mean_decision_look` or `time_distribution`: `summarise_power_object()` in
   `R/graph_optimal.R` prints only the fields it knows. Nothing fails; showing them means
   editing an existing file, a P5 decision. Record 4.7 amended.
3. `man/trial_success.Rd` was stale on the branch: `d9a3f8e` rewrote the roxygen of the
   fixed-sample `trial_success()` without `devtools::document()`. Regenerated as a side effect;
   documentation only.
4. The local `main` ref (`8b19c63`) is behind `origin/main` (`324b7ca`); the identity script
   uses the upstream commit.
5. The `hyp_weight > 1` penalty branch of `create_obj_func()` is unreachable through a normal
   `graph_constraint` (the complement weight goes negative first); the test reaches it with a
   fully fixed `hyp_constraint`, and the `trans_matrix > 1` branch with a fully fixed row.
6. `@param K` of `trial_success_gsd()` still says `K` is "checked against the p-values at first
   use"; Sol item 8 calls this imprecise. One-line P5 fix.

### Not verified

No `R CMD check` (roxygen examples were run by hand with `run.donttest = TRUE`, clean); no
full `test_package()`, by instruction; the `max(time) > K` abort has no test because the kernel
cannot trigger it; parallel closure checked at 2 threads only; the two appendix scripts were not
re-run (P0 gated them; the Figure 3b band reproducing is indirect confirmation).

### P4 review findings and their disposition

The `gsd-reviewer` agent ran against the working tree on the installed build, with its own
fixtures (scripts in the session scratchpad, `s1_items.R` to `s16_claim3.R`). Cleared, each by
running code: record section 8 items 10, 11, 12 and 16 (bitwise agreement with hand
computations from the kernel's `time` matrix; serial, 2- and 4-thread closures `identical()`;
subsampling keeps three dimensions at every `m in {1, 3}`, `K in {1, 2}` and samples from all
rows); brief claims 1 to 5 re-run on the P4 substrate; the gain/p-value contract in all six
cases; the Sol item 5 guard at top level and in named, unnamed and mixed lists, with fixed-sample
gains unaffected; `alpha` above (abort), equal (`identical()` to the default) and below
(`identical()` `local_power` and `time_distribution` to re-transforming the same raw array at the
smaller level, as 4.1 step 4 predicts); `NA` and `NaN` refused before the kernel;
`control_prepare_dims()` `identical()` to `control_prepare()` for default, user-modified and
`verbose = "detail"` controls; pruning removes only what does not lower the gain, respects fixed
weights, edges and a `power_constraint`, and is `identical()` to `prune_graph()` on a K = 1
slice; decision-time semantics under `look_back = c(TRUE, FALSE)` `identical()` to hand tables,
rows of `time_distribution` summing to exactly 1; the discount-table warning in eleven cases;
the Figure 3b block re-run (0 failures, max gap 2.30e-4 reproduced), its R formula shown to be
the paper's `Egain()`, its gain curve shown to vary by 0.007 to 0.049 over the grid (so the gate
cannot pass by constancy), and the compiled gain checked against the R formula at the reviewer's
own grid points (relative differences below 1.0e-13); a line-by-line comparison of the five new
files against their originals with no undisclosed divergence; and an end-to-end check of the
`dim<-` reshape on a hand-built 4 x 2 x 2 array with known decision times. The reviewer judged
`dev/gsd_identity_check.R` to do what it claims (two libraries, two processes, the working tree
rather than HEAD) and, having shown independently that two `trial_success()` objects built in
separate processes differ only in `func`, accepted the "every component except `trial_success`"
fallback as a correct reading of the 4.10 gate.

Findings:

1. **Does not clear the gate.** `calc_power_pvals_gsd()` checked `m` only for GSD gains, so a
   fixed-sample `trial_success()` gain with the wrong `m` in `custom_power` was handed a
   narrower `rejected` matrix and its compiled function read out of bounds (reproduced:
   `trial_success(r1 + r2 + r3)` on a 2-hypothesis object returned -15239942 with warnings and
   no abort). The same hole is pre-existing in `calc_power_pvals()`, but P4 documents that it
   accepts fixed-sample gains, so the "never feed the compiled function a matrix it does not
   fit" contract of the brief applies. **Fixed**: the `m` check now runs for every
   `multigrain_trial_success` entry; test added.
2. The `[0, 1]` range check on discount tables had no tolerance, so `d = c(1, 1 + 1e-12)`
   warned with a message naming a value that prints as 1. **Fixed**: same
   `sqrt(.Machine$double.eps)` tolerance as the monotonicity check; roxygen and tests updated.
3. Roxygen said one minus the `never` column "is `local_power`"; the difference is 1.1e-16.
   **Fixed**: reworded to "up to floating-point rounding".
4. No range check on the transformed array: a planted -1 was silently a rejection and a planted
   2 silently "never". Not required by the record. **Added** to `.gsd_kernel_matrix()` next to
   the `anyNA()` assertion (decision 14 below).
5. The Figure 3b gate has low discriminating power by construction (59 to 114 of 201 grid points
   inside the tolerance at N = 3e4): it shows the package reproduces the paper's gain surface
   near its flat maximum, not the optimal weights. **Recorded** in record 6 P4.
6. The identity script exempted the whole `trial_success` component where only `func` can
   differ. **Tightened**: `m`, `objective`, `cpp_code` and `body(func)` are now compared too.

Also from the review: the reviewer's own test runs rewrote the three CRLF snapshot files to LF
(restored by the orchestrator with `git checkout -- tests/testthat/_snaps`); probing the
deserialised compiled closure's environment with `identical()` segfaulted R (pre-existing Rcpp
behaviour, not P4's).

Re-gate after the fixes (orchestrator, second clean `clean_dll()` + `install(quick = TRUE)`,
`NOT_CRAN=true`): all twelve files 0 fail, 0 warn, 0 skip, with the counts in the status table
(`test-trial_success_gsd.R` 133, `test-objective_function_gsd.R` 18, `test-calc_power_gsd.R` 37
after the added tests); the Figure 3b block ran again with the same 21 gaps (max 2.30e-4); the
identity script, re-run by the implementer after tightening, reports `trial_success$m`,
`$objective`, `$cpp_code` and `body(func)` all `identical()` between `main` and the branch, with
only `$func` (the external pointer) differing. The reviewer's reproduction of finding 1 now
aborts naming `custom_power$total`.

14. **Values outside [0, 1] in the transformed array abort** in `.gsd_kernel_matrix()` (review
    finding 4; record 10 item 11 asked for `NA` only).

### Call-flow document

`dev/callflows/graph_optimise_gsd.md` (1283 lines), written by `gsd-documenter` against the
working tree in the format of the P0 and P3 documents: notation, the call chain at a glance,
one section per function for all 22 functions in the five new files, the reused fixed-sample
helpers described at their point of call, and closing sections on the guard hunks and the
discount-table warning. The documenter verified by running code on the installed build that
`.sample_pvals_gsd()` keeps three dimensions, `.gsd_kernel_matrix()` gives N by m*K,
`calc_power_pvals_gsd()` returns every documented field with the stated shapes and dimnames,
`.gsd_check_alpha()` and `check_trial_success_gsd()` behave as documented, and
`.eval_custom_power_gsd()` dispatches by class. It reports no roxygen-code disagreement and
confirms that the GSD copies differ from their fixed-sample originals only in the kernel call,
the matrix handed to the gain, the `K` argument and the `as.integer()` coercion in
`control_prepare_dims()`. Two observations: `calc_power_pvals_gsd()` always uses the serial
kernel (the reporting path is called once; deliberate), and `local_power_by_analysis[, K]`
carries dimnames where `local_power` does not (cosmetic).

## Environment

R 4.6.1 (Windows 11), gsDesign 3.11.0, graphicalMCP 0.3.0 (installed this session from CRAN),
mvtnorm 1.4.2, Rcpp 1.1.2, RcppParallel 6.2.1, testthat 3.3.2, devtools 2.5.2, lintr 3.4.0,
Rtools g++ 14.3.0. In the P0 to P2 session the `.claude/agents/*.md` definitions were not
registered as agent types, so each agent ran as a general-purpose Opus agent with the role text
embedded in its delegation. In the P3 session (2026-09-17) `gsd-reviewer` and `gsd-documenter`
were registered and ran as defined.

## Stopping rules (2026-10-13)

Branch `gsd-build-stopping`, from `gsd-build`. Efficacy stopping rules and the stopping-look
gain symbol `s` (user-testing suggestions 1 and 2, `dev/review/gsd_user_testing/README.md`
section 8; design record section 4.11). One implementation agent per phase; the orchestrator
re-ran every gate from `devtools::clean_dll(); devtools::install(quick = TRUE)` with
`NOT_CRAN=true` before committing. Phase order Q0, Q1, Q3, Q2, Q4, Q5. Counts below are
0 fail, 0 warn, 0 skip unless noted.

### Status

| Phase | Commit | Gate (re-run before acceptance) | Review verdict |
|---|---|---|---|
| Prerequisite | (none) | `test-graph_shortcut_gsd.R` 64 on the freshly installed toolchain, before any code | n/a |
| Q0 Kernel | `9ba7bfd` | `test-graph_shortcut_gsd.R` 135; `test-RcppExports.R` 65 (2 pre-existing gMCPLite skips); unchanged: `test-calc_power_gsd.R` 37, `test-objective_function_gsd.R` 18, `test-optimisation_gsd.R` 71 | orchestrator read the diff and re-ran the gate from a clean build; no separate review agent this session |
| Q1 Rule object | `69fe6f4` | `test-stopping_rule_gsd.R` 107; unchanged: `test-trial_success_gsd.R` 133 | orchestrator read the diff and re-ran the gate from a clean build; no separate review agent this session |
| Q3 Gain symbol `s` | `8d8f9d2` | `test-trial_success_gsd.R` 178; `test-trial_success.R` 117; consumers unchanged | orchestrator read the diff and re-ran the gate from a clean build; no separate review agent this session |
| Q2 Consumers | `e5dd92f` | `test-check_gsd.R` 56, `test-calc_power_gsd.R` 65, `test-objective_function_gsd.R` 21, `test-post_optim_processing_gsd.R` 21, `test-optimisation_gsd.R` 79 (Figure 3b block ran), `test-gsd_guards.R` 11, `test-graph_optimal.R` 27, `test-calc_power.R` 53, `test-optimisation.R` 91, `test-stopping_rule_gsd.R` 107, `test-trial_success_gsd.R` 178, `test-graph_shortcut_gsd.R` 135, `test-RcppExports.R` 65 (2 skips); `dev/gsd_identity_check.R`: every component identical except `trial_success$func` | orchestrator read the diff and re-ran the gate from a clean build; no separate review agent this session |
| Q4 Display | `f3b3ed7` | `test-graph_optimal.R` 33, `test-optimisation_gsd.R` 79, `test-optimisation.R` 91, `test-plot_graph_optimal.R` 44; existing `graph_optimal_example` snapshots byte-identical | orchestrator read the diff and re-ran the gate from a clean build; no separate review agent this session |
| Q5 Documentation | last commit on the branch | `devtools::document()` clean, a second run a no-op; `pkgdown::check_pkgdown()` "No problems found" (run on a scratch copy with the `.qmd` articles renamed `.Rmd`: the `quarto` R package is not installed, see below); no code changed | orchestrator read the diff and re-ran the gate from a clean build; no separate review agent this session |

### What was built

- **Q0, kernel** (`src/graph_shortcut_gsd.cpp`, `R/RcppExports.R`, `src/RcppExports.cpp`,
  `R/check_gsd.R`). `graph_shortcut_gsd()` and `graph_shortcut_gsd_parallel()` take a trailing
  `Nullable<LogicalMatrix> stop_rule = R_NilValue`, validated by a shared `resolve_stop_rule()`
  (`NULL` or zero columns: no rule; otherwise `2^m` rows, `K - 1` columns, and `m <= 30` so
  `1 << m` cannot overflow), and return `stop_look` beside `rejected` and `time`. After the
  cascade at an interim the trial stops if the table entry for the cumulative rejection
  bitmask is `TRUE`; a trial that rejects every hypothesis stops at that look with or without a
  rule; otherwise `stop_look` is K. `.gsd_stop_table(conditions, m, K)` in `R/check_gsd.R`
  evaluates the conditions in R on every rejection set (row `mask + 1`, hypothesis i at bit
  i - 1), aborts on a non-logical, non-scalar or `NA` result naming the interim and the set,
  caps `m` at 16 (about 2.5 s to build at `m = 16`) and returns a `2^m` by 0 matrix at `K = 1`.
  Tests: no-rule output identical to the committed kernel's saved output
  (`tests/testthat/data/gsd_kernel_reference.rds`, 18 KB); K = 1 identity with
  `graph_shortcut()`; with a rule, the time matrix identical to the decision-time edit of the
  user-testing prototype and to `gsd_ref_test()`, a reference that really stops trials, on
  every cell; parallel identical to serial. Verbatim copies of the user-testing reference are
  in `tests/testthat/helper-gsd_stopping.R`.
- **Q1, rule object** (`R/stopping_rule_gsd.R`, `man/stopping_rule_gsd.Rd`, `NAMESPACE`,
  `tests/testthat/test-stopping_rule_gsd.R`, `tests/testthat/_snaps/stopping_rule_gsd.md`).
  Exported `stopping_rule_gsd()`, `print()` and `summary()` methods; internal
  `is_stopping_rule_gsd()`, `check_stopping_rule_gsd()`, `.gsd_check_stopping_rule_dims()` and
  the capture, parsing and validation helpers (call flow:
  `dev/callflows/stopping_rule_gsd.md`). Accepts an expression, a string, `list(...)`, a
  literal `c("...")` of strings, or an injected list or character vector. `document()` also
  repaired a missing brace in `man/trial_success.Rd` and moved `Config/roxygen2/version` in
  `DESCRIPTION` from 8.1.0 to 8.1.1.
- **Q3, gain symbol `s`** (`R/trial_success_gsd.R`, `man/trial_success_gsd.Rd`,
  `tests/testthat/test-trial_success_gsd.R`, `tests/testthat/_snaps/trial_success_gsd.md`).
  `s` in comparisons, arithmetic and `d(s)`; reserved as a table name; the compiled signature
  is `double powerFunc(IntegerMatrix t, IntegerVector s = IntegerVector::create())` (Rcpp's
  wrapper shows the default as `as.integer(c())`), with a length check on `s` only when the
  gain uses it; `uses_s` on the object and in `print()`/`summary()`; new `.gsd_gain_ast()` so
  `uses_s` is read from the tree that is translated. The signature line changed in every
  `cpp_code` snapshot.
- **Q2, consumers** (`R/optimisation_gsd.R`, `R/objective_function_gsd.R`,
  `R/post_optim_processing_gsd.R`, `R/calc_power_gsd.R`, `R/check_gsd.R`, `R/graph_optimal.R`,
  their `man/` pages and tests, new `tests/testthat/test-check_gsd.R`). `stopping_rule = NULL`
  on `graph_optimise_gsd()` and `calc_power_pvals_gsd()`; the rule is settled and tabulated
  once on the public path and the table passed to the closure, the drivers, the pruning and the
  reporter; the gain is called with `(time, stop_look)` everywhere; a plain `custom_power`
  function receives `(time_row, stop_look)`; new fields `stop_distribution` and
  `mean_stop_look`; the result carries `stopping_rule` when supplied. `.gsd_check_graph_dims()`
  checks the graph's `m` against the p-values in R first (review fix). The reporter is split
  into the public `calc_power_pvals_gsd()` and the internal `.calc_power_pvals_gsd_core()`.
  `test-check_gsd.R` tests every helper in `R/check_gsd.R` directly.
- **Q4, display** (`R/graph_optimal.R`, `tests/testthat/test-graph_optimal.R`,
  `tests/testthat/_snaps/graph_optimal.md`). `summarise_power_object()` prints power by
  analysis, the decision-time distribution, the stopping-look distribution and the mean
  stopping look when present; `print()` lists the rule's conditions and `summary()` summarises
  it after the gain. Bug fixed: `power_object$local_power` partially matched
  `local_power_by_analysis` on a list without `local_power`; now `[[`.
- **Q5, documentation.** `R/stopping_rule_gsd.R` (one roxygen sentence on the 16-hypothesis
  cap) and `man/stopping_rule_gsd.Rd`; `_pkgdown.yml` (a "Group sequential designs" reference
  section; `trial_success_gsd` moved out of "Inputs"); `vignettes/articles/trial_success.qmd`
  ("Stopping for efficacy", a pointer to `s` under "Decision times", and a warning callout on a
  rule without a value for stopping); `NEWS.md` (two new-functionality bullets, one bug fix);
  `dev/gsd_design_record.md` (new 4.11, amendments to 4.3, 4.6, 4.7 and 4.10, section 10 item
  14); `dev/review/gsd_user_testing/README.md` (section 8 items 1 and 2 marked closed);
  `dev/callflows/graph_optimise_gsd.md` and `dev/callflows/trial_success_gsd.md` (updated),
  `dev/callflows/stopping_rule_gsd.md` (new); this entry.

### Decisions applied (user decisions, 2026-10-13)

- **A. `s` is a second argument of the compiled gain**, not an appended column of the time
  matrix. The time matrix stays N by m, every reported field and every gain that does not use
  `s` is unchanged, and no N by (m + 1) matrix is built per evaluation. Consequence: every
  `cpp_code` snapshot changed (signature line), and every call site passes two arguments.
- **B. A trial that rejects every hypothesis stops at that analysis, with or without a rule**;
  otherwise S = K. Consequence: without a rule `stop_distribution` is not all on analysis K
  and `d(s)` is not the constant `d(K)` in those trials (3.3% of trials in the article's
  example).
- **C. A plain `custom_power` function always receives `(time_row, stop_look)`.** A
  one-argument function written for the P4 contract now errors; there is no arity sniffing.
- **D. Toolchain installed in-session, no `air`.** Formatting by hand; see "Environment"
  below.

### Where the implementation diverged from the brief

- Phase order Q0, Q1, Q3, Q2, Q4, Q5: Q3 (the gain signature) before Q2 (the consumers), so the
  gain call sites were rewired once.
- The reporter is split into a validating `calc_power_pvals_gsd()` and an internal
  `.calc_power_pvals_gsd_core()`, so the stopping table is built once per public call;
  pruning and the drivers no longer re-run `is_graph_valid()` through the public function.
  `graph_optimise_gsd()` still calls the public function once for the final graph, so the
  table is built twice per optimisation.
- `.gsd_check_graph_dims()` uses an explicit `cli::qty()` per number so that `m = 1` reads
  "hypothesis".
- The Q2 tests needed shapes the brief did not give: test (c) uses
  `trial_success(r1 + r2 + r3)` (its `m` must match the p-values) and test (d)
  `d(s) + 0 * (r1 + r2 + r3)` (a gain must name at least one `r<i>` or `t<i>`).
- Q0: the `r1 && (r2 || r3)` fixture never fires at look 1 (13 trials reject H1 there), so the
  test also asserts at look 2 (475 trials). The kernel refuses `m > 30`. Reference outputs of
  the committed kernel are saved in `tests/testthat/data/gsd_kernel_reference.rds` and the
  user-testing reference is copied verbatim into `tests/testthat/helper-gsd_stopping.R`.
- Q1: the `r<i>` pattern is `^r[1-9][0-9]*$` (`r0`, `r01` refused); `list(r1)` has
  `n_interim = 1` and is not recycled; `TRUE` alone is refused; `K = 1` returns an empty
  condition list.
- Q3: a gain using only `s` is refused because at least one `r<i>`/`t<i>` is still required
  (spec S4); `d(t1) + d(t2) + d(t3) + 0.1 * (s == 1)` is the supported shape.
- Q5: a roxygen sentence on the 16-hypothesis cap was added to `stopping_rule_gsd()` (the cap
  was enforced but not documented); a one-paragraph `[Rev 2026-10-13]` note was added to
  record section 4.10 (a second exception to "no existing file is edited":
  `R/graph_optimal.R`), besides the amendments the brief named.
- No `air`: formatting by hand; `lintr` clean on every touched file apart from pre-existing
  `backport_linter` warnings on `%||%`.

### What was not verified

No `R CMD check`; no `test_package()`; the vignette was not rendered (the new chunks are
`eval: false`, and their printed output was produced by running the code with
`devtools::load_all()` on this branch); the kernel's NA-as-TRUE path for a hand-built table has
no test (unreachable through the public API, since `.gsd_stop_table()` never emits `NA`); only
serial, 2- and 4-thread parallel closures were tested; `pkgdown::check_pkgdown()` could not
run on the repository itself (it needs the `quarto` R package to read the `.qmd` articles)
and was run on a copy with the articles renamed; `dev/gsd_identity_check.R` was last run after
Q2, not after Q4.

### What the brief did not anticipate

1. An NA in a hand-built stopping table reads as `TRUE` in the kernel; `.gsd_stop_table()`
   never emits one. Errors raised inside a condition are not wrapped.
2. The `quarto` R package is needed for `pkgdown::check_pkgdown()` on this repository.
3. `NEWS.md` has no entry for the group sequential family itself (the wording proposed in
   record section 9 was never added), so the new bullets refer to functions the development
   NEWS does not introduce. Left for the user; out of the Q5 brief.
4. `mean_decision_look` is still not printed by `summary()` (record section 10 item 14).
5. The bare-table-name error of `trial_success_gsd()` still says a table "must be applied to a
   decision time, e.g. `d(t1)`" without mentioning `d(s)`; R code, not changed in Q5.
6. `devtools::document()` must run under `LC_ALL=C.UTF-8`: under the POSIX locale it corrupts
   `man/multigrain-package.Rd`.

### Environment

Ubuntu 24.04, R 4.3.3; packages installed from CRAN source this session (Posit Package Manager
binaries were blocked by the proxy; `libuv1-dev` and `libtiff-dev` added with apt); roxygen2
8.1.1; `cyclocomp` installed so `lintr` runs fully; `air` unavailable; pkgdown installed,
`quarto` (R package) not installed.

### Withdrawal of `s` (2026-10-13)

**Decision (user, the same day).** The gain symbol `s` added in Q3 (`8d8f9d2`) is withdrawn:
valuing the stop inside the gain was judged too complex. Everything else on the branch stays:
the rule object, the kernel's `stop_look`, `stop_distribution`, `mean_stop_look`, the display
and the graph-dimension check. Phase R1; design record 4.11 `[Rev 2026-10-13b]` and section 10
item 15. Decisions A and C above are superseded.

- The compiled gain is `double powerFunc(IntegerMatrix t)` again, called as `func(time)` by the
  objective closure, both optimiser drivers and the reporter; no `uses_s`, no length check.
- The grammar is that of 4.6: a bare `s` is an un-injected symbol (the generic error with the
  `!!` hint) and a table may be named `s`. A plain `custom_power` function is called as
  `f(time_row)`.
- The stopping rule still refuses `s`, now through the generic unknown-symbol message.
- `R/trial_success_gsd.R`, `man/trial_success_gsd.Rd`, `tests/testthat/test-trial_success_gsd.R`
  (apart from one new test) and `tests/testthat/_snaps/trial_success_gsd.md` are as at
  `69fe6f4`; `dev/callflows/trial_success_gsd.md` is as at `69fe6f4` except the comments of the
  example C++, which now match the generated code.

Files touched: `R/trial_success_gsd.R`, `R/calc_power_gsd.R`, `R/objective_function_gsd.R`,
`R/optimisation_gsd.R`, `R/stopping_rule_gsd.R`; `man/trial_success_gsd.Rd`,
`man/calc_power_pvals_gsd.Rd`, `man/stopping_rule_gsd.Rd`; `tests/testthat/test-trial_success_gsd.R`
(the `s` tests removed; one new test that `s` is an un-injected symbol and a valid table name),
`tests/testthat/_snaps/trial_success_gsd.md`, `tests/testthat/test-calc_power_gsd.R` (plain
function back to one argument; the `d(s)` test replaced by a discounted gain checked by hand on
the stopped and unstopped time matrices), `tests/testthat/test-objective_function_gsd.R`,
`tests/testthat/test-optimisation_gsd.R`, `tests/testthat/test-post_optim_processing_gsd.R`,
`tests/testthat/test-stopping_rule_gsd.R`; `NEWS.md`; `vignettes/articles/trial_success.qmd`
(the `s` pointer and the `s` paragraph removed, the gain of the end-to-end example is
`d(t1) + d(t2) + d(t3)`, the warning callout replaced by the note "A stopping rule changes which
graph is optimal"); `dev/gsd_design_record.md`; `dev/callflows/trial_success_gsd.md`,
`dev/callflows/graph_optimise_gsd.md`, `dev/callflows/stopping_rule_gsd.md`; this entry.

Gate (orchestrator, clean build, NOT_CRAN=true): `test-trial_success_gsd.R` 138, `test-trial_success.R` 117, `test-calc_power_gsd.R` 65, `test-objective_function_gsd.R` 21, `test-post_optim_processing_gsd.R` 21, `test-optimisation_gsd.R` 79, `test-stopping_rule_gsd.R` 107, `test-check_gsd.R` 56, `test-graph_optimal.R` 33, `test-gsd_guards.R` 11, `test-graph_shortcut_gsd.R` 135, all 0 fail, 0 skip; `R/trial_success_gsd.R`, `man/trial_success_gsd.Rd` and the gain snapshots identical to `gsd-build`; `devtools::document()` a no-op afterwards. The prompt for the follow-on spending-function job is `dev/implement_hsd_opt.md`.

## Spending-parameter optimisation (2026-10-14)

Branch `gsd-build-spending`, from `gsd-build` at `3e55fb5`. Joint optimisation of
Hwang-Shih-DeCani spending parameters with the graph (`dev/implement_hsd_opt.md`; design
`dev/hsd_opt_design_revised.md`; build plan `dev/hsd_opt_build_plan.md`; design record section
4.12). Phase D ran straight through by user decision and fixed its defaults by the plan's own
rules (`dev/review/hsd_opt_discovery.md`, "Decisions"); one Opus implementation agent per phase,
reporting rather than committing; the orchestrator re-ran every gate from `devtools::clean_dll()`
and a clean `devtools::install(quick = TRUE)` (or `load_all()` while a measurement battery held
the installed library) with `NOT_CRAN=true` before committing. Phase order D-a, P1, D-b, P2, P3,
P4, P5. Counts below are 0 fail, 0 warn, 0 skip.

### Status

| Phase | Commit | Gate (re-run before acceptance) | Review verdict |
|---|---|---|---|
| Prerequisite | `2cc39d2` | baselines `test-transform_pvalues_gsd.R` 66, `test-optimisation_gsd.R` 79, `test-graph_optimal.R` 33, `test-check_gsd.R` 56, `test-objective_function_gsd.R` 21; `dev/gsd_identity_check.R` green | n/a |
| D-a (d01, d03) | `c3cf03f` | scripts run end to end; no node of `seq(-10, 2, 0.25)` on four schedules or either corner fails well-ordering or invertibility; column copy 8 % of a kernel call | orchestrator read the scripts and outputs |
| P1 spec class | `2e37fce` | `test-transform_pvalues_gsd_opt.R` 168; baselines unchanged; lintr; `document()` no-op; identity green | orchestrator read the diff and re-ran the gate from a clean recompile |
| D-b (d02, d04, d05, report) | `849ba71` | scripts run end to end; report complete; warm-start gate PASS on M0, M1 and the section 2 design; no pause condition | orchestrator read the Decisions and the gap tables |
| P2 objective and bank | `6f9e7e9` | `test-objective_function_gsd_opt.R` 812 (node oracle `identical()` at every node and node pair, with a stopping rule, with two threads); baselines unchanged; no `<<-`; lintr; identity green | orchestrator read the closure, verified the copy claim with `tracemem()` itself, re-ran the gate |
| P3 driver | `74f0425` | `test-optimisation_gsd_opt.R` 115; `test-optimisation_gsd.R` 79 (Figure 3b ran), `test-graph_optimal.R` 33 with its snapshot byte-identical, `test-choose_graph.R` 11, `test-optimisation.R` 91, the P1/P2 files; lintr; `document()` no-op; identity green; `gate_p3_benchmarks.R` at nsim 1e5: six checks PASS | orchestrator read the diff, fixed two things (below), re-ran the gate from a clean install |
| P4 display | `db1ce3f` | `test-graph_optimal_hsd_opt.R` 20 in its own snapshot file; `test-graph_optimal.md` byte-identical; `test-plot_graph_optimal.R` 44; lintr; identity green | orchestrator read the diff and re-ran the gate |
| P5 documentation | `996824b` | every GSD and fixed-sample test file at its count; lintr; `document()` no-op; `pkgdown::check_pkgdown()` "No problems found" on a scratch copy with the `.qmd` articles renamed `.Rmd`; identity green | orchestrator read NEWS, 4.12, the call flow and the lay note |
| Verification | `594cddd`, `4680d9f` | `dev/review/hsd_opt/gsd_fixed_path_identity.R` PASS (untagged transform `identical()`; result identical in every component but the compiled pointer; no new field on a fixed-path result); `dev/review/hsd_opt/acceptance_m3.R` PASS (section 2 design, shared parameter, stopping rule, nsim 1e4: spec 19 s, optimiser 42 s, power `identical()` to `calc_power_pvals_gsd()` on the resolved object) | orchestrator |

### What was built

- **P1, spec class** (`R/transform_pvalues_gsd_opt.R`, touchpoints in `R/transform_pvalues_gsd.R`
  and `R/check_gsd.R`, `tests/testthat/helper-hsd_opt.R`, `test-transform_pvalues_gsd_opt.R`
  and its snapshot). `sf_hsd_opt()` (exported) tags a hypothesis; `transform_pvalues_gsd()`
  hands a tagged `spending` to `.transform_pvalues_gsd_spec()`, which transforms fixed and
  short-circuited hypotheses as today, builds one boundary table per (schedule, lattice node)
  with the well-ordering check and a new invertibility check, and returns a
  `multigrain_pvals_gsd_spec` (raw draws under `raw`, never `pvals`; no inheritance).
  `resolve_pvalues_gsd()` (exported) inverts the stored tables into a `multigrain_pvals_gsd`
  `identical()` to the fixed path's. `.gsd_spend()` aborts if a tagged function ever reaches
  the fixed transform; `check_pvals_gsd()` refuses a spec with a hint.
- **P2, objective and bank** (`R/objective_function_gsd_opt.R`, `test-objective_function_gsd_opt.R`
  and its snapshot). `create_obj_func_gsd_opt()`: parameter genes appended to the graph genes,
  clamped and mapped to nodes; the decode-and-penalise block of `create_obj_func_gsd()`
  duplicated verbatim; per-node slices copied into the working kernel matrix only for groups
  whose node moved; an archive of the best graph per node. All mutable state in one environment
  created inside the factory; the matrix is updated by take-out/put-back.
- **P3, driver** (`R/optimisation_gsd_opt.R`, dispatch in `R/optimisation_gsd.R`, `spending =`
  on `graph_optimal()`, a progress bar in `.hsd_node_tables()`, `test-optimisation_gsd_opt.R`,
  `dev/review/hsd_opt/gate_p3_benchmarks.R` with its output). Stage G: `GA::ga()` with the fixed
  driver's options on the bank, popSize `min(max(40 (m + d), 200), 500)` unless user-set. Stage
  L: `idx_local` drawn once and pinned; COBYLA on the graph at the GA's node; a warm-started
  profile along each axis from three candidate starts, with the paired batch-means difference
  and standard error; rounds for d >= 2; stale axes recomputed once against the final optimum.
  `choose_graph()`, `prune_graph_gsd()` and `calc_power_pvals_gsd()` unchanged on the resolved
  object.
- **P4, display** (`R/graph_optimal.R` two guarded calls, display helpers appended to
  `R/optimisation_gsd_opt.R`, `test-graph_optimal_hsd_opt.R` and its own snapshot file).
- **P5, documentation** (roxygen in four `R/` files and the regenerated `man/` pages,
  `_pkgdown.yml`, `NEWS.md`, `vignettes/articles/trial_success.qmd` "Optimising the spending
  parameter", `dev/gsd_design_record.md` 4.12 / 4.10 / 4.5 / section 10 item 16,
  `dev/callflows/graph_optimise_gsd_opt.md` new and `graph_optimise_gsd.md` amended,
  `dev/hsd_opt_explained.md`, `dev/review/hsd_opt/article_example.R` with its output).
- **Phase D** (`dev/review/hsd_opt/d01`-`d05`, `d_helpers.R`, `dev/review/hsd_opt_discovery.md`):
  table cost 0.71-0.88 ms per `gsBound1()` call at K = 2; eager tables 16-19 s per default
  spec; copy 8 % of a kernel call; warm-start gaps 2.1e-5 / 2.3e-5 / 3.7e-4; the parameter
  gene costs 0.90-0.98 of the generations; lattice loss 2.1e-5 at step 0.25; E_C 23-128.

### Decisions applied (user decisions, 2026-10-14)

- **A. Phase D runs straight through**, resolving its defaults by the plan's rules: step 0.25;
  full 21-node profile; popSize constant 40; Stage L as designed; a verbosity-gated progress
  indicator on table building; the same `local_opt` at profile nodes. Consequence: no pause
  between discovery and P1-P5; the numbers are in the discovery report and record 4.12.
- **B. Implementers on Opus, one per phase, reporting not committing**; the orchestrator gates
  and commits. Consequence: every phase has a brief under `dev/review/hsd_opt/briefs/`.
- **C. No `<<-` anywhere.** Mutable closure state lives in one explicit environment created
  inside the factory and is updated with ordinary `<-`. Consequence: P2 measured that
  `state$mat[, cols] <- v` inside a function copies the matrix (the d03 probe was at top
  level), so the matrix is updated by take-out/put-back; a test greps every `_opt` file for
  `<<-`.

### Where the implementation diverged from the design

The eleven reconciliations of `dev/hsd_opt_build_plan.md` section 1 (the two-argument wrapper
per node and the invertibility check; genes appended and the decode block duplicated; `raw` not
`pvals`; dispatch before `check_pvals_gsd()`; popSize adjusted in the driver, not in
`control_prepare_dims()`; `idx_local` pinned and step 1 as `.hsd_refine_node()`; labels;
display insertion points; warm starts never from `.build_start_matrix()[1, ]`; `spending =` on
`graph_optimal()`; `choose_graph()` unchanged), plus: no `allow_opt` flag (an internal abort in
`.gsd_spend()` instead); `spending_recycled` passed from `transform_pvalues_gsd()`; a progress
bar after two seconds instead of a progress step (a step at default verbosity leaked into
snapshots); an axis whose own sweep improved the optimum is recomputed once against the final
optimum with no refinement at the optimum's own node, so S, the profile and the reported gain
agree exactly at theta* (the orchestrator's fix after P4 found a 5e-4 gap on the fixture);
`spending_step` and `spending_group` attached beyond the design's seven fields (the
orchestrator's fix so the display need not infer the hypothesis-to-parameter map).

### What was not verified

- The article does not render (the `quarto` R package and binary are absent); its output is
  pasted from `dev/review/hsd_opt/article_example.R`.
- `num_threads >= 2` on the new path is exercised by the P2 oracle, the acceptance run and the
  benchmarks, not by the P3 tests. `density = 2` (every second node) is kept in the code but
  not exposed or tested. The bank memory guard is exercised with synthetic byte counts only.
- d >= 3 and tagged designs with more than two active hypotheses were not run.
- The benchmark's section 2 run took 8.0 min against the design's 3.5-6.5 min projection
  (within the 2x pause threshold); the article's M1 run took 77 s plus 20 s for the spec.

### What the brief did not anticipate

1. The "fewer than two distinct positive boundaries" abort lives in `.gsd_invert()`, not in
   the table builder, so the spec builder needs its own invertibility check.
2. A direct assignment through an environment copies a matrix inside a function; only the
   take-out/put-back form (or `<<-`, now forbidden) does not.
3. The popSize rule `40 (m + d)` is inert on every benchmark design (m + d <= 5 gives 200);
   it bites first at m = 4, d = 2 (240), which the tests assert.
4. The fixed local driver re-subsamples internally, so Stage L cannot call it on pinned rows;
   the oracle compares against a by-hand `nloptr` call and against the driver under the same
   seed.
5. Optimiser noise on the section 2 design is about 4e-4, close to the 5e-4 plateau tolerance;
   and the default `bounds = c(-4, 1)` put that design's optimum on the lower bound (the gain
   is still rising below -4). `print()` does not warn when the optimum sits on a bound.
6. The installed library went stale between phases while agents used `load_all()`: the first
   acceptance run showed no display block until the tree was reinstalled. The final gate
   reinstalls before every check.
7. The user-testing 7.4 script uses nsim 2e5; Phase D used 1e5 and re-derived the M0/M1
   references on its own draws. The 7.5 effect reproduces qualitatively (gamma moves from -4
   to -0.25, gain 0.6802 to 0.6930; with H2 fixed: sfLDOF 0.6709, sfLDPocock 0.6905 on the
   article's draws).

### Environment

R 4.3.3, gsDesign 3.11.0, GA 3.2.5, nloptr 2.2.1, testthat 3.3.2, roxygen2 8.1.1, lintr 3.4.0,
pkgdown 2.0.7; no `quarto`, no `air`; 4 cores, 15 GB. `document()` under `LC_ALL=C.UTF-8`.
The session's worker restarted twice before the build began; every phase was committed as soon
as its gate passed and every brief and report lives on the branch.
