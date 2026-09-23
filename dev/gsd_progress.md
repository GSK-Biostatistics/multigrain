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
