# Opus prompt — implementation pass: group sequential extension of `multigrain`

> Fresh Claude Code session in a clone of `GSK-Biostatistics/multigrain`, on a new branch
> `gsd-build` off `main`. The design record `dev/gsd_design_record.md` was produced by a
> separate agent in an earlier session and is the specification for this work.
>
> Before starting the session: copy `debug/gsd_design_record.md` to `dev/gsd_design_record.md`
> and commit it. `debug/` is git-ignored and a fresh clone cannot see it; `dev/` is tracked and
> build-ignored.
>
> **Status 2026-09-16.** P0 (`e4ec646`), P1 (`96dba9f`) and P2 (`d0b6556`) are done on
> `gsd-build`, each gate re-run by the orchestrator and each phase adversarially reviewed; see
> `dev/gsd_progress.md`. Remaining: P3, P4, P5. Lessons from those phases are folded in below
> and marked `[Rev 2026-09-16]`.
>
> **Status 2026-09-17.** P3 (`0908406`) is done; see `dev/gsd_progress.md`. Remaining: P4,
> P5. Lessons marked `[Rev 2026-09-17]`.

---

## Role

You are the implementation agent. You are executing a plan someone else wrote. Your job is to
implement it faithfully, not to improve it.

## Read first

- `dev/gsd_design_record.md` — the specification. Read all of it, including the appendix.
  Run the two appendix scripts before writing any code and confirm you get the numbers quoted
  in the record. If you do not, stop and report; do not proceed on a machine where the
  reference numbers do not reproduce.
- `NEWS.md`, for house conventions on naming, argument style and news entries. (There is no
  `CLAUDE.md` in this repository `[Rev 2026-09-16]`. House style as practised: 4-space
  indentation, rlang standalone checks from `R/check_types.R`, `rlang::check_dots_empty()`,
  optional arguments after `...` and named, `cli::cli_abort()` for errors.)
- The P0 to P2 files, which later phases build on: `R/transform_pvalues_gsd.R` (the
  `multigrain_pvals_gsd` object and its fields `pvals`, `nsim`, `m`, `K`, `alpha`, `info_frac`,
  `look_back`, `spending`, `tables`; the internal `.gsd_looks()`), `R/sim_pvals_gsd.R`,
  `src/graph_shortcut_gsd.cpp` (internal kernels returning `list(rejected, time)`),
  `tests/testthat/helper-gsd_reference.R` (`gsd_reference_direct()`, `gsd_pvals_matrix()`).
- Every file the roadmap says you will touch, before changing any of them:
  `src/graph_shortcut.cpp`, `R/objective_function.R`, `R/optimisation.R`, `R/trial_success.R`,
  `R/post_optim_processing.R`, `R/calc_power.R`, `R/sim_pvals.R`, `R/control_prepare.R`,
  `R/check_types.R`, `DESCRIPTION`.

## The rules that matter most

**LOCKED decisions are binding.** If you think a locked decision is wrong, stop and say so.
Do not implement something different, and do not implement it "with a small improvement".

**OPEN decisions are not yours to close.** If the record leaves a decision open, or the plan is
ambiguous at the point you need it, stop and ask. Guessing here produces work that has to be
thrown away.

**Gate by gate.** Work through phases P0 to P5 in order. At the end of each phase, stop, report
what you did and what the gate check showed, and wait. Do not run ahead.

**Report divergence.** If what you find in the repo contradicts what the record assumed, say so
before working around it.

**The `_gsd` family is additive.** Do not modify any existing exported function or the existing
kernel `graph_shortcut()`. Where the record offers "share a helper or duplicate", duplicate.
Branch `plan-for-sparsity` is rewriting `create_obj_func()` and `.try_prune()`; you must not
create merge conflicts there.

## Phases and gates

The record is authoritative; this is the summary.

- **P0 Transform. DONE (`e4ec646`).** `transform_pvalues_gsd()` and the per-hypothesis boundary
  tables built with `gsDesign::gsBound1()`. Gate: Maurer–Bretz (2013) Table 1 boundaries and
  Table 2 repeated p-values reproduced to four significant figures against the record's Appendix A
  `gsDesign` values (the paper's Table 2 digits match only to about 3.5 s.f., see record 4.1
  `[Rev 2026-09-16]`); the monotonicity assertion fires on a deliberately non-monotone spending
  function; grid inverse within 1e-5 of `uniroot()` at 200 random p-values; agreement with
  `graphicalMCP::repeated_p()` and `sequential_p()` within 1e-5 **absolute** at K = 2
  (`skip_if_not_installed("graphicalMCP")`, but the package is installed so the tests run).
- **P1 Simulator. DONE (`96dba9f`).** `simulate_pvalues_gsd()`. Gate: empirical correlations
  match the canonical joint model within 5e-3 at N = 1e6; matured columns identical to the
  maturity column. Takes `power_nominal` as `simulate_pvalues()` does (record item 1).
- **P2 Kernel. DONE (`d0b6556`).** `graph_shortcut_gsd()` and `graph_shortcut_gsd_parallel()`.
  Gate: with K = 1 the `rejected` output is `identical()` to `graph_shortcut()`; exact match to
  the R reference implementation in the record's appendix on N = 400 trials with mixed
  `look_back`; parallel output identical to serial at 1, 2, 4 and 8 threads; agreement with
  `graphicalMCP::graph_test_shortcut_gsd()` on the Maurer–Bretz case study.
- **P3 Gain. DONE (`0908406`).** `trial_success_gsd()` in `R/trial_success_gsd.R`. Gate: manuscript
  Example 5 and the three supplement forms compile and evaluate correctly on hand-built time
  matrices; `d(0)` is 0; snapshot tests of the generated C++. `[Rev 2026-09-17]` The parser
  could not copy the fixed-sample placeholder substitution (`%AND%`/`%OR%` cannot sit next to
  `==` in R); it parses natively, so precedence is R's. See record 4.6, "Parsing and
  precedence". The object's fields are `func`, `m`, `K` (may be `NULL`), `objective`,
  `cpp_code`, `tables`; class `c("multigrain_trial_success_gsd", "multigrain_trial_success")`;
  `is_trial_success_gsd()` is the internal predicate. The compiled `func` takes the kernel's
  `time` matrix (`IntegerMatrix`, 0 = never) and nothing else.
- **P4 Optimiser and post-processing.** `graph_optimise_gsd()`, `create_obj_func_gsd()`,
  `prune_graph_gsd()`, `calc_power_pvals_gsd()`. Gate: manuscript Figure 3b (optimal `w_PFS`
  against the OS/PFS value ratio at δ = 1, 0.75, 0.5) reproduced within 0.02 at N = 1e5;
  pruning never lowers the gain; `graph_optimise()` on `main` and on your branch give
  `identical()` results for the same seed and inputs.
- **P5 Documentation.** roxygen for every new export, pkgdown reference group, a GSD article
  under `vignettes/articles/`, the `NEWS.md` entry, `DESCRIPTION` (`Imports: gsDesign`;
  `Suggests: graphicalMCP (>= 0.3.0)`).

## Constraints

- **Complete files, not diffs**, unless a patch is explicitly requested.
- **Testing:** `testthat::test_file()` on the specific files affected. Never
  `testthat::test_package()`; that is reserved for CI on the PR. `[Rev 2026-09-16]` Set
  `Sys.setenv(NOT_CRAN = "true")` first, or every `skip_on_cran()` gate (P1's N = 1e6 block,
  P4's Figure 3b) silently skips when run from `Rscript`. Tests that reach the internal kernels
  need `devtools::load_all()` or `test_file(..., package = "multigrain")`. Bash on the Windows
  build machine truncates commands over about 8k characters: write scripts to files and run them
  with `Rscript`. If `testthat` rewrites `tests/testthat/_snaps/*.md` with CRLF endings when you
  run existing test files, revert them before committing, and `git add` your own paths
  explicitly rather than `git add -A`.
- **No large runs.** m ≤ 4 hypotheses, nsim ≤ 1e4. The one exception is the P4 Figure 3b
  reproduction (m = 2, K = 2, N = 1e5, a grid over a single parameter). Do not run `zhan_m5`,
  `cvot_m6` or `split_m8`. No HPC, no Slurm.
- If anything under `src/` changes, recompile and reload properly before testing; a stale
  shared object will make a broken build look like a passing one.
- Run `lintr` on changed files. Fix findings surgically; do not reformat surrounding code.
- Regenerate documentation from roxygen rather than hand-editing anything in `man/`.
- Two pre-existing issues you will notice. Do not fix them here; mention them in the report and
  do not replicate them in the `_gsd` code: `pvals[sample(nsim), ]` in `R/optimisation.R`
  permutes the first `nsim` rows rather than sampling `nsim` of all rows; unary minus in
  `parse_and_transform()` (`R/trial_success.R`) indexes a second argument that does not exist.

## Verification

Treat your own reasoning about what the existing code does as a hypothesis and check it by
running it. This applies with particular force to anything involving `graphicalMCP`, whose
conventions for matured hypotheses and `NA`-padded looks differ from the record's, and where
confident-but-wrong claims have been made before. The conventions established in P0 and P2 are
listed in record section 4.9 (`[Rev 2026-09-16]`): wrap spending functions so they return a
plain numeric vector; no column names on `p`; `decision_at` is written for tested-not-rejected
hypotheses; the oracle uses `<=`; `repeated_p()` is absolute-1e-6 and not reproducible at
K >= 3; our repeated p-values are capped at 1 above `alpha`.

`[Rev 2026-09-17]` For P4, from P3: when the gain object carries tables, its compiled function
indexes a C array with the time value, unchecked. Assert `gain$K == pvals$K` whenever `gain$K`
is not `NULL` (in `create_obj_func_gsd()` and `calc_power_pvals_gsd()`), and never call
`gain$func()` on anything but the kernel's `time` matrix. A `multigrain_trial_success_gsd`
object passes `check_trial_success()` because it inherits the fixed-sample class; P4 must
distinguish the two with `is_trial_success_gsd()` where it matters (a fixed-sample gain
expects a `LogicalMatrix` of rejections, not the time matrix). When writing cli messages that
mention comparison operators, interpolate them as values (`{.code {ops}}`): a literal `<`
inside cli markup is read as an internal delimiter and the message fails to format.

`[Rev 2026-09-16]` For P4: the kernel silently treats `NA`/`NaN` as "never reject"; put a single
`anyNA()` assertion in `create_obj_func_gsd()`. The no-change gate of record 4.10 must allow
`src/RcppExports.cpp` as well as `R/RcppExports.R`, `NAMESPACE` and `DESCRIPTION`. Subsampling
the transformed array must use `drop = FALSE` on the first dimension; note that subsetting the
raw simulator output drops its `info_frac` attribute, so subsample the `multigrain_pvals_gsd`
object, not the array. The record's Appendix B `P <- lapply(...)` block is not the joint model
of 4.7 and must not be used as a reference for the simulator.

The claims most likely to be quietly false, and which you must demonstrate by running code:

1. K = 1 is bit-identical to the fixed-sample kernel, including rows with p below 1e-12.
2. A hypothesis whose data are complete at an interim analysis can still be rejected at a later
   analysis when recycling gives it more alpha, with decision time equal to that later analysis
   (Example 5, τ_PFS = 2). To reproduce this with `graphicalMCP`, the matured hypothesis must be
   run with `look_back = TRUE`.
3. A row with p = 5e-11 and an allocation of 2.5e-11 (hypothesis weight 1e-4 times edge 1e-5
   times α) does not reject.
4. A hypothesis with weight exactly 0 never rejects, whatever its transformed p-value.
5. Evaluating the same (w, G) twice on the same transformed array gives the same gain to machine
   precision.

The record's adversarial-checks section lists what the review agent will attempt. Run them
yourself first.

## Git

- Commit at gate boundaries, not in one lump at the end.
- **Do not open the PR.** An adversarial review runs against this branch first. Push `gsd-build`
  and stop there.

## Report at the end

Write to `dev/gsd_report.md` (P0 to P2 are already reported in `dev/gsd_progress.md`; continue
that file or start `dev/gsd_report.md` for P3 onward):

- What changed, file by file.
- Where the implementation diverged from the record, and why.
- What you could not verify, and what would be needed to verify it.
- Anything you noticed that the record did not anticipate.
