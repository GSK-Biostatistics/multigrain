# Prompt: joint optimisation of a Hwang-Shih-DeCani spending parameter

> For a fresh Claude Code session in a clone of `multigrain` on branch `gsd-build-hsc-opt`,
> created from the head of `gsd-build-stopping` (`git fetch origin gsd-build-stopping &&
> git checkout gsd-build-stopping && git checkout -b gsd-build-hsc-opt`). Develop and push on
> `gsd-build-hsc-opt`; do not open a pull request. This prompt is deliberately light on
> specification: the first phase is discovery, and it ends with a report and a decision by
> the user before anything is built.

## Setting

- Toolchain: see `dev/cloud_setup.md` for the container setup script (Ubuntu apt binaries
  plus a few CRAN source installs). Before any code, confirm from a clean build
  (`devtools::clean_dll(); devtools::install(quick = TRUE); devtools::load_all()`, with
  `Sys.setenv(NOT_CRAN = "true")`) that `tests/testthat/test-transform_pvalues_gsd.R` and
  `tests/testthat/test-optimisation_gsd.R` pass. If they do not, stop and report.
- House conventions: `dev/implement_gsd.md` ("Read first", "Constraints"), plus the lessons
  in `dev/gsd_progress.md` ("Stopping rules" entry: `devtools::document()` under
  `LC_ALL=C.UTF-8`, `lintr::lint()` on every touched file, snapshot line endings, never
  `testthat::test_package()`, m <= 4 and nsim <= 1e4 in tests, recompile after any `src/`
  change, `air` may be unavailable).
- Read first: `dev/gsd_design_record.md` sections 4.1 (the transform), 4.5 (why the search
  is unchanged, and the paragraph "What would change under joint design optimisation"),
  4.7 and 4.11; `dev/review/gsd_user_testing/README.md` section 7.5 and suggestion 4
  (the spending function of the secondary endpoint moved the gain more than the graph
  did); `R/transform_pvalues_gsd.R` (`.gsd_boundary_table()`, `.gsd_invert()`,
  `.gsd_spending_list()`); `R/optimisation_gsd.R` and `R/objective_function_gsd.R`;
  `dev/callflows/graph_optimise_gsd.md`.
- Orchestration as before: the session plans, gates and reviews; implementation is
  delegated to one agent per phase, each told to report rather than commit; every gate is
  re-run from a clean build before a commit.

## Goal (binding)

Optionally optimise, jointly with the graph, a single parameter gamma of the
Hwang-Shih-DeCani spending family (`gsDesign::sfHSD(alpha, t, param = gamma)`):

- Optional: the existing fixed-spending path must be byte-identical to today. The identity
  gates of `dev/gsd_identity_check.R` (fixed-sample) and the GSD test files stand.
- One gamma, shared by the hypotheses whose spending function depends on it. Which
  hypotheses those are is expressed through the `spending` argument: a spending function of
  `(alpha, t, gamma)`, or a list of `m` functions mixing those with ordinary fixed
  `(alpha, t)` functions, so one gamma can apply to every hypothesis or only to, say, the
  secondary endpoints while the primary keeps O'Brien-Fleming. A small helper that returns
  the HSD function in that form is welcome.
- Bounds: default `c(-4, 1)` (gamma = -4 is O'Brien-Fleming-like, gamma = 1 Pocock-like).
  The user may set other bounds, with lower >= -10, upper <= 2 and lower < upper; anything
  else aborts with a message naming the limits.
- The result reports the chosen gamma and enough of the gain-versus-gamma relationship for
  the user to judge how flat it is, and the spending function actually used is visible on
  the result and in `print()`/`summary()`.
- Stopping rules (`stopping_rule_gsd()`) and every reporter field keep working; the gamma
  search must not touch the kernel.
- Error control is unaffected (the transform at any gamma is a valid repeated-p-value
  transform; the well-ordering check of 4.1 step 3 must still run at every gamma tried).

The API names, where the search lives (an argument of `graph_optimise_gsd()`, a wrapper, or
something else) and the internal structure are NOT specified: they are the output of the
discovery phase.

## Phase D: discovery (report, then stop and wait for the user's decision)

The user's first instinct was that gamma should simply be another real-valued parameter of
the evolutionary strategy, alongside the weights and transition matrix. The obstacle is
cost, not principle: a change of gamma changes the repeated p-values, so every objective
evaluation would need a transform. Measured on 2026-10-13 (4 cores, nsim = 1e5, m = 3):
one transform costs about 3.2 s at K = 2 and 4.6 s at K = 3 with `grid_size = 1024`,
0.76 s at 256 and 0.35 s at 128, against about 13 ms per kernel call. A default
`graph_optimise_gsd()` run makes thousands of objective evaluations. Measure this yourself
and answer:

1. **Where the transform's time goes.** Split the cost per gamma into the boundary tables
   (`gsDesign::gsBound1()` once per grid level per hypothesis, independent of nsim) and the
   interpolation (`.gsd_invert()`, linear in nsim), at grid sizes 128, 256, 512, 1024, at
   nsim 1e4 and 1e5, K = 2 and 3, m = 2 and 3. Report the interpolation error against the
   1024 grid at the coarser sizes (the record's Appendix B method: max relative error at
   200 random p-values) so that a coarser grid can be judged.
2. **How many evaluations the optimiser makes.** Instrument `create_obj_func_gsd()` (a
   counter in the closure, in a scratch copy) for a default run at m = 3, K = 2, nsim 1e5,
   and for `global_search = FALSE`. Report the GA count, the COBYLA count and the pruning
   count separately.
3. **Candidate designs**, each with a cost estimate from 1 and 2 and a list of what it
   changes:
   - (a) gamma as a real parameter of the GA and COBYLA, with a transform inside the
     objective. Only feasible if the per-evaluation transform can be made cheap. Explore:
     memoising boundary tables keyed on gamma (COBYLA re-evaluates nearby gammas, the GA
     does not); a coarse grid; precomputing boundary tables on a gamma grid once and
     interpolating the table in gamma (measure the error of that interpolation against the
     exact transform, and whether well-ordering survives it); caching the transformed array
     for the most recent gammas. State the per-evaluation cost you achieve and the total
     run time it implies.
   - (b) an outer loop over a gamma grid (say 11 points from lower to upper), one transform
     and one full inner `graph_optimise_gsd()` per point, pick the best, return the profile.
     This is the record's own extension path (section 4.5) and the user testing's
     suggestion 4. Optionally a one-dimensional refinement around the best point, with the
     inner optimiser warm-started from the best graph. State the total run time.
   - (c) a hybrid: the GA on the grid-cached arrays (gamma snapped to the nearest cached
     value, so the GA sees a piecewise-constant objective in gamma, which it tolerates),
     then COBYLA on the continuous gamma with a memoised transform.
   For each: how common random numbers hold (the raw draws are fixed, so gamma changes only
   the transform, section 4.5), how the gain behaves in gamma (piecewise-constant on fixed
   draws, so a derivative-free method is required), and what the result object should
   carry.
4. **Interactions to check and report**: the well-ordering check at extreme gammas within
   [-10, 2] for the information fractions the record uses; `K` inference of discount tables
   and `.gsd_check_gain_dims()` (should be gamma-independent); the stopping-rule table
   (gamma-independent by construction); `control_prepare_dims()` and subsampling (the
   subsample index must be drawn once and reused across gammas, or the objective is not
   comparable between gammas); `print()`/`summary()` of `multigrain_pvals_gsd` (the
   spending label must show the gamma used).
5. **Recommendation.** One design, with numbers, the API you propose (names, where the raw
   p-values and the transform settings enter, what the result carries), the phases and
   gates you propose in the style of `dev/implement_gsd.md`, and what you would test
   against what oracle (for example: at a fixed gamma the joint path must reproduce
   `graph_optimise_gsd()` on `transform_pvalues_gsd()` at that gamma exactly, given the
   same seed and subsample).

Write the report to `dev/review/hsd_opt_discovery.md` with the scripts that produced the
numbers under `dev/review/hsd_opt/`, commit it on `gsd-build-hsc-opt`, and stop.

## After the decision

Implement in phases with gates, one agent per phase, as in `dev/implement_gsd.md`. The
documentation deliverables are the usual: roxygen and `devtools::document()`, the pkgdown
reference section, a subsection of `vignettes/articles/trial_success.qmd` or a new article,
`NEWS.md`, design record section 4.12 (decision, rationale, alternatives, with the discovery
numbers), a call-flow document, and an entry in `dev/gsd_progress.md` in its existing format.

## Report

Append to `dev/gsd_progress.md`: what changed file by file, where the implementation diverged
from the decision and why, what could not be verified, and anything the discovery phase did
not anticipate.
