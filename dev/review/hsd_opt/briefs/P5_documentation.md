# Brief P5: documentation, design record, call flow, layperson explanation

Branch `gsd-build-spending`, repo `/home/user/multigrain-spending`. You implement phase P5 of
`dev/hsd_opt_build_plan.md`. P1–P4 are committed; no code changes in this phase except
roxygen comments and the one P1 default noted below. **Report, do not commit.** Complete files.


## House rule: no `<<-` (user decision, 2026-10-14)

Do not use `<<-` anywhere. Any state a closure must mutate (a working matrix, an archive,
a counter, a cache) lives in one explicit environment created inside the factory body
(`state <- new.env(parent = emptyenv()); state$x <- …`) and is updated with ordinary `<-`
(`state$x[i] <- v`). d03 showed that assignment through an environment never copies a
matrix, the same as `<<-`; a matrix must still be built inside the factory and never bound
to a second name. Add a test that greps your new R files for `<<-` and expects none.

## Read first (~60 minutes)

1. `dev/hsd_opt_build_plan.md` §4 "P5" (the deliverable list is binding), §1 (every
   reconciliation — they are the "Where the implementation diverged" material).
2. `dev/hsd_opt_design_revised.md` in full (you write record §4.12 from it and from what was
   built), `dev/review/hsd_opt_discovery.md` in full (the `[measured]` numbers).
3. `R/transform_pvalues_gsd_opt.R`, `R/objective_function_gsd_opt.R`,
   `R/optimisation_gsd_opt.R`, the touched lines of `R/transform_pvalues_gsd.R`,
   `R/optimisation_gsd.R`, `R/graph_optimal.R`, `R/check_gsd.R`; the four `_opt` test files.
4. Conventions: `dev/gsd_design_record.md` §4.10–4.11 and §10 (the pattern for a new 4.12,
   a `[Rev]` paragraph and item 16); `dev/callflows/graph_optimise_gsd.md` (format);
   `dev/callflows/stopping_rule_gsd.md` (a shorter example); `NEWS.md` top; `_pkgdown.yml`;
   `vignettes/articles/trial_success.qmd` (chunk style `#| label:`, `#| eval: false`, pasted
   `#>` output); `dev/review/gsd_user_testing/README.md` §7.4–7.5 and §8 (the running
   example and its numbers); `dev/gsd_progress.md` "Stopping rules" entry (format only — the
   orchestrator writes the new progress entry, not you).

## Facts from the build that the documents must reflect (read before writing)

- Result fields actually attached by the driver (`R/optimisation_gsd_opt.R`, the `spending =
  list(...)` in `.graph_optimise_gsd_opt()`): `spending_param`, `hyp_spending_param`, `spending`,
  `spending_bounds`, `spending_step`, `spending_group`, `spending_profile` (columns `param, value,
  gain, gain_fixed_graph, diff, se, n_eval`, attribute `graphs`), `spending_search` (`param, value,
  visits, best_fitness_subsample`), `profile_converged`.
- Stage L details that differ from the design text: step 1 is `.hsd_refine_node()` on rows
  `idx_local` pinned once (the fixed local driver re-subsamples internally, so it is not called);
  an axis whose own sweep improved the optimum is recomputed once against the final optimum with
  **no refinement at the optimum's own node**, so `gain == gain_fixed_graph` and `diff == se == 0`
  exactly at theta*; `profile_converged` is FALSE if a recomputed axis shows `diff < 0` or `d >= 2`
  rounds did not settle in 3.
- Table building shows a `cli::cli_progress_bar()` (not a step) after 2 s, verbosity-gated;
  a step at default verbosity leaked into snapshots.
- The objective closure keeps its working matrix, archive and last node in one environment
  created inside the factory (no `<<-`); the matrix is updated by take-out/put-back because a
  direct `state$mat[, cols] <-` inside a function copies it (P2 measured; d03's probe was top-level).
- Benchmarks at nsim 1e5 (`dev/review/hsd_opt/gate_p3_benchmarks_output.txt`): M1 joint 0.693038
  at gamma −0.25 (fixed spending at the start node −1.5: 0.691238); README §2 design joint 2.43998 at
  gamma −4 (on the lower bound; fixed at −1.5: 2.40569), wall 8.0 min against the design's 3.5–6.5
  min projection; d = 2 (M1, H1 and H2 tagged separately, bounds −2..0, step 0.5, nsim 2e4) joint
  0.694014 at (−2, 0) equal to the 25-pair product-grid maximum.
- Phase D decisions and numbers: `dev/review/hsd_opt_discovery.md` "Decisions" — quote them.
- Reconciliations the record's 4.12 must state as departures from the design document: the plan's
  §1 table (11 rows) plus the three above.
- `R/graph_optimal.R`'s `spending = NULL` argument and the two guarded display calls; the dispatch
  in `R/optimisation_gsd.R`; the handoff and `.gsd_spend()` guard in `R/transform_pvalues_gsd.R`;
  the hint in `R/check_gsd.R` — these are the §4.10 touchpoints.

## Deliverables

1. **roxygen** on `sf_hsd_opt()`, `resolve_pvalues_gsd()`, `print.multigrain_pvals_gsd_spec`
   (if documented), `transform_pvalues_gsd()` (`@param spending` gains the tagged case and the
   spec return; `@returns`), `graph_optimise_gsd()` (`@param pvals` accepts a spec; `@returns`
   lists the new fields; a paragraph on the search: lattice, GA with parameter genes, warm-
   started profile, what `spending_profile` and the plateau mean; the `choose_graph()` note —
   if the GA winner is chosen, the profile is still Stage L's; the winner's-curse caveat and
   the fresh-draw check with `resolve_pvalues_gsd()`; the `grid_size = 512` error figure from
   the discovery report). `@examples` with `\donttest{}` for anything slow. Then
   `LC_ALL=C.UTF-8 Rscript -e 'devtools::document()'` twice (second a no-op).
2. **`_pkgdown.yml`**: "Group sequential designs" gains `sf_hsd_opt` and `resolve_pvalues_gsd`
   after `transform_pvalues_gsd`. Run `pkgdown::check_pkgdown()` on a scratch copy of the repo
   with the `.qmd` articles renamed `.Rmd` (the `quarto` R package is not installed) and report
   "No problems found" or what it found.
3. **Article** `vignettes/articles/trial_success.qmd`: new `### Optimising the spending
   parameter` after "Stopping for efficacy". Running example: README §7.5 design under M1
   (`sfLDOF` on H1, `sf_hsd_opt()` on H2, `stopping_rule_gsd(r1)`, `nsim = 1e5`, seed
   20261007). Produce the output with a script `dev/review/hsd_opt/article_example.R`
   (`Rscript`, `num_threads = 2`, once; keep its output file) and paste it as `#>` lines in
   `#| eval: false` chunks. Show: building the spec and its `print()`, the optimiser call, the
   result's `print()`, the profile table from `summary()`, and one paragraph reading the
   plateau. Tie the numbers to §7.5 (gain 0.6707 → 0.6903 for a fixed Pocock H2; what the
   optimiser finds).
4. **`NEWS.md`** `## New functionality`: three bullets (`sf_hsd_opt()`, `resolve_pvalues_gsd()`,
   `graph_optimise_gsd()` on a spec and its result fields), house wording.
5. **`dev/gsd_design_record.md`**: **§4.12 Spending-parameter optimisation `[Rev 2026-10-14]`**
   inserted before `## 5.` in the 4.11 pattern: **Decision.** LOCKED (user decision,
   2026-10-14, `dev/implement_hsd_opt.md`); italic-labelled bullets *Tagging*, *Spec object and
   eager tables*, *Lattice and snapping*, *Search (Stage G, Stage L)*, *Reporting*, *Dispatch
   and touchpoints*, *Error control*, *Defaults (from Phase D)* — each stating what was built,
   with the measured numbers; **Rationale.**; **Alternatives.** (blending, penalties, fixed-graph
   sweep as headline, outer product grid, transform inside the objective, inheriting the
   class, transformed array on the result). A `[Rev 2026-10-14]` paragraph in §4.10 listing the
   touchpoints on existing files (`transform_pvalues_gsd.R`, `check_gsd.R`,
   `optimisation_gsd.R`, `graph_optimal.R`) and that the identity gate still holds. §10: item
   16 — the `start_graph` row-1 bug (open, reported) and anything the discovery report left
   open. In §4.5, append one sentence to the "What would change under joint design
   optimisation" paragraph pointing to 4.12.
6. **`dev/callflows/graph_optimise_gsd_opt.md`** in house format: "# Call flow:
   `graph_optimise_gsd()` on a spending spec and its helpers"; intro prose naming the three new
   files, the four touchpoints, the fixed helpers reused unchanged (`.gsd_boundary_table()`,
   `.gsd_invert()`, `create_obj_func_gsd()` for Stage L, `choose_graph()`,
   `prune_graph_gsd()`, `calc_power_pvals_gsd()`, `graph_optimal()`), and the kernel's
   document; a notation paragraph (m, K, d, M_r, the spec fields, the gene vector layout);
   "## The call chain at a glance" numbered (handoff → spec builder; dispatch; validation and
   control; Stage G; Stage L step 1; profile sweep and rounds; `choose_graph()`; resolve,
   prune, power, `graph_optimal()`; display); `---`; one `## fn() -- description` per function
   in execution order — `sf_hsd_opt()`, `.transform_pvalues_gsd_spec()`, `.hsd_groups()`,
   `.hsd_lattice()`, `.hsd_node_tables()` (+ `.hsd_check_invertible()` as a `###`),
   `.gsd_invert_hyp()`, `resolve_pvalues_gsd()`, `.hsd_labels()`, `.graph_optimise_gsd_opt()`,
   `.hsd_bank()`, `create_obj_func_gsd_opt()` (+ `.hsd_node_index()`),
   `.graph_optimise_ga_gsd_opt()`, `.graph_optimise_local_gsd_opt()`, `.hsd_refine_node()`,
   `.hsd_encode_graph()`, `.hsd_profile_axis()`, `.hsd_block_se()`, `.hsd_plateau()`,
   `.hsd_cat_print()`/`.hsd_cat_summary()` — each with `### Purpose`, `### Inputs and output`,
   `### How it works`, `### Why it is done this way`, `### Edge cases`. Verify every claim
   against the code (quote line numbers sparingly; they go stale). Add a `[Rev 2026-10-14]`
   paragraph to `dev/callflows/graph_optimise_gsd.md` after its intro: the dispatch on a spec
   and a pointer to the new document.
7. **`dev/hsd_opt_explained.md`** — a layperson's explanation, plain English, about 150 lines,
   at most one formula, no code beyond function names. Sections: *What a spending function is*
   (an error budget of 2.5 % and a schedule for how much may be used at each look; why looking
   early costs something); *What the HSD parameter does* (−4: save nearly everything for the
   final analysis, 1: spend generously early; at the halfway look about 12 % vs 62 % of the
   budget); *Why it became a first-order choice* (README 7.5: once a trial can stop after the
   primary endpoint succeeds, letting the secondary spend early raised the design's value from
   0.6707 to 0.6903 — 14 times the gap between the best graph and the textbook Holm graph —
   and moved the best weight on the primary from 0.32 to 0.46); *What "optimising it jointly"
   means* (the graph and the parameter are chosen together because each changes the other's
   best value); *How the search works* (a short list of candidate values, every one checked in
   advance to be a valid test; a global search that tries graph and parameter together; then
   the graph re-tuned at every candidate to draw the profile); *How to read the result* (the
   chosen value, the profile table, the two plateau intervals and what a wide plateau means);
   *What it does not change* (the error-rate guarantee, the simulation engine, any result that
   does not use `sf_hsd_opt()`); *Limits* (chosen on the same simulated trials that score it,
   so check on fresh draws; grid granularity; run time). Use the article's actual numbers for
   the running example.
8. **Baseline for the orchestrator's progress entry**: in your report, list file by file what
   changed in P5 and the gate counts.

## Gate
Clean build by the orchestrator; your own checks: `document()` twice (no-op);
`pkgdown::check_pkgdown()` on the scratch copy; `test_file()` on the four `_opt` test files and
`test-graph_optimal.R` (counts unchanged); `lintr::lint()` on every R file you touched
(roxygen edits can trip `line_length_linter`); no hand edits under `man/`; the article renders
is **not** required (quarto absent) — but every pasted output must come from a real run.

## Report
Files; gate; the pkgdown result; the article example's wall time and headline numbers;
divergences; not verified. Do not commit.
