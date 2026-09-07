# Implementation report: `graph_simplify()`

Branch `plan-for-sparsity`, implementing `dev/sparsity_design_record.md`. This
records what changed, where I departed from the record or the plan and why, what
I could not verify, and what the record did not anticipate. It is written for
the adversarial review that runs next.

Environment: R 4.6.1 (Windows 11), Rtools 4.5, `pkgload::load_all()` from the
working tree. `lintr` 3.4.0, `devtools` 2.5.2 and `cyclocomp` were installed to
run the gates (none were present).

---

## 1. What changed, file by file

### `R/objective_function.R`

* `.lexico()` — the scalar rule: feasible graphs score `u - edge_price *
  n_edges`, infeasible ones `u - edge_price * (n_free + 1)`.
* `.trial_success_range()` — exact min and max of the compiled measure over all
  `2^m` rejection patterns.
* `.make_lexico_scorer()` — **not in the record**; see §2.4.
* `create_obj_func()` gains `gain_tolerance`, `ref_graph`, `u_range`. When
  `gain_tolerance` is `NULL` the penalty offset is exactly `0` and every branch
  returns what it did before.

### `R/mutation_helpers.R`

* `.g_param_rows()` — row index of each free transition parameter.
* `.zeroing_row_targets()` — **not named in the record**; see §2.5.
* `.make_cauchy_mutation_multi()` gains `p_zero`, `param_rows`, `row_target`.
  The previous closure is kept verbatim as an inner function and returned
  unwrapped when `p_zero = 0`. The comment block carries §3.6's two paragraphs
  as the record requires.

### `R/optimisation_start.R`

* `.encode_graph()` — encodes a graph and guards derived entries that should be
  zero so they decode to `5e-6` (`5e-5` for a derived weight).
* `.build_simplify_seeds()` — reference first, then one-edge-removal
  neighbours, then `.build_start_matrix()`'s seeds, then the stage-1 population
  when its width matches; deduplicated and truncated to `popSize`.

### `R/optimisation.R`

* `.graph_optimise_ga()` gains `suggestions`, `p_zero`, `objective_args`, and
  returns `ga_objective`.
* `.graph_optimise_local()` gains `objective_args`, and returns
  `local_objective`.
* A guard after `utils::modifyList(immutable_global_args, global_opts)` restores
  `mutation` and `suggestions` when the simplification is using them; see §3.3.
* `graph_optimise()` passes `alpha` to the constructor. That is its only change.

### `R/choose_graph.R`

Compares `ga_objective` / `local_objective`, falling back with `%||%` to the
trial-success values so results that predate those elements still work — which
is what every existing `choose_graph()` test supplies.

### `R/post_optim_processing.R`

* `.prune_edges_best_first()` — evaluates every removable edge on the full
  sample, keeps candidates that stay at or above the threshold and reduce the
  edge count, removes the highest-gain one, repeats.
* `prune_graph()` gains `threshold` and returns `prune_loss`.

### `R/graph_optimal.R`

* `new_graph_optimal()` and `graph_optimal()` gain `alpha` and `sparsity`, both
  `NULL` by default and appended at the end of the list.
* `summarise_sparsity()` and its hooks in `print()` and `summary()`.

### `R/utils_verbosity.R`

* `.resolve_verbose()` — the coercion `graph_optimise()` performs inline. See
  §2.1.

### `R/graph_simplify.R` (new)

* `graph_simplify()`, `.simplify_result()`, `.trial_success_is_live()`, and the
  constant `.simplify_p_zero`.

### Everything else

* `NAMESPACE`, `man/graph_simplify.Rd` — regenerated with
  `devtools::document()`.
* `_pkgdown.yml` — `graph_simplify` under **Optimisation**.
* `vignettes/articles/get-started.Rmd` — a "Simplifying the graph" section.
* `NEWS.md` — the record's §10 wording. No bug fixes are part of this change.
* `data/graph_optimal_example.rda` — `alpha` and `sparsity` added in place.
* Tests: additions to `test-objective_function.R`, `test-mutation_helpers.R`,
  `test-optimisation_start.R`, `test-post_optim_processing.R`,
  `test-graph_optimal.R` (and its snapshot file), plus a new
  `test-graph_simplify.R`.

---

## 2. Where the implementation diverged, and why

### 2.1 `.resolve_verbose()` did not exist

The record's §5.6 calls it. It is not in the package; `graph_optimise()` and
`trial_success()` each inline the same seven lines. I added it to
`R/utils_verbosity.R` and call it **only** from `graph_simplify()`, leaving both
existing bodies untouched so "`graph_optimise()` changes in two ways only" holds
literally. Extracting the duplication is a tidy-up for a separate change.

### 2.2 `D` renamed to `edge_price`

`.lintr` runs `object_name_linter(styles = "snake_case")` with exceptions only
for `^G`, `G$`, `^gMCP`. `D` is rejected. `.lexico()`'s argument is
`edge_price`, documented as "written `D` in the design record".

### 2.3 `gc`, `ts` and `source` renamed inside `graph_simplify()`

§5.6 uses all three. `object_overwrite_linter` rejects `gc` (base) and `ts`
(stats); `source` trips both that and `undesirable_function_linter`. They are
`constraints`, `trial_success` and `result_source`. The **`sparsity$source`
element keeps its documented name** — only the local variable changed.

### 2.4 `.make_lexico_scorer()` was added

§3.7 requires `ga_objective` / `local_objective` to be the lexicographic score
**on the full sample**, but §5 gives no code for it: `create_obj_func()`'s
closure computes its threshold from the p-values *it* captured, which are the
subsample. `.make_lexico_scorer()` builds the full-sample threshold from the
reference and returns a scoring function; with an empty `objective_args` it is
the identity on trial success, so `choose_graph()`'s arithmetic is unchanged for
`graph_optimise()`.

Cost: one extra full-sample shortcut evaluation in each of the two optimisers,
to get `u_ref` on the full sample. `graph_simplify()` has already computed that
number and could pass it down; I left the two functions self-contained to keep
their signatures independent of the caller. Worth revisiting if the full sample
is large.

### 2.5 The mutation factory gained `row_target`

§5.4's sketch rescales a row to `1 - 5e-6`, and its own closing note says that
for rows with pinned non-zero entries the target must be `1 - 5e-6` **minus the
fixed sum**. That needs the constraint, which the sketched signature does not
receive. I added `.zeroing_row_targets(trans_constraint)` and a `row_target`
argument. Adversarial check A11 is covered by a test.

### 2.6 `opt_source` for the fallback

§3.1 says `"reference"`; §5.6's sketch would produce `"simplify:reference"`. I
followed §3.1, the specification.

### 2.7 `start_graph` seeds

§3.1's LOCKED signature has no `start_graph`, but §3.5 lists user-supplied start
graphs among the seeds. `.build_simplify_seeds()` receives the **reference
object's own stored `start_graph`**, which is the only one in scope.

### 2.8 The `shortcut()` helper in `create_obj_func()`

My plan said to keep the `if (use_parallel)` block inline. I used the record's
§5.2 helper instead, because the construction-time `u_ref` evaluation is a
second call site for the same branch. §5.2 explicitly permits either. Confirmed
with the user mid-implementation. Identity is proven by test, and a benchmark
(2000 evaluations x 5 reps on 1e4 x 4 p-values) put helper and inline at 1.800 s
vs 1.780 s median with fully overlapping ranges — noise.

### 2.9 Two guards the record does not specify

* **Serialisation** (§4.1) — `graph_simplify()` aborts early with an actionable
  message if the trial success function cannot be evaluated.
* **Invalid search output** — if `choose_graph()` returns a graph that is not
  valid, pruning would abort inside `calc_power_pvals()`. `graph_simplify()`
  skips pruning and falls back to the reference, which §3.2 guarantees is
  feasible. Without this the documented fallback would not actually hold in that
  case.

Both follow §3.11's pattern of handling degenerate cases; neither changes the
search.

### 2.10 The example dataset was patched, not regenerated

Agreed with the user before implementing. `data-raw/graph_optimal_example.R`
needs `2^20 x 6` simulated p-values plus a full genetic algorithm — over this
task's run limit — and calls `cran_cores()`, which only exists as a test helper.
Loading the object, adding `alpha = 0.025` (the default it was built with) and a
`NULL` `sparsity`, and re-saving with `save(compress = "bzip2", version = 2)`
achieves exactly what regenerating would, with no re-optimisation and no RNG
drift, so every existing snapshot of that object still holds. The script needs
no change; re-running it would now produce both elements on its own.

### 2.11 One existing test expectation changed

`test-graph_optimal.R`'s `new_graph_optimal` test pinned the exact element
names. It now includes `alpha` and `sparsity`. That is the change §3.10 calls
for; no snapshot was regenerated, and the new snapshots in
`_snaps/graph_optimal.md` are purely additive.

---

## 3. What the design record did not anticipate

### 3.1 A compiled trial success function does not survive serialisation

**This is the most consequential finding.** `trial_success()` compiles its
measure with Rcpp, and the resulting function does not survive `saveRDS()` /
`save()`. A reloaded object errors with `NULL value passed as symbol address`
the first time the measure is evaluated. Verified three ways: the shipped
`graph_optimal_example`, the article's cached `.rds` graphs, and a freshly built
`trial_success()` object through a `saveRDS`/`readRDS` round trip.

It is pre-existing and package-wide, and was never noticed because
`graph_optimise()` builds the measure in-session and `print()`/`plot()` only
read `$objective`. But it bears directly on `graph_simplify()`, whose premise —
"`pvals` is supplied again because the object does not store it" (§3.1) —
implies the object may well come from an earlier session. **As things stand,
a saved optimised graph cannot be simplified.**

What I did: `graph_simplify()` checks up front and aborts with the exact rebuild
command, interpolating the stored expression:

```
The trial success function of `graph_optimal` cannot be evaluated.
i Compiled trial success functions do not survive being saved and reloaded, so
  an object read back from disk needs its measure rebuilt before
  `graph_simplify()` can use it: `graph_optimal$trial_success <-
  trial_success(0.25 * (2 * (r1 && r2) + r1 * r3 + r2 * r4))`.
```

I verified that remedy works. Note it only works when the expression is typed
as a **literal**: `trial_success()` captures with `rlang::enexpr()`, so a
variable holding the string fails ("Symbol `s` ... is not a rejection
indicator"), while a literal expression, a literal string, and
`rlang::inject(trial_success(!!str2lang(s)))` all work.

**Left for the record's owner:** `graph_simplify()` could recover silently by
recompiling from `$objective` via the injection form. That would make the
function work on saved objects transparently, which is arguably what §3.1
intends. I did not do it because it is a design decision, not an implementation
detail. The alternative — fixing serialisation properly in `trial_success()` —
is a separate bug fix in the O2/O3 category.

### 3.2 `prune_loss` can be negative

§3.9 calls it "the exact loss of `U` across the removals pruning accepted".
Best-first pruning picks the highest-gain feasible candidate each round, so it
climbs whenever it can, and the "loss" is then negative. On the fixture, at a
threshold equal to the starting gain it removed 8 edges and *raised* the measure
from 0.6649 to 0.6769. I left the value exact rather than clamping it, and
documented that it can be negative. It is a diagnostic; `sparsity$gain_loss` is
computed separately against the reference and is the number the cap is about.

### 3.3 `modifyList` lets user control settings override the search

`ga_args <- utils::modifyList(immutable_global_args, global_opts)` puts
`global_opts` second, so a user's `control_global(mutation = ...)` or
`suggestions = ...` silently replaces the zeroing mutation or the warm start —
the two things that let stage 2 drop an edge at all. I restore both after the
merge when they are in use. `graph_optimise()`'s behaviour is untouched, since
it passes `p_zero = 0` and `suggestions = NULL`.

### 3.4 The subsample is a permutation, not a subsample

Both internal optimisers do `pvals_sampled <- pvals[sample(nsim), ]`.
`sample(nsim)` is a permutation of `1:nsim`, so this takes the **first `nsim`
rows in random order**, not a random subsample of all rows. §2 of the record
says "a random subsample of `nsim_global` trials". Nothing in the design depends
on it — the threshold is still computed on the same rows the candidates are
judged on — but the record's description of existing behaviour is wrong, and
anyone reasoning about subsample noise should know. **Not fixed: out of scope.**

### 3.5 The article's cached graphs are too stale to use

`vignettes/articles/data/*.rds` hold objects whose `constraints` element still
carries the pre-0.3.0 class `graph_constraint` (not
`multigrain_graph_constraint`), on top of the dead measure of §3.1. They cannot
be passed to `graph_simplify()`. The new article section therefore re-optimises
briefly on a 2e4-row subsample with `global_search = FALSE` and says so in the
text. Refreshing those caches is a separate job.

### 3.6 Best-first and fixed-order agree only at the right threshold

§7's step-5 gate says "a threshold equal to the current gain reproduces the
current 7-edge result". "Current gain" has to mean the gain of the graph *going
in*. With the threshold set instead to the gain the fixed-order prune *achieved*,
best-first removes nothing at all: it is greedy one edge at a time, and no
single removal from the dense graph reaches that higher value, whereas the
fixed-order sweep gets there through a sequence of steps each measured against a
running best. This is the difference between the record's appendix check 7,
which uses a running-best acceptance rule, and §5.5's absolute threshold. With
the threshold at the input gain the two agree on edge count (4) and on gain
(0.6769 exactly); they need not reach the same matrix, because several graphs of
that size tie on the fixture sample.

### 3.7 Two testing traps worth recording

* `expect_equal(x, 5e-6, tolerance = 1e-12)` fails on a value correct to 3e-17
  absolute, because the tolerance is **relative**: 3e-17 / 5e-6 = 6e-12. Guarded
  derived entries must be checked with an absolute comparison.
* The guard makes a row's parameters sum to `1 - 5e-6`, so after thresholding a
  guarded row sums to `0.999995`, not `1`. `param_to_solution()` renormalises it
  later. Any assertion of "rows sum to one" on an encoded seed is wrong.

---

## 4. Verification

### 4.1 `graph_optimise()` is unchanged — the claim most likely to be false

Verified by running code, not by a passing test. A script
(`scratchpad/a1_unchanged.R`) runs `graph_optimise()` over **100 cases**: seeds
1-5 x `m` in {2,3,4} x free and constrained constraints x `global_search` on and
off x `num_threads` 1 and 2, on seeded fixtures with `nsim = 1e4`. It records
the graph, `$power`, `$solution`, `$control`, `$start_graph`, the GA's solution,
fitness, population and iteration count, the nloptr solution, objective and
iteration count, `names()`, and `.Random.seed` afterwards. The same script ran
against a `git worktree` of `main` and against the branch, and the two `.rds`
files were compared element by element with `identical()`.

```
cases: main 100  branch 100
--- differences outside alpha/sparsity/names --- NONE
--- .Random.seed differences --- 0 of 100 cases differ
added: alpha, sparsity      removed: (none)
branch alpha values: 0.025  branch sparsity all NULL: TRUE
VERDICT: PASS
```

This covers adversarial check A1.

### 4.2 Per-file gates

Every file was run with `NOT_CRAN=true` so the `skip_on_cran()` snapshot tests
actually execute; the first run of `test-objective_function.R` skipped four of
them until I noticed.

| File | Result |
|---|---|
| `test-objective_function.R` | 351 pass, 0 fail, 0 skip |
| `test-mutation_helpers.R` | all pass |
| `test-optimisation_start.R` | all pass |
| `test-post_optim_processing.R` | all pass |
| `test-optimisation.R` | 86 pass |
| `test-choose_graph.R` | 11 pass |
| `test-graph_optimal.R` | 41 pass |
| `test-plot_graph_optimal.R`, `test-calc_power.R` | all pass |
| `test-graph_simplify.R` (new) | 128 pass |

No existing snapshot changed content. `_snaps/graph_optimal.md` gained three new
entries and nothing else; the other snapshot files were only rewritten with
different line endings by testthat, and their blob hashes are unchanged.

Notable evidence beyond the gates:

* The disabled objective path is compared against a **verbatim copy of the
  pre-change closure** over 200 random encodings at `m` in {3,4}, with
  `expect_identical`, including `NA`, `NaN`, out-of-range and near-zero vectors
  so every early-return branch is exercised rather than only the shortcut path.
* `p_zero = 0` is compared against a verbatim copy of the previous mutation
  closure over 50 calls **and** on the resulting `.Random.seed`.
* The cap invariant is asserted with `expect_gte`, no tolerance, across seeds
  1-5 and `gain_tolerance` in {0, 1e-3, 5e-3, 1e-2, 1}.
* `$power$trial_success` is compared with `expect_identical` against a fresh
  `calc_power_pvals()` on the returned graph.
* `num_threads = 2` gives an identical result to `num_threads = 1` at the same
  seed.

### 4.3 Lint

`lintr::lint_package()`: **229 lints on the branch, all
`object_usage_linter`**; `main` has **184, also all `object_usage_linter`**.
Zero lints of any other class on the branch. That class is "no visible global
function definition" for internal functions referenced across files, which
`codetools` cannot resolve — pre-existing and dominant on `main`.

Five lints that only `lint_package()` surfaces (per-file linting missed them,
and I truncated one output and missed two more) were fixed; all were in the
tests added here.

### 4.4 Documentation and article

`devtools::document()` runs clean. `main` produced no "could not resolve link"
warnings, and my first pass produced ten, all from `[...]` links to internal
`@noRd` topics; they are backticks now, matching the package's existing style.

`man/graph_simplify.Rd` contains the `\subsection{Edge removal during the global
search}` heading and every `sparsity` field.

The article was verified by `knitr::knit()`, which executes every chunk.
Pandoc is not installed in this environment, so **HTML rendering was not
verified** — only chunk execution. The new section runs clean and produces real
output: 6 edges to 5 (free 4 to 3, two being pinned by the constraint), losing
0.05% of a 0.5% budget.

---

## 5. What I could not verify

* **HTML rendering of the article** — no pandoc. Chunks execute; the pkgdown
  build is untested.
* **`R CMD check` with tests.** Run as `--no-tests --no-manual
  --no-build-vignettes`, because the check's test phase is `test_check()`, i.e.
  the whole suite, which the brief reserves for CI. The per-file gates above
  cover the changed code. (The check also has to run in a short output path:
  under the session scratchpad the extracted
  `_snaps/plot_graph_constraint/*.svg` paths exceed Windows `MAX_PATH` and
  extraction fails before checking starts.)
* **Open item O1.** `.simplify_p_zero` is 0.2 with an equal split, the record's
  proposal, implemented on the user's instruction and flagged in the code as
  open. The evidence O1 asks for — the fraction of stage-2 generations in which
  the elite's edge count falls, for `p_zero` in {0.1, 0.2, 0.4} at `m = 4` — was
  **not** gathered. Nothing else depends on the value; re-tuning is a one-line
  change.
* **A12 (pruning cost at scale).** Best-first is O(E) shortcut evaluations per
  removal and at most `E_max` removals, all on the full sample. Not measured at
  `m = 8`, `nsim = 1e6`; that is outside the run limit. On a 1e6-row problem
  this is the dominant cost of stage 2 and deserves a measurement before the
  feature is used in anger.
* **A14 (the Nelder-Mead step inside `GA::ga()`).** Not instrumented. The
  argument that it cannot create an edge inside a feasible support is the same
  one that applies to COBYLA — crossing `1e-5` costs `edge_price` — but I did
  not log zero counts before and after `optim` steps.
* **O4 (halving `run`).** Implemented as specified; no evidence gathered on
  whether half is the right fraction.
* **Runtime claim of §4.2** (1.3 to 1.6 times a single optimisation) — not
  measured.

## 6. Suggested follow-ups, in priority order

1. Decide whether `graph_simplify()` should recompile a dead trial success
   function automatically (§3.1). As it stands the feature does not work on a
   saved graph, which limits it considerably.
2. Fix `trial_success()` serialisation properly, so the whole package survives
   `saveRDS()`. Separate bug fix.
3. Refresh `vignettes/articles/data/*.rds` — pre-0.3.0 objects (§3.5).
4. Close O1 with the measurement it asks for.
5. Correct the record's §2 description of `pvals[sample(nsim), ]` (§3.4).
6. Extract the duplicated `verbose` coercion in `graph_optimise()` and
   `trial_success()` to use `.resolve_verbose()`.
