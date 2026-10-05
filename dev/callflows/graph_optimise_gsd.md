# Call flow: `graph_optimise_gsd()` and its helpers

This document walks through every function involved in optimising a
graph-based multiple testing procedure for a group sequential design, in the
order the code executes them. The functions live in five new files:
`R/check_gsd.R` (shared validation), `R/objective_function_gsd.R` (the
objective closure), `R/optimisation_gsd.R` (the exported entry point and the
two optimiser drivers), `R/calc_power_gsd.R` (power summary and custom
measures), and `R/post_optim_processing_gsd.R` (greedy pruning). Fixed-sample
helpers reused unchanged are `split_theta()`, `recover_full_weights()` and
`recover_full_trans_matrix()` from `R/objective_function.R`;
`param_to_solution()`, `repair_graph()`, `.redistribute_mass()` and
`.marginal_violated()` from `R/post_optim_processing.R`; `choose_graph()`
from `R/choose_graph.R`; `.build_start_matrix()`, `.is_default_start_graph()`,
`.validate_start_graphs()` and `create_start_params()` from
`R/optimisation_start.R` and `R/starting_values.R`; `graph_optimal()` from
`R/graph_optimal.R`; `.cauchy_population` and `.make_cauchy_mutation_multi()`
from `R/optimisation.R` and `R/mutation_helpers.R`; `default_control()`,
`adjust_nsim_local()` and `adjust_nsim_global()` from `R/control_prepare.R`
and `R/control_nsim.R`; `is_graph_valid()` from `R/utils.R`; and
`graph_constraint_get_names()` from the graph-constraint file. The group
sequential kernel, `graph_shortcut_gsd()` and `graph_shortcut_gsd_parallel()`
in `src/graph_shortcut_gsd.cpp`, is described at the point of call but is
documented in `dev/gsd_progress.md` and design record section 4.3.

Notation follows the design record and the earlier call flows: *m* hypotheses,
*K* analyses ("looks"). The `multigrain_pvals_gsd` object holds the
N-by-m-by-K array of repeated (or, per hypothesis, sequential) p-values, plus
the scalars `nsim`, `m`, `K`, `alpha` and the per-hypothesis metadata; it is
produced by `transform_pvalues_gsd()` and documented in the P0 call flow.
The kernel returns a list with `rejected` (N-by-m logical) and `time`
(N-by-m integer, where 0 means "never rejected" and a positive value is the
analysis at which the hypothesis was declared rejected). The gain is a
`multigrain_trial_success_gsd` object whose compiled function takes the
integer time matrix and returns the mean utility; it is produced by
`trial_success_gsd()` and documented in the P3 call flow. The encoded
parameter vector `x` that the optimisers work on is a numeric vector of
length `(number of free hypothesis weights - 1) + (number of free transition
entries - one per row that has at least one free cell)`; it lives in [0, 1]
and is decoded by `split_theta()` and the two `recover_full_*` helpers.


## The call chain at a glance

1. `graph_optimise_gsd()` validates its three main inputs with
   `check_pvals_gsd()`, `check_graph_constraint()` and
   `check_trial_success_gsd()` from `R/check_gsd.R`, resolves the testing
   level with `.gsd_check_alpha()`, asserts dimension agreement between the
   gain and the p-values with `.gsd_check_gain_dims()`, normalises the
   control object with `control_prepare_dims()`, and then follows the same
   global-then-local-then-prune-then-report flow as `graph_optimise()`.

2. The global search (`.graph_optimise_ga_gsd()`) subsamples the p-value
   object with `.sample_pvals_gsd()`, builds the objective closure with
   `create_obj_func_gsd()` (which calls `.gsd_kernel_matrix()` to reshape the
   array), runs `GA::ga()`, decodes the best individual with
   `param_to_solution()`, and, if the decoded graph is valid, re-evaluates it
   on the full sample with `calc_power_pvals_gsd()` and on the full sample
   through the kernel and the gain directly.

3. The local search (`.graph_optimise_local_gsd()`) follows the same
   subsample-then-closure pattern, runs `nloptr::nloptr()`, decodes and
   repairs with `param_to_solution()` then `repair_graph()`, and re-evaluates.

4. `choose_graph()` picks the better of the two solutions (reused unchanged).

5. `prune_graph_gsd()` drives `prune_hyp_weights_gsd()` then
   `prune_edges_gsd()`, each calling `.try_prune_gsd()` which evaluates
   candidates with `calc_power_pvals_gsd()` and accepts only if the gain does
   not decrease and marginal constraints are met (`.marginal_violated()`,
   reused unchanged).

6. `calc_power_pvals_gsd()` runs the kernel once via `.gsd_kernel_matrix()`
   and `graph_shortcut_gsd()`, computes the fixed-sample power metrics from
   the rejection matrix, the group sequential metrics (power by analysis,
   mean decision look, time distribution) from the time matrix, and evaluates
   each `custom_power` entry through `.eval_custom_power_gsd()`, which
   dispatches on the entry's class.

7. `graph_optimal()` wraps the pruned graph, the power summary and the
   optimiser outputs into the returned `multigrain_graph_optimal` object
   (reused unchanged).


---


## `is_pvals_gsd()` -- class predicate

### Purpose

A one-line predicate testing whether its argument is a `multigrain_pvals_gsd`
object. Called by `check_pvals_gsd()` and wherever GSD-specific dispatch is
needed.

### Inputs and output

Takes any R object and returns a scalar logical.

### How it works

Returns `inherits(x, "multigrain_pvals_gsd")`.


---


## `check_pvals_gsd()` -- type check for the p-value object

### Purpose

Validates that an argument is a `multigrain_pvals_gsd` object, or optionally
`NULL`. Called at the top of `graph_optimise_gsd()` and
`calc_power_pvals_gsd()`.

### Inputs and output

- `pvals`: the value to check.
- `arg`: argument name for the error message, defaulting to the caller's
  expression.
- `call`: error call context.
- `allow_null`: whether `NULL` is accepted.

Returns invisible `NULL` on success; aborts otherwise with
`rlang::stop_input_type()`.

### How it works

If `pvals` is not missing and passes `is_pvals_gsd()`, the function returns
immediately. If `allow_null` is `TRUE` and `pvals` is `NULL`, it also returns.
In all other cases it aborts with a message naming the expected type and
suggesting `transform_pvalues_gsd()`.


---


## `check_trial_success_gsd()` -- type check for the gain

### Purpose

Validates that the gain is a `multigrain_trial_success_gsd` object. A
fixed-sample `multigrain_trial_success` object is refused with a specific
message explaining the mismatch: its compiled function expects a logical
rejection matrix, and the GSD consumers only ever have decision times. Called
by `graph_optimise_gsd()`.

### Inputs and output

Takes `trial_success`, `arg` and `call`; returns invisible `NULL` or aborts.

### How it works

If the value passes `is_trial_success_gsd()`, the function returns. If the
value passes `is_trial_success()` (the fixed-sample parent class) but not
`is_trial_success_gsd()`, it aborts with a three-part message: (1) the gain
was created with `trial_success()` and scores rejection indicators only,
(2) a group sequential gain is evaluated on decision times, and (3) the user
should use `trial_success_gsd()`, noting that a rejection-only gain such as
`r1 + r2` can be written there too. In all other cases it aborts with the
generic type message.

### Why it is done this way

Design record section 4.10 `[Rev 2026-09-18]` explains the need: because
`multigrain_trial_success_gsd` inherits from `multigrain_trial_success`, the
fixed-sample `check_trial_success()` would accept a GSD gain, and
`graph_optimise()` would then hand its compiled function a logical rejection
matrix. Rcpp silently coerces `TRUE`/`FALSE` to the integers 1/0, so every
rejection is scored as an analysis-1 rejection, returning a plausible but
wrong number (measured: 0.583 on the true decision times versus 0.667 on the
logical matrix). The dual guard -- here for the gain passed to GSD consumers,
and in `graph_optimise()` and `.auto_name_custom_power()` for the gain passed
to fixed-sample consumers -- closes both directions.


---


## `.gsd_check_gain_dims()` -- dimension agreement between gain and p-values

### Purpose

Asserts that the gain and the p-value object agree on the number of
hypotheses, and, when the gain fixes a number of analyses (because it carries
discount tables), on the number of analyses as well. Called by
`graph_optimise_gsd()` (once, before the first evaluation) and by
`calc_power_pvals_gsd()` (once per compiled gain entry in `custom_power`,
whether GSD or fixed-sample).

### Inputs and output

- `trial_success`: a `multigrain_trial_success_gsd` or
  `multigrain_trial_success` object. Both have an `m` field; only the GSD
  variant has a `K` field.
- `pvals`: a `multigrain_pvals_gsd` object.
- `arg`, `call`: error context.

Returns invisible `NULL` on success.

### How it works

First it checks `trial_success$m != pvals$m`. If they disagree, it aborts
naming both values. Then, only if `trial_success$K` is not `NULL`, it checks
`trial_success$K != pvals$K`. A fixed-sample gain has no `K` field, so the
`is.null(trial_success$K)` test is true and the branch is skipped. If a GSD
gain's K disagrees, the abort message names the number of discount tables and
their lengths, because the K mismatch is most commonly caused by a discount
table of the wrong length.

### Why it is done this way

Design record section 10 item 6 explains that a gain carrying discount tables
indexes a C array with the decision time, so a time above its K reads out of
bounds. The K assertion here, together with the `!anyNA()` assertion in
`.gsd_kernel_matrix()` and the `max(time) <= K` assertion in
`calc_power_pvals_gsd()`, guarantees that the kernel's output is safe to feed
to the compiled function.


---


## `.gsd_check_alpha()` -- resolve the testing level

### Purpose

Resolves the `alpha` argument: `NULL` means the level the boundary tables
were built up to (stored as `pvals$alpha`); an explicit value is validated
against bounds. Called by `graph_optimise_gsd()` and
`calc_power_pvals_gsd()`.

### Inputs and output

- `alpha`: `NULL` or a numeric scalar.
- `pvals`: a `multigrain_pvals_gsd` object.
- `call`: error context.

Returns the resolved numeric alpha.

### How it works

If `alpha` is `NULL`, the function returns `pvals$alpha` immediately.
Otherwise it calls `rlang::check_number_decimal(alpha, min = 0, max =
pvals$alpha)` to check the type and range, then enforces the strict bounds
`alpha > 0` and `alpha <= pvals$alpha` with a custom abort. The abort message
explains that repeated p-values are capped at 1 above `pvals$alpha` (design
record section 4.1 step 4), so no allocation beyond that level could ever
reject.

### Edge cases

A value of exactly 0 is rejected by the strict `alpha > 0` check. A value
equal to `pvals$alpha` passes. A value of 0.05 when p-values were built at
0.025 is rejected. The `rlang::check_number_decimal` call rejects `NA`,
non-numeric and non-scalar inputs before the custom bounds are checked.


---


## `.gsd_kernel_matrix()` -- reshape the p-value array for the kernel

### Purpose

Reshapes the three-dimensional N-by-m-by-K array stored inside a
`multigrain_pvals_gsd` object into the two-dimensional N-by-(m*K) matrix that
the C++ kernel expects, and asserts the absence of `NA` and out-of-range
values. Called by every function that invokes the kernel:
`create_obj_func_gsd()` (once, at closure creation time),
`calc_power_pvals_gsd()` and the re-evaluation blocks of
`.graph_optimise_ga_gsd()` and `.graph_optimise_local_gsd()`.

### Inputs and output

- `pvals`: a `multigrain_pvals_gsd` object.
- `call`: error context.

Returns an N-by-(m*K) numeric matrix.

### How it works

The function extracts `pvals$pvals` (the numeric array), sets `dim(values)
<- c(pvals$nsim, pvals$m * pvals$K)`, and runs two checks. The `dim<-`
assignment does not copy or rearrange data: R's column-major layout means that
column `(k - 1) * m + i` of the reshaped matrix is hypothesis `i` at analysis
`k`, which is exactly the layout the kernel expects (design record section 4.1,
"Output" paragraph).

First, if `anyNA(values)` is true, the function aborts with a message
explaining that the kernel would treat `NA` as never rejectable (because
`NA < x` is false in C++), and suggesting that the user check how the p-value
object was built, since `transform_pvalues_gsd()` never emits `NA`.

Second, if `any(values < 0 | values > 1)`, the function aborts with the
message "The transformed p-value array contains values outside [0, 1]". The
two bullets explain the consequences: a negative value would be read as a
rejection at any positive allocation (because `-0.5 < w * alpha` is true for
any positive allocation), and a value above 1 as never rejectable (because
the kernel compares with strict `<` and the largest possible allocation is
`alpha`). The hint bullet notes that `transform_pvalues_gsd()` emits repeated
p-values in [`1e-14`, 1].

### Why it is done this way

Design record section 10 item 11 specifies a single `anyNA()` assertion on
the public path. The range check is a P4 decision that goes beyond the
record: the comment block in the code explains that the array walk is already
paid for, and that a planted negative (for example `-1`) would silently be
treated as a rejection at any positive allocation, which the adversarial
review confirmed. Both assertions live here so that both consumers
(`create_obj_func_gsd()` and `calc_power_pvals_gsd()`) get them without
repeating them, and the objective closure itself carries no check.


---


## `create_obj_func_gsd()` -- the objective closure

### Purpose

Constructs a closure that the GA and COBYLA optimisers call as
`function(x)`, where `x` is the encoded parameter vector. It is the group
sequential twin of `create_obj_func()` in `R/objective_function.R`. The
only differences from the fixed-sample factory are: the kernel called
(`graph_shortcut_gsd` or `graph_shortcut_gsd_parallel` instead of
`graph_shortcut` or `graph_shortcut_parallel`), the extra `K` argument
passed to both kernel variants, and the matrix handed to the gain
(`res$time` instead of `rej_matrix`). Called by `.graph_optimise_ga_gsd()`
and `.graph_optimise_local_gsd()`.

### Inputs and output

- `m`: number of hypotheses (integer).
- `power_criterion`: the compiled gain function (`trial_success$func`),
  which takes an integer matrix of decision times and returns a scalar.
- `hyp_constraint`: the hypothesis weight constraint vector from the
  graph constraint (NA marks free positions).
- `trans_constraint`: the transition matrix constraint (NA marks free
  cells).
- `pvals`: the reshaped N-by-(m*K) numeric matrix from
  `.gsd_kernel_matrix()`. The caller has already asserted no `NA`.
- `K`: the number of analyses (integer).
- `alpha`: the overall one-sided significance level (default 0.025).
- `num_threads`: number of threads (default 1).

Returns a function of one argument `x`.

### How it works

All captured variables are forced with `force()` to prevent lazy-evaluation
surprises in the closure. The `use_parallel` flag is set once from
`num_threads >= 2L`.

When the returned closure is called with a parameter vector `x`:

1. `split_theta(x, hyp_constraint)` splits the vector into hypothesis-weight
   parameters and transition-matrix parameters. `recover_full_weights()` and
   `recover_full_trans_matrix()` expand them to the full weight vector and
   matrix, using the constraint templates to place fixed values and derive
   the last free entry in each row by complement.

2. Five penalty branches handle degenerate or out-of-range candidates, in the
   same order and with the same values as the fixed-sample closure: if either
   the weights or the matrix contain `NA`, return `-1e6`; if any weight is
   negative, return the sum of negative weights; if any matrix entry is
   negative, return the sum of negative entries; if any weight exceeds 1,
   return minus the sum of the offending weights; if any matrix entry exceeds
   1, return minus the sum. These return values guide the optimiser away from
   infeasible regions without a hard constraint.

3. The snap-to-zero rules zero hypothesis weights below `1e-4` and transition
   entries below `1e-5`. These constants match the fixed-sample closure exactly
   and exist so that the kernel does not waste cascade iterations on allocations
   too small to reject any realistic p-value (the smallest repeated p-value is
   floored at `1e-14`, and `1e-4 * 0.025 = 2.5e-6` is a meaningful allocation,
   but `1e-5 * 0.025 = 2.5e-7` is not).

4. The group sequential kernel is called. If `use_parallel` is true,
   `graph_shortcut_gsd_parallel()` is called with `pvals`, `alpha`, the weight
   vector `w`, the transition matrix `G`, `K`, `num_threads` and `grain_size =
   -1L` (auto-tuned); otherwise `graph_shortcut_gsd()` is called without the
   threading arguments. Both have the C++ signature
   `List graph_shortcut_gsd(NumericMatrix pvals, double alpha, NumericVector w,
   NumericMatrix G, int K)` and return a list with `rejected` (N-by-m logical)
   and `time` (N-by-m integer). The `time` entry records the analysis at which
   each hypothesis was declared rejected, or 0 for never.

5. `power_criterion(res$time)` evaluates the gain on the integer time matrix
   and returns a scalar. This is the single line that differs from the
   fixed-sample closure, which calls `power_criterion(rej_matrix)`.

### Where it differs from `create_obj_func()`

The two closures are identical except for three things. First,
`create_obj_func_gsd()` takes `K` as a parameter, forces it, and passes it to
the kernel. Second, the kernel call is `graph_shortcut_gsd` (or its parallel
variant) instead of `graph_shortcut`. Third, the gain is evaluated on
`res$time` (the integer time matrix) instead of on `rej_matrix` (the logical
rejection matrix). The decode-and-penalise block, the snap-to-zero constants,
the parallel dispatch logic, and the `force()` calls are character-for-character
identical.


---


## `graph_optimise_gsd()` -- exported entry point

### Purpose

The sole exported function for optimising a graph under a group sequential
design. It is the counterpart of `graph_optimise()` for fixed-sample designs.
The user supplies a `multigrain_pvals_gsd` object, a graph constraint and a
`multigrain_trial_success_gsd` gain, and the function returns a
`multigrain_graph_optimal` object whose hypothesis weights and transition
matrix maximise the expected gain. Called by the user; calls every function
documented in this file.

### Inputs and output

- `pvals`: a `multigrain_pvals_gsd` object of repeated (or sequential)
  p-values, as created by `transform_pvalues_gsd()`.
- `graph_constraint`: a `multigrain_graph_constraint` object specifying which
  hypothesis weights and transition matrix entries are fixed and which are
  free.
- `trial_success`: a `multigrain_trial_success_gsd` object defining the gain.
  A fixed-sample `multigrain_trial_success` is refused.
- `...`: enforced empty by `rlang::check_dots_empty()`.
- `alpha`: `NULL` (default) to use the level of `pvals`, or a positive scalar
  at most `pvals$alpha`.
- `start_graph`: a list of starting graphs for the GA; each element has
  `hyp_weight` and `trans_matrix`. `NULL` weights and matrix trigger the
  default Bonferroni-Holm start.
- `global_search`: logical, whether to run the genetic algorithm before the
  local optimiser. Default `TRUE`.
- `num_threads`: number of threads for the C++ kernel. Default 1.
- `control`: a `multigrain_control` object for tuning the GA and COBYLA.
- `verbose`: `"info"`, `"detail"` or `"silent"`.

Returns a `multigrain_graph_optimal` object with `hyp_weight`, `trans_matrix`,
`constraints`, `trial_success`, `power` (the `calc_power_pvals_gsd()` output
for the pruned graph), `solution`, `global_search`, `control`,
`global_output`, `local_output` and `start_graph`.

### How it works

1. The three main inputs are validated: `check_pvals_gsd(pvals)`,
   `check_graph_constraint(graph_constraint)` and
   `check_trial_success_gsd(trial_success)`. The last of these refuses a
   fixed-sample gain with a pointer to `trial_success_gsd()`.

2. The dots, `num_threads`, `control` and `global_search` are validated
   (`rlang::check_dots_empty()`, `rlang::check_number_whole(num_threads, min =
   1)`, `check_control(control)`, `check_logical(global_search, allow_na =
   FALSE)`).

3. The `verbose` argument is normalised: `TRUE` becomes `"info"`, `FALSE`
   becomes `"silent"`, and a string is matched against the package's
   `verbosity_levels`.

4. `.gsd_check_alpha(alpha, pvals)` resolves the testing level.

5. `.gsd_check_gain_dims(trial_success, pvals)` asserts that the gain and
   p-values agree on `m` and, when the gain has a non-`NULL` `K`, on `K`.

6. `.validate_start_graphs(start_graph, m = pvals$m)` checks the starting
   graphs (reused unchanged from the fixed-sample path).

7. `control_prepare_dims(control, nsim = pvals$nsim, m = pvals$m, verbose =
   verbose)` calibrates the control object from the p-value object's
   dimensions.

8. If `global_search` is `TRUE`, `.graph_optimise_ga_gsd()` is called. Its
   best encoded parameter vector is clamped to [0, 1] and saved as
   `x0_for_local`. Otherwise, if the user supplied a non-default start graph,
   `.build_start_matrix()` extracts its first row as `x0_for_local`.

9. `.graph_optimise_local_gsd()` is called, using `x0_for_local` as the
   starting point for COBYLA (or `NULL`, in which case
   `create_start_params()` provides a default).

10. `choose_graph(ga_result, local_result)` picks the solution with the
    higher trial-success value (reused unchanged).

11. `prune_graph_gsd()` greedily removes small weights and edges from the
    best graph, accepting a removal only if the gain does not decrease.

12. `calc_power_pvals_gsd()` evaluates the pruned graph on the full sample
    and reports the power summary, with the gain passed as the sole
    `custom_power` entry.

13. Hypothesis names from `graph_constraint_get_names(graph_constraint)` are
    applied to the weight vector and the transition matrix dimnames.

14. `graph_optimal()` wraps everything into the returned
    `multigrain_graph_optimal` object.

### Where it differs from `graph_optimise()`

The flow is step-for-step the same as the fixed-sample `graph_optimise()`.
The differences are: (a) the input is a `multigrain_pvals_gsd` instead of a
numeric matrix; (b) `check_trial_success_gsd()` replaces
`check_trial_success()`, so a fixed-sample gain is refused; (c) `alpha`
defaults to `NULL` and is resolved by `.gsd_check_alpha()` rather than
defaulting to 0.025; (d) `.gsd_check_gain_dims()` has no fixed-sample
counterpart (the fixed path checks `trial_success$m != ncol(pvals)` inline);
(e) `control_prepare_dims()` replaces `control_prepare()`, taking `nsim` and
`m` as scalars rather than reading `dim(pvals)` from a matrix; (f) every
optimiser driver, pruning function and power function is the `_gsd` variant;
(g) `graph_optimise()` does not call `.gsd_check_alpha()` or
`.gsd_check_gain_dims()`.


---


## `control_prepare_dims()` -- calibrate the control object from dimensions

### Purpose

Group sequential counterpart of `control_prepare()` from
`R/control_prepare.R`. Calibrates the control object (simulation counts, GA
population size, verbosity) from the number of simulated trials and the
number of hypotheses, passed as scalars rather than read from a matrix's
dimensions. Called by `graph_optimise_gsd()`.

### Inputs and output

- `ctrl`: a `multigrain_control` object.
- `nsim`: number of simulated trials (whole number, at least 1).
- `m`: number of hypotheses (whole number, at least 1).
- `verbose`: a verbosity string.
- `call`: error context.

Returns a modified `multigrain_control`.

### How it works

1. Input validation: `check_control(ctrl)`,
   `rlang::check_number_whole(nsim, min = 1)`,
   `rlang::check_number_whole(m, min = 1)`, `rlang::check_string(verbose)`.

2. Both `nsim` and `m` are coerced with `as.integer()`. This matters because
   `control_prepare()` reads `nrow(pvals)` and `ncol(pvals)`, which R stores
   as integers, while the `multigrain_pvals_gsd` object stores `nsim` and `m`
   as the values returned by `dim()` (also integer), but the function's
   signature accepts any whole number. The `as.integer()` ensures the two
   paths agree bit for bit (design record section 6 P4).

3. A default control is built with `default_control()`. Its `nsim_local` is
   set to `nsim`, its `nsim_global` to `min(5e4L, nsim)`, and its
   `global_opt$popSize` to `min(max(40L * m, 200L), 500L)`. This last formula
   means: at least 200 individuals, scaling up with the number of hypotheses
   (40 per hypothesis), but capped at 500.

4. The user's `ctrl` is adjusted: `adjust_nsim_local()` and
   `adjust_nsim_global()` cap user-supplied simulation counts at `nsim` and
   warn if capping occurs.

5. If `verbose` is `"detail"`, the default control's COBYLA print level is set
   to 1 and the GA monitor to `TRUE`.

6. The user's `nsim_local` and `nsim_global` are filled in from the defaults
   via `%||%`.

7. `purrr::list_modify()` merges the user's `global_opt` and `local_opt`
   settings on top of the defaults, so that any user-specified field overrides
   the default while unspecified fields keep their default values.

### Where it differs from `control_prepare()`

The body is almost identical. The only difference is that `control_prepare()`
reads `nrow(pvals)` and `ncol(pvals)` from a matrix argument, while
`control_prepare_dims()` takes `nsim` and `m` as explicit parameters and
coerces them with `as.integer()`. Every default value, every adjustment call,
and every merge is the same.


---


## `.sample_pvals_gsd()` -- subsample the p-value object

### Purpose

Draws a random subsample of simulated trials from a `multigrain_pvals_gsd`
object, preserving its three-dimensional array structure. Called by
`.graph_optimise_ga_gsd()` (for the GA's smaller sample) and
`.graph_optimise_local_gsd()` (for COBYLA's sample).

### Inputs and output

- `pvals`: a `multigrain_pvals_gsd` object.
- `nsim`: the desired subsample size (a positive integer no larger than
  `pvals$nsim`).

Returns a `multigrain_pvals_gsd` object with `pvals$pvals` and `pvals$nsim`
updated.

### How it works

`sample.int(pvals$nsim, size = nsim, replace = FALSE)` draws `nsim` row
indices without replacement. The array is subsetted as `pvals$pvals[idx, , ,
drop = FALSE]`, and `pvals$nsim` is overwritten with `nsim`. The `drop =
FALSE` is critical: without it, a design with `K = 1` would reduce the
three-dimensional array `nsim-by-m-by-1` to a two-dimensional matrix, and the
subsequent `dim<-` reshape in `.gsd_kernel_matrix()` would silently produce
the wrong column layout (design record section 8 item 12).

### Where it differs from the fixed-sample subsampling

The fixed-sample path uses `pvals[sample.int(nrow(pvals), nsim), ]` on a
plain numeric matrix (in `.sample_pvals_rows()` on `main`). The GSD path
subsets a three-dimensional array with `drop = FALSE` and updates the `nsim`
field on the S3 object. The random-draw logic (`sample.int` without
replacement) is the same.


---


## `.graph_optimise_ga_gsd()` -- the genetic algorithm driver

### Purpose

Runs the global optimisation stage using `GA::ga()`. It is the group
sequential twin of `.graph_optimise_ga()` in `R/optimisation.R`. Called by
`graph_optimise_gsd()` when `global_search` is `TRUE`.

### Inputs and output

- `pvals`: the full `multigrain_pvals_gsd` object.
- `graph_constraint`, `trial_success`: as in the parent.
- `nsim`: the GA subsample size (from `control$nsim_global`).
- `global_opts`: the GA options (from `control$global_opt`).
- `alpha`, `num_threads`: passed through.
- `start_graph`: starting graphs for the GA's suggestions.
- `verbose`: verbosity string.

Returns a list with `ga_hyp_weight`, `ga_trans_matrix`,
`ga_trial_success` (scalar gain on the full sample, or `NULL`),
`ga_subset_power` (the `calc_power_pvals_gsd()` result on the subsample, or
`NULL`), `is_graph_valid` (logical) and `ga_output` (the `GA::ga` object).

### How it works

1. A progress step is emitted if not silent.

2. `.build_start_matrix(graph_constraint, start_graph)` encodes each starting
   graph into a row of the suggestion matrix `x0` (reused unchanged).

3. `.sample_pvals_gsd(pvals, nsim)` draws the subsample.

4. The objective closure is built: `create_obj_func_gsd()` receives
   `trial_success$m`, `trial_success$func`, the constraint templates, `alpha`,
   `.gsd_kernel_matrix(pvals_sampled)` (the reshaped subsample), `pvals$K`
   and `num_threads`. Note that `K` comes from the original `pvals` object,
   not from `pvals_sampled`; both have the same `K` since subsampling does not
   change the number of analyses.

5. `GA::ga()` is called via `do.call()` with the immutable arguments (type
   `"real-valued"`, the fitness closure, lower and upper bounds of 0 and 1,
   the Cauchy population initialiser `.cauchy_population`, the Cauchy mutation
   operator from `.make_cauchy_mutation_multi(p_param_mutate = 0.1, scale =
   1.0)`, `optim = TRUE` for intermittent Nelder-Mead, and the suggestions
   `x0`), merged with `global_opts` via `modifyList` so the user's settings
   override.

6. The best individual `ga_res@solution[1, ]` is decoded by
   `param_to_solution(best_raw, graph_constraint, process = TRUE)`, which
   snaps small weights to 0, epsilon values to 0.001, and normalises rows.

7. `is_graph_valid(sol$hyp_weight, sol$trans_matrix)` checks the decoded
   graph.

8. If valid, two re-evaluations happen on the full sample.
   `calc_power_pvals_gsd(pvals_sampled, ...)` evaluates on the subsample with
   the gain as `custom_power`. Then the kernel is called directly on the full
   sample via `graph_shortcut_gsd(pvals = .gsd_kernel_matrix(pvals), alpha =
   alpha, w = sol$hyp_weight, G = sol$trans_matrix, K = pvals$K)`, and
   `trial_success$func(time_mat)` evaluates the gain on the full sample's
   time matrix. The full-sample gain is stored as `ga_trial_success`.

### Where it differs from `.graph_optimise_ga()`

Three differences. First, `.sample_pvals_gsd()` replaces the fixed-sample
row-subsetting. Second, `create_obj_func_gsd()` replaces `create_obj_func()`.
Third, the re-evaluation uses `graph_shortcut_gsd()` and evaluates
`trial_success$func(time_mat)` instead of calling `graph_shortcut()` and
passing the rejection matrix. The starting-value logic, the GA arguments,
`param_to_solution()`, `is_graph_valid()` and the overall control flow are
identical.


---


## `.graph_optimise_local_gsd()` -- the COBYLA driver

### Purpose

Runs the local optimisation stage using `nloptr::nloptr()` (COBYLA by
default). It is the group sequential twin of `.graph_optimise_local()`.
Called by `graph_optimise_gsd()`.

### Inputs and output

- `pvals`, `graph_constraint`, `trial_success`: as above.
- `local_opts`: the COBYLA options (from `control$local_opt`).
- `alpha`, `num_threads`, `nsim`: as above.
- `x0`: the starting point, typically the GA's best individual clamped to
  [0, 1], or `NULL`.
- `verbose`: verbosity string.

Returns a list with `local_hyp_weight`, `local_trans_matrix`,
`local_trial_success`, `local_subset_power`, `is_graph_valid` and
`local_output` (the `nloptr` result).

### How it works

1. If `x0` is `NULL`, `create_start_params(graph_constraint)` provides a
   default starting point (reused unchanged).

2. `.sample_pvals_gsd(pvals, nsim)` draws the subsample for COBYLA.

3. `create_obj_func_gsd()` builds the closure on the subsampled, reshaped
   matrix.

4. Because `nloptr` minimises by convention, the closure is wrapped:
   `nlopt_obj_func <- function(x) -obj_fun(x)`.

5. `nloptr::nloptr()` is called with `x0`, the negated objective, lower and
   upper bounds of 0 and 1, and `local_opts`.

6. `param_to_solution(nlopt_result$solution, graph_constraint, process =
   TRUE)` decodes the solution. Then `repair_graph(sol$hyp_weight,
   sol$trans_matrix, graph_constraint)` projects it onto the feasible region
   (clamping to [0, 1], zeroing the diagonal, pinning fixed elements, and
   normalising). `repair_graph()` is reused unchanged from the fixed-sample
   path and is called here but not in the GA driver, because COBYLA's
   tolerance-based termination can leave tiny constraint violations.

7. `is_graph_valid()` checks the decoded and repaired graph.

8. If valid, re-evaluation proceeds identically to the GA driver: subsample
   power via `calc_power_pvals_gsd(pvals_sampled, ...)`, then a direct kernel
   call on the full sample followed by `trial_success$func(time_mat)`.

### Where it differs from `.graph_optimise_local()`

The same three differences as the GA driver: `_gsd` subsampling, `_gsd`
closure, and `_gsd` kernel in the re-evaluation. The `repair_graph()` call,
the COBYLA setup and the flow are identical.


---


## `calc_power_pvals_gsd()` -- power summary for a group sequential graph

### Purpose

Evaluates a graph on a `multigrain_pvals_gsd` object and reports per-hypothesis
power, per-analysis power, the decision-time distribution, and any
user-supplied custom power measures. It is the group sequential counterpart of
`calc_power_pvals()`. Called by `graph_optimise_gsd()` (for the final pruned
graph), by `.graph_optimise_ga_gsd()` and `.graph_optimise_local_gsd()` (for
the re-evaluation on the subsample), by `.try_prune_gsd()` (during greedy
pruning), and directly by the user.

### Inputs and output

- `pvals`: a `multigrain_pvals_gsd` object.
- `hyp_weight`: numeric vector of hypothesis weights (length m).
- `trans_matrix`: numeric m-by-m transition matrix.
- `...`: enforced empty.
- `alpha`: `NULL` or a scalar; resolved by `.gsd_check_alpha()`.
- `custom_power`: `NULL`, a single measure, or a list of measures. Three
  kinds are accepted: a `multigrain_trial_success_gsd`, a fixed-sample
  `multigrain_trial_success`, or a plain R function.
- `sum_to_one_constraint`: logical, whether transition rows must sum to one.
- `call`: error context.

Returns a list with:
- `local_power`: the proportion of simulations in which each hypothesis is
  rejected at any analysis (length-m numeric).
- `local_power_by_analysis`: an m-by-K matrix whose entry (i, k) is the
  proportion of simulations in which hypothesis i was rejected at or before
  analysis k. Its last column equals `local_power`.
- `exp_rejections`: the expected number of rejections.
- `disj_power`: probability of rejecting at least one hypothesis.
- `conj_power`: probability of rejecting all hypotheses.
- `mean_decision_look`: for each hypothesis, the mean analysis index over the
  simulations in which it was rejected; `NA` when it was never rejected.
- `time_distribution`: an m-by-(K+1) matrix of proportions, with a `"never"`
  column (the proportion of simulations with decision time 0) followed by one
  column per analysis. Rows sum to one, and one minus the `"never"` column
  equals `local_power` up to floating-point rounding.
- One entry per element of `custom_power`, named by the list name.

### How it works

1. `check_pvals_gsd(pvals)` validates the p-value object.
   `check_double(hyp_weight)` and `check_double_matrix(trans_matrix)` validate
   the graph. `rlang::check_dots_empty()` enforces no extra arguments.

2. `.gsd_check_alpha(alpha, pvals)` resolves the level.

3. `is_graph_valid(hyp_weight, trans_matrix, sum_to_one_constraint)` checks
   the graph; the function aborts if it is invalid.

4. `.auto_name_custom_power_gsd(custom_power)` normalises `custom_power` to
   a named list. Then, for each element that passes `is_trial_success()` --
   which covers both GSD and fixed-sample compiled gains, since the GSD class
   inherits from the fixed-sample class -- `.gsd_check_gain_dims()` asserts
   that the gain's `m` matches `pvals$m` and, when the gain has a non-`NULL`
   `K`, that `K` matches `pvals$K`. A fixed-sample gain carries no `K` field,
   so the `K` branch is skipped naturally for it. The `m` check is essential
   for both classes: a fixed-sample `trial_success(r1 + r2 + r3)` passed
   to a two-hypothesis p-value object would be handed a two-column `rejected`
   matrix, and the compiled function would read column index 2 (zero-based)
   past the end of each row, returning garbage without an abort.

5. The kernel is called: `graph_shortcut_gsd(pvals =
   .gsd_kernel_matrix(pvals, call = call), alpha = alpha, w = hyp_weight, G =
   trans_matrix, K = pvals$K)`. The serial kernel is always used here (no
   parallel dispatch), because this is the reporting path, not the hot loop.

6. A defensive assertion checks `max(res$time) > pvals$K`. This cannot
   happen with the kernel of design record section 4.3, but a decision time
   above K would be an out-of-bounds read in a gain that carries discount
   tables. If triggered, the function aborts.

7. The fixed-sample metrics are computed from `res$rejected`:
   `colMeans(res$rejected)` for `local_power`, `sum(res$rejected) /
   nrow(res$rejected)` for `exp_rejections`, and the row-sum-based
   disjunctive and conjunctive powers.

8. The group sequential metrics are computed from `res$time`:
   `.gsd_power_by_analysis()`, `.gsd_mean_decision_look()` and
   `.gsd_time_distribution()`.

9. `.eval_custom_power_gsd(custom_power, res$rejected, res$time, call)` 
   evaluates each custom measure.

10. All results are assembled by `c(list(...), .eval_custom_power_gsd(...))`.

### Where it differs from `calc_power_pvals()`

`calc_power_pvals()` calls `graph_shortcut()` (the fixed-sample kernel),
returns only the four fixed-sample metrics, and evaluates custom power via
`.eval_custom_power()` which passes only the rejection matrix.
`calc_power_pvals_gsd()` calls `graph_shortcut_gsd()`, adds the five
group-sequential-specific outputs (`local_power_by_analysis`,
`mean_decision_look`, `time_distribution`, and implicitly the custom measure
dispatch), asserts `max(time) <= K`, and evaluates custom power via
`.eval_custom_power_gsd()` which dispatches on class.


---


## `.gsd_power_by_analysis()` -- cumulative power by analysis

### Purpose

Computes the m-by-K matrix of cumulative local power: entry (i, k) is the
proportion of simulations in which hypothesis i was rejected at or before
analysis k. Called by `calc_power_pvals_gsd()`.

### Inputs and output

- `time`: the N-by-m integer matrix of decision times from the kernel.
- `m`: number of hypotheses.
- `n_look`: number of analyses (K).

Returns an m-by-K numeric matrix with dimnames `"H1", "H2", ...` on rows and
`"analysis 1", "analysis 2", ...` on columns.

### How it works

The output matrix is initialised to `NA_real_`. For each hypothesis `i` and
each analysis `k`, the entry is `mean(tau >= 1L & tau <= k)`, where `tau` is
the i-th column of `time`. The condition `tau >= 1L` excludes
never-rejected trials (where `tau` is 0), and `tau <= k` selects trials
rejected at or before analysis k. This makes the last column equal to the
proportion rejected at any analysis, which is `local_power`.


---


## `.gsd_mean_decision_look()` -- mean decision analysis index

### Purpose

For each hypothesis, computes the mean analysis index over the simulations in
which it was rejected, or `NA` when it was never rejected. Called by
`calc_power_pvals_gsd()`.

### Inputs and output

- `time`: the N-by-m integer matrix of decision times.

Returns a numeric vector of length m.

### How it works

Uses `vapply` over columns. For column `i`, `tau <- time[, i]` extracts the
decision times. `rejected <- tau > 0L` identifies trials where the hypothesis
was rejected. If no trial rejected it, `NA_real_` is returned. Otherwise,
`mean(tau[rejected])` gives the average analysis index among the rejected
trials.

### Why it is done this way

Design record section 4.4 `[Rev 2026-09-18]` explains that this quantity is
`mean_decision_look` rather than a calendar-time mean: averaging look indices
and converting afterwards would be wrong for unequally spaced analyses.
Calendar-time metadata is not stored (section 10 item 12, deferred by user
decision).


---


## `.gsd_time_distribution()` -- proportions of decision times

### Purpose

Computes the m-by-(K+1) matrix of decision-time proportions. Called by
`calc_power_pvals_gsd()`.

### Inputs and output

- `time`: the N-by-m integer matrix of decision times.
- `m`: number of hypotheses.
- `n_look`: number of analyses (K).

Returns an m-by-(K+1) numeric matrix with dimnames `"H1", "H2", ...` on rows
and `"never", "analysis 1", ..., "analysis K"` on columns.

### How it works

The output matrix is initialised to `NA_real_`. For each hypothesis `i` and
each decision time `k` from 0 to K, the entry at column `k + 1` is
`mean(tau == k)`. This gives the proportion of simulations with that exact
decision time. The "never" column (k = 0) equals one minus `local_power` up
to floating-point rounding (both are computed as `mean()` of different boolean
conditions on the same integer column), and the remaining columns partition
the rejections across the K analyses. Rows sum to one because every simulation
falls into exactly one category.


---


## `.auto_name_custom_power_gsd()` -- normalise the custom_power argument

### Purpose

Normalises the `custom_power` argument into a named list, validating each
element and assigning default names to unnamed positions. Called by
`calc_power_pvals_gsd()`.

### Inputs and output

- `x`: `NULL`, a single function or trial-success object, or a list of them.
- `call`: error context.

Returns a named list (possibly empty).

### How it works

If `x` is `NULL`, an empty list is returned. If `x` is a single function or
passes `is_trial_success()` (which includes both fixed-sample and GSD gains,
since the latter inherits from the former), it is wrapped as
`list(custom_power = x)`. If `x` is not a list, the function aborts. Each
list element is checked: it must be a function, a
`multigrain_trial_success_gsd` or a `multigrain_trial_success`. Unnamed or
blank-named positions receive names `"func1"`, `"func2"`, etc.

### Where it differs from `.auto_name_custom_power()`

The fixed-sample `.auto_name_custom_power()` additionally checks for
`is_trial_success_gsd()` on each element and aborts if found (the
fixed-consumer guard of design record section 4.10). The GSD version does
not need that guard because a GSD gain is perfectly valid in a GSD consumer.
It does, however, accept fixed-sample `multigrain_trial_success` objects
alongside GSD ones, because `.eval_custom_power_gsd()` dispatches on class.


---


## `.eval_custom_power_gsd()` -- dispatch custom measures by class

### Purpose

Evaluates each entry in a normalised `custom_power` list, choosing the right
matrix to pass based on the entry's class. Called by
`calc_power_pvals_gsd()`.

### Inputs and output

- `custom_power`: a named list (from `.auto_name_custom_power_gsd()`).
- `rejected`: the N-by-m logical matrix from the kernel.
- `time`: the N-by-m integer matrix from the kernel.
- `call`: error context.

Returns a named list of scalars, one per entry.

### How it works

The function applies `lapply` over the list. For each item, it dispatches on
three cases, in this order:

1. `is_trial_success_gsd(item)`: evaluated as `item$func(time)`. The compiled
   function takes the integer time matrix and returns the mean utility.

2. `is_trial_success(item)` (but not GSD, since the GSD check was first):
   evaluated as `item$func(rejected)`. The fixed-sample compiled function
   takes the logical rejection matrix.

3. `is.function(item)`: a plain R function, called once per trial on the
   trial's integer decision-time row. The results are averaged with `mean(
   vapply(..., numeric(1L)))`. This is the slow path, as each trial is
   evaluated separately in R.

4. Anything else aborts.

### Why it is done this way

Design record section 4.7 `[Rev 2026-09-18]` specifies that dispatch on the
gain type is explicit (Sol item 5). A `multigrain_trial_success_gsd` entry
gets the time matrix (so discount tables work); a fixed-sample entry gets the
rejection matrix (so a rejection-only gain reports the same number as
`calc_power_pvals()` would); and a plain function gets the time row (since
time is the richer output). The `is_trial_success_gsd` check must come before
`is_trial_success` because the GSD class inherits from the fixed-sample
class.

### Where it differs from `.eval_custom_power()`

The fixed-sample `.eval_custom_power()` has only two branches: compiled
(`is_trial_success(item)`, passing the rejection matrix) and plain function
(passing the rejection row). The GSD version adds a third branch for the GSD
gain class, and both the compiled and plain-function branches receive the
time matrix (the compiled GSD gain) or the time row (the plain function)
rather than the rejection matrix or row.


---


## `prune_graph_gsd()` -- greedy pruning orchestrator

### Purpose

Removes small hypothesis weights and transition-matrix edges from the
optimised graph, accepting each removal only if the gain does not decrease.
It is the group sequential twin of `prune_graph()` in
`R/post_optim_processing.R`. Called by `graph_optimise_gsd()`.

### Inputs and output

- `pvals`: the full `multigrain_pvals_gsd` object.
- `hyp_weight`, `trans_matrix`: the optimised graph.
- `trial_success`: the gain.
- `graph_constraint`: the constraint object.
- `alpha`: the testing level (default 0.025).
- `gamma`: the threshold below which a weight or edge is considered for
  removal (default 1, meaning every non-zero entry is a candidate).
- `power_constraint`: optional marginal-power thresholds.
- `verbose`: verbosity string.

Returns a list with `hyp_weight` and `trans_matrix`.

### How it works

1. `verbose` is matched and a progress step is emitted if not silent.

2. Constraint metadata is derived once: `fixed_w` (which weights are
   constrained), `fixed_edge` (which transition entries are constrained),
   and `tolerance` (the sum-to-one tolerance from the graph constraint).

3. `prune_hyp_weights_gsd()` is called with the graph, the p-values, the
   gain, and the constraint metadata. It returns a possibly-pruned weight
   vector, the unchanged transition matrix, and the current best gain.

4. `prune_edges_gsd()` is called with the (possibly pruned) weights, the
   transition matrix, and the best gain from the weight-pruning stage. It
   returns the unchanged weights, the possibly-pruned matrix, and the final
   best gain.

### Where it differs from `prune_graph()`

The body is identical except that `prune_hyp_weights_gsd` and
`prune_edges_gsd` are called instead of `prune_hyp_weights` and
`prune_edges`. The constraint metadata extraction, the two-stage flow, and
the return value are the same.


---


## `prune_hyp_weights_gsd()` -- prune small hypothesis weights

### Purpose

Attempts to zero out each small hypothesis weight by redistributing its mass,
accepting the removal if the gain does not decrease and marginal constraints
are met. Called by `prune_graph_gsd()`.

### Inputs and output

- `pvals`, `hyp_weight`, `trans_matrix`, `trial_success`, `fixed_w`,
  `alpha`, `gamma`, `power_constraint`, `tolerance`: as described in
  `prune_graph_gsd()`.

Returns a list with `hyp_weight` (pruned), `trans_matrix` (unchanged) and
`power_best`.

### How it works

1. `calc_power_pvals_gsd()` evaluates the current graph and extracts
   `power_all$custom_power` as the baseline `power_best`.

2. `constrained_idx` is set to the indices of non-`NA` entries in
   `power_constraint`, or an empty integer vector if there are no marginal
   constraints.

3. A loop iterates from hypothesis `m` down to 1. For each hypothesis `i`,
   three skip conditions are checked: the weight is already 0 or at least
   `gamma`; the weight is fixed by the constraint; or the hypothesis has a
   marginal constraint but is not among the constrained indices.

4. `.redistribute_mass(hyp_weight, drop_idx = i, fixed_idx = fixed_w,
   tolerance = tolerance)` zeros the weight and distributes its mass
   proportionally among the free recipients (or uniformly if all free
   recipients are zero). This function is reused unchanged from the
   fixed-sample path.

5. `.try_prune_gsd()` evaluates the candidate graph. If accepted (gain did
   not decrease, no marginal violation), `hyp_weight` and `power_best` are
   updated.

### Where it differs from `prune_hyp_weights()`

The only difference is that `.try_prune_gsd()` is called instead of
`.try_prune()`. The loop structure, the skip conditions, the redistribution
logic, and the acceptance rule are identical.


---


## `prune_edges_gsd()` -- prune small transition edges

### Purpose

Attempts to zero out each small transition-matrix edge by redistributing its
mass within its row, accepting the removal if the gain does not decrease.
Called by `prune_graph_gsd()`.

### Inputs and output

- `pvals`, `hyp_weight`, `trans_matrix`, `trial_success`, `fixed_edge`,
  `power_best`, `alpha`, `gamma`, `power_constraint`, `tolerance`: as above.

Returns a list with `hyp_weight` (unchanged), `trans_matrix` (pruned) and
`power_best`.

### How it works

1. `G_best` is initialised to the input transition matrix.

2. A double loop iterates: the outer loop runs target column `j` from `m`
   down to 1, and the inner loop runs source row `i` from 1 to `m`. This
   means edges into the same hypothesis are considered together. For each
   cell `(i, j)`:
   - `i == j` is skipped (the diagonal is always 0).
   - The edge is skipped if it is already 0 or at least `gamma`.
   - The edge is skipped if it is fixed by the constraint.

3. The fixed indices for row `i` are `which(fixed_edge[i, ])`, which includes
   the diagonal (always fixed in the constraint template). A candidate
   matrix `G_candidate` is built by calling `.redistribute_mass()` on
   `G_best[i, ]` with `drop_idx = j`, distributing the mass among the free
   cells in that row.

4. `.try_prune_gsd()` evaluates the candidate. If accepted, `G_best` and
   `power_best` are updated.

### Where it differs from `prune_edges()`

The only difference is `.try_prune_gsd()` instead of `.try_prune()`. The
double loop, the skip conditions, the per-row redistribution and the
acceptance rule are identical.


---


## `.try_prune_gsd()` -- evaluate a pruning candidate

### Purpose

Evaluates a candidate graph and accepts it if the gain does not decrease and
marginal constraints are satisfied. Called by `prune_hyp_weights_gsd()` and
`prune_edges_gsd()`.

### Inputs and output

- `pvals`: the full `multigrain_pvals_gsd` object.
- `hyp_weight`, `trans_matrix`: the candidate graph.
- `trial_success`: the gain.
- `power_best`: the current best gain (scalar).
- `constrained_idx`, `power_constraint`: marginal constraints.
- `alpha`: the testing level (default 0.025).

Returns a list with `hyp_weight`, `trans_matrix`, `power_best` and `accepted`
(logical).

### How it works

1. `calc_power_pvals_gsd()` evaluates the candidate graph with the gain
   passed as `custom_power = trial_success`.

2. `.marginal_violated()` checks whether any constrained hypothesis's
   `local_power` falls below its threshold (reused unchanged).

3. A default output is built with `accepted = FALSE`. If the marginal
   constraints are satisfied and `power_all$custom_power >= power_best`, the
   output's `power_best` is updated and `accepted` is set to `TRUE`.

### Where it differs from `.try_prune()`

The only difference is `calc_power_pvals_gsd()` instead of
`calc_power_pvals()`. The acceptance rule (`!violated && gain >= best`) and
the return structure are identical.


---


## The guard in the fixed-sample consumers

Two hunks were added to existing files to prevent a
`multigrain_trial_success_gsd` gain from being silently misused in the
fixed-sample pipeline.

### `graph_optimise()` in `R/optimisation.R`

After the existing `check_trial_success(trial_success)` call (which accepts
the GSD gain because it inherits from `multigrain_trial_success`), a new
block checks `is_trial_success_gsd(trial_success)`. If true, it aborts with
a three-part message: (1) the gain was created with `trial_success_gsd()` and
needs decision times, (2) `graph_optimise()` supplies only fixed-sample
rejection indicators so every rejection would be scored as an analysis-1
rejection, and (3) the user should use `graph_optimise_gsd()` with
`transform_pvalues_gsd()`.

### `.auto_name_custom_power()` in `R/calc_power.R`

A local helper `abort_gsd_gain(what)` is defined, which aborts with the same
three-part message structure, parameterised by `what` (a cli-formatted label
for the offending argument). The function checks the input at two points:
when `x` is a single object and passes `is_trial_success_gsd()`, and when
iterating over list elements, for each element that passes
`is_trial_success_gsd()`. In the list case, the label names the element
(by its list name if present, or by position).

### Why these guards exist

Design record section 4.10 `[Rev 2026-09-18]` explains: because the GSD
class inherits from the fixed-sample class, `check_trial_success()` accepts
a GSD gain. `graph_optimise()` would then pass the logical rejection matrix
to the compiled function, which Rcpp coerces to integers 1/0. The compiled
function's `double(t(i, idx) > 0)` sees every TRUE as 1 (analysis 1), so
every rejection is scored as an analysis-1 rejection. For a gain with
discount tables, this returns a plausible but wrong value (measured: 0.667
instead of 0.583). The guards add no code path for fixed-sample gains and
do not affect the identity gate; they are kept in their own commit so the
sparsity branch can rebase past them.


---


## The discount-table warning in `R/trial_success_gsd.R`

Design record section 6 P3 `[Rev 2026-09-18]` (Sol item 9) specifies that
`.gsd_gain_tables()` calls `.gsd_gain_warn_tables()` after validating the
discount tables. For each table, `.gsd_gain_table_bullets()` builds a
character vector of warning bullets. A local `tol <- sqrt(.Machine$double.eps)`
is computed once and shared by both checks:

- Values outside [0, 1] beyond the tolerance (the condition is `tab < -tol |
  tab > 1 + tol`): the offending values are listed, with a note that the
  manuscript defines discount multipliers on [0, 1]. Because of the tolerance,
  a value such as `1 + 1e-12` does not trigger the warning.
- Increases from one analysis to a later one (the condition is `diff(tab) >
  tol`): each offending pair of positions and values is listed, with a note
  that a discount table normally assigns no more value to a later decision. If
  the first value is 0, an extra bullet warns that `trial_success_gsd()`
  supplies `d(0) = 0` for "never rejected" automatically, since
  `c(0, 1, 0.75)` is the common shifted-by-one mistake.

Both checks use the same `sqrt(.Machine$double.eps)` tolerance so that a
table built by arithmetic (for example a sequence computed with `cumprod()`)
is not flagged for rounding noise. These are warnings, not errors: the
package does not forbid non-standard value tables. The rest of
`R/trial_success_gsd.R` is documented in the P3 call flow.
