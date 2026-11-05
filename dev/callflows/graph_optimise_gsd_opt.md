# Call flow: `graph_optimise_gsd()` on a spending spec and its helpers

This document walks through every function involved in optimising
Hwang–Shih–DeCani (HSD) spending parameters jointly with the graph of a group
sequential design, in the order the code executes them (design record section
4.12; design `dev/hsd_opt_design_revised.md`; build plan
`dev/hsd_opt_build_plan.md`). The functions live in three new files:
`R/transform_pvalues_gsd_opt.R` (the tag `sf_hsd_opt()`, the spec builder and
its tables, `resolve_pvalues_gsd()`, labels and the spec's `print()`),
`R/objective_function_gsd_opt.R` (the slice bank and the objective closure of
the global stage) and `R/optimisation_gsd_opt.R` (the driver, the two stages,
the profile helpers and the display helpers).

Four existing files gain guarded lines, each behind a predicate that is `FALSE`
on the fixed path (record 4.10 `[Rev 2026-10-14]`):

- `R/transform_pvalues_gsd.R`: the handoff to the spec builder in
  `transform_pvalues_gsd()` and an internal abort in `.gsd_spend()`;
- `R/check_gsd.R`: the spec hint in `check_pvals_gsd()`;
- `R/optimisation_gsd.R`: the dispatch as the first statement of
  `graph_optimise_gsd()`;
- `R/graph_optimal.R`: the `spending = NULL` argument of `graph_optimal()` and
  one guarded display call each in `print()` and `summary()`.

Fixed helpers reused unchanged:

- from `R/transform_pvalues_gsd.R`: `.gsd_boundary_table()` (one table per
  node, with its well-ordering check), `.gsd_invert()`, `.gsd_looks()`,
  `.gsd_check_spending()`, `.gsd_transform_hyp()` (fixed and short-circuited
  hypotheses), `.gsd_warn_matured()`, `.gsd_spending_labels()` and
  `new_pvals_gsd()`;
- from `R/check_gsd.R` and `R/stopping_rule_gsd.R`: `.gsd_check_alpha()`,
  `.gsd_check_graph_dims()`, `.gsd_check_gain_dims()`,
  `.gsd_check_stopping_rule_dims()`, `.gsd_stop_table()` and
  `.gsd_kernel_matrix()`;
- `create_obj_func_gsd()`, the closure every Stage L COBYLA run optimises;
- `control_prepare_dims()`, `.build_start_matrix()`, `.validate_start_graphs()`,
  `create_start_params()`, `closest_graph_to_constraints()`,
  `param_to_solution()`, `repair_graph()`, `split_theta()`,
  `recover_full_weights()`, `recover_full_trans_matrix()`,
  `.cauchy_population` and `.make_cauchy_mutation_multi()`;
- `choose_graph()`, `prune_graph_gsd()`, `calc_power_pvals_gsd()`,
  `.calc_power_pvals_gsd_core()`, `is_graph_valid()` and `graph_optimal()`.

The kernel, `graph_shortcut_gsd()` and `graph_shortcut_gsd_parallel()` in
`src/graph_shortcut_gsd.cpp`, is unchanged and is described in
`dev/callflows/graph_optimise_gsd.md` and design record section 4.3.

Notation follows the design record and the fixed call flows: *m* hypotheses,
*K* analyses, *N* (`nsim`) simulated trials, *d* the number of optimised
parameters. Parameter *r* has a lattice Γ_r of *M_r* nodes: the multiples of
`step` inside `bounds`, plus both bounds (21 nodes for the default
`bounds = c(-4, 1)`, `step = 0.25`). The spec (`multigrain_pvals_gsd_spec`)
holds:

- `raw`, the N-by-m-by-K raw p-values;
- `nsim`, `m`, `K`, `alpha`, `info_frac`, `look_back` and `grid_size`, as on
  a `multigrain_pvals_gsd`;
- `fixed` and `fixed_tables`, the transform of each fixed or short-circuited
  hypothesis (`NULL` for an active one);
- `param`, a list of `names`, `bounds` (*d*-by-2), `step`, `start` (lattice
  values) and `lattice` (a named list of node vectors);
- `group`, the parameter of each hypothesis, an integer in 1..*d*, or `NA`
  for a fixed or short-circuited one;
- `schedule_key`, the information-schedule key of each active hypothesis;
- `node_tables`, keyed `"<schedule>|<node>"`;
- `spending`, the labels; `short_circuit`; and `spending_label`.

The gene vector of the global stage is
`x = c(graph genes [p_graph], u_1, ..., u_d)`, where
`p_graph = ncol(.build_start_matrix(graph_constraint, start_graph))` is the
fixed path's encoding and each `u_r` lies in [0, 1]. Stage L's COBYLA sees
only the graph genes.


## The call chain at a glance

1. **Handoff to the spec builder.** `transform_pvalues_gsd()` validates its
   inputs as before, records `spending_recycled <- is.function(spending)` and
   builds the list with `.gsd_spending_list()`. If `.has_hsd_opt(spending)`
   holds, it returns `.transform_pvalues_gsd_spec()`, which:
   - calls `.hsd_groups()` and `.hsd_lattice()`;
   - transforms the fixed and short-circuited hypotheses with
     `.gsd_transform_hyp()`;
   - builds every node table with `.hsd_node_tables()` (which calls
     `.hsd_check_invertible()`);
   - labels the hypotheses with `.hsd_labels()`;
   - wraps everything with `new_pvals_gsd_spec()`.

   With no parameter left after short-circuits, it returns a plain
   `multigrain_pvals_gsd` instead.

2. **Dispatch.** The first statement of `graph_optimise_gsd()` sends a spec
   to `.graph_optimise_gsd_opt()` after `rlang::check_dots_empty()`.
   `check_pvals_gsd()` never sees the spec on this path; anywhere else it
   refuses it with a hint.

3. **Validation and control.** `.graph_optimise_gsd_opt()` runs the fixed
   path's checks against the spec's `nsim`, `m`, `K` and `alpha`, builds the
   stop table and calls `control_prepare_dims()`. It then sets `popSize` to
   `min(max(40 (m + d), 200), 500)` unless the user set it.

4. **Stage G.** `.graph_optimise_ga_gsd_opt()` builds the start matrix, draws
   `idx_ga`, builds the slice bank with `.hsd_bank()` and the closure with
   `create_obj_func_gsd_opt()` (using `.hsd_node_index()`), and runs
   `GA::ga()`. It decodes the winner, re-scores it on the full sample at its
   own nodes (via `resolve_pvalues_gsd()`), and returns the archive.

5. **Stage L step 1.** `.graph_optimise_local_gsd_opt()` draws `idx_local`
   once and runs `.hsd_refine_node()` (COBYLA through
   `create_obj_func_gsd()`) at the GA's parameter values. The start is the
   GA's genes, or, when `global_search = FALSE`, the user's start graph
   (`.hsd_start_genes()`, `.hsd_encode_graph()`) at the start node.

6. **Profile sweep and rounds.** `.hsd_profile()` calls
   `.hsd_profile_axis()` for each axis. That function walks the axis's nodes
   outward and calls `.hsd_profile_node()` at each: three candidate graphs,
   then `.hsd_refine_node()` from the best. The optimum moves when an axis
   finds a better gain. For *d* ≥ 2 the sweep repeats in rounds (at most 3).
   Stale axes are recomputed once against the final optimum. Paired
   differences use `.hsd_block_gains()` and `.hsd_block_diff()`, and the
   table is built by `.hsd_profile_table()`.

7. **`choose_graph()`** (unchanged) keeps the better of the two stages'
   full-sample gains. The spending values of the stage chosen ride alongside.

8. **Resolve, prune, power, `graph_optimal()`.** `resolve_pvalues_gsd(spec,
   theta_star)` is passed to the unchanged `prune_graph_gsd()` and
   `calc_power_pvals_gsd()` with the user's stopping rule. `graph_optimal()`
   is called with `spending = list(...)`, which carries the nine spending
   fields; `spending_search` comes from `.hsd_search_table()`.

9. **Display.** `print()` of the result calls `.hsd_cat_print()` (with
   `.hsd_spending_notes()`, `.hsd_infer_groups()`, `.hsd_fmt_lattice()`,
   `.hsd_plateau()`, `.hsd_fmt_interval()`), and `summary()` calls
   `.hsd_cat_summary()`. Both are guarded by `!is.null(x$spending_param)`.


---


## `sf_hsd_opt()` -- the tagged spending function

### Purpose

Marks a hypothesis whose HSD parameter is to be chosen by the optimiser. It
is the only new user-facing input; everything else follows from finding one
in `spending`.

### Inputs and output

- `id`: `NULL`, or one non-missing string or number; normalised with
  `as.character()`. It names the parameter and groups hypotheses.
- `bounds`: finite numeric of length 2, lower ≥ −10, upper ≤ 2, lower <
  upper (default `c(-4, 1)`).
- `step`: one of `0.5`, `0.25`, `0.125` (default `0.25`).
- `start`: `NULL` or one finite number inside `bounds`.

Returns `function(alpha, t, param) gsDesign::sfHSD(alpha, t, param)$spend`,
of class `c("multigrain_sf_hsd_opt", "function")`, with attributes `id`,
`bounds`, `step` and `start` (normalised values).

### How it works

1. `.hsd_check_id()`, `.hsd_check_bounds()`, `.hsd_check_step()` and
   `.hsd_check_start()` validate and normalise each argument. Each abort
   names the limit it enforces (`hsd_opt_limits = c(-10, 2)`,
   `hsd_opt_steps`).
2. The function is created and its environment is set to the package
   namespace.
3. The attributes and the class are set.

### Why it is done this way

`gsDesign` spending functions already have three formals, with
`param = NULL`, so arity cannot distinguish an optimised parameter from a
fixed one; the class is the tag (`is_sf_hsd_opt()`). Without the environment
reset, the closure environment would hold the arguments and two calls with
equal arguments would not be `identical()`. `.hsd_groups()` compares group
members with `identical()` on their attributes. Binary-fraction steps make
`k * step` exact, so nodes never drift past a bound and coincide across
parameters (tables are shared by key).

### Edge cases

`start` outside `bounds` aborts; `step = 0.1` aborts with the reason ("binary
fractions keep every lattice node exact"). `bounds` with lower = upper aborts.


---


## `.transform_pvalues_gsd_spec()` -- the spec builder

### Purpose

Does everything that does not depend on the parameter values, once: it
validates, transforms the fixed hypotheses, short-circuits the hypotheses
that spending cannot affect, and builds every node table up front. Called
only from the handoff in `transform_pvalues_gsd()`.

### Inputs and output

- `pvals`: the raw array, already checked by the caller.
- `info_frac`, `look_back`: already normalised by `.gsd_info_frac()` and
  `.gsd_look_back()`.
- `spending`: the length-*m* list from `.gsd_spending_list()`.
- `alpha`, `grid_size`, `spending_label`: as in the caller.
- `recycled`: `TRUE` when the user passed a single function.

Returns a `multigrain_pvals_gsd_spec` (fields in the notation paragraph).
When no parameter remains after short-circuits, it returns a
`multigrain_pvals_gsd` instead.

### How it works

1. `.hsd_groups()` gives the group map, the parameter names, the bounds
   matrix, the steps and the user starts.
2. `.hsd_lattice()` builds each lattice. Each start node is the user's
   `start`, or the bounds' midpoint, snapped by `.hsd_snap()` (nearest node,
   ties to the lower).
3. `.hsd_fixed_labels()` labels the hypotheses: the list names if all are
   set, else the deparsed arguments of a literal `list(...)` call recovered
   by `.hsd_spending_expr()`, else `.gsd_spending_labels()`.
4. Each untagged entry is checked with `.gsd_check_spending()`.
5. Each hypothesis is then handled in turn:
   - a tagged hypothesis with more than one look, or one look below full
     information, is left for the search; `.gsd_warn_matured()` gives the
     fixed path's warning on ignored post-maturity values;
   - a tagged hypothesis with a single look at full information is
     short-circuited. It is transformed with `.hsd_node_fn(f, start node)`
     (after `.gsd_check_spending()`) and labelled
     `sfHSD(param = <start>)`. Its boundary is the level whatever the
     spending function, so the choice of node does not matter;
   - every fixed and short-circuited hypothesis is transformed with the
     unchanged `.gsd_transform_hyp()` and stored in `fixed` and
     `fixed_tables`.
6. Short-circuited hypotheses are removed from `group`; empty groups
   disappear and the rest are renumbered 1..*d*.
7. If *d* = 0, a `cli::cli_inform()` message explains why, and the stored
   transforms are assembled into a `new_pvals_gsd()` object.
8. Otherwise `.hsd_param_names()` names the remaining groups (an `id` keeps
   its name; unnamed parameters are renumbered). `.hsd_node_tables()` then
   builds the tables, and `new_pvals_gsd_spec()` wraps everything, with the
   spec labels from `.hsd_labels()`.

### Why it is done this way

Building every node's table here discharges "the well-ordering check runs at
every value tried": the search can only try nodes, and every node has passed
before the search exists. It also lets one spec serve several gains and
stopping rules on the same draws. No random numbers are drawn, so the RNG
stream of `graph_optimise_gsd()` is the same whether or not the spec was
built in the same session.

### Edge cases

An all-short-circuited spec returns an object `identical()` to the fixed
path's (P1 test). An untagged call never enters this function. Two groups
whose `id`s give the same name abort in `.hsd_param_names()`.


---


## `.hsd_groups()` -- the group map

### Purpose

Turns the tagged entries of a spending list into parameters.

### Inputs and output

`spending` (length *m*), `m`, `recycled`. Returns a list with:

- `group`: length *m*, `NA` untagged, else 1..*d* in order of first
  appearance;
- `names`, from `.hsd_param_names()`;
- `ids`: a list, `NULL` where none was given;
- `bounds`: a *d*-by-2 matrix with dimnames;
- `step` and `start` (`NA` for a `NULL` start).

### How it works

Each tagged entry gets a key: `"id:<id>"` if it has an `id`, `"recycled"` if
the whole argument was one function, else `"hyp:<i>"`. `match()` against the
unique keys gives the group numbers. For each group, every member must have
`identical()` `bounds`, `step` and `start` to the first member; otherwise it
aborts, naming the members and the field.

### Why it is done this way

Recycling a single function has always meant "the same function for every
hypothesis" in `.gsd_spending_list()`, so a single recycled tag is one shared
parameter. In a list of *m* entries, however, `rep(list(sf_hsd_opt()), m)`
without `id`s is *m* separate parameters. That is why the caller passes
`recycled`, recorded before the list is built.

### Edge cases

No tagged entry gives *d* = 0 (not reached from the handoff). A group whose
members differ only in `start` aborts.


---


## `.hsd_lattice()` -- the nodes of one parameter

### Purpose

Gives the sorted lattice of one parameter.

### Inputs and output

`bounds` (length 2) and `step`. Returns a sorted numeric vector.

### How it works

`step * seq(ceiling(lower / step), floor(upper / step))` gives the multiples
inside the bounds. `union()` with the two bounds is then sorted.

### Why it is done this way

Both bounds are always nodes, so a lattice never has fewer than two nodes
(lower < upper is enforced), and a user's off-grid bound such as −3.9 is
still searched. The multiples are exact because `step` is a binary fraction.

### Edge cases

`bounds = c(-3.9, -3.8)` with `step = 0.25` has no inner multiple; the
lattice is the two bounds.


---


## `.hsd_node_tables()` -- one boundary table per (schedule, node)

### Purpose

Builds and checks every table the search can use. Each table is checked for
well ordering (inside `.gsd_boundary_table()`) and for invertibility
(`.hsd_check_invertible()`).

### Inputs and output

`spending`, `info_frac`, `group` (after short-circuits), `lattice` (named),
`names`, `alpha`, `grid_size`. Returns `list(tables, schedule_key)`:

- `tables` is keyed by `.hsd_table_key()` = `"<schedule>|<node>"`, where the
  schedule is the comma-joined information fractions of the hypothesis's
  distinct looks;
- `schedule_key` has length *m*, `NA` for fixed hypotheses.

### How it works

1. If verbosity is not `"silent"`, a `cli::cli_progress_bar()` is created
   with the format "Building the boundary table of {building}".
2. For each active hypothesis, `.gsd_looks()` gives its distinct looks and
   their schedule key. Then, for each node of its group's lattice:
   - `.hsd_node_fn(f, node)` = `function(a, t) f(a, t, param = node)` is
     checked with `.gsd_check_spending()`;
   - if the key has already been built (another hypothesis on the same
     schedule and node), the build is skipped;
   - otherwise the progress label is updated, and `.gsd_boundary_table()` is
     called inside `rlang::try_fetch()`. An error is re-raised as
     "Hypothesis i is not well ordered at <name> = <node>", with the hint
     "Narrow `bounds` of `sf_hsd_opt()`" and the original error as parent;
   - `.hsd_check_invertible()` runs on the new table.

### Why it is done this way

- `alpha` and `grid_size` are constant within a spec, so the schedule alone
  keys a table. The tables are deduplicated across hypotheses, and a P1 test
  checks them `identical()` to per-hypothesis builds.
- `.gsd_spend()` calls the spending function with two arguments, so the tag
  is wrapped per node (build plan, reconciliation 1).
- The bar is a progress *bar*, not a `cli_progress_step()` per table
  (discovery decision 5 asked for a step). A bar appears only once the build
  has run `cli.progress_show_after` seconds (2 by default), so the test
  suite's small specs print nothing. A step at default verbosity leaked into
  the snapshots. Only the tagged path reaches this function, so fixed-path
  output is unchanged.

### Edge cases

The progress label variable is assigned only for the bar's format string (a
`nolint` comment marks it). A spending function that fails
`.gsd_check_spending()` at one node aborts with the fixed path's message.

### `.hsd_check_invertible()` -- build-time invertibility

`.gsd_invert()` aborts on a column with fewer than two distinct positive
boundaries, and `.gsd_boundary_table()` does not check this (build plan,
"Facts verified during planning"). This helper applies the same test to every
column at build time: `sum(!duplicated(col) & col > 0) >= 2`. A failure
aborts, naming the hypothesis, the node, the look and the "narrow `bounds`"
hint, so the abort can never fire inside the search. Phase D showed every
node of `seq(-10, 2, 0.25)` passing on four schedules and two corners, with
at least 1024 distinct boundaries (discovery D3), so the check is a guarantee
rather than an expected failure.


---


## `.gsd_invert_hyp()` -- invert a prebuilt table

### Purpose

Steps 4–7 of `.gsd_transform_hyp()` on a table that already exists. It is
used by `resolve_pvalues_gsd()` (full sample) and `.hsd_bank()` (GA rows).

### Inputs and output

- `raw`: an n-by-K matrix of raw p-values of one hypothesis;
- `t_row`: its information fractions;
- `tab`: a node table;
- `look_back`: logical.

Returns `list(values, table)`: the n-by-K repeated (or sequential) p-values,
and the table entry as the fixed path stores it (`looks`, `maturity`, plus
the table's fields).

### How it works

1. `.gsd_looks(t_row)` gives the distinct looks and the maturity look.
2. Each distinct look is inverted with `.gsd_invert(tab$bounds[, l],
   grid = tab$grid, p = raw[, look])`. Looks without data stay 1.
3. Looks after maturity copy the maturity column forward.
4. Under look-back, a running `pmin()` across looks is applied.

### Why it is done this way

The step is duplicated, not shared, so `.gsd_transform_hyp()` and with it the
fixed path are untouched (house rule: share or duplicate, duplicate). The P1
test checks `resolve_pvalues_gsd()` `identical()` to the fixed transform at
every node, which is what makes the duplication safe. One inversion costs
8–12 ms per hypothesis on 1e5 rows (discovery D1b), against ≈ 0.9 s for the
table.

### Edge cases

The post-maturity warning is not repeated here; the spec builder gave it once.


---


## `resolve_pvalues_gsd()` -- a spec at chosen nodes

### Purpose

Returns the `multigrain_pvals_gsd` that the fixed transform would give with
the tagged hypotheses fixed at the chosen values, from the stored tables. It
is exported. It is used by Stage G (re-scoring the winner), by Stage L (every
node array, via `.hsd_node_array()`), by the final evaluation, and by users
for fresh-draw checks.

### Inputs and output

- `spec`: a `multigrain_pvals_gsd_spec`; anything else gets
  `rlang::stop_input_type()`.
- `spending_param`: one finite number per parameter. If named, the names
  must be the parameter names, in any order.

Returns `new_pvals_gsd()` with labels `sfHSD(param = <value>)` for active
hypotheses.

### How it works

1. `.hsd_resolve_nodes()` checks the length, finiteness and names, reorders
   by name, and maps each value to the node within `1e-12`. Otherwise
   `.hsd_abort_off_lattice()` aborts, naming the nearest node(s), the
   lattice's ends and its step.
2. An array of ones is filled hypothesis by hypothesis:
   - fixed and short-circuited hypotheses take `spec$fixed[[i]]`;
   - active ones take `.gsd_invert_hyp()` on the table keyed by their
     schedule and node, and their `tables` entry is replaced.
3. Labels come from `.hsd_labels(..., values = nodes)`.

### Why it is done this way

Refusing off-lattice values keeps "every value tried is a checked node"
literally true. The result carries no reference to the spec, so the fixed
consumers (`calc_power_pvals_gsd()`, `prune_graph_gsd()`) take it unchanged.

### Edge cases

A named vector with a wrong or duplicated name aborts. A value within 1e-12
of a node resolves to the exact node value.


---


## `.hsd_labels()` -- spending labels

### Purpose

Gives the per-hypothesis spending labels for the spec and for a resolved
object.

### Inputs and output

`fixed_labels`, `group`, `names`, `bounds`, optional `values`. Returns a
character vector of length *m*.

### How it works

Fixed and short-circuited entries keep `fixed_labels`. An active entry reads
`sfHSD(param = <name> in [lo, hi])` on the spec (`values = NULL`) and
`sfHSD(param = <value>)` when values are given, formatted with
`format(digits = 15)`.

### Why it is done this way

The fixed path's labels are deparse-based and would print the tag's body
(build plan, reconciliation 7). `.hsd_fixed_labels()` names the fixed entries
as the user wrote them (list names, or the arguments of a literal
`list(...)`), so a mixed list reads `gsDesign::sfLDOF` next to
`sfHSD(param = -0.25)`.

### Edge cases

When the `spending` argument is not a literal `list(...)` of *m* entries (a
variable, or a list built elsewhere) and its elements are unnamed,
`.hsd_fixed_labels()` falls back to the fixed path's `.gsd_spending_labels()`.
`.hsd_spending_expr()` returns `NULL` if the call frame cannot be read.


---


## `.graph_optimise_gsd_opt()` -- the driver

### Purpose

The body of `graph_optimise_gsd()` for a spec: validation, control, the two
stages, the choice, the final evaluation and the result.

### Inputs and output

The arguments of `graph_optimise_gsd()` with `spec` in place of `pvals`,
plus `call` (the error context). Returns a `multigrain_graph_optimal` with the
fixed path's fields, plus:

- `spending_param`, `hyp_spending_param` and `spending`;
- `spending_bounds`, `spending_step` and `spending_group`;
- `spending_profile`, `spending_search` and `profile_converged`.

### How it works

1. The fixed path's checks run in its order:
   - `check_graph_constraint()`, `check_trial_success_gsd()`,
     `check_stopping_rule_gsd()`, `rlang::check_number_whole(num_threads)`,
     `check_control()`, `check_logical(global_search)`;
   - verbosity is normalised;
   - `.gsd_check_alpha(alpha, spec)` (the spec carries `alpha`),
     `.gsd_check_graph_dims()`, `.gsd_check_gain_dims()` and
     `.validate_start_graphs()`;
   - the stop table is built by `.gsd_stop_table(
     .gsd_check_stopping_rule_dims(rule, spec), m, K)`.

   They read only `nsim`, `m`, `K` and `alpha`, which the spec shares with
   the p-values it resolves to.
2. Settings: *d* = `length(spec$param$names)` and `p_graph = ncol(
   .build_start_matrix(graph_constraint, start_graph))`. The user's `popSize`
   is recorded before `control_prepare_dims(control, nsim = spec$nsim, m =
   spec$m)`; if it was unset, `popSize` becomes `min(max(40 (m + d), 200),
   500)`.
3. With `global_search`, Stage G runs. The local start is then the GA
   winner's graph genes clamped to [0, 1], at `ga_spending_param`. Without
   it, the start is `.hsd_start_genes(graph_constraint, start_graph)` at
   `spec$param$start`.
4. Stage L runs with the GA archive (or `NULL`).
5. `choose_graph(ga_result, local_result)` picks the source. `theta_star` is
   the chosen source's spending values, re-snapped and named by
   `.hsd_resolve_nodes()`.
6. `pv_star <- resolve_pvalues_gsd(spec, theta_star)` is passed to
   `prune_graph_gsd(..., gamma = 1, stop_rule)`, then
   `calc_power_pvals_gsd(..., custom_power = list(trial_success =
   trial_success), stopping_rule = stopping_rule)`.
7. Names come from `graph_constraint_get_names()`. `hyp_spending_param` is
   `theta_star[spec$group]`.
8. `graph_optimal()` is called with the fixed path's arguments and
   `spending = list(...)` holding the nine fields; `spending_search` comes
   from `.hsd_search_table()`.

### Why it is done this way

- `control_prepare_dims()` is not edited (build plan, reconciliation 5); the
  population override counts each parameter as one more hypothesis. The rule
  is inert at *m* + *d* ≤ 5 (it gives 200).
- `choose_graph()` compares full-sample gains, and each stage's winner is
  scored at its own nodes, so the comparison is on gains (reconciliation 11).
  When the GA wins, the profile is still Stage L's, which started at the GA's
  nodes; roxygen says so.
- The final evaluation is on the same tables that were checked (record 4.12,
  *Error control*), so the P3 test
  `identical(res$power, calc_power_pvals_gsd(resolve_pvalues_gsd(spec,
  res$spending_param), ...))` holds.

### Edge cases

`global_search = FALSE` gives `solution$graph_valid["global"] = NA` and an
all-`NA` `spending_search`. A user `control_global(popSize = 31)` survives
(P3 test).

### `.hsd_start_genes()` and `.hsd_search_table()`

`.hsd_start_genes()` returns `create_start_params(graph_constraint)` (the
default Holm-type start) when the first `start_graph` is empty. Otherwise it
fills a missing weight vector with ones and a missing matrix with
`1 / (m - 1)` off the diagonal, as `.build_start_matrix()` fills user rows,
and encodes the graph with `.hsd_encode_graph()`. Unlike the fixed path,
which always starts COBYLA from row 1 of `.build_start_matrix()` (the Holm
seed; record 10 item 16(a)), the user's graph is honoured.

`.hsd_search_table()` lists, per parameter and node, the archive's `visits`
and `best` (`-Inf` becomes `NA`). Without an archive these are all `NA`.


---


## `.hsd_bank()` -- the slice bank of Stage G

### Purpose

Precomputes the repeated p-values of every active hypothesis at every node of
its lattice, on the GA rows, so that a GA evaluation never inverts a table.

### Inputs and output

`spec` and `idx` (the GA rows, in order). Returns `list(slices, hyp, idx,
bytes)`: `slices[[a]][[q]]` is the `length(idx)`-by-K matrix of hypothesis
`hyp[[a]]` at node *q*.

### How it works

1. `.hsd_bank_bytes()` = `8 K length(idx) Σ_active M_g(i)` is checked by
   `.check_bank_memory()`: a warning above 2 GB, an abort above
   `max_bank_gb = 4`, each naming three remedies (a coarser `step`, narrower
   `bounds`, a smaller `nsim_global`).
2. For each active hypothesis the GA rows of its raw p-values are taken once.
   Each node's table is then inverted with `.gsd_invert_hyp()`.

### Why it is done this way

Only the GA needs random access across nodes within one evaluation; Stage L
resolves one node at a time and holds no bank (design answer to question 2).
At the defaults (*M* = 21, *K* = 2, N_GA = 5e4, two active hypotheses) the
bank is 33.6 MB (discovery D5).

### Edge cases

The memory guard is tested on synthetic byte counts only (not at real sizes).


---


## `create_obj_func_gsd_opt()` -- the objective of the joint search

### Purpose

Twin of `create_obj_func_gsd()` for the gene vector with parameter genes
appended. Built once per Stage G run.

### Inputs and output

- `spec` and `bank`;
- `power_criterion`, `hyp_constraint`, `trans_constraint`, `alpha`,
  `num_threads` and `stop_rule`, as for `create_obj_func_gsd()`;
- `p_graph`;
- `start_param`, the lattice values at which the working matrix is first
  filled (default `spec$param$start`).

Returns `function(x)`, which gives the gain or the penalty of an infeasible
graph. It carries the attributes `archive` (the state environment), `d` and
`p_graph`.

### How it works

1. In the factory, all arguments are `force()`d. Then it computes:
   - for each active hypothesis, its kernel columns
     `(seq_len(K) - 1) m + i`;
   - the members of each parameter;
   - the start nodes, via `.hsd_resolve_nodes()`.
2. One environment, `state <- new.env(parent = emptyenv())`, holds:
   - the archive (`.hsd_archive_init()`): per parameter, `best` (−Inf),
     integer `visits` and an *M_r*-by-`p_graph` `genes` matrix;
   - `last_q`, the node of each parameter the matrix currently holds;
   - `mat`, the working kernel matrix. It is built once by
     `.hsd_working_matrix()`, the `.gsd_kernel_matrix()` of
     `resolve_pvalues_gsd(spec, start_param)` on rows `idx`, and never bound
     to a second name.
3. Each call `fn(x)`:
   - takes the graph genes `xg <- x[seq_len(p_graph)]` and maps the
     parameter genes with `.hsd_node_index()`;
   - decodes and penalises the graph genes with the block copied verbatim
     from `create_obj_func_gsd()` (`split_theta()`, `recover_full_*()`,
     the `-1e6` / negative / above-one penalties and the small-value
     snapping);
   - overwrites the columns of every hypothesis whose parameter's node
     differs from `last_q` with that node's bank slices, through
     `.hsd_set_columns()`, and updates `last_q`;
   - calls `graph_shortcut_gsd()` or its parallel twin with `state$mat`
     (exactly as the fixed closure does) and applies `power_criterion` to
     `res$time`;
   - updates the archive with `.hsd_archive_update()`: visits + 1 at each
     parameter's node, and `best` and `genes` replaced on a strict
     improvement. Only evaluations that reach the kernel count; penalty
     returns do not.

### Why it is done this way

- The decode block is duplicated, not shared (house rule), and the P2 node
  oracle checks the closure `identical()` to `create_obj_func_gsd()` on
  `resolve_pvalues_gsd(spec, node)[idx]` for 20 gene vectors at every node.
- Mutable state lives in one explicit environment created in the factory
  body, with no `<<-` (house rule of 2026-10-14; a test greps the new files).
- `.hsd_set_columns()` takes the matrix out of the environment
  (`mat <- state$mat; state$mat <- NULL`), assigns the columns and puts it
  back. Then the local name is the only reference and `[<-` modifies in
  place. A direct `state$mat[, j] <- v` inside a function copies the whole
  matrix on every call. P2 measured this; d03 had probed only top-level
  assignment, which does not copy. A `tracemem()` test confirms no copy per
  evaluation.
- The column copy costs about 8 % of a kernel call at N_GA = 5e4
  (discovery D5), so no C++ column copier was added.

### Edge cases

Genes outside [0, 1] are clamped to the end nodes. An evaluation that leaves
every node unchanged copies nothing.

### `.hsd_node_index()` -- genes to nodes

`.hsd_node_index(u, M)` clamps `u` to [0, 1] and returns
`min(M, 1 + floor(u M))`, vectorised over `u` and `M`. Each node owns an
interval of width 1/M of the gene, so every node has the same prior mass.
`u = 1` maps to node *M* (P2 tests: −0.3 → 1, 1.7 → *M*, 1 → *M*). The GA
suggestions place the start node's gene at the centre of its interval,
`(q_start − 0.5) / M`.


---


## `.graph_optimise_ga_gsd_opt()` -- Stage G

### Purpose

Twin of `.graph_optimise_ga_gsd()` with one gene per parameter.

### Inputs and output

`spec`, `graph_constraint`, `trial_success`, `alpha`, `nsim`
(`control$nsim_global`), `global_opts`, `num_threads`, `start_graph`,
`verbose`, `stop_rule` and `p_graph`. Returns the fixed driver's list
(`ga_hyp_weight`, `ga_trans_matrix`, `ga_trial_success`, `ga_subset_power`,
`is_graph_valid`, `ga_output`) plus:

- `ga_spending_param`, the named lattice values of the winner;
- `archive`, a plain list of `best`, `visits` and `genes`;
- `idx_ga`.

### How it works

1. A progress step "Running global optimization" is shown unless silent.
2. `x0 <- .build_start_matrix(graph_constraint, start_graph)`, which draws
   nothing.
3. `idx_ga <- sample.int(spec$nsim, size = nsim)`, the draw
   `.sample_pvals_gsd()` makes.
4. The closure is built: `create_obj_func_gsd_opt(spec, bank = .hsd_bank(spec,
   idx_ga), ...)`.
5. The suggestions are `x0` with `(q_start − 0.5) / M_r` appended per
   parameter. The fixed driver's GA arguments are used (`"real-valued"`,
   bounds `rep(0/1, p_graph + d)`, `.cauchy_population`,
   `.make_cauchy_mutation_multi(0.1, 1)`, `optim = TRUE`), merged with
   `global_opts` by `utils::modifyList()`. Then `GA::ga()` runs.
6. The archive is copied out as a plain list. The working matrix is dropped,
   along with the GA object's `@call` slot (which holds the closure and with
   it the bank; `graph_optimal()` drops that slot anyway).
7. The winner `@solution[1, ]` is sliced. Its parameter genes are mapped by
   `.hsd_node_index()` to `theta_ga`, and its graph genes are decoded by
   `param_to_solution(process = TRUE)` without `repair_graph()`, exactly as
   the fixed GA driver does.
8. If the graph is valid, it is scored on `.hsd_subset_pvals(
   resolve_pvalues_gsd(spec, theta_ga), idx_ga)` with
   `.calc_power_pvals_gsd_core()`, and on the full sample by
   `graph_shortcut_gsd()` + `trial_success$func(res$time)`. The latter is
   `ga_trial_success`.

### Why it is done this way

The RNG order of the fixed driver is kept (start matrix, subsample, GA), so
the same seed gives the same result twice (P3 test (c)). Clamping rather than
penalising out-of-range parameter genes keeps the individual's graph genes
meaningful. Phase D measured the gene's cost at 0.90–0.98× the generations of
a GA without it and 1.1–1.2× the wall time. Every gene run visited all 21
nodes and found the profile's argmax node (M0) or a node inside its plateau
(M1; discovery D8).

### Edge cases

An invalid winner gives `NULL` gains and lets `choose_graph()` take the local
result. The archive's values are on `idx_ga` and are never compared with
Stage L's gains.


---


## `.graph_optimise_local_gsd_opt()` -- Stage L

### Purpose

Refines the graph at the GA's parameter values and then traces the
warm-started profile along each parameter, which is also the local
parameter search.

### Inputs and output

As Stage G, plus `local_opts`, `x0` (graph genes, length `p_graph`), `theta0`
(lattice values), `archive` (or `NULL`) and `density` (1 profiles every node;
2 profiles every second node plus the bounds; not exposed, default 1 by
discovery decision 2). Returns the fixed local driver's list
(`local_hyp_weight`, `local_trans_matrix`, `local_trial_success`,
`local_subset_power`, `is_graph_valid`, `local_output` = step 1's `nloptr`
result) plus `local_spending_param`, `spending_profile`,
`profile_converged` and `idx_local`.

### How it works

1. A progress step is shown unless silent. `idx_local <- sample.int(N,
   nsim_local)` is drawn once.
2. `.hsd_local_context()` packs into one list:
   - the spec, the graph constraint, the gain and `alpha`;
   - `num_threads`, `stop_rule` and `local_opts`;
   - `idx = idx_local`, the archive and `density`;
   - `n_block = 100`.
3. Step 1. The node array `.hsd_node_array(ctx, theta0)` is
   `resolve_pvalues_gsd()` restricted to `idx_local` in that order. Its
   kernel matrix is built, and `.hsd_refine_node(x0, obj = .hsd_objective(ctx,
   kmat), kmat, ctx)` runs. The result becomes `current`: theta, graph, gain,
   block gains and `version = 0`.
4. Step 2. A progress step "Profiling spending parameters: …" is shown, then
   `.hsd_profile(current, ctx)` runs.
5. If the final graph is valid, it is scored on the pinned rows
   (`.calc_power_pvals_gsd_core()`) and on the full sample (`local_trial_success`).

### Why it is done this way

The fixed local driver draws its own subsample inside, so calling it would
break the pairing of Stage L. Step 1 is therefore `.hsd_refine_node()`, the
same `nloptr` call on the pinned rows (build plan, reconciliation 6). The P3
oracle (a) compares it with a by-hand `nloptr::nloptr()` call. At the default
`nsim_local = N` the draw is still a permutation, so the order must be pinned
for the block means to be paired.

### Edge cases

With `global_search = FALSE`, `x0` is the encoded user start graph (or the
default) and `theta0` the start node; step 2 is then the whole parameter
search.

### `.hsd_profile()` -- rounds and the final recompute

1. Up to 3 rounds for *d* ≥ 2, one for *d* = 1. In each round, for each axis
   *r*:
   - `version_before` is recorded and `rows[[r]] <- .hsd_profile_axis(r,
     current, ctx)` is computed;
   - if the axis's best row beats `current$gain` (strictly), `current` takes
     its value, graph, gain and blocks, and `version` increments. `moved_any`
     is set only if the node changed;
   - the axis is stamped with `version_before`.

   The rounds stop when no axis moved.
2. `converged <- d == 1 || !moved_any`.
3. Every axis whose stamp differs from the final `version` is recomputed once
   with `final = TRUE`. That includes a *d* = 1 axis whose own sweep improved
   the optimum, and earlier axes after a later one moved. So every reported
   axis passes through the final optimum.
4. `.hsd_profile_table(rows, names, current$blocks)` builds the data frame.
   Any `diff < 0` sets `converged <- FALSE`; the sign is kept, not clipped.

The recompute can expose a node better than the optimum, because the
recomputed axis starts from a graph different from the one it was swept
with. The build plan asked that the sign be kept honest, and
`profile_converged` is how the display learns of it.


---


## `.hsd_profile_axis()` and `.hsd_profile_node()` -- one axis, one node

### Purpose

`.hsd_profile_axis(r, current, ctx, final)` gives the profile rows of axis
*r* in lattice order. `.hsd_profile_node(r, q, current, prev, ctx, refine)`
gives one row.

### Inputs and output

A row is a list of `value`, `gain`, `gain_fixed_graph`, `n_eval`,
`hyp_weight`, `trans_matrix` and `blocks` (block gains on the pinned rows).

### How it works

1. The axis starts at the current node `at` (with `density = 2`, the nodes
   at even distance from `at` plus both bounds).
2. The row at `at` is computed first, with `refine = !final`. The walks then
   go upward and downward from `at`, each row warm-started from the previous
   row in its direction (`prev`).
3. At node *q*, `theta_q` is the current theta with axis *r* set to the
   node's value. The node's kernel matrix is built on the pinned rows. Up to
   three candidates are evaluated with `.hsd_eval_graph()` (the kernel plus
   the gain; an invalid graph scores `-Inf`):
   - `star`, the current optimum's graph, listed first so that it wins ties;
   - `prev`, the previous row's graph;
   - `archive`, the GA's best graph at that node, from
     `.hsd_archive_graph()`, decoded as the GA winner is (`NULL` if
     unvisited, absent or invalid).
4. With `refine = FALSE`, the row is the star graph's gain in both columns,
   with `n_eval = 0`.
5. Otherwise `.hsd_refine_node()` starts from the best candidate, encoded by
   `.hsd_encode_graph()`. The row is the best of the candidates and the
   refined graph, and `gain_fixed_graph` is the star graph's gain. The row's
   blocks are those of the winner.

### Why it is done this way

- Taking the best of the candidates and the refined graph gives
  `gain >= gain_fixed_graph` at every row (P3 test (d)), that is,
  S ≤ P̃ (design 2.5).
- Warm-starting from the neighbour follows the optimum's graph along the
  axis. The archive start covers a far node whose best graph lies in another
  basin. Phase D measured the warm-started profile within 2.1e-5 (M0),
  2.3e-5 (M1) and 3.7e-4 (§2 design) of the best of three full optimisations
  inside the plateau (discovery D4), without the archive start, which
  P3 added.
- In the final pass the optimum's own node is not refined, so there the star
  graph is the optimum by definition: `gain == gain_fixed_graph`, and
  `diff == se == 0` exactly. This makes the fixed-graph slice and the profile
  agree at θ*, as design 2.5 requires. A refinement there could find a
  slightly different graph and break the identity.

### Edge cases

`n_eval` per node was 23–35 at one graph gene and 83–128 at five (discovery
decision 6), at most 2.6 % of `maxeval`, so the budget is the main pass's
`local_opt`.


---


## `.hsd_refine_node()` -- COBYLA on the graph genes at one node

### Purpose

The `nloptr` call of the fixed local driver, on a given closure and start.

### Inputs and output

`x0` (graph genes), `obj` (from `.hsd_objective(ctx, kmat)`, i.e. the
unchanged `create_obj_func_gsd()` on `kmat`), `kmat` and `ctx`. Returns
`x`, `hyp_weight`, `trans_matrix`, `gain`, `time`, `blocks`, `n_eval` and
`nloptr_output`.

### How it works

1. A counter in its own environment wraps the negated objective.
2. `nloptr::nloptr(x0, eval_f, lb = 0, ub = 1, opts = ctx$local_opts)` runs.
3. The solution is decoded with `param_to_solution(process = TRUE)` and then
   `repair_graph()`, as the fixed local driver does.
4. The repaired graph is evaluated on `kmat` with `.hsd_eval_graph()`, and
   `.hsd_block_gains()` gives its block gains.

### Why it is done this way

Reusing the fixed closure and options keeps the oracle: the P3 test (a)
checks the result `identical()` to a by-hand `nloptr()` call with the same
closure, `x0`, options and seed. COBYLA never sees the parameter genes, which
are piecewise constant (design 2.4).

### Edge cases

The count includes every objective call, penalties included.


---


## `.hsd_encode_graph()` -- warm-start encoding

### Purpose

Encodes a decoded graph `(w, G)` as COBYLA start genes.

### Inputs and output

`gc`, `w` and `G`. Returns a numeric vector of length `p_graph`.

### How it works

`closest_graph_to_constraints(gc, w, G)` projects the graph onto the
constraint, and `create_start_params(gc, w0, G0, sum_to_one_constraint =
FALSE)` encodes it. This is the user-row branch of `.build_start_matrix()`.

### Why it is done this way

`.build_start_matrix(...)[1, ]` is always the Holm seed, whatever the user
passed (`R/optimisation_start.R`, `seeds[[1]] <-
create_start_params(graph_constraint)`). That is the pre-existing bug of
record 10 item 16(a), so it cannot serve as an encoder. The P3 test (a)
checks that `param_to_solution(.hsd_encode_graph(gc, w, G))` reproduces
`(w, G)` for three graphs.

### Edge cases

A graph that violates the constraint is projected first, so the encoded start
decodes to the projection, not to the input.


---


## `.hsd_block_se()` -- paired standard errors

### Purpose

Gives the batch-means paired difference and standard error of two
evaluations on the same rows (design 1.4).

### Inputs and output

`.hsd_block_se(time_a, time_b, func, B = 100)` returns `list(diff, se)`.

### How it works

`.hsd_block_gains()` splits the rows into `min(B, n)` contiguous blocks
(`floor(seq(0, n, length.out = B + 1))` cut points) and applies the compiled
gain to each block's time matrix. `.hsd_block_diff(a, b)` returns
`mean(a - b)` and `sd(a - b) / sqrt(B)`.

### Why it is done this way

No change to `trial_success_gsd()` is needed, and 100 calls on N/100 rows
cost far less than one COBYLA run. The driver keeps each row's block gains
(the `blocks` field) and calls `.hsd_block_diff()` against the optimum's
blocks in `.hsd_profile_table()`. `.hsd_block_se()` is the two-matrix
convenience form, used by the tests.

### Edge cases

With equal block sizes, `diff` equals the full-sample difference up to
summation order (unequal blocks weight rows slightly differently). At the
optimum's own row the blocks are identical, so `diff` and `se` are exactly 0.


---


## `.hsd_plateau()` -- the plateau interval

### Purpose

Gives the connected run of nodes around a value whose shortfall is within a
tolerance (design 1.5; the user testing's plateau convention).

### Inputs and output

`.hsd_plateau(value, diff, tol = 5e-4, at)` returns `c(lower, upper)`.

### How it works

Sorts by `value`, then extends from `at` downward and upward while the
neighbour's `diff <= tol`.

### Why it is done this way

Only a *connected* run is reported: a far node that happens to come within
tolerance on the other side of a dip is not part of the plateau (P4 test with
a gap).

### Edge cases

`at` must be one of `value`. The display falls back to the row with the
smallest `diff` if the chosen value is missing.


---


## `.hsd_cat_print()` and `.hsd_cat_summary()` -- display

### Purpose

The spending blocks of `print()` and `summary()` of a
`multigrain_graph_optimal` from the optimised path. They are called from
`R/graph_optimal.R` only when `!is.null(x$spending_param)`.

### Inputs and output

The result object. Both print and return their input invisibly.

### How it works

`.hsd_cat_print(x)`:

1. Prints "Spending functions:" and one line per hypothesis (names from
   `hyp_weight`, else `H1..Hm`). The notes from `.hsd_spending_notes()` are
   `[optimised: <name>, bounds a to b, step s]` on a group's first member,
   `[optimised: <name>]` on later members, `[optimised: a or b]` when the
   group cannot be told, and `[fixed at v: one look at full information]`
   for a short-circuited hypothesis.

   The group map is `spending_group` when present; otherwise
   `.hsd_infer_groups()` enumerates the maps consistent with the values.
   Parameters whose bounds could not go on a first member get an extra line.
   `.hsd_fmt_lattice()` reads the bounds and infers the step from the node
   spacing in `spending_search` (else the profile).
2. Per parameter: "Gain profile in <name> (graph re-optimised at each
   value):" and "within 5e-4 of the maximum on [l, u]; fixed-graph plateau
   [l', u']", where:
   - the first interval is `.hsd_plateau()` on `diff` around θ*;
   - the second is `.hsd_plateau()` on `max(gain_fixed_graph) -
     gain_fixed_graph`, i.e. measured from the fixed-graph column's own
     maximum.
3. If `profile_converged` is `FALSE`: "profile did not converge[ in 3
   rounds]; the plateau is indicative".

`.hsd_cat_summary(object)` prints, per parameter, "Spending-parameter
profile (<name>):" and a table of `value`, `gain`, `gain_fixed_graph`,
`diff` (rounded to 4), `se` (2 significant digits) and `n_eval`. It then
prints "GA visited k of M nodes", or "no global search".

### Why it is done this way

The blocks go after the trial success function and before the stopping rule
(build plan, reconciliation 8). The helpers live in the `_opt` file so that
`R/graph_optimal.R` gains only two guarded calls, and the fixed
`graph_optimal.md` snapshot stays byte-identical. The new display has its
own snapshot file, `_snaps/graph_optimal_hsd_opt.md`.

### Edge cases

A hand-built result without `spending_group` still prints (the inference
path; P4 snapshots cover *d* = 2 with a short-circuited hypothesis and
`profile_converged = FALSE`).
