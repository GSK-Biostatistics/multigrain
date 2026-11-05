<!--
Provenance: written by a Claude Code session on 2026-10-14 as the recommended design for
the job in dev/implement_hsd_opt.md, against gsd-build at 8af8adb. It revises an earlier
candidate design (hsd_opt_candidate_design.md) that was supplied to the session as an
upload and is not in this repository. The session ran in plan mode, so the Phase D
measurements it describes were planned but not run; cost figures are derived from the
brief's measured timings, as the header note below says. Saved here so the design has an
unambiguous origin in the repository. The body below is the session's document verbatim.
-->

# Joint Optimisation of Hwang–Shih–DeCani Spending Parameters in `multigrain` — revised design

> Revision of the candidate design (`hsd_opt_candidate_design.md`), written against the code on
> `gsd-build` (`8af8adb`), `dev/implement_hsd_opt.md`, `dev/gsd_design_record.md` (4.1, 4.5,
> 4.7, 4.10, 4.11, section 10) and `dev/review/gsd_user_testing/README.md` (6, 7.5, 8). Only one
> earlier review was attached (the candidate itself); if the two reviews referred to in the brief
> are elsewhere, point me at them and I will reconcile. Tags: **[code]** read from the tree,
> **[measured]** a number from the brief (2026-10-13 timings), the user testing or the record,
> **[derived]** follows from stated facts or measured numbers, **[estimate]** a projection.
>
> No new measurements were taken: the session is in plan mode, which forbids running scripts,
> so the two Opus measurement batteries I launched (transform cost and HSD well-ordering over
> $[-10, 2]$; optimiser evaluation counts and warm-started COBYLA across neighbouring gammas)
> stopped at their plans. They are ready to run and are the first thing Phase D does (4.3).
> Every cost figure below is therefore derived from the brief's measured timings and the user
> testing's run times, and is tagged accordingly.

## Changes from the candidate, in one list

1. **The local stage is an outer loop with warm starts (a true profile), not a fixed-graph
   sweep.** Each lattice node on each parameter axis through the current optimum gets its own
   warm-started COBYLA graph refinement on the full sample; the argmax of that profile is the
   local answer, and the same runs *are* the flatness report. The candidate's fixed-graph sweep
   is kept only as a second, cheaper column of the same table. Argument: section 1.3, 2.5 and
   the measurements in 4; the fixed-graph slice understates flatness exactly when spending and
   graph interact, which is the case the whole feature exists for (README 7.5: w1 moved
   0.32 → 0.46 when H2's spending changed).
2. **The global stage keeps the GA with parameter genes on a slice bank**, but the bank exists
   only there, sized by `nsim_global`; nothing else holds a bank (answers question 2).
3. **Linear blending is removed entirely**, including the extended-table comparator and
   discovery script `d02`. A blended slice is not the transform at any parameter value, so it
   fails the binding requirement that every value tried is a valid, well-ordered transform.
   Snapping is the only scheme.
4. **Out-of-range parameter genes are clamped, not penalised.** A clamped gene decodes to a
   valid node and a meaningful gain; a penalty throws away the graph genes of that individual.
5. **No special mutation operator.** `.make_cauchy_mutation_multi()` already acts per gene and
   respects `[0, 1]`; the parameter gene needs nothing more. (A mask will be needed only if the
   sparsity branch's zeroing mutation lands; noted, not built.)
6. **A single recycled `sf_hsd_opt()` means one shared parameter** (not an abort). Recycling a
   single function has always meant "the same function for every hypothesis" in
   `.gsd_spending_list()` **[code]**, and one gamma for all is the brief's original case.
7. **Node tables are built eagerly in `transform_pvalues_gsd()`**, so the well-ordering check
   for every value the search can ever try runs before any search, the spec object can be
   reused across gains and stopping rules on the same draws (README 7.5 pattern), and one grid
   size serves the search and the final evaluation. The two-grid option (`G_search`) is dropped.
8. **The result does not carry the transformed array.** It carries the parameters, labels,
   profile and archive; `resolve_pvalues_gsd(spec, spending_param)` rebuilds the exact
   `multigrain_pvals_gsd` in well under a second from the cached tables (section 4 numbers).
9. **GA budget rule stated** (question 3): `popSize = min(max(40 * (m + d), 200), 500)`,
   `run` and `maxiter` unchanged; measured in Phase D, not assumed.
10. **Dispatch in `graph_optimise_gsd()` is the first statement**, before `check_pvals_gsd()`
    (the candidate placed it after, where the spec would already have been refused).
11. **Default start node is the lattice node nearest the midpoint of the bounds**, not −2: no
    constant that is only right for the default bounds.
12. **Phase D is reshaped**: `d02` becomes "warm-started profile versus full GA at far nodes"
    (the measurement the new local stage depends on); the rest is renumbered.
13. **Warm starts are encoded with `create_start_params(gc, w0, G0, sum_to_one_constraint =
    FALSE)` after `closest_graph_to_constraints()`, never with `.build_start_matrix(...)[1, ]`.**
    Found while reading: row 1 of `.build_start_matrix()` is always the default Holm seed
    (`R/optimisation_start.R:9`), so `graph_optimise_gsd()` at `R/optimisation_gsd.R:230-233`
    (and `graph_optimise()` at `R/optimisation.R:188`) start COBYLA from Holm and ignore a
    user `start_graph` whenever `global_search = FALSE`. Pre-existing in both paths; per the
    house rule it is reported, not fixed here, and not replicated in the `_opt` code.

Decisions kept, each with the sentence that keeps it, are in section 3.

## The five questions, answered

1. **Outer loop with warm starts for one or two parameters?** Yes for the local stage, for
   every $d$: a warm-started re-optimised profile along each axis is the better product
   because the fixed-graph slice understates flatness precisely when spending and graph
   interact (2.5), and the same COBYLA runs are the local search, so the profile costs nothing
   extra beyond the search itself. Not for the global stage: an outer product grid is $M^d$
   full optimisations and a coordinate-only scheme can stall when two parameters must move
   together, so the GA keeps its parameter genes (3.1). For $d = 1$ the two stages together
   are "joint GA, then a full 21-node profile"; whether the GA's gene earns its keep at
   $d = 1$ is a `d04` measurement, but the code is needed for $d \ge 2$ regardless.
2. **Slice bank in the local stage?** No. Only the GA needs random access across nodes within
   one evaluation. Stage L resolves one node at a time from the cached tables (sub-second,
   4.1) and releases it. Appendix A.2 is GA-only.
3. **GA budget with parameter genes.** `popSize = min(max(40 (m + d), 200), 500)`, everything
   else unchanged; `run = 200` is self-adjusting so generations grow only while the gene keeps
   finding improvements. Expected 1.2–1.5 × today's GA time **[estimate]**; `d04` measures
   generations and evaluations with and without the gene and may revise the 40.
4. **Default flatness report.** One table per parameter over every lattice node on the axis
   through $\boldsymbol\vartheta^*$: re-optimised profile gain, fixed-graph gain, paired
   difference to the optimum and its batch-means SE ($B = 100$); `print()` condenses it to the
   two plateau intervals at $5\times10^{-4}$ (profile and fixed-graph); `summary()` prints the
   table; the GA archive is reported as node coverage, never as a gain comparable with the
   profile (it is on the GA subsample).
5. **Section 2 corrections.** Blending must go (a blended slice is not a transform at any
   $\vartheta$); the well-ordering guarantee is discharged at spec construction and sealed by
   `resolve_pvalues_gsd()` refusing off-lattice values; the winner's-curse bias of the reported
   gain is stated (2.2); CRN for the profile requires pinning the row order because
   `nsim_local = N` still permutes (2.3); the ordering
   $S \le \tilde P \le P$ holds only if $(\mathbf w^*, \mathbf G^*)$ is among the COBYLA starts
   at every node, which the design guarantees (2.5); and the $[-10, 2]$ well-ordering claim is
   argued but must be demonstrated, not assumed (2.1).

---

## Scope and Binding Requirements

Unchanged from the candidate and the brief: optional with a byte-identical fixed path
(`dev/gsd_identity_check.R`, the two GSD test files); up to $m$ parameters in arbitrary groups;
default bounds `c(-4, 1)`, user bounds within $[-10, 2]$ with lower < upper or an abort naming
the limits; chosen parameters and a gain-versus-parameter report on the result and in
`print()`/`summary()`; stopping rules and reporter fields untouched, kernel untouched; a valid
well-ordered transform at every value tried; raw draws fixed and one subsample index per search
stage; house rules of `dev/implement_gsd.md`; discovery before implementation.

## Notation

As in the candidate, with these additions and corrections.

| Symbol | Definition |
| --- | --- |
| $d$ | Number of optimised parameters, $0 \le d \le m$; $\mathcal I_{\text{HSD}}$ the active hypotheses, $g(i)$ their group. |
| $\Gamma_r$ | Lattice for $\vartheta_r$: the multiples of $\Delta_r$ inside $[\text{lower}_r, \text{upper}_r]$ plus the two bounds; $M_r = \lvert\Gamma_r\rvert$. Default $\Delta = 0.25$, $M = 21$ on $[-4, 1]$. |
| $u_r \in [0,1]$ | GA gene for $\vartheta_r$; node index $q_r = \min(M_r, 1 + \lfloor u_r M_r \rfloor)$ after clamping $u_r$ to $[0,1]$. |
| $S_r(\vartheta_r)$ | Fixed-graph slice: $\hat U_N(\mathbf w^*, \mathbf G^*, \vartheta_r, \boldsymbol\vartheta^*_{-r})$. |
| $\tilde P_r(\vartheta_r)$ | Warm-started profile: the gain of the COBYLA refinement at $\vartheta_r$ started from the better of the neighbouring node's optimum, $(\mathbf w^*, \mathbf G^*)$, and the GA's best graph at that node. $S_r \le \tilde P_r \le P_r$ **[derived]** (section 2.5). |
| $P_r(\vartheta_r)$ | True profile $\max_{\mathbf w, \mathbf G, \boldsymbol\vartheta_{-r}} \hat U_N$. Not computed. |
| $B$ | Number of contiguous trial blocks for batch-means standard errors (100). |

---

## 1. How It Works

### 1.1 Inputs: tagging hypotheses

`sf_hsd_opt(id = NULL, bounds = c(-4, 1), step = 0.25, start = NULL)` returns a function
`function(alpha, t, param) gsDesign::sfHSD(alpha, t, param)$spend` with class
`"multigrain_sf_hsd_opt"` and attributes `id`, `bounds`, `step`, `start`. It is therefore the
"(alpha, t, gamma) function" the brief asks for, and it is detectable without arity inspection
(every `gsDesign` spending function already has three formals **[code]**).

Grouping: in a list of $m$ entries, tagged entries with the same non-`NULL` `id` (normalised
with `as.character()`) share a parameter; `id = NULL` entries each get their own; a single
tagged function passed alone is recycled and is one shared parameter. Members of a group must
have `identical()` `bounds`, `step` and `start`, else abort. `bounds` must be finite, length 2,
lower $\ge -10$, upper $\le 2$, lower < upper, else abort naming the limits. `step` must be one
of `0.5, 0.25, 0.125` (binary fractions: `k * step` is exact in floating point **[derived]**, so
nodes never drift past a bound and coincide across parameters so tables are shared). `start`
is snapped to the nearest node; `NULL` means the node nearest the midpoint of the bounds.

### 1.2 The spec object: everything that does not depend on $\boldsymbol\vartheta$, done once

When any entry of `spending` is tagged, `transform_pvalues_gsd()` returns a
`multigrain_pvals_gsd_spec` (which does **not** inherit `multigrain_pvals_gsd`). It:

1. validates everything the fixed path validates (`.gsd_info_frac()`, `.gsd_look_back()`,
   `.gsd_spending_list()`, `.gsd_check_spending()` **[code]**);
2. transforms every fixed hypothesis exactly as today (`.gsd_transform_hyp()`) and stores the
   result;
3. applies the short-circuit rule to tagged hypotheses: a hypothesis with $\lvert L_i\rvert = 1$
   and $t \ge 1$ has boundary $= $ level for any spending function **[code]**, is transformed
   as today and dropped from $\mathcal I_{\text{HSD}}$ with `hyp_spending_param[i] = NA`; a
   group left empty disappears; if $d = 0$ after this, a plain `multigrain_pvals_gsd` is
   returned with an informational message;
4. builds the lattice $\Gamma_r$ per parameter and, for every distinct (information schedule,
   node) pair among active hypotheses, one boundary table with `.gsd_boundary_table()` at the
   user's `grid_size` **[code]**, which runs the well-ordering check. A failure aborts naming
   the hypothesis and the node value, with the hint to narrow `bounds`. This is where "the
   well-ordering check runs at every value tried" is discharged: the search can only ever try
   nodes, and every node's table has passed before the search exists;
5. stores the raw array, `nsim`, `m`, `K`, `alpha`, `info_frac`, `look_back`, `grid_size`, the
   fixed transforms, the group map, lattices, tables and labels.

`resolve_pvalues_gsd(spec, spending_param)` (exported) snaps each value to its node (abort if
not on the lattice, to keep "every value tried is a node" literally true), inverts the stored
tables for the active hypotheses (steps 4–7 of record 4.1: `.gsd_invert()` per look, matured
copy-forward, `NA` looks = 1, `cummin` under look-back), merges with the fixed transforms and
returns a `multigrain_pvals_gsd` whose `spending` labels read e.g. `sfHSD(param = -2.25)`.
`print()` of the spec shows dimensions, per-hypothesis spending (fixed label, or
"HSD gamma[r], bounds, step, M nodes"), tables built and any short-circuits.

### 1.3 The search

**Stage G — global, joint, on a bank.** `.graph_optimise_ga_gsd_opt()` draws `idx_ga`
(`sample.int(nsim, nsim_global)`, once), builds the slice bank: for each active hypothesis $i$
and node $q \in \Gamma_{g(i)}$ an $N_{\text{GA}} \times K$ matrix of repeated p-values on those
rows (inversion only; tables cached). The objective closure takes
$\boldsymbol\theta = (\mathbf w, \mathbf G\text{ genes}, u_1, \dots, u_d)$, decodes the graph
genes exactly as `create_obj_func_gsd()` does **[code]**, clamps $u_r$ to $[0,1]$, maps to
$q_r$, copies the $K$ columns of slice $q_{g(i)}$ for each active $i$ into the working
$N_{\text{GA}} \times mK$ kernel matrix (fixed hypotheses' columns never change), calls the
untouched kernel and the gain, and records in an archive environment the best fitness and best
graph genes seen at each node of each axis (plus a visit count). `GA::ga()` runs with the
existing population, crossover, Cauchy mutation and Nelder–Mead options; `popSize` follows the
rule in 3.9. The winner decodes to $(\mathbf w_{\text{GA}}, \mathbf G_{\text{GA}},
\boldsymbol\vartheta_{\text{GA}})$, is re-scored on the full sample at
$\boldsymbol\vartheta_{\text{GA}}$ exactly as `.graph_optimise_ga_gsd()` re-scores today.

**Stage L — local: graph refinement, then a warm-started profile along each axis.**
`.graph_optimise_local_gsd_opt()` draws `idx_local` once (`nsim_local`, default $N$) and
keeps that row order for everything below, so every gain in this stage is paired by trial.

1. COBYLA on the graph genes at $\boldsymbol\vartheta_{\text{GA}}$ from the GA winner, with
   the current local options; decode with `param_to_solution()` + `repair_graph()`; evaluate.
   This is `.graph_optimise_local_gsd()` on `resolve(spec, ϑ_GA)` **[code]**, and the oracle
   gate (5.3 P3 (a)) says so.
2. For each axis $r = 1..d$: for each node $q$ of $\Gamma_r$, walking outward from the current
   node in both directions, with $\boldsymbol\vartheta_{-r}$ held: evaluate three candidate
   graphs on that node's array — the previous node's refined optimum, the current overall
   optimum $(\mathbf w^*, \mathbf G^*)$ and the GA archive's best graph at that node (three
   kernel calls) — start COBYLA from the best, decode/repair, evaluate. Record the gain
   $\tilde P_r(q)$ (the best of the four values, so $\tilde P_r \ge S_r$ by construction),
   the graph, the fixed-graph value $S_r(q)$, and the batch-means paired difference and SE
   against the current optimum (1.4). Move to the argmax of the axis if it beats the current
   node by more than 0 (ties keep the current node).
3. If $d \ge 2$ and any axis moved, repeat the round (profiles of earlier axes are stale once a
   later axis moves); stop when a round moves nothing, or after 3 rounds (reported as
   `profile_converged = FALSE`). Axes whose stored profile is through a stale point are
   recomputed once at the end so that every reported profile passes through
   $\boldsymbol\vartheta^*$. For $d = 1$ this is one sweep.

Stage L's winner is $(\mathbf w^*, \mathbf G^*, \boldsymbol\vartheta^*)$ with its full-sample
gain. `choose_graph()` **[code]** then compares it with the GA winner at
$\boldsymbol\vartheta_{\text{GA}}$ on the same `ga_trial_success`/`local_trial_success`
fields; it is reused unchanged and the parameter vector rides alongside.

**Final.** `prune_graph_gsd()` and `calc_power_pvals_gsd()` run on
`resolve(spec, ϑ*)` with the user's `stopping_rule`, unchanged **[code]**.

**`global_search = FALSE`.** Stage L from the start graph at the start node; step 2 is the
whole parameter search. Step 1's oracle holds at the start node.

### 1.4 Paired standard errors

For two (graph, $\boldsymbol\vartheta$) pairs evaluated on the same rows, split the rows into
$B = 100$ contiguous blocks, call `trial_success$func()` on each block's time matrix for both
(200 calls on $N/100$-row matrices, microseconds each), and report the mean of the block
differences (equal to the full-sample difference up to summation order) and
$\mathrm{sd}/\sqrt{B}$. No change to `trial_success_gsd()`.

### 1.5 Reporting

New result fields, all absent on a fixed-path result so `print()`/`summary()` are unchanged
there: `spending_param` (named numeric, length $d$; names are the `id`s, else `gamma`,
`gamma1`, …), `hyp_spending_param` (length $m$, `NA` for fixed or short-circuited),
`spending` (length-$m$ labels, e.g. `sfLDOF`, `sfHSD(param = -2.25)`), `spending_bounds`
($d \times 2$), `spending_profile` (data frame: `param`, `value`, `gain` = $\tilde P$,
`gain_fixed_graph` = $S$, `diff` and `se` versus $\boldsymbol\vartheta^*$, `n_eval`, with the
per-node graphs as an attribute), `spending_search` (data frame from the GA archive: `param`,
`value`, `visits`, `best_fitness_subsample`), `profile_converged`.

`print()` adds, after the trial success function:

```
Spending functions:
  H1: sfLDOF
  H2: sfHSD(param = -2.25)   [optimised: gamma, bounds -4 to 1, step 0.25]
  H3: sfHSD(param = -2.25)   [optimised: gamma]
Gain profile in gamma (graph re-optimised at each value):
  within 5e-4 of the maximum on [-3.25, -1.25]; fixed-graph plateau [-2.5, -2]
```

`summary()` prints `spending_profile` per parameter (value, gain, fixed-graph gain, diff, se),
and the `spending_search` coverage line ("GA visited 21 of 21 nodes"). The two plateau
intervals are the connected run of nodes around $\vartheta^*_r$ with `diff <= 5e-4`
(the user testing's plateau convention, README 1, 7.4).

### Worked example

Unchanged from the candidate (two hypotheses, two looks, M1): $h_{-4}(0.5) \approx 0.119$,
$h_1(0.5) \approx 0.622$, Pocock-like spending for $H_2$ moved $w_1$ from 0.32 to 0.46 and the
gain from 0.6707 to 0.6903 **[measured, README 7.4–7.5]**. Under this design the profile table
would show that movement directly: the graph column of `spending_profile` changes along the
axis, and the fixed-graph column falls away from $\vartheta^*$ faster than the profile column.

---

## 2. Why It Is Valid

### 2.1 FWER control of the final design

As the candidate states (Maurer and Bretz 2013, condition (2)), with two sharpenings. (i) The
set of parameter values the search can try is exactly $\prod_r \Gamma_r$, every node's table
is built and checked in `transform_pvalues_gsd()` before the search exists, and
`resolve_pvalues_gsd()` refuses an off-lattice value, so there is no path by which an
unchecked $\boldsymbol\vartheta$ reaches the kernel. (ii) The final evaluation is
`calc_power_pvals_gsd()` on `resolve(spec, ϑ*)`, the same tables, so what is reported is what
was checked. The well-ordering tolerance is the fixed path's (`sqrt(.Machine$double.eps)`
times the column maximum **[code]**; record 10 item 9d).

Across $[-10, 2]$, what can be said without running code: $h_\vartheta(t)$ is exactly linear in
the level $a$ for HSD, so the spend increment at level $10^{-14}$ is $10^{-14} h$-fraction and
never underflows for $\vartheta \in [-10, 2]$ and $t \ge 0.1$ ($h_{-10}(0.1) =
(1-e)/(1-e^{10}) \approx 7.8\times10^{-5}$, giving an increment near $8\times10^{-19}$, an
ordinary double) **[derived]**; the first-look LDOF underflow of record 4.1 step 4 (a constant
$2.75\times10^{-89}$ below $\gamma \approx 10^{-6}$) came from LDOF's spend being
exponentially small in $a$, which HSD's linear spend is not **[derived]**. So the "fewer than
two distinct positive boundaries" abort is not expected anywhere in range. This is exactly the
claim `d01` must demonstrate rather than assume (record's own rule), at every node of
`seq(-10, 2, by = 0.25)` on the schedules $(0.5, 1)$, $(1/3, 2/3, 1)$, $(0.25, 0.5, 0.75, 1)$,
$(0.6, 1)$ and the corners $\vartheta = -10$ with $t = (0.1, 1)$ and $\vartheta = 2$ with
$t = (0.9, 1)$. `gsDesign::sfHSD()` at `param = 0` must also be checked (it documents linear
spending there; whether it accepts exactly 0 is to be confirmed) **[estimate]**.

### 2.2 Search approximations and Type I error

Unchanged: the search chooses a tuple; the tuple's error control is structural. One addition
the candidate omits: selecting $\boldsymbol\vartheta^*$ on the same draws that score it adds
an optimistic bias to the reported gain (winner's curse), as graph selection already does;
$dM$ extra candidates make it marginally larger. This is a power-reporting matter, not an
error-control one; the remedy is the user testing's fresh-draw check (README 3, G), which
`resolve_pvalues_gsd()` makes a two-line exercise.

### 2.3 Common random numbers

Raw draws fixed; `idx_ga` drawn once before the bank is built; `idx_local` drawn once before
step 1 of Stage L and reused by every node of every axis and round, in the same row order (the
default `nsim_local = N` is still a `sample.int(N, N)` permutation **[code]**, so order must be
pinned explicitly or block means are not paired). The GA archive values are on `idx_ga` and are
never compared with Stage L values; `spending_search` says so in its column name.

### 2.4 Piecewise-constant objective and the optimisers

As the candidate: COBYLA sees zero change inside a parameter cell, so it never receives
parameter genes; the GA tolerates the steps. The profile sweep is a derivative-free line search
on the lattice with an inner continuous refinement, which is the right tool for a step function
in $\vartheta_r$ whose inner problem is the one COBYLA already solves.

### 2.5 Fixed-graph slice, warm-started profile, true profile

$S_r(\vartheta_r) \le \tilde P_r(\vartheta_r) \le P_r(\vartheta_r)$, with equality of all three
at $\vartheta_r^*$ (the first because COBYLA at node $q$ starts from the best of three graphs
one of which is $(\mathbf w^*, \mathbf G^*)$ and reports its best evaluated point; the second
by definition of a maximum) **[derived]**. Hence plateau($S$) $\subseteq$ plateau($\tilde P$)
$\subseteq$ plateau($P$): the reported profile plateau is still a lower bound on the true
plateau, but a much tighter one than the slice, and the table shows both so the user can see
how much re-optimising the graph recovered. $\tilde P_r$ can fall short of $P_r$ only when the
optimal graph at a far node sits in a basin that none of the three starts reaches; the GA
archive start exists for that case, and `d02` measures the gap on the designs we have (section
4).

### 2.6 What is unaffected

The stopping-rule truth table depends on rejection sets only **[code]**; `alpha` is capped at
`spec$alpha` at every node because every table runs up to that level; `K` inference for
discount tables and `.gsd_check_gain_dims()` read `spec$K` and `spec$m`; the kernel is called
through the same two entry points with the same matrix layout.

---

## 3. Decisions

Each kept decision carries the sentence that keeps it; departures say what changed.

1. **Joint optimisation of $(\mathbf w, \mathbf G, \boldsymbol\vartheta)$ with a GA on a slice
   bank as the global stage.** KEPT. Survives because the binding requirement allows $d \ge 2$,
   where an outer product grid is $M^d$ full optimisations and a coordinate-only sweep can
   stall when two parameters must move together; the GA explores the joint space for about the
   cost of one extra generation set (3.9) and the same code serves every $d$.
2. **Pre-computed slices outside the kernel; kernel untouched.** KEPT. Survives because
   per-evaluation boundary work is seconds and the kernel is a binding no-touch; the bank is a
   one-off inversion from cached tables.
3. **Snapping to lattice nodes.** KEPT, and now the only scheme (blending removed, change 3).
   Survives because it is the only way every value tried is an exact, checked transform.
4. **Tagged constructor `sf_hsd_opt()` for arbitrary groups.** KEPT, with `step` added and
   the recycling rule changed (change 6). Survives because arity cannot distinguish an
   optimised parameter from `gsDesign`'s `param = NULL` formal **[code]**, and grouping needs
   an explicit key.
5. **Binary-fraction lattice with equal-probability gene mapping.** KEPT. Survives because
   `k * 2^{-j}` is exact, nodes coincide across parameters so tables are shared, and the
   mapping gives every node the same prior mass.
6. **Local stage.** CHANGED (change 1): graph-only COBYLA at $\boldsymbol\vartheta_{\text{GA}}$,
   then a warm-started re-optimised profile along each axis with axis rounds for $d \ge 2$,
   instead of a fixed-graph sweep plus one conditional COBYLA. Cost model and the decision
   threshold are in section 4.
7. **Per-winner full-sample re-scoring before `choose_graph()`.** KEPT. Survives because the
   two winners sit at different parameter values and each must be scored under its own.
8. **Flatness report.** CHANGED (change 1): the headline is the re-optimised profile with
   batch-means paired SEs; the fixed-graph slice is a second column of the same table; the GA
   archive is reported as coverage. Batch means survive because they need no change to the
   compiled gain.
9. **Standalone spec class and separate `_opt` files.** KEPT. Survives because inheriting
   would let `check_pvals_gsd()` pass raw p-values into the kernel through
   `calc_power_pvals_gsd()` **[code]**, and editing fixed-path closures risks the identity gate.
   Change 7 adds the eager tables to the spec; change 10 fixes the dispatch position.
10. **GA budget with parameter genes.** NEW: `popSize = min(max(40 * (m + d), 200), 500)`
    (a parameter gene counts as one more hypothesis in the existing rule), `run = 200`,
    `maxiter`, `pmutation`, `poptim` unchanged; mutation operator unchanged (change 5); the
    stopping rule `run` is self-adjusting, so generations rise only if the extra gene keeps
    producing improvements. Phase D `d04` measures generations and evaluations with and without
    the gene and may revise the 40.
11. **Result object.** CHANGED (change 8): no transformed array on the result;
    `resolve_pvalues_gsd()` is the way back to a `multigrain_pvals_gsd`.

---

## 4. What Is Not Yet Known (Phase D), and what was measured today

Numbers measured today are stated here with **[measured]**; Phase D re-runs them from scripts
committed under `dev/review/hsd_opt/` into `dev/review/hsd_opt_discovery.md` on
`gsd-build-hsc-opt`, and completes the rest.

### 4.1 What the existing measurements already settle

| Quantity | Value | Source |
| --- | --- | --- |
| One full transform, $m = 3$, $N = 10^5$, $K = 2$ | 3.2 s at $G = 1024$; 0.76 s at 256; 0.35 s at 128 | brief **[measured]** |
| Same at $K = 3$ | 4.6 s at $G = 1024$ | brief **[measured]** |
| Kernel call | about 13 ms at $N = 10^5$, $m = 3$, $K = 2$ | brief **[measured]** |
| Default run, $m = 2$ / $3$ / $4$, $N = 10^5$, 2 threads | 46–107 s / 97–154 s / 129 s | README 4.5 **[measured]** |
| Fixed-graph plateau width (within $5\times10^{-4}$) in $w_1$ | 0.08–0.24 wide | README 1, 7.4 **[measured]** |
| Optimiser's own spread of the gain at the optimum | about $4\times10^{-5}$ | README 4.2 **[measured]** |

Splitting the transform cost **[derived]**: the inversion does not depend on $G$ and the
table cost is linear in $G$, so from the $K = 2$ pair $3 \cdot 1024\,c + I = 3.2$ s and
$3 \cdot 128\,c + I = 0.35$ s: $c \approx 1.06$ ms per `gsBound1()` call and $I \approx 0$
(the 256 point predicts 0.81 s against 0.76 s measured). One boundary table therefore costs
about **1.1 s at $K = 2$ and 1.5 s at $K = 3$** for $G = 1024$ (half at 512), and the
inversion of a hypothesis on $10^5$ rows is tens of milliseconds. Consequences: the eager
tables of a default spec ($M = 21$, one schedule) cost about 23 s at $K = 2$ and 32 s at
$K = 3$, doubling with a second schedule or $M = 41$; `resolve_pvalues_gsd()` is sub-second;
the GA bank ($21 \times m_{\text{act}} \times K$ inversions on $5\times10^4$ rows) is a second
or two.

GA evaluations **[derived]**: `run = 200` stops the GA only after 200 generations without
improvement, so a run has at least $200 \times \texttt{popSize} = 4\times10^4$ fitness calls
plus Nelder–Mead steps. The measured 107 s at $m = 2$ then implies about 2.5 ms per call on
$5\times10^4$ rows with 2 threads, and 117–154 s at $m = 3$ about 3 ms; both are consistent
with 13 ms at $10^5$ rows on one thread. The GA is therefore about 90 % of today's run time,
COBYLA and pruning the rest (their counts are unmeasured; `d04`).

### 4.2 Cost model for the new path **[estimate]** ($m = 3$, $K = 2$, $N = 10^5$, 2 threads)

| Stage | Cost | Basis |
| --- | --- | --- |
| Spec tables, $d = 1$, $M = 21$ | 20–25 s once, reusable | 4.1 |
| GA on the bank | 1.2–1.5 × today's 2–2.5 min | one more gene; `run` self-adjusts; column copies $\approx 0.3$ ms per call ($5\times10^4 \times 2 \times 2 \times 8$ bytes) against a 3 ms kernel |
| Stage L step 1 (COBYLA at $\boldsymbol\vartheta_{\text{GA}}$) | as today | unchanged code |
| One profile node | inversion 50 ms + 3 kernel calls 20 ms + COBYLA $E_C \times 6.5$ ms + block means | $E_C$ unmeasured; COBYLA in 5 genes shrinking a trust region from 0.25 to $5\times10^{-8}$ on a step function needs of the order of 20 shrinks $\times$ 6 evaluations when flat, more when it moves: $E_C \approx 200$–$800$ gives 1.5–5.5 s per node |
| Profile, $d = 1$ | 21 nodes: 0.5–2 min | |
| Profile, $d = 2$ | $2 \times 21$ per round, 1–3 rounds: 1.5–12 min | the wide range is why `d04` must measure before the $d \ge 2$ default `step` is fixed |
| Total, $d = 1$ | about 3.5–6.5 min against 2–2.5 min today | |

Decision rule carried into Phase D: the profile sweep is the default if a 21-node $d = 1$
profile costs at most as much as the GA stage ($\lesssim 2.5$ min here); if it costs more,
the default profile uses every second node (`step` 0.5 for the sweep, 0.25 for the lattice)
and the full sweep is `control` opt-in. Either way the fixed-graph sweep is never the headline
report (section 2.5).

### 4.3 Remaining discovery items

1. **`d01_transform_split.R`** — the two measurement batteries already planned: (a)
   `.gsd_boundary_table()` and `.gsd_invert()` timed separately at $G \in \{128, 256, 512,
   1024\}$, $N \in \{10^4, 10^5\}$, $K \in \{2, 3\}$, with the coarse-grid interpolation error
   against $G = 1024$ at 200 random p-values for $\vartheta \in \{-4, -2, 1\}$; (b) the
   well-ordering and distinct-boundary sweep of 2.1 over $[-10, 2]$, plus the identity
   $h_\vartheta(t) = \texttt{sfHSD(1, t, }\vartheta\texttt{)\$spend}$ and exact proportionality
   in $a$ at $10^{-14}$ and $0.025$. Decide whether table building needs a progress step or a
   smaller default `grid_size` for specs (threshold: tables > 25 % of a default run at
   $m = 3$, $K = 2$, `nsim = 1e5`, $d = 1$).
2. **`d02_warm_start_profile.R`** — the measurement the local stage rests on: on Models M0
   and M1 (README 7.4) and the README section-2 design, at `nsim = 1e5`, compare
   $\tilde P(\vartheta)$ from the outward warm-started sweep against a full default
   `graph_optimise_gsd()` at each of 6 nodes spanning the bounds, 3 seeds. Gate for the
   design: $\max_q [\hat U_{\text{full}}(q) - \tilde P(q)] \le 5\times10^{-4}$ at every node
   within the reported plateau; outside the plateau report the gap. If the gate fails, the
   fallback is a second start per node from the GA archive at a coarser lattice, not the
   fixed-graph sweep.
3. **`d03_memory_and_copy.R`** — `tracemem()` on the working matrix across a kernel call and
   a column assignment; add a C++ column copier only if the R copy exceeds 15 % of kernel time
   at `nsim_global = 5e4` (today's number is in 4.1).
4. **`d04_eval_counts_and_budget.R`** — first the counts the brief asks for, on the README
   section-2 design ($m = 3$, $K = 2$, `nsim = 1e5`, `sfLDOF`): wrap `trial_success$func`
   with a counter (it is called exactly once per objective evaluation **[code]**) and run the
   internal drivers in the order `graph_optimise_gsd()` runs them (GA driver, local driver,
   `choose_graph()`, `prune_graph_gsd()`, `calc_power_pvals_gsd()`; the RNG is touched only
   by `.sample_pvals_gsd()` and the GA, so the staged run reproduces the one-call run), for
   seeds 1 and 2 and for `global_search = FALSE`; report GA generations (`@iter`), GA,
   COBYLA and pruning counts and times. Then warm-started COBYLA counts $E_C$ at
   $\vartheta \in \{-2.25, -1.75, -3, -1, -4, 1\}$ from the $\vartheta = -2$ optimum (this
   feeds 4.2). Then prototype the bank closure in a scratch script; run the GA with and
   without the parameter gene on M0/M1 at `nsim = 1e5`, 5 seeds each; report generations,
   evaluations, wall time, and the best node found against the `d02` profile argmax; test
   `popSize` 200 against the 3.10 rule. Also the `tracemem()`/column-copy timing of `d03`
   can share this script's setup.
5. **`d05_grid_step_loss.R`** — at the `d02` optimum, warm-started COBYLA on a fine grid
   ($\Delta = 0.025$) within $\pm 0.5$; the coarsest $\Delta \in \{0.5, 0.25, 0.125\}$ whose
   loss is within the optimiser's own spread ($4\times10^{-5}$ **[measured, README 4.2]**)
   becomes the default; otherwise 0.25 stands.

---

## 5. Implementation Notes

### 5.1 File map and touchpoints

New files:

- `R/transform_pvalues_gsd_opt.R`: `sf_hsd_opt()`, `is_sf_hsd_opt()`,
  `.has_hsd_opt()`, `.transform_pvalues_gsd_spec()`, `new_pvals_gsd_spec()`,
  `print.multigrain_pvals_gsd_spec()`, `.hsd_groups()`, `.hsd_lattice()`,
  `.hsd_node_tables()`, `.gsd_invert_hyp()` (duplicate of steps 4–7 of
  `.gsd_transform_hyp()` on a prebuilt table — "share or duplicate: duplicate"),
  `resolve_pvalues_gsd()` (exported), `.hsd_labels()`.
- `R/objective_function_gsd_opt.R`: `create_obj_func_gsd_opt()` (bank, clamp, map, column
  copy, archive), `.hsd_bank()`, `.check_bank_memory()`.
- `R/optimisation_gsd_opt.R`: `.graph_optimise_gsd_opt()`, `.graph_optimise_ga_gsd_opt()`,
  `.graph_optimise_local_gsd_opt()`, `.hsd_profile_axis()`, `.hsd_refine_node()`,
  `.hsd_block_se()`, `.hsd_plateau()`.
- Tests: `test-transform_pvalues_gsd_opt.R`, `test-objective_function_gsd_opt.R`,
  `test-optimisation_gsd_opt.R`, `test-graph_optimal_hsd_opt.R`; dev benchmarks under
  `dev/review/hsd_opt/`.

Untouched: `src/*`, `R/objective_function_gsd.R`, `R/choose_graph.R`,
`R/post_optim_processing_gsd.R`, `R/calc_power_gsd.R`, `dev/gsd_identity_check.R`.

Guarded touchpoints on existing files (each a few lines, each behind a predicate that is
`FALSE` on the fixed path):

1. `R/transform_pvalues_gsd.R`: after `.gsd_spending_list()`, `if (.has_hsd_opt(spending))
   return(.transform_pvalues_gsd_spec(...))`. (`.gsd_spending_list()` accepts the tagged
   functions as functions already **[code]**; `.gsd_check_spending()` is then run per node
   inside the spec builder.)
2. `R/optimisation_gsd.R`: first statement of `graph_optimise_gsd()`:
   `if (inherits(pvals, "multigrain_pvals_gsd_spec")) return(.graph_optimise_gsd_opt(...))`.
3. `R/check_gsd.R`: `check_pvals_gsd()` abort branch gains a hint for a spec
   ("call resolve_pvalues_gsd() or pass the spec to graph_optimise_gsd()"). The dimension
   helpers already read `$nsim`, `$m`, `$K`, `$alpha` and work on the spec **[code]**.
4. `R/graph_optimal.R`: `graph_optimal()` accepts `spending = NULL` list of the new fields
   and attaches them only when non-`NULL` (as `stopping_rule` is handled **[code]**);
   `print()`/`summary()` blocks guarded by `!is.null(x$spending_param)`; snapshot
   `tests/testthat/_snaps/graph_optimal.md` byte-identical.

### 5.2 Mechanics worth stating

- **RNG order.** Spec construction draws nothing. Stage G: `idx_ga`, then `GA::ga()`. Stage L:
  `idx_local`, then COBYLA, then the sweeps (no draws). Reproducible under one seed; identity
  with the fixed path is required only at $d = 0$, where the handoff returns before any of this.
- **Clamping.** `u <- pmin(pmax(u, 0), 1)` before the node map; the GA's winner is clamped the
  same way.
- **Archive.** An environment in the closure: per axis a length-$M_r$ numeric `best`, integer
  `visits`, and a $M_r \times p_{\text{graph}}$ matrix `genes`; updated when
  `fitness > best[q]`. Cost: a handful of vector ops per evaluation.
- **Warm start encoding.** A decoded graph $(\mathbf w, \mathbf G)$ is encoded for COBYLA with
  `create_start_params(graph_constraint, w0 = proj$hyp_weight, G0 = proj$trans_matrix,
  sum_to_one_constraint = FALSE)` after `proj <- closest_graph_to_constraints(graph_constraint,
  w, G)` — the user-row branch of `.build_start_matrix()` **[code]** — and never with
  `.build_start_matrix(...)[1, ]`, whose first row is the Holm seed (change 13). The P3
  oracle (a) includes `param_to_solution(x0) == (w, G)` for the encoded start. Node arrays come
  from `resolve_pvalues_gsd()` subset by `idx_local` with `drop = FALSE`. The GA's suggestion
  matrix is `.build_start_matrix()` (all rows, as today) with the start node's gene
  $(q_{\text{start}} - 0.5)/M_r$ appended to every row.
- **Mutation on the parameter gene.** `.cauchy_perturb()` with scale 1 on $[0,1]$ accepts a
  proposal within $\pm 0.05$ (one cell at $M = 21$) about 11 % of the time and otherwise jumps
  almost uniformly across the range **[derived]** ($P(|X| < 0.05) / P(|X| < 0.5)$ for a
  standard Cauchy centred mid-range); fine-grained movement in $\vartheta$ comes from
  Nelder–Mead and from Stage L, so no second operator is needed.
- **Memory guard.** Bank bytes $= 8 K N_{\text{GA}} \sum_{i \in \mathcal I_{\text{HSD}}}
  M_{g(i)}$; warn above 2 GB, abort above `max_bank_gb = 4` with the three remedies
  (coarser `step`, narrower `bounds`, smaller `nsim_global`). Defaults never approach it
  (appendix A.2).
- **Eval budget for profile COBYLA.** Same `local_opt` as the main pass by default (so the
  oracle holds); `d04` decides whether a profile-specific `maxeval` is worth a control option.

### 5.3 Phases, gates, orchestration

Orchestration and house rules as in the candidate (one agent per phase reporting rather than
committing; every gate from a clean build; `test_file()` only, `NOT_CRAN=true`;
$m \le 4$, `nsim <= 1e4` in tests; `LC_ALL=C.UTF-8 devtools::document()`; `lintr` on touched
files; revert CRLF snapshot churn; `dev/gsd_identity_check.R` at every gate).

**Phase D.** Scripts `d01`–`d05` of 4.3 (the two batteries already planned by the measurement
agents run first and replace the derivations of 4.1 with `[measured]` numbers); report with
tables for every item in 4.1–4.3 and a decision on $\Delta$, the default profile density,
`popSize` and the warm-start gate; also record the `start_graph` finding (change 13) in the
report's pre-existing-issues list. *Stop and wait for the user's decision.*

**P1 — spec class, lattice, eager tables, resolve.** Gate: bounds/step/start/group validation
messages; single recycled tag = one parameter; `d = 0` handoff returns a
`multigrain_pvals_gsd` `identical()` to the fixed path's for the same inputs (both an untagged
spec and an all-short-circuited one); table deduplication `identical()` to per-hypothesis
tables; `resolve_pvalues_gsd(spec, ϑ)` `identical()` to `transform_pvalues_gsd()` with
`function(a, t) gsDesign::sfHSD(a, t, ϑ)$spend` at every node of a 5-node lattice ($m = 3$,
$K = 2$, `nsim = 2000`, one fixed, two tagged in one group); off-lattice value refused;
`check_pvals_gsd()` refuses the spec with the hint; identity check passes.

**P2 — objective closure and bank.** Gate (node oracle): for 20 random graph gene vectors and
every node, `create_obj_func_gsd_opt()(c(x, u))` is `identical()` to `create_obj_func_gsd()`
built on `resolve(spec, node)` restricted to the same `idx_ga` rows; genes outside $[0,1]$
decode to the end nodes; archive `best`/`visits`/`genes` match a by-hand replay; copy overhead
and `tracemem()` meet the `d03` threshold; identity check passes.

**P3 — driver.** Gate: (a) warm-start oracle: `.hsd_refine_node()` from a given `x0` on
`resolve(spec, node)` rows `idx_local` is `identical()` to `.graph_optimise_local_gsd()`'s
`nloptr` result on the same array and `x0` with the same seed; (b) final-evaluation identity
`identical(res$power, calc_power_pvals_gsd(resolve_pvalues_gsd(spec, res$spending_param),
res$hyp_weight, res$trans_matrix, stopping_rule = ...))`; (c) same seed twice → `identical()`
on every component except `trial_success`; (d) field assertions at $m = 4$, $d = 2$
(`length(spending_param) == 2`, values on the lattice within bounds, `NA`s in
`hyp_spending_param` for fixed and short-circuited hypotheses, profile rows $= M_1 + M_2$,
`diff == 0` and `se == 0` at $\boldsymbol\vartheta^*$, `profile >= fixed_graph` at every
row); (e) `global_search = FALSE` runs; (f) stopping rule through the whole pipeline with
`stop_distribution` present; (g) dev benchmarks (`gate_p3_benchmarks.R`): joint gain
$\ge \max_q \hat U_{\text{outer}}(q) - 5\times10^{-4}$ at $d = 1$ (`nsim = 1e5`, 21 nodes)
and at $d = 2$ coarse (`nsim = 2e4`, $5 \times 5$); joint gain $\ge$ fixed-spending gain at
the start node $- 5\times10^{-4}$; run time at $m = 3$, $K = 2$, `nsim = 1e5`, $d = 1$ against
the 4.2 projection.

**P4 — reporting.** Gate: `print()`/`summary()` snapshots for a $d = 1$ and a $d = 2$ result;
plateau intervals computed from a hand-built profile; `test-graph_optimal.R` and the
fixed-path snapshots byte-identical; identity check passes.

**P5 — documentation and record.** roxygen, `document()`, pkgdown entry, article section
(the README 7.5 pattern as the running example), `NEWS.md`, record section 4.12 with the
discovery numbers, call-flow document, `dev/gsd_progress.md` entry. Gate: clean build, `lintr`
clean, affected `test_file()` runs green, identity check green.

---

## Appendix

### A.1 Edge-case rules

As the candidate's A.1 items 1, 2 and 4, with: item 3 replaced by change 6 (a single recycled
`sf_hsd_opt()` is one shared parameter; a list element with `id = NULL` is its own); item 5
(extended tables for blending) removed; new: `resolve_pvalues_gsd()` refuses off-lattice
values; a tagged hypothesis whose group has `step` such that the lattice has fewer than 2
nodes cannot occur (both bounds are always nodes); a node table failing well-ordering aborts
the spec with the node value and the hypothesis.

### A.2 Bank memory

Only Stage G holds a bank, so the candidate's formula and table stand with
$N_{\text{GA}} = \min(5\times10^4, N)$: about 200 MB at $K = 3$, $m_{\text{act}} = 4$,
$M = 41$; 470 MB at $M = 97$ **[estimate]**. Stage L holds one $N \times mK$ array per node at a
time (released after the node) plus the two or three graphs it compares.
