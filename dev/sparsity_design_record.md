# Design record: `graph_simplify()`, budgeted edge-count reduction of an optimised graph

Status: design specification, not yet implemented. Every design element is stated, followed by the reasoning behind it and the alternatives it was chosen over. Elements are marked **LOCKED** when the reasoning is considered settled and **OPEN** when evidence is still needed; open items are collected in section 11. Claims about existing package behaviour were checked by running the code; the scripts and their output are in the appendix.

Notation: $U(\mathbf w,\mathbf G)$ is the trial-success measure (expected gain) of a graph, estimated on a fixed matrix of simulated p-values; $\psi$ is the gain function compiled by `trial_success()`; $\mathbf G$ is the transition matrix; $E(\mathbf G)$ is the number of non-zero free entries of $\mathbf G$; $m$ is the number of hypotheses; $\lambda$ is the user's cap.

## 1. Problem

Optimised graphs frequently come back with several small non-zero transition weights. Each contributes almost nothing to expected gain, but together they make the graph hard to present: a clinical reviewer has to reason about every arrow. The package already has a post-hoc step, `prune_graph()`, that removes a weight or edge when doing so does not lower the estimated gain. That step never trades gain for simplicity.

The feature specified here is a second optimisation applied to an already-optimised graph. The user states how much expected gain they are prepared to give up in total, as a fraction $\lambda$ of the gain of the graph they already have, and receives a graph with fewer edges whose loss is guaranteed to stay within that amount on the same simulated trials. Simplicity is measured by the number of edges. Family-wise error control must be untouched, and `graph_optimise()` must behave exactly as it does today.

## 2. How the package optimises today

The facts below constrain the design and were verified by reading and running the code.

`create_obj_func()` in `R/objective_function.R` is the single objective factory. Its closure decodes the parameter vector through `split_theta()`, `recover_full_weights()` and `recover_full_trans_matrix()`, returns a negative penalty for `NA` or out-of-range entries, zeroes hypothesis weights below $10^{-4}$ and transition entries below $10^{-5}$ without renormalising, runs `graph_shortcut()` (or its parallel twin) and returns the gain function's mean. The same factory serves `GA::ga()` for the global search, including its intermittent Nelder-Mead step (`optim = TRUE`), and the final `nloptr` COBYLA run.

The parameter encoding gives each row of $\mathbf G$ with $k_i$ free entries $k_i-1$ parameters; the last free entry is derived as one minus the rest. Rows therefore always sum to one and each row keeps at least one edge. `graph_constraint()` rejects a row with exactly one free entry and rejects an incomplete row whose fixed entries already sum to one (appendix check 1), so every incomplete row has $k_i\ge2$. Re-encoding a graph that has an exact zero in a derived position produced exactly zero in every case tried (appendix check 9): `sum()` accumulates in extended precision, so one minus a sum that is one in exact arithmetic rounds to zero rather than to $-10^{-17}$. Deliberate overshoots do give negative entries and a penalty (check 2).

`param_to_solution(process = TRUE)` snaps entries below $10^{-5}$ to exact zero and entries in $(10^{-5},10^{-3})$ to the $\epsilon$-edge value $0.001$, then renormalises each row with `normalise_sum()`. `prune_graph()` runs `prune_hyp_weights()` and then `prune_edges()`, each in a fixed index order, and `.try_prune()` accepts a removal if the gain does not fall. `.redistribute_mass()` moves a dropped entry's mass proportionally onto the row's other free entries and, when those are all zero, spreads it uniformly, which creates $m-2$ new edges.

`graph_optimise()` runs the GA on a random subsample of `nsim_global` trials, COBYLA on a separate random subsample of `nsim_local` trials, chooses between the two results on the full sample with `choose_graph()`, prunes on the full sample, and only then computes the reported `$power` from the returned graph. That last ordering is deliberate: version 0.2.0 fixed a bug in which `$power$trial_success` described the graph before pruning rather than the one the user saw.

The returned `multigrain_graph_optimal` stores the constraint, the trial-success object, the prepared control object and the raw GA and nloptr results. `graph_optimal()` strips only the `@call` slot of the GA object, so the GA's final population (`global_output@population`, `popSize` rows by the encoding's parameter count) is available afterwards (check 8). The object does not store `pvals` (large, by design) and does not store `alpha`. `GA::ga()` requires `ncol(suggestions)` to equal the number of variables and fails with "number of items to replace is not a multiple of replacement length" when `suggestions` has more rows than `popSize` (check 8).

The mutation used inside the GA is a package closure (`.make_cauchy_mutation_multi()` in `R/mutation_helpers.R`), and the crossover is `GA`'s default local arithmetic crossover.

## 3. Specification

### 3.1 User-facing behaviour

A new exported function:

```r
graph_simplify(
    graph_optimal,
    pvals,
    ...,
    gain_tolerance = 1e-3,
    alpha = NULL,
    global_search = TRUE,
    num_threads = 1L,
    control = NULL,
    verbose = multigrain_verbosity()
)
```

`graph_optimal` is the result of `graph_optimise()`; it is the reference graph. `pvals` is the matrix of simulated p-values, supplied again because the object does not store it; the cap is measured on the supplied matrix. `gain_tolerance` is $\lambda$: the returned graph's gain on the full sample is at least $(1-\lambda)$ times the reference graph's gain on the same sample. `alpha` defaults to the value stored on the object (see 3.10); an explicit value overrides it. `control` defaults to the object's stored control with the GA's `run` halved (see 3.5). The function returns a `multigrain_graph_optimal` so that `print()`, `summary()`, `plot()` and `calc_power_pvals()` work unchanged, with a populated `sparsity` element (3.9) and `solution$opt_source` set to `"simplify:global"`, `"simplify:local"` or `"reference"`.

`graph_optimise()` changes in two ways only: its result stores `alpha`, and its result carries a `NULL` `sparsity` element. Its search, its random-number consumption and every other element of its result are unchanged.

**Rationale.** LOCKED. Reducing an optimised graph is a second optimisation with its own inputs (a reference graph and a cap), so it is its own verb. Making it a separate function keeps `graph_optimise()` untouched, which satisfies "unchanged when not used" by construction rather than by a fixture; makes the reference explicit as an object the user can plot next to the simplified one; lets the user try several values of $\lambda$ without repeating the first optimisation; and lets the stored control object and stored GA population serve as defaults. The name follows the package's `graph_*` verbs (`graph_optimise`, `graph_random`, `graph_constraint`) and states the user's intent. Alternatives: an argument `gain_tolerance` on `graph_optimise()` (hides the two-stage nature, forces a `NULL`-means-off default and a fixture-based identity gate; rejected); `graph_sparsify()` (precise but jargon); `graph_prune()` (collides with the internal pruning step).

### 3.2 The cap: a threshold against the reference, recomputed on every sample

Let $x_{\text{ref}}$ be the reference graph and $U_{\text{ref}}$ its gain. The cap is the constraint $U(x)\ge T$ with $T=(1-\lambda)\,U_{\text{ref}}$. Every objective closure computes $T$ itself, at construction, by running the shortcut on the reference graph over the p-value sample that closure has captured. The GA closure, the COBYLA closure and the full-sample evaluations therefore each have their own $T$, all defined relative to the same reference on the same trials as the candidates they judge.

**Rationale.** LOCKED. A subsample of 50 000 trials estimates $U$ with a standard error near $4\times10^{-3}$, larger than a budget of $10^{-3}$; a single absolute threshold carried across samples would be meaningless. Judged against the reference on the same trials, the comparison is paired: only trials whose outcome differs between the two graphs contribute. Appendix check 4 measured this at $m=4$, $n_{\text{sim}}=10^4$: removing an $\epsilon$-edge changed 59 of 10 000 trials, and the paired standard error of the difference was $1.9\times10^{-4}$ against a marginal standard error of $3.4\times10^{-3}$, a factor of 18. The reference is feasible on every sample by construction ($U_{\text{ref}}\ge(1-\lambda)U_{\text{ref}}$ whenever $U_{\text{ref}}\ge0$), so every search starts with a feasible seed and the final fallback (3.8) always has a feasible graph to return. Alternatives: a per-edge price derived from the cap (section 6, alternative 1); a threshold on the full sample only (would make the GA's and COBYLA's subsample rankings inconsistent with the cap).

### 3.3 Objective

With $U_{\max}$ and $U_{\min}$ the exact extremes of $\psi$ over all $2^m$ rejection patterns (computed once from the compiled function; check 6), $E_{\max}$ the number of free off-diagonal entries, and $D=(U_{\max}-T)+1$:
$$
f(x)=\begin{cases}
U(x)-D\,E(x) & \text{if } U(x)\ge T,\\[2pt]
U(x)-D\,(E_{\max}+1) & \text{if } U(x)<T.
\end{cases}
$$
Invalid encodings return `penalty + U_min - D * (E_max + 2)`, with `penalty` the existing negative violation measure.

Reading the score: among feasible graphs, one fewer edge always wins (because $D$ exceeds the largest possible difference in $U$ between two feasible graphs), and among graphs with the same edge count the higher gain wins. Every infeasible graph sits below every feasible one, but infeasible graphs still rise with $U$, so an optimiser that finds itself infeasible is pulled back toward the boundary. Invalid encodings sit below everything. The scalar rule is one internal helper, `.lexico()`, shared by the objective closure, `choose_graph()` and pruning, so no stage can prefer a graph that another would reject.

**Rationale.** LOCKED. The problem is "fewest edges subject to a floor on gain", and a lexicographic score is its direct scalarisation. Check 11 shows the two branches doing their jobs: from the encoded reference, COBYLA kept all seven edges and raised $U$ slightly; from a start that was infeasible at the chosen $\lambda$, COBYLA re-added one edge and restored feasibility. Alternatives: a per-edge price (section 6); a penalty large enough to dominate without a threshold (then nothing stops it removing everything; a threshold in disguise).

### 3.4 What counts as an edge

A free entry counts as an edge if it is at least $10^{-5}$ after decoding, evaluated on the same thresholded matrix the shortcut receives. Every non-zero free entry counts one, whatever its size; $\epsilon$-edges of $0.001$ count. Entries pinned by `graph_constraint()` are not counted in $E$; the reported edge count of the returned graph includes them, with the free count reported alongside.

**Rationale.** LOCKED. The threshold is the one the objective already applies, so an entry that counted during the search is still an edge after `param_to_solution()` snaps it, and an entry that did not count is exactly zero in the returned graph. Zeros in the returned graph are exact: check 5 shows a parameter set to exactly zero staying exactly zero through `param_to_solution(process = TRUE)`, `repair_graph()` and `normalise_sum()`. Counting by size would reward shrinking edges rather than removing them, which does not help a reader. Pinned entries are constants across every candidate, so counting them adds a constant to $E$.

### 3.5 Warm start: seeds, encoding guard, and the GA budget

The stage-2 GA population is seeded, in this order, with: the encoded reference graph; the encoded single-edge-removal neighbours of the reference (each free non-zero edge dropped with `.redistribute_mass()`, skipping candidates that do not reduce the edge count or leave an invalid row); the usual seeds from `.build_start_matrix()` (uniform, fixed sequence, any user-supplied start graphs); and the rows of the reference object's `global_output@population` when present and of matching width. Duplicate rows are removed and the matrix is truncated to `popSize` rows.

Encoding uses `create_start_params(gc, w0, G0, sum_to_one_constraint = FALSE)` with a guard: if a row's derived entry decodes below $10^{-5}$, that row's parameters are rescaled to sum to $1-5\times10^{-6}$ minus the row's fixed entries; if the derived hypothesis weight decodes below $10^{-4}$, the weight parameters are rescaled to sum to $1-5\times10^{-5}$ minus the fixed weights. The derived entry then decodes to exactly $5\times10^{-6}$ (or $5\times10^{-5}$), positive and below its threshold, so it is not penalised and not counted (check 9: seven edges before and after, identical objective value).

The GA runs with the reference object's stored control, `run` halved (floor 1), re-prepared against the supplied `pvals` so that `nsim_global` and `nsim_local` are clamped to its row count. A control passed explicitly is used as is.

**Rationale.** LOCKED. The reference is already near the optimum, so the population should start there: the reference guarantees a feasible individual, its neighbours give the GA one-edge-removal candidates in generation zero, the stage-1 population gives diversity in the right basin, and the fixed-sequence seed is the sparsest graph there is. Truncation to `popSize` is required because `GA::ga()` errors on oversized `suggestions` (check 8); the ordering puts the most valuable seeds first. The encoding guard costs nothing and removes a failure mode that would otherwise depend on floating-point accident; check 9 found no negative derived entries from plain re-encoding, so the guard is defensive rather than a fix for a common failure. Halving `run` reflects that stage 2 starts in the right basin and that its improvements are mostly discrete (an edge removed); the stage-1 `run` of 200 no-improvement generations is a floor on cost that a warm start does not need. Alternatives: identical stage-1 settings (simplest to explain, roughly doubles total time); GA off by default (leaves topology changes to pruning alone, which cannot find removals that only become affordable after weights are re-tuned).

### 3.6 A support-changing mutation move

**What the mutation does today.** `.make_cauchy_mutation_multi()` returns a closure the GA calls on one parent per mutation. It selects each parameter with probability $0.1$ (at least one), and replaces each selected value with a draw from a Cauchy distribution centred on the current value and truncated to $[0,1]$. That is a good operator for moving weight around continuously. It is a poor operator for removing an edge, because "removed" means the decoded entry falls below $10^{-5}$, and a Cauchy draw lands in $[0,10^{-5})$ about once in $100\,000$ tries (check 3). Crossover does not help: `GA`'s default real-valued crossover takes convex combinations of two parents, so a child entry is zero only where both parents are already zero. Nelder-Mead and COBYLA are continuous methods and never aim for the threshold. With the lexicographic score the GA would therefore *prefer* sparser graphs but almost never *produce* one to compare, and the score would be inert during the global search apart from the seeds.

**What is added.** In `graph_simplify()`, each mutation call first draws a coin with probability $p_0$. If it comes up, the call performs a zeroing move instead of a Cauchy perturbation, of one of two kinds:

1. *Zero a parameter.* Pick, uniformly, one free transition parameter currently at or above $10^{-5}$ and set it to exactly zero. That removes the corresponding edge outright.
2. *Zero a derived entry.* Each row of the transition matrix has one entry that is not a parameter: it is computed as one minus the sum of the row's parameters. Setting a parameter to zero cannot remove that edge. The move instead picks a row and rescales its parameters so they sum to $1-5\times10^{-6}$. The derived entry then equals $5\times10^{-6}$: positive, so it is not penalised as a negative entry, and below $10^{-5}$, so it does not count as an edge and is snapped to exact zero on output (check 2).

If the coin does not come up, or there is nothing to zero (no parameter above the threshold, or a row whose parameters are all zero), the call falls through to the ordinary Cauchy perturbation.

**Rationale.** LOCKED. Without the move the score would have no effect on the global search beyond ranking the seeds. Without the derived-entry variant one edge per row, $m$ of the $m(m-1)$ edges in an unconstrained graph, would be unreachable by the global search and removable only by pruning at the end. The factory keeps the current closure verbatim as an inner function and returns it unchanged when $p_0=0$, which is what `graph_optimise()` requests, so `graph_optimise()` runs the same code and consumes the same random numbers as today. The value of $p_0$ (proposed $0.2$) and the equal split between the two variants are open item O1.

**Documentation requirement.** The two paragraphs above ("What the mutation does today", "What is added") are to be carried into the package when this is implemented, in three places: the comment block of `.make_cauchy_mutation_multi()` in `R/mutation_helpers.R`; a `@details` section headed "Edge removal during the global search" in the roxygen for `graph_simplify()`, written for a user who wants to know how the optimiser can reach a sparser graph; and the corresponding paragraph of the get-started article under `vignettes/articles/`. Implementation steps 2 and 7 in section 7 carry this as part of their gates.

### 3.7 Local search and choice

After the GA (or directly from the encoded reference when `global_search = FALSE`), COBYLA runs with the lexicographic objective on its own subsample. Inside a fixed support it re-tunes the weights exactly as today; it does not create an edge, because an entry crossing $10^{-5}$ costs $D$; and if its starting point is infeasible on that subsample it climbs in $U$ and may re-add an edge to become feasible (check 11). `choose_graph()` then compares the GA and COBYLA results on the full sample by the lexicographic score; the two internal optimisers return that score as `ga_objective` and `local_objective`, which equal the raw trial success when no threshold is set, so `choose_graph()`'s decisions in `graph_optimise()` are unchanged.

**Rationale.** LOCKED. The local search's job is unchanged (polish within a support); the score only prevents it from undoing the GA's removals. Comparing on the full sample with the same score keeps the choice consistent with the cap.

### 3.8 Best-first pruning and the reference fallback

`prune_hyp_weights()` runs unchanged: a weight is removed only if $U$ does not fall. `prune_edges()`, when a threshold is supplied, becomes best-first on the full sample: evaluate every remaining removable edge, keep the candidates that are feasible ($U\ge T$) and reduce the edge count, remove the one with the highest $U$, repeat until no candidate is feasible. Candidates that would leave an invalid row (no free recipient) are skipped. The exact loss across accepted removals is recorded as `prune_loss`. Finally the pruned graph is compared with the reference under the lexicographic score on the full sample; whichever scores higher is returned, and if it is the reference, `solution$opt_source` says so.

**Rationale.** LOCKED. Pruning is the one deterministic support move on the full sample; using the same score makes it the final acceptance test of the whole pipeline. Best-first spends the budget on the cheapest edges; check 10 shows it removing two edges at $\lambda=0.02$ for 61% of the budget and stopping when the next cheapest edge would exceed it. Check 7 shows best-first at zero price reaching the same graph as the current fixed-order prune, so the change of order is safe. The fallback guarantees the cap on the full sample whatever happened on the subsamples: the reference is feasible there by construction, and the returned graph is never worse than it under the score. Alternatives: fixed-order pruning with the threshold (can spend the budget on the first edge it visits); no fallback (a graph feasible on the GA subsample but marginally infeasible on the full sample could be returned).

### 3.9 Reporting

The `sparsity` element of the returned object:

| element | meaning |
|---|---|
| `gain_tolerance` | $\lambda$ as supplied |
| `reference` | list with the reference `hyp_weight` and `trans_matrix` |
| `gain_reference` | $U_{\text{ref}}$ on the supplied `pvals` |
| `gain` | $U$ of the returned graph, equal to `$power$trial_success` |
| `gain_loss` | `gain_reference - gain` |
| `gain_loss_fraction` | `gain_loss / gain_reference` |
| `budget` | $\lambda\,U_{\text{ref}}$ |
| `n_edges_reference`, `n_edges` | total non-zero entries, pinned included, before and after |
| `n_edges_free_reference`, `n_edges_free` | non-zero free entries before and after |
| `prune_loss` | exact loss of $U$ across the removals pruning accepted |
| `source` | `"global"`, `"local"` or `"reference"` |

`$power` remains the last quantity computed, from the returned `hyp_weight` and `trans_matrix` on the full sample. `print()` and `summary()` add two lines when `sparsity` is non-`NULL`, of the form:

```
Simplified from 12 edges to 7 (free: 12 -> 7)
Trial success 0.8048 -> 0.8041: loss 0.09% of reference (cap 0.1%)
```

**Rationale.** LOCKED. The reader must be able to see the unadjusted gain, the edge count, and what the simplification cost; here all three are exact on the same trials, and the reference graph travels with the result so the two can be plotted side by side. Computing `$power` last is the 0.2.0 fix and keeps the reported gain equal to the graph the user sees.

### 3.10 Changes to `graph_optimise()` and its result

The constructor `graph_optimal()` gains `alpha = NULL` and `sparsity = NULL`; `graph_optimise()` passes its `alpha` and leaves `sparsity` as `NULL`. `names()` of an optimised graph therefore gains two entries; nothing else changes. `graph_simplify()` reads `alpha` from the object and refuses to run without it unless an explicit `alpha` is passed, so an object saved by an earlier version still works.

**Rationale.** LOCKED. The significance level is part of what was optimised, and a user who forgets to pass it again would silently simplify against the wrong reference gain. Storing it is one element on the object. `pvals` is not stored because it is large and the user has it.

### 3.11 Error control and degenerate cases

Family-wise error control is unaffected. A Bonferroni-based graphical procedure controls the family-wise error rate strongly for any weights summing to at most one and any transition matrix with non-negative entries, zero diagonal and row sums at most one. Zeroing an entry and renormalising its row preserves all three, and `is_graph_valid()` still gates `calc_power_pvals()`. Only the objective changes, not the class of procedures searched.

Cases the implementation must handle:

1. **$U_{\text{ref}}\le0$ on the full sample.** The cap has no meaning; abort with an informative error.
2. **Nothing removable** ($m=2$, or every off-diagonal entry pinned). Warn, skip the search, return the reference with `source = "reference"` and a populated `sparsity`.
3. **$\lambda=0$.** Well defined: fewest edges among graphs whose gain does not fall below the reference on each sample. Cheap and useful; documented as supported.
4. **$\lambda=1$.** Every non-negative graph is feasible, so the result has one free edge per row. Useful as a test invariant.
5. **A different `pvals` than stage 1.** The cap is measured on the supplied matrix; `gain_reference` is recomputed there and may differ from the stored `$power$trial_success`. Documented; no warning.
6. **Stored GA population absent or of the wrong width** (stage 1 ran with `global_search = FALSE`, or the constraint changed). Seeds proceed without it.
7. **Gain functions with negative values.** Invalid-encoding penalties are shifted below $U_{\min}-D(E_{\max}+1)$ so they stay below every valid score.
8. **A row reduced to one non-zero free entry during pruning.** The uniform fallback raises $E$; the candidate is skipped. The path in `.redistribute_mass()` with no recipients, which returns a row summing to less than one and makes `calc_power_pvals()` abort, is reachable only by bypassing `graph_constraint()` (check 1).
9. **Sparse two-cycles** $g_{ij}=g_{ji}=1$ hit the shortcut's `denom == 0` branch in `src/graph_shortcut.cpp`, which zeroes the row; that is the correct limit of the Bretz update and already exercised by every $m=2$ graph.
10. **Unreachable nodes** (zero weight, no incoming edge). Their outgoing edges do not affect $U$; the search reduces them to one edge.

## 4. What to expect

### 4.1 Budget spent, and lumpy

Check 10 on a four-hypothesis reference with seven edges and $U_{\text{ref}}=0.8048$: at $\lambda=10^{-3}$ one edge went for $7.5\times10^{-5}$ (9% of the budget); at $\lambda=10^{-2}$ the same, because the next cheapest edge cost $9.7\times10^{-3}$, more than the remaining budget; at $\lambda=2\times10^{-2}$ two edges went for 61% of the budget. Edge costs are discrete, so the fraction of the budget used varies; the cap is the guarantee, the spend is whatever the graph allows. A per-edge price of $1.25\times10^{-3}$, the alternative in section 6, could never have removed the second edge at any of these $\lambda$.

### 4.2 Runtime

Stage 2 repeats the global search with a warm start and half the no-improvement generations, then one COBYLA run and a best-first prune of order $E^2$ shortcut evaluations. When called immediately after `graph_optimise()` the total is roughly 1.3 to 1.6 times a single optimisation; as its own call the user sees it as its own cost, and repeated calls with different $\lambda$ do not repeat stage 1.

### 4.3 Noise at the boundary

A graph judged feasible on the GA subsample may be marginally infeasible on the local subsample or the full sample. The infeasible branch of the score pulls COBYLA back (check 11), best-first pruning only ever accepts full-sample-feasible candidates, and the fallback returns the reference if nothing feasible beats it. The paired comparison keeps the margin small (3.2), but adversarial check A2 asks the reviewer to exercise it deliberately.

## 5. Specification for implementation

### 5.1 Shared helpers

```r
# R/objective_function.R

#' Lexicographic score: feasibility, then fewer edges, then higher gain
#' @noRd
.lexico <- function(u, n_edges, threshold, D, n_free) {
    if (u >= threshold) u - D * n_edges else u - D * (n_free + 1)
}

#' Exact range of the gain function over all rejection patterns
#' @noRd
.trial_success_range <- function(trial_success) {
    m <- trial_success$m
    patterns <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), m)))
    vals <- vapply(
        seq_len(nrow(patterns)),
        function(i) trial_success$func(patterns[i, , drop = FALSE]),
        numeric(1)
    )
    c(min = min(vals), max = max(vals))
}
```

### 5.2 Objective

Complete replacement for the closure factory. The decoding and penalty code is unchanged in content; the additions are the three optional arguments, the constants computed at construction, and the last lines.

```r
create_obj_func <- function(
    m,
    power_criterion,
    hyp_constraint,
    trans_constraint,
    pvals,
    alpha = 0.025,
    num_threads = 1L,
    gain_tolerance = NULL,
    ref_graph = NULL,
    u_range = NULL
) {
    force(power_criterion)
    force(hyp_constraint)
    force(trans_constraint)
    force(alpha)
    force(pvals)
    force(num_threads)

    use_parallel <- num_threads >= 2L
    shortcut <- function(w, G) {
        if (use_parallel) {
            graph_shortcut_parallel(pvals = pvals, alpha = alpha, w = w, G = G,
                num_threads = num_threads, grain_size = -1L)
        } else {
            graph_shortcut(pvals = pvals, alpha = alpha, w = w, G = G)
        }
    }

    enabled <- !is.null(gain_tolerance)
    penalty_base <- 0
    if (enabled) {
        free_mask <- is.na(trans_constraint)
        n_free <- sum(free_mask)
        u_ref <- power_criterion(shortcut(ref_graph$hyp_weight, ref_graph$trans_matrix))
        threshold <- (1 - gain_tolerance) * u_ref
        D <- (u_range[["max"]] - threshold) + 1
        penalty_base <- u_range[["min"]] - D * (n_free + 2)
    }

    function(x) {
        theta <- split_theta(x, hyp_constraint)
        hyp_weight <- recover_full_weights(theta$w_pars, hyp_constraint)
        trans_matrix <- recover_full_trans_matrix(theta$g_pars, trans_constraint)

        if (anyNA(hyp_weight) || anyNA(trans_matrix)) {
            return(penalty_base - 1e6)
        }
        if (any(hyp_weight < 0)) {
            return(penalty_base + sum(hyp_weight[hyp_weight < 0]))
        }
        if (any(trans_matrix < 0)) {
            return(penalty_base + sum(trans_matrix[trans_matrix < 0]))
        }
        if (any(hyp_weight > 1)) {
            return(penalty_base - sum(hyp_weight[hyp_weight > 1]))
        }
        if (any(trans_matrix > 1)) {
            return(penalty_base - sum(trans_matrix[trans_matrix > 1]))
        }

        hyp_weight[hyp_weight < 1e-4] <- 0
        trans_matrix[trans_matrix < 1e-5] <- 0

        u <- power_criterion(shortcut(hyp_weight, trans_matrix))
        if (!enabled) {
            return(u)
        }
        .lexico(u, sum(trans_matrix[free_mask] != 0), threshold, D, n_free)
    }
}
```

With the defaults every returned value is arithmetically the value returned today (`0 + x` and `0 - 1e6` are exact), and no extra random numbers are drawn. The implementer may keep the existing `if (use_parallel)` block inline instead of the `shortcut()` helper; the point is that the disabled path evaluates the same expressions.

### 5.3 Encoding guard and seeds

```r
# R/optimisation_start.R

#' Encode a full graph for the optimisers, guarding derived entries
#'
#' A derived entry (the last free entry of a row, or the last free weight)
#' is one minus the rest. When it should be zero, the rest is rescaled so the
#' derived entry decodes to a small positive value below the zeroing
#' threshold (5e-6 for edges, 5e-5 for weights) instead of an exact zero
#' that floating point could turn into a negative and a penalty.
#' @noRd
.encode_graph <- function(graph_constraint, hyp_weight, trans_matrix) {
    hc <- graph_constraint$hyp_constraint
    tc <- graph_constraint$trans_constraint
    x <- as.numeric(create_start_params(
        graph_constraint, w0 = hyp_weight, G0 = trans_matrix,
        sum_to_one_constraint = FALSE
    ))
    n_w <- max(sum(is.na(hc)) - 1L, 0L)
    if (n_w > 0L) {
        last_w <- max(which(is.na(hc)))
        if (hyp_weight[last_w] < 1e-4) {
            target <- 1 - 5e-5 - sum(hc, na.rm = TRUE)
            s <- sum(x[seq_len(n_w)])
            if (s > 0) x[seq_len(n_w)] <- x[seq_len(n_w)] * (target / s)
        }
    }
    k <- rowSums(is.na(tc))
    rows <- rep.int(seq_along(k), pmax(k - 1L, 0L))
    for (i in unique(rows)) {
        idx <- n_w + which(rows == i)
        last_g <- max(which(is.na(tc[i, ])))
        if (trans_matrix[i, last_g] < 1e-5) {
            target <- 1 - 5e-6 - sum(tc[i, ], na.rm = TRUE)
            s <- sum(x[idx])
            if (s > 0) x[idx] <- x[idx] * (target / s)
        }
    }
    x
}

#' Seed matrix for the simplification GA, most valuable rows first,
#' truncated to pop_size (GA::ga() errors on oversized suggestions).
#' @noRd
.build_simplify_seeds <- function(
    graph_constraint,
    ref_graph,
    pop_size,
    start_graph = NULL,
    population = NULL
) {
    tc <- graph_constraint$trans_constraint
    fixed_edge <- !is.na(tc)
    G <- ref_graph$trans_matrix
    w <- ref_graph$hyp_weight

    seeds <- list(.encode_graph(graph_constraint, w, G))
    cand <- which(G != 0 & !fixed_edge, arr.ind = TRUE)
    for (k in seq_len(nrow(cand))) {
        i <- cand[k, 1L]
        j <- cand[k, 2L]
        Gc <- G
        Gc[i, ] <- .redistribute_mass(G[i, ], drop_idx = j, fixed_idx = which(fixed_edge[i, ]))
        if (sum(Gc != 0) >= sum(G != 0)) next
        if (abs(sum(Gc[i, ]) - 1) > sqrt(.Machine$double.eps)) next
        seeds[[length(seeds) + 1L]] <- .encode_graph(graph_constraint, w, Gc)
    }
    mat <- do.call(rbind, seeds)
    mat <- rbind(mat, .build_start_matrix(graph_constraint, start_graph))
    if (!is.null(population) && ncol(population) == ncol(mat)) {
        mat <- rbind(mat, population)
    }
    mat <- unique(mat, MARGIN = 1)
    storage.mode(mat) <- "double"
    mat[seq_len(min(nrow(mat), pop_size)), , drop = FALSE]
}
```

### 5.4 Zeroing move in the mutation closure

```r
# R/mutation_helpers.R

#' Row index of each free transition parameter in the encoding
#' (NA for hypothesis-weight parameters). Mirrors the parameter order
#' used by recover_full_trans_matrix().
#' @noRd
.g_param_rows <- function(hyp_constraint, trans_constraint) {
    n_w <- max(sum(is.na(hyp_constraint)) - 1L, 0L)
    k <- rowSums(is.na(trans_constraint))
    rows <- rep.int(seq_along(k), pmax(k - 1L, 0L))
    c(rep(NA_integer_, n_w), rows)
}

.make_cauchy_mutation_multi <- function(
    p_param_mutate = 0.1,
    scale = 1.0,
    p_zero = 0,
    param_rows = NULL
) {
    force(p_param_mutate)
    force(scale)
    force(p_zero)
    force(param_rows)

    cauchy_mutation <- function(object, parent) {
        parent_vec <- as.numeric(object@population[parent, ])
        d <- length(parent_vec)
        lower <- object@lower
        upper <- object@upper

        mutate_mask <- stats::runif(d) < p_param_mutate

        if (!any(mutate_mask)) {
            mutate_mask[sample.int(d, 1L)] <- TRUE
        }

        for (j in which(mutate_mask)) {
            parent_vec[j] <- .cauchy_perturb(
                parent_vec[j],
                lower[j],
                upper[j],
                scale
            )
        }

        parent_vec
    }

    if (p_zero <= 0) {
        return(cauchy_mutation)
    }

    g_idx <- which(!is.na(param_rows))
    rows <- unique(param_rows[g_idx])

    function(object, parent) {
        if (stats::runif(1) >= p_zero) {
            return(cauchy_mutation(object, parent))
        }
        parent_vec <- as.numeric(object@population[parent, ])
        if (stats::runif(1) < 0.5) {
            candidates <- g_idx[parent_vec[g_idx] >= 1e-5]
            if (length(candidates) == 0L) {
                return(cauchy_mutation(object, parent))
            }
            j <- candidates[sample.int(length(candidates), 1L)]
            parent_vec[j] <- 0
        } else {
            i <- rows[sample.int(length(rows), 1L)]
            idx <- g_idx[param_rows[g_idx] == i]
            s <- sum(parent_vec[idx])
            if (s <= 0) {
                return(cauchy_mutation(object, parent))
            }
            parent_vec[idx] <- pmin(parent_vec[idx] * ((1 - 5e-6) / s), 1)
        }
        parent_vec
    }
}
```

The inner `cauchy_mutation` is the present closure verbatim and is what the factory returns when `p_zero` is zero. The row-rescale branch assumes no fixed entries in the row other than the diagonal; for rows with pinned non-zero entries the target is $1-5\times10^{-6}$ minus the fixed sum, which the implementer obtains from `trans_constraint` as in `.encode_graph()`.

### 5.5 Best-first pruning with a threshold

```r
# R/post_optim_processing.R: used by prune_graph() when a threshold is given.
# The existing fixed-order prune_edges() is kept unchanged otherwise.

.prune_edges_best_first <- function(
    pvals,
    hyp_weight,
    trans_matrix,
    trial_success,
    fixed_edge,
    threshold,
    alpha = 0.025,
    tolerance = sqrt(.Machine$double.eps)
) {
    u_of <- function(G) {
        calc_power_pvals(
            pvals,
            hyp_weight = hyp_weight,
            trans_matrix = G,
            alpha = alpha,
            custom_power = trial_success
        )$custom_power
    }
    G_best <- trans_matrix
    u_best <- u_of(G_best)
    loss <- 0
    n_removed <- 0L

    repeat {
        n_edges <- sum(G_best != 0)
        cand <- which(G_best != 0 & !fixed_edge, arr.ind = TRUE)
        best <- NULL
        for (k in seq_len(nrow(cand))) {
            i <- cand[k, 1L]
            j <- cand[k, 2L]
            G_try <- G_best
            G_try[i, ] <- .redistribute_mass(
                G_best[i, ],
                drop_idx = j,
                fixed_idx = which(fixed_edge[i, ]),
                tolerance = tolerance
            )
            if (sum(G_try != 0) >= n_edges) next
            if (abs(sum(G_try[i, ]) - 1) > tolerance) next
            u_try <- u_of(G_try)
            if (u_try >= threshold && (is.null(best) || u_try > best$u)) {
                best <- list(G = G_try, u = u_try)
            }
        }
        if (is.null(best)) break
        loss <- loss + (u_best - best$u)
        u_best <- best$u
        G_best <- best$G
        n_removed <- n_removed + 1L
    }

    list(
        hyp_weight = hyp_weight,
        trans_matrix = G_best,
        power_best = u_best,
        prune_loss = loss,
        n_removed = n_removed
    )
}
```

`prune_graph()` gains `threshold = NULL`, forwards it, dispatches to this function when it is non-`NULL`, and its return list gains `prune_loss` (zero otherwise). At most $E_{\max}$ iterations, each at most $E$ shortcut evaluations on the full sample.

### 5.6 `graph_simplify()`

```r
# R/graph_simplify.R

graph_simplify <- function(
    graph_optimal,
    pvals,
    ...,
    gain_tolerance = 1e-3,
    alpha = NULL,
    global_search = TRUE,
    num_threads = 1L,
    control = NULL,
    verbose = multigrain_verbosity()
) {
    check_graph_optimal(graph_optimal)
    check_double_matrix(pvals)
    rlang::check_dots_empty()
    rlang::check_number_decimal(gain_tolerance, min = 0, max = 1)
    rlang::check_number_decimal(alpha, min = 0, max = 1, allow_null = TRUE)
    check_logical(global_search, allow_na = FALSE)
    rlang::check_number_whole(num_threads, min = 1)
    check_control(control, allow_null = TRUE)
    verbose <- .resolve_verbose(verbose)   # same coercion as graph_optimise()

    gc <- graph_optimal$constraints
    ts <- graph_optimal$trial_success
    if (is.null(gc) || is.null(ts)) {
        cli::cli_abort("{.arg graph_optimal} must carry its constraint and trial-success objects.")
    }
    if (ts$m != ncol(pvals)) {
        cli::cli_abort("`pvals` must have one column per hypothesis of {.arg graph_optimal}.")
    }
    alpha <- alpha %||% graph_optimal$alpha
    if (is.null(alpha)) {
        cli::cli_abort("{.arg graph_optimal} does not store {.arg alpha}; supply it explicitly.")
    }

    ref <- list(hyp_weight = unname(graph_optimal$hyp_weight),
                trans_matrix = unname(graph_optimal$trans_matrix))
    ref_power <- calc_power_pvals(pvals, hyp_weight = ref$hyp_weight,
        trans_matrix = ref$trans_matrix, alpha = alpha,
        custom_power = list(trial_success = ts))
    u_ref <- ref_power$trial_success
    if (u_ref <= 0) {
        cli::cli_abort("The reference graph has non-positive trial success on {.arg pvals}; nothing to trade.")
    }
    u_range <- .trial_success_range(ts)
    free_mask <- is.na(gc$trans_constraint)
    n_removable <- sum(pmax(rowSums(free_mask) - 1L, 0L))

    if (is.null(control)) {
        control <- graph_optimal_get_control(graph_optimal) %||% multigrain_control()
        control$global_opt$run <- max(1L, control$global_opt$run %/% 2L)
    }
    control <- control_prepare(control, pvals = pvals, verbose = verbose)

    if (n_removable == 0L) {
        cli::cli_warn("No removable edges under this constraint; returning the reference graph.")
        return(.simplify_result(graph_optimal, ref, ref, ref_power, ref_power, u_ref,
            gain_tolerance, prune_loss = 0, source = "reference", control, alpha,
            global_search, NULL, NULL))
    }

    obj_args <- list(gain_tolerance = gain_tolerance, ref_graph = ref, u_range = u_range)
    seeds <- .build_simplify_seeds(gc, ref, pop_size = control$global_opt$popSize,
        population = if (!is.null(graph_optimal$global_output)) graph_optimal$global_output@population)

    ga_result <- NULL
    x0 <- .encode_graph(gc, ref$hyp_weight, ref$trans_matrix)
    if (global_search) {
        ga_result <- .graph_optimise_ga(pvals, gc, ts, nsim = control$nsim_global,
            global_opts = control$global_opt, alpha = alpha, num_threads = num_threads,
            verbose = verbose, suggestions = seeds, p_zero = 0.2, objective_args = obj_args)
        x0 <- pmin(pmax(ga_result$ga_output@solution[1, ], 0), 1)
    }
    local_result <- .graph_optimise_local(pvals, gc, ts, nsim = control$nsim_local,
        local_opts = control$local_opt, alpha = alpha, num_threads = num_threads,
        x0 = x0, verbose = verbose, objective_args = obj_args)

    best <- choose_graph(ga_result, local_result)   # compares *_objective
    threshold_full <- (1 - gain_tolerance) * u_ref
    pruned <- prune_graph(pvals, best$hyp_weight, best$trans_matrix, ts, gc,
        alpha = alpha, gamma = 1, threshold = threshold_full, verbose = verbose)

    # Fallback: the reference is feasible on the full sample by construction.
    score <- function(w, G) {
        u <- calc_power_pvals(pvals, hyp_weight = w, trans_matrix = G, alpha = alpha,
            custom_power = ts)$custom_power
        n_free <- sum(free_mask)
        D <- (u_range[["max"]] - threshold_full) + 1
        .lexico(u, sum(G[free_mask] != 0), threshold_full, D, n_free)
    }
    use_ref <- score(ref$hyp_weight, ref$trans_matrix) >=
        score(pruned$hyp_weight, pruned$trans_matrix)
    final <- if (use_ref) ref else pruned
    source <- if (use_ref) "reference" else best$source

    final_power <- calc_power_pvals(pvals, hyp_weight = final$hyp_weight,
        trans_matrix = final$trans_matrix, alpha = alpha,
        custom_power = list(trial_success = ts))
    .simplify_result(graph_optimal, ref, final, ref_power, final_power, u_ref,
        gain_tolerance, prune_loss = pruned$prune_loss, source = source,
        control, alpha, global_search, ga_result, local_result)
}
```

`.simplify_result()` applies the hypothesis names, assembles `sparsity` (3.9) and calls `graph_optimal()` with `solution = list(opt_source = paste0("simplify:", source), graph_valid = ...)`. `.graph_optimise_ga()` and `.graph_optimise_local()` gain the optional arguments `suggestions`, `p_zero` and `objective_args`, with defaults reproducing today's behaviour exactly, and return `ga_objective` / `local_objective` alongside the trial-success values.

## 6. Alternatives considered

1. **A single-pass per-edge price inside `graph_optimise()`.** Derive $c=\lambda U_{\max}/P_G$ with $P_G$ the number of removable edges and maximise $U-cE$ in one pass. It certifies a total cap of $\lambda U_{\max}$ without a reference and needs no second optimisation. Rejected because the price is calibrated for the worst case (every removable edge goes) and is therefore tiny: at $\lambda=10^{-3}$ and $m=6$ about $4\times10^{-5}$ per edge, which removes only edges that flip a few dozen trials per million (check 4: an $\epsilon$-edge costing $1.5\times10^{-3}$ would stay), leaves most of the budget unspent, and cannot report the cost against the dense optimum because that optimum is never computed. The two-stage design spends the budget the user stated and reports exactly what it cost.
2. **`gain_tolerance` as an argument of `graph_optimise()`** running the second stage internally. Hides a second optimisation behind an argument, forces a `NULL`-means-off default, and requires a fixture-based identity gate to show that the off path is unchanged. A separate function makes the stage explicit and leaves `graph_optimise()` alone.
3. **An uncalibrated per-edge price.** No total cap; the price's meaning depends on the scale of $\psi$. Rejected.
4. **A multi-objective GA** over $(U,-E)$ with a choice from the Pareto front afterwards. Needs a new dependency and rewrites the global search. Rejected.
5. **Stage 2 without a GA** (COBYLA and best-first pruning only). Cannot find removals that only become affordable after the weights are re-tuned, and cannot move to a different topology. Available to the user through `global_search = FALSE`; not the default.
6. **Entropy or $\ell_q$ sparsity measures; counting pinned entries.** Rejected for the reasons in 3.4.
7. **A threshold on the full sample only.** Would make the GA's and COBYLA's subsample rankings inconsistent with the cap. Rejected in favour of a threshold per sample.

## 7. Implementation plan

Each step ends with a gate that must pass before the next starts. Test files are run one at a time with `testthat::test_file()`.

1. **Shared helpers and objective** (`R/objective_function.R`). Add `.lexico()`, `.trial_success_range()` and the three optional arguments of `create_obj_func()` as in 5.1 and 5.2. Gate: `test-objective_function.R` passes unchanged; new tests: for 200 random encodings at $m\in\{3,4\}$ the closure with default arguments returns values identical (`expect_identical`) to a copy of the pre-change closure kept in the test file; with a threshold it returns `.lexico(U, E, T, D, n_free)` with $E$ computed independently from `param_to_solution(process = TRUE)`; the encoded reference scores exactly $U_{\text{ref}}-D\,E_{\text{ref}}$; every invalid encoding scores below every valid one; `.trial_success_range()` returns $(0,1)$, $(0,4)$ and $(0,1)$ for the three functions in check 6.
2. **Mutation move** (`R/mutation_helpers.R`) as in 5.4, plus `.g_param_rows()`, with the comment block carrying the explanation from 3.6. Gate: `test-mutation_helpers.R` passes; new tests: with `p_zero = 0` the factory returns a closure whose outputs under `set.seed()` are identical to the present closure's; with `p_zero = 1` every call either sets a free transition parameter to exactly 0 or scales a row to $1-5\times10^{-6}$ within `1e-12`; the row map matches `recover_full_trans_matrix()`'s parameter order for a constrained example; the comment block describes both variants.
3. **Encoding guard and seeds** (`R/optimisation_start.R`) as in 5.3. Gate: `test-optimisation_start.R` passes; new tests: encoding then decoding the reference of check 9 gives no negative entry, derived entries of exactly $5\times10^{-6}$ where the reference has zeros, and the same edge count; the seed matrix has the reference first, at most `pop_size` rows, and no duplicates; a population of the wrong width is ignored.
4. **Internal optimisers and choice** (`R/optimisation.R`, `R/choose_graph.R`). Add `suggestions`, `p_zero`, `objective_args` to `.graph_optimise_ga()` and `objective_args` to `.graph_optimise_local()`, both defaulting to today's behaviour; add `ga_objective` / `local_objective`; `choose_graph()` compares them. Gate: `test-optimisation.R` and `test-choose_graph.R` pass with existing snapshots untouched.
5. **Pruning** (`R/post_optim_processing.R`) as in 5.5; `prune_graph()` gains `threshold` and returns `prune_loss`. Gate: `test-post_optim_processing.R` passes unchanged; new tests: on the 4-hypothesis fixture in that file, a threshold equal to the current gain reproduces the current 7-edge result; a threshold never admits a candidate that raises the edge count or lowers $U$ below it; `prune_loss` equals the difference in $U$ before and after; no row is left with a sum different from one.
6. **Object and methods** (`R/graph_optimal.R`, `R/optimisation.R`). Constructor gains `alpha` and `sparsity`; `graph_optimise()` stores `alpha`; print and summary lines. Gate: `test-graph_optimal.R` and `test-optimisation.R` pass with existing snapshots untouched; new test that `names()` of a `graph_optimise()` result equals the previous names plus `alpha` and `sparsity`, with `sparsity` `NULL`; new snapshot of print and summary with a populated `sparsity`.
7. **`graph_simplify()`** (`R/graph_simplify.R`, `NAMESPACE`, `_pkgdown.yml` reference index, roxygen, article). Documentation, all in the roxygen and rendered to `man/graph_simplify.Rd`: `@param` entries stating that the cap is a fraction of the reference gain on the supplied `pvals`, that `alpha` and `control` default to the stored values with `run` halved, and that `pvals` may differ from stage 1; a `@details` section "Edge removal during the global search" carrying the explanation from 3.6 in user-facing terms; a paragraph on what to expect (4.1 to 4.3), including that with `global_search = FALSE` edges are removed only by pruning; the `@returns` entry listing every `sparsity` field from 3.9. The get-started article gains a section calling `graph_simplify()` on the example graph, plotting reference and result side by side, and repeating the edge-removal explanation in one paragraph. Gate: `devtools::document()` runs clean; `man/graph_simplify.Rd` contains the `@details` heading and every `sparsity` field; the article renders; end-to-end test at $m=4$, $n_{\text{sim}}=10^4$ with the fixture of `test-post_optim_processing.R`: for $\lambda\in\{0,10^{-3},10^{-2},1\}$ the returned graph is valid, `gain >= (1 - lambda) * gain_reference`, `n_edges <= n_edges_reference`, `$power$trial_success` is identical to a fresh `calc_power_pvals()` on the returned graph, and at $\lambda=1$ `n_edges_free` equals $m$.
8. **Regenerate `data/graph_optimal_example.rda`** with `data-raw/graph_optimal_example.R` (the object gained two elements), update `NEWS.md` (section 10), run `R CMD check`. Gate: check clean locally.

## 8. Test plan

Beyond the per-step gates above:

**`graph_optimise()` unchanged.** Every existing snapshot in `test-optimisation.R` and `test-graph_optimal.R` passes without update. A new test asserts that a `graph_optimise()` result differs from the previous object shape only by the `alpha` and `sparsity` elements. Because no code on the `graph_optimise()` path changes except the constructor call, no fixture is needed; the reviewer's check A1 covers random-number consumption.

**Cap invariant.** For seeds 1 to 5 and $\lambda\in\{0,10^{-3},5\times10^{-3},10^{-2},1\}$: `sparsity$gain >= (1 - lambda) * sparsity$gain_reference` on the supplied `pvals`, exactly, with no tolerance.

**Fallback.** A reference already at its floor (one free edge per row) and $\lambda=0$ returns the reference with `source = "reference"` and zero loss.

**Seeds.** `.build_simplify_seeds()` returns at most `popSize` rows with the encoded reference first; with the example object's population appended, `GA::ga()` runs without error.

**Inputs.** An object without `alpha` and no explicit `alpha` errors; an explicit `alpha` is used; `pvals` with the wrong number of columns errors; `pvals` with fewer rows than the stored `nsim_*` triggers the existing `control_prepare()` warnings and runs; `num_threads = 2` gives the same result as `num_threads = 1` for the same seed (check 4 verified identical rejections from the parallel shortcut).

**Degenerate cases.** `graph_constraint_free(2)` warns and returns the reference; a gain function that is never positive errors before any search.

## 9. Adversarial checks for review

- **A1.** Run `graph_optimise()` before and after the change for seeds 1 to 5, $m\in\{2,3,4\}$, constrained and free, `global_search` on and off, `num_threads` 1 and 2. Any non-identical element other than `alpha` and `sparsity`, or any difference in `.Random.seed` afterwards, fails the unchanged-behaviour requirement.
- **A2.** Boundary feasibility across samples: choose $\lambda$ so that a one-edge-removal neighbour of the reference sits within one paired standard error of the threshold on the GA subsample. Confirm that whatever the GA and COBYLA return, the final graph satisfies the cap on the full sample, and that `source` reports `"reference"` when nothing feasible beat it.
- **A3.** Seeds larger than `popSize`: a reference with many edges and a small `popSize`. `GA::ga()` must not error; the reference must still be in the population.
- **A4.** Encoding guard on constrained graphs with pinned non-zero entries in rows that also have free entries; on a reference whose derived weight is zero; on a row whose parameters sum to more than one before the guard.
- **A5.** Gain functions with negative values and with $U_{\text{ref}}\le0$: the error path, not a silent result.
- **A6.** $\lambda=0$ and $\lambda=1$: the invariants in step 7; at $\lambda=1$ the result has exactly one free edge per row and `source` is not `"reference"`.
- **A7.** `pvals` different from stage 1 (a fresh sample of the same design, and a sample of a different design): `gain_reference` is recomputed on the supplied matrix and the cap holds there.
- **A8.** Stored population of the wrong width (object optimised under a different constraint than the one passed): seeds ignore it, no error.
- **A9.** Start COBYLA from a graph infeasible at the chosen $\lambda$ (as in check 11) and confirm the returned graph is feasible or the fallback fires; never an infeasible return.
- **A10.** Stored `run` of 1: halving floors at 1, GA still runs.
- **A11.** Rows with pinned non-zero entries under the row-rescale mutation move: the rescale target must subtract the fixed sum; check that decoded entries stay non-negative.
- **A12.** Best-first pruning termination and cost: at most $E_{\max}$ iterations; measure shortcut evaluations at $m=4$ and extrapolate to $m=8$, $n_{\text{sim}}=10^6$.
- **A13.** $\epsilon$-edges: after `param_to_solution(process = TRUE)` snaps an entry to $0.001$, confirm pruning can remove it when feasible and that it is counted in `n_edges_free` when it stays.
- **A14.** The intermittent Nelder-Mead step inside `GA::ga()` receives the lexicographic fitness and does not re-create edges in feasible individuals: count zeros before and after `optim` steps on a logged run.
- **A15.** Two consecutive calls, `graph_simplify(graph_simplify(res, pvals, 1e-3), pvals, 1e-3)`: the second call's reference is the first call's result, so the cap compounds; confirm the documentation says so and that nothing breaks.

## 10. Proposed `NEWS.md` wording

Under `# multigrain (development version)`:

```
## New functionality

* New `graph_simplify()` takes an optimised graph and searches for one with
  fewer edges whose trial-success value, on the supplied p-values, is at
  least `1 - gain_tolerance` times that of the input graph. The search reuses
  the global and local optimisers with a score that ranks feasibility first,
  fewer edges second and trial success third, warm-starts from the input
  graph and the stored population, and finishes with best-first pruning; the
  input graph is returned if nothing feasible improves on it. The result is a
  `multigrain_graph_optimal` with a `sparsity` element reporting both graphs'
  trial success, the edge counts and the exact loss, shown by `print()` and
  `summary()`.
* `graph_optimise()` results now store `alpha` and carry a `sparsity`
  element (`NULL`); nothing else about `graph_optimise()` changes.
```

No bug fixes are part of this change. Two pre-existing issues found while reading the code are recorded as open items O2 and O3 and should be fixed separately.

## 11. Open items

- **O1. $p_0$ and the split between parameter zeroing and row rescaling.** Proposed $0.2$ and one half. Evidence to close: on $m=4$ examples, the fraction of stage-2 generations in which the elite's edge count decreases, for $p_0\in\{0.1,0.2,0.4\}$.
- **O2. `global_opt_power` / `local_opt_power`.** The 0.2.0 news entry says the returned object stores the pre-pruning power of the global and local solutions. The constructor in `R/graph_optimal.R` has no such elements; `.graph_optimise_ga()` and `.graph_optimise_local()` compute them as `ga_subset_power` and `local_subset_power`, and `graph_optimise()` discards them. Either reinstate them or correct the news entry. Separate change.
- **O3. Pre-existing abort in `prune_edges()`** when a row has no free recipients (check 1). Unreachable through `graph_constraint()`; fix separately if `prune_graph()` is ever exported.
- **O4. Halving `run` for stage 2.** Proposed default. Evidence to close: on realistic problems of five to eight hypotheses, the generation at which stage 2 last improved, compared with its `run` value.
- **O5. Whether `graph_optimise()` should offer a pass-through** (for example `simplify = 1e-3`) that calls `graph_simplify()` on its own result. Convenience only; decide after the function has been used.
- **O6. Warning when `pvals` differs from stage 1.** The object could store a hash of the p-value matrix and `graph_simplify()` could note when the supplied matrix differs. Decide after use.

## Appendix: verification runs

Environment: R 4.6.1 (Windows), Rtools 4.5 on the path, package loaded with `pkgload::load_all()` from the working tree at commit `8b19c63`. All runs used $m\le4$ and $n_{\text{sim}}\le10^4$. The script is reproduced first, then the output.

```r
## Verification checks for the sparsity design record.
## Constraints: m <= 4, nsim <= 1e4, no test_package(), read-only w.r.t. repo.
.libPaths(c("C:/Users/advsp/AppData/Local/R/win-library/4.6", .libPaths()))
suppressPackageStartupMessages({
    pkgload::load_all("C:/Users/advsp/Documents/repos/multigrain", quiet = TRUE)
})
options(multigrain_verbosity = "silent")
sec <- function(x) cat("\n==== ", x, " ====\n", sep = "")

## ---------------------------------------------------------------- 1
sec("1. prune_edges: dropping a row's only free edge -> row sum < 1 -> abort?")
tc_single <- rbind(c(0, 0.6, NA), c(NA, 0, NA), c(NA, NA, 0))
res_gc <- tryCatch(graph_constraint(trans_constraint = tc_single), error = function(e) conditionMessage(e))
cat("graph_constraint() with a single-NA row ->", res_gc, "\n")
tc_sum1 <- rbind(c(0, 1, NA, NA), c(NA, 0, NA, NA), c(NA, NA, 0, NA), c(NA, NA, NA, 0))
res_gc2 <- tryCatch(graph_constraint(trans_constraint = tc_sum1), error = function(e) conditionMessage(e))
cat("graph_constraint() with fixed row sum 1 + NAs ->", res_gc2, "\n")
fixed_edge3 <- !is.na(tc_single) # bypass validation to exercise prune_edges directly
set.seed(1)
pv3 <- simulate_pvalues(c(0.9, 0.8, 0.7), corr_matrix = diag(3), nsim = 5e3)
ts3 <- trial_success(r1 + r2 + r3, verbose = "silent")
w3 <- c(0.5, 0.3, 0.2)
G3 <- rbind(c(0, 0.6, 0.4), c(0.5, 0, 0.5), c(0.5, 0.5, 0))
cat("valid input graph:", is_graph_valid(w3, G3), "\n")
rm3 <- .redistribute_mass(G3[1, ], drop_idx = 3L, fixed_idx = c(1L, 2L))
cat("row after .redistribute_mass (no recipients):", rm3, " sum =", sum(rm3), "\n")
res1 <- tryCatch(
    prune_edges(
        pv3, hyp_weight = w3, trans_matrix = G3, trial_success = ts3,
        fixed_edge = fixed_edge3, power_best = 0, gamma = 1
    ),
    error = function(e) paste("ERROR:", conditionMessage(e))
)
print(res1)

## ---------------------------------------------------------------- 2
sec("2. derived entry can be a tiny negative; objective returns sum(negatives)")
gc4 <- graph_constraint_free(4)
set.seed(2)
pv4 <- simulate_pvalues(c(0.93, 0.91, 0.90, 0.85),
    corr_matrix = matrix(0.2, 4, 4) + diag(0.8, 4), nsim = 1e4)
ts4 <- trial_success(r1 && r2 && r3 && r4, verbose = "silent")
ts_custom <- trial_success(0.25 * (2 * (r1 && r2) + r1 * r3 + r2 * r4), verbose = "silent")
obj4 <- create_obj_func(4, ts4$func, gc4$hyp_constraint, gc4$trans_constraint,
    pvals = pv4)
x0 <- as.numeric(create_start_params(gc4))
cat("x0 length:", length(x0), " f(x0) =", obj4(x0), "\n")
show_derived <- function(a, b, label) {
    x <- x0
    x[4:5] <- c(a, b) # row-1 params (after the 3 weight params)
    G <- recover_full_trans_matrix(split_theta(x, gc4$hyp_constraint)$g_pars,
        gc4$trans_constraint)
    cat(sprintf("%s: params (%.17g, %.17g) -> derived G[1,4] = %.3g ; f(x) = %.3g\n",
        label, a, b, G[1, 4], obj4(x)))
}
show_derived(0.6, 0.6, "row sum 1.2 (ordinary overshoot)")
a <- 0.3; b <- (1 - a) * (1 + .Machine$double.eps)
show_derived(a, b, "row sum 1 + 1 ulp")
show_derived(0.3, 0.7 - 5e-6, "row scaled to 1 - 5e-6 (zeroing move)")

## ---------------------------------------------------------------- 3
sec("3. how often does Cauchy mutation / GA crossover create a zero?")
set.seed(4)
n <- 2e5
below <- mean(replicate(n, .cauchy_perturb(0.3, 0, 1)) < 1e-5)
below_small <- mean(replicate(n, .cauchy_perturb(0.001, 0, 1)) < 1e-5)
cat(sprintf("P(perturb(0.3) < 1e-5) = %.2e ; P(perturb(0.001) < 1e-5) = %.2e (n=%d)\n",
    below, below_small, n))
ctl <- GA::gaControl("real-valued")
cat("GA default real-valued crossover:", ctl$crossover, "\n")
print(body(GA:::gareal_laCrossover_R))

## ---------------------------------------------------------------- 4
sec("4. paired-difference SE vs marginal SE when removing one small edge (m=4, nsim=1e4)")
w <- c(0.5, 0.3, 0.2, 0)
G <- rbind(
    c(0, 0.6, 0.39, 0.01),
    c(0.5, 0, 0.5, 0),
    c(0.3, 0.3, 0, 0.4),
    c(0.5, 0.5, 0, 0)
)
stopifnot(is_graph_valid(w, G))
ts_avg <- trial_success(r1 + r2 + r3 + r4, verbose = "silent")
psi_row <- function(rej) rowSums(rej) # r1+r2+r3+r4 per simulation
rej_a <- graph_shortcut(pv4, 0.025, w, G)
G2 <- G
G2[1, ] <- .redistribute_mass(G[1, ], drop_idx = 4L, fixed_idx = 1L)
rej_b <- graph_shortcut(pv4, 0.025, w, G2)
pa <- psi_row(rej_a); pb <- psi_row(rej_b)
cat(sprintf("U(G) = %.5f, U(G - edge) = %.5f, diff = %.2e\n", mean(pa), mean(pb), mean(pa) - mean(pb)))
cat(sprintf("marginal SE(U) = %.2e ; paired SE(diff) = %.2e ; sims whose outcome changed = %d / %d\n",
    sd(pa) / sqrt(length(pa)), sd(pa - pb) / sqrt(length(pa)), sum(pa != pb), length(pa)))
# same for the 0.39 edge
G3b <- G
G3b[1, ] <- .redistribute_mass(G[1, ], drop_idx = 3L, fixed_idx = 1L)
pc <- psi_row(graph_shortcut(pv4, 0.025, w, G3b))
cat(sprintf("removing the 0.39 edge: diff = %.2e, paired SE = %.2e, changed = %d\n",
    mean(pa) - mean(pc), sd(pa - pc) / sqrt(length(pa)), sum(pa != pc)))
# an epsilon edge (0.001) on the custom (0-1 scale) objective
Ge <- G; Ge[1, ] <- c(0, 0.6, 0.399, 0.001)
Ge0 <- Ge; Ge0[1, ] <- .redistribute_mass(Ge[1, ], drop_idx = 4L, fixed_idx = 1L)
psi_c <- function(rej) vapply(seq_len(nrow(rej)), function(i) ts_custom$func(rej[i, , drop = FALSE]), numeric(1))
pe <- psi_c(graph_shortcut(pv4, 0.025, w, Ge)); pe0 <- psi_c(graph_shortcut(pv4, 0.025, w, Ge0))
cat(sprintf("epsilon edge 0.001, custom psi in [0,1]: U = %.5f, diff = %.2e, marginal SE = %.2e, paired SE = %.2e, changed = %d\n",
    mean(pe), mean(pe) - mean(pe0), sd(pe) / sqrt(length(pe)), sd(pe - pe0) / sqrt(length(pe)), sum(pe != pe0)))
# parallel backend gives the same rejections (edge count unaffected by num_threads)
cat("parallel == serial rejections:", identical(graph_shortcut_parallel(pv4, 0.025, w, G, num_threads = 2L), graph_shortcut(pv4, 0.025, w, G)), "\n")

## ---------------------------------------------------------------- 5
sec("5. exact zeros survive param_to_solution / repair_graph / normalise_sum")
x <- x0
x[4] <- 0            # row-1 first param exactly 0
sol <- param_to_solution(x, gc4, process = TRUE)
cat("param 0 -> G[1,2] identical to 0:", identical(sol$trans_matrix[1, 2], 0), "\n")
rg <- repair_graph(sol$hyp_weight, sol$trans_matrix, gc4)
cat("after repair_graph still identical 0:", identical(rg$trans_matrix[1, 2], 0), "\n")
ns <- normalise_sum(c(0, 0.3, 0.7000001), fixed_idx = 1L)
cat("normalise_sum keeps 0:", identical(ns[1], 0), " sum:", sum(ns), "\n")
cat("row with all-zero free entries after repair_graph:\n")
Gz <- sol$trans_matrix; Gz[2, ] <- 0
print(repair_graph(sol$hyp_weight, Gz, gc4)$trans_matrix[2, ])

## ---------------------------------------------------------------- 6
sec("6. U_max over all 2^m patterns, and P_G from the encoding")
u_max <- function(ts) {
    m <- ts$m
    pats <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), m)))
    vals <- vapply(seq_len(nrow(pats)), function(i) ts$func(pats[i, , drop = FALSE]), numeric(1))
    c(min = min(vals), max = max(vals))
}
print(rbind(conj = u_max(ts4), avg = u_max(ts_avg), custom = u_max(ts_custom)))
p_g <- function(gc) sum(pmax(rowSums(is.na(gc$trans_constraint)) - 1L, 0L))
gc_c <- graph_constraint(
    hyp_constraint = c(NA, NA, 0, 0),
    trans_constraint = rbind(c(0, NA, NA, 0), c(NA, 0, 0, NA), c(0, 1, 0, 0), c(1, 0, 0, 0))
)
for (g in list(free4 = gc4, constrained = gc_c)) {
    n_g_params <- length(create_start_params(g)) - max(sum(is.na(g$hyp_constraint)) - 1L, 0L)
    cat(sprintf("P_G by formula = %d ; free G params in encoding = %d ; m(m-2) = %d\n",
        p_g(g), n_g_params, 4 * 2))
}

## ---------------------------------------------------------------- 7
sec("7. toy best-first prune with acceptance dU <= c, vs current gamma = 1 prune")
ctrl <- multigrain_control() |> control_local(maxeval = 2000, print_level = 0)
ctrl <- control_prepare(ctrl, pvals = pv4)
set.seed(5)
loc <- .graph_optimise_local(pv4, gc4, ts_custom, local_opts = ctrl$local_opt, verbose = "silent")
w_raw <- loc$local_hyp_weight; G_raw <- loc$local_trans_matrix
cat("raw local optimum: edges =", sum(G_raw != 0), "\n"); print(round(G_raw, 4))
U_of <- function(w, G) calc_power_pvals(pv4, hyp_weight = w, trans_matrix = G,
    custom_power = ts_custom)$custom_power
U_raw <- U_of(w_raw, G_raw)
cur <- prune_graph(pv4, w_raw, G_raw, ts_custom, gc4, gamma = 1, verbose = "silent")
cat(sprintf("current prune (gamma=1): edges %d -> %d, U %.5f -> %.5f\n",
    sum(G_raw != 0), sum(cur$trans_matrix != 0), U_raw, U_of(cur$hyp_weight, cur$trans_matrix)))
best_first <- function(w, G, c, gc) {
    fixed_edge <- !is.na(gc$trans_constraint)
    U_best <- U_of(w, G); loss <- 0; removed <- 0L
    repeat {
        cand <- which(G != 0 & !fixed_edge, arr.ind = TRUE)
        if (!nrow(cand)) break
        trial <- lapply(seq_len(nrow(cand)), function(k) {
            i <- cand[k, 1]; j <- cand[k, 2]
            Gc <- G
            Gc[i, ] <- .redistribute_mass(G[i, ], drop_idx = j, fixed_idx = which(fixed_edge[i, ]))
            if (sum(Gc[i, ] != 0) == 0 || abs(sum(Gc[i, ]) - 1) > 1e-8) return(NULL)
            list(G = Gc, U = U_of(w, Gc), E = sum(Gc != 0))
        })
        trial <- Filter(Negate(is.null), trial)
        trial <- Filter(function(t) t$E < sum(G != 0), trial) # must reduce count
        if (!length(trial)) break
        Us <- vapply(trial, `[[`, numeric(1), "U")
        k <- which.max(Us)
        if (U_best - Us[k] > c) {
            cat(sprintf("   stop: cheapest remaining edge costs %.2e (> c = %.2e)\n", U_best - Us[k], c))
            break
        }
        loss <- loss + (U_best - Us[k]); U_best <- Us[k]; G <- trial[[k]]$G; removed <- removed + 1L
    }
    list(G = G, U = U_best, loss = loss, removed = removed)
}
umax <- u_max(ts_custom)[["max"]]
for (lam in c(0, 1e-3, 5e-3, 1e-2)) {
    cc <- lam * umax / p_g(gc4)
    bf <- best_first(cur$hyp_weight, cur$trans_matrix, cc, gc4)
    cat(sprintf("after gamma=1 prune, lambda=%.3f  c=%.2e : edges %d -> %d, removed %d, prune loss %.2e (bound lambda*Umax = %.2e), U = %.5f\n",
        lam, cc, sum(cur$trans_matrix != 0), sum(bf$G != 0), bf$removed, bf$loss, lam * umax, bf$U))
}
cat("\nbest-first from the RAW graph (weights from gamma=1 prune), c = 0 vs current fixed-order prune:\n")
bf0 <- best_first(cur$hyp_weight, G_raw, 0, gc4)
cat(sprintf("  best-first c=0: edges %d -> %d, U = %.5f ; current prune: %d edges, U = %.5f\n",
    sum(G_raw != 0), sum(bf0$G != 0), bf0$U, sum(cur$trans_matrix != 0), U_of(cur$hyp_weight, cur$trans_matrix)))
print(round(bf0$G, 4))
cat("\nDONE\n")
```

Output (warnings emitted by `is_graph_valid()` inside check 1 omitted):

```
==== 1. prune_edges: dropping a row's only free edge -> row sum < 1 -> abort? ====
graph_constraint() with a single-NA row -> At least one incomplete transition matrix row has a single optimisable
(i.e. `NA`) value.
graph_constraint() with fixed row sum 1 + NAs -> At least one incomplete transition matrix row has a sum equal to 1.
valid input graph: TRUE
row after .redistribute_mass (no recipients): 0 0.6 0  sum = 0.6
[1] "ERROR: The supplied `hyp_weight` and `trans_matrix` do not build a valid graph."

==== 2. derived entry can be a tiny negative; objective returns sum(negatives) ====
x0 length: 11  f(x0) = 0.6572
row sum 1.2 (ordinary overshoot): params (0.59999999999999998, 0.59999999999999998) -> derived G[1,4] = -0.2 ; f(x) = -0.2
row sum 1 + 1 ulp: params (0.29999999999999999, 0.70000000000000007) -> derived G[1,4] = 0 ; f(x) = 0.658
row scaled to 1 - 5e-6 (zeroing move): params (0.29999999999999999, 0.69999499999999992) -> derived G[1,4] = 5e-06 ; f(x) = 0.658

==== 3. how often does Cauchy mutation / GA crossover create a zero? ====
P(perturb(0.3) < 1e-5) = 1.00e-05 ; P(perturb(0.001) < 1e-5) = 1.00e-05 (n=200000)
GA default real-valued crossover: gareal_laCrossover
{
    parents <- object@population[parents, , drop = FALSE]
    n <- ncol(parents)
    children <- matrix(as.double(NA), nrow = 2, ncol = n)
    a <- runif(n)
    children[1, ] <- a * parents[1, ] + (1 - a) * parents[2, ]
    children[2, ] <- a * parents[2, ] + (1 - a) * parents[1, ]
    out <- list(children = children, fitness = rep(as.double(NA), 2))
    return(out)
}

==== 4. paired-difference SE vs marginal SE when removing one small edge (m=4, nsim=1e4) ====
U(G) = 3.41110, U(G - edge) = 3.39380, diff = 1.73e-02
marginal SE(U) = 9.72e-03 ; paired SE(diff) = 1.39e-03 ; sims whose outcome changed = 186 / 10000
removing the 0.39 edge: diff = 2.70e-03, paired SE = 2.38e-03, changed = 266
epsilon edge 0.001, custom psi in [0,1]: U = 0.78892, diff = 1.48e-03, marginal SE = 3.44e-03, paired SE = 1.91e-04, changed = 59
parallel == serial rejections: TRUE

==== 5. exact zeros survive param_to_solution / repair_graph / normalise_sum ====
param 0 -> G[1,2] identical to 0: TRUE
after repair_graph still identical 0: TRUE
normalise_sum keeps 0: TRUE  sum: 1
row with all-zero free entries after repair_graph:
[1] 0.3333333 0.0000000 0.3333333 0.3333333

==== 6. U_max over all 2^m patterns, and P_G from the encoding ====
       min max
conj     0   1
avg      0   4
custom   0   1
P_G by formula = 8 ; free G params in encoding = 8 ; m(m-2) = 8
P_G by formula = 2 ; free G params in encoding = 2 ; m(m-2) = 8

==== 7. toy best-first prune with acceptance dU <= c, vs current gamma = 1 prune ====
raw local optimum: edges = 12
       [,1]   [,2]   [,3]   [,4]
[1,] 0.0000 0.6410 0.3351 0.0239
[2,] 0.2834 0.0000 0.3354 0.3812
[3,] 0.5966 0.3361 0.0000 0.0673
[4,] 0.3369 0.5851 0.0780 0.0000
current prune (gamma=1): edges 12 -> 7, U 0.79788 -> 0.80483
   stop: cheapest remaining edge costs 7.50e-05 (> c = 0.00e+00)
after gamma=1 prune, lambda=0.000  c=0.00e+00 : edges 7 -> 7, removed 0, prune loss 0.00e+00 (bound lambda*Umax = 0.00e+00), U = 0.80483
   stop: cheapest remaining edge costs 9.67e-03 (> c = 1.25e-04)
after gamma=1 prune, lambda=0.001  c=1.25e-04 : edges 7 -> 6, removed 1, prune loss 7.50e-05 (bound lambda*Umax = 1.00e-03), U = 0.80475
   stop: cheapest remaining edge costs 9.67e-03 (> c = 6.25e-04)
after gamma=1 prune, lambda=0.005  c=6.25e-04 : edges 7 -> 6, removed 1, prune loss 7.50e-05 (bound lambda*Umax = 5.00e-03), U = 0.80475
   stop: cheapest remaining edge costs 9.67e-03 (> c = 1.25e-03)
after gamma=1 prune, lambda=0.010  c=1.25e-03 : edges 7 -> 6, removed 1, prune loss 7.50e-05 (bound lambda*Umax = 1.00e-02), U = 0.80475

best-first from the RAW graph (weights from gamma=1 prune), c = 0 vs current fixed-order prune:
   stop: cheapest remaining edge costs 7.50e-05 (> c = 0.00e+00)
  best-first c=0: edges 12 -> 7, U = 0.80483 ; current prune: 7 edges, U = 0.80483
       [,1]   [,2]   [,3]   [,4]
[1,] 0.0000 0.6567 0.3433 0.0000
[2,] 0.2834 0.0000 0.3354 0.3812
[3,] 0.0000 1.0000 0.0000 0.0000
[4,] 0.0000 1.0000 0.0000 0.0000

DONE
```

Two remarks on the output. In check 2 the local optimum used to seed `x0` differs from the graph used in `show_derived`, so `f(x0)` and the `f(x)` values are not comparable with each other; what matters is that the rescaled row gives a derived entry of exactly $5\times10^{-6}$ and a finite objective. In check 7 the local search was capped at 2000 evaluations and had not converged, which is why the existing prune raised $U$; the comparison of interest is best-first at zero price against the fixed-order prune, which agree.

### Additional checks 8 to 11

Same environment. These checks build the reference graph the same way as check 7 (a local optimum at 2000 evaluations, then the current prune), so the reference has seven edges and $U_{\text{ref}}=0.804825$ on the $10^4$ trials.

```r
## Additional verification checks for the graph_simplify() design.
## Constraints: m <= 4, nsim <= 1e4, no test_package(), read-only w.r.t. repo.
.libPaths(c("C:/Users/advsp/AppData/Local/R/win-library/4.6", .libPaths()))
suppressPackageStartupMessages({
    pkgload::load_all("C:/Users/advsp/Documents/repos/multigrain", quiet = TRUE)
})
options(multigrain_verbosity = "silent")
sec <- function(x) cat("\n==== ", x, " ====\n", sep = "")

## shared setup (same as checks.R sections 2 and 7)
gc4 <- graph_constraint_free(4)
set.seed(2)
pv4 <- simulate_pvalues(c(0.93, 0.91, 0.90, 0.85),
    corr_matrix = matrix(0.2, 4, 4) + diag(0.8, 4), nsim = 1e4)
ts_custom <- trial_success(0.25 * (2 * (r1 && r2) + r1 * r3 + r2 * r4), verbose = "silent")
U_of <- function(w, G) calc_power_pvals(pv4, hyp_weight = w, trans_matrix = G,
    custom_power = ts_custom)$custom_power
ctrl <- multigrain_control() |> control_local(maxeval = 2000, print_level = 0)
ctrl <- control_prepare(ctrl, pvals = pv4)
set.seed(5)
loc <- .graph_optimise_local(pv4, gc4, ts_custom, local_opts = ctrl$local_opt, verbose = "silent")
cur <- prune_graph(pv4, loc$local_hyp_weight, loc$local_trans_matrix, ts_custom, gc4,
    gamma = 1, verbose = "silent")
w_ref <- cur$hyp_weight; G_ref <- cur$trans_matrix
U_ref <- U_of(w_ref, G_ref)
cat("reference (stage-1 output): edges =", sum(G_ref != 0), " U_ref =", format(U_ref, digits = 6), "\n")
print(round(G_ref, 4)); cat("w_ref:", round(w_ref, 4), "\n")
free_mask <- is.na(gc4$trans_constraint)
n_free <- sum(free_mask)

## ---------------------------------------------------------------- 8
sec("8. stage-1 GA population is stored on the object; GA truncates suggestions")
x <- graph_optimal_example
cat("names(graph_optimal_example):", names(x), "\n")
cat("population dim:", dim(x$global_output@population), " popSize:", x$global_output@popSize, "\n")
nvars_example <- length(create_start_params(x$constraints))
cat("encoding nvars for its constraint:", nvars_example, " (matches ncol:", ncol(x$global_output@population) == nvars_example, ")\n")
set.seed(9)
sugg <- matrix(runif(12 * 3), nrow = 12, ncol = 3)
tiny_ga <- function(s, pop) suppressWarnings(GA::ga(type = "real-valued",
    fitness = function(v) -sum((v - 0.5)^2), lower = rep(0, 3), upper = rep(1, 3),
    popSize = pop, maxiter = 3, run = 3, suggestions = s, monitor = FALSE))
res_over <- tryCatch(tiny_ga(sugg, 10), error = function(e) paste("ERROR:", conditionMessage(e)))
cat("GA::ga with popSize 10 and 12 suggestion rows ->", if (is.character(res_over)) res_over else "ran", "\n")
res_ok <- tiny_ga(sugg[1:10, ], 10)
cat("GA::ga with popSize 10 and 10 suggestion rows -> ran, iter =", res_ok@iter, "\n")
res_fewer <- tiny_ga(sugg[1:3, ], 10)
cat("GA::ga with popSize 10 and 3 suggestion rows -> ran, iter =", res_fewer@iter, "\n")

## ---------------------------------------------------------------- 9
sec("9. encoding a pruned graph: derived entries, the guard, and the objective")
obj_plain <- create_obj_func(4, ts_custom$func, gc4$hyp_constraint, gc4$trans_constraint, pvals = pv4)
x_ref <- as.numeric(create_start_params(gc4, w0 = w_ref, G0 = G_ref, sum_to_one_constraint = FALSE))
decode <- function(xx) {
    th <- split_theta(xx, gc4$hyp_constraint)
    list(w = recover_full_weights(th$w_pars, gc4$hyp_constraint),
         G = recover_full_trans_matrix(th$g_pars, gc4$trans_constraint))
}
count_edges <- function(G) { G[G < 1e-5] <- 0; sum(G[free_mask] != 0) }
d0 <- decode(x_ref)
cat("derived col per row (last free col): ", apply(gc4$trans_constraint, 1, function(r) max(which(is.na(r)))), "\n")
cat("derived entries after plain re-encoding: ", format(c(d0$G[1, 4], d0$G[2, 4], d0$G[3, 4], d0$G[4, 3]), digits = 3), "\n")
cat("derived weight (w[4]):", format(d0$w[4], digits = 3), "\n")
cat("any negative entry:", any(d0$G < 0) || any(d0$w < 0), " -> objective value:", format(obj_plain(x_ref), digits = 6),
    " (U of reference on same sample:", format(U_ref, digits = 6), ")\n")
# the guard: rescale rows whose derived entry is < 1e-5 so params sum to 1 - 5e-6;
# rescale weights if the derived weight is < 1e-4 so they sum to 1 - 5e-5
guard_encode <- function(xx, gc) {
    n_w <- max(sum(is.na(gc$hyp_constraint)) - 1L, 0L)
    k <- rowSums(is.na(gc$trans_constraint))
    rows <- rep.int(seq_along(k), pmax(k - 1L, 0L))
    d <- decode(xx)
    if (n_w > 0) {
        last_free_w <- max(which(is.na(gc$hyp_constraint)))
        if (d$w[last_free_w] < 1e-4) {
            s <- sum(xx[seq_len(n_w)]) + sum(gc$hyp_constraint, na.rm = TRUE)
            if (s > 0) xx[seq_len(n_w)] <- xx[seq_len(n_w)] * ((1 - 5e-5 - sum(gc$hyp_constraint, na.rm = TRUE)) / sum(xx[seq_len(n_w)]))
        }
    }
    for (i in unique(rows)) {
        idx <- n_w + which(rows == i)
        last_free <- max(which(is.na(gc$trans_constraint[i, ])))
        if (d$G[i, last_free] < 1e-5) {
            fixed_sum <- sum(gc$trans_constraint[i, ], na.rm = TRUE)
            s <- sum(xx[idx])
            if (s > 0) xx[idx] <- xx[idx] * ((1 - 5e-6 - fixed_sum) / s)
        }
    }
    xx
}
x_ref_g <- guard_encode(x_ref, gc4)
d1 <- decode(x_ref_g)
cat("derived entries after guard:              ", format(c(d1$G[1, 4], d1$G[2, 4], d1$G[3, 4], d1$G[4, 3]), digits = 3), "\n")
cat("derived weight after guard:", format(d1$w[4], digits = 3), "\n")
cat("objective after guard:", format(obj_plain(x_ref_g), digits = 6),
    "; edges counted:", count_edges(d1$G), " (reference has", sum(G_ref != 0), ")\n")
# an adversarial case: parameters summing to exactly 1 in a row (derived exactly 0)
xx <- x_ref; xx[4:5] <- c(0.3, 0.7)
cat("row params (0.3, 0.7): derived =", format(decode(xx)$G[1, 4], digits = 3), "; objective:", format(obj_plain(xx), digits = 6), "\n")
xx[4:5] <- c(0.1, 0.9)
cat("row params (0.1, 0.9): derived =", format(decode(xx)$G[1, 4], digits = 3), "; objective:", format(obj_plain(xx), digits = 6), "\n")
xx[4:5] <- c(0.7, 0.3)
cat("row params (0.7, 0.3): derived =", format(decode(xx)$G[1, 4], digits = 3), "; objective:", format(obj_plain(xx), digits = 6), "\n")
xx[4:5] <- c(0.35, 0.65)
cat("row params (0.35, 0.65): derived =", format(decode(xx)$G[1, 4], digits = 3), "; objective:", format(obj_plain(xx), digits = 6), "\n")

## ---------------------------------------------------------------- 10
sec("10. threshold prune on the toy reference: budget spent in full")
best_first_T <- function(w, G, threshold, gc) {
    fixed_edge <- !is.na(gc$trans_constraint)
    U_best <- U_of(w, G); loss <- 0; removed <- 0L
    repeat {
        cand <- which(G != 0 & !fixed_edge, arr.ind = TRUE)
        best <- NULL
        for (k in seq_len(nrow(cand))) {
            i <- cand[k, 1]; j <- cand[k, 2]
            Gc <- G
            Gc[i, ] <- .redistribute_mass(G[i, ], drop_idx = j, fixed_idx = which(fixed_edge[i, ]))
            if (sum(Gc != 0) >= sum(G != 0)) next
            if (abs(sum(Gc[i, ]) - 1) > 1e-8) next
            u <- U_of(w, Gc)
            if (u >= threshold && (is.null(best) || u > best$u)) best <- list(G = Gc, u = u)
        }
        if (is.null(best)) break
        loss <- loss + (U_best - best$u); U_best <- best$u; G <- best$G; removed <- removed + 1L
    }
    list(G = G, U = U_best, loss = loss, removed = removed)
}
for (lam in c(1e-3, 1e-2, 2e-2)) {
    Tt <- (1 - lam) * U_ref
    r <- best_first_T(w_ref, G_ref, Tt, gc4)
    cat(sprintf("lambda=%.3f  T=%.5f : edges %d -> %d, removed %d, loss %.2e (budget %.2e, spent %.0f%%), U=%.5f, within cap: %s\n",
        lam, Tt, sum(G_ref != 0), sum(r$G != 0), r$removed, r$loss, lam * U_ref, 100 * r$loss / (lam * U_ref), r$U, r$U >= Tt))
    if (lam == 2e-2) print(round(r$G, 4))
}

## ---------------------------------------------------------------- 11
sec("11. stage-2 COBYLA within a support: lexicographic objective from the encoded reference")
u_range <- {
    pats <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), 4)))
    v <- vapply(seq_len(nrow(pats)), function(i) ts_custom$func(pats[i, , drop = FALSE]), numeric(1))
    c(min = min(v), max = max(v))
}
lam <- 1e-2
Tt <- (1 - lam) * U_ref
D <- (u_range[["max"]] - Tt) + 1
lexico <- function(u, n_edges) if (u >= Tt) u - D * n_edges else u - D * (n_free + 1)
f_lex <- function(xx) {
    d <- decode(xx)
    if (anyNA(d$w) || anyNA(d$G)) return(u_range[["min"]] - D * (n_free + 2) - 1e6)
    pen <- sum(d$w[d$w < 0]) + sum(d$G[d$G < 0]) - sum(d$w[d$w > 1] - 1) - sum(d$G[d$G > 1] - 1)
    if (pen < 0) return(u_range[["min"]] - D * (n_free + 2) + pen)
    d$w[d$w < 1e-4] <- 0; d$G[d$G < 1e-5] <- 0
    u <- ts_custom$func(graph_shortcut(pv4, 0.025, d$w, d$G))
    lexico(u, sum(d$G[free_mask] != 0))
}
cat(sprintf("T = %.5f, D = %.4f, f(reference, guarded) = %.5f, E = %d\n", Tt, D, f_lex(x_ref_g), count_edges(decode(x_ref_g)$G)))
set.seed(11)
nl <- nloptr::nloptr(x0 = x_ref_g, eval_f = function(v) -f_lex(v), lb = rep(0, length(x_ref_g)), ub = rep(1, length(x_ref_g)),
    opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 5e-8, xtol_abs = 5e-9, maxeval = 1000, print_level = 0))
sol <- param_to_solution(nl$solution, gc4, process = TRUE)
sol <- repair_graph(sol$hyp_weight, sol$trans_matrix, gc4)
U_sol <- U_of(sol$hyp_weight, sol$trans_matrix)
cat(sprintf("COBYLA: evaluations = %d, f = %.5f, edges %d -> %d, U %.5f -> %.5f, feasible (U >= T): %s\n",
    nl$iterations, -nl$objective, sum(G_ref != 0), sum(sol$trans_matrix != 0), U_ref, U_sol, U_sol >= Tt))
print(round(sol$trans_matrix, 4))
# and COBYLA started from the 2e-2 threshold-pruned graph (sparser support), same objective
r2 <- best_first_T(w_ref, G_ref, (1 - 2e-2) * U_ref, gc4)
x_sp <- guard_encode(as.numeric(create_start_params(gc4, w0 = w_ref, G0 = r2$G, sum_to_one_constraint = FALSE)), gc4)
cat(sprintf("sparse start: E = %d, U = %.5f, f = %.5f\n", count_edges(decode(x_sp)$G), U_of(w_ref, r2$G), f_lex(x_sp)))
nl2 <- nloptr::nloptr(x0 = x_sp, eval_f = function(v) -f_lex(v), lb = rep(0, length(x_sp)), ub = rep(1, length(x_sp)),
    opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 5e-8, xtol_abs = 5e-9, maxeval = 1000, print_level = 0))
sol2 <- param_to_solution(nl2$solution, gc4, process = TRUE)
sol2 <- repair_graph(sol2$hyp_weight, sol2$trans_matrix, gc4)
U_sol2 <- U_of(sol2$hyp_weight, sol2$trans_matrix)
cat(sprintf("COBYLA from sparse start: evaluations = %d, edges %d -> %d, U %.5f -> %.5f, feasible at lambda=1e-2: %s\n",
    nl2$iterations, sum(r2$G != 0), sum(sol2$trans_matrix != 0), U_of(w_ref, r2$G), U_sol2, U_sol2 >= Tt))
cat("\nDONE\n")
```

Output:

```
reference (stage-1 output): edges = 7  U_ref = 0.804825
       [,1]   [,2]   [,3]   [,4]
[1,] 0.0000 0.6567 0.3433 0.0000
[2,] 0.2834 0.0000 0.3354 0.3812
[3,] 1.0000 0.0000 0.0000 0.0000
[4,] 0.0000 1.0000 0.0000 0.0000
w_ref: 0.9951 0 0 0.0049

==== 8. stage-1 GA population is stored on the object; GA truncates suggestions ====
names(graph_optimal_example): hyp_weight trans_matrix constraints trial_success power solution global_search control global_output local_output start_graph
population dim: 200 3  popSize: 200
encoding nvars for its constraint: 3  (matches ncol: TRUE )
GA::ga with popSize 10 and 12 suggestion rows -> ERROR: number of items to replace is not a multiple of replacement length
GA::ga with popSize 10 and 10 suggestion rows -> ran, iter = 3
GA::ga with popSize 10 and 3 suggestion rows -> ran, iter = 3

==== 9. encoding a pruned graph: derived entries, the guard, and the objective ====
derived col per row (last free col):  4 4 4 3
derived entries after plain re-encoding:  0.000 0.381 0.000 0.000
derived weight (w[4]): 0.00494
any negative entry: FALSE  -> objective value: 0.804825  (U of reference on same sample: 0.804825 )
derived entries after guard:               0.000005 0.381207 0.000005 0.000005
derived weight after guard: 0.00494
objective after guard: 0.804825 ; edges counted: 7  (reference has 7 )
row params (0.3, 0.7): derived = 0 ; objective: 0.801075
row params (0.1, 0.9): derived = 0 ; objective: 0.79595
row params (0.7, 0.3): derived = 0 ; objective: 0.805125
row params (0.35, 0.65): derived = 0 ; objective: 0.8017

==== 10. threshold prune on the toy reference: budget spent in full ====
lambda=0.001  T=0.80402 : edges 7 -> 6, removed 1, loss 7.50e-05 (budget 8.05e-04, spent 9%), U=0.80475, within cap: TRUE
lambda=0.010  T=0.79678 : edges 7 -> 6, removed 1, loss 7.50e-05 (budget 8.05e-03, spent 1%), U=0.80475, within cap: TRUE
lambda=0.020  T=0.78873 : edges 7 -> 5, removed 2, loss 9.75e-03 (budget 1.61e-02, spent 61%), U=0.79507, within cap: TRUE
     [,1] [,2]  [,3]  [,4]
[1,]    0    1 0.000 0.000
[2,]    0    0 0.468 0.532
[3,]    1    0 0.000 0.000
[4,]    0    1 0.000 0.000

==== 11. stage-2 COBYLA within a support: lexicographic objective from the encoded reference ====
T = 0.79678, D = 1.2032, f(reference, guarded) = -7.61774, E = 7
COBYLA: evaluations = 118, f = -7.61759, edges 7 -> 7, U 0.80483 -> 0.80497, feasible (U >= T): TRUE
       [,1]   [,2]   [,3]   [,4]
[1,] 0.0000 0.6567 0.3433 0.0000
[2,] 0.2834 0.0000 0.3354 0.3812
[3,] 1.0000 0.0000 0.0000 0.0000
[4,] 0.0000 1.0000 0.0000 0.0000
sparse start: E = 5, U = 0.79507, f = -14.84683
COBYLA from sparse start: evaluations = 142, edges 5 -> 6, U 0.79507 -> 0.79890, feasible at lambda=1e-2: TRUE

DONE
```

Remarks on the output. In check 9 no plain re-encoding produced a negative derived entry: `sum()` accumulates in extended precision, so one minus a sum that is one in exact arithmetic rounds to exactly zero. The guard is therefore defensive; it costs nothing and removes a dependence on rounding behaviour. In check 11 the second COBYLA run starts from the five-edge graph produced at $\lambda=0.02$, which is infeasible at $\lambda=0.01$ ($U=0.79507<T=0.79678$, hence the score in the infeasible band, $-14.85$); COBYLA climbs in $U$, crosses the threshold by re-admitting one edge, and ends feasible with six edges. That is the infeasible branch behaving as designed, not an edge being created inside a feasible support: the first run, from the feasible reference, kept all seven edges.
