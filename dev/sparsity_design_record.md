# Design record: budgeted edge-count minimisation in `graph_optimise()`

Status: design specification, not yet implemented. Every design element is stated, followed by the reasoning behind it and the alternatives it was chosen over. Elements are marked **LOCKED** when the reasoning is considered settled and **OPEN** when evidence is still needed; open items are collected in section 11. Claims about existing package behaviour were checked by running the code; the script and its output are in the appendix.

Notation: $U(\mathbf w,\mathbf G)$ is the trial-success measure (expected gain) of a graph, estimated on the fixed matrix of simulated p-values; $\psi$ is the gain function compiled by `trial_success()`; $\mathbf G$ is the transition matrix; $E(\mathbf G)$ is the number of non-zero free entries of $\mathbf G$; $m$ is the number of hypotheses.

## 1. Problem

Optimised graphs frequently come back with several small non-zero transition weights. Each contributes almost nothing to expected gain, but together they make the graph hard to present: a clinical reviewer has to reason about every arrow. The package already has a post-hoc step, `prune_graph()`, that removes a weight or edge when doing so does not lower the estimated gain. That step never trades gain for simplicity.

The feature specified here lets the user state, once, how much expected gain they are prepared to give up in total for a simpler graph, and returns a graph with fewer edges whose total loss is guaranteed to stay within that amount. Simplicity is measured by the number of edges. The guarantee must hold under the design assumptions encoded in the p-value matrix, family-wise error control must be untouched, and with the feature switched off the package must behave exactly as it does today, including its consumption of random numbers.

## 2. How the package optimises today

The facts below constrain the design and were verified by reading and running the code.

`create_obj_func()` in `R/objective_function.R` is the single objective factory. Its closure decodes the parameter vector through `split_theta()`, `recover_full_weights()` and `recover_full_trans_matrix()`, returns a negative penalty for `NA` or out-of-range entries, zeroes hypothesis weights below $10^{-4}$ and transition entries below $10^{-5}$ without renormalising, runs `graph_shortcut()` (or its parallel twin) and returns the gain function's mean. The same factory serves `GA::ga()` for the global search, including the intermittent Nelder-Mead step (`optim = TRUE`), and the final `nloptr` COBYLA run.

The parameter encoding gives each row of $\mathbf G$ with $k_i$ free entries $k_i-1$ parameters; the last free entry is derived as one minus the rest. Rows therefore always sum to one, and the derived entry cannot be set to zero by moving a single parameter. `graph_constraint()` rejects a row with exactly one free entry and rejects an incomplete row whose fixed entries already sum to one (appendix check 1), so every incomplete row has $k_i\ge2$ and a positive amount of mass to distribute.

`param_to_solution(process = TRUE)` snaps entries below $10^{-5}$ to exact zero and entries in $(10^{-5},10^{-3})$ to the $\epsilon$-edge value $0.001$, then renormalises each row with `normalise_sum()`. `prune_graph()` runs `prune_hyp_weights()` and then `prune_edges()`, each in a fixed index order, and `.try_prune()` accepts a removal if the gain does not fall. `.redistribute_mass()` moves a dropped entry's mass proportionally onto the row's other free entries and, when those are all zero, spreads it uniformly, which creates $m-2$ new edges.

`graph_optimise()` runs the GA on a random subsample of `nsim_global` trials, COBYLA on a separate random subsample of `nsim_local` trials, chooses between the two results on the full sample with `choose_graph()`, prunes on the full sample, and only then computes the reported `$power` from the returned graph. That last ordering is deliberate: version 0.2.0 fixed a bug in which `$power$trial_success` described the graph before pruning rather than the one the user saw.

The mutation used inside the GA is a package closure (`.make_cauchy_mutation_multi()` in `R/mutation_helpers.R`), and the crossover is `GA`'s default local arithmetic crossover.

## 3. Specification

### 3.1 User-facing behaviour

`graph_optimise()` gains one argument, `gain_tolerance`, placed after `...` so it must be named. Its default is `NULL`, which switches the feature off. A number $\lambda\in[0,1]$ switches it on. The value is a cap on the total expected gain that may be lost, expressed as a fraction of the largest value the gain function can take. For a gain function valued in $[0,1]$, such as any probability of trial success, $\lambda$ is a cap in units of power: $\lambda=0.001$ means at most one tenth of a percentage point.

The returned `multigrain_graph_optimal` gains a `sparsity` element that is `NULL` when the feature is off and otherwise reports the edge counts, the internal price, the bound on the loss and the loss actually incurred during pruning. `print()` and `summary()` show these when present.

### 3.2 The cap becomes a fixed per-edge price

Let $U_{\max}=\max_{\mathbf r\in\{0,1\}^m}\psi(\mathbf r)$ be the largest value the gain function takes over all $2^m$ rejection patterns, and let
$$
P_G=\sum_{i=1}^{m}\max(k_i-1,\,0)
$$
be the number of edges that can be removed at all, where $k_i$ is the number of free entries in row $i$ of `trans_constraint`. Each row keeps at least one edge because the encoding forces rows to sum to one, so $P_G$ counts one fewer than the free entries in each row. For an unconstrained graph $P_G=m(m-2)$; in general it equals the number of free transition parameters in the encoding, which `recover_full_trans_matrix()` already computes (appendix check 6 confirms equality for a free and a constrained example).

The internal price per edge is
$$
c=\frac{\lambda\,U_{\max}}{P_G},
$$
and the objective maximised everywhere is
$$
f_\lambda(x)=U(x)-c\,E(x).
$$
Any optimiser that prefers a graph with $E_1$ edges over one with $E_0>E_1$ edges under $f_\lambda$ has accepted a loss in $U$ of at most $c\,(E_0-E_1)$. At most $P_G$ edges can be removed, so the total loss against any graph the search compared with is at most $cP_G=\lambda U_{\max}$. The cap therefore holds by construction, in a single pass, without knowing the unpenalised optimum.

**Rationale.** LOCKED. Three properties were required: a total cap rather than a per-edge price, a single optimisation pass, and an objective that is fixed before the search starts, because the package ranks candidates on common random numbers and a moving target would invalidate those comparisons. A cap defined relative to the unpenalised optimum $U^*$ would need $U^*$ first, and hence a second pass (section 6, alternative 1). Deriving the price from a quantity known in advance, $U_{\max}$, and from the number of edges that could possibly go, $P_G$, gives all three properties at once. The cost is conservatism: the bound assumes every removable edge goes, so at small $\lambda$ the price only removes edges that are individually very cheap. Section 4.1 quantifies this.

### 3.3 What counts as an edge

A free entry counts as an edge if it is at least $10^{-5}$ after decoding, evaluated on the same thresholded matrix the shortcut receives. Every non-zero free entry counts one, whatever its size; $\epsilon$-edges of $0.001$ count. Entries pinned by `graph_constraint()` are not counted in $E$ and not counted in $P_G$.

**Rationale.** LOCKED. The threshold $10^{-5}$ is the one the objective already applies, so an entry that counted during the search is still an edge after `param_to_solution()` snaps it (to exact zero below $10^{-5}$, to $0.001$ between $10^{-5}$ and $10^{-3}$), and an entry that did not count is exactly zero in the returned graph. Counting by size (entropy, an $\ell_q$ norm) would reward shrinking edges rather than removing them, which does not help a reader. Pinned entries are constants across every candidate, so counting them adds a constant to $E$ and changes nothing; excluding them keeps $P_G$ honest. The reported edge count of the returned graph does include pinned non-zero entries, because that is what a reader sees, and the free count is reported next to it.

Zeros in the returned graph are exact. Verified (appendix check 5): a parameter set to exactly zero stays exactly zero through `param_to_solution(process = TRUE)`, `repair_graph()` and `normalise_sum()`; the last scales free entries proportionally and adjusts the largest one, never a zero. No new snapping step is needed. `repair_graph()` fills a row whose free entries are all zero uniformly, but that state cannot arise from the encoding, because the derived entry is one minus the rest and `graph_constraint()` forbids incomplete rows whose fixed entries sum to one.

### 3.4 One objective everywhere

$f_\lambda$ is the fitness of the GA, the objective of its intermittent Nelder-Mead step, the objective of the final COBYLA run, the quantity `choose_graph()` compares on the full sample, and the acceptance criterion in `prune_graph()`. With `gain_tolerance = NULL` the price is zero and every one of these reduces to today's $U$.

**Rationale.** LOCKED. If the price applied only in a post-processing step, the search would first find a dense optimum and then remove edges without re-tuning the weights around the removal, and would spend the cap on whichever edges the removal order happened to visit first. Putting the price into the shared objective lets the global search find sparse regions, lets the local searches re-tune inside a support and refuse to re-create an edge, and makes the final acceptance test the same criterion the rest of the search used. It also removes a class of inconsistency: no stage can prefer a graph that another stage would reject.

### 3.5 A support-changing mutation move

When the price is positive, the GA's mutation closure gains a zeroing move. With probability $p_0$ per mutation call (proposed $p_0=0.2$, OPEN O3), instead of the usual Cauchy perturbation, it either sets one free transition parameter currently at or above $10^{-5}$ to exactly zero, or picks a row and rescales that row's parameters so they sum to $1-5\times10^{-6}$. The second variant puts the derived entry at $5\times10^{-6}$, below the counting threshold and strictly positive.

**Rationale.** LOCKED. Verified (appendix check 3): a truncated Cauchy perturbation lands below $10^{-5}$ with probability about $10^{-5}$ per mutated parameter, and the default crossover forms convex combinations of two parents, so a child entry is zero only where both parents are zero. Nelder-Mead and COBYLA are continuous methods and cross the threshold only by accident. The objective would reward fewer edges but nothing in the search would create them. The derived entry needs its own variant because it is not a parameter: an overshoot of the row's parameters gives a negative entry and a penalty (appendix check 2 shows $(0.6,0.6)$ returning $-0.2$), whereas rescaling to $1-5\times10^{-6}$ gives exactly $5\times10^{-6}$ and a normal evaluation. The move is added only when the price is positive, so the closure returned when the feature is off is the present one and draws exactly the random numbers it draws today.

### 3.6 Best-first pruning

When the price is positive, `prune_edges()` evaluates every remaining removable edge, removes the one whose loss in $U$ on the full sample is smallest, and repeats until the cheapest removal would cost more than $c$. Candidates that do not reduce the edge count are skipped: this covers the uniform fallback in `.redistribute_mass()`, which would replace one edge by $m-2$. Hypothesis-weight pruning is unchanged, since weights are not edges, and a weight still goes only if $U$ does not fall. The exact loss across accepted removals is recorded as `prune_loss`.

**Rationale.** LOCKED. The acceptance test `loss <= c` is the statement $f_\lambda(\text{candidate})\ge f_\lambda(\text{current})$ for a candidate with one fewer edge, so pruning uses the same objective as everything before it and nothing is paid twice: the price is a fixed number, not a running total. Best-first order spends the price on the edges that cost least; the present fixed index order could spend it on the first edge it happens to visit. Verified (appendix check 7): best-first at zero price from a 12-edge local optimum reaches the same 7-edge graph and the same $U$ as the current fixed-order prune, so the change of order is safe on that example. With the feature off the fixed-order loop runs unchanged.

### 3.7 Ties and Monte Carlo noise

Two graphs with the same edge count compare on $U$, as today. Whether an edge goes is decided by comparing $U$ before and after its removal on the same simulated trials, so only trials whose outcome changes contribute to the noise of that decision. Verified (appendix check 4): removing an $\epsilon$-edge of $0.001$ at $m=4$, $n_{\text{sim}}=10^4$ changed 59 of 10 000 trials; the paired standard error of the difference was $1.9\times10^{-4}$ against a marginal standard error of $3.4\times10^{-3}$ for $U$ itself, a factor of 18. At $10^6$ trials the paired error is about $2\times10^{-5}$. Because $c$ is a fixed number rather than a threshold on a sample estimate, nothing about the cap depends on which subsample a stage uses.

### 3.8 Reporting

The `sparsity` element:

| element | meaning |
|---|---|
| `gain_tolerance` | $\lambda$ as supplied |
| `edge_price` | $c$ |
| `u_max` | $U_{\max}$ |
| `n_removable` | $P_G$ |
| `n_edges` | `sum(trans_matrix != 0)` on the returned graph, pinned entries included |
| `n_edges_free` | non-zero free entries on the returned graph |
| `loss_bound` | $\lambda U_{\max}$ |
| `prune_loss` | exact loss of $U$ across the removals pruning accepted, on the full sample |

`$power` stays the last quantity computed, from the returned `hyp_weight` and `trans_matrix` on the full sample. The distance to the dense optimum is not reported because it is never computed; a user who wants it runs once more with `gain_tolerance = NULL` and compares.

**Rationale.** LOCKED. The reader must be able to see the unadjusted expected gain, the edge count, and what the simplification cost. The first two are direct. The third is available in two forms: the bound, which is guaranteed, and the prune-stage loss, which is exact. Reporting a comparison with the dense optimum would require the second pass this design avoids.

### 3.9 API placement

`gain_tolerance` is a direct argument of `graph_optimise()`, not a member of `multigrain_control`. `NULL` is off; a number is on; `0` is a meaningful setting (best-first removal at zero price, which removes only edges whose removal does not lower $U$ on the full sample). The documented recommendation is $10^{-3}$, with a note that on graphs of six or more hypotheses the cap begins to bind around $5\times10^{-3}$.

**Rationale.** LOCKED. The value changes what "optimal" means, as `trial_success` and `alpha` do; the control object holds tuning of the optimisers. The package moved `global_search` out of the control object and back to a direct argument in 0.3.0 for the same reason. A separate logical switch plus a numeric default would be two arguments for one idea.

### 3.10 Error control

Unchanged. A Bonferroni-based graphical procedure controls the family-wise error rate strongly for any weights summing to at most one and any transition matrix with non-negative entries, zero diagonal and row sums at most one. Zeroing an entry and renormalising its row preserves all three, and `is_graph_valid()` still gates `calc_power_pvals()`. Only the objective changes, not the class of procedures searched.

Degenerate cases the implementation must handle:

1. **Nothing removable** ($P_G=0$, for example $m=2$ or a fully pinned `trans_constraint`). `graph_optimise()` warns that `gain_tolerance` has no effect, sets the price to zero and still returns a `sparsity` element.
2. **$U_{\max}\le0$.** The gain function never rewards anything; abort with an informative error before optimising.
3. **Gain functions with negative values.** The existing invalid-encoding penalties (`-1e6` for `NA`, otherwise the sum of the violating entries) assume $U\ge0$; a valid graph with $U<0$ already ranks below an encoding with a $-10^{-17}$ violation. With a price, valid values extend down to $-\lambda U_{\max}$, so when the price is positive the penalties are shifted by $-(1+\lambda U_{\max})$ to stay below every valid value. Off-behaviour is untouched.
4. **A row reduced to one non-zero free entry during pruning.** Trying to remove it triggers the uniform fallback, which raises $E$; the candidate is skipped. The path in `.redistribute_mass()` with no recipients at all, which returns a row summing to less than one and makes `calc_power_pvals()` abort, is reachable only by bypassing `graph_constraint()` (appendix check 1 does so deliberately).
5. **Sparse two-cycles** $g_{ij}=g_{ji}=1$ hit the shortcut's `denom == 0` branch in `src/graph_shortcut.cpp`, which zeroes the row; that is the correct limit of the Bretz update and already exercised by every $m=2$ graph.
6. **Unreachable nodes** (zero weight, no incoming edge). Their outgoing edges do not affect $U$; at any positive price they are pruned to one edge, which is the interpretable result.

## 4. What to expect from the price

### 4.1 Conservatism

The bound holds even if every removable edge goes, so at small $\lambda$ the per-edge price is small. With $U_{\max}=1$:

| $m$ | $P_G$ | $c$ at $\lambda=10^{-3}$ | $c$ at $\lambda=5\times10^{-3}$ |
|---|---|---|---|
| 4 | 8 | $1.25\times10^{-4}$ | $6.25\times10^{-4}$ |
| 6 | 24 | $4.2\times10^{-5}$ | $2.1\times10^{-4}$ |
| 8 | 48 | $2.1\times10^{-5}$ | $1.0\times10^{-4}$ |

In appendix check 4 an $\epsilon$-edge of $0.001$ on a probability-valued gain cost $1.5\times10^{-3}$ at $m=4$, more than $c$ at $\lambda=10^{-3}$, so it would stay. In check 7, after the existing prune had already reduced a 12-edge local optimum to 7 edges, the price removed one more edge costing $7.5\times10^{-5}$, and the next cheapest cost $9.7\times10^{-3}$, so $\lambda\in\{10^{-3},5\times10^{-3},10^{-2}\}$ all returned the same 6-edge graph. The feature therefore removes edges that are individually negligible and leaves the rest; a user who wants more simplification raises $\lambda$ and the bound rises with it. Open item O1 asks for the distribution of per-edge losses on realistic problems, which would settle whether the recommended default should be $10^{-3}$ or $5\times10^{-3}$.

### 4.2 Landscape

$E$ is piecewise constant in the parameters, so $f_\lambda$ is $U$ with a downward step of $c$ each time an entry crosses $10^{-5}$. Inside a fixed support the landscape is exactly today's, so Cauchy mutation, crossover, Nelder-Mead and COBYLA behave as they do now there. The local searches re-tune weights inside a support and, because $f_\lambda$ drops when an entry re-enters, do not re-create edges. Support changes come from the zeroing move during the global search and from best-first pruning at the end. With `global_search = FALSE` only pruning changes the support: COBYLA from a dense start does not become sparse on its own, and the documentation says so.

### 4.3 Cost

One pass. The extra work is $2^m$ evaluations of the compiled gain function for $U_{\max}$ (at most 4096 rows for $m\le12$, no random numbers), the mutation move (cheaper than a Cauchy perturbation), and best-first pruning at up to $E^2/2$ shortcut evaluations on the full sample, which is seconds to a couple of minutes at $10^6$ trials.

## 5. Specification for implementation

### 5.1 Price and counts

```r
# R/objective_function.R (new helpers)

#' Number of removable free edges under a graph constraint
#' @noRd
.n_removable_edges <- function(trans_constraint) {
    k <- rowSums(is.na(trans_constraint))
    sum(pmax(k - 1L, 0L))
}

#' Exact maximum of the gain function over all rejection patterns
#' @noRd
.trial_success_max <- function(trial_success) {
    m <- trial_success$m
    patterns <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), m)))
    vals <- vapply(
        seq_len(nrow(patterns)),
        function(i) trial_success$func(patterns[i, , drop = FALSE]),
        numeric(1)
    )
    max(vals)
}

#' Per-edge price implied by a total cap
#' @noRd
.edge_price <- function(gain_tolerance, trial_success, graph_constraint) {
    if (is.null(gain_tolerance)) {
        return(0)
    }
    n_removable <- .n_removable_edges(graph_constraint$trans_constraint)
    if (n_removable == 0L) {
        return(0)
    }
    u_max <- .trial_success_max(trial_success)
    if (u_max <= 0) {
        cli::cli_abort(
            "{.arg trial_success} never takes a positive value; \\
            {.arg gain_tolerance} cannot be applied."
        )
    }
    gain_tolerance * u_max / n_removable
}
```

### 5.2 Objective

Complete replacement for the closure factory. The decoding and penalty code is unchanged in content; the additions are the `edge_price` argument, the free-entry mask, the penalty shift, and the last four lines.

```r
create_obj_func <- function(
    m,
    power_criterion,
    hyp_constraint,
    trans_constraint,
    pvals,
    alpha = 0.025,
    num_threads = 1L,
    edge_price = 0
) {
    force(power_criterion)
    force(hyp_constraint)
    force(trans_constraint)
    force(alpha)
    force(pvals)
    force(num_threads)
    force(edge_price)

    use_parallel <- num_threads >= 2L
    free_mask <- is.na(trans_constraint)
    n_free <- sum(free_mask)
    # Shift for invalid encodings so they stay below every valid U - c * E.
    # Zero when the feature is off, so existing penalties are returned as is.
    penalty_shift <- if (edge_price > 0) 1 + edge_price * n_free else 0

    function(x) {
        theta <- split_theta(x, hyp_constraint)
        hyp_weight <- recover_full_weights(theta$w_pars, hyp_constraint)
        trans_matrix <- recover_full_trans_matrix(
            theta$g_pars,
            trans_constraint
        )

        if (anyNA(hyp_weight) || anyNA(trans_matrix)) {
            return(-1e6 - penalty_shift)
        }
        if (any(hyp_weight < 0)) {
            return(sum(hyp_weight[hyp_weight < 0]) - penalty_shift)
        }
        if (any(trans_matrix < 0)) {
            return(sum(trans_matrix[trans_matrix < 0]) - penalty_shift)
        }
        if (any(hyp_weight > 1)) {
            return(-sum(hyp_weight[hyp_weight > 1]) - penalty_shift)
        }
        if (any(trans_matrix > 1)) {
            return(-sum(trans_matrix[trans_matrix > 1]) - penalty_shift)
        }

        hyp_weight[hyp_weight < 1e-4] <- 0
        trans_matrix[trans_matrix < 1e-5] <- 0

        rej_matrix <- if (use_parallel) {
            graph_shortcut_parallel(
                pvals = pvals,
                alpha = alpha,
                w = hyp_weight,
                G = trans_matrix,
                num_threads = num_threads,
                grain_size = -1L
            )
        } else {
            graph_shortcut(
                pvals = pvals,
                alpha = alpha,
                w = hyp_weight,
                G = trans_matrix
            )
        }

        u <- power_criterion(rej_matrix)
        if (edge_price > 0) {
            u - edge_price * sum(trans_matrix[free_mask] != 0)
        } else {
            u
        }
    }
}
```

With `edge_price = 0` every returned value is arithmetically the value returned today (`x - 0` is exact for finite `x`), and no extra random numbers are drawn.

### 5.3 Objective on the full sample, and `choose_graph()`

`.graph_optimise_ga()` and `.graph_optimise_local()` keep `ga_trial_success` and `local_trial_success` (raw $U$ on the full sample) and add `ga_objective` and `local_objective`, computed as $U-cE$ from the processed graph so that $E$ uses the same zeroed entries the user will see. `choose_graph()` compares the objectives. With a zero price they equal the trial-success values, so its decisions are unchanged.

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
        # Half the zeroing moves target a parameter, half a derived entry.
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

The inner `cauchy_mutation` is the present closure verbatim and is what the factory returns when `p_zero` is zero. The row rescaling relies on `recover_full_trans_matrix()` computing the derived entry as `1 - sum(G[i, ])`; appendix check 2 shows the result is exactly $5\times10^{-6}$. When the row's parameters sum to more than one before the move the rescaling moves mass off the derived entry and onto the explicit ones; the `pmin` clamp keeps each within its bound and the objective's penalty handles the rare remainder.

### 5.5 Best-first pruning

```r
# R/post_optim_processing.R: used by prune_graph() when edge_price > 0.
# The existing fixed-order prune_edges() is kept unchanged for edge_price == 0.

.prune_edges_best_first <- function(
    pvals,
    hyp_weight,
    trans_matrix,
    trial_success,
    fixed_edge,
    edge_price,
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
            # Skip candidates that do not reduce the edge count (uniform
            # fallback) or that leave an invalid row (no recipients).
            if (sum(G_try != 0) >= n_edges) next
            if (abs(sum(G_try[i, ]) - 1) > tolerance) next
            u_try <- u_of(G_try)
            if (is.null(best) || u_try > best$u) {
                best <- list(G = G_try, u = u_try)
            }
        }
        if (is.null(best) || u_best - best$u > edge_price) break
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

`prune_graph()` gains `edge_price = 0`, forwards it, dispatches to this function when it is positive, and its return list gains `prune_loss` (zero when off).

### 5.6 `graph_optimise()` and the returned object

In `graph_optimise()`: validate `gain_tolerance` with `rlang::check_number_decimal(gain_tolerance, min = 0, max = 1, allow_null = TRUE)`; compute `edge_price <- .edge_price(gain_tolerance, trial_success, graph_constraint)`; warn once if `gain_tolerance` is non-`NULL` and the price is zero because nothing is removable; pass `edge_price` to `.graph_optimise_ga()`, `.graph_optimise_local()` and `prune_graph()`; in `.graph_optimise_ga()`, build the mutation closure with `p_zero = if (edge_price > 0) 0.2 else 0` and `param_rows = .g_param_rows(...)`; after `final_power`, assemble `sparsity` when `gain_tolerance` is non-`NULL`. `graph_optimal()` and `new_graph_optimal()` gain `sparsity = NULL`, appended after `start_graph`, so `names(result)` gains one entry even when the feature is off. `print.multigrain_graph_optimal()` adds, when `sparsity` is non-`NULL`:

```
Edges: 6 (5 free of 8 removable); edge price 1.25e-04
Total loss of trial success bounded by 1e-03 (prune stage: 7.5e-05)
```

`summary()` prints the same two lines under the power block.

## 6. Alternatives considered

1. **A threshold against the unpenalised optimum in two stages.** Run the present pipeline, take its graph as reference $x_{\text{ref}}$, then re-run the local search and pruning (optionally the GA) with a lexicographic objective ranking feasibility $U\ge(1-\lambda)U_{\text{ref}}$ first, fewer edges second and $U$ third, with the threshold recomputed on each stage's sample. This certifies the cap against the dense optimum and can spend the whole budget. Rejected: it needs a second pass, three thresholds threaded through three samples, and a fallback rule, and it roughly doubles runtime when the GA is repeated. The single-pass price gives the same guarantee against a slightly different reference ($U_{\max}$ rather than $U^*$) at no extra cost.
2. **An uncalibrated per-edge price.** One pass, but no total cap, and a price whose meaning depends on the scale of $\psi$. Rejected.
3. **A single-pass price calibrated on a cheap lower bound of $U^*$**, $c=\lambda U_{\text{seed}}/P_G$. Certifies the cap against $U^*$ rather than $U_{\max}$, but is smaller still and needs a seed evaluation before the objective exists. Rejected as strictly more conservative for no practical gain.
4. **A multi-objective GA** over $(U,-E)$ with a choice from the Pareto front afterwards. Needs a new dependency and rewrites the global search. Rejected.
5. **Counting pinned entries.** Adds a constant, changes nothing in the search, confuses the report. Rejected.
6. **Entropy or $\ell_q$ sparsity measures.** Reward shrinking edges rather than removing them. Rejected.
7. **A lexicographic scalarisation with a large edge weight** so that fewer edges always win. Something must then stop it from removing everything, which is a threshold in disguise. Superseded by the price.
8. **Placing the setting on `multigrain_control`.** Rejected for the reason in 3.9.

## 7. Implementation plan

Each step ends with a gate that must pass before the next starts. Test files are run one at a time with `testthat::test_file()`.

1. **Helpers and objective** (`R/objective_function.R`). Add `.n_removable_edges()`, `.trial_success_max()`, `.edge_price()` and the `edge_price` argument as in 5.1 and 5.2. Gate: `test-objective_function.R` passes unchanged; new tests: for 200 random encodings at $m\in\{3,4\}$ the closure with `edge_price = 0` returns values identical (`expect_identical`) to a copy of the pre-change closure kept in the test file; with a positive price it returns $U-cE$ with $E$ computed independently from `param_to_solution(process = TRUE)`; every invalid encoding scores below $-\lambda U_{\max}$; `.n_removable_edges()` equals the parameter count of `create_start_params()` for three constraints; `.trial_success_max()` returns 1, 4 and 1 for the three functions in appendix check 6.
2. **Mutation move** (`R/mutation_helpers.R`) as in 5.4, plus `.g_param_rows()`. Gate: `test-mutation_helpers.R` passes; new tests: with `p_zero = 0` the factory returns a closure whose outputs under `set.seed()` are identical to the present closure's; with `p_zero = 1` every call either sets a free transition parameter to exactly 0 or scales a row to $1-5\times10^{-6}$ within `1e-12`; the row map matches `recover_full_trans_matrix()`'s parameter order for a constrained example.
3. **Pruning** (`R/post_optim_processing.R`) as in 5.5; `prune_graph()` gains `edge_price` and returns `prune_loss`. Gate: `test-post_optim_processing.R` passes unchanged; new tests: on the 4-hypothesis fixture in that file, a zero price reproduces the current 7-edge result of the `gamma = 1` prune; a positive price never accepts a candidate that raises the edge count; `prune_loss <= edge_price * n_removed`; no row is left with a sum different from one.
4. **Pipeline** (`R/optimisation.R`, `R/choose_graph.R`). Thread `edge_price` through `.graph_optimise_ga()` and `.graph_optimise_local()`; add `ga_objective` and `local_objective`; `choose_graph()` compares them. Gate: `test-optimisation.R` and `test-choose_graph.R` pass with existing snapshots untouched.
5. **Object and methods** (`R/graph_optimal.R`). Add `sparsity`; print and summary lines. Gate: `test-graph_optimal.R` passes; new snapshot with a non-`NULL` `sparsity`.
6. **User-facing argument** (`R/optimisation.R`, roxygen). `gain_tolerance = NULL`; validation; warning for $P_G=0$; error for $U_{\max}\le0$; documentation of the semantics, the conservatism, the recommended $10^{-3}$ and the note about `global_search = FALSE`. Gate: end-to-end test at $m=4$, $n_{\text{sim}}=10^4$, `gain_tolerance = 5e-3`: returned graph valid, `sparsity$n_edges` at most the count from a `NULL` run with the same seed, `$power$trial_success` identical to a fresh `calc_power_pvals()` on the returned graph, and `sparsity$prune_loss <= sparsity$loss_bound`.
7. **Bit-identity gate** (section 8). Gate: passes on the branch against a fixture generated before any change.
8. **Regenerate `data/graph_optimal_example.rda`** with `data-raw/graph_optimal_example.R` (the object gained an element), run `devtools::document()`, update `NEWS.md` (section 10). Gate: `R CMD check` clean locally.

## 8. Test plan

Beyond the per-step gates above:

**Bit-identical when disabled.** Before any change, run the recipe below and save the result and the RNG state. After the change, the test re-runs the recipe and compares.

```r
# tests/testthat/data/make_sparsity_baseline.R (run before the change; commits the RDS)
set.seed(20260906)
pvals <- multigrain::simulate_pvalues(
    c(0.93, 0.91, 0.90, 0.85),
    corr_matrix = matrix(0.2, 4, 4) + diag(0.8, 4),
    nsim = 1e4
)
ts <- multigrain::trial_success(0.25 * (2 * (r1 && r2) + r1 * r3 + r2 * r4), verbose = "silent")
ctrl <- multigrain::multigrain_control() |>
    multigrain::control_global(maxiter = 30, popSize = 40, run = 10) |>
    multigrain::control_local(maxeval = 500)
set.seed(1)
res <- multigrain::graph_optimise(
    pvals = pvals,
    graph_constraint = multigrain::graph_constraint_free(4),
    trial_success = ts,
    control = ctrl,
    verbose = "silent"
)
saveRDS(
    list(result = res, seed_after = .Random.seed, versions = c(
        R = R.version.string, GA = as.character(packageVersion("GA")),
        nloptr = as.character(packageVersion("nloptr"))
    )),
    "tests/testthat/data/sparsity_baseline.rds"
)
```

```r
# tests/testthat/test-sparsity_baseline.R
test_that("graph_optimise() is bit-identical to the baseline when gain_tolerance is NULL", {
    skip_on_cran()
    base <- readRDS(test_path("data", "sparsity_baseline.rds"))
    # same recipe as make_sparsity_baseline.R
    ...
    res <- graph_optimise(..., control = ctrl, verbose = "silent")
    keep <- setdiff(names(res), "sparsity")
    expect_identical(unclass(res)[keep], unclass(base$result)[keep])
    expect_identical(.Random.seed, base$seed_after)
    expect_null(res$sparsity)
})
```

The test is skipped on CRAN because nloptr and GA binaries may differ across platforms in the last bits; it is the gate for step 7 locally and in the Ubuntu CI job. The fixture records the R, GA and nloptr versions so a mismatch is diagnosable.

**Feature tests.** In addition to the per-step gates: `gain_tolerance = 0` runs, returns a non-`NULL` `sparsity`, and never lowers `$power$trial_success` below the `NULL` run's value with the same seed; `num_threads = 2` gives the same `sparsity$n_edges` as `num_threads = 1` for the same seed (appendix check 4 confirms the parallel shortcut returns identical rejections); a constrained example with pinned non-zero edges reports `n_edges > n_edges_free`; `gain_tolerance = 1e-3` on `graph_constraint_free(2)` warns that nothing is removable.

## 9. Adversarial checks for review

- **A1.** Run `graph_optimise()` before and after the change for seeds 1 to 5, $m\in\{2,3,4\}$, constrained and free, `global_search` on and off, `num_threads` 1 and 2, `gain_tolerance = NULL`. Any non-identical element other than `sparsity`, or any difference in `.Random.seed` afterwards, fails the off-by-default requirement.
- **A2.** Use a gain function with negative values, such as `r1 - r2`, and a price. Confirm no invalid encoding ever outranks a valid graph in the GA's final population.
- **A3.** Constant gain function and one that is never positive: confirm the error path in `.edge_price()`, not a silent zero or `NaN`.
- **A4.** $m=2$ and a fully pinned `trans_constraint`: confirm the warning and that `sparsity$n_removable` is 0.
- **A5.** Construct a GA population where the zeroing move rescales a row whose parameters sum to more than one. After the move the row sums to $1-5\times10^{-6}$ unless clamped; confirm no parameter exceeds its bound and the objective is finite.
- **A6.** After `param_to_solution(process = TRUE)` snaps an entry to $0.001$, confirm pruning can remove it when its loss is below $c$, and that it is counted in `n_edges_free` when it stays.
- **A7.** Try to make best-first pruning loop for ever: every accepted step reduces the edge count by one, so at most $P_G$ iterations; verify with a `gain_tolerance = 1` run.
- **A8.** Confirm `choose_graph()` picks the GA graph when the local one is invalid and a price is set, and that the comparison uses the objective, not raw $U$.
- **A9.** The pre-existing abort when `.redistribute_mass()` has no recipients (appendix check 1): confirm it is unreachable through `graph_constraint()` and that the validity skip in best-first pruning makes it unreachable there as well.
- **A10.** Noise floor at scale: at $m=8$, $n_{\text{sim}}=10^6$, $\lambda=10^{-3}$, $c\approx2\times10^{-5}$ is close to the paired standard error for an edge that flips 0.1% of trials. Check whether pruning's accept/reject at the margin flips between two independent p-value matrices; if it does, the documentation must say that $\lambda$ below $5\times10^{-3}$ at $m\ge8$ is at the noise floor. This is the design's weakest point.
- **A11.** `global_search = FALSE` with a price: confirm the result differs from the `NULL` run only through pruning, and that the documentation says so.
- **A12.** Confirm the intermittent Nelder-Mead step inside `GA::ga()` receives the priced fitness (it calls the same `fitness` slot) and does not re-create edges: count zeros before and after `optim` steps on a logged run.

## 10. Proposed `NEWS.md` wording

Under `# multigrain (development version)`:

```
## New functionality

* `graph_optimise()` gains `gain_tolerance`, a total cap on the trial-success
  value the optimiser may give up in exchange for a graph with fewer edges.
  The cap is expressed as a fraction of the largest value the trial-success
  function can take, so for probability-valued functions it is a cap in units
  of power. Internally it becomes a fixed per-edge price
  `gain_tolerance * max(psi) / n_removable` that is part of the objective used
  by the global search, the local search and the pruning step, so the total
  loss is bounded by the cap however many edges are removed. When
  `gain_tolerance` is `NULL` (the default) behaviour is unchanged.
* When `gain_tolerance` is set, the returned `multigrain_graph_optimal` gains a
  `sparsity` element reporting the edge count, the price, the bound on the
  total loss and the exact loss incurred during pruning; `print()` and
  `summary()` show these.
* When `gain_tolerance` is set, edge pruning is best-first: the edge whose
  removal costs least is removed first.
```

No bug fixes are part of this change. Two pre-existing issues found while reading the code are recorded as open items O2 and O4 and should be fixed separately so that the off-by-default gate stays clean.

## 11. Open items

- **O1. Conservatism of the price.** At $\lambda=10^{-3}$ and $m\ge6$ the price removes only edges that flip a few dozen trials per million. Evidence that would close it: on realistic problems of five to eight hypotheses, the distribution of per-edge losses of the edges left after the existing prune. If most lie between $10^{-4}$ and $10^{-3}$, the documented recommendation should be $5\times10^{-3}$, or the two-stage threshold (section 6, alternative 1) becomes worth its cost.
- **O2. `global_opt_power` / `local_opt_power`.** The 0.2.0 news entry says the returned object stores the pre-pruning power of the global and local solutions. The constructor in `R/graph_optimal.R` has no such elements; `.graph_optimise_ga()` and `.graph_optimise_local()` compute them as `ga_subset_power` and `local_subset_power`, and `graph_optimise()` discards them. Either reinstate them or correct the news entry. Separate change.
- **O3. $p_0$ and the split between parameter zeroing and row rescaling.** Proposed $0.2$ and one half. Evidence to close: on $m=4$ examples, the fraction of GA generations in which the elite's edge count decreases, for $p_0\in\{0.1,0.2,0.4\}$.
- **O4. Pre-existing abort in `prune_edges()`** when a row has no free recipients (appendix check 1). Unreachable through `graph_constraint()`; fix separately if `prune_graph()` is ever exported.
- **O5. Whether `gain_tolerance = 0` is documented as a supported mode.** It is well defined (best-first, zero price) and cheap; giving it a name is a documentation decision.

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
