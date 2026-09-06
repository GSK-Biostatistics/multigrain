# Design record: budgeted edge-count minimisation in `graph_optimise()`

Status: design only. Nothing under `R/`, `src/`, `tests/`, `man/` or `NAMESPACE` changes in this record. A later session implements; a third reviews adversarially.

Format used here (there were no earlier design records in the repository, so this one sets the convention): every decision is stated, the alternatives considered are listed with the reason each was rejected, and the decision is marked **LOCKED** or **OPEN**. Numbered claims about existing code were verified by running it; the script and its output are in the appendix.

## 1. Purpose

Optimised graphs come back with many small non-zero transition weights that contribute almost nothing to expected gain but make the graph hard to present in a protocol. `prune_graph()` already removes weights and edges when the gain does not fall. This design adds an explicit, capped willingness to give up a small amount of gain for a simpler graph: the user states a single total $\lambda$ they are prepared to lose, and the optimiser returns a graph with fewer edges whose total loss is bounded by that amount.

Throughout, $U(\mathbf w,\mathbf G)$ is the trial-success measure (expected gain) estimated on the fixed p-value matrix, $\mathbf G$ the transition matrix, $E(\mathbf G)=\lVert\mathbf G\rVert_0$ restricted to free entries, and $m$ the number of hypotheses.

## 2. Where this record departs from the brief

The brief was written before the repository was read. Five points differ and are flagged rather than worked around.

1. The brief asks for "the locked-decision format already used in this project's design records". There are none in the repository (`dev/` holds scripts, diagrams and old vignettes). This record defines the format.
2. The brief says `$global_opt_power` and `$local_opt_power` are populated on a `multigrain_graph_optimal`, and `NEWS.md` (0.2.0, bug fixes) says the same. The constructor in `R/graph_optimal.R` has no such slots. `.graph_optimise_ga()` and `.graph_optimise_local()` compute `ga_subset_power` and `local_subset_power`, and `graph_optimise()` discards them. See OPEN item O2.
3. The brief's problem statement is a threshold against $U^*$, "the maximum achievable expected gain", with $\lambda$ a fraction of it. $U^*$ is not available until an unpenalised optimisation has finished, so that statement needs a second pass. After discussion the maintainer chose a single-pass realisation of the same cap: a per-edge price derived from $\lambda$ and the number of removable edges, so that the total loss is bounded by $\lambda$ times the largest value the gain function can take. Decisions 3 and 7 below are rewritten accordingly. The threshold formulation is recorded under rejected alternatives.
4. `create_obj_func*` is a single function, `create_obj_func()` in `R/objective_function.R`, not a family.
5. "Default $10^{-3}$" and "off by default" cannot both hold for one numeric argument. Resolved as `NULL` = off, with $10^{-3}$ the documented recommended value.

## 3. Decisions

**D1. Sparsity measure is edge count, $E=\lVert\mathbf G\rVert_0$ over free entries.** LOCKED (brief).

**D2. Every non-zero free edge counts one, including $\epsilon$-edges.** LOCKED (brief). Numerically, an entry counts if it survives the same `< 1e-5` zeroing that the objective already applies before calling the shortcut (`R/objective_function.R`, line 67). The $\epsilon$-edge snap in `param_to_solution(process = TRUE)` maps $(10^{-5},10^{-3})$ to $0.001$, so an entry counted during the search is still an edge in the returned graph, and an entry not counted is exactly zero there.

**D3. A single total cap $\lambda$, realised as a derived per-edge price.** LOCKED (rewritten with the maintainer). The user supplies $\lambda\in[0,1]$. Internally,
$$
c=\frac{\lambda\,U_{\max}}{P_G},\qquad
U_{\max}=\max_{\mathbf r\in\{0,1\}^m}\psi(\mathbf r),\qquad
P_G=\sum_{i=1}^{m}\max(k_i-1,\,0),
$$
where $\psi$ is the gain function compiled by `trial_success()` and $k_i$ is the number of `NA` entries in row $i$ of `trans_constraint`. $P_G$ is the number of edges that can be removed at all: every row keeps at least one edge, because the encoding forces each row to sum to exactly one (the last free entry is derived as one minus the rest, `recover_full_trans_matrix()` lines 194–204). For an unconstrained graph $P_G=m(m-2)$. Verified: the formula equals the encoding's free transition-parameter count for `graph_constraint_free(4)` (8) and for the constrained example in `tests/testthat/data/test_data.R` (2); see appendix check 6. Since at most $P_G$ edges can go and each costs at most $c$ against any graph the search compared it with, the total loss is bounded by $cP_G=\lambda U_{\max}$. For probability-valued $\psi$, $U_{\max}=1$ and $\lambda$ reads directly as a cap in units of power.

Alternatives: absolute per-edge price (no total cap; rejected by the brief); price in Monte Carlo standard-error units (rejected by the brief); a threshold $(1-\lambda)U^*$ (needs a second pass; rejected by the maintainer on runtime, see section 6).

**D4. The price lives in the general objective used everywhere.** LOCKED (brief). `create_obj_func()` returns $f_\lambda(x)=U(x)-c\,E(x)$ when a price is set; the same closure factory serves `GA::ga()` (including its intermittent Nelder-Mead step, `optim = TRUE`), the final `nloptr` COBYLA run, and the same scalar is used by `choose_graph()` and by the acceptance test in `prune_graph()`.

**D5. Off by default; identical behaviour and RNG consumption when off.** LOCKED (brief). With `gain_tolerance = NULL` no new code path executes: `create_obj_func()` returns its present closure, the mutation closure is the present one, prune runs its present fixed-order loop, and no random draws are added. The gate is in section 8.

**D6. Package feature, not a pipeline experiment.** LOCKED (brief).

**D7. Symbol $\lambda$ is the user's total cap; $c$ is the internal price.** LOCKED (rewritten). The user-facing argument is `gain_tolerance`. $\epsilon$, $\tau$, $\delta$ remain taken.

**D8. Encoding unchanged.** LOCKED (brief). The zeroing move in D11 works within the existing parameterisation.

**D9. Pinned entries are not counted and not removable.** LOCKED. Entries fixed by `graph_constraint()` are constants across all candidates, so counting them would add a constant to $E$; they are excluded from $E$ and from $P_G$. The reported edge count of the returned graph includes them, because that is what a reader sees, and the free count is reported alongside.

**D10. Ties resolve by higher $U$ on the common sample.** LOCKED. Two graphs with the same $E$ compare on $U$ exactly as today. Because all candidates in a run are evaluated on the same fixed p-value matrix, the decision to drop an edge is a paired comparison, and only simulated trials whose outcome changes contribute to its noise. Verified at $m=4$, $n_{\text{sim}}=10^4$ (appendix check 4): removing an $\epsilon$-edge of $0.001$ changed 59 of 10 000 trials; the paired standard error of the difference was $1.9\times10^{-4}$ against a marginal standard error of $3.4\times10^{-3}$ for $U$ itself, a ratio of 18. At $10^6$ trials the paired error scales to about $2\times10^{-5}$.

**D11. The GA gains a support-changing mutation move, active only when a price is set.** LOCKED. Verified (appendix check 3): a truncated Cauchy perturbation lands below $10^{-5}$ with probability about $10^{-5}$ per mutated parameter, and the default real-valued crossover `gareal_laCrossover` forms convex combinations of two parents, so a child entry is zero only if both parents are zero there. Nelder-Mead and COBYLA are continuous optimisers and do not cross the threshold on purpose. Without a dedicated move the GA would see the price but could not act on it. The move, with probability $p_0$ per mutation call (proposed $p_0=0.2$, OPEN O3): pick one free transition parameter currently at or above $10^{-5}$ uniformly at random and set it to exactly zero; or, with probability proportional to the number of rows, pick a row and rescale its parameters so that they sum to $1-5\times10^{-6}$, which puts the derived entry at $5\times10^{-6}$, below the threshold and strictly positive. Verified (appendix check 2): parameters $(0.3,\,0.7-5\times10^{-6})$ give a derived entry of exactly $5\times10^{-6}$ and the closure evaluates normally; an ordinary overshoot such as $(0.6,0.6)$ gives $-0.2$ and the closure returns $-0.2$ as its penalty, and a one-ulp overshoot rounds to exactly zero because `sum()` accumulates in extended precision. The move needs a map from parameter index to row, computed once from `trans_constraint`.

**D12. `prune_graph()` stays and becomes best-first when a price is set.** LOCKED. `.try_prune()` accepts on $f_\lambda$ rather than on raw $U$; with no price this is the existing test. With a price, `prune_edges()` evaluates every remaining removable edge, removes the one whose loss in $U$ is smallest, and repeats until the cheapest removal would cost more than $c$. Candidates that do not reduce $E$ (the uniform fallback in `.redistribute_mass()` when a row's other free entries are all zero) are skipped. Hypothesis-weight pruning is unchanged: weights are not edges (D1), and a weight is still removed only if $U$ does not fall. Because every stage pays the same fixed price, nothing is spent twice.

**D13. API: `graph_optimise(..., gain_tolerance = NULL)`.** LOCKED (agreed with the maintainer). Named-only, after `...`, checked with `rlang::check_number_decimal(min = 0, max = 1, allow_null = TRUE)`. `NULL` is off. A numeric value switches the feature on; `0` is a meaningful setting (best-first removal at zero price, which removes only edges whose removal does not reduce $U$ on the full sample). Rationale: the value changes what "optimal" means, like `trial_success` and `alpha`, and the package removed a control-object switch (`control_global_search()`) in 0.3.0 in favour of a direct argument. Alternatives: `control_sparsity()` on `multigrain_control` (tuning, not problem definition; rejected); a separate logical switch plus a numeric default (two arguments for one idea; rejected).

**D14. Reporting.** LOCKED. A new `sparsity` element on `multigrain_graph_optimal`, `NULL` when off, otherwise a list:

| element | meaning |
|---|---|
| `gain_tolerance` | $\lambda$ as supplied |
| `edge_price` | $c$ |
| `u_max` | $U_{\max}$ |
| `n_removable` | $P_G$ |
| `n_edges` | `sum(trans_matrix != 0)` on the returned graph, pinned entries included |
| `n_edges_free` | non-zero free entries on the returned graph |
| `loss_bound` | $\lambda U_{\max}$ |
| `prune_loss` | exact loss of $U$ across the removals prune accepted, on the full sample |

`$power` remains the last quantity computed, from the returned `hyp_weight` and `trans_matrix` on the full sample, which is the 0.2.0 fix and keeps the reported trial success equal to the graph the user sees. `print()` and `summary()` each gain two lines when `sparsity` is non-`NULL`. The distance to the dense optimum is not reported because it is never computed; the documented way to obtain it is a second run with `gain_tolerance = NULL`.

**D15. Error control is unaffected.** LOCKED. A Bonferroni-based graphical procedure controls the family-wise error rate strongly for any $\mathbf w$ with $\sum w_i\le1$ and any $\mathbf G$ with non-negative entries, zero diagonal and row sums at most one. Zeroing an entry and renormalising the row keeps all three properties, and `is_graph_valid()` still gates `calc_power_pvals()`. The search space has not changed, only the objective.

## 4. Answers to the brief's questions

**(a) Making the cap operational.** No reference value is needed. $U_{\max}$ is exact and computed once from the $2^m$ rejection patterns by calling the compiled `trial_success$func` on one-row matrices (verified: conjunctive 1, `r1 + r2 + r3 + r4` 4, the custom example 1). $P_G$ is read off `trans_constraint`. $c$ is fixed before the first evaluation, so the objective is stationary within every optimiser run and across runs, and common-random-number ranking holds. One caveat is honest and important: the price is conservative. It is calibrated so that the cap holds even if every removable edge goes, so at small $\lambda$ it only removes edges that are cheap individually. The table gives $c$ for $U_{\max}=1$.

| $m$ | $P_G$ | $c$ at $\lambda=10^{-3}$ | $c$ at $\lambda=5\times10^{-3}$ |
|---|---|---|---|
| 4 | 8 | $1.25\times10^{-4}$ | $6.25\times10^{-4}$ |
| 6 | 24 | $4.2\times10^{-5}$ | $2.1\times10^{-4}$ |
| 8 | 48 | $2.1\times10^{-5}$ | $1.0\times10^{-4}$ |

In appendix check 4 an $\epsilon$-edge of 0.001 on a probability-valued gain cost $1.5\times10^{-3}$ at $m=4$, more than $c$ at $\lambda=10^{-3}$, so it would stay. In check 7, after the existing prune had already taken a 12-edge local optimum to 7 edges, the price removed one further edge costing $7.5\times10^{-5}$ and the next cheapest cost $9.7\times10^{-3}$, so $\lambda\in\{10^{-3},5\times10^{-3},10^{-2}\}$ all gave the same 6-edge graph. The recommended default in documentation is $10^{-3}$, with the note that on graphs of six or more hypotheses $5\times10^{-3}$ is where the cap starts to bind.

**(b) Landscape.** $E$ is piecewise constant in the parameters, so $f_\lambda$ is $U$ with a downward step of $c$ each time an entry crosses $10^{-5}$. Within a fixed support the landscape is exactly today's. The GA's Cauchy mutation and crossover, its Nelder-Mead step, and the final COBYLA run therefore behave as now inside a support; the intermittent and final local searches re-tune weights within a support and, because $f_\lambda$ drops when an entry re-enters, do not re-create an edge. Support changes come from the zeroing move (D11) during the global search and from best-first prune (D12) at the end. When `global_search = FALSE` only prune changes the support; COBYLA from a dense start does not become sparse on its own, and the record says so in the documentation text.

**(c) Numerical zero.** An entry is an edge if it is at least $10^{-5}$ after decoding, evaluated on the same thresholded matrix the shortcut receives. `normalise_sum()`'s tolerance concerns the row sum, not individual entries, and its proportional scaling maps an exact zero to an exact zero; the anchor it adjusts is the largest free entry, never a zero. Verified (appendix check 5): a parameter set to exactly 0 is exactly 0 after `param_to_solution(process = TRUE)`, after `repair_graph()`, and after `normalise_sum()`. The returned graph is therefore snapped to exact zeros in three places that already exist: `param_to_solution()` (entries below $10^{-5}$), `repair_graph()` (clamp), and `.redistribute_mass()` (`vec[drop_idx] <- 0`). No new snapping step is needed. One existing behaviour to keep in mind: `repair_graph()` fills a row whose free entries are all zero uniformly (check 5 shows $(1/3,0,1/3,1/3)$); this cannot arise from the encoding, because the derived entry equals one minus the rest, and `graph_constraint()` rejects incomplete rows whose fixed entries already sum to one (check 1).

**(d) Pinned entries.** Not counted, not removable (D9). Both counts are reported.

**(e) Ties and noise.** D10. The acceptance test in prune compares $U$ before and after a removal on the same trials; only changed trials contribute. The search inside the GA compares individuals on one fixed subsample; the price is the same number on every sample, unlike a threshold, so nothing about the cap depends on which subsample is in use.

**(f) `prune_graph()`.** D12. It is complementary: it is the one deterministic support move on the full sample, and it uses the same objective as everything before it. It is not spending the budget twice, because the budget is a fixed price and not a running total. Verified (check 7): best-first at zero price from a raw 12-edge local optimum reached the same 7-edge graph and the same $U$ as the current fixed-order prune, so switching the order when a price is set is safe on this example.

**(g) Reporting.** D14. The historical mismatch (NEWS 0.2.0) came from computing `$power` before pruning; the implementation plan keeps `final_power` as the last call in `graph_optimise()` and adds a test that `$power$trial_success` equals a fresh `calc_power_pvals()` on the returned graph with a price set.

**(h) API.** D13. Shape of the returned object: one new element `sparsity` (a list or `NULL`) appended after `start_graph`. `names(result)` therefore changes even when the feature is off; the bit-identity gate compares every pre-existing element and `.Random.seed`, and `data/graph_optimal_example.rda` is regenerated.

**(i) Error control and degenerate cases.** D15. Cases the implementer must handle:

1. A row with no removable free edge. `graph_constraint()` rejects rows with exactly one `NA` and incomplete rows whose fixed entries sum to one (check 1), so through the public constructor every incomplete row has $k_i\ge2$. Inside prune, a row can be reduced to one non-zero free entry; trying to remove it triggers the uniform fallback, which raises $E$, and the candidate is skipped. The `.redistribute_mass()` path with no recipients, which returns a row summing to less than one and makes `calc_power_pvals()` abort with "do not build a valid graph", is reachable only by bypassing validation (check 1 does so on purpose). Left as is; see adversarial check A9.
2. $P_G=0$, for example $m=2$ or a fully pinned `trans_constraint`. There is nothing to remove; `graph_optimise()` warns that `gain_tolerance` has no effect, sets the price to 0 and still returns a `sparsity` element.
3. $U_{\max}\le0$. The gain function never rewards anything; abort with an informative error before optimisation.
4. $\psi$ taking negative values. The existing invalid-encoding penalties (`-1e6` for `NA`, the sum of negative entries otherwise) assume $U\ge0$; a valid graph with $U<0$ already ranks below an encoding with a $-10^{-17}$ violation. With a price, valid graphs range down to $-cP_G=-\lambda U_{\max}$, so when a price is set the penalties are shifted by $-(1+\lambda U_{\max})$ to stay below every valid value. Behaviour when off is untouched.
5. Sparse two-cycles $g_{ij}=g_{ji}=1$ hit the shortcut's `denom == 0` branch (`src/graph_shortcut.cpp`, lines 132–148), which zeroes the row; that is the correct limit of the Bretz update and already exercised by every $m=2$ graph.
6. Unreachable nodes (zero weight, no incoming edge). Their outgoing edges do not affect $U$; at any positive price they are pruned to a single edge, which is the interpretable outcome.

**(j) Against the simplest alternative.** The simplest thing that could plausibly work is the current pipeline followed by `prune_graph()` with acceptance "loss at most $c$". It certifies the same cap, costs nothing extra, and is the last step of the recommended design. It is not enough on its own for three reasons: it violates D4; it cannot find removals that only become affordable after the weights are re-tuned, which the GA with the zeroing move and the intermittent local search can; and with the fixed index order it can spend the price on an edge that happened to come first when a cheaper one existed. The recommended design is that simple thing plus the price wired into the one objective, a best-first order, and a mutation move; there is no second pass, no reference graph, and no threshold.

## 5. Recommended design, stated for implementation

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

`expand.grid` uses no random numbers. For $m\le12$ the pattern matrix is at most 4096 rows.

### 5.2 Objective

Complete replacement for the closure factory. The decoding and penalty code is unchanged in content; the only additions are the `edge_price` argument, the free-entry mask, and the last four lines.

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

With `edge_price = 0` every returned value is arithmetically the value returned today (`x - 0` is exact in IEEE arithmetic for finite `x`), and no extra random numbers are drawn.

### 5.3 Objective on the full sample, and `choose_graph()`

`.graph_optimise_ga()` and `.graph_optimise_local()` currently compute `ga_trial_success` and `local_trial_success` on the full `pvals`. They keep those and add `ga_objective` and `local_objective`, computed as $U-c\,E$ from the processed graph (`param_to_solution(process = TRUE)` output, so $E$ uses the same zeroed entries the user will see). `choose_graph()` compares the objectives; when the price is zero they equal the trial-success values, so its decisions are unchanged.

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
            parent_vec[idx] <- parent_vec[idx] * ((1 - 5e-6) / s)
        }
        parent_vec
    }
}
```

The inner `cauchy_mutation` is the present closure verbatim, and it is what the factory returns when `p_zero` is zero, so the disabled path draws exactly the random numbers it draws today. Scaling a row to $1-5\times10^{-6}$ relies on `recover_full_trans_matrix()` computing the derived entry as `1 - sum(G[i, ])`; check 2 shows the result is $5\times10^{-6}$, not a negative number, and the closure then zeroes it. Row parameters in the encoding are bounded by `upper = 1`, and the rescaled values cannot exceed the originals when $s\ge1-5\times10^{-6}$; when $s<1-5\times10^{-6}$ the move scales up, which is intentional (it moves mass off the derived entry onto the explicit ones) and stays within bounds only if no single parameter exceeds one afterwards; the implementer clamps to `upper` and lets the objective's penalty handle the rare remainder.

### 5.5 Best-first prune

```r
# R/post_optim_processing.R: replacement for prune_edges() when a price is set.
# The existing fixed-order loop is kept unchanged and used when edge_price == 0.

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

`prune_graph()` gains `edge_price = 0`, forwards it, and dispatches to this function when it is positive; its return list gains `prune_loss` (0 when off). Acceptance `u_best - best$u <= edge_price` is the statement $f_\lambda(\text{candidate})\ge f_\lambda(\text{current})$ for a candidate with one fewer edge. Cost is at most $E^2/2$ shortcut evaluations on the full sample.

### 5.6 `graph_optimise()` and the returned object

In `graph_optimise()`: validate `gain_tolerance`; compute `edge_price <- .edge_price(gain_tolerance, trial_success, graph_constraint)`; warn once if `gain_tolerance` is non-`NULL` and the price is 0 because `n_removable == 0`; pass `edge_price` to `.graph_optimise_ga()`, `.graph_optimise_local()` and `prune_graph()`; in `.graph_optimise_ga()`, build the mutation closure with `p_zero = if (edge_price > 0) 0.2 else 0` and `param_rows = .g_param_rows(...)`; after `final_power`, assemble `sparsity` when `gain_tolerance` is non-`NULL`. `graph_optimal()` and `new_graph_optimal()` gain `sparsity = NULL`. `print.multigrain_graph_optimal()` adds, when non-`NULL`:

```
Edges: 6 (5 free of 8 removable); edge price 1.25e-04
Total loss of trial success bounded by 1e-03 (prune stage: 7.5e-05)
```

`summary()` prints the same two lines under the power block.

## 6. Rejected alternatives

1. **Threshold against $U^*$ in two stages** (the brief's literal statement): run the present pipeline, take its graph as reference, then re-run COBYLA and prune (optionally the GA) with a lexicographic objective that ranks feasibility, then $-E$, then $U$, against $T=(1-\lambda)U_{\text{ref}}$ recomputed on each sample. It certifies the cap against the dense optimum and can spend the whole budget. Rejected by the maintainer: a second pass, three thresholds to thread, a fallback rule, and roughly double runtime when the GA is repeated.
2. **Uncalibrated per-edge price.** One pass, but no total cap and a price whose meaning depends on the scale of $\psi$. Rejected by the brief and the maintainer.
3. **Single-pass price calibrated on a cheap lower bound of $U^*$**, $c=\lambda U_{\text{seed}}/P_G$. Certifies the cap against $U^*$ rather than $U_{\max}$, but is smaller still and needs a seed evaluation before the objective is built. Rejected as strictly more conservative than D3 for no practical gain.
4. **Multi-objective GA** (Pareto front over $(U,-E)$, choose afterwards). Needs a new dependency and rewrites the global search. Rejected.
5. **Counting all entries including pinned ones.** Adds a constant; changes nothing in the search and confuses the report. Rejected (D9).
6. **Entropy or $\ell_q$ sparsity.** Rejected by the brief.
7. **Lexicographic scalarisation with a large edge weight** $D>U_{\max}-U_{\min}$ so that fewer edges always win. This is a threshold design in disguise (something must stop it from removing everything) and was superseded by D3.

## 7. Implementation plan

Each step ends with a gate that must pass before the next starts. Test files are run one at a time with `testthat::test_file()`.

1. **Helpers and objective** (`R/objective_function.R`). Add `.n_removable_edges()`, `.trial_success_max()`, `.edge_price()`, and the `edge_price` argument as in 5.1 and 5.2. Gate: `test-objective_function.R` passes unchanged; new tests: for 200 random encodings at $m\in\{3,4\}$ the closure with `edge_price = 0` returns values identical (`expect_identical`) to a copy of the pre-change closure kept in the test file; with a positive price it returns $U-cE$ with $E$ computed independently from `param_to_solution(process = TRUE)`; every invalid encoding scores below $-\lambda U_{\max}$; `.n_removable_edges()` equals the parameter count of `create_start_params()` for three constraints; `.trial_success_max()` returns 1, 4, 1 for the three functions in check 6.
2. **Mutation move** (`R/mutation_helpers.R`) as in 5.4, plus `.g_param_rows()`. Gate: `test-mutation_helpers.R` passes; new tests: with `p_zero = 0` the factory returns a closure whose outputs under `set.seed()` are identical to the present closure's; with `p_zero = 1` every call either sets a free transition parameter to exactly 0 or scales a row to $1-5\times10^{-6}$ within `1e-12`; the row map matches `recover_full_trans_matrix()`'s parameter order for a constrained example.
3. **Prune** (`R/post_optim_processing.R`) as in 5.5; `prune_graph()` gains `edge_price` and returns `prune_loss`. Gate: `test-post_optim_processing.R` passes unchanged; new tests: on the 4-hypothesis fixture in that file, a price of 0 reproduces the current 7-edge result of the `gamma = 1` prune; a positive price never accepts a candidate that raises the edge count; `prune_loss <= edge_price * n_removed`; a row is never left with sum different from one.
4. **Pipeline** (`R/optimisation.R`, `R/choose_graph.R`). Thread `edge_price` through `.graph_optimise_ga()` and `.graph_optimise_local()`; add `ga_objective` and `local_objective`; `choose_graph()` compares them. Gate: `test-optimisation.R` and `test-choose_graph.R` pass with existing snapshots untouched.
5. **Object and methods** (`R/graph_optimal.R`). Add `sparsity`; print and summary lines. Gate: `test-graph_optimal.R` passes; new snapshot with a non-`NULL` `sparsity`.
6. **User-facing argument** (`R/optimisation.R`, roxygen). `gain_tolerance = NULL`; validation; warning for $P_G=0$; error for $U_{\max}\le0$; documentation of the semantics, the conservatism, the recommended $10^{-3}$ and the note about `global_search = FALSE`. Gate: end-to-end test at $m=4$, $n_{\text{sim}}=10^4$, `gain_tolerance = 5e-3`: returned graph valid, `sparsity$n_edges` at most the count from a `NULL` run with the same seed, `$power$trial_success` identical to a fresh `calc_power_pvals()` on the returned graph, and `sparsity$prune_loss <= sparsity$loss_bound`.
7. **Bit-identity gate** (section 8). Gate: passes on the branch against a fixture generated on `main`.
8. **Regenerate `data/graph_optimal_example.rda`** with `data-raw/graph_optimal_example.R` (object shape changed), run `devtools::document()`, update `NEWS.md` (section 10), update `_pkgdown.yml` only if a new exported symbol is added (none is planned). Gate: `R CMD check` clean locally with `--as-cran` off.

## 8. Test plan

Beyond the per-step gates above:

**Bit-identical when disabled.** On `main`, before any change, run the recipe below and save the result and the RNG state. On the branch, the test re-runs the recipe and compares.

```r
# tests/testthat/data/make_sparsity_baseline.R (run on main; commits the RDS)
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
test_that("graph_optimise() is bit-identical to main when gain_tolerance is NULL", {
    skip_on_cran()
    base <- readRDS(test_path("data", "sparsity_baseline.rds"))
    # same recipe as make_sparsity_baseline.R, with global_search = TRUE
    ...
    res <- graph_optimise(..., control = ctrl, verbose = "silent")
    keep <- setdiff(names(res), "sparsity")
    expect_identical(unclass(res)[keep], unclass(base$result)[keep])
    expect_identical(.Random.seed, base$seed_after)
    expect_null(res$sparsity)
})
```

The test is skipped on CRAN because nloptr and GA binaries may differ across platforms in the last bits; it is the gate for step 7 locally and in the Ubuntu CI job. The fixture records the R, GA and nloptr versions so a mismatch is diagnosable.

**Feature tests.** In addition to the per-step gates: `gain_tolerance = 0` runs, returns a non-`NULL` `sparsity`, and never lowers `$power$trial_success` below the `NULL` run's value with the same seed (best-first at zero price); `num_threads = 2` gives the same `sparsity$n_edges` as `num_threads = 1` for the same seed (check 4 verified that the parallel shortcut returns identical rejections); a constrained example with pinned non-zero edges reports `n_edges > n_edges_free`; `gain_tolerance = 1e-3` on `graph_constraint_free(2)` warns that nothing is removable.

## 9. Adversarial checks for the review agent

- **A1.** Run `graph_optimise()` on `main` and on the branch for seeds 1–5, $m\in\{2,3,4\}$, constrained and free, `global_search` on and off, `num_threads` 1 and 2, `gain_tolerance = NULL`. Any non-identical element other than `sparsity`, or any difference in `.Random.seed` afterwards, fails D5.
- **A2.** Use a gain function with negative values, such as `r1 - r2`, and a price. Confirm no invalid encoding ever outranks a valid graph in the GA's final population.
- **A3.** Constant gain function (`U_max = U_min`) and one that is never positive: confirm the error path in `.edge_price()`, not a silent price of zero or `NaN`.
- **A4.** $m=2$ and fully pinned `trans_constraint`: confirm the warning and that `sparsity$n_removable` is 0.
- **A5.** Construct a GA population where the zeroing move rescales a row whose parameters sum to more than 1 (so the derived entry is negative before the move). After the move the row sums to $1-5\times10^{-6}$; confirm no parameter exceeds `upper` and the objective is finite.
- **A6.** After `param_to_solution(process = TRUE)` snaps an entry to 0.001, confirm the price can remove it in prune when its loss is below $c$, and that it is counted in `n_edges_free` when it stays.
- **A7.** Try to make best-first prune loop for ever: every accepted step reduces the edge count by one, so at most $P_G$ iterations; verify with a `gain_tolerance = 1` run.
- **A8.** Confirm `choose_graph()` picks the GA graph when the local one is invalid and a price is set, and that the comparison uses the objective, not raw $U$.
- **A9.** The pre-existing abort when `.redistribute_mass()` has no recipients (check 1): confirm it is unreachable through `graph_constraint()` and that the best-first loop's validity skip makes it unreachable through prune with a price as well.
- **A10.** Paired-noise floor at scale: at $m=8$, $n_{\text{sim}}=10^6$, $\lambda=10^{-3}$, $c\approx2\times10^{-5}$ is close to the paired standard error for an edge that flips 0.1% of trials. Check whether prune's accept/reject at the margin flips between two independent p-value matrices; if it does, the documentation must say that $\lambda$ below $5\times10^{-3}$ at $m\ge8$ is at the noise floor. (This is the design's weakest point.)
- **A11.** `global_search = FALSE` with a price: confirm the result differs from the `NULL` run only through prune, and that the documentation says so.
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

## Bug fixes

* None in this change. The `global_opt_power` / `local_opt_power` elements
  described in the 0.2.0 news entry are not present on the returned object;
  see the open item in `dev/sparsity_design_record.md`.
```

The second bullet under bug fixes is a placeholder to be resolved by O2 before release, not shipped as written.

## 11. Open items

- **O1. Conservatism of the price.** At $\lambda=10^{-3}$ and $m\ge6$ the price removes only edges that flip a few dozen trials per million. Evidence that would close it: on the three benchmark problems the maintainers already use (not run here), the distribution of per-edge losses of the edges left after the existing prune. If most sit between $10^{-4}$ and $10^{-3}$, the documentation should recommend $5\times10^{-3}$ or the maintainers may revisit D3 with the two-stage threshold as the alternative.
- **O2. `global_opt_power` / `local_opt_power`.** Either reinstate them on the object (the values are already computed) or correct the 0.2.0 NEWS entry. Separate change; not part of the bit-identity gate.
- **O3. $p_0$ and the split between parameter-zeroing and row-rescaling.** Proposed $0.2$ and $1/2$. Evidence to close: on the $m=4$ examples, fraction of GA generations in which the elite's edge count decreases, for $p_0\in\{0.1,0.2,0.4\}$.
- **O4. Pre-existing abort in `prune_edges()` when a row has no free recipients.** Unreachable via `graph_constraint()`; fix separately if `prune_graph()` is ever exported.
- **O5. Whether `gain_tolerance = 0` should be documented as a supported mode.** It is well defined (best-first, zero price) and cheap; a user-facing name for it is a documentation decision.

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
