# Review: `graph_simplify()` on `plan-for-sparsity`

Reviewed at commit `2ee9241` against `main` at `8b19c63`.

Environment:

```text
R 4.6.0 (2026-04-24 ucrt)
multigrain 0.3.0
GA 3.2.5
NOT_CRAN=true
```

The package loaded with `pkgload::load_all()` without rebuilding. The requested
baseline:

```r
Sys.setenv(NOT_CRAN = "true")
pkgload::load_all(quiet = TRUE)
testthat::test_file(
    "tests/testthat/test-graph_simplify.R",
    reporter = "summary"
)
```

produced:

```text
graph_simplify: ....................................................................................................................................

== DONE =======================================================================
```

## Verdict

I found one PR-blocking defect, one non-blocking defect already identified in
the implementation report, one design gap, and one reproducibility gap in the
implementation report's evidence.

| # | Classification | Severity | Blocks PR | Finding |
|---|---|---:|:---:|---|
| 1 | Defect | Medium | Yes | `graph_simplify()` silently discards user-supplied `mutation` and `suggestions` in `control_global()` |
| 2 | Defect | Medium | No | GA and COBYLA do not draw random subsamples; they use permuted, overlapping prefixes |
| 3 | Design gap | Low | No | A gain is printed as a negative "loss" |
| 4 | Documentation gap | Low | No | The A1 evidence script named by the implementation report is not present or recoverable from Git |

The scorer/fallback scale, constrained edge counting, Nelder-Mead support
argument, serialisation guard for real compiled objects, example-data alpha,
and the requested remaining adversarial cases held up in the checks below.

## Findings

### 1. Explicit GA `mutation` and `suggestions` are silently discarded

**Classification:** defect  
**Severity:** medium  
**Blocks the PR:** yes

`control_global()` is documented as accepting `GA::ga()` options directly, and
the `graph_simplify()` documentation says that an explicitly supplied control
object "is used as it is". That is not true for `mutation` and `suggestions`.

In `.graph_optimise_ga()`, `global_opts` initially wins the `modifyList()`
merge. The implementation then restores the internal mutation whenever
`p_zero > 0`, and restores the internal suggestions whenever the separate
`suggestions` argument is non-`NULL`. `graph_simplify()` always passes
`p_zero = .simplify_p_zero` (`0.2`) and always passes its generated seed
matrix. Consequently, both user settings are always replaced and no condition
is signalled.

This is the mirror image of report section 3.3: the guard protects the
algorithm, but turns two documented expert controls into silently ignored
inputs. Either the inputs need to be honored, or `graph_simplify()` needs to
reject/reserve them explicitly; silently accepting them is the defect.

Reproduction:

```r
Sys.setenv(NOT_CRAN = "true")
pkgload::load_all(quiet = TRUE)

set.seed(20260907)
corr <- matrix(0.2, 4, 4)
diag(corr) <- 1
pvals <- stats::pnorm(
    mvtnorm::rmvnorm(
        1e4,
        mean = calc_ncp(c(0.93, 0.91, 0.89, 0.86)),
        sigma = corr
    ),
    lower.tail = FALSE
)
ts <- trial_success(
    0.25 * (2 * (r1 && r2) + r1 * r3 + r2 * r4),
    verbose = "silent"
)
gc <- graph_constraint_free(4)
w <- rep(0.25, 4)
G <- matrix(1 / 3, 4, 4)
diag(G) <- 0
ctrl <- multigrain_control() |>
    control_global(run = 2, maxiter = 2, popSize = 12) |>
    control_local(maxeval = 20)
ref <- graph_optimal(
    w, G,
    constraints = gc,
    trial_success = ts,
    power = calc_power_pvals(
        pvals, w, G,
        custom_power = list(trial_success = ts)
    ),
    global_search = FALSE,
    control = ctrl,
    alpha = 0.025
)

x_ref <- .encode_graph(gc, w, G)
user_mutation <- function(object, parent) rep(0.123, length(x_ref))
user_suggestions <- matrix(0.234, nrow = 1, ncol = length(x_ref))
user_control <- ctrl |>
    control_global(
        mutation = user_mutation,
        suggestions = user_suggestions
    )

captured <- NULL
trace(
    "ga",
    where = asNamespace("GA"),
    tracer = quote({
        assign(
            "captured",
            list(mutation = mutation, suggestions = suggestions),
            envir = .GlobalEnv
        )
        stop("captured")
    }),
    print = FALSE
)
warnings <- character()
tryCatch(
    withCallingHandlers(
        graph_simplify(
            ref, pvals,
            gain_tolerance = 0.01,
            control = user_control,
            verbose = "silent"
        ),
        warning = function(cnd) {
            warnings <<- c(warnings, conditionMessage(cnd))
            invokeRestart("muffleWarning")
        }
    ),
    error = function(cnd) {
        stopifnot(grepl("captured", conditionMessage(cnd), fixed = TRUE))
    }
)
untrace("ga", where = asNamespace("GA"))

cat("kept mutation:", identical(captured$mutation, user_mutation), "\n")
cat(
    "kept suggestions:",
    identical(captured$suggestions, user_suggestions),
    "\n"
)
cat("warning:", length(warnings) > 0, "\n")
```

Output:

```text
kept mutation: FALSE
kept suggestions: FALSE
warning: FALSE
```

The guard conditions themselves are precise. Direct calls to
`.graph_optimise_ga()` gave:

```text
neither internal feature active: kept mutation TRUE, kept suggestions TRUE
only p_zero active:              kept mutation FALSE, kept suggestions TRUE
only internal seeds active:      kept mutation TRUE, kept suggestions FALSE
both active:                     kept mutation FALSE, kept suggestions FALSE
```

Thus ordinary `graph_optimise()` keeps these user options because it passes
`p_zero = 0` and no separate suggestions; `graph_simplify()` discards both.

### 2. The two optimiser samples are correlated prefixes, not independent random samples

**Classification:** defect  
**Severity:** medium  
**Blocks the PR:** no

Both optimiser functions use:

```r
pvals_sampled <- pvals[sample(nsim), ]
```

`sample(nsim)` permutes `1:nsim`; it never selects a row above `nsim`. When
`nsim_global == nsim_local`, GA and COBYLA therefore use the same rows in
different orders. When the sizes differ, the smaller sample is a subset of the
larger prefix. Their sampling noise is not independent.

Reproduction:

```r
n_total <- 10000L
n_equal <- 2000L
set.seed(81); ga_idx <- sample(n_equal)
set.seed(82); local_idx <- sample(n_equal)
set.seed(83); small_idx <- sample(1000L)
set.seed(84); large_idx <- sample(2000L)

cat(max(ga_idx), max(local_idx), "\n")
cat(identical(sort(ga_idx), sort(local_idx)), "\n")
cat(identical(ga_idx, local_idx), "\n")
cat(all(small_idx %in% large_idx), "\n")
```

Output:

```text
2000 2000
TRUE
FALSE
TRUE
```

Order can still change the last bits of a floating-point reduction. With seed
`2026`, 10,000 four-hypothesis trials, and
`trial_success(0.1*r1 + 0.2*r2 + 0.3*r3 + 0.4*r4)`, two permutations of the
same first 2,000 rows produced:

```text
equal-row-set thresholds identical: FALSE
equal-row-set threshold absolute difference: 4.4408920985006262e-16
```

The design record's statement that each closure computes its own threshold
remains mechanically true. Its statistical implication does not: equal-sized
closures have the same empirical sample and almost the same threshold, and
unequal-sized closures use nested samples. The paired reference-versus-
candidate comparison within each closure is still sound, and the final
full-sample fallback still enforces the public cap. What is lost is the claimed
independent sampling noise between GA and COBYLA. This can affect search
diversity and exact boundary decisions, but did not invalidate any returned
graph in the checks below.

This behavior predates this branch, which is why I do not treat it as a blocker
for this PR, but it is contrary to the record and to the `control_nsim_*`
documentation's "random sample" wording.

### 3. Improvements print as a negative "loss"

**Classification:** design gap  
**Severity:** low  
**Blocks the PR:** no

In the reproduced case, pruning improves on the reference itself, so
`prune_loss`, `gain_loss`, and `gain_loss_fraction` are all negative and
coherent. More generally `prune_loss` is measured from the graph entering
pruning while `gain_loss` is measured from the reference, so their signs need
not match. The print and summary methods use `gain_loss_fraction` but always
label it as a "loss", producing text such as `loss -1.59%`.
That follows the record's prescribed signed field and literal output template;
the gap is that the record did not specify gain-aware wording when the signed
"loss" is negative.

Reproduction, seed `5`:

```r
set.seed(5)
corr <- matrix(0.2, 4, 4)
diag(corr) <- 1
pvals <- stats::pnorm(
    mvtnorm::rmvnorm(
        1e4,
        mean = calc_ncp(c(0.93, 0.91, 0.90, 0.85)),
        sigma = corr
    ),
    lower.tail = FALSE
)
ts <- trial_success(r1 && r2 && r3 && r4, verbose = "silent")
w <- c(0.5, 0.5, 0, 0)
G <- rbind(
    c(0, 0.641, 0.335, 0.024),
    c(0.283, 0, 0.335, 0.382),
    c(0.597, 0.336, 0, 0.067),
    c(0.337, 0.585, 0.078, 0)
)
gc <- graph_constraint_free(4)
u_ref <- calc_power_pvals(
    pvals, w, G, custom_power = ts
)$custom_power
pruned <- .prune_edges_best_first(
    pvals, w, G, ts,
    fixed_edge = !is.na(gc$trans_constraint),
    threshold = u_ref
)

sp <- list(
    gain_reference = u_ref,
    gain = pruned$power_best,
    gain_loss = u_ref - pruned$power_best,
    gain_loss_fraction = (u_ref - pruned$power_best) / u_ref,
    gain_tolerance = 0,
    n_edges_reference = sum(G != 0),
    n_edges = sum(pruned$trans_matrix != 0),
    n_edges_free_reference = sum(G != 0),
    n_edges_free = sum(pruned$trans_matrix != 0),
    prune_loss = pruned$prune_loss
)
summarise_sparsity(sp)
cat("prune_loss:", sp$prune_loss, "\n")
cat("gain_loss:", sp$gain_loss, "\n")
cat("gain_loss_fraction:", sp$gain_loss_fraction, "\n")
```

Output:

```text
Simplified from 12 edges to 5 (free: 12 -> 5)
Trial success 0.6605 -> 0.6710: loss -1.59% of reference (cap 0.0%)
prune_loss: -0.0105
gain_loss: -0.0105
gain_loss_fraction: -0.0158970476911
```

`summary()` prints the same sparsity line. In an actual reference fallback with
no accepted pruning, all three values are zero and the line is coherent:

```text
fallback source: reference
fallback prune_loss: 0
fallback gain_loss: 0
fallback gain_loss_fraction: 0
Trial success 0.8048 -> 0.8048: loss 0.00% of reference (cap 0.0%)
```

### 4. The A1 evidence script is unavailable in a fresh clone

**Classification:** documentation gap  
**Severity:** low  
**Blocks the PR:** no

The implementation report cites `scratchpad/a1_unchanged.R`, but the file is
not in the working tree, ignored files, or any Git ref. Its stated Cartesian
grid contains `5 * 3 * 2 * 2 * 2 = 120` cases, while the narrative and output
report 100. The sweep therefore cannot be inspected, rerun, or reconciled from
the report alone.

Reproduction:

```powershell
"working tree exists: $(Test-Path -LiteralPath 'scratchpad\a1_unchanged.R')"
$history = git log --all --format='%h' -- scratchpad/a1_unchanged.R
"commits containing path: $(@($history).Count)"
```

Output:

```text
working tree exists: False
commits containing path: 0
```

I did not recreate the claimed sweep. The smaller independent spot-check below
held.

## Checked and held up

### Full-sample scorer, optimiser scores, and fallback use the same scale

`calc_power_pvals()` and `.make_lexico_scorer()` both use serial
`graph_shortcut()` for their full-sample reference calculation.
`num_threads` affects the optimiser's subsample objective, not these
full-sample comparisons.

The following check used seed `20260907`, 10,000 rows, the reference graph, and
all 12 single-edge-removal neighbors:

```r
set.seed(20260907)
corr <- matrix(0.2, 4, 4)
diag(corr) <- 1
pvals <- stats::pnorm(
    mvtnorm::rmvnorm(
        1e4,
        mean = calc_ncp(c(0.93, 0.91, 0.89, 0.86)),
        sigma = corr
    ),
    lower.tail = FALSE
)
ts <- trial_success(
    0.25 * (2 * (r1 && r2) + r1 * r3 + r2 * r4),
    verbose = "silent"
)
gc <- graph_constraint_free(4)
w_ref <- rep(0.25, 4)
G_ref <- matrix(1 / 3, 4, 4)
diag(G_ref) <- 0

u_calc <- calc_power_pvals(
    pvals, w_ref, G_ref,
    alpha = 0.025,
    custom_power = list(trial_success = ts)
)$trial_success
rej_serial <- graph_shortcut(pvals, 0.025, w_ref, G_ref)
rej_parallel <- graph_shortcut_parallel(
    pvals, 0.025, w_ref, G_ref,
    num_threads = 2L,
    grain_size = -1L
)
u_short <- ts$func(rej_serial)
u_parallel <- ts$func(rej_parallel)

args <- list(
    gain_tolerance = 0.01,
    ref_graph = list(hyp_weight = w_ref, trans_matrix = G_ref),
    u_range = .trial_success_range(ts)
)
scorer <- .make_lexico_scorer(
    pvals, 0.025, ts, gc$trans_constraint, args
)
serial_objective <- do.call(
    create_obj_func,
    c(list(
        m = 4,
        power_criterion = ts$func,
        hyp_constraint = gc$hyp_constraint,
        trans_constraint = gc$trans_constraint,
        pvals = pvals,
        alpha = 0.025,
        num_threads = 1L
    ), args)
)
parallel_objective <- do.call(
    create_obj_func,
    c(list(
        m = 4,
        power_criterion = ts$func,
        hyp_constraint = gc$hyp_constraint,
        trans_constraint = gc$trans_constraint,
        pvals = pvals,
        alpha = 0.025,
        num_threads = 2L
    ), args)
)

free <- is.na(gc$trans_constraint)
threshold <- 0.99 * u_calc
edge_price <- (.trial_success_range(ts)[["max"]] - threshold) + 1
fallback_score <- function(G) {
    u <- calc_power_pvals(
        pvals, w_ref, G,
        alpha = 0.025,
        custom_power = ts
    )$custom_power
    .lexico(u, sum(G[free] != 0), threshold, edge_price, sum(free))
}

candidates <- list(reference = G_ref)
for (k in seq_len(nrow(which(G_ref != 0, arr.ind = TRUE)))) {
    ij <- which(G_ref != 0, arr.ind = TRUE)[k, ]
    G_try <- G_ref
    G_try[ij[1], ] <- .redistribute_mass(
        G_ref[ij[1], ],
        drop_idx = ij[2],
        fixed_idx = which(!is.na(gc$trans_constraint[ij[1], ]))
    )
    candidates[[paste(ij, collapse = "_")]] <- G_try
}

routes <- lapply(candidates, function(G) {
    u <- ts$func(graph_shortcut(pvals, 0.025, w_ref, G))
    x <- .encode_graph(gc, w_ref, G)
    c(
        scorer = scorer(u, G),
        fallback = fallback_score(G),
        serial = serial_objective(x),
        parallel = parallel_objective(x)
    )
})
route_matrix <- do.call(rbind, routes)
route_identical <- function(a, b) {
    all(vapply(seq_len(nrow(route_matrix)), function(i) {
        identical(
            unname(route_matrix[i, a]),
            unname(route_matrix[i, b])
        )
    }, logical(1)))
}

cat(
    "serial rejection matrix identical to parallel:",
    identical(rej_serial, rej_parallel), "\n"
)
cat(
    "u_ref calc_power identical to direct shortcut:",
    identical(u_calc, u_short), "\n"
)
cat(
    "u_ref direct serial identical to parallel:",
    identical(u_short, u_parallel), "\n"
)
cat(
    "scorer identical to fallback for all candidates:",
    route_identical("scorer", "fallback"), "\n"
)
cat(
    "scorer identical to serial objective for all candidates:",
    route_identical("scorer", "serial"), "\n"
)
cat(
    "serial objective identical to parallel objective for all candidates:",
    route_identical("serial", "parallel"), "\n"
)
cat(
    "max absolute scorer/fallback difference:",
    max(abs(route_matrix[, "scorer"] - route_matrix[, "fallback"])),
    "\n"
)
cat(
    "candidate ranking identical:",
    identical(
        order(route_matrix[, "scorer"], decreasing = TRUE),
        order(route_matrix[, "fallback"], decreasing = TRUE)
    ),
    "\n"
)

losses <- vapply(candidates[-1], function(G) {
    u_calc - ts$func(graph_shortcut(pvals, 0.025, w_ref, G))
}, numeric(1))
boundary_name <- names(losses)[which(losses > 0)[1]]
boundary_G <- candidates[[boundary_name]]
boundary_u <- ts$func(graph_shortcut(
    pvals, 0.025, w_ref, boundary_G
))
boundary_lambda <- 1 - boundary_u / u_calc
boundary_args <- list(
    gain_tolerance = boundary_lambda,
    ref_graph = list(hyp_weight = w_ref, trans_matrix = G_ref),
    u_range = .trial_success_range(ts)
)
boundary_scorer <- .make_lexico_scorer(
    pvals, 0.025, ts, gc$trans_constraint, boundary_args
)
boundary_threshold <- (1 - boundary_lambda) * u_calc
boundary_price <- (
    .trial_success_range(ts)[["max"]] - boundary_threshold
) + 1
boundary_fallback <- .lexico(
    boundary_u,
    sum(boundary_G[free] != 0),
    boundary_threshold,
    boundary_price,
    sum(free)
)
cat("boundary candidate:", boundary_name, "\n")
cat("boundary lambda:", format(boundary_lambda, digits = 17), "\n")
cat(
    "candidate u equals reconstructed threshold:",
    identical(boundary_u, boundary_threshold), "\n"
)
cat(
    "boundary scorer identical to fallback:",
    identical(
        boundary_scorer(boundary_u, boundary_G),
        boundary_fallback
    ),
    "\n"
)
```

Output:

```text
serial rejection matrix identical to parallel: TRUE
u_ref calc_power identical to direct shortcut: TRUE
u_ref direct serial identical to parallel: TRUE
scorer identical to fallback for all candidates: TRUE
scorer identical to serial objective for all candidates: TRUE
serial objective identical to parallel objective for all candidates: TRUE
max absolute scorer/fallback difference: 0
candidate ranking identical: TRUE
boundary candidate: 2_1
boundary lambda: 0.0016210031148687909
candidate u equals reconstructed threshold: TRUE
boundary scorer identical to fallback: TRUE
```

I found no route by which the independently constructed scorer changes the
graph selected by `choose_graph()` or by the final fallback.

### Total-edge filters are equivalent to free-edge filters here

Pinned entries are excluded from redistribution and remain unchanged. Their
count is therefore a constant offset:

```text
total_edges(G) = free_edges(G) + pinned_nonzero_edges
```

so the sign of the before/after edge-count change is identical for total and
free counts.

The concrete constrained graph below includes two pinned non-zero entries. One
removal activates more free recipients and increases both counts; another
removal decreases both:

```r
tc <- matrix(
    c(
        0, NA, NA, NA,
        0, 0, 1, 0,
        NA, NA, 0, NA,
        1, 0, 0, 0
    ),
    nrow = 4,
    byrow = TRUE
)
fixed <- !is.na(tc)
G <- matrix(
    c(
        0, 1, 0, 0,
        0, 0, 1, 0,
        0.5, 0, 0, 0.5,
        1, 0, 0, 0
    ),
    nrow = 4,
    byrow = TRUE
)
count_edges <- function(x) {
    c(total = sum(x != 0), free = sum(x[!fixed] != 0))
}
drop <- function(G, i, j) {
    out <- G
    out[i, ] <- .redistribute_mass(
        G[i, ],
        drop_idx = j,
        fixed_idx = which(fixed[i, ])
    )
    out
}
count_edges(G)
count_edges(drop(G, 1, 2))
count_edges(drop(G, 3, 1))

supports <- Filter(
    length,
    lapply(0:7, function(mask) which(as.logical(intToBits(mask)[1:3])))
)
deltas_match <- fixed_match <- logical()
for (s1 in supports) {
    for (s3 in supports) {
        G_case <- G
        G_case[1, c(2, 3, 4)] <- 0
        G_case[1, c(2, 3, 4)[s1]] <- 1 / length(s1)
        G_case[3, c(1, 2, 4)] <- 0
        G_case[3, c(1, 2, 4)[s3]] <- 1 / length(s3)
        cand <- which(G_case != 0 & !fixed, arr.ind = TRUE)
        for (k in seq_len(nrow(cand))) {
            G_try <- drop(G_case, cand[k, 1], cand[k, 2])
            delta <- count_edges(G_try) - count_edges(G_case)
            deltas_match <- c(deltas_match, delta[1] == delta[2])
            fixed_match <- c(
                fixed_match,
                identical(G_try[fixed], G_case[fixed])
            )
        }
    }
}
stopifnot(all(deltas_match), all(fixed_match))
```

Output:

```text
pinned non-zero entries: 2
activation case before: 5/3 after: 6/4
reduction case before: 5/3 after: 4/2
```

I also enumerated all non-empty support combinations for the two free rows,
168 candidate removals in total:

```text
fixed entries unchanged in every comparison: TRUE
total and free edge deltas equal in every comparison: TRUE
filter decisions equal in every comparison: TRUE
```

The filters correctly skip a "removal" when redistribution creates at least as
many free edges as it removes. I found no constrained case where counting total
entries instead of free entries changes the decision.

### A14: Nelder-Mead did not recreate edges in feasible elites

I forced `poptim = 1` and `pressel = 1` so each GA generation selected the
current elite for its `stats::optim()` call. A temporary trace recorded the
input and returned parameter vectors. Each vector was decoded with the same
thresholding and objective closure used by the GA.

Core instrumentation:

```r
set.seed(314159)
corr <- matrix(0.2, 4, 4)
diag(corr) <- 1
pvals <- stats::pnorm(
    mvtnorm::rmvnorm(
        4000,
        mean = calc_ncp(c(0.93, 0.91, 0.89, 0.86)),
        sigma = corr
    ),
    lower.tail = FALSE
)
ts <- trial_success(
    0.25 * (2 * (r1 && r2) + r1 * r3 + r2 * r4),
    verbose = "silent"
)
gc <- graph_constraint_free(4)
w_ref <- rep(0.25, 4)
G_ref <- matrix(1 / 3, 4, 4)
diag(G_ref) <- 0
objective_args <- list(
    gain_tolerance = 0.01,
    ref_graph = list(hyp_weight = w_ref, trans_matrix = G_ref),
    u_range = .trial_success_range(ts)
)
seeds <- .build_simplify_seeds(
    gc,
    objective_args$ref_graph,
    pop_size = 30L
)
global_opts <- default_control()$global_opt
global_opts$popSize <- 30L
global_opts$maxiter <- 8L
global_opts$run <- 8L
global_opts$monitor <- FALSE
global_opts$optim <- TRUE

diagnose_par <- function(par, fn) {
    theta <- split_theta(par, gc$hyp_constraint)
    w <- recover_full_weights(theta$w_pars, gc$hyp_constraint)
    G <- recover_full_trans_matrix(theta$g_pars, gc$trans_constraint)
    if (
        anyNA(w) || anyNA(G) ||
            any(w < 0 | w > 1) ||
            any(G < 0 | G > 1)
    ) {
        return(list(
            feasible = FALSE,
            edges = NA_integer_,
            support = rep(NA, sum(is.na(gc$trans_constraint))),
            fitness = fn(par)
        ))
    }
    w[w < 1e-4] <- 0
    G[G < 1e-5] <- 0
    env <- environment(fn)
    u <- env$power_criterion(env$shortcut(w, G))
    list(
        feasible = u >= env$threshold,
        edges = sum(G[is.na(gc$trans_constraint)] != 0),
        support = G[is.na(gc$trans_constraint)] != 0,
        fitness = fn(par)
    )
}

.review_optim_entries <- list()
.review_optim_logs <- list()
trace(
    "optim",
    where = asNamespace("stats"),
    tracer = quote({
        frames <- sys.frames()
        ga_frames <- which(vapply(frames, function(frame) {
            exists("Fitness", envir = frame, inherits = FALSE) &&
                exists("Pop", envir = frame, inherits = FALSE) &&
                exists("i", envir = frame, inherits = FALSE)
        }, logical(1)))
        ga_frame <- frames[[tail(ga_frames, 1L)]]
        fitness_values <- get("Fitness", envir = ga_frame)
        selected <- get("i", envir = ga_frame)
        entries <- get(".review_optim_entries", envir = .GlobalEnv)
        entries[[length(entries) + 1L]] <- list(
            par = par,
            fn = fn,
            selected_is_elite = selected %in%
                which(fitness_values == max(fitness_values, na.rm = TRUE))
        )
        assign(".review_optim_entries", entries, envir = .GlobalEnv)
    }),
    exit = quote({
        entries <- get(".review_optim_entries", envir = .GlobalEnv)
        entry <- entries[[length(entries)]]
        logs <- get(".review_optim_logs", envir = .GlobalEnv)
        logs[[length(logs) + 1L]] <- list(
            before = entry,
            after_par = returnValue()$par,
            after_fitness = entry$fn(returnValue()$par)
        )
        assign(".review_optim_logs", logs, envir = .GlobalEnv)
    }),
    print = FALSE
)

global_opts$optimArgs <- list(
    method = "Nelder-Mead",
    poptim = 1,
    pressel = 1,
    control = list(fnscale = -1, maxit = c(20L, 20L))
)
set.seed(271828)
ga_result <- .graph_optimise_ga(
    pvals = pvals,
    graph_constraint = gc,
    trial_success = ts,
    nsim = 2000L,
    global_opts = global_opts,
    alpha = 0.025,
    num_threads = 1L,
    verbose = "silent",
    suggestions = seeds,
    p_zero = 0.2,
    objective_args = objective_args
)
untrace("optim", where = asNamespace("stats"))

rows <- lapply(.review_optim_logs, function(item) {
    before <- diagnose_par(item$before$par, item$before$fn)
    after <- diagnose_par(item$after_par, item$before$fn)
    c(
        elite = item$before$selected_is_elite,
        before_feasible = before$feasible,
        after_feasible = after$feasible,
        before_edges = before$edges,
        after_edges = after$edges,
        support_changed = !identical(before$support, after$support),
        fitness_decreased = after$fitness < before$fitness
    )
})
do.call(rbind, rows)
```

The run used seed `314159` for the 4,000-row p-value fixture and seed `271828`
for the GA:

```text
optim steps logged: 9
all selected individuals were elite: TRUE
feasible-to-feasible steps: 9
free-edge increases among feasible-to-feasible steps: 0
support changes among feasible-to-feasible steps: 0
fitness decreases: 0

before_edges -> after_edges
6 -> 6
6 -> 6
6 -> 6
6 -> 6
6 -> 6
6 -> 6
6 -> 6
5 -> 5
5 -> 5
```

This supports the report's argument for the tested run. It is evidence about
the observed Nelder-Mead steps, not a proof for every objective or optimiser
configuration.

### The serialisation guard catches real saved objects and gives a working command

Three valid live objectives were accepted, and the same three compiled
objects were rejected after `saveRDS()`/`readRDS()`:

```r
objectives <- list(
    quote(r1 + r2 + r3 + r4),
    quote(r1 && r2 && r3 && r4),
    quote(0.25 * (2 * (r1 && r2) + r1 * r3 + r2 * r4))
)
live <- vapply(objectives, function(expr) {
    x <- rlang::inject(trial_success(!!expr, verbose = "silent"))
    .trial_success_is_live(x)
}, logical(1))
dead <- vapply(objectives, function(expr) {
    x <- rlang::inject(trial_success(!!expr, verbose = "silent"))
    path <- tempfile(fileext = ".rds")
    saveRDS(x, path)
    on.exit(unlink(path), add = TRUE)
    .trial_success_is_live(readRDS(path))
}, logical(1))
cat(live, "\n")
cat(dead, "\n")
```

Output:

```text
TRUE TRUE TRUE
FALSE FALSE FALSE
```

For the non-trivial expression, the early error included:

```text
graph_optimal$trial_success <-
trial_success(0.25 * (2 * (r1 && r2) + r1 * r3 + r2 * r4))
```

Evaluating that expression rebuilt a live measure, after which
`graph_simplify()` completed.

The round-trip and rebuild check was:

```r
# `ref` is the four-hypothesis object constructed in finding 1.
path <- tempfile(fileext = ".rds")
saveRDS(ref, path)
dead_ref <- readRDS(path)
tryCatch(
    graph_simplify(
        dead_ref, pvals,
        global_search = FALSE,
        verbose = "silent"
    ),
    error = function(cnd) cat(conditionMessage(cnd), "\n")
)
rebuild <- paste0(
    "trial_success(",
    dead_ref$trial_success$objective,
    ", verbose = \"silent\")"
)
dead_ref$trial_success <- eval(parse(text = rebuild))
stopifnot(.trial_success_is_live(dead_ref$trial_success))
stopifnot(inherits(
    graph_simplify(
        dead_ref, pvals,
        gain_tolerance = 0,
        global_search = FALSE,
        control = ctrl,
        verbose = "silent"
    ),
    "multigrain_graph_optimal"
))
unlink(path)
```

The helper is intentionally only a liveness probe. A synthetic object whose
function succeeds for the one-row all-`TRUE` probe but fails for other matrix
shapes passes the guard and then exposes its own error:

```r
synthetic <- ts
synthetic$func <- function(x) {
    if (nrow(x) == 1L && all(x)) return(1)
    stop("synthetic late failure")
}
.trial_success_is_live(synthetic)
```

```text
TRUE
shape-sensitive failure after probe: synthetic late failure
```

I do not classify that as a defect: it requires replacing the compiled
function inside an otherwise valid-looking object. No genuinely serialised
compiled function slipped through, and no genuine live object was rejected.

### The patched example data uses the generation script's alpha

The generation call in `data-raw/graph_optimal_example.R` does not pass
`alpha`; the `graph_optimise()` default at the commit that introduced the
dataset was `0.025`. The data-raw script has not changed since that commit, and
the patched object stores `0.025`.

Reproduction:

```r
data("graph_optimal_example", package = "multigrain")
cat("stored alpha:", graph_optimal_example$alpha, "\n")

source_lines <- readLines("data-raw/graph_optimal_example.R")
call_text <- paste(source_lines[45:54], collapse = "\n")
cat("generation call supplies alpha:", grepl("alpha\\s*=", call_text), "\n")

old <- system2(
    "git",
    c("show", "1c028a7:R/optimisation.R"),
    stdout = TRUE
)
i <- grep("^graph_optimise <- function", old)
grep("alpha\\s*=", old[i:(i + 20)], value = TRUE)

system2(
    "git",
    c(
        "diff", "--quiet", "1c028a7", "HEAD", "--",
        "data-raw/graph_optimal_example.R"
    )
)
```

Output:

```text
stored alpha: 0.025
generation call supplies alpha: FALSE
alpha = 0.025,
data-raw script unchanged since dataset introduction: TRUE
```

This verifies the patch without regenerating the `2^20`-row dataset.

### Shared fixture for the remaining adversarial checks

The following fixture is the one used for A2, A7, A9, A15, and the requested
degenerate cases:

```r
set.seed(2)
pvals <- simulate_pvalues(
    c(0.93, 0.91, 0.90, 0.85),
    corr_matrix = matrix(0.2, 4, 4) + diag(0.8, 4),
    nsim = 1e4
)
ts <- trial_success(
    0.25 * (2 * (r1 && r2) + r1 * r3 + r2 * r4),
    verbose = "silent"
)
gc <- graph_constraint_free(4)
U_of <- function(w, G, pv = pvals) {
    calc_power_pvals(
        pv, w, G,
        custom_power = ts
    )$custom_power
}

stage1_control <- multigrain_control() |>
    control_local(maxeval = 1000, print_level = 0)
prepared <- control_prepare(stage1_control, pvals, verbose = "silent")
set.seed(5)
local_ref <- .graph_optimise_local(
    pvals, gc, ts,
    local_opts = prepared$local_opt,
    nsim = nrow(pvals),
    verbose = "silent"
)
pruned_ref <- prune_graph(
    pvals,
    local_ref$local_hyp_weight,
    local_ref$local_trans_matrix,
    ts,
    gc,
    gamma = 1,
    verbose = "silent"
)
w_ref <- pruned_ref$hyp_weight
G_ref <- pruned_ref$trans_matrix
u_ref <- U_of(w_ref, G_ref)
ref_object <- graph_optimal(
    w_ref, G_ref,
    constraints = gc,
    trial_success = ts,
    power = calc_power_pvals(
        pvals, w_ref, G_ref,
        custom_power = list(trial_success = ts)
    ),
    global_search = FALSE,
    control = stage1_control,
    alpha = 0.025
)
small_control <- multigrain_control() |>
    control_global(run = 3, maxiter = 5, popSize = 20) |>
    control_local(maxeval = 250, print_level = 0)
```

### A2: a subsample-boundary candidate did not escape the full-sample cap

Using the design-record fixture (p-values seed `2`, stage-1 local search seed
`5`), I enumerated single-edge removals and selected one whose gain was exactly
at its 2,000-row prefix threshold but below the corresponding full-sample
threshold. I then ran `graph_simplify()` with seed `7001` and both optimiser
sample sizes set to 2,000.

The selection and execution code was:

```r
row_score <- function(rej) {
    0.25 * (
        2 * (rej[, 1] & rej[, 2]) +
            rej[, 1] * rej[, 3] +
            rej[, 2] * rej[, 4]
    )
}
ref_sub <- row_score(graph_shortcut(
    pvals[1:2000, ], 0.025, w_ref, G_ref
))
ref_full <- row_score(graph_shortcut(pvals, 0.025, w_ref, G_ref))

G_candidate <- G_ref
G_candidate[1, ] <- .redistribute_mass(
    G_ref[1, ],
    drop_idx = 3,
    fixed_idx = which(!is.na(gc$trans_constraint[1, ]))
)
candidate_sub <- row_score(graph_shortcut(
    pvals[1:2000, ], 0.025, w_ref, G_candidate
))
candidate_full <- row_score(graph_shortcut(
    pvals, 0.025, w_ref, G_candidate
))

lambda <- 1 - mean(candidate_sub) / mean(ref_sub)
subsample_margin <- mean(candidate_sub) - (1 - lambda) * mean(ref_sub)
paired_se <- stats::sd(candidate_sub - ref_sub) / sqrt(2000)
full_margin <- mean(candidate_full) - (1 - lambda) * mean(ref_full)

boundary_control <- small_control |>
    control_nsim_global(2000) |>
    control_nsim_local(2000)
set.seed(7001)
result <- graph_simplify(
    ref_object,
    pvals,
    gain_tolerance = lambda,
    control = boundary_control,
    verbose = "silent"
)
stopifnot(
    result$sparsity$gain >=
        (1 - lambda) * result$sparsity$gain_reference
)
```

Output:

```text
boundary removal: 1 -> 3
subsample rows: 2000
lambda: 0.010429638854296375
subsample margin: 0
paired standard error: 0.00246522683727
candidate full-sample margin: -0.00135596590909
returned source: local
returned full-sample cap holds: TRUE
```

The tested neighbor would have violated the full-sample cap, but it was not
returned. The selected local result passed the full-sample comparison.

### A3: oversized seed sets are truncated with the reference first

Using a dense four-hypothesis reference produced 14 unique seeds. With
`popSize = 10`, tracing the arguments at `GA::ga()` showed:

```text
available unique seeds before truncation: 14
configured popSize: 10
suggestions passed to GA: 10
reference is first suggestion: TRUE
GA completed: TRUE
```

The core assertions were:

```r
dense_w <- rep(0.25, 4)
dense_G <- matrix(1 / 3, 4, 4)
diag(dense_G) <- 0
dense_ref <- list(hyp_weight = dense_w, trans_matrix = dense_G)
all_seeds <- .build_simplify_seeds(gc, dense_ref, pop_size = 1000L)
passed <- .build_simplify_seeds(gc, dense_ref, pop_size = 10L)
stopifnot(nrow(all_seeds) > 10L)
stopifnot(nrow(passed) == 10L)
stopifnot(identical(
    unname(passed[1, ]),
    unname(.encode_graph(gc, dense_w, dense_G))
))

dense_object <- graph_optimal(
    dense_w, dense_G,
    constraints = gc,
    trial_success = ts,
    power = calc_power_pvals(
        pvals, dense_w, dense_G,
        custom_power = list(trial_success = ts)
    ),
    global_search = FALSE,
    control = stage1_control,
    alpha = 0.025
)
captured_suggestions <- NULL
trace(
    "ga",
    where = asNamespace("GA"),
    tracer = quote(assign(
        "captured_suggestions",
        suggestions,
        envir = .GlobalEnv
    )),
    print = FALSE
)
set.seed(7002)
result <- graph_simplify(
    dense_object,
    pvals,
    gain_tolerance = 0.01,
    control = multigrain_control() |>
        control_global(run = 2, maxiter = 2, popSize = 10) |>
        control_local(maxeval = 50, print_level = 0),
    verbose = "silent"
)
untrace("ga", where = asNamespace("GA"))
stopifnot(nrow(captured_suggestions) == 10L)
stopifnot(identical(
    unname(captured_suggestions[1, ]),
    unname(.encode_graph(gc, dense_w, dense_G))
))
stopifnot(is_ga(result$global_output))
```

### A7: the reference and cap are recomputed on supplied p-values

Fresh same-design p-values used seed `99`; different-design p-values used seed
`100`. Simplification seeds were `7003` and `7004`.

```r
set.seed(99)
same_design <- simulate_pvalues(
    c(0.93, 0.91, 0.90, 0.85),
    corr_matrix = matrix(0.2, 4, 4) + diag(0.8, 4),
    nsim = 1e4
)
set.seed(100)
different_design <- simulate_pvalues(
    c(0.65, 0.95, 0.75, 0.55),
    corr_matrix = diag(4),
    nsim = 1e4
)

samples <- list(
    same_design = same_design,
    different_design = different_design
)
seeds <- c(same_design = 7003, different_design = 7004)
for (name in names(samples)) {
    pv <- samples[[name]]
    set.seed(seeds[[name]])
    result <- graph_simplify(
        ref_object,
        pv,
        gain_tolerance = 0.01,
        global_search = FALSE,
        control = small_control,
        verbose = "silent"
    )
    direct <- calc_power_pvals(
        pv, w_ref, G_ref, custom_power = ts
    )$custom_power
    stopifnot(identical(result$sparsity$gain_reference, direct))
    stopifnot(result$sparsity$gain >= 0.99 * direct)
}
```

Output:

```text
same_design stored stage1 gain differs: TRUE recomputed exactly: TRUE cap holds: TRUE
different_design stored stage1 gain differs: TRUE recomputed exactly: TRUE cap holds: TRUE
```

### A9: COBYLA recovered from an infeasible sparse start

I reproduced the design record's construction: prune the reference at
`lambda = 0.02`, then use that five-edge graph as the start for a
`lambda = 0.01` local objective. The p-values used seed `2`, the stage-1 local
search seed `5`, and this COBYLA call seed `11`.

```r
loose <- .prune_edges_best_first(
    pvals, w_ref, G_ref, ts,
    fixed_edge = !is.na(gc$trans_constraint),
    threshold = 0.98 * u_ref
)
x_sparse <- .encode_graph(gc, w_ref, loose$trans_matrix)
set.seed(11)
local <- .graph_optimise_local(
    pvals, gc, ts,
    local_opts = control_prepare(
        multigrain_control() |>
            control_local(maxeval = 1000, print_level = 0),
        pvals,
        verbose = "silent"
    )$local_opt,
    nsim = nrow(pvals),
    x0 = x_sparse,
    verbose = "silent",
    objective_args = list(
        gain_tolerance = 0.01,
        ref_graph = list(hyp_weight = w_ref, trans_matrix = G_ref),
        u_range = .trial_success_range(ts)
    )
)
```

Output:

```text
start edges: 5 start gain: 0.795075 threshold: 0.79677675 start feasible: FALSE
COBYLA result edges: 6 gain: 0.7989 feasible: TRUE
```

### A15: nested calls compound the cap and the documentation says so

With seed `7005` for each call:

```r
set.seed(7005)
once <- graph_simplify(
    ref_object, pvals,
    gain_tolerance = 0.01,
    global_search = FALSE,
    control = small_control,
    verbose = "silent"
)
set.seed(7005)
twice <- graph_simplify(
    once, pvals,
    gain_tolerance = 0.01,
    global_search = FALSE,
    control = small_control,
    verbose = "silent"
)
stopifnot(identical(
    twice$sparsity$gain_reference,
    once$sparsity$gain
))
stopifnot(twice$sparsity$gain >= 0.99 * once$sparsity$gain)
```

Output:

```text
second reference equals first result: TRUE
second cap holds: TRUE
documentation says cap compounds: TRUE
```

The wording appears in `R/graph_simplify.R` and generated
`man/graph_simplify.Rd`.

### Design-record cases 3.11.1, 3.11.2, and 3.11.4

Independent calls produced:

```text
case 1 non-positive gain error: TRUE
case 2 warning emitted: TRUE
case 2 returns exact reference: TRUE
case 2 source/prune_loss: reference / 0
case 4 free edges: 4
case 4 one free edge per row: TRUE
```

The essential calls were:

```r
# Case 1
zero_ts <- trial_success(
    0 * r1 + 0 * r2 + 0 * r3 + 0 * r4,
    verbose = "silent"
)
zero_object <- ref_object
zero_object$trial_success <- zero_ts
tryCatch(
    graph_simplify(zero_object, pvals, verbose = "silent"),
    error = conditionMessage
)

# Case 2
set.seed(3)
pvals2 <- simulate_pvalues(
    c(0.9, 0.8),
    corr_matrix = diag(2),
    nsim = 2000
)
G2 <- matrix(c(0, 1, 1, 0), 2, 2, byrow = TRUE)
ref2 <- graph_optimal(
    c(0.5, 0.5), G2,
    constraints = graph_constraint_free(2),
    trial_success = trial_success(r1 + r2, verbose = "silent"),
    alpha = 0.025
)
graph_simplify(ref2, pvals2, verbose = "silent")

# Case 4
case4 <- graph_simplify(
    ref_object, pvals,
    gain_tolerance = 1,
    global_search = FALSE,
    control = small_control,
    verbose = "silent"
)
rowSums(unname(case4$trans_matrix) != 0)
```

The `lambda = 1` result had exactly one free edge in each row.

### `graph_optimise()` spot-checks matched `main`

Because the cited A1 script was unavailable, I ran only three cases rather
than repeating the claimed sweep:

| Case | `m` | Constraint | Global | Threads |
|---|---:|---|:---:|---:|
| `m2_local_serial` | 2 | free | no | 1 |
| `m3_global_constrained` | 3 | pinned rows/weights | yes | 1 |
| `m4_global_parallel` | 4 | free | yes | 2 |

Each worktree ran the same script, with p-value seeds `11`, `12`, and `13` and
optimisation seeds `101`, `102`, and `103`. Before comparison I removed only
the branch's documented new `alpha` and `sparsity` fields and the compiled
function pointer, which cannot be compared meaningfully after serialisation.
The rest of each complete result object was compared with `identical()`,
including every GA slot, every nloptr field, constraints, class and attributes.
Warnings, messages, console output, and `.Random.seed` were compared
separately.

Runner:

```r
args <- commandArgs(trailingOnly = TRUE)
setwd(args[[1]])
Sys.setenv(NOT_CRAN = "true")
pkgload::load_all(quiet = TRUE)

make_pvals <- function(seed, power) {
    withr::with_seed(seed, {
        m <- length(power)
        corr <- matrix(0.2, m, m)
        diag(corr) <- 1
        stats::pnorm(
            mvtnorm::rmvnorm(
                1000,
                mean = calc_ncp(power),
                sigma = corr
            ),
            lower.tail = FALSE
        )
    })
}
normalise_result <- function(x) {
    x$alpha <- NULL
    x$sparsity <- NULL
    x$trial_success$func <- NULL
    x
}
run_case <- function(seed, pvals, constraint, global_search, threads) {
    objective <- switch(
        as.character(ncol(pvals)),
        "2" = quote(r1 + r2),
        "3" = quote(r1 + r2 + r3 + 0.5 * (r1 && r2)),
        "4" = quote(r1 + r2 + r3 + r4 + 0.5 * (r1 && r2))
    )
    ts <- rlang::inject(trial_success(!!objective, verbose = "silent"))
    ctrl <- multigrain_control() |>
        control_global(run = 3, maxiter = 5, popSize = 20) |>
        control_local(maxeval = 100)
    warnings <- messages <- character()
    set.seed(seed)
    console <- capture.output(
        result <- withCallingHandlers(
            graph_optimise(
                pvals, constraint, ts,
                control = ctrl,
                global_search = global_search,
                num_threads = threads,
                verbose = "silent"
            ),
            warning = function(cnd) {
                warnings <<- c(warnings, conditionMessage(cnd))
                invokeRestart("muffleWarning")
            },
            message = function(cnd) {
                messages <<- c(messages, conditionMessage(cnd))
                invokeRestart("muffleMessage")
            }
        )
    )
    list(
        result = normalise_result(result),
        result_names = names(result),
        result_class = class(result),
        warnings = warnings,
        messages = messages,
        console = console,
        random_seed = .Random.seed
    )
}

gc3 <- graph_constraint(
    hyp_constraint = c(NA, NA, 0),
    trans_constraint = matrix(
        c(0, NA, NA, NA, 0, NA, 1, 0, 0),
        nrow = 3,
        byrow = TRUE
    )
)
results <- list(
    m2_local_serial = run_case(
        101, make_pvals(11, c(0.9, 0.8)),
        graph_constraint_free(2), FALSE, 1L
    ),
    m3_global_constrained = run_case(
        102, make_pvals(12, c(0.9, 0.82, 0.75)),
        gc3, TRUE, 1L
    ),
    m4_global_parallel = run_case(
        103, make_pvals(13, c(0.93, 0.9, 0.85, 0.8)),
        graph_constraint_free(4), TRUE, 2L
    )
)
saveRDS(results, args[[2]], version = 2)
```

Save that block as `a1_spotcheck.R`, then run:

```powershell
$branch = (Get-Location).Path
$main = 'C:\temp\multigrain-main'
git worktree add --detach $main main
Rscript "$branch\a1_spotcheck.R" $branch "$branch\a1_branch.rds"
Rscript "$branch\a1_spotcheck.R" $main "$branch\a1_main.rds"
```

It was run once in each worktree, then compared with:

```r
main_results <- readRDS("a1_main.rds")
branch_results <- readRDS("a1_branch.rds")

normalise_result <- function(x) {
    x$alpha <- NULL
    x$sparsity <- NULL
    x$trial_success$func <- NULL
    x
}

for (name in names(main_results)) {
    m <- main_results[[name]]
    b <- branch_results[[name]]
    added <- setdiff(b$result_names, m$result_names)
    removed <- setdiff(m$result_names, b$result_names)
    cat("\n", name, "\n", sep = "")
    cat("normalised result identical:",
        identical(m$result, b$result), "\n")
    cat("random seed identical:",
        identical(m$random_seed, b$random_seed), "\n")
    cat("warnings identical:",
        identical(m$warnings, b$warnings), "\n")
    cat("messages identical:",
        identical(m$messages, b$messages), "\n")
    cat("console identical:",
        identical(m$console, b$console), "\n")
    cat("class identical:",
        identical(m$result_class, b$result_class), "\n")
    cat("added names:", paste(added, collapse = ", "), "\n")
    cat("removed names:", paste(removed, collapse = ", "), "\n")
}
```

Output for all three cases:

```text
normalised result identical: TRUE
random seed identical: TRUE
warnings identical: TRUE
messages identical: TRUE
console identical: TRUE
class identical: TRUE
added names: alpha, sparsity
removed names:
```

The report's description of A1 does not say it recorded top-level
`constraints`, `trial_success`, or `global_search`; most GA slots; most nloptr
fields; or emitted conditions/output. Those are observable paths that its
listed element-by-element comparison could miss. The three spot-checks above
included them and found no difference.

## Not tested

- **The claimed 100-case A1 sweep itself.** Its script is unavailable, so I
  could not judge its exact fixture construction, comparison code, or raw
  output. Settling this requires adding the cited script and its two result
  files, or rerunning it from a preserved copy. The three spot-checks above are
  not a substitute for 100 cases.
- **HTML rendering, full `R CMD check`, A12 at `m = 8`/`nsim = 1e6`, O1
  mutation-rate tuning, O4 `run` tuning, and the runtime multiplier.** These
  remain outside the requested review scope or explicit run limits, matching
  the implementation report's exclusions.
- **A proof of Nelder-Mead behavior for all objectives.** The requested logged
  `m = 4` run held for all nine optimiser steps. A proof would require reasoning
  about every simplex trajectory and tie case, beyond the requested
  instrumentation.
