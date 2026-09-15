# Design record: group sequential extension of `multigrain`

Status: design specification, not yet implemented. It supersedes `debug/docs/stale-design.md` (April 2026) in full, including the three decisions that document called governing; section 2 says why. Every design element is stated, followed by the reasoning behind it and the alternatives it was chosen over. Elements are marked **LOCKED** when the reasoning is considered settled and **OPEN** when evidence or a decision is still needed; open items are collected in section 10. Claims about `gsDesign`, `graphicalMCP` and the existing package were checked by running code; the scripts and their output are in the appendices.

Part I is written for a developer who knows multiplicity and graphical procedures but not group sequential methodology. Part II is the specification.

Notation: $m$ hypotheses, $K$ analyses ("looks"). $p_{i,k}$ is the nominal one-sided p-value of $H_i$ at look $k$; $t_{i,k}\in(0,1]$ its information fraction there; $\phi_i(a,t)$ its alpha-spending function (cumulative level spent by information fraction $t$ when the hypothesis holds level $a$). $\alpha^*_{i,k}(a)$ is the nominal boundary: the largest p-value at look $k$ that rejects when $H_i$ holds level $a$. $w_i(I)$ is the current weight of $H_i$ given the set $I$ of hypotheses not yet rejected, so its level is $w_i(I)\alpha$. $p^r_{i,k}$ is the repeated p-value and $p^s_{i,k}=\min_{l\le k}p^r_{i,l}$ the sequential p-value. $\tau_i\in\{1,\dots,K,\infty\}$ is the look at which $H_i$ is first declared rejected. $\psi(\mathbf r,\boldsymbol\tau)$ is the gain function; $U(\mathbf w,\mathbf G)$ its expectation, estimated on a fixed set of $N$ simulated trials.

---

# Part I. The GSD extension in plain language

## What a group sequential design changes

In the fixed-sample package each simulated trial produces one p-value per hypothesis and the graph decides which hypotheses are rejected. In a group sequential design the trial is analysed several times as data accumulate, say at one third, two thirds and all of the planned information. So each hypothesis now has one p-value *per analysis*, and a simulated trial is an $m\times K$ table of p-values rather than a row of $m$.

Testing the same hypothesis repeatedly at the same level inflates the type I error, so the local level $w_i\alpha$ that a hypothesis holds has to be spread across the $K$ looks. An *alpha-spending function* says how much of it may be used by each look. From the spending function, the information fractions, and the known correlation between a hypothesis's test statistics at successive looks, one can compute a *nominal boundary* for each look: the p-value threshold that spends exactly the right amount. The boundary is smaller than $w_i\alpha$ at every look before the last, and it depends on $w_i\alpha$ in a monotone way: more alpha, larger boundary.

Maurer and Bretz (2013) showed that the graphical procedure carries over almost unchanged. At each analysis you run the familiar cascade (reject, recycle alpha along the edges, retest) but compare each p-value with its boundary instead of with $w_i\alpha$. The updated graph is carried into the next analysis, and a hypothesis that has been rejected stays rejected. Because the gain function cares when a claim is established (an interim rejection is worth more than a final one), the outcome of a trial is now a rejection vector **and** a vector of decision times, one per hypothesis, with "never" as a possible value.

## What made the old design complicated

The boundary a hypothesis is compared against depends on its *current* weight, which changes every time alpha is recycled, and computing a boundary is expensive (a multivariate normal probability and a root-find). The stale design therefore tabulated boundaries over a grid of allocated alpha before the search and interpolated inside the kernel each time a hypothesis was tested. It worked on paper, but it put a table lookup on the hottest line of the C++ kernel, threaded spending functions and information fractions all the way into the optimiser, and treated look-back as an extra kernel feature to defer.

## The key idea: convert every p-value once, then reuse the fixed-sample kernel

Because the boundary is monotone in the allocated level, the question "is $p_{i,k}$ below the boundary for level $a$?" has an equivalent form: "what is the smallest level $a$ at which $p_{i,k}$ would cross its boundary?" That number is the *repeated p-value* $p^r_{i,k}$ (Maurer and Bretz, section 3.3). It depends only on the p-value, the spending function and the information fractions. It does **not** depend on the graph. And the group sequential rejection rule "$p_{i,k}$ below the boundary at level $w_i(I)\alpha$" is exactly the fixed-sample rule "$p^r_{i,k}\le w_i(I)\alpha$".

So the whole extension is:

1. Before optimisation, convert every simulated p-value at every analysis into its repeated p-value. This is a one-time transform of the $N\times m\times K$ array.
2. To evaluate a graph on one trial, run the existing fixed-sample cascade on the analysis-1 column of that trial, keep the updated graph, run it again on the analysis-2 column, and so on. Record the pass at which each hypothesis fell.

The kernel never sees a boundary, a spending function or an information fraction. The only kernel change is an outer loop over analyses and a second output matrix of decision times. The optimiser (encoding, genetic algorithm, local search, pruning) is untouched; each evaluation costs at most $K$ times what it costs today.

A useful analogy: adjusted p-values let you compare every hypothesis against a single $\alpha$ instead of a per-hypothesis threshold. Repeated p-values do the same for looks: they express each look's evidence in fixed-sample currency.

This was checked, not assumed. Inverting `gsDesign` boundaries reproduces the repeated p-values that Maurer and Bretz published (computed with ADDPLAN) to four significant figures, and the existing cascade on those values reproduces their worked example (Appendix A). A direct implementation that computes boundaries on the fly agreed with transform-then-cascade on every one of 1200 hypothesis-by-trial decisions and decision times, in four look-back configurations (Appendix B).

## Look-back comes for free

A variant of the procedure (`graphicalMCP` calls it `look_back`; Zhao et al. 2025) lets a hypothesis that receives extra alpha at a later analysis be rejected on the strength of its evidence at *any* earlier analysis. In converted currency that is simply the running minimum of $p^r_{i,1},\dots,p^r_{i,k}$ (the *sequential p-value*): one `cummin` along a row, switchable per hypothesis. No kernel change. It is off by default, matching Algorithm 1 of Maurer and Bretz and the `graphicalMCP` default.

## Two conventions the implementer must know

- **Endpoints that finish early.** If an endpoint's data are complete at analysis $k$ (information fraction 1), its p-value and its boundary do not change afterwards, but it can still be rejected at a later analysis if recycling gives it more alpha. That is exactly how PFS behaves in the manuscript's Example 5. In converted currency: copy $p^r$ forward from the analysis at which the endpoint matured. Never ask `gsDesign` for a boundary with a repeated information fraction of 1: the zero spending increment makes it return an unrejectable boundary (Appendix B, first block).
- **Endpoints not yet tested at an early analysis.** Set $p^r=1$ ("cannot reject") for that analysis. No data mask is needed.

## How the conversion is computed

For each hypothesis, take a log-spaced grid of allocated levels from $10^{-14}$ up to $\alpha$. At each grid level, `gsDesign::gsBound1()` returns the nominal boundary at every look (deterministic quadrature, 2.5 to 5 ms per call on this machine). That gives a monotone table of boundary against level for each look, which is inverted by interpolation to map any p-value to its level. Interpolation error is about $10^{-6}$ relative with 1024 grid points; converting a million trials takes about a second. If a table is not monotone the user's spending function is not "well ordered", which means the group sequential graphical procedure itself is invalid, and the transform stops with an error.

## Gain functions with time

`trial_success_gsd()` speaks the same language as `trial_success()` plus three additions: symbols `t1 … tm` for the analysis at which each hypothesis was rejected (0 = never); comparisons such as `t1 == 1`; and discount tables written as `d(t1)`, where `d` is a numeric vector of length $K$ supplied alongside the expression and `d(0)` is always 0. Example 5 becomes `v_pfs * d(t1) + v_os * d(t2)` with `d = c(1, 0.75)`.

## The pipeline

```
simulate_pvalues_gsd()        raw p-values, N x m x K, canonical joint normal model
        |
transform_pvalues_gsd()       repeated (or sequential) p-values + boundary tables
        |                     <- the only place spending functions and gsDesign appear
graph_optimise_gsd()          same GA / COBYLA / pruning as today; objective =
        |                     graph_shortcut_gsd() (K passes) + trial_success_gsd()
calc_power_pvals_gsd()        power per hypothesis and per analysis, decision-time
        |                     distribution, expected gain
multigrain_graph_optimal      same object class as today
```

---

# Part II. Design record

## 1. Problem

Extend `multigrain` so that a candidate graph $(\mathbf w,\mathbf G)$ can be evaluated under a group sequential design with a time-dependent gain $\psi(\mathbf r,\boldsymbol\tau)$ (manuscript section 2.3.4, Example 5, supplement S2.3 and S3), and optimised with the existing search. Spending functions, information fractions and analysis times are fixed inputs; the decision variables are the graph. The user's three questions were whether boundary evaluation belongs in the hot path at all, whether inverting the monotone boundary map and transforming the p-value array once reduces the problem to the fixed-sample kernel, and whether look-back deserves deferral or falls out of a better formulation. The answers are no, yes, and it falls out.

## 2. Why the stale design is superseded

**Stale decision 1.1 (boundary grid interpolated in the kernel).** It solved the wrong problem. It kept the comparison "p-value against boundary at the current level" and made the boundary cheap to fetch. The right move is to invert the boundary map once so the comparison becomes "repeated p-value against the current level", which the existing kernel already performs. Its supporting argument (a per-candidate closure cache would take about 71 hours at $m=8$) is true and irrelevant, because neither the cache nor the grid is needed in the search. The `has_data` mask and the `bounds_grid` kernel argument disappear with it.

**Stale decision 1.3 (`look_back = FALSE` only; look-back deferred as a 20 to 30 percent kernel cost).** Under the transform, look-back is a running minimum along the analysis axis of the transformed array, applied per hypothesis before the search. It costs nothing in the kernel and nothing in the search. Deferral was a consequence of decision 1.1, not a property of the problem.

**Stale decision 1.2 (extend the `trial_success()` parser).** Reasonable on its own; the user has chosen a parallel `_gsd` family instead (section 4.7), which keeps the fixed-sample surface untouched while branch `plan-for-sparsity` is changing the same files.

Two further stale assumptions do not survive: that the kernel needs the spending functions (it does not), and that `graphicalMCP::gs_boundaries()` should be the boundary engine (it uses `mvtnorm::pmvnorm` with the randomised Genz–Bretz algorithm and a per-value `uniroot`, which is fine for one trial and infeasible for $N\cdot m\cdot K$ values; `gsDesign::gsBound1()` is deterministic quadrature and one call yields all looks).

## 3. How the package evaluates and searches today

The facts below constrain the design; each was verified by reading the current tree (line numbers checked on 2026-09-15).

`create_obj_func()` in `R/objective_function.R:24-89` is the single objective factory. Its closure decodes the parameter vector through `split_theta()`, `recover_full_weights()` and `recover_full_trans_matrix()`, returns a negative penalty for `NA` or out-of-range entries, zeroes hypothesis weights below $10^{-4}$ and transition entries below $10^{-5}$ (`:66-67`), runs `graph_shortcut()` or `graph_shortcut_parallel()` (`:69-85`) and returns `power_criterion(rej_matrix)`.

`graph_shortcut()` in `src/graph_shortcut.cpp:62-166` runs, per trial, a `while(true)` cascade: scan hypotheses in index order, reject the first with `cur_p[i] < cur_a[i]` (`:110`, strict), propagate `cur_a`, update the working graph with the Bretz formula, zero the rejected row and column, set `cur_a[rej] = 0` (`:161`), repeat until no rejection. The parallel worker is the same code per chunk of trials. Output is an $N\times m$ logical matrix.

`graph_optimise()` in `R/optimisation.R:88-239` requires a double matrix (`:105`), takes $m$ from `ncol(pvals)` (`:124`, `:133`), calibrates control on the matrix (`:136`; `control_prepare()` in `R/control_prepare.R:73-86`), runs the GA on `pvals[sample(nsim), ]` (`:283`) and COBYLA on another such subsample (`:403`), re-evaluates both on the full matrix with `graph_shortcut()` (`:332`, `:452`), chooses with `choose_graph()`, prunes with `prune_graph()` and reports `calc_power_pvals()` on the pruned graph. `prune_graph()` (`R/post_optim_processing.R:481-532`) drives `.try_prune()` (`:218-255`), which calls `calc_power_pvals()` (`R/calc_power.R:114-164`; requires a double matrix at `:124`, calls `graph_shortcut` at `:147`, evaluates custom measures at `:212-230` by passing the rejection matrix as the single argument).

`trial_success()` in `R/trial_success.R` accepts only `r<digit>` symbols and the operators `+ - * / ( && ||` (`:150-192`), counts $m$ from `r[0-9]+` matches (`:314-353`), transforms the expression with a two-type ("bool", "real") system in `parse_and_transform()` (`:382-471`; unknown calls fall to "real" at `:436`) and compiles `double powerFunc(LogicalMatrix x)` (`:218`).

`check_double_matrix()` in `R/check_types.R:41-66` rejects arrays and S3 objects. `summary.multigrain_graph_optimal()` prints `local_power` as a vector (`R/graph_optimal.R:127-130`).

Branch `plan-for-sparsity` rewrites `create_obj_func()`, `.try_prune()`, `choose_graph()` and `R/mutation_helpers.R`. Anything in this design that touches those files will conflict.

## 4. Specification

### 4.1 The transform: repeated p-values from `gsDesign` boundary tables

**Inputs.** A raw array of nominal p-values $P$ of dimension $N\times m\times K$; an $m\times K$ matrix $T$ of information fractions with `NA` where a hypothesis has no data at a look; a spending function per hypothesis (one function recycled, or a list of $m$); the overall level $\alpha$; a logical `look_back` of length 1 or $m$; a grid size $G$.

**Per hypothesis $i$.** Let $D_i=\{k: T_{i,k}\ne\mathrm{NA}\}$ (non-empty, information fractions positive and non-decreasing over $D_i$). Let $k^{\text{mat}}_i$ be the first look in $D_i$ with $T_{i,k}\ge1$, or the last look in $D_i$ if none. Let $L_i=\{k\in D_i: k\le k^{\text{mat}}_i\}$ and $t_i=T_{i,L_i}$; these are the looks that carry distinct information.

1. *Short-circuit.* If $|L_i|=1$ and $t_i\ge1$, the boundary at level $a$ is $a$ itself, so $p^r_{i,k}=p_{i,k_0}$ exactly for every $k\ge k_0$, where $k_0$ is that look. No table is built. This covers $K=1$ and the PFS endpoint of Example 5, and is what makes $K=1$ bit-identical to the fixed-sample kernel (`gsBound1` returns $0.012499999999999968$ for level $0.0125$; Appendix B).
2. *Table.* Otherwise take $\gamma_g=\exp\{\log\gamma_{\min}+(g-1)(\log\alpha-\log\gamma_{\min})/(G-1)\}$, $g=1..G$, with $\gamma_{\min}=10^{-14}$. For each $\gamma_g$ compute the per-look spend increments $\Delta_l=\phi_i(\gamma_g,t_{i,l})-\phi_i(\gamma_g,t_{i,l-1})$ and the boundaries
   `b <- gsDesign::gsBound1(theta = 0, I = t_i, a = rep(-20, length(t_i)), probhi = Delta)$b`, `B[g, ] <- pnorm(b, lower.tail = FALSE)`. One call gives all looks in $L_i$ because boundaries at look $l$ depend only on looks $\le l$.
3. *Well-ordering check.* Assert every column of $B$ is non-decreasing in $g$. This is condition (2) of Maurer and Bretz; if it fails the sequentially rejective procedure is not valid for that spending function, and the transform aborts naming the hypothesis. All seven `gsDesign` spending families checked (`sfLDOF`, `sfLDPocock`, `sfHSD` at $\pm4$, `sfPower` at 0.5 and 3, `sfLinear`) pass.
4. *Inverse.* For $k\in L_i$, with $l$ its position in $L_i$: drop duplicated and non-positive entries of $B_{\cdot,l}$ (the first-look LDOF boundary underflows to a constant $2.75\times10^{-89}$ below $\gamma\approx10^{-6}$), then $p^r_{i,k}=\exp\{\mathrm{approx}(\log B_{\cdot,l},\log\gamma,\ \log p_{i,k})\}$. If $p_{i,k}\ge B_{G,l}$ set $p^r=1$ (never rejectable, since no allocation exceeds $\alpha$). If $p_{i,k}$ is below the smallest positive table entry set $p^r=\gamma_{\min}$.
5. *Matured looks.* For $k>k^{\text{mat}}_i$: $p^r_{i,k}=p^r_{i,k^{\text{mat}}_i}$. Any non-`NA` raw p-value the user supplies after maturity is ignored; if it differs from the maturity p-value the function warns.
6. *Looks without data.* For $k\notin D_i$ with $k<k^{\text{mat}}_i$: $p^r_{i,k}=1$.
7. *Look-back.* If `look_back[i]`, replace row $i$ of the transformed array by its running minimum along $k$. With step 6 this reproduces the `graphicalMCP` carry-forward of the last sequential p-value across a look without data.

**Output.** An S3 object `multigrain_pvals_gsd` holding the transformed array (numeric, `dim = c(N, m, K)`), `m`, `K`, `alpha`, `info_frac`, `look_back`, a label per spending function, and the tables $B_i$ with their grids (for printing and for reporting the boundaries of a final graph). The kernel receives the array with `dim <- c(N, m * K)`; column $(k-1)m+i$ is hypothesis $i$ at look $k$, which R's column-major layout gives without copying.

**Rationale.** LOCKED. Correctness rests on section 3.3 of Maurer and Bretz: under well-ordered boundaries $p_{i,k}\le\alpha^*_{i,k}(a)\iff p^r_{i,k}\le a$, and $p^r$ does not involve the graph. Appendix A reproduces their Table 1 boundaries and Table 2 repeated p-values from `gsBound1` to four significant figures and reproduces the rejections in their case study. Appendix B shows the whole chain (transform, then $K$ passes of the fixed-sample cascade) agreeing with a direct implementation on 0 of 1200 decisions and 0 of 1200 decision times in each of four look-back configurations, including one hypothesis maturing at look 2 and one at look 1.

The floor is $\gamma_{\min}$ rather than 0 because the kernel's snaps leave allocations as small as $10^{-4}\cdot10^{-5}\cdot\alpha=2.5\times10^{-11}$, and about 0.9 percent of simulated rows at noncentrality 4 have $p<10^{-10}$; a 0 floor would reject those rows at any positive allocation even when the true $p^r$ exceeds it. Appendix B shows the case $p=5\times10^{-11}$, allocation $2.5\times10^{-11}$: the old rule rejects, the new rule does not. The residual effect of the $\gamma_{\min}$ floor is a conservative non-rejection when the true $p^r$ lies below $\gamma_{\min}$ and the allocation lies between the two, which cannot happen for allocations above $10^{-14}$.

Grid size: with levels on $[10^{-14},\alpha]$ and log-log linear interpolation, the maximum relative error against `uniroot` at 200 random p-values was $1.0\times10^{-5}$ ($G=256$), $4.0\times10^{-6}$ ($G=512$) and $1.5\times10^{-6}$ ($G=1024$). Monte Carlo noise in $U$ at $N=10^5$ is of order $10^{-3}$; an interpolation error only matters when a repeated p-value lies within its relative error of a threshold, which at $10^{-6}$ affects of order one row in a million. Default $G=1024$; OPEN whether 512 suffices (section 10).

Spending function interface: any `function(alpha, t)` returning a numeric vector of cumulative spend or an object with a `$spend` element, evaluated at the full vector $t_i$. `gsDesign` functions work directly (`sfLDOF`, `sfLDPocock`, `sfPower`; parametrised ones as `function(a, t) gsDesign::sfHSD(a, t, param = -2)`), and so do those of `graphicalMCP`. The transform checks $\phi_i(a,1)=a$ within tolerance (a spending function that does not spend everything by full information would silently waste alpha).

**Alternatives.** Direct root-finding per value as `graphicalMCP::repeated_p()` does: one `uniroot` over `pmvnorm` per value, tens of milliseconds each, infeasible at $N\cdot m\cdot K$. Boundaries from `mvtnorm` in-house: `pmvnorm` is randomised quasi-Monte Carlo above three dimensions, so the transform would not be reproducible bit-for-bit for four or more looks, and a root-finder would have to be written and maintained. A linear grid in $\gamma$ (the stale design): wastes points where boundaries are large and is coarse where they are tiny; the log grid is uniform in relative error. Clamp to 0 below the table: rejected for the reason above.

### 4.2 Semantics: repeated by default, sequential per hypothesis

**Decision.** LOCKED (user decision). Default `look_back = FALSE` for every hypothesis: a hypothesis is rejected at look $k$ only on the strength of its look-$k$ evidence, as in Algorithm 1 of Maurer and Bretz and the `graphicalMCP` default. `look_back` may be a logical vector of length $m$; a `TRUE` entry switches that hypothesis to sequential p-values.

**Rationale.** Both semantics strongly control the FWER (section 4.9). The sequential variant gives uniformly more rejections and earlier decision times, hence higher expected gain for any monotone $\psi$, but is less established with regulators and is not what Algorithm 1 does. The manuscript's wording in section 2.3.4 invokes the sequential framework only for the matured-endpoint case, where the two semantics coincide (the p-value is the same at every later look). Making the switch per hypothesis costs nothing and mirrors `graphicalMCP`. The protocol must state which semantics the trial will use; the package does not decide that.

### 4.3 The kernel: the fixed-sample cascade with an outer loop over looks

**Decision.** LOCKED. New file `src/graph_shortcut_gsd.cpp` with `graph_shortcut_gsd(pvals, alpha, w, G, K)` and `graph_shortcut_gsd_parallel(pvals, alpha, w, G, K, num_threads, grain_size)`. `pvals` is the $N\times mK$ matrix described in 4.1. Per trial: initialise `cur_a = w * alpha`, `cur_G = G`, a rejected flag per hypothesis, `sumrej = 0`; then for `k = 0..K-1`: load `cur_p[i] = pvals(set, k*m + i)`, run the existing cascade verbatim (scan in index order, reject the first unrejected $i$ with `cur_p[i] < cur_a[i]`, propagate, update the graph, zero row and column, set `cur_a[rej] = 0`, repeat), recording `time(set, i) = k + 1` and `rejected(set, i) = TRUE` at each rejection; stop early when `sumrej == m`. Return a list with `rejected` ($N\times m$ logical) and `time` ($N\times m$ integer, 0 = never). The parallel worker is the same code per chunk; the auto-tuned grain size uses `K * m^3` operations per row in place of `m^3`.

**Rationale.** The graph state after the last rejection at look $k$ is exactly the state Maurer and Bretz carry into look $k+1$; a rejected hypothesis stays rejected because its weight is zero and its flag is set. Strict `<` is kept from the fixed kernel: it keeps a hypothesis with weight exactly zero from rejecting a floored $p^r$, and it makes the $K=1$ case identical rather than merely equal. (`graphicalMCP` uses `<=`; the two differ on a set of measure zero.) The kernel takes no spending function, information fraction or mask; every convention of section 4.1 is already encoded in the array.

**Alternatives.** Interpolated boundary lookup in the scan (stale design): rejected in section 2. A separate kernel call per look with graph state passed back to R: $K$ round trips per candidate and per-look allocation of output; rejected. A `<=` comparison to match `graphicalMCP`: rejected for the two reasons above.

### 4.4 Decision time

**Decision.** LOCKED. $\tau_i$ is the look index at which the kernel processed the rejection of $H_i$. Under both semantics this equals the manuscript's definition (the first analysis at which the procedure declares $H_i$ rejected) and `decision_at` in `graphicalMCP`. In Example 5 a PFS rejection that becomes possible only when OS recycles alpha at the final analysis gets $\tau_{\text{PFS}}=2$, as the manuscript specifies. The `first_rejected_at` quantity of `graphicalMCP` (the earliest look whose boundary the hypothesis crossed, under look-back) is a reporting quantity, is not needed by the gain, and can be recovered from the transformed array if a user asks.

Calendar times $s_k$ enter only through the user's discount tables (section 4.6); the package stores them, if supplied, as labels.

### 4.5 The search is unchanged

**Decision.** LOCKED. Encoding of $(\mathbf w,\mathbf G)$, constraint handling, Cauchy population and mutation, GA with intermittent Nelder–Mead, COBYLA refinement, `choose_graph()`, greedy pruning: all reused as they are. Only the objective closure changes: it calls `graph_shortcut_gsd()` and hands the `time` matrix to the compiled gain.

**Rationale.** The objective is still a piecewise-constant function of $(\mathbf w,\mathbf G)$ on fixed draws, of the same class as today, with more pieces. Common random numbers still apply (the transformed array is fixed). Per-evaluation cost is at most $K$ times the fixed-sample cost, and the transform is amortised over the whole run. Nothing about the landscape asks for a different optimiser.

**What would change under joint design optimisation** (out of scope; user decision). If spending parameters or interim timings became decision variables, the transform would move inside the objective (one boundary table per design, seconds each) and a change of timing would change the correlation structure and hence require re-simulation on fixed standard normals, as the sample-size minimisation already does. The natural extension path within this design is one transformed array per candidate design and an outer loop over a discrete set of designs.

### 4.6 The time-dependent gain function

**Decision.** LOCKED (language); OPEN (function name, section 10). `trial_success_gsd(objective, ..., K = NULL)` compiles an expression over:

- `r<i>`: rejection indicator of $H_i$ (bool), defined as `t<i> > 0`;
- `t<i>`: the look at which $H_i$ was rejected, an integer in $0..K$ with 0 meaning never;
- numeric literals and values injected with `!!`, as today;
- `+ - * /`, `&& || and or`, parentheses, as today;
- comparisons `== != < <= > >=` between a `t<i>`, a literal or an arithmetic expression; the result is bool;
- discount tables: any name passed through `...` as a numeric vector of length $K$ may be applied as `name(t<i>)`; it compiles to a lookup with index 0 mapped to 0.

The generated C++ has signature `double f(IntegerMatrix t)`; `r<i>` becomes `double(t(i, idx) > 0)`, `t<i>` becomes `double(t(i, idx))`, a comparison becomes `double(A op B)`, and each table becomes `static const double name_tab[] = {0.0, v1, ..., vK}` indexed by `t(i, idx)`. $m$ is the largest index over `r` and `t` symbols; $K$ is the common length of the tables, or the `K` argument, or unset until first use, when it is checked against the p-value object. The object has class `c("multigrain_trial_success_gsd", "multigrain_trial_success")` so the existing constructor checks, `print()` and `summary()` apply.

Examples the parser must accept (values injected with `!!` where they are variables):

- Example 5: `trial_success_gsd(!!v_pfs * d(t1) + !!v_os * d(t2), d = c(1, 0.75))`.
- Supplement, dual primaries with a same-look bonus: `!!nu1 * d(t1) + !!nu2 * d(t2) + !!nu12 * ((t1 == 1 && t2 == 1) + !!delta12 * (t1 == 2 && t2 == 2))`, with `d = c(1, delta)`.
- Supplement, co-primaries at the same look: `!!nu12 * ((t1 == 1 && t2 == 1) + !!delta12 * (t1 == 2 && t2 == 2))`.
- Supplement, PFS (`t1`) and OS (`t2`): `(t2 == 1) * (!!a1 + !!b1 * (t1 == 1)) + (t2 == 2) * (!!a2 + !!b2 * (t1 == 1) + !!c2 * (t1 == 2)) + (t2 == 0) * (t1 == 1) * !!b3`.

**Rationale.** Comparisons alone can express every gain in the manuscript and supplement; discount tables are sugar for the common step-function case and make Example 5 read as it is written. Using the time matrix as the sole input keeps one compiled signature and makes `r<i>` a derived quantity, so the two symbols can never disagree. A separate function and class, rather than extending `trial_success()`, is the user's decision and also avoids `R/trial_success.R`, which the sparsity branch touches. The parser is a copy of the fixed-sample one with three additions; the changes needed at each site are listed in section 6, step P3.

**Alternatives.** A user-supplied R function on `(rejected, time)`: cannot be compiled, and would be evaluated per row in R inside the GA; retained only as the slow path of `custom_power` in `calc_power_pvals_gsd()`, as today. A calendar-time argument on the gain object: unnecessary, the table already maps look index to value. Replacing the `t<i> = 0` convention by `NA` or `K + 1`: `0` indexes the table naturally and matches the stale design.

### 4.7 API surface: a parallel `_gsd` family

**Decision.** LOCKED (shape, user decision); names OPEN.

- `simulate_pvalues_gsd(power_nominal, ..., alpha = 0.025, corr_matrix, info_frac, nsim = 1e5)`. `power_nominal` and `corr_matrix` refer to the final-look statistics as in `simulate_pvalues()`; `info_frac` is a length-$K$ vector or an $m\times K$ matrix with `NA` for looks without data. Simulates the canonical joint model of supplement S2.3, $E[Z_{i,k}]=\Delta_i\sqrt{t_{i,k}}$ and $\mathrm{Corr}(Z_{i,k},Z_{j,l})=\rho_{ij}\sqrt{\min(t_{i,k},t_{j,l})/\max(t_{i,k},t_{j,l})}$, over the distinct information levels only, then copies matured columns forward so that the covariance is never singular. Returns the raw array with `info_frac` attached.
- `transform_pvalues_gsd(pvals, ..., info_frac = NULL, spending, alpha = 0.025, look_back = FALSE, grid_size = 1024L)`: section 4.1; exported so that users with their own simulator, and the test suite, can call it directly. `info_frac` defaults to the attribute set by the simulator.
- `trial_success_gsd()`: section 4.6.
- `graph_optimise_gsd(pvals, graph_constraint, trial_success, ..., alpha, start_graph, global_search, num_threads, control, verbose)`: the same signature and flow as `graph_optimise()`, with `pvals` a `multigrain_pvals_gsd`. Internally: `create_obj_func_gsd()` (a copy of the decode-and-penalise block calling the GSD kernel and passing `time` to the gain), `.graph_optimise_ga_gsd()` and `.graph_optimise_local_gsd()` (subsampling on the first array dimension with `drop = FALSE`; re-evaluation with the GSD kernel), `choose_graph()` unchanged, `prune_graph_gsd()` and `calc_power_pvals_gsd()`. Control calibration reads $N$ and $m$ from the object through a new internal `control_prepare_dims()` rather than from a matrix. Returns a `multigrain_graph_optimal`.
- `calc_power_pvals_gsd(pvals, hyp_weight, trans_matrix, ..., alpha, custom_power, sum_to_one_constraint)`: returns `local_power` (probability of rejecting each hypothesis at any look), `local_power_by_analysis` ($m\times K$, cumulative by look), `exp_rejections`, `disj_power`, `conj_power`, `mean_decision_time` (over rejected trials), `time_distribution` ($m\times(K+1)$, proportions at "never" and looks $1..K$), and one entry per `custom_power` measure. A `multigrain_trial_success_gsd` entry is evaluated by its compiled function on the time matrix; a plain R function receives the integer time row.
- `prune_graph_gsd()`: the two greedy loops of `prune_graph()` with an acceptance test that evaluates `calc_power_pvals_gsd()`.
- `print()` and `summary()` methods for `multigrain_pvals_gsd` and `multigrain_trial_success_gsd`.

**Rationale.** The family is additive: no existing exported function, internal helper or kernel changes, so `graph_optimise()` on `main` and on the feature branch are the same code, and the gate in section 6 (identical results for the same seed) is a check on that fact rather than a fixture. The cost is duplicated plumbing (objective closure, two optimiser drivers, pruning acceptance); section 10 records the follow-up to merge the two paths once the sparsity branch has landed. Returning the same result class means `plot()`, `print()`, `summary()` and `graph_optimal_get_control()` work unchanged; `summary()` will print the by-look power matrix acceptably.

**Alternatives.** Extend `simulate_pvalues()`, `trial_success()`, `graph_optimise()` and `calc_power_pvals()` in place with $K=1$ as the special case (one code path; rejected by the user for now, and it would collide with the sparsity branch). A single `_gsd` entry point that takes raw p-values and spending functions (hides the transform; users could not inspect repeated p-values or reuse one simulation under several spending choices).

### 4.8 Dependencies

**Decision.** LOCKED (user decision). `Imports: gsDesign` for `gsBound1()` and its spending functions. `Suggests: graphicalMCP (>= 0.3.0)` as the test oracle only (`repeated_p()`, `sequential_p()`, `graph_test_shortcut_gsd()`, all exported since 0.3.0, released to CRAN on 27 August 2026). `mvtnorm` is already imported and is used by the simulator.

**Rationale.** `gsBound1()` is deterministic, one call returns all looks, and it is the standard implementation in the field, so a reviewer can reproduce every boundary with `gsDesign::gsDesign()` (Appendix A cross-checks the two). `graphicalMCP` computes one root per value with a randomised integrand; correct for a trial, not for a transform, and its conventions for matured hypotheses differ (section 4.9).

### 4.9 Error control and conventions

The group sequential graphical procedure strongly controls the FWER at $\alpha$ when the boundaries are well ordered (Maurer and Bretz, section 3.1). The transform changes nothing about which hypotheses are rejected: it is a re-expression of the same test. Three conventions need their own justification.

*Matured endpoints keep their last p-value and boundary.* Once $t_{i,k}=1$ the statistic is fixed and the spending function has spent the whole level, so for any level $a$ the event "reject $H_i$ at some look" is $\{p_{i,k^{\text{mat}}}\le\alpha^*_{i,k^{\text{mat}}}(a)\}$, whose null probability is exactly $a$. Re-testing the same p-value at a later look against the boundary for a larger level $a'$ is the same event at level $a'$; the Bonferroni bound $\sum_{i\in J}w_i(J)\alpha$ on each intersection is unchanged, and consonance follows from monotonicity of the event in $a$. This is the manuscript's Example 5 convention. `graphicalMCP` instead pads later looks with `NA` and, with `look_back = FALSE`, never re-tests a matured hypothesis; with `look_back = TRUE` it carries the sequential p-value forward, which for a matured hypothesis is the same number. So the oracle reproduces Example 5 only with `look_back = TRUE` on the matured hypothesis, and every oracle test must say which convention it exercises.

*Look-back.* The union over looks $l\le k$ of $\{p_{i,l}\le\alpha^*_{i,l}(a)\}$ has null probability $\phi_i(a,t_{i,k})\le a$ by construction of the boundaries, so the per-hypothesis bound and the closed-test argument hold, and the event is monotone in $a$, so the shortcut remains valid. Mixing semantics across hypotheses is fine because the bound is per hypothesis.

*Looks without data.* $p^r=1$ is "cannot reject", which is what not testing means.

The strict `<` and the $\gamma_{\min}$ floor can only withhold a rejection, never add one, so neither affects error control.

### 4.10 What does not change, and how that is verified

No file under `R/` or `src/` that exists today is edited except `R/RcppExports.R` (regenerated) and `DESCRIPTION`; new code lives in new files. The gate is an end-to-end run of `graph_optimise()` with a fixed seed on `main` and on the branch, compared with `identical()` on the returned object, plus `identical()` on `.Random.seed` afterwards.

## 5. Rejected alternatives (summary)

| Alternative | Reason rejected |
|---|---|
| Boundary grid interpolated inside the kernel (stale 1.1) | Solves the wrong problem; the inversion removes boundaries from the search entirely |
| Per-candidate closure cache of exact boundaries | Same, and $m\cdot2^{m-1}$ boundary sets per candidate |
| `graphicalMCP::repeated_p()` as the transform | One `uniroot` over randomised `pmvnorm` per value; infeasible at $N\cdot m\cdot K$ |
| In-house boundaries via `mvtnorm` | Not bit-reproducible above three looks; a root-finder to maintain; no gain over `gsDesign` |
| `look_back = TRUE` as default | Not Algorithm 1; less established with regulators; user decision |
| Deferring look-back | It is one `cummin` |
| A `has_data` mask passed to the kernel | Encoded as $p^r=1$ |
| `<=` in the kernel | Breaks $K=1$ identity; lets a zero weight reject a floored value |
| Clamp to 0 below the table | Spurious rejections at tiny recycled allocations (Appendix B) |
| Extending the fixed-sample functions in place | User decision; conflicts with the sparsity branch |
| Joint optimisation of spending or timing | Out of scope; extension path noted in 4.5 |

## 6. Implementation plan

Each step names the files it creates, the gate that closes it, and an effort estimate. Steps are sequential; nothing in a later step is needed to close an earlier gate. Only new files are created except where stated. The companion brief for the implementation agent is `debug/implement_gsd.md`.

**P0. Transform.** (2 days)
Create `R/transform_pvalues_gsd.R`: exported `transform_pvalues_gsd()`; internal `.gsd_boundary_table(t, spending, alpha, grid_size)` (steps 2 and 3 of 4.1), `.gsd_invert(table, grid, p, floor)` (step 4), `.gsd_check_spending(spending, alpha)` (the $\phi(a,1)=a$ check); the `multigrain_pvals_gsd` constructor, `print()` and `summary()`. Add `gsDesign` to `Imports` and `graphicalMCP (>= 0.3.0)` to `Suggests` in `DESCRIPTION`.
Gate: `tests/testthat/test-transform_pvalues_gsd.R` passes; it asserts Table 1 boundaries and Table 2 repeated p-values of Maurer and Bretz to four significant figures (Appendix A values), that a deliberately non-monotone spending function aborts, that the grid inverse is within $10^{-5}$ of `uniroot` at 200 random p-values, that the first-look short-circuit returns the raw column `identical()`, and agreement with `graphicalMCP::repeated_p()` and `sequential_p()` within $10^{-5}$ on random inputs (`skip_if_not_installed`).

**P1. Simulator.** (1 day)
Create `R/sim_pvals_gsd.R`: exported `simulate_pvalues_gsd()` building the canonical covariance over distinct information levels and copying matured columns.
Gate: `test-sim_pvals_gsd.R`: dimensions and the `info_frac` attribute; empirical within-hypothesis correlation $\sqrt{t_k/t_l}$ and cross-hypothesis correlation $\rho\sqrt{\min/\max}$ within $5\times10^{-3}$ at $N=10^6$ (`skip_on_cran`); matured columns identical.

**P2. Kernel.** (3 days)
Create `src/graph_shortcut_gsd.cpp` with the serial and parallel functions of 4.3; regenerate `R/RcppExports.R`. Create `tests/testthat/helper-gsd_reference.R` containing the direct boundary-on-the-fly R implementation from Appendix B (the `direct()` function and its helpers).
Gate: `test-graph_shortcut_gsd.R`: with $K=1$, `rejected` is `identical()` to `graph_shortcut()` on the same matrix and `time` is 1 wherever rejected; exact agreement with the R reference on $N=400$ trials with `look_back` = FFF, TTT, TFT, FTF (both `rejected` and `time`); parallel output `identical()` to serial at 1, 2, 4 and 8 threads; a hypothesis with weight exactly 0 and $p^r=\gamma_{\min}$ is not rejected; `graphicalMCP::graph_test_shortcut_gsd()` agreement on the Maurer and Bretz case study (no matured hypotheses, `look_back = FALSE`) and on an Example 5 configuration run in the oracle with `look_back = c(TRUE, FALSE)`.

**P3. Gain.** (2 days)
Create `R/trial_success_gsd.R`: exported `trial_success_gsd()`; its own copies of the validator, index counter and transformer with these changes relative to `R/trial_success.R`: the allowed-operator list at `:153` gains `== != < <= > >=` and the names of supplied tables; the symbol rule at `:174` becomes `^[rt]\d+$`; the index counter at `:326` matches `[rt][0-9]+`; `parse_and_transform()` gains a comparison branch returning `double(A op B)` of type `bool`, a `t<i>` branch returning `double(t(i, idx))` of type `real`, a table-call branch returning `name_tab[t(i, idx)]` of type `real`, and `r<i>` becomes `double(t(i, idx) > 0)`; the generated code at `:218` becomes `double f(IntegerMatrix t)` preceded by one `static const double name_tab[] = {0.0, ...}` per table; unary minus is handled by checking `length(parts) == 2` rather than indexing a second argument. Class `c("multigrain_trial_success_gsd", "multigrain_trial_success")`; fields `func`, `m`, `K`, `objective`, `cpp_code`, `tables`.
Gate: `test-trial_success_gsd.R`: the four expressions in 4.6 compile and, on hand-built time matrices, return the value computed in R; `d(0)` is 0; `K` is inferred from tables and checked against an explicit `K`; an unknown symbol, an unsupported operator, and a table of the wrong length each error with the intended message; `cpp_code` snapshots.

**P4. Optimiser and post-processing.** (3 days)
Create `R/objective_function_gsd.R` (`create_obj_func_gsd()`), `R/optimisation_gsd.R` (`graph_optimise_gsd()`, `.graph_optimise_ga_gsd()`, `.graph_optimise_local_gsd()`, `control_prepare_dims()`), `R/calc_power_gsd.R` (`calc_power_pvals_gsd()`, `.eval_custom_power_gsd()`), `R/post_optim_processing_gsd.R` (`prune_graph_gsd()`, `.try_prune_gsd()`). Reuse unchanged: `split_theta()`, `recover_full_weights()`, `recover_full_trans_matrix()`, `param_to_solution()`, `repair_graph()`, `.build_start_matrix()`, `create_start_params()`, `choose_graph()`, `graph_optimal()`.
Gate: `test-optimisation_gsd.R` runs the whole pipeline at $m=3$, $K=2$, $N=2000$ and returns a valid graph with a populated `power` element; the manuscript's Figure 3b is reproduced as a `skip_on_cran` test: optimal $w_{\text{PFS}}$ on a grid at $N=10^5$ for $\delta\in\{1,0.75,0.5\}$ and value ratios in the figure's range, within 0.02 of the reference values (OPEN item 7 names the reference); `test-post_optim_processing_gsd.R` shows pruning never lowers the gain and respects fixed entries; `test-calc_power_gsd.R` checks shapes, that `local_power` equals the last column of `local_power_by_analysis`, that rows of `time_distribution` sum to one, and both kinds of `custom_power`. The no-change gate of 4.10 is run as a script (`dev/gsd_identity_check.R`) against an installed `main` and an installed branch build, and its output is pasted into the report.

**P5. Documentation and release notes.** (1 to 2 days)
roxygen for every new export with a runnable example built on a two-hypothesis, two-look design; a "Group sequential designs" group in `_pkgdown.yml`; an article `vignettes/articles/group-sequential.Rmd` that walks Example 5 end to end; the `NEWS.md` entry of section 9.

Total: 12 to 13 days. Cross-cutting: timings come from installed builds (`devtools::install(quick = TRUE)`), never `load_all()`; `testthat::test_file()` only; no runs above $m=4$, $N=10^4$ outside the two named exceptions.

## 7. Test plan

| File | What it asserts |
|---|---|
| `test-transform_pvalues_gsd.R` | Table 1 and Table 2 reproduction; monotonicity abort; grid accuracy; short-circuit identity; matured copy-forward; `NA` before first look gives 1; `cummin` under `look_back`; spend-at-1 check; information fractions above 1 accepted; oracle agreement with `graphicalMCP` |
| `test-sim_pvals_gsd.R` | Dimensions, attribute, canonical correlations, matured columns |
| `test-graph_shortcut_gsd.R` | $K=1$ identity; R-reference equivalence with mixed look-back; parallel equals serial; zero weight never rejects; floor value never rejects at zero weight; oracle agreement on the case study and on Example 5 (oracle `look_back = c(TRUE, FALSE)`) |
| `test-trial_success_gsd.R` | The four example gains; `d(0)`; `K` inference and mismatch error; symbol and operator errors; snapshots |
| `test-optimisation_gsd.R` | End-to-end small run; Figure 3b (`skip_on_cran`); subsampling keeps a 3-D array; result class and `power` fields; `print()`/`summary()` succeed |
| `test-calc_power_gsd.R` | Output shapes and internal consistency; custom measures of both kinds; evaluating the same graph twice gives identical values |
| `test-post_optim_processing_gsd.R` | Pruning monotone in gain; fixed entries respected |
| `dev/gsd_identity_check.R` | `graph_optimise()` identical on `main` and branch for a fixed seed; `.Random.seed` identical afterwards |

## 8. Adversarial checks for the review agent

1. Example 5 configuration where PFS is below its level only after OS recycles alpha at the final look: the kernel must report $\tau_{\text{PFS}}=2$; the oracle reproduces this only with `look_back = c(TRUE, FALSE)`.
2. $K=1$ with $t=1$: GSD kernel `rejected` `identical()` to `graph_shortcut()`, including rows with $p<10^{-12}$; no `all.equal`, no tolerance.
3. A row with $p=5\times10^{-11}$ at a full-information look and an allocation of $2.5\times10^{-11}$ (weight $10^{-4}$, edge $10^{-5}$): must not reject.
4. Hypothesis weight exactly 0 and $p^r=\gamma_{\min}$: no rejection.
5. A spending function that is not well ordered (for example one that switches family with the level): the transform must abort.
6. `NA` at look 1, data at looks 2 and 3, `look_back = FALSE`: the hypothesis is inactive at look 1 and look 3 uses look-3 evidence only.
7. Data at looks 1 and 3, `NA` at 2, `look_back = TRUE`: the look-1 sequential p-value is carried into look 2.
8. A hypothesis with data at a single interim look ($t<1$) and nothing after: its repeated p-value exceeds its raw p-value (the boundary is the partial spend, not the level).
9. Parallel versus serial kernel `identical()` at 1, 2 and 4 threads with mixed look-back.
10. Gains `(t1 == 1) && r2` and `d(t1)` with `d = c(1, 0.75)`: compile, and `d(0)` is 0.
11. Evaluate the same $(\mathbf w,\mathbf G)$ twice on the same transformed array: identical gain to machine precision.
12. GA and COBYLA subsampling on the 3-D array keeps `dim` (a `drop = TRUE` slip silently turns $N\times m\times 1$ into a matrix).
13. `graph_optimise()` on `main` and on the branch with the same seed: `identical()` objects.
14. Non-`NA` raw p-values after maturity that differ from the maturity p-value: warning, and the maturity value is used.
15. A user spending function that returns cumulative spend as a plain numeric (no `$spend`): accepted.
16. `print()` and `summary()` of a GSD `multigrain_graph_optimal` and of a `multigrain_pvals_gsd` with $m=6$, $K=3$ and `NA` padding.

## 9. Proposed `NEWS.md` wording

```
# multigrain (development version)

## New functionality

* Group sequential designs are supported through a new family of functions.
  `simulate_pvalues_gsd()` simulates p-values at several analyses under the
  canonical joint normal model; `transform_pvalues_gsd()` converts them to
  repeated (or, per hypothesis, sequential) p-values using alpha-spending
  boundaries computed with gsDesign; `trial_success_gsd()` defines a gain
  function of rejections and decision times (symbols `t1, t2, ...`,
  comparisons, and discount tables such as `d(t1)`); `graph_optimise_gsd()`
  optimises a graph for a group sequential design; `calc_power_pvals_gsd()`
  reports power by hypothesis and by analysis together with the decision-time
  distribution. The fixed-sample functions are unchanged.
* New dependency: gsDesign (Imports). graphicalMCP (>= 0.3.0) is used in the
  test suite only (Suggests).
```

No bug-fix entry: nothing in the fixed-sample path changes. The two pre-existing issues noticed during this design (row subsampling permutes the first `nsim` rows rather than sampling from all rows; unary minus in the gain parser) should become GitHub issues.

## 10. Open items

1. **Function names.** `transform_pvalues_gsd()` versus `gsd_pvalues()`; `trial_success_gsd()`; whether `simulate_pvalues_gsd()` takes `power_nominal` (consistent with `simulate_pvalues()`) or noncentrality parameters directly (consistent with the manuscript's tables). Closed by the user before P0.
2. **`grid_size` default.** 1024 (error $1.5\times10^{-6}$) or 512 ($4\times10^{-6}$); both are far below Monte Carlo noise. Closed by choosing; no further evidence needed.
3. **Record location.** `debug/` is git-ignored. Copy this record to `dev/gsd_design_record.md` and commit it before the implementation session, or the implementation agent cannot read it.
4. **Merging the two code paths.** Once `plan-for-sparsity` has landed, fold `create_obj_func_gsd()`, the two optimiser drivers and `.try_prune_gsd()` back into the fixed-sample versions with a kernel-dispatch argument, so that the fixed-sample case is $K=1$ of one implementation. Not part of this roadmap.
5. **Reporting for protocols.** Whether `calc_power_pvals_gsd()` should also return the nominal boundaries of the final graph at each weight in its closure (as `graphicalMCP` does with `verbose = TRUE`). The tables on the `multigrain_pvals_gsd` object make this cheap; needs a user decision on scope.
6. **`K` on the gain object.** Whether `K` must be supplied when no discount table fixes it, or may be left unset until first use. The record allows unset; a stricter rule is a one-line change.
7. **Reference values for the Figure 3b gate.** The manuscript reports the curves graphically. The article's reproduction scripts (supplied with the paper as a data file) should provide the numbers; if not, the gate compares against the grid run in this record's own reproduction, to be added to Appendix B before P4 starts.
8. **Per-look local power constraints in pruning.** The fixed-sample `power_constraint` is per hypothesis. Under a GSD a constraint could be per (hypothesis, look). Default: per hypothesis over any look; per-look constraints deferred.

---

# Appendix A. Reproduction of Maurer and Bretz (2013) Tables 1 and 2

Script run on 2026-09-15 with R 4.6.1 and gsDesign 3.11.0. Table 1 lists nominal boundaries $\alpha^*_{k}(\gamma)$ for the O'Brien–Fleming-type spending function at information fractions 1/3, 2/3, 1; Table 2 lists the repeated p-values computed with ADDPLAN for the case-study p-values. The paper reports four significant figures.

```r
# check1_tables.R -- Does inverting gsDesign boundaries reproduce Maurer & Bretz
# (2013) Tables 1 and 2, and does the fixed-sample cascade on repeated p-values
# reproduce their case-study rejections?
suppressPackageStartupMessages(library(gsDesign))
options(digits = 6)

# Forward map: nominal p-value boundaries alpha*_k(gamma) for a one-sided
# alpha-spending design at information fractions t, spending function sf.
nom_bounds <- function(gamma, t, sf = sfLDOF) {
  K <- length(t)
  inc <- diff(c(0, sf(gamma, t)$spend))                     # per-analysis spend
  b <- gsBound1(theta = 0, I = t, a = rep(-20, K), probhi = inc)$b
  pnorm(b, lower.tail = FALSE)
}

t3 <- c(1/3, 2/3, 1)
cat("== Maurer-Bretz Table 1 (LDOF, t = 1/3, 2/3, 1) ==\n")
cat("gamma=0.0125  :", format(nom_bounds(0.0125, t3), digits = 4),  "  paper: 0.00002 0.0022 -\n")
cat("gamma=0.01875 :", format(nom_bounds(0.01875, t3), digits = 4), "  paper k=2: 0.004\n")
cat("gamma=0.00625 :", format(nom_bounds(0.00625, t3), digits = 4), "  paper k=2: 0.0008\n")
cat("gamma=0.025   :", format(nom_bounds(0.025, t3), digits = 4),   "  paper k=2: 0.006, k=3: 0.02313\n")
d <- gsDesign(k = 3, test.type = 1, alpha = 0.0125, timing = t3, sfu = sfLDOF)
cat("cross-check gsDesign():", format(pnorm(d$upper$bound, lower.tail = FALSE), digits = 4), "\n")

# Inverse: repeated p-value = the gamma solving alpha*_k(gamma) = p
rep_p <- function(p, k, t, sf = sfLDOF) {
  uniroot(function(g) nom_bounds(g, t, sf)[k] - p, c(1e-12, 0.999), tol = 1e-10)$root
}
cat("\n== Maurer-Bretz Table 2 (repeated p-values, ADDPLAN) ==\n")
p1 <- c(0.0062, 0.017, 0.009, 0.13); p2 <- c(0.0002, 0.0035, 0.002, 0.06)
r1 <- sapply(p1, rep_p, k = 1, t = t3); r2 <- sapply(p2, rep_p, k = 2, t = t3)
cat("k=1:", format(r1, digits = 4), "  paper: 0.1141 0.1683 0.1316 0.382\n")
cat("k=2:", format(r2, digits = 4), "  paper: 0.0024 0.0172 0.0117 0.1285\n")

cat("\n== Fixed-sample cascade on the k=2 repeated p-values (Figure 1 graph) ==\n")
w <- c(0.5, 0.5, 0, 0)
G <- rbind(c(0, .5, .5, 0), c(.5, 0, 0, .5), c(0, 1, 0, 0), c(1, 0, 0, 0))
alpha <- 0.025; p <- r2; rej <- rep(FALSE, 4); order <- integer(0)
repeat {
  j <- which(!rej & p <= w * alpha)[1]; if (is.na(j)) break
  rej[j] <- TRUE; order <- c(order, j)
  wn <- w + w[j] * G[j, ]; wn[j] <- 0
  Gn <- G
  for (l in 1:4) for (k in 1:4)
    Gn[l, k] <- if (l == k || rej[l] || rej[k]) 0 else (G[l, k] + G[l, j] * G[j, k]) / (1 - G[l, j] * G[j, l])
  w <- wn; G <- Gn
}
cat("rejected:", rej, " order:", order, "  paper: H1, H2, H3 rejected; H4 retained\n")
```

Output:

```
== Maurer-Bretz Table 1 (LDOF, t = 1/3, 2/3, 1) ==
gamma=0.0125  : 1.517e-05 2.215e-03 1.180e-02   paper: 0.00002 0.0022 -
gamma=0.01875 : 4.679e-05 3.976e-03 1.750e-02   paper k=2: 0.004
gamma=0.00625 : 2.179e-06 8.105e-04 5.986e-03   paper k=2: 0.0008
gamma=0.025   : 0.0001035 0.0060122 0.0231281   paper k=2: 0.006, k=3: 0.02313
cross-check gsDesign(): 1.517e-05 2.215e-03 1.180e-02 

== Maurer-Bretz Table 2 (repeated p-values, ADDPLAN) ==
k=1: 0.1141 0.1682 0.1315 0.3820   paper: 0.1141 0.1683 0.1316 0.382
k=2: 0.002393 0.017161 0.011649 0.128460   paper: 0.0024 0.0172 0.0117 0.1285

== Fixed-sample cascade on the k=2 repeated p-values (Figure 1 graph) ==
rejected: TRUE TRUE TRUE FALSE  order: 1 2 3   paper: H1, H2, H3 rejected; H4 retained
```

# Appendix B. End-to-end equivalence, edge cases, clamp safety and cost

Script run on 2026-09-15 with R 4.6.1 and gsDesign 3.11.0. The `direct()` function is the Maurer and Bretz procedure with boundaries computed at the current weight on the fly (memoised); `transform()` and `kern()` are the design of sections 4.1 and 4.3 written in R. Three hypotheses, three looks: $H_1$ with LDOF spending at $(1/3, 2/3, 1)$, $H_2$ with HSD($-2$) spending maturing at look 2, $H_3$ mature at look 1. The four `look_back` configurations mix semantics across hypotheses. The clamp-safety block reproduces the red-team case of section 4.1.

```r
# check2_equivalence.R -- Is "transform once + K passes of the fixed-sample
# cascade" identical to a direct boundary-on-the-fly implementation of the
# Maurer-Bretz procedure, for mixed per-hypothesis look-back and matured
# hypotheses?  Also: gsBound1 edge cases, grid accuracy, clamp safety, cost.
suppressPackageStartupMessages(library(gsDesign)); options(digits = 6, warn = -1)
alpha <- 0.025; g_min <- 1e-14
nom_bounds <- function(gamma, t, sf = sfLDOF) {
  K <- length(t); inc <- diff(c(0, sf(gamma, t)$spend))
  b <- gsBound1(theta = 0, I = t, a = rep(-20, K), probhi = inc)$b
  pnorm(b, lower.tail = FALSE)
}

cat("== gsBound1 edge cases ==\n")
cat("duplicate t=c(0.7,1,1):", format(nom_bounds(0.025, c(0.7, 1, 1)), digits = 4), "  <- third look unrejectable: never pass repeated t=1\n")
cat("t > 1, t=c(0.5,1.2):   ", format(nom_bounds(0.025, c(0.5, 1.2)), digits = 4), "\n")
cat("K=1, t=1, gamma=0.0125:", format(nom_bounds(0.0125, 1), digits = 17), " (identical to gamma?", identical(nom_bounds(0.0125, 1), 0.0125), ")\n")
for (g in c(1e-4, 1e-8, 1e-14)) cat(sprintf("gamma=%g, t=(1/3,2/3,1): %s\n", g, paste(format(nom_bounds(g, c(1/3, 2/3, 1)), digits = 3), collapse = " ")))

cat("\n== Grid size vs inverse accuracy (t=(0.7,1), levels in [1e-14, alpha], log-log linear) ==\n")
t2 <- c(0.7, 1)
rep_p <- function(p, k, t) uniroot(function(g) nom_bounds(g, t)[k] - p, c(1e-15, 0.999), tol = 1e-13)$root
set.seed(2); ptest <- 10^runif(200, -9, log10(0.03))
ex <- lapply(1:2, function(k) sapply(ptest, rep_p, k = k, t = t2))
for (Gn in c(256, 512, 1024)) {
  grid <- exp(seq(log(g_min), log(alpha), length.out = Gn)); fwd <- do.call(rbind, lapply(grid, nom_bounds, t = t2))
  for (k in 1:2) { f <- fwd[, k]; ok <- !duplicated(f) & f > 0
    ap <- exp(approx(log(f[ok]), log(grid[ok]), xout = log(ptest), rule = 1)$y); sel <- ex[[k]] <= alpha & !is.na(ap)
    cat(sprintf("G=%4d k=%d: n=%3d  max rel err %.2e\n", Gn, k, sum(sel), max(abs(ex[[k]] - ap)[sel] / ex[[k]][sel]))) } }

cat("\n== End-to-end equivalence: direct vs transform-once + cascade (m=3, K=3, N=400) ==\n")
set.seed(11); N <- 400; m <- 3; K <- 3
tmat <- rbind(c(1/3, 2/3, 1), c(0.5, 1, 1), c(1, 1, 1))     # H2 matures at k=2, H3 at k=1
sfs  <- list(sfLDOF, function(a, t) sfHSD(a, t, param = -2), sfLDOF)
w0 <- c(0.5, 0.3, 0.2); G0 <- rbind(c(0, .5, .5), c(.5, 0, .5), c(.5, .5, 0)); Delta <- c(2.6, 2.4, 2.9)
mature_k <- apply(tmat, 1, function(t) which(t >= 1)[1]); eff_l <- function(i, k) min(k, mature_k[i])
cache <- new.env()
bound <- function(i, l, a) { key <- sprintf("%d_%d_%.15g", i, l, a); v <- cache[[key]]
  if (is.null(v)) { v <- nom_bounds(a, tmat[i, 1:l], sfs[[i]])[l]; assign(key, v, envir = cache) }; v }
P <- lapply(1:m, function(i) { kk <- mature_k[i]; t <- tmat[i, 1:kk]
  S <- outer(t, t, function(a, b) sqrt(pmin(a, b) / pmax(a, b)))
  Z <- mvtnorm::rmvnorm(N, mean = Delta[i] * sqrt(t), sigma = S); p <- pnorm(Z, lower.tail = FALSE)
  cbind(p, matrix(p[, kk], N, K - kk)) })                    # matured p carried forward
upd <- function(w, G, i, rej) { wn <- w + w[i] * G[i, ]; wn[i] <- 0; Gn <- matrix(0, m, m)
  for (l in 1:m) for (j in 1:m) if (l != j && !rej[l] && !rej[j]) { d <- 1 - G[l, i] * G[i, l]; Gn[l, j] <- if (d > 0) (G[l, j] + G[l, i] * G[i, j]) / d else 0 }
  list(w = wn, G = Gn) }
# Direct Maurer-Bretz procedure: boundaries computed on the fly at the current weight
direct <- function(pk, lb) { w <- w0; G <- G0; rej <- rep(FALSE, m); tm <- integer(m)
  for (k in 1:K) repeat { hit <- NA
    for (i in which(!rej)) { a <- w[i] * alpha; if (a <= 0) next
      ls <- if (lb[i]) 1:eff_l(i, k) else eff_l(i, k)
      if (any(pk[i, ls] < sapply(ls, function(l) bound(i, l, a)))) { hit <- i; break } }
    if (is.na(hit)) break; rej[hit] <- TRUE; tm[hit] <- k; u <- upd(w, G, hit, rej); w <- u$w; G <- u$G }
  c(rej, tm) }
# Transform: grid tables per hypothesis, inverse with floor clamp and first-look short-circuit
grids <- lapply(1:m, function(i) { g <- exp(seq(log(g_min), log(alpha), length.out = 1024)); kk <- mature_k[i]
  list(g = g, fwd = do.call(rbind, lapply(g, function(gm) nom_bounds(gm, tmat[i, 1:kk], sfs[[i]])))) })
inv <- function(i, l, p) {
  if (l == 1 && tmat[i, 1] >= 1) return(p)                    # short-circuit: boundary is the level itself
  g <- grids[[i]]$g; f <- grids[[i]]$fwd[, l]; ok <- !duplicated(f) & f > 0
  out <- exp(approx(log(f[ok]), log(g[ok]), xout = log(p), rule = 1)$y)
  out[p >= max(f)] <- 1; out[p < min(f[ok])] <- g_min; out }
transform <- function(lb) { A <- array(NA_real_, c(N, m, K))
  for (i in 1:m) for (k in 1:K) { l <- eff_l(i, k); A[, i, k] <- inv(i, l, P[[i]][, l]) }
  for (i in which(lb)) A[, i, ] <- t(apply(A[, i, ], 1, cummin)); A }
kern <- function(pr) { w <- w0; G <- G0; rej <- rep(FALSE, m); tm <- integer(m)
  for (k in 1:K) repeat { i <- which(!rej & pr[, k] < w * alpha)[1]; if (is.na(i)) break
    rej[i] <- TRUE; tm[i] <- k; u <- upd(w, G, i, rej); w <- u$w; G <- u$G }
  c(rej, tm) }
for (lb in list(c(F, F, F), c(T, T, T), c(T, F, T), c(F, T, F))) {
  A <- transform(lb)
  D  <- t(sapply(1:N, function(n) direct(rbind(P[[1]][n, ], P[[2]][n, ], P[[3]][n, ]), lb)))
  Tr <- t(sapply(1:N, function(n) kern(A[n, , ])))
  cat(sprintf("look_back=%s: rejection mismatches %d/%d, time mismatches %d/%d; local power %s; mean decision time %s\n",
    paste(substr(as.character(lb), 1, 1), collapse = ""), sum(D[, 1:m] != Tr[, 1:m]), N * m, sum(D[, -(1:m)] != Tr[, -(1:m)]), N * m,
    paste(format(colMeans(Tr[, 1:m]), digits = 3), collapse = " "),
    paste(format(colMeans(replace(Tr[, -(1:m)], Tr[, -(1:m)] == 0, NA), na.rm = TRUE), digits = 3), collapse = " "))) }
cat("matured-at-first-look hypothesis: transformed column identical() to raw p:", identical(transform(c(F, F, F))[, 3, 1], P[[3]][, 1]), "\n")

cat("\n== Clamp safety (red-team check): tiny allocation 1e-4 * 1e-5 * alpha after the kernel's snaps ==\n")
a_tiny <- 1e-4 * 1e-5 * alpha; p_row <- 5e-11                 # true p^r at a t=1 look ~ p = 5e-11 > a_tiny
old_floor <- 1e-10; g_old <- exp(seq(log(old_floor), log(alpha), length.out = 1024)); f_old <- sapply(g_old, function(gm) nom_bounds(gm, t2)[2])
pr_old <- if (p_row < min(f_old)) 0 else NA                    # old rule: below table -> 0
f_new <- grids[[1]]$fwd[, 3]; ok <- !duplicated(f_new) & f_new > 0
pr_new <- exp(approx(log(f_new[ok]), log(grids[[1]]$g[ok]), xout = log(p_row), rule = 1)$y)
cat(sprintf("allocation %.2e, p = %.0e: old rule p^r = %g -> rejects (wrong: %s); new rule p^r = %.2e -> rejects %s\n",
  a_tiny, p_row, pr_old, pr_old < a_tiny, pr_new, pr_new < a_tiny))
set.seed(3); cat(sprintf("fraction of p < 1e-10 at NCP 4: %.4f\n", mean(pnorm(rnorm(1e6, 4), lower.tail = FALSE) < 1e-10)))

cat("\n== Cost ==\n")
tm <- system.time(fw <- do.call(rbind, lapply(exp(seq(log(g_min), log(alpha), length.out = 1024)), nom_bounds, t = c(1/3, 2/3, 1))))
cat(sprintf("boundary table, K=3, 1024 levels: %.2f s (%.2f ms per gsBound1 call)\n", tm[["elapsed"]], 1000 * tm[["elapsed"]] / 1024))
g <- grids[[1]]; p <- runif(1e6); ok <- !duplicated(g$fwd[, 3])
tm <- system.time(for (j in 1:12) exp(approx(log(g$fwd[ok, 3]), log(g$g[ok]), xout = log(p), rule = 2)$y))
cat(sprintf("inverse interpolation, 12 slices x 1e6 values: %.2f s\n", tm[["elapsed"]]))
```

Output:

```
== gsBound1 edge cases ==
duplicate t=c(0.7,1,1): 7.384e-03 2.275e-02 2.754e-89   <- third look unrejectable: never pass repeated t=1
t > 1, t=c(0.5,1.2):    0.001525 0.024333 
K=1, t=1, gamma=0.0125: 0.012499999999999968  (identical to gamma? FALSE )
gamma=0.0001, t=(1/3,2/3,1): 1.60e-11 1.89e-06 9.93e-05
gamma=1e-08, t=(1/3,2/3,1): 2.75e-89 2.24e-12 1.00e-08
gamma=1e-14, t=(1/3,2/3,1): 2.75e-89 2.75e-89 9.99e-15

== Grid size vs inverse accuracy (t=(0.7,1), levels in [1e-14, alpha], log-log linear) ==
G= 256 k=1: n=178  max rel err 1.04e-05
G= 256 k=2: n=199  max rel err 1.70e-05
G= 512 k=1: n=178  max rel err 2.58e-06
G= 512 k=2: n=199  max rel err 4.01e-06
G=1024 k=1: n=178  max rel err 6.72e-07
G=1024 k=2: n=199  max rel err 1.48e-06

== End-to-end equivalence: direct vs transform-once + cascade (m=3, K=3, N=400) ==
look_back=FFF: rejection mismatches 0/1200, time mismatches 0/1200; local power 0.690 0.568 0.730; mean decision time 2.58 1.87 1.15
look_back=TTT: rejection mismatches 0/1200, time mismatches 0/1200; local power 0.69 0.57 0.73; mean decision time 2.58 1.87 1.15
look_back=TFT: rejection mismatches 0/1200, time mismatches 0/1200; local power 0.690 0.568 0.730; mean decision time 2.58 1.87 1.15
look_back=FTF: rejection mismatches 0/1200, time mismatches 0/1200; local power 0.69 0.57 0.73; mean decision time 2.58 1.87 1.15
matured-at-first-look hypothesis: transformed column identical() to raw p: TRUE 

== Clamp safety (red-team check): tiny allocation 1e-4 * 1e-5 * alpha after the kernel's snaps ==
allocation 2.50e-11, p = 5e-11: old rule p^r = 0 -> rejects (wrong: TRUE); new rule p^r = 5.00e-11 -> rejects FALSE
fraction of p < 1e-10 at NCP 4: 0.0092

== Cost ==
boundary table, K=3, 1024 levels: 2.58 s (2.52 ms per gsBound1 call)
inverse interpolation, 12 slices x 1e6 values: 0.89 s
```

# Appendix C. Red-team findings and their disposition

An independent review of the transform-once architecture was run against the current tree before this record was finalised. Findings and what was done with them:

1. Repeated information fractions of 1 passed to `gsBound1` yield an unrejectable boundary at later looks (Appendix B, first block). Adopted as the matured copy-forward rule, 4.1 step 5.
2. `gsBound1` at a single look with $t=1$ returns the level with about $10^{-16}$ error, so $K=1$ would not be bit-identical through the table. Adopted as the short-circuit, 4.1 step 1.
3. A floor of 0 below the table rejects rows with $p<10^{-10}$ at any positive allocation, and allocations of $2.5\times10^{-11}$ exist after the kernel's snaps. Adopted: floor at $\gamma_{\min}=10^{-14}$, 4.1 step 4, verified in Appendix B.
4. The transform assumes well-ordered boundaries. Adopted as the monotonicity assertion, 4.1 step 3.
5. `graphicalMCP` with `NA` padding and `look_back = FALSE` does not reproduce Example 5. Adopted as the oracle note in 4.9 and in the test plan.
6. The list of places where the fixed-sample code assumes an $N\times m$ double matrix (section 3) and the parser sites that need changing (section 6, P3) were taken from the review and verified line by line.
7. Two pre-existing issues outside this work (row subsampling permutes the first `nsim` rows; unary minus in the gain parser) are recorded in section 9 for GitHub issues.

Confirmed to hold, no change needed: mixed per-hypothesis look-back with recycling; zero-initial-weight hypotheses receiving alpha later; information fractions above 1; decision-time semantics under both look-back settings.
