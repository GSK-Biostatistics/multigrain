# How multigrain Applies a Multiple-Testing Graph

*Developer notes on `graph_shortcut()` and `graph_shortcut_parallel()` (`src/graph_shortcut.cpp`). Written for a statistician who is comfortable with simulation and power, but new to multiple testing. Last updated 22 October 2026, branch `fast-shortcut`.*

**Contents**

1. [What the function does](#1-what-the-function-does)
2. [Background: multiple testing with graphs](#2-background-multiple-testing-with-graphs)
3. [A worked example](#3-a-worked-example)
4. [The sequential kernel](#4-the-sequential-kernel)
5. [The table kernel](#5-the-table-kernel)
6. [Choosing a kernel and running in parallel](#6-choosing-a-kernel-and-running-in-parallel)
7. [Floating point: keeping the two kernels identical](#7-floating-point-keeping-the-two-kernels-identical)
8. [Testing and working on the code](#8-testing-and-working-on-the-code)
9. [Glossary and references](#9-glossary-and-references)

---

## 1. What the Function Does

multigrain designs multiple-testing strategies for clinical trials and estimates their power by simulation. The simulation produces a matrix of p-values, one row per simulated trial and one column per hypothesis. `graph_shortcut()` takes that matrix and a testing strategy, and returns which hypotheses each simulated trial rejects:

```r
rej <- graph_shortcut(pvals, alpha, w, G)  # internal; call as multigrain:::graph_shortcut
# pvals: N x m numeric matrix (N simulated trials, m hypotheses)
# alpha: overall significance level, e.g. 0.025
# w:     length-m vector of initial weights
# G:     m x m transition matrix
# rej:   N x m logical matrix, TRUE where the trial rejects that hypothesis
```

Everything else is computed from `rej`:

* **Power per hypothesis:** `colMeans(rej)`.
* **Disjunctive power:** the probability of rejecting at least one hypothesis.
* **Conjunctive power:** the probability of rejecting all of them.
* **User-defined success criteria,** built with `trial_success()`.

It is called in three places:

| Caller | File | Purpose |
|---|---|---|
| `calc_power_pvals()` | `R/calc_power.R` | Power of a given graph |
| The objective closure from `create_obj_func()` | `R/objective_function.R` | Scores each candidate graph during optimisation, many thousands of times per run |
| `graph_optimise()` and its helpers | `R/optimisation.R` | Re-scores the chosen graph on the full simulated sample |

The optimiser calls it so often that its speed matters, and much of this note explains how it is made fast without changing its answers.

`graph_shortcut_parallel()` does the same job across several threads. Its result is identical whatever the number of threads.

---

## 2. Background: Multiple Testing with Graphs

### 2.1 The Problem

A trial often tests several hypotheses: two doses, or a primary and secondary endpoint. If each is tested at $\alpha = 0.025$, the chance of at least one false rejection is larger than 0.025. Regulators require control of the **familywise error rate** (FWER): the probability of rejecting at least one true null hypothesis, whichever hypotheses happen to be true, must be at most $\alpha$.

The simplest fix is **Bonferroni**: test each of $m$ hypotheses at $\alpha/m$. **Weighted Bonferroni** splits $\alpha$ unevenly. Hypothesis $H_j$ is tested at $\alpha w_j$, with weights $w_j \ge 0$ and $\sum_j w_j \le 1$. Both control the FWER, but they waste power: once $H_1$ is rejected, its share of $\alpha$ is never used again.

### 2.2 Graphs Recycle Alpha

The graphical approach (Bretz et al. 2009) recycles that share. A strategy is a directed graph:

* **Nodes** are the hypotheses $H_1, \dots, H_m$. Each starts with a weight $w_j$, so it is tested at level $\alpha w_j$.
* **Edges** carry transition weights $g_{jk}$, the share of $H_j$'s level passed to $H_k$ if $H_j$ is rejected. These form the $m \times m$ **transition matrix** $\mathbf{G}$.
* **Rules for a valid graph:** the diagonal is zero, every entry is non-negative, and each row sums to at most 1.

Testing proceeds as follows (Algorithm 1 in Bretz et al. 2009):

1. Find any unrejected $H_j$ with $p_j < \alpha w_j$. If there is none, stop.
2. Reject $H_j$.
3. Pass its weight along its edges. For every other unrejected $H_k$:

$$w_k \leftarrow w_k + w_j\, g_{jk}.$$

4. Remove $H_j$ and update the edges between the remaining nodes. For every pair $k \ne l$ of remaining nodes:

$$g_{kl} \leftarrow \frac{g_{kl} + g_{kj}\, g_{jl}}{1 - g_{kj}\, g_{jk}}.$$

   The numerator adds the indirect route $k \to j \to l$ to the direct route $k \to l$. The denominator removes the weight that would otherwise bounce back and forth between $k$ and $j$. If the denominator is 0, $H_k$ and $H_j$ pass everything to each other, and row $k$ is set to zero.
5. Go back to step 1.

This procedure controls the FWER. It is a shortcut for a **closed test** of weighted Bonferroni tests, and Section 2.3 explains why that matters for the code.

### 2.3 The One Fact the Code Relies On

Call $S$ the set of hypotheses not yet rejected. A key result (Bretz et al. 2009) is that **the weights after the updates depend only on $S$, not on the order in which the other hypotheses were rejected.** So each set $S$ has well-defined **local weights** $w_j(S)$ and **local critical values**

$$c_j(S) = \alpha\, w_j(S), \qquad j \in S.$$

The final set of rejections is also the same whichever qualifying hypothesis is rejected first.

Two consequences shape the implementation:

* **Order does not matter,** so the kernel can always pick the first qualifying hypothesis by index (Section 4), or reject all qualifying hypotheses at once (Section 5).
* **Weights depend only on $S$,** so the critical values for every possible $S$ can be computed once per graph and looked up, instead of recomputed in every trial (Section 5).

One more property is used repeatedly: **local weights only grow as $S$ shrinks.** If $S' \subseteq S$, then $w_j(S') \ge w_j(S)$. Removing a hypothesis only ever passes weight *to* the others. The same holds for the transition weights $g_{kl}$, because the update above has a non-negative numerator and a denominator at most 1. This is why a hypothesis that can be rejected stays rejectable as testing proceeds.

---

## 3. A Worked Example

The example has three hypotheses and $\alpha = 0.025$. $H_1$ and $H_2$ share the initial alpha. $H_3$ starts with nothing, so it is tested only if alpha reaches it. Each node passes half its weight to each of the others.

$$\mathbf{w} = (0.5,\ 0.5,\ 0), \qquad \mathbf{G} = \begin{pmatrix} 0 & 0.5 & 0.5 \\ 0.5 & 0 & 0.5 \\ 0.5 & 0.5 & 0 \end{pmatrix}$$

**Trial with $\mathbf{p} = (0.010,\ 0.018,\ 0.015)$:**

| Step | Unrejected $S$ | Critical values $c_j(S)$ | Decision |
|---|---|---|---|
| 1 | $\{H_1, H_2, H_3\}$ | $(0.0125,\ 0.0125,\ 0)$ | $p_1 = 0.010 < 0.0125$: reject $H_1$ |
| 2 | $\{H_2, H_3\}$ | $(-,\ 0.01875,\ 0.00625)$ | $p_2 = 0.018 < 0.01875$: reject $H_2$ |
| 3 | $\{H_3\}$ | $(-,\ -,\ 0.025)$ | $p_3 = 0.015 < 0.025$: reject $H_3$ |

The step 2 numbers come from the update rules:

* **Weights after removing $H_1$:**
  - $w_2 = 0.5 + 0.5 \times 0.5 = 0.75$;
  - $w_3 = 0 + 0.5 \times 0.5 = 0.25$.
* **Graph after removing $H_1$:** $g_{23} = (0.5 + 0.5 \times 0.5)/(1 - 0.5 \times 0.5) = 1$, and by symmetry $g_{32} = 1$.

Removing $H_2$ then gives $w_3 = 0.25 + 0.75 \times 1 = 1$. All three hypotheses are rejected, although $p_2 = 0.018$ and $p_3 = 0.015$ were both too large for their starting levels.

With $\mathbf{p} = (0.010,\ 0.020,\ 0.015)$, step 2 fails ($0.020 \ge 0.01875$), and only $H_1$ is rejected.

### 3.1 All Local Critical Values for This Graph

These are exactly what `graph_shortcut_table(0.025, w, G)` returns. Row `mask + 1` holds set $S$, where bit $j-1$ of `mask` is set when $H_j \in S$.

| mask | $S$ | $c_1$ | $c_2$ | $c_3$ |
|---:|---|---:|---:|---:|
| 0 | $\emptyset$ | 0 | 0 | 0 |
| 1 | $\{H_1\}$ | 0.025 | 0 | 0 |
| 2 | $\{H_2\}$ | 0 | 0.025 | 0 |
| 3 | $\{H_1, H_2\}$ | 0.0125 | 0.0125 | 0 |
| 4 | $\{H_3\}$ | 0 | 0 | 0.025 |
| 5 | $\{H_1, H_3\}$ | 0.01875 | 0 | 0.00625 |
| 6 | $\{H_2, H_3\}$ | 0 | 0.01875 | 0.00625 |
| 7 | $\{H_1, H_2, H_3\}$ | 0.0125 | 0.0125 | 0 |

Three things are worth noticing:

* **Hypotheses outside $S$ have level 0.**
* **Levels only grow** as you move to smaller sets.
* **Some sets are "unreachable" in this example.** $\{H_1, H_2\}$ would need $H_3$ rejected first, and $H_3$ starts at level 0. The table is still filled in for every set.

---

## 4. The Sequential Kernel

The sequential kernel is the direct implementation of Section 2.2, and the reference that everything else must match. It lives in `sequential_trial()`. That one function serves the serial path, the parallel worker (`SequentialShortcutWorker`), and the table kernel's fallback (Section 7).

For each trial (one row of `pvals`), it does the following:

1. **Copy the starting state** into scratch buffers:
   - the levels `cur_a[j]` $= \alpha w_j$, precomputed once per call as `init_a`;
   - the graph `cur_G`;
   - the trial's p-values `cur_p`.
2. **Scan** $j = 1, \dots, m$ and take the **first** $j$ with `cur_p[j] < cur_a[j]`. If there is none, the trial is finished.
3. **Record** the rejection in the output, and stop if every hypothesis is now rejected.
4. **Pass on the level:** `cur_a[k] += cur_a[j] * G[j, k]` for every $k$.
5. **Update the graph** with the formula in Section 2.2, writing into a second buffer `new_G` and swapping the two. The update is skipped when only one hypothesis remains, because no further transitions can be used.
6. **Zero row and column $j$** of the graph and set `cur_a[j] = 0`, then go back to step 2.

Some implementation details matter for correctness:

* **The code works with levels, not weights.** It stores levels $\alpha w_j$ ("alpha space") rather than weights, so step 2 is a single comparison. The update in step 4 is the same rule multiplied through by $\alpha$.
* **Strict `<`.** A p-value equal to its level is not rejected. In particular, a hypothesis with level 0 can never be rejected, even if $p = 0$.
* **Degenerate pairs.** If $1 - g_{kj}g_{jk} = 0$, row $k$ of the new graph is set to zero instead of dividing by zero.
* **Layout.** Matrices are column-major (R's layout): entry $(i, j)$ of an $m \times m$ matrix is `G[i + j * m]`. Each trial's p-values are copied into a small dense buffer first, because reading along a row of an $N \times m$ column-major matrix jumps through memory.

**Cost:** $O(m^2)$ to copy the graph, plus $O(m^2)$ per rejection for the graph update. A trial that rejects $r$ hypotheses costs about $O(m^2 (1 + r))$.

---

## 5. The Table Kernel

### 5.1 Idea

By Section 2.3, the critical values depend only on the unrejected set $S$. There are $2^m$ possible sets, and for the $m$ values multigrain usually works with (up to about 12), $2^m$ is far smaller than the number of simulated trials. So:

1. **Once per call,** compute $c_j(S)$ for all $2^m$ sets and store them in a table.
2. **For each trial,** walk through the table instead of updating a graph.

This follows `graph_test_shortcut_fast()` in the graphicalMCP package (Xi and Brockmann, 2024), but is an independent C++ implementation.

### 5.2 Storing Sets as Bitmasks

A set $S$ is stored as an integer whose bit $j-1$ is 1 when $H_j \in S$, as in the table in Section 3.1. The table is a flat array of $2^m \times m$ doubles, with row `S` holding $c_1(S), \dots, c_m(S)$. At the maximum $m = 20$ that is $2^{20} \times 20$ doubles, or 160 MB. That is one reason for the hard limit `TABLE_HARD_MAX_M = 20`; the other is that it keeps the 32-bit mask valid.

### 5.3 Building the Table (`CritTableBuilder`)

The builder visits every set once, depth-first, starting from the full set.

* **Removal order.** It removes hypotheses in increasing index order, so each set is reached from exactly one parent. From a set whose last removed index is `last`, it only removes indices above `last`.
* **Same arithmetic.** Each step applies the same level and graph updates as the sequential kernel (Section 4).
* **Memory.** It keeps one copy of the levels and the graph per depth, so the working memory is $O(m^3)$ rather than one graph per set.
* **Skipped work.** The graph update is skipped for sets that will not be expanded further, because their graph is never used. The exception is three-member sets, whose graph is needed for the check in Section 7.2.

The total cost is $O(2^m m^2)$ per call, independent of $N$.

### 5.4 Walking the Table (`walk_crit_table()`)

For each trial:

1. Start at the full set: `S = (1 << m) - 1`.
2. Look up row `S`, and mark **every** $j \in S$ with $p_j < c_j(S)$.
3. If none qualify, stop. Otherwise remove them all from `S` and go back to step 2.
4. The rejected hypotheses are the bits of the full set that are no longer in `S`.

Rejecting all qualifying hypotheses at once is safe because levels only grow as $S$ shrinks. A hypothesis that qualifies now still qualifies after the others are removed, so the final set is the same as with one-at-a-time rejection. Each step is $O(m)$, with no graph copy and no division.

In the example of Section 3, take $\mathbf{p} = (0.010,\ 0.012,\ 0.030)$:

* **At the full set,** both $H_1$ and $H_2$ qualify ($0.010$ and $0.012 < 0.0125$), so the walk removes both in one step and moves to mask 4.
* **At $\{H_3\}$,** $p_3 = 0.030 \ge 0.025$, so the walk stops.

The sequential kernel reaches the same answer in two steps.

---

## 6. Choosing a Kernel and Running in Parallel

### 6.1 Which Kernel Runs

The table is worth building only when there are enough trials to pay for it. Under the default `"auto"`, the table kernel is used when both of these hold:

* $m \le 20$ (or a lower cap set by an option);
* $2^m \le \kappa N$, where $\kappa$ is `TABLE_KAPPA` $= 2$.

The value of $\kappa$ came from serial benchmarks. Otherwise the sequential kernel runs. The rule deliberately ignores the number of threads.

Two internal, undocumented options override this. They are read once per call:

| Option | Values | Effect |
|---|---|---|
| `multigrain.shortcut_kernel` | `"auto"` (default), `"sequential"`, `"table"` | Forces a kernel. Use this for tests and benchmarks |
| `multigrain.shortcut_table_max_m` | whole number in 0–20 (default 20) | Caps $m$ for the table under `"auto"` |

`multigrain:::graph_shortcut_kernel(n_trials, m)` reports which kernel `"auto"` would choose.

Since both kernels give identical output (Section 7), the choice only ever affects speed.

### 6.2 Parallel Version

`graph_shortcut_parallel()` splits the trials into chunks with `RcppParallel::parallelFor`:

* **Read-only shared inputs.** The p-values, $\mathbf{G}$, the starting levels and (for the table kernel) the table are shared by all threads. The table is built once on the main thread before the parallel loop.
* **No write conflicts.** Each trial writes only its own row of the output.
* **Per-chunk scratch.** Each chunk allocates its own scratch buffers.
* **No R API calls inside workers.** Raw pointers to R's memory are taken on the main thread first, which is the rule for RcppParallel code.
* **Grain size.** The chunk size is tuned automatically from a target amount of work per chunk; `grain_size` overrides it.

---

## 7. Floating Point: Keeping the Two Kernels Identical

### 7.1 The Problem

In exact arithmetic the two kernels give the same answer, by Section 2.3. In floating point they need not:

* The table builder computes $c_j(S)$ along one removal order (increasing index).
* The sequential kernel computes it along the order in which that trial rejected hypotheses.

The two results can differ in the last bit or two. That only matters when a p-value lies within that rounding error of a critical value, which is rare with continuous simulated p-values but common with rounded ones.

Here is an example found in code review, with $\alpha = 0.025$, $\mathbf{w} = (0.1, 0.4, 0.5)$, every off-diagonal $g = 0.5$, and $\mathbf{p} = (0.005, 0, 0.025)$:

* The sequential kernel computes $c_3(\{H_3\})$ as exactly `0.025` and does not reject $H_3$, since `0.025 < 0.025` is false.
* The table computes `0.025000000000000005` and does reject $H_3$.

Before the fix, under `"auto"`, one copy of this row gave one answer and four copies gave the other, because four rows switched to the table kernel.

### 7.2 The Fix: Detect Near-Ties and Defer

The sequential kernel is the reference, so the table kernel defers to it whenever a decision is too close to call.

* **Near-tie flag.** During the walk, a trial is flagged if any p-value it compares lies within a relative $10^{-6}$ of the critical value: $\lvert p_j - c_j(S)\rvert < 10^{-6}\, c_j(S)$ (`TIE_REL_TOL`). That trial is then rerun with `sequential_trial()`.
  - **Rate.** Only about 6 in a million trials are rerun with continuous p-values; with p-values rounded to 3 dp, up to a few percent are.
  - **Why the walk's own comparisons are enough.** Levels only grow as $S$ shrinks. If no comparison the walk makes is close, the sequential kernel cannot keep a hypothesis the walk rejects, or reject one the walk keeps.
* **Ill-conditioned graphs.** The rounding difference between the kernels grows as $\varepsilon/d$, where $d = 1 - g_{kj}g_{jk}$ is the smallest graph-update denominator. When two hypotheses pass almost all their weight to each other, $d$ is tiny. The fixed tolerance then no longer covers the difference, so the whole call runs on the sequential kernel instead (`MIN_DENOM` $= 10^{-6}$). The builder finds $d$ like this:
  - It checks pairs only within three-member sets $\{i, k, l\}$, because denominators only fall as $S$ shrinks and the sequential kernel only updates the graph while three or more hypotheses remain.
  - It counts a pair only if $g_{il}$ or $g_{kl}$ is non-zero. A pair that sends nothing outside itself only ever divides exact zeros, so its denominator cannot cause a difference however small it is. This keeps common optimiser graphs on the fast path: there, small edges are zeroed and a near-degenerate pair is often left with nothing else in its rows.
  - The value can be inspected as `attr(multigrain:::graph_shortcut_table(alpha, w, G), "min_denom")`.

**The result.** The output is identical to the sequential kernel's for every graph, every $N$ and every thread count. A development sweep of 79 million trials over random, sparse, round-number and near-degenerate graphs found no differences. Without the fix it found 3,422. The scripts are in `dev/review/shortcut_table/`, which is untracked.

**The cost.** The extra check makes the table walk roughly 10% slower in low-power settings, where most trials stop after one step, and no slower at high power.

---

## 8. Testing and Working on the Code

### 8.1 Where the Tests Are

| File | What it covers |
|---|---|
| `tests/testthat/test-RcppExports.R` | Agreement with `gMCPLite::graphTest()`, an independent implementation of the graphical procedure. This is the external check that the sequential kernel is right |
| `tests/testthat/test-graph_shortcut_table.R` | Table kernel against the sequential kernel: random graphs, edge cases, ties, options, and the table builder against `calc_local_weights()` |
| `tests/testthat/helper-graph_shortcut.R` | `shortcut_all_kernels()` (every kernel through every backend), `expect_identical_kernels()`, and graph generators |
| `tests/testthat/helper-pvals.R` | The shared p-value fixtures |

`calc_local_weights(w, G)` (C++, exported) computes the local weights of every intersection in the closed-testing layout. It is the slow, general route to the same numbers the table builder produces.

### 8.2 Rules for Changing This Code

* **Do not change the sequential kernel's arithmetic.** Not even an equivalent reordering is allowed. It is the reference for the table kernel and for the `gMCPLite` parity tests. If you must touch it, save its output on a set of graphs and p-values before the change and check that it is `identical()` afterwards. `dev/review/shortcut_table/d1_reference.R` does this.
* **Any change to the table kernel must keep the kernels identical.** Test it with forced kernels:

  ```r
  withr::with_options(
    list(multigrain.shortcut_kernel = "table"),
    multigrain:::graph_shortcut(pvals, alpha, w, G)
  )
  ```

  Include rounded p-values and round-number graphs, because that is where ties happen.
* **Run only the affected test files,** with `testthat::test_file()`. The full suite runs in CI.
* **Keep `testthat` cases small** (at most $10^4$ rows per call). Large sweeps belong in dev scripts.
* **Compare speeds carefully.** Only compare timings made with the same `num_threads` on the same machine, and include a same-code control.

---

## 9. Glossary and References

| Term | Meaning |
|---|---|
| FWER | Familywise error rate: the probability of at least one false rejection |
| Weight $w_j$ | Share of $\alpha$ currently allocated to $H_j$ |
| Transition weight $g_{jk}$ | Share of $H_j$'s weight passed to $H_k$ when $H_j$ is rejected |
| Local weight $w_j(S)$, critical value $c_j(S)$ | $H_j$'s weight, and its level $\alpha w_j(S)$, once exactly the hypotheses outside $S$ have been rejected |
| Closed test | Rejects $H_j$ if every intersection hypothesis containing it is rejected; the graphical procedure is a shortcut for one |
| Sequential kernel | One-at-a-time implementation that updates the graph after each rejection; the reference |
| Table kernel | Precomputes $c_j(S)$ for every $S$ and walks the table |
| Near-tie | A p-value within a relative $10^{-6}$ of a critical value the walk compares it with |

**References:**
* Bretz F, Maurer W, Brannath W, Posch M (2009). A graphical approach to sequentially rejective multiple test procedures. *Statistics in Medicine* 28:586–604. This is the algorithm and the order-independence result.
* Bretz F, Posch M, Glimm E, Klinglmueller F, Maurer W, Rohmeyer K (2011). Graphical approaches for multiple comparison procedures using weighted Bonferroni, Simes, or parametric tests. *Biometrical Journal* 53:894–913.
* Xi D, Brockmann E (2024). graphicalMCP: Graphical Multiple Comparison Procedures. R package. Its `graph_test_shortcut_fast()` is the source of the table idea.
* gMCPLite (R package): the independent implementation used as the test oracle.
