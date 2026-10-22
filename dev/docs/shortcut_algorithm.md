# Evaluating Graphical Multiple-Testing Procedures in `multigrain`

*Developer documentation for `graph_shortcut()` and `graph_shortcut_parallel()` ([`src/graph_shortcut.cpp`](https://github.com/GSK-Biostatistics/multigrain/blob/main/src/graph_shortcut.cpp)). Written for statisticians and software developers familiar with multiple testing and simulation.*

**Contents**

1. [Overview](#1-overview)
2. [Statistical Background: Graphical Procedures and Closed Testing](#2-statistical-background-graphical-procedures-and-closed-testing)
3. [Worked Example ($m = 3$)](#3-worked-example-m--3)
4. [The Sequential Method (`sequential_trial()`)](#4-the-sequential-method-sequential_trial)
5. [The Lookup-Table Method (`CritTableBuilder` and `walk_crit_table()`)](#5-the-lookup-table-method-crittablebuilder-and-walk_crit_table)
6. [Method Selection and Parallel Execution](#6-method-selection-and-parallel-execution)
7. [Floating-Point Arithmetic and Exact Agreement Between Methods](#7-floating-point-arithmetic-and-exact-agreement-between-methods)
8. [Verification, Testing, and Maintenance Guidelines](#8-verification-testing-and-maintenance-guidelines)
9. [Notation Summary and References](#9-notation-summary-and-references)

---

## 1. Overview

`multigrain` designs and optimises graphical multiple-testing procedures for clinical trials, estimating their operating characteristics by Monte Carlo simulation. A simulation generates an $N \times m$ matrix of $p$-values `pvals`, where rows index simulated trials ($1, \dots, N$) and columns index null hypotheses $H_1, \dots, H_m$.

Given `pvals` and a graphical testing procedure specified by an overall significance level $\alpha$, an initial weight vector $\mathbf{w}$, and a transition matrix $\mathbf{G}$, `graph_shortcut()` determines which hypotheses are rejected in each simulated trial:

```r
rej <- multigrain:::graph_shortcut(pvals, alpha, w, G)
# pvals: N x m numeric matrix of simulated p-values (N trials, m hypotheses)
# alpha: scalar familywise significance level (e.g. 0.025)
# w:     numeric vector of length m containing initial hypothesis weights
# G:     m x m numeric matrix of transition weights
# rej:   N x m logical matrix; rej[i, j] is TRUE if trial i rejects H_j
```

`graph_shortcut_parallel()` provides a multithreaded implementation using `RcppParallel` and returns the exact same logical matrix regardless of the number of threads.

All simulation-based power summaries in the package are computed directly from `rej`:

* **Marginal power** for each hypothesis: `colMeans(rej)`.
* **Disjunctive power** (probability of rejecting at least one hypothesis).
* **Conjunctive power** (probability of rejecting all hypotheses).
* **Custom trial-success criteria** defined via `trial_success()`.

Within the package, `graph_shortcut()` and `graph_shortcut_parallel()` are called in three places:

| Caller | Source file | Role |
|---|---|---|
| `calc_power_pvals()` | [`R/calc_power.R`](https://github.com/GSK-Biostatistics/multigrain/blob/main/R/calc_power.R) | Evaluates power metrics for a single user-supplied graph |
| Objective closure in `create_obj_func()` | [`R/objective_function.R`](https://github.com/GSK-Biostatistics/multigrain/blob/main/R/objective_function.R) | Evaluates candidate graphs repeatedly during numerical optimisation |
| `graph_optimise()` and helpers | [`R/optimisation.R`](https://github.com/GSK-Biostatistics/multigrain/blob/main/R/optimisation.R) | Evaluates the optimised graph on the full matrix of simulated $p$-values |

Because the optimiser may evaluate thousands of candidate graphs across $N = 10^4$ to $10^5$ simulated trials in a single run, `src/graph_shortcut.cpp` provides two mathematically equivalent evaluation methods (referred to in the C++ source and internal options as `"sequential"` and `"table"` kernels):

1. **The sequential method** (Section 4), which applies the Bretz et al. (2009) graph-update equations step-by-step within each simulated trial.
2. **The lookup-table method** (Section 5), which precomputes the local critical values for all $2^m - 1$ non-empty intersection hypotheses once per graph and evaluates each trial by table lookup.

Both methods are constructed to return identical rejection matrices for every input (Section 7), so choosing between them affects execution time only.

---

## 2. Statistical Background: Graphical Procedures and Closed Testing

### 2.1 Specification of a Graphical Procedure

Consider testing $m$ null hypotheses $H_1, \dots, H_m$ while controlling the familywise error rate (FWER) strongly at level $\alpha$. In the graphical framework of Bretz et al. (2009), a multiple-testing procedure is represented by a directed graph with nodes $H_1, \dots, H_m$ parameterised by:

* An **initial weight vector** $\mathbf{w} = (w_1, \dots, w_m)$ satisfying $w_j \ge 0$ for all $j$ and $\sum_{j=1}^m w_j \le 1$. Initially, $H_j$ is allocated local significance level $\alpha w_j$.
* An $m \times m$ **transition matrix** $\mathbf{G} = (g_{jk})$ satisfying $g_{jj} = 0$, $g_{jk} \ge 0$ for all $j, k$, and $\sum_{k=1}^m g_{jk} \le 1$ for every row $j$. The entry $g_{jk}$ specifies the proportion of $H_j$'s weight that is transferred to $H_k$ if $H_j$ is rejected.

### 2.2 Sequential Graph Updates (Bretz et al., 2009, Algorithm 1)

For a single observed vector of $p$-values $\mathbf{p} = (p_1, \dots, p_m)$, let $S \subseteq \{1, \dots, m\}$ denote the index set of currently unrejected hypotheses, initialised to $S = \{1, \dots, m\}$. Testing proceeds iteratively:

1. **Test against current thresholds:** Identify any unrejected hypothesis $j \in S$ with $p_j < \alpha w_j$. (Following `gMCPLite::graphTest()`, `multigrain` uses strict inequality `<`; see Section 4.) If no such $j$ exists, stop and retain all hypotheses remaining in $S$.
2. **Reject $H_j$ and redistribute its weight:** Remove $j$ from $S$ ($S \leftarrow S \setminus \{j\}$) and update the weights of all remaining hypotheses $k \in S$:
   $$w_k \leftarrow w_k + w_j g_{jk}.$$
3. **Update the transition matrix:** For every ordered pair of distinct remaining hypotheses $k, l \in S$ ($k \ne l$), update the transition weight from $k$ to $l$:
   $$g_{kl} \leftarrow \frac{g_{kl} + g_{kj} g_{jl}}{1 - g_{kj} g_{jk}},$$
    provided the denominator $1 - g_{kj} g_{jk} > 0$. If $g_{kj} g_{jk} = 1$ (which occurs when $H_k$ and $H_j$ form a closed cycle transferring $100\%$ of their weight to each other), set $g_{kl} \leftarrow 0$ for all $l \in S$.
4. **Repeat:** Return to step 1 with the updated weights and transition matrix on $S$.

In the transition update of step 3, the numerator $g_{kl} + g_{kj} g_{jl}$ combines the direct edge $k \to l$ with the two-step path $k \to j \to l$ passing through the newly rejected node $j$. The denominator $1 - g_{kj} g_{jk} = (\sum_{r=0}^\infty (g_{kj} g_{jk})^r)^{-1}$ renormalises for weight cycling between $k$ and $j$.

### 2.3 Connection to Closed Testing, Path Independence, and Monotonicity

Algorithm 1 is a consonant shortcut for a closed testing procedure based on weighted Bonferroni intersection tests. Understanding why this shortcut works explains both the sequential and lookup-table implementations in `multigrain`.

Recall that a closed testing procedure defines a local level-$\alpha$ test for every non-empty intersection hypothesis $H_S = \bigcap_{j \in S} H_j$, where $S \subseteq \{1, \dots, m\}$, and rejects an elementary hypothesis $H_j$ if and only if $H_S$ is rejected for every $S \ni j$. For a weighted Bonferroni closed test, each intersection $H_S$ has a vector of local weights $\{w_j(S) : j \in S\}$ satisfying $\sum_{j \in S} w_j(S) \le 1$, and $H_S$ is rejected if $p_j < c_j(S)$ for at least one $j \in S$, where

$$c_j(S) = \alpha w_j(S), \qquad j \in S$$

are the **local critical values** for subset $S$.

Bretz et al. (2009, Theorem 1 and Appendix) prove two fundamental properties of the update equations in Section 2.2:

1. **Path independence (uniqueness of $w_j(S)$):** Starting from $(\mathbf{w}, \mathbf{G})$ on $\{1, \dots, m\}$ and sequentially removing the hypotheses in $\{1, \dots, m\} \setminus S$ via the update equations yields the exact same local weights $w_j(S)$ and transition weights $g_{kl}(S)$ regardless of the order in which the hypotheses outside $S$ were removed. Thus, $w_j(S)$ and $c_j(S) = \alpha w_j(S)$ are well-defined functions of the subset $S$ alone.
2. **Monotonicity and consonance:** Because all weights and transition entries are non-negative and $1 - g_{kj} g_{jk} \le 1$, removing a hypothesis never decreases the weights or transition entries among the remaining hypotheses:
   $$S' \subseteq S \implies w_j(S') \ge w_j(S) \quad \text{and} \quad c_j(S') \ge c_j(S) \qquad \text{for all } j \in S'.$$
   Monotonicity guarantees **consonance**: whenever an intersection hypothesis $H_S$ is rejected because $p_j < c_j(S)$ for some $j \in S$, every smaller intersection $H_{S'}$ with $j \in S' \subset S$ has $p_j < c_j(S) \le c_j(S')$ and is therefore rejected automatically. Consequently, testing all $2^m - 1$ intersection hypotheses reduces to tracking the shrinking set $S$ of currently unrejected hypotheses.

These two properties have direct computational consequences:

* **Invariance to rejection order within a trial:** If several hypotheses in $S$ simultaneously satisfy $p_j < c_j(S)$, rejecting the lowest-indexed qualifying hypothesis first (as the sequential method does) and rejecting all qualifying hypotheses in $S$ simultaneously in a single batch (as the lookup-table method does) both lead to the same final set of rejected hypotheses.
* **Shared critical values across simulated trials:** Because $c_j(S)$ depends only on the subset $S$ and the graph $(\alpha, \mathbf{w}, \mathbf{G})$—and not on the $p$-values—the critical values for all $2^m - 1$ non-empty subsets $S$ can be computed once before looping over the $N$ simulated trials.

---

## 3. Worked Example ($m = 3$)

Consider $m = 3$ hypotheses tested at $\alpha = 0.025$ with initial weights and transition matrix

$$\mathbf{w} = (0.5,\ 0.5,\ 0), \qquad \mathbf{G} = \begin{pmatrix} 0 & 0.5 & 0.5 \\ 0.5 & 0 & 0.5 \\ 0.5 & 0.5 & 0 \end{pmatrix}.$$

Here $H_1$ and $H_2$ split the initial $\alpha$ equally ($c_1(\{1,2,3\}) = c_2(\{1,2,3\}) = 0.0125$), while $H_3$ has initial weight $0$ and can only be rejected after at least one of $H_1$ or $H_2$ is rejected.

Suppose a simulated trial produces $\mathbf{p} = (0.010,\ 0.018,\ 0.015)$:

| Step | Unrejected set $S$ | Local critical values $(c_1(S), c_2(S), c_3(S))$ | Comparison and decision |
|---:|---|---|---|
| 1 | $\{H_1, H_2, H_3\}$ | $(0.0125,\ 0.0125,\ 0)$ | $p_1 = 0.010 < 0.0125 \implies$ reject $H_1$ |
| 2 | $\{H_2, H_3\}$ | $(0,\ 0.01875,\ 0.00625)$ | $p_2 = 0.018 < 0.01875 \implies$ reject $H_2$ |
| 3 | $\{H_3\}$ | $(0,\ 0,\ 0.025)$ | $p_3 = 0.015 < 0.025 \implies$ reject $H_3$ |

After $H_1$ is rejected in step 1, the updated weights and transition entries on $S = \{H_2, H_3\}$ are obtained from the formulas in Section 2.2:

* $w_2(\{2, 3\}) = 0.5 + 0.5 \times 0.5 = 0.75$, giving $c_2(\{2, 3\}) = 0.025 \times 0.75 = 0.01875$.
* $w_3(\{2, 3\}) = 0 + 0.5 \times 0.5 = 0.25$, giving $c_3(\{2, 3\}) = 0.025 \times 0.25 = 0.00625$.
* $g_{23}(\{2, 3\}) = \frac{0.5 + 0.5 \times 0.5}{1 - 0.5 \times 0.5} = \frac{0.75}{0.75} = 1$, and by symmetry $g_{32}(\{2, 3\}) = 1$.

When $H_2$ is rejected in step 2, its entire weight is passed to $H_3$, yielding $w_3(\{3\}) = 0.25 + 0.75 \times 1 = 1$ and $c_3(\{3\}) = 0.025$. All three hypotheses are rejected even though neither $p_2 = 0.018$ nor $p_3 = 0.015$ met its initial threshold. By contrast, if $\mathbf{p} = (0.010,\ 0.020,\ 0.015)$, step 2 fails to reject either $H_2$ ($0.020 \ge 0.01875$) or $H_3$ ($0.015 \ge 0.00625$), so the procedure stops with only $H_1$ rejected.

### 3.1 Complete Table of Local Critical Values

The internal helper `multigrain:::graph_shortcut_table(0.025, w, G)` returns the $2^3 \times 3$ matrix of local critical values $c_j(S)$ for all subsets $S \subseteq \{H_1, H_2, H_3\}$. As explained in Section 5.2, each subset $S$ is indexed by an integer `mask` $= \sum_{j \in S} 2^{j-1} \in \{0, \dots, 7\}$, stored in row `mask + 1` of the R matrix:

| `mask` | R row | Subset $S$ | $c_1(S)$ | $c_2(S)$ | $c_3(S)$ |
|---:|---:|---|---:|---:|---:|
| 0 | 1 | $\emptyset$ | 0 | 0 | 0 |
| 1 | 2 | $\{H_1\}$ | 0.025 | 0 | 0 |
| 2 | 3 | $\{H_2\}$ | 0 | 0.025 | 0 |
| 3 | 4 | $\{H_1, H_2\}$ | 0.0125 | 0.0125 | 0 |
| 4 | 5 | $\{H_3\}$ | 0 | 0 | 0.025 |
| 5 | 6 | $\{H_1, H_3\}$ | 0.01875 | 0 | 0.00625 |
| 6 | 7 | $\{H_2, H_3\}$ | 0 | 0.01875 | 0.00625 |
| 7 | 8 | $\{H_1, H_2, H_3\}$ | 0.0125 | 0.0125 | 0 |

By convention, $c_j(S) = 0$ whenever $H_j \notin S$. Note that $c_j(S') \ge c_j(S)$ whenever $S' \subseteq S$, and that the table includes all $2^m$ subsets even when a subset cannot be reached during testing from the initial weights (for instance, reaching $S = \{H_1, H_2\}$ would require rejecting $H_3$ first, which is impossible because $c_3(\{H_1, H_2, H_3\}) = 0$).

---

## 4. The Sequential Method (`sequential_trial()`)

The sequential method is the direct per-trial implementation of Algorithm 1 in Section 2.2 and serves as the reference implementation against which the lookup-table method is validated. It is implemented in [`sequential_trial()`](https://github.com/GSK-Biostatistics/multigrain/blob/main/src/graph_shortcut.cpp#L246-L322), which is called by the single-threaded loop in `graph_shortcut()`, the multithreaded worker `SequentialShortcutWorker`, and the near-tie fallback in `walk_crit_table()`.

Before iterating over trials, the caller precomputes the initial critical values `init_a[i] = w[i] * alpha` ($i = 0, \dots, m - 1$ in 0-based C++ indexing) and allocates reusable working vectors (`SequentialScratch`: `cur_a`, `cur_G`, `new_G`, and `cur_p`) so that no dynamic memory allocation occurs inside the trial loop. For each trial row `set`, `sequential_trial()` executes the following steps:

1. **Initialise trial state:** Copy `init_a` into `cur_a`, copy the $m \times m$ transition matrix `G_ptr` into `cur_G`, and extract the $m$ $p$-values for row `set` from `p_ptr` into the contiguous vector `cur_p`. Because R stores matrices in column-major order, entry `(set, i)` of the $N \times m$ matrix `pvals` sits at memory offset `set + i * N`; copying the row once into `cur_p` avoids strided memory reads during the inner rejection loop.
2. **Find the first qualifying hypothesis:** Scan indices $i = 0, \dots, m - 1$ in order and select the first index `rej` satisfying `cur_p[i] < cur_a[i]`. If no index satisfies this condition, terminate the trial.
3. **Record the rejection:** Set `h_ptr[set + rej * N] = 1` in the output matrix and increment the trial's rejection counter `sumrej`. If `sumrej == m`, all hypotheses have been rejected and the trial terminates immediately.
4. **Update local critical values:** Redistribute `cur_a[rej]` to the remaining hypotheses directly on the $\alpha$-scale:
   ```cpp
   for (int i = 0; i < m; ++i) {
     cur_a[i] += cur_a[rej] * cur_G[rej + i * m];
   }
   ```
5. **Update the transition matrix (if at least two unrejected hypotheses remain):** When `sumrej < m - 1`, compute the updated transition weights into `new_G`. For each row $i$, compute `denom = 1.0 - gir * gri`, where `gir = cur_G[i + rej * m]` and `gri = cur_G[rej + i * m]`. If `denom != 0.0`, hoist the reciprocal `inv_denom = 1.0 / denom` outside the inner loop over columns $j$ and set `new_G[i + j * m] = (gij + gir * grj) * inv_denom` for $j \ne i$ (with zeros on the diagonal). If `denom == 0.0`, set row $i$ of `new_G` to zero. Swap the internal buffer pointers of `cur_G` and `new_G` via `std::swap(cur_G, new_G)` in $O(1)$ time, and zero out row `rej` and column `rej` of `cur_G`. (When `sumrej == m - 1`, only one unrejected hypothesis remains, so the transition matrix will never be used again and step 5 is skipped.)
6. **Zero the rejected threshold and repeat:** Set `cur_a[rej] = 0.0` (ensuring hypothesis `rej` cannot satisfy `cur_p[i] < cur_a[i]` in subsequent iterations) and return to step 2.

Two numerical conventions in `sequential_trial()` are essential for exact reproducibility:

* **Operating on the $\alpha$-scale (`cur_a`):** Multiplying the weight update $w_k \leftarrow w_k + w_j g_{jk}$ through by $\alpha$ gives $c_k \leftarrow c_k + c_j g_{jk}$. Storing $c_k = \alpha w_k$ directly avoids multiplying by $\alpha$ on every comparison in step 2.
* **Strict inequality (`cur_p[i] < cur_a[i]`):** A hypothesis with $p_i = c_i(S)$ is not rejected. In particular, an unrejected hypothesis with $c_i(S) = 0$ and $p_i = 0$ is not rejected, and once a hypothesis is rejected and its threshold is set to `cur_a[rej] = 0.0`, it can never be selected again even if its $p$-value is `0.0`.

**Computational complexity:** Copying the initial $m \times m$ transition matrix costs $O(m^2)$ per trial, and each rejection before the penultimate one triggers an $O(m^2)$ transition matrix update. A trial that rejects $r$ hypotheses therefore requires $O(m^2(1 + r))$ operations, and evaluating $N$ trials requires $O(N m^2 (1 + \bar{r}))$ operations, where $\bar{r}$ is the mean number of rejections per trial.

---

## 5. The Lookup-Table Method (`CritTableBuilder` and `walk_crit_table()`)

### 5.1 Precomputing Intersection Critical Values

In typical clinical trial simulations, the number of hypotheses $m$ is modest (commonly $m \in \{2, \dots, 12\}$), whereas the number of simulated trials $N$ is large ($10^4$ to $10^5$). Under the sequential method, the same graph updates are recomputed independently in every trial.

Because the local critical values $c_j(S) = \alpha w_j(S)$ depend only on the subset $S \subseteq \{1, \dots, m\}$ of unrejected hypotheses (Section 2.3), the lookup-table method splits the computation into two stages:

1. **Table construction (once per call):** Compute $c_j(S)$ for all $2^m$ subsets $S \subseteq \{1, \dots, m\}$ and store them in a flat $2^m \times m$ array `crit`.
2. **Per-trial table lookup (for each of the $N$ trials):** Starting at $S = \{1, \dots, m\}$, read the row of `crit` corresponding to $S$, reject all hypotheses $j \in S$ with $p_j < c_j(S)$, remove them from $S$, and repeat until no further hypotheses in $S$ can be rejected.

This two-stage strategy was introduced by Xi and Brockmann (2024) in `graphicalMCP::graph_test_shortcut_fast()`; `src/graph_shortcut.cpp` provides a standalone C++ implementation designed to agree bit-for-bit with `sequential_trial()`.

### 5.2 Indexing Subsets $S \subseteq \{1, \dots, m\}$ with Binary Integers (Bitmasks)

To look up the critical values for a subset $S$ in $O(1)$ time without hashing or pointer indirection, each subset $S \subseteq \{1, \dots, m\}$ is represented as an unsigned 32-bit integer (`std::uint32_t`) whose binary digits (bits) indicate which hypotheses belong to $S$. In 0-based C++ indexing, bit $i$ (corresponding to $2^i$) is `1` if $H_{i+1} \in S$ and `0` if $H_{i+1} \notin S$:

$$\text{mask}(S) = \sum_{j \in S} 2^{j-1} \in \{0, 1, \dots, 2^m - 1\}.$$

For $m = 3$, this bijection maps the $2^3 = 8$ subsets to the integers $0, \dots, 7$ shown in Section 3.1:

| Subset $S$ | Binary representation ($H_3 H_2 H_1$) | Integer `mask` |
|---|---|---:|
| $\emptyset$ | `000` | $0$ |
| $\{H_1\}$ | `001` | $1$ |
| $\{H_2\}$ | `010` | $2$ |
| $\{H_1, H_2\}$ | `011` | $1 + 2 = 3$ |
| $\{H_3\}$ | `100` | $4$ |
| $\{H_1, H_3\}$ | `101` | $1 + 4 = 5$ |
| $\{H_2, H_3\}$ | `110` | $2 + 4 = 6$ |
| $\{H_1, H_2, H_3\}$ | `111` | $1 + 2 + 4 = 7$ |

Because `mask(S)` takes every integer value from $0$ to $2^m - 1$ exactly once, the critical values $c_1(S), \dots, c_m(S)$ are stored contiguously in row `mask(S)` of the row-major array `crit` (at memory offset `S * m`). Moreover, standard set operations on $S$ reduce to single CPU bitwise instructions:

| Set operation | C++ expression (`src/graph_shortcut.cpp`) | Explanation | Equivalent R expression |
|---|---|---|---|
| Full set $\{1, \dots, m\}$ | `(std::uint32_t(1) << m) - 1` | Shift `1` left by $m$ bits ($2^m$) and subtract $1$, setting the lowest $m$ bits to `1` | `2^m - 1` |
| Test membership $H_{i+1} \in S$ | `(S >> i) & 1u` | Shift $S$ right by $i$ bits and extract bit $0$ | `bitwAnd(S, 2^i) != 0` |
| Remove single hypothesis $S \setminus \{H_{k+1}\}$ | `S & ~(std::uint32_t(1) << k)` | Clear bit $k$ using bitwise AND with the bitwise complement of $2^k$ | `bitwAnd(S, bitwNot(2^k))` |
| Remove rejected batch $S \setminus R$ | `S &= ~rej` | Clear all bits that are set in the rejected bitmask `rej` | `bitwAnd(S, bitwNot(rej))` |
| Test whether $S = \emptyset$ | `S == 0` | All bits are zero | `S == 0` |
| Rejected complement $\{1, \dots, m\} \setminus S$ | `full & ~S` | Bits in the full set that are no longer in $S$ | `bitwAnd(2^m - 1, bitwNot(S))` |

When inspecting the table from R via `tab <- multigrain:::graph_shortcut_table(alpha, w, G)`, row `mask + 1` contains the critical values for subset `mask` (since R uses 1-based row indexing):

```r
mask <- sum(2^(c(1, 3) - 1))  # S = {H1, H3} -> mask = 5
tab[mask + 1, ]               # returns c_1(S), ..., c_m(S)
```

**Memory footprint and upper bound on $m$:** The table stores $2^m \times m$ double-precision floats ($8$ bytes each). At $m = 12$, this requires $4096 \times 12 \times 8 \approx 393\text{ KB}$; at the hard limit `TABLE_HARD_MAX_M = 20`, it requires $2^{20} \times 20 \times 8 = 160\text{ MB}$. Bounding $m \le 20$ limits memory consumption and keeps the 32-bit bitmask shift `1u << m` well within the range of `std::uint32_t`.

### 5.3 Constructing the Table via Depth-First Traversal (`CritTableBuilder`)

Rather than computing each subset $S$ from scratch, [`CritTableBuilder`](https://github.com/GSK-Biostatistics/multigrain/blob/main/src/graph_shortcut.cpp#L375-L470) constructs the $2^m$ rows of `crit` via a depth-first traversal of the subset lattice starting from the full set $S = \{1, \dots, m\}$ (`depth = 0`).

To ensure that every non-empty subset $S$ is visited once along a unique path, hypotheses are removed in **strictly increasing index order**. If the most recently removed hypothesis leading to the current subset $S$ had 0-based index `last` (with `last = -1` at the root), `dfs(S, depth, last)` only branches by removing indices $k \in \{\text{last} + 1, \dots, m - 1\}$. Consequently, every proper subset $S$ is generated from the unique parent subset $S \cup \{k_{\max}\}$, where $k_{\max}$ is the largest index absent from $S$.

At each node `(S, depth, last)` of the recursion:

1. The current local critical values `as_[depth][0..m-1]` are copied into row `S` of `crit`.
2. For each $k > \text{last}$ whose removal leaves a non-empty child subset `child = S & ~(1u << k)`, the child's critical values `as_[depth + 1]` and updated transition matrix `Gs_[depth + 1]` are computed from `as_[depth]` and `Gs_[depth]` using the exact same floating-point formulas as `sequential_trial()`.
3. The transition matrix update into `Gs_[depth + 1]` is omitted when `child` is a leaf of the recursion tree (`k == m - 1`) and does not have three members (`m - depth - 1 != 3`), because `child` will not generate further recursive branches and its transition matrix is only needed at three-member subsets to evaluate the conditioning check described in Section 7.2.

Because the recursion has maximum depth $m$, `CritTableBuilder` allocates only $m + 1$ working level vectors (`as_`) and $m + 1$ working $m \times m$ matrices (`Gs_`), requiring $O(m^3)$ working memory in addition to the $O(2^m m)$ output table. Visiting all $2^m - 1$ non-empty subsets takes $O(2^m m^2)$ operations in total, independent of the number of trials $N$.

### 5.4 Evaluating a Trial by Table Lookup (`walk_crit_table()`)

Once `crit` has been populated, [`walk_crit_table()`](https://github.com/GSK-Biostatistics/multigrain/blob/main/src/graph_shortcut.cpp#L492-L535) evaluates each trial `set` without performing any weight or matrix arithmetic:

1. Extract the trial's $p$-values into `cur_p[0..m-1]` and initialise the unrejected set bitmask to `S = full` ($2^m - 1$).
2. Point `row` to `crit + S * m`, which holds $c_1(S), \dots, c_m(S)$.
3. Loop over $i = 0, \dots, m - 1$ and accumulate a bitmask `rej` of all hypotheses in $S$ that satisfy $p_{i+1} < c_{i+1}(S)$:
   ```cpp
   rej |= static_cast<std::uint32_t>(((S >> i) & 1u) & (cur_p[i] < row[i])) << i;
   ```
   (At the same time, check whether $|p_{i+1} - c_{i+1}(S)| < 10^{-6} c_{i+1}(S)$ to detect floating-point near-ties; see Section 7.2.)
4. If `rej == 0`, no further hypotheses in $S$ can be rejected; exit the `while (S)` loop. Otherwise, remove all newly rejected hypotheses simultaneously via `S &= ~rej` and return to step 2.
5. After the loop terminates, the rejected set is the bitmask complement `R = full & ~S`. For each $i$ with bit $i$ set in `R`, write `h_ptr[set + i * N] = 1`.

As shown in Section 2.3, rejecting all qualifying hypotheses in `rej` simultaneously in step 4 is guaranteed by monotonicity ($c_j(S') \ge c_j(S)$ for $S' \subseteq S$) to reach the same final unrejected set $S$ as rejecting them one at a time.

To see the difference in work per trial, consider the $m = 3$ graph from Section 3 evaluated on a trial with $\mathbf{p} = (0.010,\ 0.012,\ 0.030)$:

* **Sequential method (`sequential_trial()`):**
  - *Iteration 1 ($S = \{H_1, H_2, H_3\}$):* Scans thresholds $(0.0125, 0.0125, 0)$, stops at the first qualifying hypothesis $H_1$ ($0.010 < 0.0125$), rejects $H_1$, updates the critical-value vector, and performs a $3 \times 3$ transition matrix update (including floating-point divisions).
  - *Iteration 2 ($S = \{H_2, H_3\}$):* Scans updated thresholds $(0, 0.01875, 0.00625)$, rejects $H_2$ ($0.012 < 0.01875$), and updates the critical-value vector.
  - *Iteration 3 ($S = \{H_3\}$):* Tests $p_3 = 0.030$ against $c_3(\{H_3\}) = 0.025$, fails to reject, and terminates.
* **Lookup-table method (`walk_crit_table()`):**
  - *Iteration 1 (`S = 7`, binary `111`, $\{H_1, H_2, H_3\}$):* Reads row $7$ of `crit` $(0.0125, 0.0125, 0)$, sees that both $H_1$ and $H_2$ pass (`rej = 3`, binary `011`), and clears both bits at once (`S = 7 & ~3 = 4`, binary `100`).
  - *Iteration 2 (`S = 4`, binary `100`, $\{H_3\}$):* Reads row $4$ of `crit` $(0, 0, 0.025)$, sees that $p_3 = 0.030 \ge 0.025$ (`rej = 0`), and terminates with rejected mask `R = 7 & ~4 = 3` ($\{H_1, H_2\}$).

When hypotheses only qualify cascade-fashion one after another (as in Section 3 with $\mathbf{p} = (0.010, 0.018, 0.015)$, where `walk_crit_table()` visits rows `7 -> 6 -> 4 -> 0`), the lookup-table method takes the same number of iterations as the sequential method, but each iteration performs $m$ comparisons against precomputed table entries rather than an $O(m^2)$ graph update.

---

## 6. Method Selection and Parallel Execution

### 6.1 Automatic Selection Rule and Diagnostic Options

Constructing the lookup table costs $O(2^m m^2)$ once per call, whereas the sequential method costs $O(N m^2 (1 + \bar{r}))$ across $N$ trials. Precomputing the table is therefore advantageous when $2^m$ is not too large relative to $N$.

Under the default `"auto"` setting, [`use_table_kernel()`](https://github.com/GSK-Biostatistics/multigrain/blob/main/src/graph_shortcut.cpp#L168-L184) selects the lookup-table method when both of the following conditions hold:

1. $1 \le m \le m_{\max}$, where $m_{\max} \le 20$ defaults to `TABLE_HARD_MAX_M = 20`.
2. $2^m \le \kappa N$, where $\kappa =$ `TABLE_KAPPA` $= 2.0$.

The threshold $\kappa = 2.0$ was calibrated from single-threaded benchmarks across varying numbers of hypotheses and rejection rates (where the lookup-table method was faster in every scenario up to $2^m / N = 1.64$, and the sequential method first became faster in low-power scenarios at $2^m / N = 2.62$). To keep behaviour predictable, the selection rule depends only on $m$ and $N$, and not on `num_threads`.

Two internal R options—read once on the main thread at the start of each call—allow tests and benchmarks to override the automatic rule:

| R option | Valid values | Effect |
|---|---|---|
| `multigrain.shortcut_kernel` | `"auto"` (default), `"sequential"`, `"table"` | Forces either the sequential method (`"sequential"`) or the lookup-table method (`"table"`, valid for $m \le 20$) regardless of $N$. |
| `multigrain.shortcut_table_max_m` | Integer in $0, \dots, 20$ (default `20`) | Lowers the maximum $m$ for which `"auto"` will select the lookup-table method. |

The internal helper `multigrain:::graph_shortcut_kernel(n_trials, m)` returns `"table"` or `"sequential"` to report which method `graph_shortcut()` will choose under the currently active options.

### 6.2 Multithreading with `RcppParallel`

Because trials are independent, `graph_shortcut_parallel()` partitions the trial rows $0, \dots, N - 1$ across threads using `RcppParallel::parallelFor` and Intel TBB work-stealing workers (`SequentialShortcutWorker` and `TableShortcutWorker`):

* **Shared read-only data:** Raw pointers to `pvals`, `G`, `init_a`, and (when using the lookup-table method) `crit` are obtained on the main R thread before entering `parallelFor`. When the lookup-table method is used, `build_crit_table()` builds `crit` once on the main thread and all worker threads read from it concurrently.
* **Disjoint output writes and thread-local scratch space:** Worker threads never call the R C API. Each trial `set` writes only to entries `set + i * N` ($i = 0, \dots, m - 1$) of the pre-allocated output matrix `ans`, avoiding write contention. Scratch buffers (`SequentialScratch` or `cur_p`) are allocated locally inside `Worker::operator()(begin, end)`, once per dispatched chunk of trials.
* **Chunk size (`grain_size`) tuning:** When `grain_size = -1` (the default), [`choose_grain_size()`](https://github.com/GSK-Biostatistics/multigrain/blob/main/src/graph_shortcut.cpp#L189-L219) chooses the minimum chunk size so that each chunk performs roughly $10^5$ inner-loop operations (estimating $m^3$ operations per trial for the sequential method and $m^2$ operations per trial for the lookup-table method), subject to creating at least $4 \times \text{num\_threads}$ chunks for load balancing and clamping the chunk size to $[1, 10000]$.

---

## 7. Floating-Point Arithmetic and Exact Agreement Between Methods

### 7.1 Why Removal Order Affects Floating-Point Critical Values

In exact real arithmetic, the path-independence theorem of Bretz et al. (2009) guarantees that the local critical value $c_j(S)$ is identical regardless of the order in which the hypotheses outside $S$ were rejected (Section 2.3). In IEEE 754 double-precision floating-point arithmetic, however, addition and division are not associative, so computing $c_j(S)$ along two different removal orders can produce values that differ in their least significant bits (a few units in the last place, or ulps).

Specifically:

* **`CritTableBuilder`** always reaches a subset $S$ by removing the hypotheses in $\{1, \dots, m\} \setminus S$ in **increasing index order**.
* **`sequential_trial()`** reaches $S$ by removing hypotheses in the **rejection order dictated by that trial's $p$-values**.

Whenever a trial's rejection order removes a higher-indexed hypothesis before a lower-indexed one, the table entry `crit[S * m + j]` and the value `cur_a[j]` computed by `sequential_trial()` follow different sequences of floating-point operations.

For almost all continuous simulated $p$-values, a discrepancy of a few ulps (roughly $10^{-16} c_j(S)$) has no effect on the comparison $p_j < c_j(S)$. However, if a $p$-value happens to fall on or within rounding error of $c_j(S)$—which occurs naturally when $p$-values or graph weights are rounded to a decimal grid—the two methods can make opposite rejection decisions.

A concrete $m = 3$ example illustrates this phenomenon. Take $\alpha = 0.025$, initial weights $\mathbf{w} = (0.1,\ 0.4,\ 0.5)$, off-diagonal transition weights $g_{jk} = 0.5$ for all $j \ne k$, and trial $p$-values $\mathbf{p} = (0.005,\ 0,\ 0.025)$:

1. At $S = \{H_1, H_2, H_3\}$, the initial critical values are $(0.0025,\ 0.010,\ 0.0125)$. Only $H_2$ qualifies ($p_2 = 0 < 0.010$), so the sequential method rejects $H_2$ first, then rejects $H_1$ at $S = \{H_1, H_3\}$ (where $c_1(\{H_1, H_3\}) = 0.0075 > 0.005$), reaching $S = \{H_3\}$ via the removal order $H_2 \to H_1$. Along this order, `sequential_trial()` computes $c_3(\{H_3\}) =$ `0.025` (exact IEEE 754 double). Since $p_3 = 0.025$, the strict test `0.025 < 0.025` evaluates to `false`, and $H_3$ is **not** rejected.
2. By contrast, `CritTableBuilder` reaches $S = \{H_3\}$ via the increasing-index removal order $H_1 \to H_2$. Along that order, the intermediate floating-point operations produce $c_3(\{H_3\}) =$ `0.025000000000000005` ($1$ ulp above `0.025`). Without a safeguard, `0.025 < 0.025000000000000005` evaluates to `true`, so a naive table lookup would reject $H_3$. Worse still, under `"auto"`, passing $N = 1$ copy of this row would run the sequential method ($2^3 > 2 \times 1$) and retain $H_3$, whereas passing $N = 4$ copies would switch to the lookup-table method ($2^3 \le 2 \times 4$) and reject $H_3$.

### 7.2 Ensuring Bit-for-Bit Identical Rejections

To guarantee that `graph_shortcut()` and `graph_shortcut_parallel()` return the exact same rejections as `sequential_trial()` for every graph, every $N$, and every thread count, the lookup-table implementation incorporates two safeguards: a per-trial near-tie check during table lookup and a conditioning check on the transition matrix denominators during table construction.

#### 1. Per-trial near-tie detection and sequential re-evaluation (`TIE_REL_TOL = 1e-6`)

During `walk_crit_table()`, whenever $p_j$ is compared against the table entry $c_j(S) = \text{row}[i]$, the loop also checks whether $p_j$ lies within a relative tolerance $\tau = 10^{-6}$ (`TIE_REL_TOL`) of $c_j(S)$:

$$|p_j - c_j(S)| < \tau \, c_j(S), \qquad \tau = 10^{-6}.$$

(Because the inequality is strict, a zero critical value $c_j(S) = 0$ never triggers the flag, even when $p_j = 0$; hypotheses outside $S$ or with zero weight have exact $0.0$ in both methods.)

If any comparison at any subset $S$ visited during the trial satisfies this inequality, `walk_crit_table()` immediately aborts the table lookup for that trial and re-evaluates that single trial row using `sequential_trial()`.

* **Why checking only the subsets visited by `walk_crit_table()` is sufficient:** Suppose no comparison made by `walk_crit_table()` lies within $\tau c_j(S)$ of its threshold, and that the relative discrepancy between `CritTableBuilder` and `sequential_trial()` at every subset is bounded by $\tau$. Then at every step, any hypothesis rejected by `walk_crit_table()` satisfies $p_j \le (1 - \tau) c_j^{\text{table}}(S) < c_j^{\text{seq}}(S)$, and by monotonicity of $c_j^{\text{seq}}(\cdot)$ as the unrejected set shrinks, `sequential_trial()` must also reject all of those hypotheses. Conversely, at the terminal unrejected set $S_{\text{final}}$ where `walk_crit_table()` stops, every remaining $j \in S_{\text{final}}$ satisfies $p_j \ge (1 + \tau) c_j^{\text{table}}(S_{\text{final}}) > c_j^{\text{seq}}(S_{\text{final}})$, so `sequential_trial()` must also stop at $S_{\text{final}}$.
* **Frequency of re-evaluation:** For continuous $p$-values, only about $6$ in $10^6$ trials trigger the near-tie check, so the amortised cost of calling `sequential_trial()` on flagged trials is negligible. For $p$-values rounded to $3$ decimal places on round-number graphs, up to a few percent of trials are re-evaluated. Adding the near-tie check inside the inner loop increases lookup time by roughly $8\text{--}12\%$ in low-power scenarios (where most trials stop after one step) and has negligible impact at high power.

#### 2. Conditioning check for near-degenerate cycles (`MIN_DENOM = 1e-6`)

The near-tie tolerance $\tau = 10^{-6}$ is sufficient only as long as the relative floating-point discrepancy between different removal orders is strictly smaller than $10^{-6}$. In the transition update

$$g_{kl} \leftarrow \frac{g_{kl} + g_{kj} g_{jl}}{1 - g_{kj} g_{jk}},$$

dividing by a small denominator $d = 1 - g_{kj} g_{jk}$ amplifies prior rounding differences by a factor proportional to $1 / d$. When two hypotheses $H_i$ and $H_k$ form a nearly closed cycle ($g_{ik} g_{ki} \approx 1$, so $d \ll 1$), the relative difference between critical values computed along different removal orders scales empirically as $O(m \varepsilon_{\text{mach}} / d)$, where $\varepsilon_{\text{mach}} = 2^{-52} \approx 2.22 \times 10^{-16}$. If $d$ is arbitrarily small (e.g. $d = 10^{-12}$), the relative discrepancy between removal orders can exceed $\tau = 10^{-6}$.

To prevent an ill-conditioned graph from exceeding $\tau$, `CritTableBuilder` tracks the smallest relevant update denominator `min_denom` encountered while building the table, and `graph_shortcut()` / `graph_shortcut_parallel()` only use the table if `min_denom >= MIN_DENOM` ($10^{-6}$); otherwise the entire call falls back to `sequential_trial()`. At $d \ge 10^{-6}$ and $m \le 20$, the worst-case relative discrepancy between the two methods is on the order of $5 \times 10^{-9}$ (and in a sweep of $79.2$ million trials across $m = 2, \dots, 12$, never exceeded $2.2 \times 10^{-10}$), remaining more than three orders of magnitude inside $\tau = 10^{-6}$.

The way `CritTableBuilder` computes `min_denom` relies on two structural properties of the graph updates:

1. **Checking three-hypothesis subsets $\{i, k, l\}$ suffices:** Because transition weights $g_{ik}(S)$ only grow as $S$ shrinks, $g_{ik}(S) g_{ki}(S)$ is largest—and the denominator $1 - g_{ik}(S) g_{ki}(S)$ is smallest—at the smallest subsets where a graph update actually occurs. Since `sequential_trial()` skips the graph update once only two unrejected hypotheses remain (`sumrej < m - 1`), the smallest subsets at which a denominator $1 - g_{ik} g_{ki}$ is ever divided by are the three-member subsets $S = \{i, k, l\}$.
2. **Only pairs $(i, k)$ with positive outgoing weight to the third hypothesis $l$ can cause a discrepancy:** At a three-member subset $\{i, k, l\}$, removing $k$ divides the numerator $g_{il} + g_{ik} g_{kl}$ by $1 - g_{ik} g_{ki}$, whereas removing $i$ divides $g_{kl} + g_{ki} g_{il}$ by $1 - g_{ik} g_{ki}$. Because every transition entry is a sum of non-negative products, any entry that is zero in exact arithmetic evaluates to exact floating-point `0.0` under every removal order. Therefore, if both $g_{il} = 0$ and $g_{kl} = 0$ (meaning the pair $\{i, k\}$ sends no weight to $l$), both numerators are exact `0.0`, and dividing `0.0` by any non-zero denominator—or setting the row to `0.0` when the denominator rounds to `0.0`—yields exact `0.0` in both methods. Furthermore, since transition entries are non-decreasing as $S$ shrinks, $g_{il}(S) = g_{kl}(S) = 0$ at $\{i, k, l\}$ implies $g_{il}(S') = g_{kl}(S') = 0$ at every larger superset $S' \supset \{i, k, l\}$.

Accordingly, at each three-member subset $\{i, k, l\}$, [`CritTableBuilder::dfs()`](https://github.com/GSK-Biostatistics/multigrain/blob/main/src/graph_shortcut.cpp#L412-L426) includes the denominator $1 - g_{ik} g_{ki}$ in `min_denom` if and only if $g_{il} \ne 0$ or $g_{kl} \ne 0$. Filtering out pairs with no outgoing weight to $l$ is important in practice: during numerical optimisation, `create_obj_func()` truncates very small edge weights to zero, frequently leaving two hypotheses that pass $100\%$ of their weight to each other ($g_{ik} = g_{ki} = 1$) and $0\%$ elsewhere. Because $g_{il} = g_{kl} = 0$ for every third hypothesis $l$, such graphs are completely unaffected by rounding in $1 - g_{ik} g_{ki}$ and remain on the fast lookup-table path.

For diagnostic and testing purposes, `multigrain:::graph_shortcut_table(alpha, w, G)` attaches `min_denom` as the `"min_denom"` attribute of the returned matrix.

---

## 8. Verification, Testing, and Maintenance Guidelines

### 8.1 Test Suite Structure

The correctness of `graph_shortcut()` and `graph_shortcut_parallel()` is verified by two layers of tests: external validation of the sequential method against an independent R implementation (`gMCPLite::graphTest()`), and internal equivalence tests verifying that the lookup-table method matches the sequential method across serial and parallel backends.

| File | Coverage |
|---|---|
| [`tests/testthat/test-RcppExports.R`](https://github.com/GSK-Biostatistics/multigrain/blob/main/tests/testthat/test-RcppExports.R) | Compares `graph_shortcut()` and `graph_shortcut_parallel()` against `gMCPLite::graphTest()` across standard clinical trial graphs, edge cases, and multi-thread configurations. |
| [`tests/testthat/test-graph_shortcut_table.R`](https://github.com/GSK-Biostatistics/multigrain/blob/main/tests/testthat/test-graph_shortcut_table.R) | Verifies exact agreement between `"sequential"`, `"table"`, and `"auto"` (both serial and parallel) across random dense/sparse graphs, boundary cases ($p = 0$, $p = 1$, zero weights), exact and near-tie $p$-values, rounded $p$-values on grid graphs, near-degenerate pairs around `MIN_DENOM`, and `CritTableBuilder` vs. `calc_local_weights()`. |
| [`tests/testthat/helper-graph_shortcut.R`](https://github.com/GSK-Biostatistics/multigrain/blob/main/tests/testthat/helper-graph_shortcut.R) | Defines `shortcut_all_kernels()` (runs `"sequential"`, `"table"`, and `"auto"` through both serial and parallel entry points), `expect_identical_kernels()`, and random/grid graph generators (`make_grid_graph()`). |
| [`tests/testthat/helper-pvals.R`](https://github.com/GSK-Biostatistics/multigrain/blob/main/tests/testthat/helper-pvals.R) | Provides shared simulated $p$-value matrices (`pvals_fixture` and `pvals_high_power`). |

Note that `calc_local_weights(w, G)` (exported from C++) computes the local weights $w_j(S)$ for all $2^m - 1$ intersection hypotheses in the closed-testing row order used by `apply_ctp()`, providing an independent check on the entries of `graph_shortcut_table()`.

### 8.2 Guidelines for Modifying `src/graph_shortcut.cpp`

1. **Treat `sequential_trial()`'s floating-point arithmetic as frozen:** Because `sequential_trial()` is the reference oracle for both the `gMCPLite` parity tests and the lookup-table equivalence tests, do not reorder its floating-point operations (for example, do not replace `(gij + gir * grj) * inv_denom` with division inside the loop or fused multiply-add instructions that alter rounding). If `sequential_trial()` must be modified, capture its output across a diverse reference suite before the change and verify `identical()` output afterwards.
2. **Preserve exact equivalence between `"sequential"` and `"table"`:** Any modification to `CritTableBuilder` or `walk_crit_table()` must maintain bit-for-bit identical logical output with `sequential_trial()`, tested via `expect_identical_kernels()` or explicit option overrides:
   ```r
   withr::with_options(
     list(multigrain.shortcut_kernel = "table"),
     multigrain:::graph_shortcut(pvals, alpha, w, G)
   )
   ```
   Always include test cases with rounded $p$-values (e.g. `round(pvals, 3)`) and rational-grid graphs (`make_grid_graph()`), where exact ties $p_j = c_j(S)$ occur frequently.
3. **Run targeted test files during development:** Execute `testthat::test_file("tests/testthat/test-graph_shortcut_table.R")` and `testthat::test_file("tests/testthat/test-RcppExports.R")` rather than running the entire test suite on every iteration, and keep individual `testthat` matrices small ($\le 10^4$ rows).

---

## 9. Notation Summary and References

### 9.1 Notation and Constants

| Symbol / Identifier | Description |
|---|---|
| $m$ | Number of null hypotheses $H_1, \dots, H_m$ (`ncol(pvals)`) |
| $N$ | Number of simulated trials (`nrow(pvals)`) |
| $\alpha$ | Overall familywise significance level |
| $\mathbf{w} = (w_1, \dots, w_m)$ | Vector of initial hypothesis weights ($w_j \ge 0$, $\sum_{j=1}^m w_j \le 1$) |
| $\mathbf{G} = (g_{jk})$ | $m \times m$ transition matrix ($g_{jj} = 0$, $g_{jk} \ge 0$, $\sum_{k=1}^m g_{jk} \le 1$) |
| $S \subseteq \{1, \dots, m\}$ | Set of currently unrejected hypotheses (also denotes the index set of the intersection hypothesis $H_S = \bigcap_{j \in S} H_j$) |
| $w_j(S)$ | Local weight of hypothesis $H_j$ ($j \in S$) after all hypotheses outside $S$ have been rejected |
| $c_j(S) = \alpha w_j(S)$ | Local critical value (rejection threshold) for $H_j$ when the unrejected set is $S$ |
| $\text{mask}(S) = \sum_{j \in S} 2^{j-1}$ | 32-bit integer bitmask encoding subset $S \subseteq \{1, \dots, m\}$, used as the 0-based row index into `crit` |
| `TABLE_HARD_MAX_M` ($20$) | Maximum number of hypotheses $m$ supported by the lookup-table method ($2^{20} \times 20 \times 8\text{ B} = 160\text{ MB}$) |
| `TABLE_KAPPA` ($2.0$) | Automatic selection constant: `"auto"` uses the lookup-table method when $m \le 20$ and $2^m \le 2 N$ |
| `TIE_REL_TOL` ($\tau = 10^{-6}$) | Relative near-tie tolerance: if $|p_j - c_j(S)| < 10^{-6} c_j(S)$ during table lookup, that trial is re-evaluated with `sequential_trial()` |
| `MIN_DENOM` ($10^{-6}$) | Minimum conditioning threshold on $1 - g_{ik} g_{ki}$ across active pairs in three-member subsets; below this threshold, the call falls back to the sequential method |

### 9.2 References

* Bretz F, Maurer W, Brannath W, Posch M (2009). A graphical approach to sequentially rejective multiple test procedures. *Statistics in Medicine* 28(4):586–604. (Introduces the graphical testing framework, Algorithm 1, and the path-independence and consonance results.)
* Bretz F, Posch M, Glimm E, Klinglmueller F, Maurer W, Rohmeyer K (2011). Graphical approaches for multiple comparison procedures using weighted Bonferroni, Simes, or parametric tests. *Biometrical Journal* 53(6):894–913.
* Xi D, Brockmann E (2024). *graphicalMCP: Graphical Multiple Comparison Procedures*. R package. (`graph_test_shortcut_fast()` is the reference for precomputing intersection critical values into a lookup table.)
* *gMCPLite* (R package): Independent implementation (`gMCPLite::graphTest()`) used as the external test oracle in `test-RcppExports.R`.
