# The Split–Merge Mutation Algorithm in `multigrain`

*Developer documentation for `.make_split_merge_mutation()` and its helpers ([`R/mutation_helpers.R`](https://github.com/GSK-Biostatistics/multigrain/blob/main/R/mutation_helpers.R)). Written for statistical software engineers familiar with multiple testing and simulation.*

**Contents**

1. [Overview](#1-overview)
2. [Background: Simplexes, Parameter Groups, and Encoding](#2-background-simplexes-parameter-groups-and-encoding)
3. [How the Split–Merge Algorithm Works](#3-how-the-splitmerge-algorithm-works)
4. [Worked Example ($m = 4$)](#4-worked-example-m--4)
5. [Implementation Details in `R/mutation_helpers.R`](#5-implementation-details-in-rmutation_helpersr)
6. [Testing, Notation Summary, and References](#6-testing-notation-summary-and-references)

---

## 1. Overview

When designing a graphical multiple-testing procedure (Bretz et al., 2009), `graph_optimise()` searches for an initial hypothesis weight vector $\mathbf{w}$ and a transition matrix $\mathbf{G}$ that maximise a user-defined trial-success criterion across simulated clinical trials (Spiers et al., 2026).

Because the simulated objective function is a step function of the rejection thresholds—and can have multiple local optima—`graph_optimise()` begins with a global search using a genetic algorithm (`GA::ga()`) before polishing the best candidate graph with local derivative-free optimisation (`nloptr::nloptr()` using COBYLA) and pruning redundant weights and edges (`prune_graph()`).

Inside the genetic algorithm, `multigrain` maintains a population of candidate graphs and evolves them across generations:

* Crossover between different graphs is disabled (`pcrossover = 0`), and every member of the population is mutated in every generation (`pmutation = 1`).
* Elitism preserves the best candidate graphs found so far, and periodic Nelder–Mead local search (`optim = TRUE`, `poptim = 0.2`) refines promising candidates.
* Exploration of the space of valid graphs is driven by the custom **split–merge mutation operator** created by `.make_split_merge_mutation()`:

```r
mutation_fn <- multigrain:::.make_split_merge_mutation(
  hyp_constraint   = graph_constraint$hyp_constraint,
  trans_constraint = graph_constraint$trans_constraint,
  p_del            = 0.45,
  scale            = 0.2
)
# Returns a closure function(object, parent) passed to GA::ga(mutation = ...)
```

Think of the mutation operator as the genetic-algorithm analogue of a proposal distribution in Metropolis–Hastings Markov chain Monte Carlo (MCMC): given a parent graph in the current population, it proposes a modified child graph to be evaluated on the simulated $p$-value matrix via `graph_shortcut()`.

At a high level, a graphical procedure consists of up to $m + 1$ independent parameter groups—the vector of initial hypothesis weights $\mathbf{w}$ and the $m$ rows of the transition matrix $\mathbf{G}$—each of which sums to $1$. On each call, the split–merge operator selects one group and applies one of two moves:

1. **A merge move** (with probability $p_{\text{del}} = 0.45$), which picks two active entries in the chosen group and folds the smaller entry into the larger one, setting the smaller entry to zero. This removes an edge or initial weight, simplifying the graph.
2. **A transfer move** (with probability $1 - p_{\text{del}} = 0.55$, or whenever no group has two active entries to merge), which picks one active entry and one other free entry in the same group and shifts a random fraction of their combined mass between them. When the second entry is already positive, this tunes the balance between the two entries; when the second entry is zero, it splits the first entry and activates a previously unused edge or hypothesis weight.

Within the package, the mutation operator interacts with three components:

| Component | Source file | Role |
|---|---|---|
| `.graph_optimise_ga()` | [`R/optimisation.R`](https://github.com/GSK-Biostatistics/multigrain/blob/main/R/optimisation.R) | Instantiates `.make_split_merge_mutation()` and passes it to `GA::ga()` |
| `.make_split_merge_mutation()` and helpers | [`R/mutation_helpers.R`](https://github.com/GSK-Biostatistics/multigrain/blob/main/R/mutation_helpers.R) | Decodes a parent's parameter groups, applies a merge or transfer move, and encodes the child vector |
| `create_obj_func()`, `split_theta()`, `recover_full_weights()`, `recover_full_trans_matrix()` | [`R/objective_function.R`](https://github.com/GSK-Biostatistics/multigrain/blob/main/R/objective_function.R) | Decodes the parameter vector into $(\mathbf{w}, \mathbf{G})$, zeroes out below-threshold entries, and evaluates trial success |

---

## 2. Background: Simplexes, Parameter Groups, and Encoding

### 2.1 Why Optimise on Simplexes? (Exhaustive Graphs and Fewer Decision Variables)

In the general graphical framework of Bretz et al. (2009), the initial weights $\mathbf{w} = (w_1, \dots, w_m)$ and each row $\mathbf{g}_i = (g_{i1}, \dots, g_{im})$ of the transition matrix $\mathbf{G}$ are required only to be non-negative and sum to *at most* $1$:

$$\sum_{j=1}^m w_j \le 1 \qquad \text{and} \qquad \sum_{\substack{k=1 \\ k \ne i}}^m g_{ik} \le 1 \quad \text{for each } i = 1, \dots, m.$$

However, any weight left unallocated is simply discarded: allocating less than $100\%$ of the initial significance level $\alpha$, or transferring less than $100\%$ of a rejected hypothesis's weight to the remaining hypotheses, can never increase the power to reject any hypothesis. During optimisation, `multigrain` therefore restricts attention to *exhaustive* graphs in which the initial weights and every transition row sum to $1$:

$$\sum_{j=1}^m w_j = 1 \qquad \text{and} \qquad \sum_{\substack{k=1 \\ k \ne i}}^m g_{ik} = 1 \quad \text{for each } i = 1, \dots, m.$$

In mathematics, a vector of $k$ non-negative numbers $\mathbf{s} = (s_1, \dots, s_k)$ that sum to $1$ ($s_j \ge 0$ and $\sum_{j=1}^k s_j = 1$) is called a point on the **unit simplex**. Geometrically, for $k = 2$ the unit simplex is the line segment from $(1, 0)$ to $(0, 1)$, and for $k = 3$ it is the triangle with vertices $(1, 0, 0)$, $(0, 1, 0)$, and $(0, 0, 1)$. For practical purposes, "a point on the simplex" simply means a vector of non-negative proportions (or **shares**) that add up to $1$.

Enforcing the sum-to-one constraint on the initial weights and on every transition row has two major benefits for numerical optimisation:

1. **No wasted significance level:** Every candidate graph evaluated by the optimiser uses the full significance level $\alpha$ and redistributes all rejected weight.
2. **Fewer decision variables to optimise:** Whenever $k$ free weights must add up to a known total, choosing the first $k - 1$ weights automatically determines the $k$-th weight as whatever is left over. This removes one decision variable from the initial weight vector and one from each of the $m$ transition rows. For an unconstrained graph with $m$ hypotheses, the number of decision variables drops from $m + m(m - 1) = m^2$ to $(m - 1) + m(m - 2) = m^2 - m - 1$ (for example, from $16$ to $11$ when $m = 4$).

### 2.2 Parameter Groups, Constraints, and the Residual Encoding

A user can fix specific initial weights or transition edges via `graph_constraint()` (for instance, fixing a weight to $0$ to forbid an edge, or fixing $w_1 = 0.5$). In `hyp_constraint` and `trans_constraint`, fixed entries hold their numeric value and free entries are marked as `NA` (with the diagonal of `trans_constraint` always fixed to $0$).

This structure divides a graph into independent **groups** of free entries—one group for the free entries of $\mathbf{w}$, and one group for the free entries of each row of $\mathbf{G}$:

* Suppose a group has $k$ free entries and its fixed entries sum to $f$. Then the **free mass** available to be divided among those $k$ free entries is
  $$c = 1 - f.$$
* If $k = 0$ or $k = 1$ (or if $c = 0$), the group has no degrees of freedom: its entries are completely determined by the constraints and require no optimisation parameters.
* If $k \ge 2$, the optimiser stores only the first $k - 1$ free entries $(p_1, \dots, p_{k-1})$ in the parameter vector $\mathbf{x}$, and the decoder functions `recover_full_weights()` and `recover_full_trans_matrix()` in [`R/objective_function.R`](https://github.com/GSK-Biostatistics/multigrain/blob/main/R/objective_function.R) compute the $k$-th free entry as the **derived (residual) entry**:
  $$p_k = c - \sum_{j=1}^{k-1} p_j.$$

Dividing a group's $k$ free entries $(p_1, \dots, p_k)$ by its free mass $c$ converts them into normalised **shares** on the unit simplex:

$$\mathbf{s} = (s_1, \dots, s_k) = \left(\frac{p_1}{c}, \dots, \frac{p_k}{c}\right), \qquad \text{where } s_j \ge 0 \text{ and } \sum_{j=1}^k s_j = 1.$$

### 2.3 Why Mutating Simplex-Constrained Parameters Requires Care

Because `GA::ga()` passes the flat parameter vector $\mathbf{x}$ (containing $k - 1$ stored numbers per group) to the mutation function, it is tempting to mutate $\mathbf{x}$ directly. Earlier versions of `multigrain` used an independent coordinate-wise Cauchy perturbation on the stored entries of $\mathbf{x}$ (archived in [`dev/archived_code/cauchy_mutation.R`](https://github.com/GSK-Biostatistics/multigrain/blob/main/dev/archived_code/cauchy_mutation.R)). However, mutating parameters subject to sum-to-one constraints introduces four subtle pitfalls that the split–merge algorithm is designed to avoid:

1. **Violating the sum-to-one constraint (invalid child graphs):** If stored parameters $p_1, \dots, p_{k-1}$ in $[0, 1]$ are perturbed independently, their sum $\sum_{j=1}^{k-1} p_j$ frequently exceeds $c$. When that happens, the derived $k$-th entry $p_k = c - \sum_{j=1}^{k-1} p_j$ becomes negative, producing an invalid graph that `create_obj_func()` rejects with a penalty without ever running `graph_shortcut()`.
2. **Positional bias toward the derived $k$-th entry:** Because only $k - 1$ entries are stored explicitly and the $k$-th entry absorbs whatever mass is left over, any operation applied directly to the stored coordinates treats the $k$-th hypothesis differently from the first $k - 1$ hypotheses. For example, setting a stored parameter $p_j$ to zero automatically transfers all of its mass into the derived entry $p_k$. Over many generations, this causes weight to accumulate on whichever hypothesis happens to be listed last in the group, making the optimisation result sensitive to how the user numbered the hypotheses.
3. **Wasted evaluations from boundary clipping (duplicate offspring):** A natural way to preserve the free mass $c$ is to pick two entries $(i, j)$ in a group and shift a random amount of weight between them while keeping their sum $t = s_i + s_j$ fixed. However, if a group already has several zero entries (sparse edges) and out-of-bounds steps are simply clipped back to the interval $[0, 1]$, any pair where one entry is already $0$ sits on the boundary of $[0, 1]$. Half of all symmetric random steps will push further past the boundary, get clipped back to their starting value, and return a child graph that is identical to its parent. Evaluating an unchanged child graph runs a full Monte Carlo simulation of `graph_shortcut()` for no new information (Carvalho Pinto and Doerr, 2018).
4. **Balancing graph sparsity with continuous weight tuning:** Clinical trial graphs often benefit from *sparsity*—setting unhelpful transition edges or starting weights cleanly to zero so that rejected $\alpha$ is concentrated on the most promising pathways. Continuous random perturbations are effective at tuning the exact split between positive weights, but rarely land on exact zero to remove an edge.

The split–merge operator resolves all four issues by:

* **Decoding to simplex shares $\mathbf{s}$ before mutating** (and choosing entries uniformly across all $k$ positions, including the derived $k$-th entry), which eliminates positional bias and guarantees that every child graph satisfies the sum-to-one constraints.
* **Separating edge removal from weight tuning**, using a dedicated **merge** move to zero out entries and a **transfer** move with boundary reflection (rather than clipping) so that no mutation ever returns a copy of its parent.

---

## 3. How the Split–Merge Algorithm Works

The split–merge operator adapts two classic ideas from MCMC on simplexes to evolutionary graph optimisation: dimension-changing **split and merge moves** on simplex weights (Richardson and Green, 1997) and **pairwise mass-redistribution steps** that adjust two simplex entries while holding the remaining entries fixed (Smith, 2014).

### 3.1 Active vs. Inactive Entries

When `create_obj_func()` in [`R/objective_function.R`](https://github.com/GSK-Biostatistics/multigrain/blob/main/R/objective_function.R) evaluates a candidate graph, it zeroes out tiny weights and edges before calling `graph_shortcut()`:

* Initial hypothesis weights below $\tau_w = 10^{-4}$ are set to $0$.
* Transition matrix entries below $\tau_g = 10^{-5}$ are set to $0$.

To stay aligned with what the objective function actually sees, `.sm_active()` classifies a free entry $j$ in a group as **active** if its actual weight $c s_j$ is at or above the group's threshold $\tau$ ($\tau_w = 10^{-4}$ for initial weights, $\tau_g = 10^{-5}$ for transition rows), and **inactive** (empty) otherwise.

### 3.2 Choosing Between Merge and Transfer

On each call to the mutation closure created by `.make_split_merge_mutation()`, the algorithm decodes every mutable group (every group with $k \ge 2$ free entries and free mass $c > 0$) into its vector of simplex shares $\mathbf{s} = (s_1, \dots, s_k)$. It then draws $U \sim \text{Uniform}(0, 1)$:

1. **If $U < p_{\text{del}}$ (default $p_{\text{del}} = 0.45$):** Attempt a **merge move**. The algorithm identifies all mutable groups that currently have at least two active entries. If at least one such group exists, it picks one uniformly at random and merges two of its active entries (Section 3.3).
2. **Otherwise (or if no group has at least two active entries):** Perform a **transfer move**. The algorithm picks any mutable group uniformly at random and transfers weight between two of its entries (Section 3.4).

Only the single chosen group is modified; all other groups in the parent graph are left unchanged.

### 3.3 The Merge Move (`.sm_merge()`)

Within the selected group (which has at least two active entries):

1. Pick two distinct active entries uniformly at random.
2. Identify which of the two chosen entries has the smaller share (`small`) and which has the larger share (`large`, breaking ties in favour of the second entry drawn).
3. Add the smaller entry's share to the larger entry, and set the smaller entry to zero:
   $$s_{\text{large}} \leftarrow s_{\text{large}} + s_{\text{small}}, \qquad s_{\text{small}} \leftarrow 0.$$

Because $s_{\text{large}} + s_{\text{small}}$ is unchanged and all other entries in $\mathbf{s}$ are untouched, the shares still sum to $1$, and the number of active entries in the group decreases by $1$.

**Why fold the smaller entry into the larger?** If a mutation suddenly deleted a large weight (for example, $w_1 = 0.8$), the resulting child graph would behave very differently from its parent and would almost always have lower power, causing the genetic algorithm to reject it. By always absorbing the *smaller* of the two randomly chosen active entries into the larger, the merge move preferentially prunes small or moderate weights, producing a gentler structural step that selection is much more likely to retain.

### 3.4 The Transfer Move (`.sm_transfer()`)

Within the selected group of $k \ge 2$ free entries:

1. **Pick a donor $i$:** Choose $i$ uniformly at random from the group's active entries (or, in the rare case that rounding left no entry above the threshold, choose `which.max(s)`).
2. **Pick a recipient $j$:** Choose $j$ uniformly at random from the remaining $k - 1$ free entries in the group—regardless of whether $j$ is currently active or inactive.
3. **Compute the pair total and current split:** Let
   $$t = s_i + s_j > 0 \qquad \text{and} \qquad u = \frac{s_i}{t} \in (0, 1]$$
   be the combined share of the pair and the fraction of that combined share currently held by the donor $i$.
4. **Perturb the split with a reflected Cauchy step:** Draw a standard Cauchy random variable $C \sim \text{Cauchy}(0, 1)$, perturb $u$ by `scale * C` (with default `scale` $\sigma = 0.2$), and reflect the result into the unit interval via `.reflect01()`:
   $$u_{\text{new}} = \operatorname{reflect}_{[0,1]}(u + \sigma C).$$
   (If $u_{\text{new}}$ is non-finite or lands on exact floating-point `0` or `1`, redraw $C$ until $0 < u_{\text{new}} < 1$.)
5. **Update the pair's shares:**
   $$s_i \leftarrow t \, u_{\text{new}}, \qquad s_j \leftarrow t - s_i.$$

A Cauchy perturbation is used rather than a Gaussian perturbation because the Cauchy distribution has a sharp peak near zero and heavy tails: most draws make small local adjustments to $u$, while occasional larger draws jump across the interval $(0, 1)$ to explore distant weight allocations.

Depending on whether the randomly chosen recipient $j$ was already active or inactive, the transfer move performs one of two roles:

* **When the recipient $j$ is already active ($s_j > 0$, so $0 < u < 1$): Weight tuning.** Both $s_i$ and $s_j$ start positive and end strictly positive ($0 < u_{\text{new}} < 1$). The move adjusts the relative balance between entries $i$ and $j$ without deleting either entry or changing any other entry in the group.
* **When the recipient $j$ is inactive ($s_j = 0$, so $u = 1$): Splitting (activating an edge).** Because $s_j = 0$, the donor starts with $100\%$ of the pair's mass ($u = 1$). After reflection into $(0, 1)$, the donor retains fraction $u_{\text{new}} \in (0, 1)$ of its original mass and passes the remaining fraction $1 - u_{\text{new}} > 0$ to the previously empty recipient $j$. This splits a single weight into two positive weights, bringing a dormant edge or hypothesis weight back into the graph.

### 3.5 Why Reflection (`.reflect01()`) Instead of Clipping?

In step 4 of a transfer move, the raw step $v = u + \sigma C$ can fall outside $[0, 1]$ (either $v < 0$ or $v > 1$). Rather than clipping $v$ to $0$ or $1$, `.reflect01(v)` folds $v$ back and forth between the walls at $0$ and $1$ like a billiard ball bouncing between two cushions:

```r
.reflect01 <- function(v) {
  v <- v %% 2
  ifelse(v > 1, 2 - v, v)
}
```

To see how this works:

* **Inside $[0, 1]$:** If $v = 0.65$, `v %% 2` is $0.65 \le 1$, so `.reflect01(0.65) = 0.65` (unchanged).
* **Just above $1$:** If $v = 1.25$ (overshooting $1$ by $0.25$), `v %% 2` is $1.25 > 1$, so $2 - 1.25 = 0.75$ (bounced $0.25$ inward from $1$).
* **Just below $0$:** If $v = -0.25$ (overshooting $0$ by $0.25$), R's modulo operator returns `(-0.25) %% 2 = 1.75 > 1`, so $2 - 1.75 = 0.25$ (bounced $0.25$ inward from $0$).
* **Larger steps:** Any real number $v$ is mapped periodically onto $[0, 1]$ by taking $v \bmod 2$ and reflecting values in $(1, 2]$ across $1$.

Using reflection instead of clipping achieves two properties simultaneously:

1. **Zero duplicate offspring when splitting ($u = 1$):** When the recipient $j$ is empty ($u = 1$), a positive Cauchy draw $C > 0$ gives $1 + \sigma C > 1$. Clipping would clamp this back to $1$, leaving $s_i = t$ and $s_j = 0$ unchanged ($50\%$ of all draws wasted). Reflection bounces $1 + \sigma C$ inward to $1 - \sigma C$ (for small steps), so positive and negative Cauchy draws both move mass from the donor $i$ into the empty recipient $j$.
2. **Transfers never accidentally delete an entry:** Because $u + \sigma C$ is a continuous random variable, the reflected value $u_{\text{new}}$ lands in the open interior $(0, 1)$ with probability $1$ (and the `repeat` loop guards against rare floating-point underflow to exact `0` or `1`). Consequently, a transfer never zeroes out either entry; deletions occur *only* through the explicit merge move.

### 3.6 How Merge and Transfer Balance Graph Sparsity

Because deletions happen only via merges and activations happen only when a transfer picks an inactive recipient, the balance between sparse and dense graphs is governed directly by the merge probability $p_{\text{del}}$.

Consider a group with $k$ free entries, of which $k_{\text{act}}$ are currently active ($2 \le k_{\text{act}} < k$):

* With probability $p_{\text{del}}$, a **merge** reduces the number of active entries by $1$ ($\Delta k_{\text{act}} = -1$).
* With probability $1 - p_{\text{del}}$, a **transfer** picks one active donor and chooses a recipient uniformly among the remaining $k - 1$ entries. Since $k - k_{\text{act}}$ of those $k - 1$ entries are inactive, the transfer activates a new entry ($\Delta k_{\text{act}} = +1$) with probability $(k - k_{\text{act}}) / (k - 1)$, and tunes two already-active entries ($\Delta k_{\text{act}} = 0$) otherwise.

The expected change in the number of active entries on a single mutation in that group (in the absence of selection) is therefore:

$$\mathbb{E}[\Delta k_{\text{act}}] = (1 - p_{\text{del}}) \frac{k - k_{\text{act}}}{k - 1} - p_{\text{del}}.$$

Setting $\mathbb{E}[\Delta k_{\text{act}}] = 0$ gives the natural equilibrium number of active entries $k_{\text{act}}^*$ where activations and deletions balance out:

$$k_{\text{act}}^* = k - (k - 1) \frac{p_{\text{del}}}{1 - p_{\text{del}}}.$$

At `multigrain`'s default `p_del = 0.45` (where $\frac{p_{\text{del}}}{1 - p_{\text{del}}} = \frac{9}{11} \approx 0.818$), an unconstrained transition row in a $10$-hypothesis graph ($k = 9$ off-diagonal entries) has an equilibrium density of about $k_{\text{act}}^* \approx 2.5$ to $3.6$ active outgoing edges per row (accounting for the floor at $k_{\text{act}} = 1$, where merges cannot occur). This keeps the population exploring moderately sparse graphs while still proposing plenty of weight-tuning moves ($55\%$ of calls).

---

## 4. Worked Example ($m = 4$)

To see how `.make_split_merge_mutation()` operates on a concrete graph, consider $m = 4$ hypotheses with the constrained layout used in [`tests/testthat/test-mutation_helpers.R`](https://github.com/GSK-Biostatistics/multigrain/blob/main/tests/testthat/test-mutation_helpers.R):

* **Hypothesis weight constraints (`hyp_constraint`):** `(NA, NA, 0.2, NA)`. Here $w_3 = 0.2$ is fixed, leaving $k = 3$ free weights ($w_1, w_2, w_4$) with free mass $c = 1 - 0.2 = 0.8$.
* **Transition matrix constraints (`trans_constraint`):**
  $$\mathbf{G}_{\text{constraint}} = \begin{pmatrix} 0 & \text{NA} & \text{NA} & \text{NA} \\ 1 & 0 & 0 & 0 \\ \text{NA} & \text{NA} & 0 & 0.3 \\ \text{NA} & \text{NA} & \text{NA} & 0 \end{pmatrix}.$$
  Here row 2 is completely fixed ($g_{21} = 1$), row 3 has one fixed off-diagonal entry ($g_{34} = 0.3$) and two free entries ($g_{31}, g_{32}$) sharing $c = 1 - 0.3 = 0.7$, and rows 1 and 4 each have three free off-diagonal entries sharing $c = 1$.

### 4.1 Extracting the Mutable Groups (`.sm_groups()`)

Calling `.sm_groups(hyp_constraint, trans_constraint)` scans the constraints in the same order as the decoder (`split_theta()`, `recover_full_weights()`, and `recover_full_trans_matrix()`) and identifies four mutable groups occupying $7$ positions in the encoded parameter vector $\mathbf{x} \in \mathbb{R}^7$:

| Group | Free entries ($k$) | Stored positions in $\mathbf{x}$ (`latent`) | Derived $k$-th entry | Free mass $c$ | Active threshold $\tau$ |
|---|---:|---|---|---:|---:|
| 1. Initial weights $\mathbf{w}$ | $w_1, w_2, w_4$ ($k = 3$) | `1:2` ($w_1, w_2$) | $w_4 = 0.8 - (w_1 + w_2)$ | $0.8$ | $10^{-4}$ |
| 2. Transition row 1 | $g_{12}, g_{13}, g_{14}$ ($k = 3$) | `3:4` ($g_{12}, g_{13}$) | $g_{14} = 1 - (g_{12} + g_{13})$ | $1.0$ | $10^{-5}$ |
| *(Transition row 2)* | *(none, $k = 0$)* | *(none)* | *(none)* | $0.0$ | — |
| 3. Transition row 3 | $g_{31}, g_{32}$ ($k = 2$) | `5` ($g_{31}$) | $g_{32} = 0.7 - g_{31}$ | $0.7$ | $10^{-5}$ |
| 4. Transition row 4 | $g_{41}, g_{42}, g_{43}$ ($k = 3$) | `6:7` ($g_{41}, g_{42}$) | $g_{43} = 1 - (g_{41} + g_{42})$ | $1.0$ | $10^{-5}$ |

Suppose the parent individual has parameter vector

$$\mathbf{x} = (0.40,\ 0.16,\ \ 0.50,\ 0.00,\ \ 0.35,\ \ 0.20,\ 0.50).$$

### 4.2 Tracing a Merge Move on Group 1 (Initial Weights)

Suppose the mutation operator selects a **merge move** ($U < 0.45$) and chooses **Group 1** (the initial weights $w_1, w_2, w_4$, with $c = 0.8$ and `latent = 1:2`):

1. **Decode to simplex shares (`.sm_to_shares()`):**
   * Stored entries: $p_1 = x_1 = 0.40$ ($w_1$) and $p_2 = x_2 = 0.16$ ($w_2$).
   * Derived entry: $p_3 = c - (p_1 + p_2) = 0.80 - 0.56 = 0.24$ ($w_4$).
   * Divide by $c = 0.8$ to obtain the simplex shares:
     $$\mathbf{s} = \left(\frac{0.40}{0.80},\ \frac{0.16}{0.80},\ \frac{0.24}{0.80}\right) = (0.50,\ 0.20,\ 0.30).$$
   * All three entries satisfy $c s_j \ge 10^{-4}$, so all three positions $\{1, 2, 3\}$ (corresponding to $w_1, w_2, w_4$) are active and equally likely to be chosen—even though $w_4$ is not stored directly in $\mathbf{x}$.
2. **Merge two active entries (`.sm_merge()`):**
   * Suppose the random pair drawn from $\{1, 2, 3\}$ is $(2, 3)$ (corresponding to $w_2$ with share $s_2 = 0.20$ and $w_4$ with share $s_3 = 0.30$).
   * Because $s_2 = 0.20 < s_3 = 0.30$, entry $2$ is `small` and entry $3$ is `large`.
   * Entry $2$'s share is folded into entry $3$:
     $$s_3 \leftarrow 0.30 + 0.20 = 0.50, \qquad s_2 \leftarrow 0, \qquad \mathbf{s}_{\text{new}} = (0.50,\ 0.00,\ 0.50).$$
3. **Re-encode into $\mathbf{x}$ (`.sm_from_shares()`):**
   * Multiply the first $k - 1 = 2$ shares by $c = 0.8$:
     $$(p_1, p_2) = 0.8 \times (0.50, 0.00) = (0.40, 0.00).$$
   * Write these back into positions `1:2` of $\mathbf{x}$:
     $$\mathbf{x}_{\text{child}} = (0.40,\ 0.00,\ \ 0.50,\ 0.00,\ \ 0.35,\ \ 0.20,\ 0.50).$$
   * When decoded by `recover_full_weights()`, the child graph has initial weights $\mathbf{w} = (0.40,\ 0.00,\ 0.20,\ 0.40)$: hypothesis $H_2$'s initial weight has been cleanly set to zero and absorbed into $H_4$.

### 4.3 Tracing a Transfer Move on Group 2 (Transition Row 1)

Now suppose instead that the operator performs a **transfer move** on **Group 2** (transition row 1: free entries $g_{12}, g_{13}, g_{14}$ with $c = 1.0$ and `latent = 3:4`):

1. **Decode to simplex shares:**
   * Stored entries: $p_1 = x_3 = 0.50$ ($g_{12}$) and $p_2 = x_4 = 0.00$ ($g_{13}$).
   * Derived entry: $p_3 = 1.0 - (0.50 + 0.00) = 0.50$ ($g_{14}$).
   * Shares: $\mathbf{s} = (0.50,\ 0.00,\ 0.50)$. Here entries $1$ ($g_{12}$) and $3$ ($g_{14}$) are active, while entry $2$ ($g_{13}$) is inactive ($0$).
2. **Select donor and recipient:**
   * Donor $i$ is drawn uniformly from the active set $\{1, 3\}$; suppose $i = 1$ ($g_{12}$, with $s_1 = 0.50$).
   * Recipient $j$ is drawn uniformly from the remaining entries $\{2, 3\}$:
     * **Case A ($j = 3$, active recipient $g_{14}$ — Tuning):** The pair total is $t = s_1 + s_3 = 0.50 + 0.50 = 1.00$, and the donor's current fraction is $u = 0.50 / 1.00 = 0.50$. If the Cauchy draw is $C = -0.5$, then $u + 0.2 C = 0.40$, which lies inside $(0, 1)$ so $u_{\text{new}} = 0.40$. The updated shares are $s_1 = 0.40$ and $s_3 = 0.60$, giving $\mathbf{s}_{\text{new}} = (0.40,\ 0.00,\ 0.60)$ (tuning the split between $g_{12}$ and $g_{14}$ while leaving $g_{13} = 0$).
     * **Case B ($j = 2$, inactive recipient $g_{13}$ — Splitting / Activating):** The pair total is $t = s_1 + s_2 = 0.50 + 0.00 = 0.50$, and the donor starts with $100\%$ of the pair's mass ($u = 0.50 / 0.50 = 1.00$). Suppose the Cauchy draw is $C = +0.5$, pushing outward to $v = 1.00 + 0.2(0.5) = 1.10$. Instead of clipping back to $1.00$, `.reflect01(1.10)` reflects $1.10$ across $1$ to $u_{\text{new}} = 2 - 1.10 = 0.90$. The donor retains $s_1 = 0.50 \times 0.90 = 0.45$ and the previously empty recipient receives $s_2 = 0.50 - 0.45 = 0.05$, giving $\mathbf{s}_{\text{new}} = (0.45,\ 0.05,\ 0.50)$ and activating edge $H_1 \to H_3$.

---

## 5. Implementation Details in `R/mutation_helpers.R`

The implementation in [`R/mutation_helpers.R`](https://github.com/GSK-Biostatistics/multigrain/blob/main/R/mutation_helpers.R) is compact (around 200 lines of R code) and organised into small, self-contained helper functions:

| Function | Role |
|---|---|
| `.make_split_merge_mutation(hyp_constraint, trans_constraint, p_del = 0.45, scale = 0.2)` | Validates arguments, precomputes `groups` once via `.sm_groups()`, and returns the `function(object, parent)` closure required by `GA::ga()`. |
| `.sm_groups(hyp_constraint, trans_constraint)` | Walks `hyp_constraint` and the rows of `trans_constraint`, recording `list(latent, c, threshold)` for every group with $k \ge 2$ free entries and filtering to those with positive free mass $c > 0$. |
| `.sm_free_entries(x, g)` | Extracts the $k - 1$ stored parameters `p <- x[g$latent]` and appends the derived $k$-th entry `g$c - sum(p)`. |
| `.sm_to_shares(x, g)` | Divides `.sm_free_entries(x, g)` by `g$c` to obtain simplex shares $\mathbf{s}$, repairing any out-of-bounds vector via `.simplex_project()`. |
| `.simplex_project(v)` | Computes the exact Euclidean projection of any vector $\mathbf{v} \in \mathbb{R}^k$ onto the unit simplex. |
| `.sm_active(g, s)` | Returns the indices of `s` where `g$c * s >= g$threshold`. |
| `.sm_merge(g, s)` | Performs a merge move on shares `s` (Section 3.3). |
| `.sm_transfer(g, s, scale)` | Performs a transfer move on shares `s` using `.reflect01()` (Section 3.4). |
| `.reflect01(v)` | Folds any real number `v` into $[0, 1]$ via `v %% 2` and `ifelse(v > 1, 2 - v, v)` (Section 3.5). |
| `.sm_from_shares(x, g, s)` | Multiplies the first $k - 1$ shares by `g$c`, applies a floating-point safety margin, and writes them back into `x[g$latent]`. |

Two defensive numerical details in `R/mutation_helpers.R` are worth noting when reading or maintaining the code:

### 5.1 Repairing Out-of-Bounds Parents (`.simplex_project()`)

Every child produced by `.sm_merge()` or `.sm_transfer()` lies on the simplex. However, `GA::ga()` is configured with `optim = TRUE`, which periodically runs `stats::optim(..., method = "Nelder-Mead")` on individuals in the population. Because Nelder–Mead operates in the unconstrained box $[0, 1]^d$, it can leave a group's stored parameters summing to slightly more than $c$, so the derived $k$-th entry in `.sm_free_entries(x, g)` is negative when that individual is later selected as a parent.

Whenever `.sm_to_shares(x, g)` detects a negative entry (`any(v < 0)`), it repairs the group before mutating by projecting $\mathbf{v}$ onto the unit simplex via `.simplex_project(v)`:

```r
.simplex_project <- function(v) {
  u <- sort(v, decreasing = TRUE)
  rho <- max(which(u - (cumsum(u) - 1) / seq_along(u) > 0))
  pmax(v - (sum(u[seq_len(rho)]) - 1) / rho, 0)
}
```

This is the standard $O(k \log k)$ Euclidean simplex projection algorithm (Duchi et al., 2008): geometrically, it shifts all entries of $\mathbf{v}$ down (or up) by a single constant $\theta = (\sum_{j=1}^\rho u_j - 1) / \rho$ and clamps negative results to zero (`pmax(v - theta, 0)`), where $\rho$ is chosen so that the remaining positive entries sum to $1$. If $\mathbf{v}$ is already on the unit simplex, `.simplex_project(v)` returns $\mathbf{v}$ unchanged.

### 5.2 Preventing Negative Derived Entries from Floating-Point Rounding (`.sm_from_shares()`)

In exact arithmetic, if $\sum_{j=1}^k s_j = 1$ with $s_k \ge 0$, then the stored parameters $p_j = c s_j$ ($j = 1, \dots, k - 1$) satisfy $\sum_{j=1}^{k-1} p_j \le c$, so the derived $k$-th entry is non-negative.

In IEEE 754 double-precision floating-point arithmetic, however, two things can happen when $s_k = 0$ (i.e. when the derived $k$-th entry is empty):

1. Rounding each product `p <- g$c * s[-length(s)]` up by half a unit in the last place (ulp) can make `sum(p)` exceed `g$c` by $1\text{--}2$ ulps (roughly $10^{-16}$).
2. Even if `sum(p) <= g$c`, `recover_full_weights()` and `recover_full_trans_matrix()` compute the derived entry as `1 - sum(fixed_and_stored_entries)`. Because floating-point addition is not associative, summing the fixed entries and the stored parameters `p` together can still round to $1 + 10^{-16}$, leaving the decoded $k$-th entry at $-10^{-16}$. Inside `create_obj_func()`, the strict check `any(hyp_weight < 0)` or `any(trans_matrix < 0)` would then reject the child graph as invalid.

To prevent floating-point rounding from ever producing a negative derived entry, `.sm_from_shares()` keeps the sum of the stored parameters inside `g$c * (1 - 1e-14)` by shaving any tiny excess off the largest stored parameter in the group:

```r
.sm_from_shares <- function(x, g, s) {
  p <- g$c * s[-length(s)]
  for (attempt in 1:3) {
    excess <- sum(p) - g$c * (1 - 1e-14)
    if (excess <= 0) {
      break
    }
    j <- which.max(p)
    p[j] <- max(p[j] - excess, 0)
  }
  x[g$latent] <- pmin(pmax(p, 0), 1)
  x
}
```

When the derived $k$-th entry was supposed to be zero ($s_k = 0$), this leaves it with a tiny positive residual of about $10^{-14} c$. Because $10^{-14} c$ is at least nine orders of magnitude below the active truncation thresholds ($\tau_w = 10^{-4}$ and $\tau_g = 10^{-5}$), `create_obj_func()` immediately truncates it to exact `0.0` after passing the non-negativity check, and `.sm_active()` continues to treat it as inactive on subsequent mutations.

### 5.3 Tracking Latent Indices Across Zero-Mass Groups (`.sm_groups()`)

A subtle edge case occurs if a user specifies constraints where a row has $k \ge 2$ entries marked `NA`, but the fixed entries in that same row already sum to $1$ (so $c = 0$). Even though $c = 0$ means those `NA` entries must all be $0$ and cannot be mutated, the decoder `recover_full_trans_matrix()` still allocates $k - 1$ positions in the parameter vector $\mathbf{x}$ for that row.

For that reason, `.sm_groups()` first advances the cursor `pos <- pos + k - 1L` for *every* group with $k \ge 2$ free entries (so that the latent indices of all subsequent rows remain aligned with `recover_full_trans_matrix()`), and only filters out groups with $c \le 0$ at the very end via `Filter(function(g) g$c > 0, groups)`.

---

## 6. Testing, Notation Summary, and References

### 6.1 Unit Test Suite

All invariants of the split–merge operator are verified in [`tests/testthat/test-mutation_helpers.R`](https://github.com/GSK-Biostatistics/multigrain/blob/main/tests/testthat/test-mutation_helpers.R):

| Test block in `test-mutation_helpers.R` | Property verified |
|---|---|
| `.reflect01` tests | Maps arbitrary positive and negative values into $[0, 1]$, acts as the identity on $(0, 1)$, and returns `NaN` for `Inf`. |
| `.simplex_project` tests | Projects out-of-bounds vectors onto the unit simplex ($p_j \ge 0$, $\sum p_j = 1$) and leaves valid simplex points unchanged. |
| `.sm_groups` tests | Matches the latent parameter layout of `split_theta()`, `recover_full_weights()`, and `recover_full_trans_matrix()`, including when zero-mass groups ($c = 0, k \ge 2$) are skipped. |
| `.sm_merge` and `.sm_transfer` tests | Merge zeroes out the smaller of two active entries while preserving the group sum; transfer preserves the pair sum and changes active support only when activating an empty recipient. |
| Graph validity and single-group mutation | Across 300 random parents on both unconstrained (`free_4`) and constrained (`pinned_h`, `pinned_t`) graphs, every child decodes to a valid graph (`is_graph_valid() == TRUE`), preserves fixed entries, and modifies exactly one group. |
| Zero duplicate-parent rate | Across 1,000 mutations at `p_del` $\in \{0, 0.45, 1\}$, no child vector equals its parent. |
| Positional symmetry | On a symmetric 6-hypothesis group, all 6 positions—including the 6th derived entry—are mutated at equal rates under both pure transfer (`p_del = 0`) and pure merge (`p_del = 1`). |
| Invalid-parent repair | When a parent has stored parameters summing to more than $c$, mutating a group repairs that group onto the simplex. |

### 6.2 Notation and Function Summary

| Symbol / Identifier | Description |
|---|---|
| $m$ | Number of hypotheses $H_1, \dots, H_m$ |
| $\mathbf{w} = (w_1, \dots, w_m)$ | Initial hypothesis weight vector ($w_j \ge 0$, $\sum_{j=1}^m w_j = 1$) |
| $\mathbf{G} = (g_{ik})$ | $m \times m$ transition matrix ($g_{ii} = 0$, $g_{ik} \ge 0$, $\sum_{k \ne i} g_{ik} = 1$) |
| $k$ | Number of free (`NA`) entries in a single group (the initial weights or one transition row) |
| $c = 1 - \sum_{\text{fixed}} v_j$ | Free mass ($1$ minus the sum of fixed entries) shared by the $k$ free entries of a group |
| $\mathbf{s} = (s_1, \dots, s_k)$ | Decoded shares on the unit simplex ($s_j \ge 0$, $\sum_{j=1}^k s_j = 1$), where each free entry equals $c s_j$ |
| $\tau$ (`threshold`) | Minimum weight for an entry to be treated as active ($\tau_w = 10^{-4}$ for initial weights, $\tau_g = 10^{-5}$ for transition entries) |
| `p_del` ($0.45$) | Probability of attempting a merge move on a mutation call |
| `scale` ($\sigma = 0.2$) | Scale parameter of the standard Cauchy perturbation in a transfer move |

### 6.3 References

* Bretz F, Maurer W, Brannath W, Posch M (2009). A graphical approach to sequentially rejective multiple test procedures. *Statistics in Medicine* 28(4):586–604. (Introduces the graphical multiple-testing framework parameterised by $\mathbf{w}$ and $\mathbf{G}$.)
* Spiers A D V, Grayling M J, Wheeler G M, Mander A P (2026). Gain-function optimisation of graphical multiple testing procedures for confirmatory clinical trials. *arXiv:2609.19994v1*. (Describes the `multigrain` optimisation pipeline combining a genetic algorithm, local COBYLA refinement, and post-optimisation pruning.)
* Richardson S, Green P J (1997). On Bayesian analysis of mixtures with an unknown number of components. *Journal of the Royal Statistical Society: Series B* 59(4):731–792. (Introduces split and merge moves on simplex weights that change the number of active components while conserving their sum.)
* Smith A (2014). A Gibbs sampler on the $n$-simplex. *Annals of Applied Probability* 24(1):114–130. (Analyses pairwise moves on the unit simplex that redistribute mass between two randomly chosen coordinates while holding the remaining coordinates fixed.)
* Carvalho Pinto E, Doerr C (2018). Towards a more practice-aware runtime analysis of evolutionary algorithms. *arXiv:1812.00493*. (Demonstrates the efficiency gain from mutation operators that prevent offspring identical to their parent.)
* Duchi J, Shalev-Shwartz S, Singer Y, Chandra T (2008). Efficient projections onto the $\ell_1$-ball for learning in high dimensions. *Proceedings of the 25th International Conference on Machine Learning (ICML)*, 272–279. (Gives the $O(k \log k)$ sorting algorithm used by `.simplex_project()` to project an out-of-bounds vector onto the unit simplex.)
