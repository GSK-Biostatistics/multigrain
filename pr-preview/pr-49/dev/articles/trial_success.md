# Trial success measures

## Gain functions

A graphical test (see [Getting
started](https://gsk-biostatistics.github.io/multigrain/dev/articles/get-started.md))
applies a graph $(\mathbf{w},\mathbf{G})$ to simulated p-values and
produces a rejection pattern $\mathbf{r} = (r_{1},\ldots,r_{m})$ for
each trial, where $r_{i} = 1$ if $H_{i}$ is rejected. A **gain
function** $\psi(\mathbf{r})$ assigns a numerical score to each pattern
— average power, disjunctive power, a regulatory hurdle, or any
composite you design.
[`graph_optimise()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/graph_optimisation.md)
searches for the graph that maximises $E\lbrack\psi(\mathbf{r})\rbrack$.

Because FWER control holds for every valid graph, the choice among them
is free to be made on value grounds. The gain is a score on whatever
scale you choose; it does not need to be a probability.

### Fixed-sample designs

[`trial_success()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/trial_success.md)
takes an expression in `r1, r2, ..., rm` and compiles it to C++ so it
can be evaluated quickly over thousands of simulated trials.

``` r
set.seed(1)
rej <- matrix(runif(4000) < 0.6, ncol = 4)   # 1000 trials, m = 4

# Expected number of rejections (= m × average power)
trial_success(r1 + r2 + r3 + r4)$func(rej)
#> ✔ Trial success function compiled and sourced successfully.
#> [1] 2.445

# Disjunctive power: at least one rejection
trial_success(r1 || r2 || r3 || r4)$func(rej)
#> ✔ Trial success function compiled and sourced successfully.
#> [1] 0.978

# Conjunctive power: all four rejected
trial_success(r1 && r2 && r3 && r4)$func(rej)
#> ✔ Trial success function compiled and sourced successfully.
#> [1] 0.145

# Dual-endpoint hurdle: success for at least one dose
# H1, H2 = high dose (primary, secondary); H3, H4 = low dose
trial_success((r1 && r2) || (r3 && r4))$func(rej)
#> ✔ Trial success function compiled and sourced successfully.
#> [1] 0.596

# Gated secondary with incremental value
trial_success(r1 * (1 + 0.5 * r2))$func(rej)
#> ✔ Trial success function compiled and sourced successfully.
#> [1] 0.7875

# Hurdle plus weighted incremental claims
trial_success((r1 && r2) * (1 + 0.3 * r3 + 0.2 * r4))$func(rej)
#> ✔ Trial success function compiled and sourced successfully.
#> [1] 0.4749
```

The first example returns approximately 2.4, not 0.6. That is correct —
it counts the expected *number* of rejections across four hypotheses,
not a probability. Divide by $m$ to get average power. The gain function
is on whatever scale you define; it does not need to lie in
$\lbrack 0,1\rbrack$.

To inject values from your R session into the expression, use rlang’s
\[`!!`\]\[rlang::injection-operator\] operator:

``` r
w <- 2
trial_success(!!w * r1 + r2 + r3)
#> ✔ Trial success function compiled and sourced successfully.
#> <multigrain_trial_success>
#> 2 * r1 + r2 + r3
```

> **Always parenthesise around `&&` and `||`**
>
> The internal parser does not follow R’s standard operator precedence
> for `&&` and `||`. Write `(r1 && r2) || (r3 && r4)` rather than
> `r1 && r2 || r3 && r4` to ensure the expression is interpreted as you
> intend.

### Group sequential designs

In a group sequential design you test each hypothesis at $K$ planned
analyses (e.g., one interim and one final). The outcome of a simulated
trial is not just *which* hypotheses were rejected but *when* each
rejection was declared. An early rejection — allowing a programme to
proceed sooner — is typically worth more than a late one.
[`trial_success_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/trial_success_gsd.md)
lets you express this.

#### Decision times

[`trial_success_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/trial_success_gsd.md)
uses two families of symbols:

- `t1, t2, ..., tm`: the **decision time** for each hypothesis — the
  analysis (an integer from `1` to `K`) at which the procedure declared
  the rejection, or **`0` if the hypothesis was never rejected**.
- `r1, r2, ..., rm`: the usual rejection indicators. These are derived
  as `ri = (ti > 0)` inside the compiled function, so `r` and `t` are
  always consistent. You can mix both in one expression.

Decision time is the analysis at which the *procedure declared* the
rejection. This is not necessarily the analysis at which the data first
crossed a boundary. Alpha recycling can delay the declaration. For
example, suppose PFS data are fully mature at analysis 1 but PFS carries
zero initial weight ($w_{\text{PFS}} = 0$). PFS cannot reject at
analysis 1 because it has no $\alpha$ to spend. At analysis 2, OS
rejects and recycles its entire level to PFS. PFS then rejects — on its
analysis-1 data — but the decision time is $\tau_{\text{PFS}} = 2$, not
1.

#### Discount tables

You assign values to declaration analyses by passing one or more named
numeric vectors through `...`. For $K = 2$ analyses:

``` r
d = c(1, 0.75)   # full value at analysis 1, 75% at analysis 2
```

You supply only the $K$ values for analyses 1 through $K$. The value for
“never rejected” ($t_{i} = 0$) is **added automatically** as zero — you
do not include it yourself.

In the `objective` expression, apply a table to a decision time as
`d(t1)`. When $t_{1} = 0$ the result is 0; when $t_{1} = 1$ the result
is `d[1]`; when $t_{1} = 2$ the result is `d[2]`.

The same table can be shared across hypotheses, as in `d(t1) + d(t2)`,
or each hypothesis can use a different table, as in
`d_pfs(t1) + d_os(t2)`. All tables must contain $K$ values, but the
values may differ.

> **Discounts are indexed by analysis, not information fraction**
>
> A discount vector maps the shared analysis indices $1,\ldots,K$ to
> values. It does not look up or interpolate from the numeric
> information fraction. Information fractions determine the group
> sequential boundaries and may differ by hypothesis; the gain function
> receives only the analysis at which each rejection was declared.

#### Worked example: hypothesis-specific discounts

Consider a trial with two hypotheses — PFS ($H_{1}$) and OS ($H_{2}$) —
and $K = 2$ analyses. PFS is fully mature at analysis 1, while OS has
70% information at analysis 1 and reaches full information at analysis
2:

| Analysis | PFS information fraction | PFS discount | OS information fraction | OS discount |
|----------|-------------------------:|-------------:|------------------------:|------------:|
| 1        |                      1.0 |         1.00 |                     0.7 |        1.00 |
| 2        |                      1.0 |         0.50 |                     1.0 |        0.75 |

You value PFS at 0.4 and OS at 0.6. A PFS declaration delayed until
analysis 2 retains 50% of its early value, whereas an OS declaration at
analysis 2 retains 75%:

``` r
v_pfs <- 0.4
v_os  <- 0.6
gain <- trial_success_gsd(
    !!v_pfs * d_pfs(t1) + !!v_os * d_os(t2),
    d_pfs = c(1, 0.50),
    d_os  = c(1, 0.75)
)
gain
#> <multigrain_trial_success_gsd/multigrain_trial_success>
#> 0.4 * d_pfs(t1) + 0.6 * d_os(t2)
#> Analyses (K): 2
#> Discount table d_pfs(t): 1.0, 0.5
#> Discount table d_os(t): 1.00, 0.75
```

The gain for every possible combination of decision times is:

|                                     | $t_{2} = 0$ (OS never) | $t_{2} = 1$ (OS at analysis 1) | $t_{2} = 2$ (OS at analysis 2) |
|-------------------------------------|------------------------|--------------------------------|--------------------------------|
| **$t_{1} = 0$ (PFS never)**         | 0.000                  | 0.600                          | 0.450                          |
| **$t_{1} = 1$ (PFS at analysis 1)** | 0.400                  | 1.000                          | 0.850                          |
| **$t_{1} = 2$ (PFS at analysis 2)** | 0.200                  | 0.800                          | 0.650                          |

Reading the table: $(t_{1} = 0,t_{2} = 1)$ means PFS was not rejected
and OS was rejected at analysis 1. The gain is
$0.4 \times d_{PFS}(0) + 0.6 \times d_{OS}(1) = 0 + 0.6 \times 1 = 0.6$.
The entry $(t_{1} = 2,t_{2} = 2)$ — both declared at analysis 2 — gives
$0.4 \times 0.50 + 0.6 \times 0.75 = 0.65$.

For PFS, $t_{1} = 2$ still selects the second discount value even though
its information fraction was already 1 at analysis 1. The discount
reflects the later **declaration**, not a second information-fraction
lookup. This is the alpha-recycling case described above: mature PFS
data can support a decision only after OS releases alpha at the later
analysis.

#### More expressions

You can combine rejection indicators and decision times freely:

``` r
# OS value plus a bonus only if PFS was declared at the interim
trial_success_gsd(r2 + 0.5 * (r1 && t1 == 1), K = 2)

# Hurdle gain with a time-discounted base value
trial_success_gsd((r1 && r2) * d(t2), d = c(1, 0.75))

# Co-primary endpoints declared at the same analysis
trial_success_gsd(
    (t1 == 1 && t2 == 1) + 0.5 * (t1 == 2 && t2 == 2)
)
```

#### End-to-end evaluation

To evaluate the gain on simulated GSD data, you simulate p-values,
transform them to repeated p-values, and evaluate a graph on them with
[`calc_power_pvals_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/calc_power_pvals_gsd.md).
Pass the gain through `custom_power`; the result also carries the power
of each hypothesis by analysis and the distribution of its decision
time.

``` r
alpha <- 0.025
info_frac <- rbind(PFS = c(1, 1), OS = c(0.7, 1))
corr      <- rbind(c(1, 0.5), c(0.5, 1))
pow_nom   <- c(PFS = 1 - pnorm(qnorm(1 - alpha) - 4.01),
               OS  = 1 - pnorm(qnorm(1 - alpha) - 3.44))

set.seed(20260918)
sim <- simulate_pvalues_gsd(
    power_nominal = pow_nom, alpha = alpha,
    corr_matrix = corr, info_frac = info_frac,
    nsim = 20000L
)

ldof <- function(a, t) gsDesign::sfLDOF(a, t)$spend
pv <- transform_pvalues_gsd(
    sim, info_frac = info_frac,
    spending = list(function(a, t) a * (t >= 1), ldof),
    alpha = alpha, grid_size = 1024L
)

G <- rbind(c(0, 1), c(1, 0))   # full recycling between PFS and OS
power <- calc_power_pvals_gsd(
    pv,
    hyp_weight = c(0.5, 0.5),
    trans_matrix = G,
    custom_power = list(gain = gain)
)

round(power$gain, 3)
#> [1] 0.903

round(power$time_distribution, 3)
#>    never analysis 1 analysis 2
#> H1 0.028      0.966      0.006
#> H2 0.078      0.667      0.255
```

The gain can be checked by hand from the decision-time distribution:
$0.4 \times (0.966 + 0.50 \times 0.006) + 0.6 \times (0.667 + 0.75 \times 0.255) \approx 0.903$.

#### Stopping for efficacy

The gains above are evaluated on a trial that always continues: the
kernel analyses every planned look of every simulated trial, and stops
early only once every hypothesis has been rejected. If the trial will in
fact stop for efficacy at an interim analysis, a hypothesis can then be
counted as rejected at a look the trial would never have reached, so
late rejections, and every power that counts them, are overstated. In a
two-hypothesis design that stops once $H_{1}$ is rejected, a graph
reported at a gain of 0.708 with $H_{2}$ power 0.754 was worth 0.663
with $H_{2}$ power 0.614 once the stop was taken into account (section
7.4 of `dev/review/gsd_user_testing/README.md`).

A stopping rule is its own object, created with
[`stopping_rule_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/stopping_rule_gsd.md)
and passed to the optimiser and the power calculation rather than
written into the gain. It is a Boolean condition over the rejection
indicators `r1, r2, ..., rm`, evaluated after each interim analysis on
the hypotheses rejected so far. When the condition is `TRUE` the trial
stops at that analysis and nothing is rejected afterwards. The final
analysis always ends the trial, so it has no condition.

``` r
# Stop at an interim analysis once H1 and at least one of H2 and H3 are
# rejected
rule <- stopping_rule_gsd(r1 && (r2 || r3))
rule
#> <multigrain_stopping_rule_gsd>
#> Every interim: r1 && (r2 || r3)
#> Hypotheses referenced: up to r3

# With K = 3 analyses, one condition per interim analysis: stop at the first
# interim once H1 is rejected, at the second only once H2 or H3 is too
rule3 <- stopping_rule_gsd(list(r1, r1 && (r2 || r3)))
rule3
#> <multigrain_stopping_rule_gsd>
#> Interim 1: r1
#> Interim 2: r1 && (r2 || r3)
#> Hypotheses referenced: up to r3
```

A single condition, as in `rule`, applies at every interim analysis. A
list needs one condition per interim analysis, `K - 1` in all; a list of
`K` conditions is accepted with a warning and its last condition is
dropped. The rule does not store `m` or `K`: both are settled against
the p-values when the rule is used. Strings work too, with the words
`and`, `or` and `not` (`stopping_rule_gsd("r1 and (r2 or r3)")`).

Stopping only removes rejections, so the procedure keeps its error
control under any rule: a stopped trial keeps the rejections it has
declared and makes no new ones.

The gain below values each claim at its own decision look. It is
evaluated on a design with three hypotheses and $K = 2$ analyses:

``` r
set.seed(20261013)
sim3 <- simulate_pvalues_gsd(
    power_nominal = c(0.9, 0.8, 0.8),
    corr_matrix = matrix(0.5, 3, 3) + diag(0.5, 3),
    info_frac = c(0.5, 1),
    nsim = 20000L
)
pvals <- transform_pvalues_gsd(sim3, spending = gsDesign::sfLDOF)

gain <- trial_success_gsd(
    d(t1) + d(t2) + d(t3),
    d = c(1, 0.8)
)

# Optimise the graph under the rule
result <- graph_optimise_gsd(
    pvals, graph_constraint_free(3), gain,
    stopping_rule = rule
)

# Or evaluate a given graph under the rule
w <- c(0.5, 0.25, 0.25)
G <- rbind(c(0, 0.5, 0.5), c(0.5, 0, 0.5), c(0.5, 0.5, 0))
power <- calc_power_pvals_gsd(
    pvals, w, G,
    stopping_rule = rule,
    custom_power = gain
)

round(power$stop_distribution, 3)
#> analysis 1 analysis 2
#>      0.066      0.933
round(power$mean_stop_look, 3)
#> [1] 1.934
round(power$local_power, 3)
#> [1] 0.855 0.718 0.724
```

`stop_distribution` is the proportion of simulated trials that stopped
at each analysis, and `mean_stop_look` its mean. Every other field of
the result, and the gain, is computed from the stopped trials: here the
power of $H_{2}$ and $H_{3}$ is lower than the 0.732 and 0.741 the same
graph has without the rule. The optimised result records the rule, and
its [`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) show it with the
stopping-look distribution.

> **A trial that rejects every hypothesis stops in any case**
>
> Once every hypothesis is rejected there is nothing left to test, so
> the trial stops at that analysis whether or not a rule is supplied.
> Without a rule every other trial runs to the final analysis. In the
> example above, 3.3% of trials stop at the interim without the rule,
> all of them because every hypothesis was rejected there.

#### Optimising the spending parameter

A stopping rule also changes which spending function suits each
hypothesis. If the trial stops once $H_{1}$ is rejected, $H_{2}$ gets no
analysis after that stop, so spending more of $H_{2}$’s level at the
interim becomes worth more. In the user testing (sections 7.4 and 7.5 of
`dev/review/gsd_user_testing/README.md`) the design below had a best
gain of 0.6707 (at $w_{1} = 0.32$) with an O’Brien–Fleming-type function
for $H_{2}$. With a Pocock-type function for $H_{2}$ the best gain rose
to 0.6903 (at $w_{1} = 0.46$). That gain is 14 times the gap between the
best graph and Holm, so under the rule the choice of $H_{2}$’s spending
function mattered more than the graph did.

Instead of trying spending functions one at a time, you can tag a
hypothesis with
[`sf_hsd_opt()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/sf_hsd_opt.md).
Its Hwang–Shih–DeCani parameter $\gamma$ is then chosen together with
the graph. At an interim analysis at half the information, $\gamma = -4$
has spent about 12% of the level (close to O’Brien–Fleming) and
$\gamma = 1$ about 62% (close to Pocock). Given a tagged entry,
[`transform_pvalues_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/transform_pvalues_gsd.md)
returns a *spending spec* instead of transformed p-values. The spec
holds the raw p-values and a boundary table for every value the search
may try: by default $\gamma$ runs from $-4$ to $1$ in steps of $0.25$,
21 values, each checked for well ordering when the spec is built.
Building the tables is the slow part, about 20 seconds here. A spec can
be reused with other gains and stopping rules on the same draws.

``` r
set.seed(20261007)
raw <- simulate_pvalues_gsd(
    power_nominal = c(H1 = 0.90, H2 = 0.80),
    corr_matrix = matrix(c(1, 0.5, 0.5, 1), nrow = 2),
    info_frac = c(0.5, 1),
    nsim = 1e5
)

spec <- transform_pvalues_gsd(
    raw,
    spending = list(gsDesign::sfLDOF, sf_hsd_opt())
)
spec
#> <multigrain_pvals_gsd_spec>
#> 100000 simulation(s), 2 hypothes(es), 2 analys(es); alpha = 0.025
#> Information fractions (NA = no data at that analysis):
#>    analysis 1 analysis 2
#> H1        0.5          1
#> H2        0.5          1
#> Look-back: H1 = FALSE, H2 = FALSE
#> Spending:
#> H1: gsDesign::sfLDOF
#> H2: sfHSD(param = gamma in [-4, 1])
#> Spending parameters:
#> gamma: bounds [-4, 1], step 0.25, 21 nodes, start -1.5
#> Boundary tables built: 21
```

To optimise the graph and $\gamma$ together, pass the spec to
[`graph_optimise_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/graph_optimise_gsd.md)
in place of the p-values. The gain values $H_{1}$ at 0.6 and $H_{2}$ at
0.4, a claim at the final analysis is worth 0.8 of one at the interim,
and the trial stops once $H_{1}$ is rejected:

``` r
gain <- trial_success_gsd(0.6 * d(t1) + 0.4 * d(t2), d = c(1, 0.8))
rule <- stopping_rule_gsd(r1)

result <- graph_optimise_gsd(
    spec,
    graph_constraint_free(2),
    gain,
    stopping_rule = rule,
    num_threads = 2
)
result
#> <multigrain_graph_optimal>
#> Optimal graph found (given user-defined constraints on graph and computational resources):
#>
#> Hypothesis weights:
#>     H1     H2
#> 0.4447 0.5553
#>
#> Transition matrix:
#>    H1 H2
#> H1  0  1
#> H2  1  0
#>
#> Trial success function:
#> 0.6 * d(t1) + 0.4 * d(t2)
#>
#> Spending functions:
#>   H1: gsDesign::sfLDOF
#>   H2: sfHSD(param = -0.25)   [optimised: gamma, bounds -4 to 1, step 0.25]
#> Gain profile in gamma (graph re-optimised at each value):
#>   within 5e-4 of the maximum on [-0.75, 0]; fixed-graph plateau [-0.75, 0]
#>
#> Stopping rule:
#> Every interim: r1
#> Hypotheses referenced: up to r1
#>
#> Value of trial success measure:
#> 0.693
```

The search runs in two stages. First a genetic algorithm searches the
graph and $\gamma$ together. Then the graph is re-optimised at every one
of the 21 values of $\gamma$, starting each time from a neighbouring
value’s graph. The second stage traces how the best achievable gain
depends on $\gamma$, the *profile*. Most of the run (77 seconds here, on
two threads) is the genetic algorithm, as it is without
[`sf_hsd_opt()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/sf_hsd_opt.md).
The chosen value is in `result$spending_param`, the label per hypothesis
in `result$spending`, and the profile in `result$spending_profile`.
[`summary()`](https://rdrr.io/r/base/summary.html) prints the profile
(alongside the graph, the gain, the power and the constraints, which are
not shown here):

``` r
summary(result)
#> Spending-parameter profile (gamma):
#>  value   gain gain_fixed_graph   diff      se n_eval
#>  -4.00 0.6802           0.6802 0.0128 0.00035     25
#>  -3.75 0.6815           0.6815 0.0115 0.00033     27
#>  -3.50 0.6829           0.6828 0.0102 0.00034     26
#>  -3.25 0.6840           0.6839 0.0090 0.00035     24
#>  -3.00 0.6853           0.6852 0.0077 0.00034     28
#>  -2.75 0.6865           0.6864 0.0066 0.00030     30
#>  -2.50 0.6876           0.6876 0.0054 0.00030     28
#>  -2.25 0.6885           0.6885 0.0045 0.00027     24
#>  -2.00 0.6896           0.6896 0.0034 0.00026     25
#>  -1.75 0.6906           0.6906 0.0025 0.00025     25
#>  -1.50 0.6913           0.6913 0.0017 0.00023     25
#>  -1.25 0.6918           0.6917 0.0012 0.00021     28
#>  -1.00 0.6923           0.6923 0.0007 0.00019     27
#>  -0.75 0.6927           0.6927 0.0003 0.00016     27
#>  -0.50 0.6929           0.6929 0.0001 0.00012     26
#>  -0.25 0.6930           0.6930 0.0000 0.00000     25
#>   0.00 0.6928           0.6927 0.0002 0.00013     31
#>   0.25 0.6925           0.6925 0.0005 0.00016     25
#>   0.50 0.6921           0.6920 0.0009 0.00019     25
#>   0.75 0.6914           0.6913 0.0016 0.00022     26
#>   1.00 0.6906           0.6905 0.0024 0.00024     26
#> GA visited 21 of 21 nodes
```

The table has one row per value of $\gamma$. `gain` is the best gain
with the graph re-optimised at that value, and `gain_fixed_graph` is the
gain of the chosen graph kept as it is. `diff` is how far the row falls
short of the optimum, computed on the same simulated trials, and `se` is
the standard error of that difference.
[`print()`](https://rdrr.io/r/base/print.html) sums the table up as two
*plateaus*: the range of $\gamma$ around the optimum whose gain is
within $5 \times 10^{-4}$ of the best.

Here the optimum is $\gamma = -0.25$, between the two textbook choices,
with a gain of 0.6930. Both plateaus run from $-0.75$ to $0$, so any
$\gamma$ in that range is as good as the optimum to within the precision
of the simulation. Below $-1$ the gain falls steadily: at $\gamma = -4$,
close to O’Brien–Fleming, it is 0.0128 lower, about 37 standard errors.
Above $0.25$ the gain falls more slowly. The two gain columns barely
differ because the best graph hardly moves with $\gamma$ in this design.
Every row’s re-optimised weight on $H_{1}$ lies between 0.40 and 0.46,
and the plateau of $w_{1}$ under this rule is wide. When the best graph
does move with the parameter, the fixed-graph plateau is the narrower of
the two, and the re-optimised one is the one to read.

On the same draws, optimising the graph with $H_{2}$’s spending fixed
gives the user testing’s comparison (these runs are in the example
script):

| $H_{2}$ spending                                                                                                      | $w_{1}$ |   Gain |
|-----------------------------------------------------------------------------------------------------------------------|--------:|-------:|
| `sfLDOF` (O’Brien–Fleming type)                                                                                       |   0.292 | 0.6709 |
| `sfLDPocock` (Pocock type)                                                                                            |   0.549 | 0.6905 |
| [`sf_hsd_opt()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/sf_hsd_opt.md), chosen $\gamma = -0.25$ |   0.445 | 0.6930 |

These reproduce the 0.6707 → 0.6903 of the user testing. Choosing
$\gamma$ adds a further 0.0025 over a fixed Pocock function, which is
small next to the 0.0196 gained by moving away from O’Brien–Fleming. So
the main point of the earlier finding stands: under this rule, $H_{2}$
should spend early. The profile adds that $H_{2}$ need not spend quite
as early as Pocock, and that anywhere in $\lbrack-0.75,0\rbrack$ will
do.

The parameter and the graph were chosen on the same simulated trials
that score them, so the reported gain is optimistic, as it is for an
optimised graph alone. Check the chosen design on fresh draws by
building a spec from them with the same spending functions and resolving
it at the chosen value with
[`resolve_pvalues_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/resolve_pvalues_gsd.md):

``` r
set.seed(20261008)
raw_new <- simulate_pvalues_gsd(
    power_nominal = c(H1 = 0.90, H2 = 0.80),
    corr_matrix = matrix(c(1, 0.5, 0.5, 1), nrow = 2),
    info_frac = c(0.5, 1),
    nsim = 1e5
)
spec_new <- transform_pvalues_gsd(
    raw_new,
    spending = list(gsDesign::sfLDOF, sf_hsd_opt())
)
pv_new <- resolve_pvalues_gsd(spec_new, result$spending_param)

check <- calc_power_pvals_gsd(
    pv_new,
    hyp_weight = result$hyp_weight,
    trans_matrix = result$trans_matrix,
    stopping_rule = rule,
    custom_power = list(gain = gain)
)
round(check$gain, 4)
#> [1] 0.6936
```

On these fresh draws the design scores 0.6936, in line with the 0.6930
reported: here the optimism is smaller than the simulation noise. Only
the 21 values of the lattice can be resolved, since only their tables
were built and checked.
[`sf_hsd_opt()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/sf_hsd_opt.md)
takes `bounds` (within $-10$ to $2$) and `step` (`0.5`, `0.25` or
`0.125`) to change the lattice. If the optimum sits on a bound, the gain
may still be rising beyond it. Hypotheses tagged with the same `id`
share one parameter.

## Things that go wrong quietly

The gain-function constructors validate your expression, but some
mistakes compile without error and produce silently wrong results.

> **Do not include the zero in your discount table**
>
> The constructor adds a leading zero (for “never rejected”)
> automatically. If you include it yourself, `K` is inferred as one more
> than you intended and every claim is mispriced. The increase from `0`
> to the next value triggers a warning that names the table, but the
> table still compiles.
>
> ``` r
> # WRONG — K is inferred as 3; d(1) = 0, d(2) = 1
> trial_success_gsd(d(t1), d = c(0, 1, 0.75))
>
> # RIGHT — K = 2; d(1) = 1, d(2) = 0.75
> trial_success_gsd(d(t1), d = c(1, 0.75))
> ```

> **Do not use `!!` in discount table arguments**
>
> rlang’s `!!` operator is for the `objective` expression only. In
> ordinary arguments, `!!` is double negation: `!!0.75` evaluates to
> `TRUE`, which R coerces to `1`. The discount is silently lost.
>
> ``` r
> delta <- 0.75
>
> # WRONG — !!delta becomes TRUE, so d = c(1, 1)
> trial_success_gsd(d(t1), d = c(1, !!delta))
>
> # RIGHT — just use the variable directly
> trial_success_gsd(d(t1), d = c(1, delta))
> ```

> **Do not pass a GSD gain to the fixed-sample optimiser**
>
> [`graph_optimise()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/graph_optimisation.md)
> and
> [`calc_power_pvals()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/calc_power_pvals.md)
> supply a *logical* rejection matrix, not decision times, so a gain
> created with
> [`trial_success_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/trial_success_gsd.md)
> would score every rejection as though it happened at analysis 1. Both
> functions refuse such an object with an error that points to
> [`graph_optimise_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/graph_optimise_gsd.md)
> and
> [`calc_power_pvals_gsd()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/calc_power_pvals_gsd.md).

> **A stopping rule changes which graph is optimal**
>
> Under a stopping rule, stopping costs the hypotheses not yet rejected
> their later analyses, and the package does not value the stop itself:
> the gain counts only the claims made before the trial stopped. The
> optimiser responds by making the stop less likely, which can move the
> optimal graph away from the primary hypothesis. In the user testing
> (section 7.4 of `dev/review/gsd_user_testing/README.md`), with a rule
> that stops once $H_{1}$ is rejected, the optimal weight on $H_{1}$
> fell from 0.98 without the rule to 0.32 with it. Optimise the graph
> with and without the rule, compare the two, and decide what the design
> is for before adopting either.

> **Two rules of thumb**
>
> 1.  **Always parenthesise** expressions around `&&` and `||` to avoid
>     precedence surprises.
> 2.  **The gain is on your own scale**, not $\lbrack 0,1\rbrack$.
>     Normalise it yourself for reporting; only relative values affect
>     the optimiser.
