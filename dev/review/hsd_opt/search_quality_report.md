# Spending-parameter search: merge check and search quality

**Status: final, 2026-11-05 02:00.** All runs finished with no failures.

Run with R 4.6.1, GA 3.2.5, nloptr 2.2.1, gsDesign 3.11.0, on `gsd-build-spending` at `e103474`.

## Summary

1. **The merge did not weaken the search.** Before and after `e103474`, the spending search is bit-for-bit
   identical, and the mutation has not changed.
2. **You are right that the search can fail to find good graphs, but the weakness is in the graph search, not the
   lattice.** It usually gets the spending parameter right or close, then lands on graphs of varying quality. With one
   parameter per hypothesis on H4xi, 29 of 45 runs fell below a gain inside their own search space, by up to 0.0038.
3. **The fixed-node GA has the same weakness, so this is not specific to the joint search.** At m = 4 the joint search is
   within the fixed-node GA's own seed noise at 5 of the 6 nodes where they can be compared. At m = 6 the same node gives
   gains up to 0.0069 apart across seeds with the fixed-node GA, and the joint search does at least as well.
4. **The GA matters less as m grows.** Without it, gains are up to 0.004 to 0.005 lower at m = 4, from any start. At
   m = 6 the one run without it was level with the GA runs.
5. **Starting in the middle is not established as a cause.** None of 14 comparisons of start or initial population reaches
   2 SE, and the direction changes between designs. The strongest hint in its favour, on H4xi, shrank from 1.9 to 1.2 SE
   with 15 seeds per arm. One lead remains: on that design the default start has more very poor runs (Part 4).
6. **Repeated runs give very different graphs with almost the same gain**, for the fixed-node GA as much as for the joint
   search. By eye that could look like a failure to find good graphs.
7. **Merging `split-merge` would break the GSD search**, with no code conflict from git. See Part 2.

**Your own evidence is still unexplained in detail.** I have not seen your design. Sending me the design and what you
compared it with is the most direct next step.

### Correction to what I said during the session

While the runs were in progress I said the H4 benchmark supported your observation, and that the joint search was
about 40 times less repeatable than the fixed-node GA at node -0.50. **Both claims were wrong.** They rested on two
fixed-node seeds that happened to agree (0.61230, 0.61233). With 5 seeds the fixed-node GA's standard deviation at that
node is 0.00050, against 0.00046 for the joint search.

## Part 1: Did the merge change the search?

No. Three builds were installed into private libraries from `git archive` exports: `5585caf` (immediately before
the merge), `e103474` (HEAD) and `e7a0f27` (14 Oct, before the lint commits and the earlier main merge).

| Test | Result |
|---|---|
| Mutation operator, before vs after | Same code; 10,000 replayed children per design identical |
| Functions changed by the merge | Added `graph_shortcut_kernel`, `graph_shortcut_table` and their native symbols; changed `multigrain_verbosity`; none of the 19 search functions checked |
| HSD search, before vs after the merge | 10 of 10 runs identical, including the spec, final RNG state and GA population |
| HSD search, 14 Oct build vs now | 10 of 10 runs identical |
| Same build run twice (negative control) | Identical |
| `p_param_mutate = 0.2` (positive control) | Differs from GA generation 2, so the test can detect a mutation change |
| Recorded baselines | Reproduced: S3 seed 7 gain 2.3986 at -4; D2 seed 1 gain 0.694014 at (-2, 0) |
| Fixed-sample `graph_optimise()`, before vs after, m = 3, 4, 6 | 9 of 9 runs identical, although the merge changed its kernel |
| Table kernel vs sequential kernel | 9 of 9 runs identical |

Comparisons used `identical(num.eq = FALSE)`, so doubles were compared bit for bit.

Designs: S3 is `dev/review/hsd_opt/acceptance_m3.R` (m = 3, nsim 1e4, shared parameter, stopping rule). M1 is the
m = 2 design of `gate_p3_benchmarks.R` at nsim 2e4. D2 is its two-parameter version with a coarse lattice. Seeds
were 1, 7 and 101, with 2 threads, plus 1 thread for S3 seed 7.

**What this does not cover:** a result seen on another build, or a comparison between runs with different seeds.

## Part 2: How has the mutation changed?

It has not changed on this branch. `R/mutation_helpers.R` is the same blob (`8fb295a`) at the merge base, before
the merge, on main and at HEAD. The search uses `.make_cauchy_mutation_multi(0.1, 1.0)`.

What the operator does to a spending gene:

| Design | Genes (graph + spending) | Chance a call changes a spending gene |
|---|---|---|
| S3 | 5 + 1 | 18% |
| M1 | 1 + 1 | 50% |
| D2 | 1 + 2 | 34 to 35% each |

A changed gene jumps about 0.30 of its range on average, about 6 nodes on a 21-node lattice, so this is a global
move and not a local one.

### Where it does change: `split-merge`

`split-merge` (tip `ff9d907`) is not in main. It replaces the operator and has no GSD code. A `git merge-tree` dry
run into this branch shows:

- Conflicts only in `.gitignore` and `NEWS.md`.
- No definition of `.make_cauchy_mutation_multi`, while `R/optimisation_gsd.R:572` and
  `R/optimisation_gsd_opt.R:380` still call it. Read from the merged tree, not run. It should fail with "could not
  find function" whenever `global_search = TRUE` on the GSD path.
- `pcrossover = 0`, `pmutation = 1` in `default_control()`, which reach the GSD paths through
  `control_prepare_dims()`.

`gh` was not logged in, so I could not check whether `split-merge` has an open PR.

### Would that operator hurt the spending search?

`.sm_groups()` covers graph genes only, so the split-merge move never changes a spending gene. I expected that to
do harm. On these designs it did not. Three designs, 10 seeds, paired on the same draws, set through
`control_global()` with no source edited:

| Design | Arm | Nodes visited by GA | Final picked the best node | Mean final gain | Paired difference from `cauchy` (SE) |
|---|---|---|---|---|---|
| S3 | `cauchy` | 21.0 of 21 | 10 of 10 | 2.398930 | |
| S3 | `sm_frozen` | 19.4 of 21 | 10 of 10 | 2.398750 | -0.000180 (0.000342) |
| S3 | `sm_redraw` | 21.0 of 21 | 10 of 10 | 2.398630 | -0.000300 (0.000230) |
| M1 | `cauchy` | 21.0 of 21 | 10 of 10 | 0.694978 | |
| M1 | `sm_frozen` | 20.1 of 21 | 10 of 10 | 0.694978 | 0 (all tied) |
| M1 | `sm_redraw` | 21.0 of 21 | 10 of 10 | 0.694978 | 0 (all tied) |
| D2 | `cauchy` | 10.0 of 10 | 10 of 10 | 0.694014 | |
| D2 | `sm_frozen` | 10.0 of 10 | 10 of 10 | 0.694014 | 0.000000 (0.000000) |
| D2 | `sm_redraw` | 10.0 of 10 | 10 of 10 | 0.694011 | -0.000003 (0.000003) |

`sm_frozen` is `split-merge` as written with its defaults. `sm_redraw` adds a 20% uniform redraw of each spending
gene. The S3 95% intervals are (-0.00095, +0.00059) and (-0.00082, +0.00022).

**This does not clear `split-merge` for larger designs.** These designs are small, M1 and D2 tie almost everywhere,
and the initial population favours S3's optimum (see Part 3).

## Part 3: Does the joint search find good graphs at m = 4, and does the start matter?

### What "starts in the middle" means in the code

It happens in three separate places:

1. **Seeded graphs.** `sf_hsd_opt(start = NULL)` uses the node nearest the midpoint of `bounds`, -1.5 by default.
   Every seeded graph in the GA gets that node (`R/optimisation_gsd_opt.R:363-372`).
2. **The initial population.** `esch_population_rcpp()` is far from uniform on a spending gene, although the comment
   at `R/objective_function_gsd_opt.R:6-8` gives each node equal mass only "under a uniform gene". Measured on
   200,000 draws, 21-node lattice, where uniform would be 4.8% per node:

   | Node value | -4 | -3.75 | -3.5 | -2 | -1.75 | -1.5 | -1.25 | -1 | 0 | 0.5 | 1 |
   |---|---|---|---|---|---|---|---|---|---|---|---|
   | Share of genes | 16.2% | 11.4% | 7.2% | 1.0% | 0.8% | 8.9% | 14.1% | 9.1% | 1.9% | 1.1% | 0.7% |

3. **The local sweep without a GA.** With `global_search = FALSE`, the sweep walks outward from the start node. With
   a GA, it walks outward from the GA's node instead.

### The test

All three designs have m = 4, nsim 1e4, one parameter shared by all hypotheses, the default 21-node lattice, powers
(0.9, 0.85, 0.7, 0.6), block correlation 0.5 within pairs and 0.2 across, and no stopping rule. All use the same draws.

- **H4:** `0.4 d(t1) + 0.3 d(t2) + 0.2 d(t3) + 0.1 d(t4)`, `d = c(1, 0.7)`.
- **H4x:** the same with `d = c(1, 0.4)`, which puts the optimum at the thin top of the lattice.
- **H4g:** `(r1 || r2) * (0.6 + 0.2 d(t3) + 0.2 d(t4))`, `d = c(1, 0.7)`, a gated gain.

The benchmark is a full `graph_optimise_gsd()` at fixed nodes: every second node with 2 seeds, then 5 seeds at the
nodes the joint search chose most. The joint search ran 5 seeds in each of four arms: default start, low start, high
start, and default start with uniform initial spending genes. `global_search = FALSE` ran once from each start.

### Does it choose the right node?

Mostly, and where it misses, little is lost.

| Design | Best fixed-node mean | Nodes the joint search chose (20 runs) |
|---|---|---|
| H4 | -0.25 (0.612374) | -0.25 x8, -0.50 x10, -0.75 x2 |
| H4x | 1.00 (0.508430) | 1.00 x20 |
| H4g | -3.75 (0.708495) | -4 x4, -3.75 x3, -3.5 x7, -3.25 x2, -3 x2, -2.75 x1, -2.5 x1 |

The H4g fixed-node means are nearly flat in this region: 0.708439, 0.708495, 0.708248 and 0.707912 at -4, -3.75, -3.5
and -3.25. H4x is the strongest evidence against a start problem: its optimum is at the node with the fewest initial
genes (0.7%), and all 20 runs found it, including those started at the other end of the lattice.

### Does it find good graphs at that node?

Fixed-node GA against the joint runs that chose the same node, where each side has at least 3 runs:

| Design | Node | Fixed-node GA, 5 seeds: mean (SD) | Joint search: mean (SD), runs | Fixed minus joint (SE) | Rank-sum p |
|---|---|---|---|---|---|
| H4 | -0.50 | 0.612033 (0.000495) | 0.611722 (0.000458), 10 | +0.000311 (0.000265) | 0.165 |
| H4 | -0.25 | 0.612374 (0.000106) | 0.612445 (0.000205), 8 | -0.000071 (0.000087) | 0.354 |
| H4x | 1.00 | 0.508430 (0.000146) | 0.508283 (0.000509), 20 | +0.000148 (0.000131) | 0.974 |
| H4g | -4.00 | 0.708439 (0.000297) | 0.708388 (0.000283), 4 | +0.000051 (0.000194) | 0.905 |
| H4g | -3.75 | 0.708495 (0.000604) | 0.708687 (0.000311), 3 | -0.000192 (0.000324) | 0.786 |
| H4g | -3.50 | 0.708248 (0.000257) | 0.707656 (0.000213), 7 | +0.000592 (0.000140) | 0.003 |

Only H4g at -3.50 clearly favours the fixed-node GA. Treat that one with care: these rows condition on the node the
joint search chose, and a run tends to stop at a node when it failed to find something better elsewhere.

Without that conditioning, comparing all 20 joint runs with the fixed-node GA at its best node:

| Design | Fixed-node GA at best node, mean | Joint, default start, mean | Joint, all 20 runs, mean | Average shortfall, all runs | Highest gain, fixed vs joint |
|---|---|---|---|---|---|
| H4 | 0.612374 | 0.611808 | 0.612027 | 0.00035 | 0.612506 vs 0.612654 |
| H4x | 0.508430 | 0.508518 | 0.508283 | 0.00015 | 0.508568 vs 0.508882 |
| H4g | 0.708495 | 0.707605 | 0.707907 | 0.00059 | 0.708958 vs 0.709044 |

The best node was picked after seeing the fixed-node means, which flatters the benchmark. The single highest gain came
from a joint run on all three designs.

**Cost:** the fixed-node GA took about 80 s per node and the joint search 104 to 144 s per run. A full outer loop over
21 nodes would take about 12 to 17 times as long. These times are rough, because processes shared the machine.

### Does the start matter?

Not consistently. Mean final gain against the default start, 5 seeds each:

| Design | Default start, mean | Low start minus default (SE) | High start minus default (SE) | Uniform spending genes minus default (SE) |
|---|---|---|---|---|
| H4 | 0.611808 | +0.000476 (0.000246) | +0.000349 (0.000342) | +0.000052 (0.000320) |
| H4x | 0.508518 | -0.000083 (0.000235) | -0.000303 (0.000293) | -0.000554 (0.000385) |
| H4g | 0.707605 | +0.000282 (0.000329) | +0.000532 (0.000327) | +0.000395 (0.000384) |

None of the nine differences reaches 2 SE, and the sign changes between designs. A small effect is not excluded: five
seeds per arm cannot rule out differences of about 0.0005 to 0.001.

### `global_search = FALSE` is where graphs are poor

| Design | Start at -4 | Start at -1.5 | Start at 1 | Joint search with GA, 20-run mean |
|---|---|---|---|---|
| H4 | 0.609020 at -0.25 | 0.608140 at -0.50 | 0.608730 at -0.50 | 0.612027 |
| H4x | 0.504040 at 1.00 | 0.504570 at 1.00 | 0.503170 at 1.00 | 0.508283 |
| H4g | 0.707500 at -4.00 | 0.707404 at -4.00 | 0.707460 at -4.00 | 0.707907 |

From any start, the sweep reaches a good node, but the graph is poorer by up to about 0.004 to 0.005 on H4 and H4x. The
start changes the gain by at most 0.0014, with no pattern. So the GA, not the start, is what supplies the good graph.
Each cell is one run.

### Other findings

- **Different graphs, nearly the same gain.** The mean absolute difference per graph entry between two runs is 0.15 to
  0.16 on H4 and H4x and about 0.10 on H4g. That is the same within fixed-node runs, within joint runs and between the
  two: 0.1527, 0.1544, 0.1698 on H4 at -0.50. The landscape has many near-equivalent graphs, so comparing graphs by eye
  will exaggerate any difference.
- **The GA decides the result.** The local stage moved the node in 1, 0 and 3 of 20 runs on H4, H4x and H4g, and otherwise
  added a median of about 0.00003. When it did move the node, on H4g, it added a median of 0.001.
- **The profile runs low.** Its gains at other nodes are typically 0.0005 to 0.0018 below the best of two fixed-node runs.
  Taking the best of two flatters the benchmark, so this is an upper estimate. Read the plateaus as the lower bound the
  design says they are.

## Part 4: Larger designs

### H6: m = 6, one shared parameter

Powers (0.9, 0.85, 0.8, 0.7, 0.6, 0.5), block correlation as before, nsim 1e4, gain
`0.3 d(t1) + 0.25 d(t2) + 0.2 d(t3) + 0.1 d(t4) + 0.1 d(t5) + 0.05 d(t6)`, `d = c(1, 0.7)`.

| Search | Runs | Mean final gain | SD | Range |
|---|---|---|---|---|
| Fixed-node GA at -2 | 2 | 0.554753 | | 0.553040 to 0.556465 |
| Fixed-node GA at -1 | 5 | 0.554507 | 0.00280 | 0.551965 to 0.557928 |
| Fixed-node GA at -0.5 | 5 | 0.555566 | 0.00276 | 0.551303 to 0.558234 |
| Fixed-node GA at 0 | 5 | 0.556448 | 0.00089 | 0.555548 to 0.557581 |
| Fixed-node GA at 1 | 2 | 0.550857 | | 0.549417 to 0.552296 |
| Joint search, default start | 5 | 0.557058 | 0.00095 | 0.55588 to 0.55823 |
| Joint search, low start | 5 | 0.557010 | 0.00194 | 0.55437 to 0.55899 |
| No GA, default start | 1 | 0.557426 | | |

- **The joint search is not the weak link.** Its mean is above the mean of every fixed node, and it chose -0.50 six times
  and -0.75 four times.
- **The graph search as a whole is unreliable at this size.** Seed noise at a single node (SD up to 0.0028, range up to
  0.0069) is larger than the differences between the node means tested at 5 seeds (at most 0.0019), so a single GA run
  cannot tell those nodes apart.
- **The GA adds little here.** The single run without it is level with the GA runs, and the local stage added a median of
  0.0012 to the GA's result, against about 0.00003 at m = 4. So the GA is further from converged when it stops.
- **The true optimum is unknown.** The highest gain seen, 0.558995, is only a lower bound on it, so the real shortfall could
  be larger than these numbers suggest.

### H4i and H4xi: one spending parameter per hypothesis

These use the H4 and H4x gains, with four parameters and a search space of 21^4 lattice points. Sharing one value across
hypotheses is inside that space, so a sound search should at least match the shared-parameter fixed-node GA at its best
node, on the same draws. That is the reference below, with 5 seeds.

| Design | Reference, mean (SD) | Arm | Seeds | Mean final gain (SD) | Range | Minus reference | Runs below reference |
|---|---|---|---|---|---|---|---|
| H4i | 0.612374 (0.000106) | Default start | 15 | 0.615256 (0.000411) | 0.614549 to 0.615725 | +0.002882 | 0 |
| H4i | | High start | 5 | 0.614956 (0.000347) | 0.614563 to 0.615385 | +0.002582 | 0 |
| H4i | | Uniform initial genes | 15 | 0.614851 (0.000784) | 0.613107 to 0.615786 | +0.002477 | 0 |
| H4xi | 0.508430 (0.000146) | Default start | 15 | 0.507767 (0.001215) | 0.504630 to 0.508962 | -0.000663 | 10 |
| H4xi | | High start | 15 | 0.508206 (0.000360) | 0.507320 to 0.508646 | -0.000225 | 10 |
| H4xi | | Uniform initial genes | 15 | 0.508159 (0.000485) | 0.507346 to 0.508802 | -0.000271 | 9 |

- **H4i works.** All 35 runs beat the shared-parameter gain by using the extra freedom. The default-start runs agree
  roughly on the pattern: H1 0.25 to 0.75, H2 -0.25 to -1.25, H3 -1.75 to -2.75, H4 -3 to -4. The reported gain is
  optimistic, as the parameters and graph are chosen on the same draws.
- **H4xi is the failure.** The extra freedom should not lose anything, yet 29 of 45 runs fall below the reference, by up to
  0.0038. Since the best run reached 0.508962, every one of those runs missed a better design that was inside its search
  space.
- **That failure is in the graph, not the lattice.** All 45 H4xi runs put H1 and H2 at the upper bound or close to it.
  Twenty-one runs chose exactly (1, 1, 1, -0.25), and their gains range from 0.506134 to 0.508836. Same spending
  parameters, a range of 0.0027.
- **Convergence:** `profile_converged` was `FALSE` in 11 of 80 runs.

### Is the start involved? 15 seeds on the strongest hint

Across all designs, 14 comparisons of start or initial population were made against the default start, and none reaches
2 SE. H4xi is the design where your hypothesis should bite hardest, because its optimum is in a corner of the lattice where
the initial population is thinnest. At 5 seeds the uniform arm looked better by +0.000949 (1.9 SE), so I re-ran the
comparison with 10 more seeds per arm:

| Design | Comparison with the default start | Difference (SE) | 95% CI | Welch p | Rank-sum p |
|---|---|---|---|---|---|
| H4xi | High start | +0.000439 (0.000327) | -0.00025 to +0.00113 | 0.198 | 0.838 |
| H4xi | Uniform initial genes | +0.000393 (0.000338) | -0.00032 to +0.00110 | 0.260 | 0.806 |
| H4i | Uniform initial genes | -0.000405 (0.000229) | -0.00088 to +0.00007 | 0.091 | 0.161 |

- **The hint shrank.** It fell to 1.2 SE on H4xi, the medians barely differ, and one very poor default-start run (seed 14,
  0.50463) accounts for much of the gap: without it, the uniform arm's advantage is about 0.00017.
- **The same change goes the other way on H4i**, where the uniform arm is lower.
- **Any effect on the mean is small.** Even at the top of its interval, about 0.0011, it is well below the 0.0027 range seen
  at one fixed set of spending parameters.
- **One lead remains: the spread.** On H4xi the default start's SD is 0.00122, against 0.00036 and 0.00049, and it has all
  4 of the runs below 0.507 (4 of 15, against 0 of 30). On H4i it reverses: 0.00041 for the default against 0.00078 for
  uniform. So the initial population may affect how often a run goes badly, in a way that depends on the design. I noticed
  this after seeing the data and chose the 0.507 cut-off afterwards, so it is a lead and not a finding.

### What I would do with this today

Run several seeds and keep the best. It needs no code change, and blocks of 5 default-start runs were repeatable:

| Design | Best of each block of 5 | Single-run SD |
|---|---|---|
| H4i | 0.615646, 0.615725, 0.615685 | 0.00041 |
| H4xi | 0.508602, 0.508836, 0.508962, all above the reference | 0.00122 |

On H6, with one block per start, the bests were 0.55823 and 0.55899, against a single-run SD of 0.00095 to 0.00194, so
expect less repeatability there. Keeping the best costs 5 times the computing and increases the optimism in the reported
gain, so check the chosen design on fresh draws.

## What has not been tested

Your design, constrained graphs, more than two analyses, stopping rules at m >= 4, nsim above 2e4, and any change to the
GA's budget, such as a larger `popSize` or `run`.

## Other observations (none caused by the merge)

- **Seed noise can look like a regression.** On identical code and draws, final gain varied by 0.0019 on S3 over 10
  seeds, and by 0.0056 on the fixed-sample path at m = 6 over 3 seeds (0.7334 to 0.7391). Compare builds on the same
  seed and draws, and use more than 2 seeds when judging repeatability.
- **Verbosity rename.** `options(multigrain_verbosity = ...)` is now silently ignored, for example at
  `dev/review/hsd_opt/gate_p3_benchmarks.R:59`. It affects console output only.
- **Local builds.** The default library and `src/multigrain.dll` are 22 Oct builds with no GSD code, so
  `library(multigrain)` cannot currently run `sf_hsd_opt()`.

## Files

This report is the only file committed. The scripts, logs and raw results are in the session scratchpad, which is
temporary:
`C:\Users\advsp\AppData\Local\Temp\claude\f6d42908bacbf4a5\f4b4e7c9-65a9-4292-9d26-e19527986f13\scratchpad\merge_regression\`
