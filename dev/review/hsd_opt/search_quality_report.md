# Spending-parameter search: merge check and search quality

**Status: INTERIM, written 2026-11-04 at 23:10.** Parts 1 and 2 are complete. Part 3, the question of whether
the joint search finds good graphs, was still running when this was written and has no results yet.

Run with R 4.6.1, GA 3.2.5, nloptr 2.2.1, gsDesign 3.11.0, on `gsd-build-spending` at `e103474`.

## Summary

1. **The merge did not weaken the search.** Before and after `e103474`, the spending search is bit-for-bit
   identical, and the mutation has not changed.
2. **That does not show the search is good.** It shows only that the search did not get worse. Until now, every
   benchmark of it has used m <= 3, so its quality on larger designs has never been measured.
3. **Whether it finds good graphs on larger designs is being tested** in Part 3, including whether starting in
   the middle of the lattice is the cause.
4. **Merging `split-merge` would break the GSD search**, with no code conflict from git. See Part 2.

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

## Part 3: Does the joint search find good graphs, and does the start matter?

**In progress. No results yet.**

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

   So the population bunches at the lower bound and just above the middle, and is thin at the top of the lattice.
3. **The local sweep without a GA.** With `global_search = FALSE`, the sweep walks outward from the start node. With
   a GA, it walks outward from the GA's node instead.

### The test

Each design is run at fixed nodes with a full `graph_optimise_gsd()` as the benchmark, on every second lattice node
with 2 seeds. Against that:

| Arm | Change | Seeds |
|---|---|---|
| `joint_mid` | Default start, the lattice midpoint | 5 |
| `joint_low` | Started at the lower bound | 5 |
| `joint_high` | Started at the upper bound | 5 |
| `joint_unif` | Default start, spending genes uniform at initialisation | 5 |
| `local_mid` | `global_search = FALSE` from the default start | 1 |

Designs, all m = 4, nsim 1e4, one parameter shared by all hypotheses, default 21-node lattice, powers
(0.9, 0.85, 0.7, 0.6), block correlation 0.5 within pairs and 0.2 across:

- **H4:** `0.4 d(t1) + 0.3 d(t2) + 0.2 d(t3) + 0.1 d(t4)`, `d = c(1, 0.7)`.
- **H4x:** the same with `d = c(1, 0.4)`, to push the optimum towards aggressive spending.
- **H4g:** `(r1 || r2) * (0.6 + 0.2 d(t3) + 0.2 d(t4))`, `d = c(1, 0.7)`, a gated gain.

The results will separate the two ways the search can fail: choosing the wrong node, or choosing the right node
with a poor graph.

### Results

Pending.

## Other observations (none caused by the merge)

- **Seed noise can look like a regression.** On identical code and draws, final gain varied by 0.0019 on S3 over 10
  seeds, and by 0.0056 on the fixed-sample path at m = 6 over 3 seeds (0.7334 to 0.7391). Compare builds on the same
  seed and draws.
- **Verbosity rename.** `options(multigrain_verbosity = ...)` is now silently ignored, for example at
  `dev/review/hsd_opt/gate_p3_benchmarks.R:59`. It affects console output only.
- **Local builds.** The default library and `src/multigrain.dll` are 22 Oct builds with no GSD code, so
  `library(multigrain)` cannot currently run `sf_hsd_opt()`.

## Files

This report is not committed. The scripts, logs and raw results are in the session scratchpad, which is temporary:
`C:\Users\advsp\AppData\Local\Temp\claude\f6d42908bacbf4a5\f4b4e7c9-65a9-4292-9d26-e19527986f13\scratchpad\merge_regression\`
