# User testing of the group sequential extension

Tester: Niamh Fitzgerald. Date: 2026-10-05. Branch tested: `gsd-build` at `6e5e509`; the two commits of this pull request sit on top of it.

This folder is a first outside run of the `_gsd` family: a 2-analysis trial with 2 hypotheses, the same with 3 hypotheses, an independent check of both, the features those two do not touch (look-back, an endpoint that finishes early, an endpoint with no interim data, mixed spending functions, a third analysis, a fourth hypothesis, constraints), the two case studies of the graphicalMCP vignette, the package's own gates, and a prototype for the stopping-criteria gap. Nothing under `R/` or `src/` is changed. Unless a source is named, the numbers below come from the scripts in this folder: the `*_output.txt` files are the logs of the runs reported here and the `.rds` files hold the same results at full precision.

Claims are labelled **[proved]** when they follow from the structure of the algorithm, and **[simulated]** when they rest on Monte Carlo runs on one design.

## Summary

1. **It works as documented.** Both optimisations ran first time from the roxygen examples and the design record. They took 68 to 107 s (2 hypotheses) and 117 to 154 s (3 hypotheses) at `nsim = 1e5` with default control and 2 threads, against the roughly 5 minutes I was told to expect.
2. **The answers are right, as far as I can check them.** An independent reference (Maurer and Bretz Algorithm 1 on raw p-values, boundaries by quadrature, no gsDesign, no transform, no C++) agrees with transform plus kernel on every one of 500,000 decision times on the optimisation draws, and on all but 1 of 1,000,000 on fresh draws. The optimiser matches a grid search (2 hypotheses) and beats a 20,000-graph random search (3 hypotheses). FWER of the optimised graphs is at or below 0.025 within Monte Carlo error in all 10 null configurations.
3. **The package's own gates hold on a second machine.** All 444 expectations in the 9 GSD test files pass on R 4.3.3 (the record was developed on R 4.6.1) with no skips; the full suite gives 1622 expectations with no failures; and `dev/gsd_identity_check.R` finds the fixed-sample `graph_optimise()` identical on `main` and on the branch (section 4.6).
4. **The rest of the feature set behaves.** Both case studies of the graphicalMCP vignette are reproduced decision for decision. Example 5 of the manuscript is reproduced to within Monte Carlo error in all 21 cells once the paper's boundary is used. Look-back, an endpoint complete at the interim, an endpoint with no interim data, mixed spending functions, three analyses, a constrained graph and four hypotheses all ran and checked out (section 4).
5. **Documentation is behind the code in four places**, one of them a wrong number in the vignette (section 5). The second commit of this pull request fixes the vignette; drop it if P5 already covers it.
6. **Discounting interacts with the spending function.** Under O'Brien-Fleming-type spending the interim boundary is strongly convex in the level a hypothesis holds, so a discount on late claims pushes the optimum towards concentrating alpha, and in a symmetric design it makes the gain two-peaked with a trough around the Holm point. Under Pocock-type spending neither happens (section 6).
7. **Stopping criteria.** An efficacy stopping rule that depends only on which hypotheses have been rejected is a deterministic edit of the decision-time matrix, so it needs no kernel or transform change, and small cases can already be written in today's gain grammar. On a 2-hypothesis example, ignoring a "stop when the primary is rejected" rule overstates the value of the design by about 7%, moves the optimal weight on the primary from 0.98 to 0.32, and makes the choice of spending function for the secondary worth more than the choice of graph (section 7). A non-binding futility rule can be applied as a mask before the search. Suggestions for the package are collected in section 8.

## Files

| File | What it does |
|---|---|
| `00_reference.R` | Independent reference for K = 2: simulator by independent increments, boundaries by one-dimensional quadrature, Algorithm 1 in plain R with per-hypothesis information fractions and optional stopping rules |
| `01_two_hyp_two_stage.R` | 2 hypotheses, 2 analyses: simulate, transform, gain with a discount, optimise, compare with simple graphs and a grid |
| `02_three_hyp_two_stage.R` | The same with 3 hypotheses and a hurdle on the primary; two seeds; random-search check |
| `03_verify.R` | Checks A to G on both optimised designs |
| `04_stopping_rules.R` | Stopping-rule prototype |
| `05_discount_sensitivity.R` | Discount and spending function against the optimal split of alpha |
| `06_graphicalmcp_vignette.R` | The two case studies of the graphicalMCP group sequential vignette through multigrain |
| `07_features.R` | Manuscript Example 5, look-back, and a mixed design with an endpoint that has no interim data |
| `08_scaling_constraints.R` | Three analyses, a constrained graph, four hypotheses, reproducibility |

## Environment

R 4.3.3 on Linux, 2 cores. multigrain 0.3.0 installed from the branch (`R CMD INSTALL`, not `load_all()`), gsDesign 3.11.0, graphicalMCP 0.3.0, gMCPLite 0.1.7, rlang 1.3.0, RcppParallel 6.2.1, ggplot2 3.4.4. gsDesign was built without its `gt` and `r2rtf` imports, which only affects its table-formatting functions; `gsBound1()` and the spending functions are untouched.

## 1. Two hypotheses, two analyses

Design: H1 nominal power 0.80, H2 0.90, correlation 0.5, analyses at information fractions 0.5 and 1, `sfLDOF` for both, one-sided alpha 0.025, `nsim = 1e5`. H1 is the more valuable claim and the less well powered one.

Gain: `0.6 * d(t1) + 0.4 * d(t2)` with `d = c(1, 0.8)`, so a claim made at the final analysis keeps 80% of its value. The undiscounted gain `0.6 * r1 + 0.4 * r2` was optimised as well.

| Graph | w1 | Discounted gain | Undiscounted gain | Power H1 | Power H2 | Interim H1 | Interim H2 |
|---|---:|---:|---:|---:|---:|---:|---:|
| Optimised for the discounted gain | 0.149 | 0.6814 | 0.8160 | 0.767 | 0.890 | 0.089 | 0.224 |
| Optimised for the undiscounted gain | 0.241 | 0.6807 | 0.8162 | 0.771 | 0.885 | 0.093 | 0.208 |
| Holm | 0.500 | 0.6783 | 0.8154 | 0.779 | 0.869 | 0.109 | 0.161 |
| Fixed sequence H1 then H2 | 1 | 0.6501 | 0.7794 | 0.799 | 0.750 | 0.163 | 0.087 |
| Fixed sequence H2 then H1 | 0 | 0.6782 | 0.8094 | 0.750 | 0.898 | 0.087 | 0.252 |

Both optimised graphs recycle fully (transition weights of 1). With two hypotheses the only free parameter is w1, so the optimiser can be checked against a grid of step 0.01:

- discounted: grid maximum 0.68135 at w1 = 0.15, optimiser 0.68135 at w1 = 0.149; every w1 in [0.05, 0.21] is within 5e-4 of the maximum;
- undiscounted: grid maximum 0.81620 at w1 = 0.24, optimiser 0.81624 at w1 = 0.241; plateau [0.15, 0.44].

The discount moves the optimum from 0.24 to 0.15, towards the hypothesis more likely to be rejected at the interim. The gain is flat near its maximum, as the record found for Figure 3b, so the weight is resolved only to the width of the plateau.

Time: simulate 0.1 s, transform 1.6 s, optimise 107 s (discounted) and 68 s (undiscounted).

## 2. Three hypotheses, two analyses

Design: nominal powers 0.90, 0.80, 0.70, pairwise correlation 0.5, information fractions 0.5 and 1, `sfLDOF`, `nsim = 1e5`.

Gain: `0.5 * d(t1) + r1 * (0.3 * d(t2) + 0.2 * d(t3))` with `d = c(1, 0.8)`. The secondaries have value only if the primary is rejected.

| Graph | Discounted gain | Power H1 | Power H2 | Power H3 | Interim H1 | Interim H2 | Interim H3 |
|---|---:|---:|---:|---:|---:|---:|---:|
| Optimised, seed 1 | 0.6670 | 0.898 | 0.745 | 0.608 | 0.251 | 0.084 | 0.036 |
| Optimised, seed 2 | 0.6670 | 0.898 | 0.744 | 0.608 | 0.251 | 0.084 | 0.036 |
| Optimised for the undiscounted gain | 0.6669 | 0.898 | 0.743 | 0.610 | 0.251 | 0.083 | 0.036 |
| Fixed sequence H1, H2, H3 | 0.6660 | 0.898 | 0.751 | 0.592 | 0.251 | 0.088 | 0.036 |
| H1 gate, then Holm on H2 and H3 | 0.6624 | 0.898 | 0.712 | 0.635 | 0.251 | 0.059 | 0.045 |
| Holm | 0.6179 | 0.837 | 0.736 | 0.649 | 0.115 | 0.073 | 0.051 |

The optimised graph puts all initial weight on H1 and sends 93% of its level to H2 and 7% to H3:

```
w = (1, 0, 0)        H1    H2    H3
                 H1   0  0.93  0.07
                 H2   1  0     0
                 H3   1  0     0
```

The edges from H2 and H3 back to H1 look odd, since H1 is always rejected first, but the graph update turns them into full recycling between H2 and H3 once H1 is gone: (0 + 1 x 0.07) / (1 - 1 x 0.93) = 1. A reader of the printed graph could miss that.

The two seeds agree on the gain to five decimals (0.66696) and on the split to 0.004 (0.930 against 0.926). A random search over 20,000 graphs on the same draws reached 0.66521 at best, below the optimiser. The discount hardly changes this graph (0.93 against 0.91 for the undiscounted optimum), because the hurdle already forces everything through H1.

Time: simulate 0.1 s, transform 2.6 s, optimise 143 s, 117 s and 154 s; the random search took 84 s.

## 3. Verification (`03_verify.R`)

| Check | 2 hypotheses | 3 hypotheses |
|---|---|---|
| A. Reference on the raw p-values against transform plus kernel, all 1e5 trials | 0 of 200,000 decision times differ | 0 of 300,000 |
| B. Compiled gain against the same gain in plain R | difference 6e-13 | difference 6e-13 |
| C. Fresh draws from the independent simulator, scored by the reference (2e5 trials) | 0.68009 (s.e. 0.00064) against 0.68135 in sample, z = 1.1 | 0.66595 (s.e. 0.00061) against 0.66696, z = 1.0 |
| C. The same fresh draws through transform plus kernel | 1 of 400,000 decision times differs | 0 of 600,000 |
| D. `graphicalMCP::graph_test_shortcut_gsd()` on 300 trials | 0 disagreements in 600 reject or retain decisions and in the decision times of 494 rejections | 0 in 900 and 662 |
| E. FWER of the optimised graph, every set of true nulls, 1e6 trials each (s.e. 0.00016) | 0.0238 to 0.0247 | 0.0242 to 0.0251 |
| F. Power fields consistent with the kernel output and with each other | all 5 ok | all 5 ok |
| G. Fresh package draws (1e6): optimised minus best simple graph, paired | +0.0030 against Holm and against H2 then H1 (z about 30) | +0.0009 against the fixed sequence (z = 26) |

The single disagreement in C is the kind the record predicts: a repeated p-value within the interpolation error of its threshold, of order one decision in a million.

G shows the optimiser's advantage survives out of sample, so it is not an artefact of fitting to the optimisation draws. It is small in both designs: 0.4% of the gain over the best simple graph with 2 hypotheses and 0.1% with 3. Against Holm with 3 hypotheses it is 7.9%, which is the hurdle and not the search.

All of section 3 is **[simulated]** on two designs with K = 2 and a common information fraction. Section 4 extends the checks to an endpoint that is complete at the interim, an endpoint with no interim data, unequal information fractions, look-back and a third analysis.

## 4. Other features, a third analysis and a fourth hypothesis

### 4.1 The graphicalMCP vignette (`06_graphicalmcp_vignette.R`)

The two case studies of the graphicalMCP 0.3.0 vignette "Group Sequential Design with Graphical Approaches", each run as one observed trial through `transform_pvalues_gsd()` and the kernel and compared with `graph_test_shortcut_gsd()`.

| Case | Setting | Rejections | Decision times |
|---|---|---|---|
| Diabetes trial of Maurer and Bretz: 4 hypotheses, analyses at 1/3 and 2/3 | look-back off | agree: H1, H2, H3 at analysis 2 | agree |
| Oncology trial: 6 hypotheses with 3, 3, 2, 2, 1 and 1 analyses | look-back on, the vignette's p-values | agree: H5 at analysis 1, H1 and H3 at analysis 2 | agree |
| The same | look-back on, the "look-back makes a difference" p-values | agree: H1 to H5 | agree with `decision_at` |
| The same | those p-values, look-back off | agree: H1, H2, H3 | agree |

Repeated p-values for the diabetes trial match Table 2 of the paper to 8.5e-05, the ADDPLAN against gsDesign residual noted in the record.

Two points for a reader coming from that vignette. First, in the third row graphicalMCP attributes H2, H4 and H5 to analysis 1 (`first_rejected_at`, the earliest boundary crossed), whereas multigrain's `t<i>` is 2, the analysis at which the rejection could be declared; for a gain the second is the right one. Second, graphicalMCP wants `NA` after an endpoint's last analysis and multigrain wants its information fraction carried forward at 1; the script converts one to the other.

### 4.2 Example 5 of the manuscript (`07_features.R`, part A)

PFS complete at the interim, OS at information fractions 0.7 and 1, correlation 0.5, inputs and reference values from `tests/testthat/data/gsd_example5_reference.rds`.

- All 21 (r, delta) cells at `nsim = 1e5` with `sfLDOF`: using the paper's optimal weight in place of this run's grid optimum loses at most 2.4e-04 of gain, the flatness the record reports.
- The expected gain at the paper's weight is below the paper's by 0.0002 to 0.0078, more as delta falls. That is the boundary: the paper used the classical O'Brien-Fleming boundary (0.0082 at the interim) where `sfLDOF` gives 0.0074.
- With the paper's boundary supplied as a spending function the shortfall is 0.0003 to 0.0016 at `nsim = 1e5`. At `nsim = 1e6` it is gone: this run minus the paper, over the 21 cells, runs from +0.0001 to +0.0003 on one seed and from -0.0001 to +0.0004 on another, against a standard error of roughly 0.0001 to 0.0003 for a gain at that size. The package reproduces the paper's expected gains to within Monte Carlo error. **[simulated]**
- The optimiser at three cells returns w1 = 0.875, 0.627 and 0.442 against grid optima 0.875, 0.665 and 0.455, each at a gain within 4e-05 of the grid maximum, in 46 to 47 s.
- Reference against transform plus kernel at w = (0.2, 0.8): 0 of 200,000 decision times differ. PFS is declared at the final analysis, on data that were complete at the interim, in 2.2% of trials. That is the recycling case the record describes.

Supplying the paper's boundary did not work as the record suggests; see section 5 item 11.

### 4.3 Look-back (`07_features.R`, part B)

The three-hypothesis design of section 2 with sequential in place of repeated p-values.

| Graph | Gain, look-back off | Look-back on H2 and H3 | Look-back on all | Trials changed | Rejections lost | Rejections gained | Decided later |
|---|---:|---:|---:|---:|---:|---:|---:|
| Optimised graph of section 2 | 0.66696 | 0.66715 | 0.66715 | 0.08% | 0 | 97 | 0 |
| Holm | 0.61792 | 0.61808 | 0.61811 | 0.08% | 0 | 91 | 0 |

Look-back never removed or delayed a rejection, as it should not. With two analyses it changes fewer than 1 trial in 1,000. FWER under the global null (Holm, 1e6 trials, s.e. about 0.00015) is 0.02226 with and without it.

### 4.4 Mixed spending, unequal information fractions, no interim data (`07_features.R`, part C)

H1 with `sfLDOF` at (0.5, 1), H2 with Pocock-type spending at (0.6, 1), H3 with no data at the interim; powers, correlation and gain as in section 2.

- The simulated correlations between test statistics match rho x sqrt(min t / max t): 0.460 against 0.456, 0.354 against 0.354, 0.501 against 0.500 and 0.385 against 0.387.
- H3 is `NA` at the interim in the raw array, 1 after the transform, and is never rejected there.
- Optimised in 97 s: gain 0.6552, w = (0.952, 0.048, 0), with H1 passing 83% of its level to H2 and 17% to H3.
- Reference against transform plus kernel: 0 of 300,000 decision times differ on the optimised graph and 0 of 300,000 on Holm.

### 4.5 Three analyses, constraints, four hypotheses, reproducibility (`08_scaling_constraints.R`)

| Run | Result | Check | Time |
|---|---|---|---:|
| 2 hypotheses, analyses at 1/3, 2/3 and 1, discount `c(1, 0.9, 0.8)` | w1 = 0.180, gain 0.69231 | grid maximum 0.69230 at w1 = 0.15, plateau [0.09, 0.28]; graphicalMCP on 200 trials: 0 disagreements in 400 reject or retain decisions and in the decision times of 346 rejections | 105 s |
| Section 2 design with H1 fixed as gatekeeper and H2, H3 fixed to recycle to each other | gain 0.66697, H1 split 0.927 and 0.073 | fixed entries respected; the free optimum of section 2 was 0.66696 | 110 s |
| 4 hypotheses, 2 analyses, additive gain with values 0.4, 0.3, 0.2, 0.1 and `d = c(1, 0.8)` | gain 0.6657 | Holm 0.6308, fixed sequence 0.6618 | 129 s |

The timed optimiser runs at `nsim = 1e5` on 2 threads took 46 to 107 s with 2 hypotheses, 97 to 154 s with 3 and 129 s with 4. The same seed gives an identical result when run twice, and 1 thread gives the same result as 2.

### 4.6 The package's own gates

- `dev/gsd_identity_check.R`, run unmodified: `.Random.seed`, `hyp_weight`, `trans_matrix`, `power`, `solution` and `global_output@solution` are `identical()` between `main` (`324b7ca`) and the branch. Only `trial_success$func` differs, as the record says it must.
- Full test suite with gMCPLite, vdiffr and svglite installed: 42 files, 1622 expectations, 0 failures, 0 errors, 3 skips (two for removed functions, one Linux-only). The 30 warnings are ggplot2 deprecation notices from the plot tests on this machine's older ggplot2.
- Without those three packages the same run gives at least 32 errors, not skips (section 5 item 12).

## 5. Things noticed while using it

Numbered roughly by how much they would cost a new user.

1. **Wrong number in the vignette.** `vignettes/articles/trial_success.qmd`, end-to-end chunk: the printed `#> [1] 0.5833` is not what the code gives. Run as written (with the `gain` defined earlier on the page and the stated seed) it returns 0.9028 through the internal kernel and through `calc_power_pvals_gsd()`. With nominal powers near 0.98 and 0.93 a gain near 0.9 is what one would expect; 0.583 looks like a value carried over from the hand-built matrix of the P4 review. The chunks are `eval: false`, so nothing catches this.
2. **The same page says `calc_power_pvals_gsd()` does not exist yet** and routes the reader through `multigrain:::graph_shortcut_gsd()`. It exists and gives the same value.
3. **The same page says a leading zero in a discount table gives no warning.** It warns, as the roxygen says.
4. **The same page says `graph_optimise()` accepts a GSD gain.** It now aborts with a clear message.
5. **`vignettes/articles/trial_success.rmarkdown` is committed.** It looks like a Quarto intermediate of the `.qmd` and carries the same stale text. Suggest deleting it and adding `*.rmarkdown` to `.gitignore`. Not touched here.
6. **`print()` and `summary()` of the result hide the GSD output** (`local_power_by_analysis`, `mean_decision_look`, `time_distribution`). The record notes this (its section 4.7). From a user's side it is the first thing missing: interim power is why one runs a group sequential design.
7. **`spending = gsDesign::sfHSD` fails with R's own "argument "param" is missing, with no default".** The wrapper form is in the help, but the error does not point to it.
8. **Information fractions that never reach 1 are accepted silently**, for example `c(0.4, 0.8)`. That is by design (record, section 8 item 8), but the hypothesis then never spends its whole level and nothing says so. A line in `summary()` would do.
9. **Error messages are otherwise very good.** Raw array into the optimiser, fixed-sample gain into the GSD optimiser, K mismatch, m mismatch, unnamed `spending`, `alpha` above the table, decreasing information fractions: each aborts with a message that says what to do.
10. **Grammar limits that matter for section 7.** A table takes only a bare `t<i>`, so `d(pmax(t1, t2))` is refused; there is no negation and no min or max. Comparisons between two times (`t2 <= t1`) do work.

11. **The classical O'Brien-Fleming boundary cannot be supplied the way Appendix B of the record suggests.** `function(a, t) cumsum(gsDesign::gsDesign(k = length(t), test.type = 1, alpha = a, sfu = "OF", timing = t)$upper$spend)` aborts the transform twice over. First, the check that a spending function spends its whole level calls it with a single information fraction of 1, and `gsDesign(k = 1)` refuses ("input timing of interim analyses must be increasing strictly between 0 and 1"). Second, with that case guarded, `gsDesign()` fails with "f() values at end points not of opposite sign" at 130 of the 1024 table levels, all between 4.0e-14 and 2.8e-11. Guarding the single-look case and falling back to `sfLDOF` at the failing levels works (transform in 13 s) and reproduces the record's boundaries of 0.008197 and 0.022321; `07_features.R` has the wrapper. Either the table builder could catch a failing level and say which one, or the workaround could be documented.
12. **Tests that need a Suggests package error when it is missing.** Without gMCPLite, vdiffr and svglite the full suite gives at least 32 errors (the reporter stops listing them after ten) in six fixed-sample test files (`test-RcppExports.R`, `test-calc_power.R`, `test-objective_function.R`, `test-optimisation_start.R`, `test-plot_graph_constraint.R`, `test-plot_graph_optimal.R`); the ten that are printed all read "there is no package called 'gMCPLite'". With the packages installed all of them pass. `skip_if_not_installed()` is used in the GSD tests and in parts of `test-RcppExports.R` but not in these. Not a GSD matter; noted since it is what a new contributor sees first.
13. **An unnamed list of spending functions prints deparsed code in `summary()`**, for example `list(sfLDOF, pocock, sfLDOF)[[2]]`. The record has this (section 10 item 9c); a named list prints cleanly.

Items 1 to 4 are fixed in the second commit of this pull request.

## 6. Discount, spending function and the shape of the gain (`05_discount_sensitivity.R`)

Two hypotheses, information fractions 0.5 and 1, correlation 0.5, `nsim = 2e5`, full recycling, grid of step 0.02 over w1. `delta` is the value kept by a final-analysis claim. The first row is the design of section 1 on a different set of draws and a coarser grid; its optima (0.28 and 0.12) differ from those of section 1 (0.24 and 0.15) by less than the width of the plateau.

| Design | Spending | w1 at delta = 1 | 0.8 | 0.5 | 0.2 |
|---|---|---:|---:|---:|---:|
| Power 0.80, 0.90; value 0.6, 0.4 | LDOF | 0.28 | 0.12 | 0.02 | 0.00 |
| Power 0.80, 0.90; value 0.6, 0.4 | Pocock | 0.22 | 0.22 | 0.36 | 0.36 |
| Power 0.85, 0.85; value 0.5, 0.5 | LDOF | 0.50 | 0.50 | 0.12 | 0.02 |
| Power 0.85, 0.85; value 0.5, 0.5 | Pocock | 0.54 | 0.54 | 0.46 | 0.52 |

Under LDOF a stronger discount drives the optimum to a corner. Under Pocock it barely moves (all four Pocock optima sit inside each other's plateaux). The reason is the interim boundary as a function of the level held:

| Level | 0.025 | 0.0125 | 0.00625 |
|---|---:|---:|---:|
| LDOF interim boundary | 0.001525 | 0.000412 | 0.000110 |
| Pocock interim boundary | 0.015503 | 0.007751 | 0.003876 |

Halving the level divides the LDOF interim boundary by 3.7 and the Pocock one by 2. Splitting alpha between hypotheses therefore costs much more than proportionally at an LDOF interim, and a gain that rewards interim claims prefers to keep the level in one place.

In the symmetric LDOF design at `delta = 0.5` the profile over w1 is two-peaked: 0.4771 at w1 = 0.10, 0.4756 at the Holm point w1 = 0.50, 0.4769 at w1 = 0.90. The dip is small but not noise: the paired difference between w1 = 0.10 and w1 = 0.50 on the same draws is 0.00144 with standard error 0.00027 (z = 5.3). The default start graph is Bonferroni-Holm, which sits in the trough between the peaks (the profile is flat to about 1e-4 from w1 = 0.36 to 0.60, and by the symmetry of the design w1 = 0.5 is a stationary point). Both `global_search = FALSE` and the default search left it (w1 = 0.92 and 0.81, gains 0.47702 and 0.47685, both on the upper peak), so I have no failure to report, only the observation that time discounting can make this objective non-concave in the weights where the undiscounted one is not. **[simulated]**

## 7. Stopping criteria (`04_stopping_rules.R`)

### 7.1 What the gap is

The kernel analyses every look of every trial; it leaves the look loop early only when all m hypotheses are rejected. So `t<i>` is the decision time in a trial that always continues, and three things follow:

- a hypothesis can be counted as rejected at a look the trial would never have reached;
- `d(t<i>)` assumes each claim can be acted on at its own decision look, whether or not the trial is still running;
- the package has no trial-level time, so a gain cannot say what stopping early is worth.

Expected utility needs three ingredients: which claims are made, when each is realised, and when the trial ends. Today the gain sees the first, a version of the second, and not the third.

### 7.2 Error control

Stopping early can only remove rejections from the always-continue procedure, so the Maurer and Bretz procedure keeps its FWER control under any stopping rule, data dependent or not. **[proved]**: the rejection set of the stopped procedure is a subset of the rejection set of the full one.

What stopping does not allow is giving a hypothesis its full level at the look where the trial stops. Hung, Wang and O'Neill (2007) showed that testing the secondary at level alpha whenever the primary is significant inflates the FWER, and Glimm, Maurer and Bretz (2010) bounded the inflation. In the graphical group sequential procedure each hypothesis has its own spending function, which is what makes it valid; the consequence is that the only lever for a hypothesis the trial may leave behind is that spending function, and the package treats it as a fixed input.

### 7.3 An efficacy stopping rule is an edit of the decision-time matrix

Take any rule of the form "stop after the first look k at which C holds", where C depends only on the set of hypotheses rejected at looks up to k. Let S be that look (K if C never holds). The decision times of the stopped trial are those of the always-continue trial with every time above S set to 0.

**[proved]**: the kernel state after look k (levels, working graph, rejection flags) depends only on the repeated p-values at looks up to k, so rejections at looks up to S are the same whether or not the trial continues, and a stopped trial has no later ones. Look-back and matured endpoints are already encoded in the transformed array, so the argument covers them.

**[simulated]**: three implementations of "stop once H1 is rejected" agree on a 2-hypothesis design (powers 0.90 and 0.80, correlation 0.5, `nsim = 2e5`), for two spending choices:

- (a) the kernel's time matrix edited in R as above;
- (b) the rule written in today's grammar, `0.6 * d(t1) + 0.4 * ((t2 == 1) + 0.8 * (t2 == 2 && t1 != 1))`;
- (c) the reference of `00_reference.R` stopping for real on raw p-values.

(a) and (b) give the same gain to within 2e-12, and (a) and (c) agree on all 200,000 decision times checked (the first 50,000 trials under each spending choice).

So efficacy stopping needs no change to the kernel or the transform. It can be a post-step on the time matrix, or equivalently a `break` in the kernel's loop over looks, which also saves the later passes. graphicalMCP treats a stopped trial the same way: in its diabetes case study the trial stops after the second of three planned analyses, and only those two analyses are passed to the test (reproduced in section 4.1).

### 7.4 What a stopping rule does to the answer

Same design, three value models:

- **M0**: each claim is worth its value at its own decision look and the trial always continues. This is the gain of sections 1 and 2 and of Example 5.
- **M1**: the same, but the trial stops once H1 is rejected, so H2 cannot be rejected afterwards.
- **M2**: the trial stops once H1 is rejected and every claim is realised when the trial ends, so the discount applies to the stopping look and not to each decision look.

Optimal w1 on a grid of step 0.02 with full recycling:

| Spending (H1, H2) | Model | w1 | Plateau (within 5e-4) | Gain | Holm | Fixed sequence | Power H1 | Power H2 | P(stop at interim) |
|---|---|---:|---|---:|---:|---:|---:|---:|---:|
| LDOF, LDOF | M0 | 0.98 | [0.94, 0.98] | 0.7081 | 0.6937 | 0.7074 | 0.896 | 0.754 | 0 |
| LDOF, LDOF | M1 | 0.32 | [0.22, 0.40] | 0.6707 | 0.6693 | 0.6619 | 0.854 | 0.735 | 0.130 |
| LDOF, LDOF | M2 | 0.44 | [0.30, 0.52] | 0.6671 | 0.6667 | 0.6619 | 0.865 | 0.712 | 0.151 |
| LDOF, Pocock | M0 | 0.88 | [0.78, 0.94] | 0.7017 | 0.6981 | 0.6976 | 0.890 | 0.707 | 0 |
| LDOF, Pocock | M1 | 0.46 | [0.40, 0.64] | 0.6903 | 0.6902 | 0.6812 | 0.861 | 0.704 | 0.195 |
| LDOF, Pocock | M2 | 0.98 | [0.92, 1.00] | 0.6816 | 0.6755 | 0.6812 | 0.896 | 0.649 | 0.248 |

All **[simulated]**, one design.

**Valuation.** With LDOF for both, the graph optimised under M0 (w1 = 0.98) is reported at 0.7081. If the trial in fact stops when H1 is rejected, that graph is worth 0.6633: the reported value is about 7% too high. Its H2 power is reported as 0.754 and is 0.614.

**Selection.** The best graph under M1 has w1 = 0.32 and scores 0.6707, so choosing the graph while ignoring the rule costs 0.0074, about 1.1% of the gain. The argmax moves a long way and the regret is modest, the same pattern as the correlation work, though here the regret is above the 0.8% benchmark used there.

**Why the optimum moves away from the primary.** Under M1 nothing rewards stopping except the undiscounted H1 claim, and stopping costs H2 its final analysis. The optimiser responds by making an interim rejection of H1 less likely (0.25 down to 0.13). Whether that is the right behaviour depends on what stopping is worth. Adding a saving c for each trial that stops at the interim, `+ c * (t1 == 1)` on a scale where the claims total 1:

| Saving c | 0 | 0.02 | 0.05 | 0.10 | 0.20 |
|---|---:|---:|---:|---:|---:|
| Optimal w1, LDOF and LDOF | 0.32 | 0.32 | 0.44 | 0.94 | 0.98 |
| Optimal w1, LDOF and Pocock | 0.46 | 0.54 | 0.56 | 0.82 | 0.92 |

The optimal weight on the primary runs from 0.32 to 0.98 across plausible values of c. A stopping rule without a value for stopping gives the optimiser the wrong incentive.

### 7.5 The spending function becomes a first-order choice

Under M0, giving H2 a Pocock-type function lowers the best achievable gain (0.7081 to 0.7017). Under M1 it raises it (0.6707 to 0.6903, about 3%), and that is 14 times the gap between the M1 optimum and Holm under LDOF (0.0014). Under a stopping rule, choosing the secondary's spending function matters more than optimising the graph.

At the fixed sequence w1 = 1, the setting of Tamhane, Mehta and Liu (2010), H2 power under the stopping rule is 0.606 with LDOF and 0.644 with Pocock, which agrees with their conclusion that an O'Brien-Fleming boundary for the primary with a Pocock boundary for the secondary does best. Without the stopping rule the order reverses (0.748 against 0.695). These four powers are in the saved profiles of `04_stopping_rules_results.rds`, not in the log.

The record puts joint optimisation of spending out of scope and names "one transformed array per candidate design and an outer loop over a discrete set of designs" as the extension path (its section 4.5). With stopping I think that outer loop is needed for the answer to be useful. It is cheap: one transform per spending assignment (1.6 to 2.6 s here) on the same raw draws, so common random numbers hold across the assignments.

### 7.6 Futility

A non-binding futility rule depends on the raw statistics and not on the graph, and it leaves the efficacy boundaries unchanged. It can therefore be applied once, in `transform_pvalues_gsd()`, which is the only function that sees raw p-values: for each trial, set every repeated p-value after its futility look to 1. The kernel and the search are untouched. **[proved]** by the same argument as 7.3.

**[simulated]** in part 6 of `04_stopping_rules.R`, with the rule "stop after the interim if the raw interim p-value of H1 exceeds 0.3", which stops 3.9% of trials. Masking the transformed array agrees with the reference stopping those trials for real on every decision time checked: 100,000 for each of the two spending choices, with the futility rule alone and combined with the efficacy rule of 7.3 (400,000 in all). In this example the rule lowers the best gain by 0.009 to 0.011 and leaves the optimal w1 where it was under both M0 and M1.

Two consequences. The time a futility rule saves does not depend on the graph, so on its own a value for that time adds a constant to the objective and cannot move the argmax; it matters only for reporting the value of the design. Combined with an efficacy rule the two meet in S, the earlier of the two stopping looks, which has to be computed per trial. The rejections a futility rule removes do depend on the graph, so the gain changes and in principle the optimum can move, though it did not in the example above. A binding rule would change the efficacy boundaries and would have to enter the boundary tables; I would leave it out.

### 7.7 A possible shape for it

1. A stopping rule as its own object, a logical expression in `r<i>` parsed by the existing gain parser, for example `stop_when_gsd(r1)` or `stop_when_gsd(r1 && r2)`, passed to `graph_optimise_gsd()` and `calc_power_pvals_gsd()` and not buried in the gain. If it lives only in the gain, as in (b) above, every reported power field ignores it: in the end-to-end run of `04_stopping_rules.R` the result object reports H2 power 0.722 where the power under the rule is 0.695.
2. Applied after the kernel (or as a `break` inside it): compute S per trial, zero the later times, and append S as column m + 1 of the time matrix so the compiled gain keeps its single `IntegerMatrix` signature.
3. A new gain symbol `s` for the stopping look, usable in tables and comparisons: `d(s)` for trial-level timing (M2), `cost(s)` for the cost of running to look s. Without it M2 needs one comparison term per combination of decision times, which grows as (K + 1)^m.
4. `calc_power_pvals_gsd()` reports the distribution of S and its mean alongside `time_distribution`. The calendar-time metadata of open item 12 becomes more useful once S exists, since expected duration is the natural summary.
5. Optional `futility` argument on `transform_pvalues_gsd()` as in 7.6, with the per-trial futility look stored on the object so that S is the earlier of the two.
6. An outer loop, or just a documented pattern, over a short list of spending assignments.

### 7.8 Questions

1. Which trials is this for? If the trial continues after the primary succeeds so that the other endpoints can mature (Example 5, and one of the two cases in Glimm et al.), M0 is right and there is no gap. If it stops, M1 or M2 applies.
2. Is a claim worth its value at its own decision look, or when the trial reports? The first needs something like an interim filing while the trial continues.
3. Should stopping have a value of its own? Section 7.4 says the optimal graph depends on it heavily.
4. Is a discrete choice of spending function per hypothesis in scope once stopping is?
5. Is futility wanted in the first version?

### Limits of this section

Two hypotheses, two looks, one set of powers, one correlation, one efficacy rule and one futility rule. The direction of each effect has a mechanism behind it, but the sizes are for this design only.

## 8. Suggestions, in the order I would take them

1. **A stopping rule as its own input, and a stopping-look symbol in the gain** (7.7). Without it the gain cannot say what stopping is worth and the reported power ignores the rule.
2. **Show the group sequential output in `print()` and `summary()`**: power by analysis and the decision-time distribution. They are computed already and are the first thing a user of a group sequential design looks for.
3. **Report the nominal boundaries of the final graph** (the record's open item 5). After optimising, the next thing a protocol needs is the boundary at each analysis for each weight a hypothesis can hold, as graphicalMCP prints with `verbose = TRUE`. The tables on the `multigrain_pvals_gsd` object already hold what is needed.
4. **A documented loop over spending assignments**, reusing one set of raw draws (7.5). With a stopping rule this choice moved the gain more than the graph did.
5. **Calendar times on the p-value object** (open item 12), so that the expected stopping time and a discount tied to dates can be reported. Worth more once item 1 exists.
6. **One end-to-end article** (P5) that walks Example 5 through simulate, transform, optimise and report, with a short table mapping graphicalMCP's conventions to multigrain's (section 4.1).
7. **Small things**: section 5 items 5, 7, 8, 11, 12 and 13.

## Reproducing

From the repository root, with the branch installed:

```
Rscript dev/review/gsd_user_testing/01_two_hyp_two_stage.R
Rscript dev/review/gsd_user_testing/02_three_hyp_two_stage.R
Rscript dev/review/gsd_user_testing/03_verify.R        # needs the results of 01 and 02
Rscript dev/review/gsd_user_testing/04_stopping_rules.R
Rscript dev/review/gsd_user_testing/05_discount_sensitivity.R
Rscript dev/review/gsd_user_testing/06_graphicalmcp_vignette.R
Rscript dev/review/gsd_user_testing/07_features.R        # needs the results of 02
Rscript dev/review/gsd_user_testing/08_scaling_constraints.R   # needs the results of 02
```

`GSD_TEST_THREADS` sets the number of threads for 01, 02, 04, 05, 07 and 08 (default 2); 03 and 06 use the serial kernel. 03, 04 and 07 source `00_reference.R`; 06 needs graphicalMCP; check D of 03 and part D of 08 use it if it is installed. Total run time here was about 40 minutes.

## References

- Glimm E, Maurer W, Bretz F (2010). Hierarchical testing of multiple endpoints in group-sequential trials. *Statistics in Medicine* 29(2):219-228. doi:10.1002/sim.3748
- Hung HMJ, Wang SJ, O'Neill R (2007). Statistical considerations for testing multiple endpoints in group sequential or adaptive clinical trials. *Journal of Biopharmaceutical Statistics* 17(6):1201-1210. doi:10.1080/10543400701645405
- Maurer W, Bretz F (2013). Multiple testing in group sequential trials using graphical approaches. *Statistics in Biopharmaceutical Research* 5(4):311-320.
- Tamhane AC, Mehta CR, Liu L (2010). Testing a primary and a secondary endpoint in a group sequential design. *Biometrics* 66(4):1174-1184. doi:10.1111/j.1541-0420.2010.01402.x
