# Exercises: seeing under the hood of `transform_pvalues_gsd()`

Companion to `dev/callflows/transform_pvalues_gsd.md` (the prose) and
`dev/diagrams/out/gsd-workflow.excalidraw` (the picture). The callflow document
tells you *what* each function does. These exercises are for building the
intuition to explain *why* it works, by doing each step by hand on a toy problem
and then checking it against the package.

Run each script from the repo root with `Rscript dev/review/NN_name.R`, or
source it interactively and stop at each `# >>>` marker to answer the question
in the comment before running on. Each script prints what it computes and ends
with a `stopifnot()` block so you know the claim held on your machine.

Every script starts with `devtools::load_all()`. They need `gsDesign` and
`mvtnorm` (both already installed).

## The one idea to hold on to

A group sequential boundary at look *k* is a function of the alpha allocated to
the hypothesis: give it more alpha, the boundary rises. The graphical procedure
keeps changing that allocation as it recycles. Rather than recompute the boundary
every time the allocation changes, ask the inverse question once: **"what is the
smallest allocation at which this p-value would have crossed?"** That number is
the repeated p-value. Comparing it with the current allocation is the same
decision as comparing the raw p-value with the boundary at that allocation, and
it is a plain fixed-sample comparison. Everything in the transform is machinery
for computing that inverse fast and safely.

## The exercises

| # | Script | What you should be able to say afterwards |
|---|--------|-------------------------------------------|
| 01 | `01_one_look.R` | Why, with a single look at full information, the boundary *is* the level; why the code short-circuits instead of trusting `gsBound1()`. |
| 02 | `02_forward_map_by_hand.R` | How a spending function plus information fractions becomes nominal boundaries. You compute the look-1 boundary with `qnorm()` and the look-2 boundary with a bivariate normal, then match `gsBound1()`. This is the math the rest of the chain rests on. |
| 03 | `03_boundary_table.R` | What `.gsd_boundary_table()` builds: a monotone table of boundary against level, one column per look, on a log grid. You plot it, see the look-1 underflow, and break the well-ordering check with a deliberately bad spending function. |
| 04 | `04_invert.R` | How `.gsd_invert()` turns a p-value into its level by log-log interpolation; what the two clamps (1 above the table, `1e-14` below) do and why the lower clamp is not 0. You measure interpolation error against `uniroot()`. |
| 05 | `05_transform_by_hand.R` | The whole chain reproduced in ~30 lines of plain R, compared to `transform_pvalues_gsd()` and to Maurer and Bretz (2013) Table 2. |
| 06 | `06_maturity_and_missing.R` | The two conventions: `NA` before data means "cannot reject" (`p^r = 1`); after maturity the value is copied forward. You build the PFS/OS example and watch both happen. |
| 07 | `07_look_back.R` | Sequential p-values as a running minimum, and a case where `look_back = TRUE` changes the decision. |
| 08 | `08_why_it_works.R` | The equivalence theorem checked numerically: for random allocations *w* and random p-values, `p <= boundary(w * alpha)` and `p^r <= w * alpha` agree on every case. Then a type I error simulation under the global null. |

Suggested order is numeric. 02 and 04 are the mathematical core; 08 is the
payoff. 03, 06 and 07 are mostly about conventions and can be skimmed on a
first pass.

## Ideas for a vignette

Collected in `09_vignette_ideas.md`.
