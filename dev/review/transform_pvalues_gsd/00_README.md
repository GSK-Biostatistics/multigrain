# Exercises: seeing under the hood of `transform_pvalues_gsd()`

Companion to `dev/callflows/transform_pvalues_gsd.md` (the prose, which
documents what each function does and why) and
`dev/diagrams/out/gsd-workflow.excalidraw` (the picture). These exercises are
for building the intuition to *explain* the method, by doing each step by hand
on a toy problem and then checking it against the package.


## Prerequisites

Every script starts with `devtools::load_all()` and requires `gsDesign` and
`mvtnorm` (exercise 02 uses `mvtnorm` directly; the rest need only `gsDesign`).
Run from the repo root:

```
Rscript dev/review/transform_pvalues_gsd/NN_name.R
```

Or source interactively and stop at each `# >>>` marker to answer the question
in the comment before running on. Each script prints what it computes and ends
with a `stopifnot()` block so every claim is machine-checked.


## The one idea to hold on to

A group sequential boundary at look *k* is a function of the alpha allocated to
the hypothesis: more alpha, larger boundary. The graphical procedure keeps
changing the allocation as it recycles. Rather than recompute the boundary every
time, ask the inverse question once: **"what is the smallest allocation at which
this p-value would have crossed?"** That number is the repeated p-value.
Comparing it with the current allocation is the same decision as comparing the
raw p-value with the boundary at that allocation -- a plain fixed-sample
comparison. Everything in the transform is machinery for computing that inverse
fast and safely.


## The exercises

**First pass (the mathematical core):** 02, 04, 05, 08. These four cover the
forward map, the inversion, the full chain, and the proof that it works.

**Full pass:** all eight in numeric order. Exercises 01, 03, 06, and 07 cover
edge cases and conventions that matter for a correct implementation but are not
the central idea.


### 01 -- `01_one_look.R`: the short-circuit

*Question:* why does a single look at full information need no table?

With one look, the boundary at level *a* is *a* itself, so the repeated p-value
equals the raw p-value. The script shows that `gsBound1()` returns
0.012499999999999968 for level 0.0125 -- not bit-identical to 0.0125 (the
difference is about 2e-17). That is why the code short-circuits: it copies the
raw value through without calling `gsBound1()`, so that the K = 1 path is
bit-identical to the fixed-sample kernel. At the end you will have confirmed
that `identical(raw, tr$pvals)` holds and no boundary table is built.


### 02 -- `02_forward_map_by_hand.R`: the forward map

*Question:* given a spending function, information fractions, and an allocated
level, how do you get the nominal boundary at each look?

You compute the look-1 boundary with `qnorm()` (one-dimensional) and the look-2
boundary with a bivariate normal root-find via `pmvnorm()` and `uniroot()`,
matching `gsBound1()` to within 1e-8. The script then shows boundaries at four
levels (0.001, 0.005, 0.0125, 0.025) to demonstrate monotonicity -- the
property the entire transform depends on. At the end you have a table you can
read by eye: "your look-2 p-value is 0.02; read down the look-2 column; the
boundary first exceeds 0.02 between gamma = 0.0125 and 0.025; that is
approximately the repeated p-value." Exercise 04 makes this lookup precise.


### 03 -- `03_boundary_table.R`: what `.gsd_boundary_table()` builds

*Question:* what does the boundary table look like, and what goes wrong when the
spending function is not well ordered?

The script builds a 128-row table on a log-spaced grid from 1e-14 to 0.025 and
plots both columns on log-log axes. You will see the look-1 underflow: at grid
levels below about 1e-6, the LDOF look-1 boundary saturates at a constant
2.75e-89 (numerical underflow from `gsBound1`). Exercise 04 shows why this does
not matter. The script then feeds a deliberately non-monotone spending function
(one whose look-1 spend *decreases* as the level increases) to
`.gsd_boundary_table()` and confirms it aborts with a "well ordered" error.

*Convention:* the well-ordering check is condition (2) of Maurer and Bretz
(2013). If the boundary is not monotone in the level, the group sequential
graphical procedure itself is invalid.


### 04 -- `04_invert.R`: log-log interpolation and the two clamps

*Question:* how does `.gsd_invert()` turn a p-value into its repeated p-value,
and how accurate is it?

The script inverts several look-2 p-values via the table, checks them against
`uniroot()` (the exact answer), and sweeps over grid sizes to show the
interpolation error trend (40 random p-values, look 2, LDOF):

```
G =  128  max relative error = 5.6e-05
G =  256  max relative error = 9.3e-06
G =  512  max relative error = 3.0e-06
G = 1024  max relative error = 5.3e-07
```

The error roughly quarters each time G doubles (consistent with second-order
interpolation on the log scale). The design record (section 4.1) quotes the
same trend on 200 p-values: about 1e-5 at G = 256, 4e-6 at G = 512, and 1.5e-6
at G = 1024. Two clamp rules handle values outside the table:

- **Above the table** (p larger than the boundary at alpha): the repeated
  p-value is 1, meaning no allocation up to alpha could reject it.
- **Below the table** (p smaller than the smallest positive boundary): the
  repeated p-value is floored at `gsd_grid_min = 1e-14`, not zero. The design
  record (section 4.1) explains why: the kernel's weight-snapping produces
  allocations as small as ~2.5e-11; a zero floor would make any tiny p-value
  rejectable at any positive allocation, causing spurious rejections.

An `NA` raw p-value maps to repeated p = 1 ("cannot reject").


### 05 -- `05_transform_by_hand.R`: the whole chain in 30 lines

*Question:* can you reproduce `transform_pvalues_gsd()` from scratch?

The script writes a `my_transform()` function in about 15 lines of plain R
(build the grid, compute boundaries at each level with `gsBound1`, drop
underflow entries, interpolate on the log scale, clamp). It matches the package
to within 1e-12 on 200 simulated three-look trials.

The script then reproduces Maurer and Bretz (2013) Table 2: four hypotheses,
three looks at information fractions 1/3, 2/3, 1, Lan-DeMets O'Brien-Fleming
spending. **Surprise:** every look-1 repeated p-value and H4 at look 2 come
back as 1, while the paper reports values like 0.1141 and 0.1285. The reason:
the table only runs up to alpha = 0.025, and a repeated p-value above alpha
means the hypothesis cannot be rejected at any allocation the graph could give
it, so the package reports 1. The paper reports the mathematically exact inverse
over (0, 1); for the decision they are equivalent (0.114 > 0.025 and 1 > 0.025
both mean "not rejected"). Widening the grid to alpha = 0.5 recovers the
paper's values to within 2e-4, which is the ADDPLAN-versus-gsDesign numerical
residual (design record, section 4.1).


### 06 -- `06_maturity_and_missing.R`: endpoints that start late or finish early

*Question:* what does the transform do with ragged data -- PFS that matures at
the first analysis, OS that does not start until look 2?

Two conventions (design record, "Two conventions the implementer must know"):

- **No data yet** (info_frac is `NA`): repeated p = 1 ("cannot reject").
- **Matured** (info_frac has reached 1): the repeated p-value is copied forward
  to all later looks, so the hypothesis can still be picked up if recycling
  gives it more alpha.

The script builds a PFS/OS example: PFS is complete at look 1, OS accrues
through three looks. PFS's repeated p-value equals its raw p-value at every
look (single matured look triggers the short-circuit; value is then copied
forward). The script also confirms that supplying *different* raw p-values
after maturity triggers a warning, and the maturity value is used regardless.


### 07 -- `07_look_back.R`: sequential p-values

*Question:* what does `look_back = TRUE` do, and when does it change a
decision?

The sequential p-value at look *k* is the running minimum of the repeated
p-values through look *k*: one `cummin` per hypothesis. The script confirms
this with `all.equal(seq_p, cummin(rep_p))`.

**Surprise:** the script constructs a specific two-hypothesis example where the
decision flips. H1 has strong evidence at look 2 (repeated p about 0.015) but
weaker evidence at look 3 (repeated p about 0.032). H2 is rejected at look 3
and recycles its alpha to H1, which now holds the full 0.025. Without
look-back, H1 is judged on its look-3 evidence (0.032 > 0.025: not rejected).
With look-back, H1 is judged on its best evidence so far
(min(1, 0.015, 0.032) = 0.015 < 0.025: rejected). The `stopifnot()` block
confirms: without look-back `rejected = (FALSE, TRUE)`; with look-back
`rejected = (TRUE, TRUE)`. The decision flips only because the hand-picked
raw p-values put H1's look-2 repeated p-value between 0.0125 and 0.025.


### 08 -- `08_why_it_works.R`: the equivalence theorem and FWER simulation

*Question:* does "raw p below the boundary at allocation *a*" really give the
same answer as "repeated p below *a*"?

The script checks this on 300 random p-values at each of 3 looks, crossed with
20 random allocations: 18,000 decisions compared, zero disagreements. This is
the numerical version of the theorem (Maurer and Bretz, section 3.3): because
the boundary is monotone, the set of allocations that reject is an upper set,
and membership in that set is equivalent to exceeding the repeated p-value.

Then a type I error simulation: 200,000 trials under the global null, three
looks, LDOF spending. Testing the repeated p-values against alpha = 0.025 gives
a rejection rate within Monte Carlo noise of 0.025 (the spending design's
guarantee). Testing the *raw* p-values against 0.025 at every look gives about
0.054 -- the cost of repeated testing that the transform eliminates.


## Two things that will surprise you

1. **Repeated p-values above alpha are reported as 1** (exercise 05). The table
   only covers levels up to alpha; above that, no allocation the graph could
   give would reject the hypothesis, so the package reports 1. Maurer and
   Bretz's Table 2 shows the exact inverse (e.g. 0.114), but 0.114 > 0.025 and
   1 > 0.025 produce the same decision.

2. **The look-back example flips the decision only for specific raw p-values**
   (exercise 07). The effect requires H1's best repeated p-value to fall between
   its initial allocation (0.0125) and the full level (0.025). If it were below
   0.0125, H1 would be rejected either way; if above 0.025, not even look-back
   would help.


## Where to go next

- `dev/callflows/transform_pvalues_gsd.md` documents what every function in the
  transform does, argument by argument, branch by branch.
- `09_vignette_ideas.md` collects notes for a user-facing vignette, including
  which toy examples earned their place and which details to leave out.
