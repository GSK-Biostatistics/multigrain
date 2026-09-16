# Vignette ideas: explaining repeated p-values

Notes toward a user-facing vignette, drawn from what the exercises made
visible. Working title: **"One conversion, then business as usual: how
multigrain handles group sequential designs"**.

## The narrative arc

The vignette should tell one story, not list features. The story is:

1. **The problem.** In a group sequential trial each hypothesis is tested at
   several looks. The threshold at each look depends on how much alpha the
   hypothesis holds, and the graph keeps changing that as it recycles. Naively
   you would recompute a boundary (a multivariate normal root-find) every time
   the graph moves. That is slow and it tangles the spending machinery into the
   optimiser.

2. **The trick.** Ask the inverse question once: *what is the smallest alpha at
   which this p-value would have crossed?* Call it the repeated p-value. Now
   "p below the boundary at the current allocation" is exactly "repeated p
   below the current allocation", and the graph can run unchanged.

3. **Why it is legitimate.** One line of mathematics (the boundary is monotone
   in the allocation) plus one citation (Maurer and Bretz 2013, section 3.3).
   Exercise 08's zero-disagreement check is the numerical version of this.

4. **What the user actually does.** `simulate_pvalues_gsd()` then
   `transform_pvalues_gsd()` then the same optimiser as before.

5. **The two things that surprise people.** Repeated p-values above alpha come
   back as 1 (exercise 05), and a matured endpoint keeps its value at later
   looks rather than dropping out (exercise 06).

## Toy examples that earned their place

Each of these produced an "oh, I see" moment in the exercises. Any vignette
should reuse at least the first three.

### A. The table you can read by eye (exercise 02, step 5)

Four levels, two looks, LDOF:

```
              look 1    look 2
gamma = 0.001  3.3e-06  0.00100
gamma = 0.005  7.2e-05  0.00498
gamma = 0.0125 4.1e-04  0.01236
gamma = 0.025  1.5e-03  0.02450
```

Then: "your look-2 p-value is 0.02. Read down the look-2 column. The boundary
first exceeds 0.02 somewhere between gamma = 0.0125 and 0.025. So you would
need an allocation of about 0.020 to reject it. That is its repeated p-value."
This is the entire method in one table and one sentence. The interpolation is
just doing this lookup precisely.

### B. Maurer and Bretz Table 2 (exercise 05)

Four hypotheses, raw p at look 1 = 0.0062 but repeated p = 0.114. Readers who
know O'Brien-Fleming will recognise instantly why: at a third of the
information, almost no alpha has been spent, so a raw 0.006 is nowhere near
enough. This example also sets up the "why is it 1?" surprise, and the answer
(the package caps at alpha because nothing above alpha can ever matter)
teaches what the number means.

### C. The naive-versus-correct FWER (exercise 08, steps 2 and 3)

Under the global null with three looks: testing the raw p against 0.025 at
every look gives FWER ~0.054; testing the repeated p against 0.025 gives
~0.025. Two numbers, one plot. This is the "so what" for a reader who does not
care about the internals.

### D. PFS and OS (exercise 06)

The Example 5 shape: PFS mature at look 1, OS accruing through look 3. Show the
`summary()` output with its "Nominal boundaries at the full level" block: PFS
has 0.025 at every look (no sequential penalty, nothing to invert), OS has
0.0004, 0.007, 0.023. Then explain the copy-forward: PFS keeps its repeated
p-value alive so that if OS is rejected at look 3 and recycles alpha, PFS can
still be picked up. This is the convention people get wrong when they roll
their own.

### E. Look-back flips a decision (exercise 07, step 2)

Two hypotheses, one edge, hand-picked p-values where H1 is rejected with
look-back and not without. Small enough to trace on paper. Good for the
"which semantics does your protocol use?" paragraph.

## Figures worth making

- **Boundary curves on log-log axes** (exercise 03's plot). One panel per
  spending family (LDOF, Pocock, HSD -4, HSD +4), three looks each. Readers see
  monotonicity, see how Pocock's look-1 curve sits much higher than LDOF's, and
  see the look-K curve hug the diagonal. The transform is "read this graph
  horizontally".

- **Repeated versus raw p-value**, look by look, for a few thousand simulated
  trials. Points above the diagonal at early looks; on the diagonal at the last
  look. Quantifies the early-look penalty.

- **Interpolation error versus grid size** (exercise 04, step 3). One small
  panel justifying the default `grid_size = 1024`.

- **Decision-time distribution** with and without look-back, from the same
  simulated trials. Shows that look-back moves decisions earlier, never later.

## Sections to keep short or drop

- The floor at `1e-14` rather than 0. One sentence and a footnote to the design
  record. Users never see it.
- The underflow filter in `.gsd_invert()`. Internal.
- The short-circuit and the `0.012499999999999968` story (exercise 01). Fun
  but belongs in developer notes, not a user vignette.
- The well-ordering check. One sentence: "if your spending function is not
  monotone in the level the method is invalid and the transform will say so."

## Open questions before writing

- Where does the vignette sit relative to the existing fixed-sample vignettes
  (`dev/vignettes/`)? It should assume the reader has seen the graph and
  `graph_optimise()` already.
- Should the "why it is legitimate" section carry the one-line proof or just
  the citation? The proof is genuinely one line (the set of allocations that
  reject is an upper set), and readers who want it are the ones who will
  otherwise distrust the method.
- Whether to include a "reproduce it yourself" appendix with `my_transform()`
  from exercise 05. It is 15 lines and removes all mystery. I lean yes.
