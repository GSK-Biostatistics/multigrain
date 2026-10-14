<!--
Provenance: written on 2026-10-14 in phase P5 of dev/hsd_opt_build_plan.md, as the
plain-English companion to design record section 4.12. The numbers of the running
example come from dev/review/hsd_opt/article_example.R (output kept beside it),
the same run as the "Optimising the spending parameter" section of the
trial_success article.
-->

# Choosing a spending function by optimisation, explained

This note explains, without formulas or code, what the spending-parameter
optimisation in `multigrain` does, why it was added and how to read what it
reports. The technical record is section 4.12 of `dev/gsd_design_record.md`.

## What a spending function is

A confirmatory trial is allowed a small chance of a false claim, usually
2.5 %. Think of that as an error budget. A trial with an interim analysis looks
at the data more than once, and every look is another chance to be fooled by
noise. So the budget has to be shared between the looks.

A spending function is the schedule for that sharing: it says how much of the
2.5 % may have been used by the time a given share of the trial's information
has come in. By the final analysis the whole budget must have been spent.

Looking early costs something either way. Spending a lot at the interim makes
an early success easier, but leaves less for the final analysis, where most
of the information is. Spending almost nothing early keeps the final analysis
nearly as strong as in a trial with no interim, but then an early success
needs overwhelming evidence. The classical O'Brien–Fleming boundary takes the
second view, and the Pocock boundary leans towards the first.

## What the HSD parameter does

The Hwang–Shih–DeCani family of spending functions covers that whole range
with one number, its parameter, written gamma. Low values save the budget for
the end and high values spend it early. In formula form, for a share t of
the information, the fraction of the budget spent so far is

> (1 − e^(−gamma·t)) / (1 − e^(−gamma)),

which is just t when gamma is 0 (spending in proportion to information).

The search uses −4 to 1 by default:

- at gamma = −4, the halfway look has used about **12 %** of the budget,
  much like O'Brien–Fleming;
- at gamma = 1 it has used about **62 %**, much like Pocock.

Each hypothesis in a design can have its own spending function, so a trial
with a primary and a secondary endpoint can spend cautiously on one and
generously on the other.

## Why it became a first-order choice

For a trial that always runs to the end, the spending function of a
secondary endpoint is a detail. The package's user testing (section 7.5 of
`dev/review/gsd_user_testing/README.md`) found that it stops being a detail
once the trial **stops as soon as the primary endpoint succeeds**. After that
stop the secondary has no further analysis. Whatever it has not spent by
then, it never gets to use.

In that example, letting the secondary spend early raised the value of the
best design from **0.6707 to 0.6903**. That gain was **14 times** the gap
between the best graph and the textbook Holm graph. In other words, choosing
the secondary's spending function mattered more than optimising the graph. It
also moved the best weight on the primary from **0.32 to 0.46**, so the graph
and the spending function cannot be chosen separately.

## What "optimising it jointly" means

The *graph* of a design decides how the error budget is split between
hypotheses at the start, and how a hypothesis's share passes on to the others
once it succeeds. The spending functions decide how each share is spread over
time. The two interact: the best split depends on how each hypothesis spends,
and the best way to spend depends on the split.

Optimising jointly means choosing the graph and the spending parameter
together, as one design, rather than fixing one and tuning the other. The
user tags each hypothesis whose parameter should be chosen (with
`sf_hsd_opt()`). Hypotheses can share one parameter or have their own, and
the rest keep the spending functions the user gave them.

## How the search works

1. **A short list of candidate values, each checked in advance.** The
   parameter is not searched continuously. Instead there is a short list of
   values: by default −4 to 1 in steps of 0.25, 21 values in all. Before any
   search begins, `transform_pvalues_gsd()` works out the analysis boundaries
   for every value on the list and checks that each gives a valid test. The
   search can only ever try values from the list, so it can never land on an
   unchecked one. This step takes about 20 seconds and is done once; the
   result (a *spec*) can be reused.
2. **A global search over graph and parameter together.** A genetic
   algorithm, the same one `multigrain` already uses for graphs, now also
   carries the parameter as one more "gene". It explores combinations of
   graph and parameter on a subsample of the simulated trials.
3. **The graph re-tuned at every candidate value.** The search then visits
   each of the 21 values in turn, walking outward from the best one found so
   far. At each value it re-optimises the graph, starting from the best graph
   of the neighbouring value. The value with the best result is the answer.
   The 21 results are also the *profile*: how good the design can be at each
   value of the parameter. All of these comparisons use the same simulated
   trials, so the differences between values are not blurred by simulation
   luck.

With two or more parameters the third step goes through them one at a time
and repeats until nothing moves, at most three rounds.

## How to read the result

The running example is the design from the user testing: two hypotheses, an
interim at half the information, and a trial that stops once the primary
succeeds. The secondary's parameter is optimised and the primary keeps an
O'Brien–Fleming-type function. The results below come from 100,000 simulated
trials (vignette `trial_success`, "Optimising the spending parameter").

- **The chosen value.** gamma = **−0.25**, with a design value (gain) of
  **0.6930** and a weight of 0.445 on the primary. This lies between the two
  textbook choices.
- **The profile table** (`summary()` of the result) has one row per candidate
  value. It gives the best gain with the graph re-tuned there, and the gain
  of the chosen graph left as it is. It also gives how far short of the
  optimum the row falls, with a standard error for that shortfall. Here the
  gain falls steadily below −1: at −4 it is 0.0128 lower, about 37 standard
  errors. Above zero it falls more gently.
- **The two plateau intervals** (`print()` of the result) give the range of
  values whose gain is within 0.0005 of the best. The first uses the re-tuned
  graphs, the second the chosen graph kept fixed. Here both run from **−0.75
  to 0**.

A wide plateau means the exact value hardly matters: any value in it is as
good as the optimum to within the precision of the simulation. In practice
the plateau is the more useful answer, since it lets a team pick a value it
can explain. If the re-tuned plateau is clearly wider than the fixed-graph
one, the best graph changes with the parameter, and the graph should be
re-optimised if a different value is adopted.

On the same simulated trials, fixing the secondary's spending at
O'Brien–Fleming type gave 0.6709, and at Pocock type 0.6905. These match the
user testing's 0.6707 and 0.6903. The optimised value adds a further 0.0025
over Pocock. The main lesson stands (spend early on the secondary), and the
optimiser adds that it need not spend quite as early as Pocock.

## What it does not change

- **The error-rate guarantee.** Every value the search can try was checked to
  give a valid group sequential test before the search started. The reported
  design is evaluated on exactly those checked boundaries. Choosing the
  parameter is like choosing the graph: whichever is chosen, the error rate
  is controlled.
- **The simulation engine.** The core routine that runs the graphical test
  on each simulated trial is untouched, as are the gain functions and
  stopping rules.
- **Any result that does not use `sf_hsd_opt()`.** Without a tagged
  hypothesis, every function behaves exactly as before; this is checked by
  comparing results before and after the change.

## Limits

- **Chosen on the same trials that score it.** The parameter and the graph
  are picked because they did best on these simulated trials, so the
  reported value is slightly optimistic. Check the design on freshly
  simulated trials: build a new spec from them and use
  `resolve_pvalues_gsd()` to fix the chosen value. In the running example the
  fresh trials gave 0.6936 against 0.6930 reported, so the optimism was
  smaller than the simulation noise. That will not always be so.
- **Granularity.** Only the values on the list are tried. A step of 0.25 lost
  about 0.00002 against a ten-times finer list in the package's tests, which
  is noise level. If the best value sits at an end of the range, the true
  optimum may lie beyond it. On one three-hypothesis design it lay below −4,
  so widen the bounds (anywhere from −10 to 2 is allowed) and check.
- **Run time.** Checking the 21 values takes about 20 seconds per spec. The
  search itself took about 1.3 minutes on the two-hypothesis example and
  about 8 minutes on a three-hypothesis design, with 100,000 simulated trials
  on two threads. That is 1.4 to 3.5 times a search with the spending fixed.
