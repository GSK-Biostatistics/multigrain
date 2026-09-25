# Follow-up fixes for `graph_simplify()`

This records the changes made after the adversarial review in
`dev/sparsity_review.md`.

## Decisions

### Stage-2 `mutation` and `suggestions`

`graph_simplify()` reserves both `GA::ga()` options while its global search is
active. The stage-2 design depends on a support-changing mutation and on a
reference-first seed population; allowing either to be replaced would remove
those guarantees.

An explicit `control` that sets either option now fails before the GA runs and
names the reserved settings. If the default control inherited from the
stage-1 result contains them, `graph_simplify()` warns and removes them. With
`global_search = FALSE`, no GA runs, so they are left alone. `graph_optimise()`
continues to honour both settings unchanged.

### Restoring trial-success functions

The recovery logic lives in `R/trial_success.R`, not solely in
`graph_simplify()`, because `$objective`, `$m`, and `$func` form one
`multigrain_trial_success` object and the helper is reusable by future callers.

A live function is verified semantically rather than only by comparing
formatted strings. For `m <= 12`, it and the stored expression are evaluated on
every one-row Boolean rejection pattern. This is exact for the requested
`m <= 4` scope, works for legacy objects, does not recompile a live function,
and leaves existing live objects unchanged. Above that limit, the stored
generated C++ source must match the objective and bounded deterministic
boundary, singleton, complement and alternating patterns are compared, avoiding
an exponential cost in ordinary use.

A missing or dead function is rebuilt from `$objective`. The restored measure
is installed on the local graph object and is carried by the returned result.
Missing or invalid objectives, dimension mismatches, failed rebuilds, and live
function/objective disagreements all abort explicitly.

## Other changes

* `print()` and `summary()` retain signed sparsity fields but display a negative
  `gain_loss_fraction` as `gain X% over reference`.
* `dev/sparsity_report.md` no longer cites the missing
  `scratchpad/a1_unchanged.R` or claims its unreconciled 100-case result. It
  records that the described grid contains 120 cases and points to the three
  reproducible spot-checks in the committed review.
* `NEWS.md`, roxygen documentation, and targeted tests were updated for the new
  behavior.

## Deliberately not changed

The existing optimiser subsampling expressions were not touched. Review
finding 2 is being handled in a separate pull request.

The existing parser semantics for unparenthesised mixtures of AND and OR were
also preserved. The restoration helper uses the same parsed expression tree as
the compiler rather than redefining that separate behavior.

## Validation

Validation was limited to the affected files, with `NOT_CRAN=true` and
`testthat::test_file()` as required:

```powershell
$env:NOT_CRAN='true'
Rscript -e "pkgload::load_all(); testthat::test_file('tests/testthat/test-trial_success.R')"
# 135 pass, 0 fail, 0 warning, 0 skip

Rscript -e "pkgload::load_all(); testthat::test_file('tests/testthat/test-graph_simplify.R')"
# 144 pass, 0 fail, 0 warning, 0 skip

Rscript -e "pkgload::load_all(); testthat::test_file('tests/testthat/test-graph_optimal.R')"
# 44 pass, 0 fail, 0 warning, 0 skip
```

The final diff was also checked for whitespace errors and for changes to the
two existing `pvals[sample(nsim), ]` expressions.
