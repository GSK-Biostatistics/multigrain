# Changelog

## multigrain (development version)

- New
  [`family_spec()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/family_spec.md)
  and `families` argument of
  [`graph_optimise()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/graph_optimisation.md)
  and
  [`calc_power_pvals()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/calc_power_pvals.md):
  optimise graphs whose nodes are families of hypotheses (Maurer and
  Bretz, 2014). Hypotheses within a family are tested with Hommel’s,
  Hochberg’s or Holm’s procedure, and a family passes its weight on only
  when all of its hypotheses are rejected. With `families`,
  [`calc_power_pvals()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/calc_power_pvals.md)
  also returns `family_power`. See the new article “Families of
  hypotheses”.

## multigrain 1.0.0

CRAN release: 2026-10-18

- Initial release.
- Added a NEWS.md file to track package changes.
