# Changelog

## multigrain (development version)

- The graphical shortcut used to score candidate graphs (in
  [`graph_optimise()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/graph_optimisation.md),
  [`calc_power_pvals()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/calc_power_pvals.md)
  and pruning) now precomputes the local critical values for every
  intersection once per graph and looks them up in each simulated trial,
  instead of updating the graph after every rejection. This is faster
  for small and moderate numbers of hypotheses and gives the same
  rejections; the previous algorithm is still used for large numbers of
  hypotheses. There is no user-visible API change. The approach follows
  `graph_test_shortcut_fast()` in the graphicalMCP package (Xi and
  Brockmann, 2024).

## multigrain 1.0.0

CRAN release: 2026-10-18

- Initial release.
- Added a NEWS.md file to track package changes.
