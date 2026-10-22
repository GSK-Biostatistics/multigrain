# multigrain (development version)

## New functionality

## Minor improvements and bug fixes

* The global search in `graph_optimise()` now uses a split-merge mutation
  instead of the per-parameter truncated Cauchy mutation. Each mutation works
  on one decoded weight vector or transition row: it either merges two edges
  (the smaller gives its weight to the larger) or moves weight between two
  edges by a reflected Cauchy step. Every mutated graph is valid and differs
  from its parent, and no edge is favoured by where it is stored. Along with
  this, the genetic algorithm's defaults are now `pcrossover = 0` and
  `pmutation = 1`. In benchmarks the new search is as good as the old on small
  problems and better on larger ones, at the same cost. Results for a given
  seed differ from multigrain 1.0.0.

* The graphical shortcut used to score candidate graphs (in `graph_optimise()`,
  `calc_power_pvals()` and pruning) now pre-computes the local critical values
  for every intersection once per graph and looks them up in each simulated
  trial, instead of updating the graph after every rejection. This is faster
  for small and moderate numbers of hypotheses and gives the same rejections;
  the previous algorithm is still used for large numbers of hypotheses. There
  is no user-visible API change. The approach follows
  `graph_test_shortcut_fast()` in the graphicalMCP package (Xi and Brockmann,
  2024).

# multigrain 1.0.0

* Initial release.
* Added a NEWS.md file to track package changes.
