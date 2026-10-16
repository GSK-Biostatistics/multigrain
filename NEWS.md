# multigrain (development version)

* The graphical shortcut used to score candidate graphs (in `graph_optimise()`,
  `calc_power_pvals()` and pruning) now precomputes the local critical values
  for every intersection once per graph and looks them up in each simulated
  trial, instead of updating the graph after every rejection. This is faster
  for small and moderate numbers of hypotheses and gives the same rejections;
  the previous algorithm is still used for large numbers of hypotheses. There
  is no user-visible API change.

# multigrain 1.0.0

* Initial release.
* Added a NEWS.md file to track package changes.
