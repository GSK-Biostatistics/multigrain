# multigrain (development version)

* New `family_spec()` and `families` argument of `graph_optimise()` and
  `calc_power_pvals()`: optimise graphs whose nodes are families of hypotheses
  (Maurer and Bretz, 2014). Hypotheses within a family are tested with
  Hommel's, Hochberg's or Holm's procedure, and a family passes its weight on
  only when all of its hypotheses are rejected. With `families`,
  `calc_power_pvals()` also returns `family_power`. See the new article
  "Families of hypotheses".

# multigrain 1.0.0

* Initial release.
* Added a NEWS.md file to track package changes.
