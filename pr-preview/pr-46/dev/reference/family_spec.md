# Create a *family specification* for graphs of families of hypotheses

A *family specification* groups hypotheses into families. When it is
passed to
[`graph_optimise()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/graph_optimisation.md)
or
[`calc_power_pvals()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/calc_power_pvals.md),
each node of the graph is a family rather than a single hypothesis, and
the graph is tested with the procedure of Maurer and Bretz (2014):

## Usage

``` r
family_spec(members, ..., test = "hommel")
```

## Arguments

- members:

  A list of integer vectors, one per family, giving the hypotheses
  (columns of the p-value matrix) in each family. Together they must
  contain each of `1, ..., m` exactly once. If the list is named, the
  names are used as family names; otherwise families are named `"F1"`,
  `"F2"`, and so on.

- ...:

  These dots are for future extensions and must be empty.

- test:

  The within-family test: `"hommel"` (default), `"hochberg"` or
  `"holm"`. Either a single test used for every family, or one test per
  family.

## Value

A `multigrain_family_spec` object containing:

- `members`: a named list with the (sorted) hypotheses in each family.

- `test`: a named character vector with the within-family test of each
  family.

## Details

- within a family, the hypotheses are tested with `test` at the family's
  local significance level;

- a family passes its weight on only when every hypothesis in it is
  rejected.

Families are combined with a weighted Bonferroni test, as the hypotheses
of a standard graph are. Within each family, `test` is applied to that
family's p-values alone:

- `"hommel"`: the Hommel (1988) procedure, i.e. the closed test with
  Simes intersection tests.

- `"hochberg"`: the Hochberg (1988) step-up procedure.

- `"holm"`: the Holm (1979) step-down procedure.

Hommel and Hochberg control the family-wise error rate when the test
statistics within a family are positively dependent, for example
one-sided tests based on multivariate normal statistics with
non-negative correlations. This cannot be checked from the p-values; use
`"holm"` for a family where negative dependence is possible. For a
family with a single hypothesis every test gives the same result.

The hypotheses are the columns of the p-value matrix. In
[`trial_success()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/trial_success.md),
`r1, r2, ...` still refer to these individual hypotheses; the success of
a whole family is written as a conjunction, for example `r2 && r3`.

## References

Maurer, W. and Bretz, F. (2014). A note on testing families of
hypotheses using graphical procedures. *Statistics in Medicine*, 33(30),
5340–5346. [doi:10.1002/sim.6307](https://doi.org/10.1002/sim.6307)

Hommel, G. (1988). A stagewise rejective multiple test procedure based
on a modified Bonferroni test. *Biometrika*, 75(2), 383–386.
[doi:10.1093/biomet/75.2.383](https://doi.org/10.1093/biomet/75.2.383)

Hochberg, Y. (1988). A sharper Bonferroni procedure for multiple tests
of significance. *Biometrika*, 75(4), 800–802.
[doi:10.1093/biomet/75.4.800](https://doi.org/10.1093/biomet/75.4.800)

Holm, S. (1979). A simple sequentially rejective multiple test
procedure. *Scandinavian Journal of Statistics*, 6(2), 65–70.

## See also

[`graph_optimise()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/graph_optimisation.md)
and
[`calc_power_pvals()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/calc_power_pvals.md),
which take the family specification through their `families` argument.

## Examples

``` r
# Five hypotheses in three families: H1 alone, H2 to H4 together and H5
# alone. Hommel is used within every family.
family_spec(list(primary = 1, secondary = 2:4, other = 5))
#> <multigrain_family_spec>
#> Hypotheses (columns of the p-value matrix) in each family:
#> primary   (hommel): 1
#> secondary (hommel): 2 3 4
#> other     (hommel): 5

# One test per family
family_spec(
  list(primary = 1, secondary = 2:4, other = 5),
  test = c("hommel", "hochberg", "holm")
)
#> <multigrain_family_spec>
#> Hypotheses (columns of the p-value matrix) in each family:
#> primary   (hommel): 1
#> secondary (hochberg): 2 3 4
#> other     (holm): 5
```
