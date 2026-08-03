# Create an unconstrained *graph constraint*

Create a graph constraint that allows for all hypotheses transitions to
be optimised.

## Usage

``` r
graph_constraint_free(num_hyp, names = "auto")
```

## Arguments

- num_hyp:

  An integer denoting the number of hypotheses. Must be greater than or
  equal to 2.

- names:

  An optional character vector containing hypotheses' names. If not
  provided it defaults to `"auto"` meaning the hypotheses will be
  automatically named `"H1"`, `"H2"`, etc.

## Value

An unconstrained `multigrain_graph_constraint` object. In
`hyp_constraint` all values are `NA` and similarly in the
`trans_constraint`, except for the diagonal which is set to `0`. If we
have only 2 hypotheses, then the matrix will have 0 on the diagonal and
the other 2 elements are set to 1.

## Examples

``` r
# Create a graph constraint object with 3 hypotheses and no constraints
graph_constraint_free(3)
#> <multigrain_graph_constraint>
#> Constraints on hypothesis weights:
#> H1 H2 H3 
#> NA NA NA 
#> 
#> Constraints on transition matrix:
#>    H1 H2 H3
#> H1  0 NA NA
#> H2 NA  0 NA
#> H3 NA NA  0

# Create a graph constraint object with 2 hypotheses results in set values in
# the transition matrix
graph_constraint_free(2)
#> <multigrain_graph_constraint>
#> Constraints on hypothesis weights:
#> H1 H2 
#> NA NA 
#> 
#> Constraints on transition matrix:
#>    H1 H2
#> H1  0  1
#> H2  1  0
```
