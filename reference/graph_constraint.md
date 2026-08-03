# Create a *graph constraint* object for optimisation procedures

A *graph constraint* object defines constraints on the hypothesis weight
vector and transition matrix for optimisation of graph-based multiple
testing procedures.

## Usage

``` r
graph_constraint(
  hyp_constraint = NULL,
  trans_constraint = NULL,
  ...,
  names = "auto",
  diagnose = FALSE,
  tolerance = sqrt(.Machine$double.eps)
)
```

## Arguments

- hyp_constraint:

  A numeric vector defining the constraints on the hypothesis weight
  vector. If `NULL`, all hypothesis weights are free parameters to be
  optimised.

- trans_constraint:

  A numeric matrix defining the constraints on the transition matrix
  between hypotheses. If `NULL`, all transition matrix weights are free
  parameters to be optimised, except for diagonal elements which remain
  set to 0.

- ...:

  These dots are for future extensions and must be empty.

- names:

  An optional character vector containing hypotheses' names. If not
  provided it defaults to `"auto"` meaning the hypotheses will be
  automatically named `"H1"`, `"H2"`, etc.

- diagnose:

  A logical value enabling detailed diagnosis. Default is `FALSE`.

- tolerance:

  numeric \>= 0. Differences smaller than `tolerance` will not be
  reported. The default value is close to `1.5e-8` - i.e.
  `sqrt(.Machine$double.eps)` (the standard R definition of "practically
  equal", as used by
  [`base::all.equal()`](https://rdrr.io/r/base/all.equal.html)).

## Value

A multigrain *graph constraint* object (an S3 list with class
`multigrain_graph_constraint`) containing:

- `hyp_constraint`: a numeric vector representing the constraints on the
  hypothesis weight vector.

- `trans_constraint`: a numeric matrix representing the constraints on
  the transition matrix. If an element is `NA`, it is a free parameter
  to be optimised by
  [`graph_optimise()`](https://gsk-biostatistics.github.io/multigrain/reference/graph_optimise.md).

## Details

The object contains a constraint on the weights (`hyp_constraint`) and a
constraint on the transition matrix (`trans_constraint`).

The *graph constraint* object is used to define constraints on both the
hypothesis weight vector and the transition matrix in graph-based
optimisation procedures. The
[`graph_optimise()`](https://gsk-biostatistics.github.io/multigrain/reference/graph_optimise.md)
function will read the graph constraints and only optimise free
parameters (specified by `NA` in `hyp_constraint` and
`trans_constraint`).

Either `hyp_constraint` or `trans_constraint` must be provided (they
can't both be `NULL` at the same time). If only one is provided, the
*graph constraint* object will allow
[`graph_optimise()`](https://gsk-biostatistics.github.io/multigrain/reference/graph_optimise.md)
to optimise any parameter (i.e., no constraints will be specified) in
the other one.

`hyp_constraint` and `trans_constraint` must refer to the same number of
hypotheses.

## References

Xi, D. and Chen, Y. (2024). Optimal weighted Bonferroni tests and their
graphical extensions. *Statistics in Medicine*, 43(3), 475–500.
<https://doi.org/10.1002/sim.9958>

## Examples

``` r
# Create a graph_constraint object with predefined weight constraints
graph_constraint(hyp_constraint = c(NA, 0.4, NA))
#> <multigrain_graph_constraint>
#> Constraints on hypothesis weights:
#>  H1  H2  H3 
#>  NA 0.4  NA 
#> 
#> Constraints on transition matrix:
#>    H1 H2 H3
#> H1  0 NA NA
#> H2 NA  0 NA
#> H3 NA NA  0
```
