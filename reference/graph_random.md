# Generate random graph

Generate a random graph, consisting of a vector of hypothesis weights
and a transition matrix. The weights and transition matrix are randomly
generated respecting any constraints provided, ensuring that the sum of
the weights equals 1 and each row of the transition matrix sums to 1.

## Usage

``` r
graph_random(m = NULL, graph_constraint = NULL, names = "auto")
```

## Arguments

- m:

  An integer representing the number of hypotheses. It defines both the
  length of the weight vector and the dimensions (`m x m`) of the
  transition matrix. Optional when `graph_constraint` is supplied
  (inferred from constraint dimensions). If both `m` and
  `graph_constraint` are supplied, they must agree.

- graph_constraint:

  An optional graph constraint object created by
  [`graph_constraint()`](https://gsk-biostatistics.github.io/multigrain/reference/graph_constraint.md).
  When supplied, fixed elements are honoured and only free (`NA`)
  positions are randomised.

- names:

  An optional character vector containing hypotheses' names. If not
  provided it defaults to `"auto"` meaning the hypotheses will be
  automatically named `"H1"`, `"H2"`, etc.

## Value

A list containing:

- `hyp_weight`: A numeric vector of length `m` representing the
  generated hypothesis weights.

- `trans_matrix`: A numeric matrix of dimension `m x m` representing the
  generated transition matrix.

## Examples

``` r
# Generate a random graph for 5 hypotheses
random_graph <- graph_random(5)

# print the weight vector
random_graph$hyp_weight
#>         H1         H2         H3         H4         H5 
#> 0.34995750 0.18347771 0.02713369 0.26768173 0.17174937 

# print the transition matrix
random_graph$trans_matrix
#>            H1        H2         H3        H4        H5
#> H1 0.00000000 0.1970907 0.02826882 0.1035610 0.6710795
#> H2 0.33323627 0.0000000 0.29235191 0.0588705 0.3155413
#> H3 0.09436979 0.2198000 0.00000000 0.2900734 0.3957568
#> H4 0.20755374 0.2837108 0.31231181 0.0000000 0.1964237
#> H5 0.22885255 0.2476563 0.19637350 0.3271177 0.0000000

# Generate a random graph respecting constraints
gc <- graph_constraint(
    hyp_constraint = c(0.5, NA, NA),
    trans_constraint = matrix(c(0, NA, NA, NA, 0, NA, NA, NA, 0), 3, 3)
)

random_graph <- graph_random(graph_constraint = gc)
random_graph
#> $hyp_weight
#>        H1        H2        H3 
#> 0.5000000 0.2642026 0.2357974 
#> 
#> $trans_matrix
#>           H1        H2        H3
#> H1 0.0000000 0.4271691 0.5728309
#> H2 0.3382072 0.0000000 0.6617928
#> H3 0.2538808 0.7461192 0.0000000
#> 
```
