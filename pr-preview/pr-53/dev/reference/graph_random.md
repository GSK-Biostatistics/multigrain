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
  [`graph_constraint()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/graph_constraint.md).
  When supplied, fixed elements are honoured and only free (`NA`)
  positions are randomised.

- names:

  An optional character vector containing hypotheses' names. If not
  provided it defaults to `"auto"` meaning the hypotheses will be
  automatically named `"H1"`, `"H2"`, and so on.

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
#>        H1        H2        H3        H4        H5 
#> 0.2270011 0.2923408 0.1849178 0.1093510 0.1863892 

# print the transition matrix
random_graph$trans_matrix
#>            H1         H2        H3         H4          H5
#> H1 0.00000000 0.47874745 0.4613150 0.05654358 0.003393991
#> H2 0.26753951 0.00000000 0.2897257 0.26253030 0.180204459
#> H3 0.09650060 0.28280445 0.0000000 0.31551864 0.305176308
#> H4 0.07820677 0.36827784 0.4195843 0.00000000 0.133931043
#> H5 0.31038277 0.03020637 0.5575341 0.10187677 0.000000000

# Generate a random graph respecting constraints
gc <- graph_constraint(
    hyp_constraint = c(0.5, NA, NA),
    trans_constraint = matrix(c(0, NA, NA, NA, 0, NA, NA, NA, 0), 3, 3)
)

random_graph <- graph_random(graph_constraint = gc)
random_graph
#> $hyp_weight
#>        H1        H2        H3 
#> 0.5000000 0.2305589 0.2694411 
#> 
#> $trans_matrix
#>            H1        H2        H3
#> H1 0.00000000 0.1296672 0.8703328
#> H2 0.06513572 0.0000000 0.9348643
#> H3 0.53474349 0.4652565 0.0000000
#> 
```
