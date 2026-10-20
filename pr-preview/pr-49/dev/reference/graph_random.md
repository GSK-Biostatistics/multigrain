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
#>         H1         H2         H3         H4         H5 
#> 0.35071387 0.06866372 0.34471708 0.17832142 0.05758390 

# print the transition matrix
random_graph$trans_matrix
#>           H1        H2         H3         H4        H5
#> H1 0.0000000 0.4591870 0.04061141 0.03885603 0.4613455
#> H2 0.1263199 0.0000000 0.23973113 0.40756689 0.2263821
#> H3 0.1590820 0.1849153 0.00000000 0.21572220 0.4402805
#> H4 0.2370370 0.1164889 0.39764040 0.00000000 0.2488337
#> H5 0.3913403 0.3328453 0.17856518 0.09724917 0.0000000

# Generate a random graph respecting constraints
gc <- graph_constraint(
    hyp_constraint = c(0.5, NA, NA),
    trans_constraint = matrix(c(0, NA, NA, NA, 0, NA, NA, NA, 0), 3, 3)
)

random_graph <- graph_random(graph_constraint = gc)
random_graph
#> $hyp_weight
#>        H1        H2        H3 
#> 0.5000000 0.2573798 0.2426202 
#> 
#> $trans_matrix
#>           H1        H2         H3
#> H1 0.0000000 0.9510283 0.04897169
#> H2 0.2635302 0.0000000 0.73646982
#> H3 0.6601752 0.3398248 0.00000000
#> 
```
