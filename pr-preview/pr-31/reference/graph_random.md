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
#> 0.1106319 0.1533899 0.1524895 0.4028438 0.1806449 

# print the transition matrix
random_graph$trans_matrix
#>           H1        H2         H3         H4         H5
#> H1 0.0000000 0.2561394 0.33795524 0.02935635 0.37654905
#> H2 0.2284497 0.0000000 0.53047896 0.03013680 0.21093456
#> H3 0.3902660 0.1488752 0.00000000 0.09658506 0.36427375
#> H4 0.4149580 0.1438325 0.40983617 0.00000000 0.03137333
#> H5 0.1979900 0.5556418 0.07636948 0.16999875 0.00000000

# Generate a random graph respecting constraints
gc <- graph_constraint(
    hyp_constraint = c(0.5, NA, NA),
    trans_constraint = matrix(c(0, NA, NA, NA, 0, NA, NA, NA, 0), 3, 3)
)

random_graph <- graph_random(graph_constraint = gc)
random_graph
#> $hyp_weight
#>        H1        H2        H3 
#> 0.5000000 0.1896911 0.3103089 
#> 
#> $trans_matrix
#>           H1        H2        H3
#> H1 0.0000000 0.6793908 0.3206092
#> H2 0.5065631 0.0000000 0.4934369
#> H3 0.7867861 0.2132139 0.0000000
#> 
```
