# Check the validity of a graph-based MTP

Check that a given weight vector and transition matrix satisfy several
conditions to determine if they form a valid graph. It verifies that the
matrix is square, the dimensions match the length of the weight vector,
the weight vector elements and matrix elements lie within the interval
\[0, 1\], and that the weight vector sums to 1. Additionally, it checks
the matrix is properly normalised, with all diagonal elements equal to 0
and each row summing to 1.

## Usage

``` r
is_graph_valid(
  hyp_weight,
  trans_matrix,
  sum_to_one_constraint = TRUE,
  tolerance = sqrt(.Machine$double.eps)
)
```

## Arguments

- hyp_weight:

  A numeric vector representing the weights of the nodes (hypotheses) in
  the graph. All elements must be in the interval \[0, 1\] and the
  vector must sum to 1.

- trans_matrix:

  A numeric matrix representing the transition matrix between nodes. All
  elements must be in the interval \[0, 1\], with diagonal elements
  equal to 0, and each row must sum to 1.

- sum_to_one_constraint:

  A logical indicating whether to allow graphs where both the transition
  matrix rows are not constrained to sum-to-one, for example in a fixed
  sequence. Defaults to `TRUE`.

- tolerance:

  numeric \>= 0. The tolerance when evaluating the sum-to-one
  constraints of the hypothesis weights and transition matrix rows. The
  default value is close to `1.5e-8` - i.e. `sqrt(.Machine$double.eps)`
  (the standard R definition of "practically equal", as used by
  [`base::all.equal()`](https://rdrr.io/r/base/all.equal.html)).

## Value

A logical value: `TRUE` if the graph is valid, otherwise `FALSE`. The
function issues warnings if any of the validity checks fail.

## Details

The function performs the following checks to determine graph validity:

- `trans_matrix` is a square matrix.

- the length of `hyp_weight` matches the number of rows and columns in
  `trans_matrix`.

- all diagonal elements of `trans_matrix` are 0.

- all elements of `hyp_weight` and `trans_matrix` lie in the interval
  \[0, 1\].

- the sum of `hyp_weight` equals 1.

- each row of `trans_matrix` sums to 1.

If any of these conditions are not satisfied, the function returns
`FALSE` and issues the corresponding warning message.

## Examples

``` r
hyp_weight <- c(0.4, 0.3, 0.3)
trans_matrix <- matrix(
                    c(0, 0.5, 0.5, 0.3, 0, 0.7, 0.6, 0.4, 0),
                    nrow = 3,
                    byrow = TRUE
                )
is_graph_valid(hyp_weight, trans_matrix)
#> [1] TRUE
```
