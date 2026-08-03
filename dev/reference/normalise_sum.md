# Normalise graph weights to sum to a target value

Ensures that a numeric vector of graph weights is rescaled so its
elements sum to a target value (default 1). This is useful when
exporting optimised graphs (e.g. to gMCP), where weights must form a
valid probability vector.

## Usage

``` r
normalise_sum(
  x,
  ...,
  fixed_idx = integer(0),
  target = 1,
  tolerance = sqrt(.Machine$double.eps)
)
```

## Arguments

- x:

  A numeric vector of weights to be normalised.

- ...:

  These dots are for future extensions and must be empty.

- fixed_idx:

  An integer vector representing the indices of elements that must not
  be modified (e.g. diagonal entries in a transition matrix, or
  positions locked by
  [`graph_constraint()`](https://gsk-biostatistics.github.io/multigrain/dev/reference/graph_constraint.md)).
  Defaults to `integer(0)` (no fixed elements).

- target:

  A number representing the desired sum for the output vector. Default
  is `1`.

- tolerance:

  numeric \>= 0. The tolerance to be used when checking whether the sum
  has converged to `target`. The default value is close to `1.5e-8` -
  i.e. `sqrt(.Machine$double.eps)` (the standard R definition of
  "practically equal", as used by
  [`base::all.equal()`](https://rdrr.io/r/base/all.equal.html)).

## Value

A numeric vector of the same length as `x`, adjusted so that `sum(x)` is
within `tolerance` of `target`.

## Details

Many downstream tools (e.g. gMCP) require hypothesis weights to sum
to 1. Direct division by `sum(x)` may leave tiny discrepancies due to
floating-point rounding. This helper fixes such issues automatically,
ensuring exported graphs are valid.

The algorithm first scales all free (non-fixed) elements proportionally
so they sum to `target - sum(x[fixed_idx])`. It then computes the
largest free element (the "anchor") as the exact complement of all other
elements (`target - sum(x[-anchor])`), absorbing any remaining rounding
residual. A single additive fallback pass handles the rare case where
the prior scaling introduced enough rounding for the complement to land
1 ULP off.

## Note

Two edge cases return early without enforcing `target`:

- If all elements of `x` are zero, the zero vector is returned unchanged
  (proportional scaling is undefined for an all-zero input).

- If every index is in `fixed_idx` (no free elements), the input is
  returned as-is — the caller is responsible for ensuring
  `sum(x) == target` when no elements may be modified.

## Examples

``` r

x <- c(0.2, 0.3, 0.500000001)
print(sum(x) == 1)
#> [1] FALSE

x_sum_to_1 <- normalise_sum(x)
print(all.equal(sum(x_sum_to_1), 1))
#> [1] TRUE

# With fixed elements: indices 1 and 3 are held constant
x2 <- c(0.25, 0.45, 0.30)
normalise_sum(x2, fixed_idx = c(1L, 3L))
#> [1] 0.25 0.45 0.30
```
