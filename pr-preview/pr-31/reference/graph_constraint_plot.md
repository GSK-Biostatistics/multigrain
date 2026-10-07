# Autoplot method for `multigrain_graph_constraint` objects

Autoplot method for `multigrain_graph_constraint` objects

## Usage

``` r
# S3 method for class 'multigrain_graph_constraint'
autoplot(object, ..., root = NULL, digits = NULL, title = NULL)

# S3 method for class 'multigrain_graph_constraint'
plot(x, ..., root = NULL, digits = NULL, title = NULL)
```

## Arguments

- object:

  A `multigrain_graph_constraint` object.

- ...:

  Additional arguments.

- root:

  A numeric vector indicating which nodes to be regarded as root in the
  tree layout. Passed down to
  [`igraph::layout_as_tree()`](https://r.igraph.org/reference/layout_as_tree.html).

- digits:

  Number of decimal places to round to (between 0 and 3). Defaults to
  `NULL` which will choose the number of digits based on how crowded the
  graph plot will be:

  - a single digit if more than 8 hypotheses

  - 2 digits if 5 to 8 hypotheses

  - 3 digits if 4 or fewer hypotheses

- title:

  An optional plot title.

- x:

  A `multigrain_graph_constraint` object.

## Value

A
[`ggraph::ggraph()`](https://ggraph.data-imaginist.com/reference/ggraph.html)
object.

## Details

Both [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
methods require an optional dependency,
[ggplot2](https://ggplot2.tidyverse.org/reference/ggplot2-package.html).

## Examples

``` r
library(ggplot2)
gc1 <- graph_constraint(
    c(NA, NA, 0, 0, 0),
    matrix(
       c(
           0, 0.8, 0.2, 0, 0,
           NA, 0, 0, NA, 0,
           NA, 0.8, 0, 0.1, NA,
           0.9, NA, NA, 0, NA,
           0, 0, 0, 1, 0
       ),
       nrow = 5,
       byrow = TRUE
   )
)

autoplot(gc1)


# If you want to control which nodes are treated as root, you can pass them
# as a numeric vector via the `root` argument. For example, we want nodes 1
# and 3 to be plotted as root:
autoplot(gc1, root = c(1, 3))


# You can control the number of decimal places to round to
autoplot(gc1, digits = 1)


# You can supply a title
autoplot(gc1, title = "My graph constraint")
```
