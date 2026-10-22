# Autoplot method for `multigrain_graph_optimal` objects

Autoplot method for `multigrain_graph_optimal` objects

## Usage

``` r
# S3 method for class 'multigrain_graph_optimal'
autoplot(object, ..., root = NULL, digits = NULL, title = NULL)

# S3 method for class 'multigrain_graph_optimal'
plot(x, ..., root = NULL, digits = NULL, title = NULL)
```

## Arguments

- object:

  A `multigrain_graph_optimal` object.

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

  A `multigrain_graph_optimal` object.

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

# plotting the example graph optimal
autoplot(graph_optimal_example)


# If you want to control which nodes are treated as root, you can pass them
# as a numeric vector via the `root` argument. For example, we do not like
# the default layout and we would like nodes 1 and 2 to be plotted as root:
autoplot(graph_optimal_example, root = c(1, 2))


# control the rounding with the `digits` argument.
autoplot(graph_optimal_example, root = c(1, 2), digits = 2)


# title
autoplot(graph_optimal_example, root = c(1, 2), title = "My graph optimal")
```
