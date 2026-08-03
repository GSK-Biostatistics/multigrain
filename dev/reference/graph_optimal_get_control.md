# Get the optimisation control

Retrieve the settings used for the graph optimisation.

## Usage

``` r
graph_optimal_get_control(graph_optimal)
```

## Arguments

- graph_optimal:

  A `multigrain_graph_optimal` object.

## Value

The `multigrain_control` object holding the optimisation settings.

## Examples

``` r
# graph_optimal_example is an example optimised graph
graph_optimal_get_control(graph_optimal_example)
#> <multigrain_control>
#> local simulations: 20000
#> global simulations: 50000
#> local optimisation:
#> • algorithm: "NLOPT_LN_COBYLA"
#> • xtol_rel: 5e-08
#> • xtol_abs: 5e-09
#> • maxeval: 5000
#> • print_level: 1
#> global optimisation:
#> • pcrossover: 0.2
#> • pmutation: 0.8
#> • maxiter: 1e+05
#> • popSize: 200
#> • run: 7
#> • monitor: TRUE
#> • optimArgs: 
#>     • method: "Nelder-Mead"
#>     • poptim: 0.2
#>     • pressel: 0.6
```
