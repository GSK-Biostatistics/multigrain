# moved the code for `ga_N_options()` and `ga_gsd_options()`

#' @rdname ga_options
#' @export
ga_N_options <- function(
    ...,
    p_crossover = 0.2,
    p_mutation = 0.8,
    max_iter = 1e5,
    pop_size = 200,
    max_stagnation = 200,
    monitor = FALSE
) {
    output <- ga_options(
        ...,
        p_crossover = p_crossover,
        p_mutation = p_mutation,
        max_iter = max_iter,
        pop_size = pop_size,
        max_stagnation = max_stagnation,
        monitor = monitor
    )

    attr(output, "optimisation_type") <- "N"

    output
}

#' @rdname ga_options
#' @export
ga_gsd_options <- function(
    ...,
    p_crossover = 0.2,
    p_mutation = 0.8,
    max_iter = 1e5,
    pop_size = 200,
    max_stagnation = 200,
    monitor = FALSE
) {
    output <- ga_options(
        ...,
        p_crossover = p_crossover,
        p_mutation = p_mutation,
        max_iter = max_iter,
        pop_size = pop_size,
        max_stagnation = max_stagnation,
        monitor = monitor
    )

    attr(output, "optimisation_type") <- "gsd"

    output
}
