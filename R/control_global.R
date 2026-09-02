#' Modify global optimisation options
#'
#' `control_global()` is for expert use only; it allows you to directly set
#' [GA::ga()] options to access features that are otherwise not available in
#' multigrain.
#'
#' @inheritParams control_local
#' @param ... <[`dynamic-dots`][rlang::dyn-dots]> Name-value pairs of [GA::ga()]
#'   options and their values.
#'
#' @returns a modified [multigrain_control].
#'
#' @export
#' @examples
#' # `control_global()` allows you to access `GA::ga()` options that are not
#' # otherwise exposed by multigrain. For example, in special cases you may
#' # need to control the population size (the number of candidate graphs
#' # evaluated at each generation).
#' # multigrain makes some informed choices, but if you're convinced you want to
#' # try other values, you can access this (and other) [GA::ga()] options:
#' multigrain_control() |>
#'     control_global(pcrossover = 0.2)
control_global <- function(.ctrl, ...) {
    check_control(.ctrl)

    .ctrl$global_opt <- modify_list(.ctrl$global_opt, ...)
    .ctrl
}
