#' Modify local optimisation options
#'
#' `control_local()` is for expert use only; it allows you to directly set
#' [nloptr::nloptr()] options to access features that are otherwise not
#' available in multigrain.
#'
#' @param .ctrl A [multigrain_control] object.
#' @param ... <[`dynamic-dots`][rlang::dyn-dots]> Name-value pairs of
#'   [nloptr::nloptr()] options and their values.
#'
#' @returns A modified [multigrain_control].
#'
#' @export
#' @examples
#' # `control_local()` allows you to access `nloptr::nloptr()` options that are
#' # not otherwise exposed by multigrain. For example, in special cases you
#' # may want to try a different local optimisation algorithm. multigrain uses
#' # "NLOPT_LN_COBYLA", but if you're convinced you want to try other
#' # algorithms, you can access this (and other) `nloptr::nloptr()` options:
#' multigrain_control() |>
#'     control_local(algorithm = "NLOPT_LN_NEWUOA")
control_local <- function(.ctrl, ...) {
    check_control(.ctrl)

    .ctrl$local_opt <- modify_list(.ctrl$local_opt, ...)
    .ctrl
}
