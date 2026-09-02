#' Enable or disable global search
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' Deprecated as we realised the user experience is much better with
#' `global_search` as a direct argument to [graph_optimise()].
#'
#' @inheritParams control_nsim_local
#' @param global_search (boolean) If `TRUE`, performs a global optimisation
#'   before the local optimisation. If unset, `TRUE` will be used.
#'
#' @returns a modified [multigrain_control].
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' # before
#' ctrl <- multigrain_control() |>
#'     control_global_search(FALSE)
#'
#' graph_optimise(..., control = ctrl)
#'
#' # now you pass global_search directly to graph_optimise()
#' graph_optimise(..., global_search = FALSE)
#' }
#'
control_global_search <- function(ctrl, global_search) {
    lifecycle::deprecate_stop(
        when = "0.3.0",
        what = I("`control_global_search()`"),
        with = I("the `global_search` `graph_optimise()` argument")
    )
}
