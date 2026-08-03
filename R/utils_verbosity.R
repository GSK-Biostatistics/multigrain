verbosity_levels <- c("info", "detail", "silent")

#' Multigrain verbosity
#'
#' `multigrain_verbosity()` returns the option named `multigrain_verbosity`,
#' which controls multigrain's verbosity. There are 3 possible levels
#' ("detail" > "info" > "silent"):
#' * "info" (default): will only show milestones / informational messages
#' highlighting the progress of the optimisation at coarse-grained level.
#' * "detail":  will show milestones and information about fine-grained
#' optimisation events.
#' * "silent": no information about the progress of the optimisation is
#' printed to the console. Errors and warnings are still thrown normally.
#'
#' If the `multigrain_verbosity` option is unset, then the "info" level will
#' be used.
#'
#' @returns A string indicating the verbosity level as set with options or
#' "info" if `multigrain_verbosity` is unset.
#'
#' @export
#' @examples
#' multigrain_verbosity()
multigrain_verbosity <- function() {
    mv_opt <- getOption("multigrain_verbosity", "info")

    if (!rlang::is_string(mv_opt) || !(mv_opt %in% verbosity_levels)) {
        cli::cli_abort(
            'Option "multigrain_verbosity" must be one of: \\
            {.or {.field {verbosity_levels}}}.'
        )
    }
    mv_opt
}
