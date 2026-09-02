#' Set number of local simulations
#'
#' @param ctrl a [multigrain_control] object.
#' @param nsim_local The number of simulations to use when evaluating the trial
#'   success function in local optimisation.
#'   - If set and lower than the number of rows in the `pvals` matrix,
#'   `graph_optimise()` will use a random sample of `nsim_local` rows from
#'   `pvals` for local optimisation.
#'   - If unset or greater than or equal to the number of rows in the `pvals`
#'   matrix, _all_ rows from `pvals` will be used.
#'
#' @returns a modified [multigrain_control].
#'
#' @export
#' @examples
#' multigrain_control() |>
#'     control_nsim_local(10000)
control_nsim_local <- function(ctrl, nsim_local) {
    check_control(ctrl)
    rlang::check_number_whole(
        nsim_local,
        min = 1
    )

    ctrl$nsim_local <- nsim_local

    ctrl
}

adjust_nsim_local <- function(ctrl, nrow_pvals, call = rlang::caller_env()) {
    if (is.null(ctrl$nsim_local)) {
        return(ctrl)
    }

    if (ctrl$nsim_local > nrow_pvals) {
        msg <- glue::glue(
            "Number of simulations for trial success calculation \\
            (`nsim_local`) is greater than the number of rows in `pvals`."
        )
        bullet <- "Setting `nsim_local` = `nrows(pvals)`."
        cli::cli_warn(
            c(msg, "*" = bullet),
            call = call
        )

        ctrl$nsim_local <- nrow_pvals
    }

    ctrl
}

#' Set number of global simulations
#'
#' @inheritParams control_nsim_local
#' @param nsim_global The number of simulations to use when evaluating the trial
#'   success function in global optimisation. If unset, the minimum between
#'   50000 and the number of sets of p-values will be used.
#'
#' @returns A modified [multigrain_control].
#'
#' @export
#' @examples
#' multigrain_control() |>
#'     control_nsim_global(10000)
control_nsim_global <- function(ctrl, nsim_global) {
    check_control(ctrl)
    rlang::check_number_whole(
        nsim_global,
        min = 1
    )

    ctrl$nsim_global <- nsim_global

    ctrl
}

adjust_nsim_global <- function(ctrl, nrow_pvals, call = rlang::caller_env()) {
    if (is.null(ctrl$nsim_global)) {
        return(ctrl)
    }

    if (ctrl$nsim_global > nrow_pvals) {
        msg <- glue::glue(
            "Number of simulations for trial success calculation within \\
            global optimisation algorithm (`nsim_global`) is greater than \\
            the number of rows in `pvals`."
        )
        bullet <- "Setting `nsim_global` = `nrows(pvals)`."

        rlang::warn(
            c(msg, "*" = bullet),
            use_cli_format = TRUE,
            call = call
        )

        ctrl$nsim_global <- nrow_pvals
    }

    ctrl
}
