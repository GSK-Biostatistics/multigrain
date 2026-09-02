#' Set parameters for graph optimisation
#'
#' @description
#' There are two steps needed to fine tune the graph optimisation with
#' multigrain:
#'
#' 1. Create a control object with `multigrain_control()`.
#' 2. Define its behaviour with `control_` functions:
#'    * [control_nsim_local()] to set the number of local simulations.
#'    * [control_nsim_global()] to set the number of global simulations.
#'    * [control_local()] to modify the local optimisation options.
#'    * [control_global()] to modify the global optimisation options.
#'
#' Any unset parameters will automatically be set before running the
#' optimisation. There are predefined defaults, but they are calibrated based on
#' the `pvals` dimensions.
#'
#' @returns an S3 object (a list) with class `multigrain_control`.
#'
#' @export
#' @examples
#' multigrain_control()
multigrain_control <- function() {
    new_multigrain_control()
}


# low-level constructor --------------------------------------------------
new_multigrain_control <- function(
    nsim_local = NULL,
    nsim_global = NULL,
    global_opt = list(),
    local_opt = list()
) {
    structure(
        list(
            nsim_local = nsim_local,
            nsim_global = nsim_global,
            local_opt = local_opt,
            global_opt = global_opt
        ),
        class = "multigrain_control"
    )
}


# helpers ----------------------------------------------------------------
is_control <- function(x) {
    inherits(x, "multigrain_control")
}

check_control <- function(
    ctrl,
    arg = rlang::caller_arg(ctrl),
    call = rlang::caller_env(),
    allow_null = FALSE
) {
    if (!missing(ctrl)) {
        if (is_control(ctrl)) {
            return(invisible(NULL))
        }

        if (allow_null && is.null(ctrl)) {
            return(invisible(NULL))
        }
    }

    rlang::stop_input_type(
        ctrl,
        "a multigrain control object",
        allow_null = allow_null,
        arg = arg,
        call = call
    )
}


# print method -----------------------------------------------------------
#' @export
print.multigrain_control <- function(x, ...) {
    if (is.null(x)) {
        return()
    }

    cli::cat_line(cli::format_inline("{.cls {class(x)}}"))

    cat_control_opt("local simulations", x$nsim_local)
    cat_control_opt("global simulations", x$nsim_global)

    bullets_with_header("local optimisation:", x$local_opt)
    global_opt_sans_optim_args <- x$global_opt[
        names(x$global_opt) != "optimArgs"
    ]
    bullets_with_header("global optimisation:", global_opt_sans_optim_args)

    global_opt_optim_args <- x$global_opt[names(x$global_opt) == "optimArgs"]
    cat_optim_args(global_opt_optim_args)

    invisible(x)
}
