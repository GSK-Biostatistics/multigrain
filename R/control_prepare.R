#' Create the default control object
#'
#' This will be used to set any values that were not modified by the user.
#'
#' @returns a [multigrain_control] object with default values
#'
#' @noRd
default_control <- function() {
    # nsim placeholders
    nsim_local <- integer()
    nsim_global <- integer()

    local_opt <- list(
        algorithm = "NLOPT_LN_COBYLA",
        xtol_rel = 5e-8,
        xtol_abs = 5e-9,
        maxeval = 5000,
        print_level = 0L
    )

    global_opt <- list(
        pcrossover = 0.2,
        pmutation = 0.8,
        maxiter = 1e5,
        popSize = 200,
        run = 200,
        monitor = FALSE,
        optimArgs = list(
            method = "Nelder-Mead",
            poptim = 0.2,
            pressel = 0.6
        )
    )

    new_multigrain_control(
        nsim_local = nsim_local,
        nsim_global = nsim_global,
        local_opt = local_opt,
        global_opt = global_opt
    )
}

#' Prepare the control object for optimisation
#'
#' Once the user has finished modifying the control object `control_prepare()`
#' injects the unspecified defaults. `control_prepare()` also calibrates some
#' of the arguments (context-aware defaults).
#'
#' @inheritParams control_nsim_local
#' @param pvals (optional) matrix of p values to use for calibration. We are
#'   interested in the number of rows and columns.
#' @inheritParams graph_optimise
#' @inheritParams rlang::args_error_context
#'
#' @returns a modified [multigrain_control]
#'
#' @noRd
#' @examples
#' # inject the defaults and calibrate the defaults based on pvals
#' pvals <- simulate_pvalues(
#'     power_nominal = c(0.9, 0.85, 0.8, 0.75),
#'     alpha = 0.025,
#'     corr_matrix = diag(4),
#'     nsim = 1e5
#' )
#' multigrain_control() |> control_prepare(pvals)
control_prepare <- function(
    ctrl,
    pvals,
    trace = FALSE,
    call = rlang::caller_env()
) {
    check_control(ctrl, call = call)
    check_double_matrix(pvals, call = call)
    check_logical(trace, allow_na = FALSE, call = call)

    default_ctrl <- default_control()

    # calibrate the default number of simulations with pvals
    default_ctrl$nsim_local <- nrow(pvals)
    default_ctrl$nsim_global <- min(5e4L, nrow(pvals))

    default_ctrl$global_opt$popSize <- min(max(40L * ncol(pvals), 200L), 500L)

    # calibrate user-supplied nsim values
    ctrl <- adjust_nsim_local(ctrl, nrow(pvals), call = call)
    ctrl <- adjust_nsim_global(ctrl, nrow(pvals), call = call)

    if (trace) {
        # adjusting the default verbosity allows the user-set values to pass
        default_ctrl$local_opt$print_level <- 1L
        default_ctrl$global_opt$monitor <- TRUE
    }

    ctrl$nsim_local <- ctrl$nsim_local %||% default_ctrl$nsim_local
    ctrl$nsim_global <- ctrl$nsim_global %||% default_ctrl$nsim_global

    ctrl$global_opt <- purrr::list_modify(
        default_ctrl$global_opt,
        !!!ctrl$global_opt
    )

    ctrl$local_opt <- purrr::list_modify(
        default_ctrl$local_opt,
        !!!ctrl$local_opt
    )

    ctrl
}
