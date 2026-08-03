is_nloptr <- function(x) {
    inherits(x, "nloptr")
}

check_nloptr <- function(
    nloptr,
    arg = rlang::caller_arg(nloptr),
    call = rlang::caller_env(),
    allow_null = FALSE
) {
    if (!missing(nloptr)) {
        if (is_nloptr(nloptr)) {
            return(invisible(NULL))
        }

        if (allow_null && is.null(nloptr)) {
            return(invisible(NULL))
        }
    }

    rlang::stop_input_type(
        nloptr,
        "a nloptr object",
        allow_null = allow_null,
        arg = arg,
        call = call
    )
}
