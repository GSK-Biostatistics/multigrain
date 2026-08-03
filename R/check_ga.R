is_ga <- function(x) {
    inherits(x, "ga")
}


check_ga <- function(
    ga,
    arg = rlang::caller_arg(ga),
    call = rlang::caller_env(),
    allow_null = FALSE
) {
    if (!missing(ga)) {
        if (is_ga(ga)) {
            return(invisible(NULL))
        }

        if (allow_null && is.null(ga)) {
            return(invisible(NULL))
        }
    }

    rlang::stop_input_type(
        ga,
        "a GA object",
        allow_null = allow_null,
        arg = arg,
        call = call
    )
}
