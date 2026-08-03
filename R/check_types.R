# these `check_()` functions complement the rlang standalone checks:
#   * `check_double()`: checks an input is a double *vector*
#     * `rlang::check_number_decimal()` checks for scalars, not vectors
#     * can handle vectors of unoptimised weights (all NAs)
#   * `check_integerish()`: checks an input is an integer-like *vector*
#   * `check_double_matrix()`: checks an input is a double matrix

check_double <- function(
    x,
    ...,
    allow_na = FALSE,
    allow_null = FALSE,
    arg = rlang::caller_arg(x),
    call = rlang::caller_env()
) {
    if (!missing(x)) {
        if (rlang::is_double(x)) {
            return(invisible(NULL))
        }
        if (allow_null && rlang::is_null(x)) {
            return(invisible(NULL))
        }
        # allow a vector of all NAs to pass (as it can be cast to double)
        # this is mostly useful for hypothesis constraints (unoptimised weights)
        if (!rlang::is_null(x) && allow_na && all(is.na(x))) {
            return(invisible(NULL))
        }
    }

    rlang::stop_input_type(
        x,
        "a double",
        ...,
        allow_na = allow_na,
        allow_null = allow_null,
        arg = arg,
        call = call
    )
}

check_double_matrix <- function(
    x,
    ...,
    allow_null = FALSE,
    arg = rlang::caller_arg(x),
    call = rlang::caller_env()
) {
    if (!missing(x)) {
        if (rlang::is_double(x) && is.matrix(x)) {
            return(invisible(NULL))
        }
        if (allow_null && rlang::is_null(x)) {
            return(invisible(NULL))
        }
    }

    rlang::stop_input_type(
        x,
        "a double matrix",
        ...,
        allow_na = FALSE,
        allow_null = allow_null,
        arg = arg,
        call = call
    )
}

check_integerish <- function(
    x,
    ...,
    allow_na = FALSE,
    allow_null = FALSE,
    arg = rlang::caller_arg(x),
    call = rlang::caller_env()
) {
    if (!missing(x)) {
        if (rlang::is_integerish(x)) {
            return(invisible(NULL))
        }
        if (allow_null && rlang::is_null(x)) {
            return(invisible(NULL))
        }
    }

    rlang::stop_input_type(
        x,
        "integer-like",
        ...,
        allow_na = allow_na,
        allow_null = allow_null,
        arg = arg,
        call = call
    )
}
