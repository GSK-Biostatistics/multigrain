check_length <- function(
    x,
    m,
    ...,
    allow_null = FALSE,
    arg = rlang::caller_arg(x),
    call = rlang::caller_env()
) {
    if (!missing(x)) {
        if (allow_null && rlang::is_null(x)) {
            return(invisible(NULL))
        }
        if (length(x) == m) {
            return(invisible(NULL))
        }
    }

    cli::cli_abort(
        "{.arg {arg}} has length {length(x)}; expected {m}.",
        call = call
    )
}

check_dim <- function(
    x,
    m,
    ...,
    allow_null = FALSE,
    arg = rlang::caller_arg(x),
    call = rlang::caller_env()
) {
    if (!missing(x)) {
        if (allow_null && rlang::is_null(x)) {
            return(invisible(NULL))
        }
        if (all(dim(x) == c(m, m))) {
            return(invisible(NULL))
        }
    }

    cli::cli_abort(
        "{.arg {arg}} has dim ({readable_dim(x)}); expected ({m} x {m}).",
        call = call
    )
}

readable_dim <- function(x) {
    glue::glue("{dim(x)[1]} x {dim(x)[2]}")
}
