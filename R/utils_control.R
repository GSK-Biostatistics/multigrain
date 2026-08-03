cat_control_opt <- function(param, value) {
    if (length(value) == 0) {
        return()
    }

    cli::cat_line(
        cli::format_inline(
            "{.strong {param}:} {value}"
        )
    )
}

cat_optim_args <- function(optim_args) {
    if (length(optim_args) == 0) {
        return()
    }

    glue::glue("{cli::col_blue('optimArgs')}: ") |>
        cli::cat_bullet()

    flat_optim_args <- purrr::flatten(optim_args)

    nms <- names(flat_optim_args)
    vals <- purrr::map_chr(flat_optim_args, as_simple)

    # cli_() functions return "semantic CLI elements"
    # cli::cat_() functions do not
    # semantic elements show up as "messages"
    # cli_fmt allows us to capture the semantic element and then call cat_line
    # on it, if needed. for example, we do not need this in the diagnosis
    # functions where we want messages, but we need it here as we want "output"
    # also see https://github.com/r-lib/cli/issues/703
    output <- cli::cli_fmt({
        cli::cli_div(
            theme = list(
                ul = list(
                    `margin-left` = 4,
                    before = ""
                )
            )
        )

        ulid <- cli::cli_ul()
        glue::glue("{cli::col_green(nms)}: {vals}") |>
            cli::cli_li()
        cli::cli_end(ulid)
    })

    cli::cat_line(output)
}
