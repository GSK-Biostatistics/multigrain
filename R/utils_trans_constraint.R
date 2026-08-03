# diagnosis helper functions ----------------------------------------------

#' Print offending transition matrix row
#'
#' Helper function used for printing rows that fail certain checks. Used by the
#' `diagnose_trans_constraint_...()` family of functions.
#'
#' @param row (numeric) row to print
#' @param x (matrix) transition matrix
#' @param type (character) check type. One of `"sum"`, `"value"` or `"na".`
#' @param ref (numeric) reference value for comparison. Either `0` or `1`.
#' @inheritParams validate_graph_constraint
#'
#' @returns `NULL` invisibly as it is called for its side-effects.
#' @noRd
print_offending_row <- function(
    row,
    x,
    type = c("sum", "value", "na"),
    ref = c(0, 1),
    tolerance = sqrt(.Machine$double.eps)
) {
    type <- rlang::arg_match(type)
    offending_rows_values <- pull_rows(x)[[row]]

    values_string <- offending_rows_values |>
        str_trunc_light() |>
        # we need to use paste() since NAs are contagious in {glue}
        paste(collapse = ", ")

    values_string <- paste0("[", values_string, "]")

    if (type == "sum") {
        cli::cli_li(
            "Row {row}: {cli::col_blue(values_string)} has a sum of \\
            {cli::col_red(sum(offending_rows_values, na.rm = TRUE))}."
        )
    }

    if (type == "value") {
        if (ref == 1) {
            # nolint start: object_usage_linter
            # {lintr} false positive: it does not pick up the use of
            # `offending_values` in the cli_li call
            offending_values <- offending_rows_values[
                !is.na(offending_rows_values) &
                    offending_rows_values - 1 > tolerance
            ]

            cli::cli_li(
                "Row {row}: {cli::col_blue(values_string)} contains at least \\
                one value greater than {.field 1}: \\
                {cli::col_red(offending_values)}."
            )
        }

        if (ref == 0) {
            offending_values <- offending_rows_values[
                !is.na(offending_rows_values) &
                    offending_rows_values < 0
            ]
            # nolint end

            cli::cli_li(
                "Row {row}: {cli::col_blue(values_string)} contains at least \\
                one value less than {.field 0}: \\
                {cli::col_red(offending_values)}."
            )
        }
    }

    if (type == "na") {
        # nolint start: object_usage_linter
        # {lintr} false positive: it does not pick up the use of `row_sum` in
        # the cli_li() call
        row_sum <- sum(offending_rows_values, na.rm = TRUE)
        # nolint end

        cli::cli_li(
            "Row {row}: {cli::col_blue(values_string)} the optimisable value \\
            is effectively equal to {cli::col_red(1 - row_sum)}."
        )
    }

    invisible(NULL)
}

#' Decompose a matrix into its rows
#'
#' Used for pulling out rows that do not match expectations, in preparation for
#' printing.
#'
#' @param x (numeric matrix)
#'
#' @returns A list of vectors each representing a matrix row.
#' @noRd
pull_rows <- function(x) {
    #  x is a transition matrix
    purrr::map(
        seq_len(nrow(x)),
        pull_row,
        x
    )
}

#' Pull/ extract a matrix row
#'
#' @param row (numeric) row to extract
#' @param x (matrix) matrix to extract row from
#'
#' @returns A vector representing the values in the given matrix row.
#' @noRd
pull_row <- function(row, x) {
    x[row, ]
}

# we do not want to depend on stringr just for stringr::str_trunc. this is a
# lighter version doing only what we need with base R

#' Truncate a vector to a given number of characters
#'
#' Unlike `stringr::str_trunc()`, `str_trunc_light()` only truncates to the
#' right.
#'
#' @param x A vector, either character or something that can be cast to
#'   character.
#' @param width A scalar integerish the desired output width/length of the
#'   string.
#' @param ellipsis A character scalar to indicate content has been removed.
#'
#' @returns A character vector of truncated values.
#'
#' @noRd
#'
#' @examples
#' str_trunc_light(0.000000434)
str_trunc_light <- function(
    x,
    width = 7,
    ellipsis = "..."
) {
    withr::local_options(list(scipen = 999))
    x <- as.character(x)

    out <- x
    is_na <- is.na(x)
    num_char <- nchar(x)
    keep <- pmax(width - nchar(ellipsis), 0L)
    needs_trunc <- !is_na & num_char > width

    out[needs_trunc] <- paste0(substr(x[needs_trunc], 1L, keep), ellipsis)

    out[is_na] <- NA_character_
    out
}

offending_rows_bullets <- function(
    rows,
    x,
    type = c("sum", "value", "na"),
    ref = c(0, 1),
    tolerance = sqrt(.Machine$double.eps)
) {
    # cli_() functions return "semantic CLI elements"
    # cli::cat_() functions do not
    # semantic elements show up as "messages"
    # cli_fmt allows us to capture the semantic element and then call cat_line
    # on it, if needed. for example, we do not need this in the diagnosis
    # functions where we want messages

    # we are happy with `offending_rows_bullets()` outputting messages, so we
    # do not need cli_fmt to capture the div
    cli::cli_div(
        theme = list(
            ul = list(
                `margin-left` = 2,
                before = ""
            )
        )
    )
    ulid <- cli::cli_ul()
    purrr::walk(
        rows,
        print_offending_row,
        x = x,
        type = type,
        ref = ref,
        tolerance = tolerance
    )
    cli::cli_end(ulid)
}
