# main assertion function -------------------------------------------------

#' Assert transition matrix constraint is valid
#'
#' `assert_trans_constr()` errors when the `trans_constraint` component of
#' the input does not meet the validation criteria.
#'
#' Wraps several granular assertion functions:
#'   * `assert_trans_constr_diagonal()`: asserts the diagonal is 0.
#'   * `assert_trans_constr_square()`: asserts the matrix is square.
#'   * `assert_trans_constr_values()`: asserts values are between 0 and 1.
#'   * `assert_trans_constr_row_na()`: asserts rows have 0 or more than 2
#'   `NA`s.
#'   * `assert_trans_constr_row_sums()`: asserts the row sums are less than
#'   or equal to 1.
#'
#' @inheritParams validate_graph_constraint
#'
#' @returns `x` invisibly.
#' @noRd
assert_trans_constr <- function(
    x,
    call = rlang::caller_env()
) {
    assert_trans_constr_diagonal(x, call = call)
    assert_trans_constr_square(x, call = call)
    assert_trans_constr_values(x, call = call)
    assert_trans_constr_row_na(x, call = call)
    assert_trans_constr_row_sums(x, call = call)

    invisible(x)
}


## granular assertion functions -------------------------------------------

assert_trans_constr_diagonal <- function(
    x,
    call = rlang::caller_env()
) {
    tc <- x$trans_constraint
    diagonal <- diag(tc)

    if (anyNA(diagonal)) {
        cli::cli_abort(
            c(
                "There must not be any {.code NA}s on the transition matrix \\
                diagonal. All values must be {.field 0}.",
                i = info_diagnose
            ),
            call = call
        )
    }

    if (any(diagonal != 0)) {
        cli::cli_abort(
            c(
                "All values on the transition matrix diagonal must be \\
                {.field 0}.",
                i = info_diagnose
            ),
            call = call
        )
    }

    invisible(x)
}

assert_trans_constr_square <- function(
    x,
    call = rlang::caller_env()
) {
    tc <- x$trans_constraint

    if (dim(tc)[1] != dim(tc)[2]) {
        cli::cli_abort(
            c(
                "The transition matrix is not square.",
                i = info_diagnose
            ),
            call = call
        )
    }

    invisible(x)
}

assert_trans_constr_values <- function(
    x,
    call = rlang::caller_env()
) {
    tc <- x$trans_constraint
    tolerance <- attr(x, "tolerance")

    if (any(tc - 1 > tolerance, na.rm = TRUE)) {
        cli::cli_abort(
            c(
                "Some transition matrix values are greater than {.field 1}.",
                i = info_diagnose
            ),
            call = call
        )
    }

    if (any(tc < 0, na.rm = TRUE)) {
        cli::cli_abort(
            c(
                "Some transition matrix values are less than {.field 0}.",
                i = info_diagnose
            ),
            call = call
        )
    }

    invisible(x)
}

assert_trans_constr_row_na <- function(x, call = rlang::caller_env()) {
    tc <- x$trans_constraint

    na_per_row <- rowSums(is.na(tc))

    if (any(na_per_row == 1)) {
        cli::cli_abort(
            c(
                "At least one {.emph incomplete} transition matrix row has a \\
                single optimisable (i.e. {.code NA}) value.",
                i = info_diagnose
            ),
            call = call
        )
    }

    invisible(x)
}

assert_trans_constr_row_sums <- function(
    x,
    call = rlang::caller_env()
) {
    tc <- x$trans_constraint
    tolerance <- attr(x, "tolerance")

    trans_constraint_row_complete <- vctrs::vec_detect_complete(tc)
    complete_rows <- which(trans_constraint_row_complete)
    incomplete_rows <- which(!trans_constraint_row_complete)

    row_sums <- rowSums(tc, na.rm = TRUE)
    sums_complete_rows <- row_sums[complete_rows]
    sums_incomplete_rows <- row_sums[incomplete_rows]

    if (any(row_sums - 1 > tolerance)) {
        cli::cli_abort(
            c(
                "The sum of transition matrix constraint rows cannot be \\
                greater than {.field 1}.",
                i = info_diagnose
            ),
            call = call
        )
    }

    if (!isTRUE(check_equal(sums_complete_rows, 1, tolerance))) {
        cli::cli_abort(
            c(
                "The sum of {.emph complete} transition matrix constraint \\
                rows must be {.field 1}.",
                i = info_diagnose
            ),
            call = call
        )
    }

    if (any(vec_check_equal(sums_incomplete_rows, 1, tolerance))) {
        cli::cli_abort(
            c(
                "At least one {.emph incomplete} transition matrix row has a \\
                sum equal to {.field 1}.",
                i = info_diagnose
            ),
            call = call
        )
    }

    invisible(x)
}


# main diagnosis function -------------------------------------------------

#' Diagnose the transition matrix constraint
#'
#' `diagnose_trans_constr()` does not error, it simply reports on the state
#' of the transition matrix constraint (`trans_constraint`) component of the
#' input relative to the validation criteria.
#'
#' Wraps individual diagnosis functions:
#'   * `diagnose_trans_constr_diagonal()`: informs if the diagonal is 0.
#'   * `diagnose_trans_constr_square()`: informs on the squareness.
#'   * `diagnose_trans_constr_gt1()`: informs on values compared to 1.
#'   * `diagnose_trans_constr_lt0()`: informs on values compared to 0.
#'   * `diagnose_trans_constr_row_na()`: informs on number of `NA`s.
#'   * `diagnose_trans_constr_row_sums()`: informs on row sums.
#'
#' @inheritParams validate_graph_constraint
#'
#' @returns `x` invisibly.
#' @noRd
diagnose_trans_constr <- function(x) {
    cli::cli_h3("Transition matrix constraint diagnosis:")

    diagnose_trans_constr_diagonal(x)
    diagnose_trans_constr_square(x)
    diagnose_trans_constr_gt1(x)
    diagnose_trans_constr_lt0(x)
    diagnose_trans_constr_row_na(x)
    diagnose_trans_constr_row_sums(x)

    invisible(x)
}


## granular diagnosis functions -------------------------------------------

diagnose_trans_constr_diagonal <- function(x) {
    tc <- x$trans_constraint

    diagonal <- diag(tc)

    if (anyNA(diagonal)) {
        cli::cli_alert_danger(
            "All values on the transition matrix diagonal must be {.field 0}. \\
            Some values are {.code NA}."
        )
    } else if (any(diagonal != 0)) {
        cli::cli_alert_danger(
            "Some values on the transition matrix diagonal are not {.field 0}."
        )
    } else {
        cli::cli_alert_success(
            "All values on the transition matrix diagonal are {.field 0}."
        )
    }

    invisible(x)
}

diagnose_trans_constr_square <- function(x) {
    tc <- x$trans_constraint

    if (dim(tc)[1] != dim(tc)[2]) {
        cli::cli_alert_danger(
            "The transition matrix is not square."
        )
    } else {
        cli::cli_alert_success(
            "The transition matrix is square."
        )
    }
    invisible(x)
}

diagnose_trans_constr_row_sums <- function(x) {
    tc <- x$trans_constraint
    tolerance <- attr(x, "tolerance")

    tc_row_complete <- vctrs::vec_detect_complete(tc)
    complete_rows <- which(tc_row_complete)
    incomplete_rows <- which(!tc_row_complete)

    row_sums_without_na <- rowSums(tc, na.rm = TRUE)
    sums_complete_rows <- row_sums_without_na[complete_rows]
    sums_incomplete_rows <- row_sums_without_na[incomplete_rows]

    # check complete rows
    if (!rlang::is_empty(complete_rows)) {
        # check the sum of the complete rows is equal to 1
        if (isTRUE(check_equal(sums_complete_rows, 1, tolerance))) {
            cli::cli_alert_success(
                "All {.emph complete} transition matrix rows sum up to \\
                {.field 1}."
            )
        } else {
            cli::cli_alert_danger(
                "At least one {.emph complete} transition matrix row does not \\
                sum up to {.field 1}."
            )

            # we can't use sums_complete_rows != 1 as we need to account for the
            # tolerance
            offending_complete_rows <- sums_complete_rows |>
                vec_check_not_equal(1, tolerance) |>
                which()
            offending_rows <- complete_rows[offending_complete_rows]

            offending_rows_bullets(
                offending_rows,
                x = tc
            )
        }
    }

    # check incomplete rows
    if (!rlang::is_empty(incomplete_rows)) {
        # check the sum of the incomplete rows is not greater than 1
        if (any(sums_incomplete_rows - 1 > tolerance)) {
            cli::cli_alert_danger(
                "At least one {.emph incomplete} transition matrix row has a \\
                sum greater than {.field 1}."
            )
            offending_incomplete_rows <- which(
                sums_incomplete_rows - 1 > tolerance
            )
            offending_rows <- incomplete_rows[offending_incomplete_rows]

            offending_rows_bullets(
                offending_rows,
                x = tc
            )
        } else {
            cli::cli_alert_success(
                "No {.emph incomplete} transition matrix rows have a sum \\
                greater than {.field 1}."
            )
        }

        # check the sum of the incomplete rows is not equal to 1
        if (any(vec_check_equal(sums_incomplete_rows, 1, tolerance))) {
            cli::cli_alert_danger(
                "At least one {.emph incomplete} transition matrix row has \\
                a sum equal to {.field 1}."
            )

            offending_complete_rows <- sums_incomplete_rows |>
                vec_check_equal(1, tolerance) |>
                which()
            offending_rows <- incomplete_rows[offending_complete_rows]

            offending_rows_bullets(
                offending_rows,
                x = tc
            )
        } else {
            cli::cli_alert_success(
                "No {.emph incomplete} transition matrix rows have a sum \\
                equal to {.field 1}."
            )
        }
    }

    invisible(x)
}

diagnose_trans_constr_gt1 <- function(x) {
    tc <- x$trans_constraint
    tolerance <- attr(x, "tolerance")

    if (any(tc - 1 > tolerance, na.rm = TRUE)) {
        cli::cli_alert_danger(
            "Some transition matrix values are greater than {.field 1}."
        )
        # find the offending rows
        offending_rows <- purrr::map_lgl(
            pull_rows(tc),
            ~ any(.x - 1 > tolerance, na.rm = TRUE)
        )

        offending_rows <- which(offending_rows)

        offending_rows_bullets(
            offending_rows,
            x = tc,
            type = "value",
            ref = 1,
            tolerance = tolerance
        )
    } else {
        cli::cli_alert_success(
            "All transition matrix values are less than or equal to {.field 1}."
        )
    }

    invisible(x)
}

diagnose_trans_constr_lt0 <- function(x) {
    tc <- x$trans_constraint

    if (any(tc < 0, na.rm = TRUE)) {
        cli::cli_alert_danger(
            "Some transition matrix values are less than {.field 0}."
        )
        # find the offending rows
        offending_rows <- purrr::map_lgl(
            pull_rows(tc),
            ~ any(.x < 0, na.rm = TRUE)
        )
        offending_rows <- which(offending_rows)

        offending_rows_bullets(
            offending_rows,
            x = tc,
            type = "value",
            ref = 0
        )
    } else {
        cli::cli_alert_success(
            "All transition matrix values are greater than or equal to \\
            {.field 0}."
        )
    }

    invisible(x)
}

diagnose_trans_constr_row_na <- function(x) {
    tc <- x$trans_constraint

    na_per_row <- rowSums(is.na(tc))

    if (any(na_per_row == 1)) {
        cli::cli_alert_danger(
            "At least one {.emph incomplete} transition matrix row has a \\
            single optimisable (i.e. {.code NA}) value."
        )
        offending_rows <- which(na_per_row == 1)

        offending_rows_bullets(
            offending_rows,
            x = tc,
            type = "na"
        )
    } else {
        cli::cli_alert_success(
            "All transition matrix rows are either complete or have at least \\
            {.field 2} values to optimise."
        )
    }

    invisible(x)
}
