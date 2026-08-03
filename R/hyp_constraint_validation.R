# main assertion function -------------------------------------------------

#' Assert hypothesis weight vector constraint is valid
#'
#' `assert_hyp_constraint()` errors when the `hyp_constraint` component of the
#' input does not meet the validation criteria.
#'
#' Wraps individual assessment functions:
#'   * `assert_hyp_constraint_values()`: asserts individual values are between 0
#'   and 1.
#'   * `assert_hyp_constraint_na()`: asserts the vector does not have a single
#'   optimisable (i.e. `NA`) value.
#'   * `assert_hyp_constraint_sum()`: asserts the sum of the vector is not
#'   greater than 1. For _complete_ vectors it must be 1, while for
#'   _incomplete_ ones it should be less than 1.
#'
#' @inheritParams validate_graph_constraint
#'
#' @returns The input `x` invisibly.
#' @noRd
assert_hyp_constraint <- function(
    x,
    call = rlang::caller_env()
) {
    assert_hyp_constraint_values(x, call = call)
    assert_hyp_constraint_na(x, call = call)
    assert_hyp_constraint_sum(x, call = call)

    invisible(x)
}


## granular assertion functions -------------------------------------------

assert_hyp_constraint_sum <- function(
    x,
    call = rlang::caller_env()
) {
    hc <- x$hyp_constraint
    tolerance <- attr(x, "tolerance")

    sum_hc <- sum(hc, na.rm = TRUE)

    hc_complete <- all(vctrs::vec_detect_complete(hc))

    if (vctrs::vec_is_empty(hc)) {
        hc_complete <- FALSE
    }

    if (sum_hc - 1 > tolerance) {
        cli::cli_abort(
            c(
                "The sum of the hypothesis weight constraint vector cannot be \\
                greater than {.field 1}. It is {cli::col_red(sum_hc)}.",
                i = info_diagnose
            ),
            call = call
        )
    }

    # we cannot use `==` for equality comparisons due to small floating point
    # differences
    if (
        hc_complete && isTRUE(check_not_equal(sum_hc, 1, tolerance = tolerance))
    ) {
        cli::cli_abort(
            c(
                "The sum of a {.emph complete} hypothesis weight constraint \\
                vector must be {.field 1}. It is {cli::col_red(sum_hc)}.",
                i = info_diagnose
            ),
            call = call
        )
    }

    if (!hc_complete && isTRUE(check_equal(sum_hc, 1, tolerance = tolerance))) {
        cli::cli_abort(
            c(
                "The sum of an {.emph incomplete} hypothesis weight \\
                constraint vector must not be equal to {.field 1}.",
                i = info_diagnose
            ),
            call = call
        )
    }

    invisible(x)
}

assert_hyp_constraint_values <- function(
    x,
    call = rlang::caller_env()
) {
    hc <- x$hyp_constraint
    tolerance <- attr(x, "tolerance")

    if (any(hc - 1 > tolerance, na.rm = TRUE)) {
        cli::cli_abort(
            c(
                "Values in the hypothesis weight constraint vector cannot be \\
                greater than {.field 1}.",
                i = info_diagnose
            ),
            call = call
        )
    }

    if (any(hc < 0, na.rm = TRUE)) {
        cli::cli_abort(
            c(
                "Values in the hypothesis weights constraint vector cannot be \\
                less than {.field 0}.",
                i = info_diagnose
            ),
            call = call
        )
    }

    invisible(x)
}

assert_hyp_constraint_na <- function(x, call = rlang::caller_env()) {
    hc <- x$hyp_constraint

    number_nas <- sum(is.na(hc))

    if (number_nas == 1) {
        cli::cli_abort(
            c(
                "An {.emph incomplete} hypothesis weight constraint vector \\
                cannot have a single optimisable (i.e. {.code NA}) value.",
                i = info_diagnose
            ),
            call = call
        )
    }

    invisible(x)
}


# main diagnosis function -------------------------------------------------

#' Diagnose the hypothesis weight constraint vector
#'
#' `diagnose_hyp_constraint()` does not error, it simply reports on the state of
#' the `hyp_constraint` component of the input relative to the validation
#' criteria.
#'
#' Wraps individual diagnosis functions:
#'   * `diagnose_hyp_constraint_gt1`: checks individual values are not greater
#'   than 1.
#'   * `diagnose_hyp_constraint_lt0()`: checks individual values are not less
#'   than 0.
#'   * `diagnose_hyp_constraint_na()`: checks the vector does not have a single
#'   optimisable (i.e. `NA`) value.
#'   * `diagnose_hyp_constraint_sum()`: checks the sum of the vector is not
#'   greater than 1. For _complete_ vectors it must be 1, while for _incomplete_
#'   ones it should be less than 1.
#'
#' @inheritParams validate_graph_constraint
#'
#' @returns The input `x` invisibly.
#' @noRd
diagnose_hyp_constraint <- function(x) {
    cli::cli_h3("Hypothesis weight constraint diagnosis:")

    diagnose_hyp_constraint_gt1(x)
    diagnose_hyp_constraint_lt0(x)
    diagnose_hyp_constraint_na(x)
    diagnose_hyp_constraint_sum(x)

    invisible(x)
}


## granular diagnosis functions -------------------------------------------

diagnose_hyp_constraint_sum <- function(x) {
    hc <- x$hyp_constraint
    tolerance <- attr(x, "tolerance")

    sum_hc <- sum(hc, na.rm = TRUE)
    hc_complete <- all(vctrs::vec_detect_complete(hc))

    if (vctrs::vec_is_empty(hc)) {
        hc_complete <- FALSE
    }

    if (!hc_complete) {
        if (sum_hc - 1 > tolerance) {
            cli::cli_alert_danger(
                "The sum of the {.emph incomplete} hypothesis weight \\
                constraint vector cannot be greater than {.field 1}. It is \\
                {cli::col_red(sum_hc)}."
            )
        } else if (isTRUE(check_equal(sum_hc, 1, tolerance = tolerance))) {
            cli::cli_alert_danger(
                "The sum of the {.emph incomplete} hypothesis weight \\
                constraint vector is {.field 1}."
            )
        } else {
            cli::cli_alert_success(
                "The sum of the {.emph incomplete} hypothesis weights \\
                constraint vector is less than {.field 1}."
            )
        }
    }

    if (hc_complete) {
        if (isTRUE(check_equal(sum_hc, 1, tolerance = tolerance))) {
            cli::cli_alert_success(
                "The sum of the {.emph complete} hypothesis weight constraint \\
                vector is {.field 1}."
            )
        } else {
            cli::cli_alert_danger(
                "The sum of the {.emph complete} hypothesis weights \\
                constraint vector is {cli::col_red(sum_hc)}. It should be \\
                {.field 1}."
            )
        }
    }

    invisible(x)
}

diagnose_hyp_constraint_gt1 <- function(x) {
    hc <- x$hyp_constraint
    tolerance <- attr(x, "tolerance")

    if (any(hc - 1 > tolerance, na.rm = TRUE)) {
        cli::cli_alert_danger(
            "Some hypothesis weight constraint values are greater than \\
            {.field 1}."
        )
    } else {
        cli::cli_alert_success(
            "All hypothesis weights constraint values are less than or equal \\
            to {.field 1}."
        )
    }

    invisible(x)
}

diagnose_hyp_constraint_lt0 <- function(x) {
    hc <- x$hyp_constraint

    if (any(hc < 0, na.rm = TRUE)) {
        cli::cli_alert_danger(
            "Some hypothesis weight constraint values are less than {.field 0}."
        )
    } else {
        cli::cli_alert_success(
            "All hypothesis weights constraint values are greater than or \\
            equal to {.field 0}."
        )
    }

    invisible(x)
}

diagnose_hyp_constraint_na <- function(x) {
    hc <- x$hyp_constraint

    number_nas <- sum(is.na(hc))

    if (number_nas == 0) {
        cli::cli_alert_info(
            "The hypothesis weights constraint vector is {.emph complete} and \\
            no weights will be optimised."
        )
    } else if (number_nas == 1) {
        cli::cli_alert_danger(
            "A single hypothesis weight constraint cannot be optimised as \\
            weights must add up to {.field 1}."
        )
    } else {
        cli::cli_alert_success(
            "The hypothesis weights constraint vector is correctly defined \\
            for optimising {cli::col_blue(number_nas)} weights."
        )
    }

    invisible(x)
}
