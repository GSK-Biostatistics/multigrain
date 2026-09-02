#' Validate a `multigrain_graph_constraint` object
#'
#' `validate_graph_constraint()` supports `graph_constraint()` to ensure its
#' output is always a valid object. It wraps the following functions:
#'   * diagnosis functions: optional, they inform both on success and failure:
#'      * `diagnose_graph_constraint_not_full()`: checks there are values to
#'      optimise (either in `hyp_constraint` or in `trans_constraint`).
#'      * `diagnose_hyp_constraint()`: runs:
#'        * `diagnose_hyp_constraint_gt1()`: informs on values greater than 1
#'        in the hypothesis weight constraint vector.
#'        * `diagnose_hyp_constraint_lt0()`: informs on values less than 1 in
#'        the hypothesis weight constraint vector.
#'        * `diagnose_hyp_constraint_na()`: informs on the number of optimisable
#'        values (i.e. `NA`s) in the hypothesis weight constraint vector.
#'        * `diagnose_hyp_constraint_sum()`: informs on the sum of the elements
#'        of the hypothesis weight constraint vector.
#'      * `diagnose_trans_constr()`: wraps:
#'        * `diagnose_trans_constr_diagonal()`: informs on the transition matrix
#'        constraint diagonal values.
#'        * `diagnose_trans_constr_square()`: informs on the squareness of the
#'        transition matrix constraint.
#'        * `diagnose_trans_constr_gt1()`: informs on values greater than 1.
#'        * `diagnose_trans_constr_lt0()`: informs on values less than 0.
#'        * `diagnose_trans_constr_row_na()`: informs on the number of `NA`s.
#'        * `diagnose_trans_constr_row_sums()`: informs on row sums.
#'      * `diagnose_hc_tc_consistency()`: informs on dimensional consistency
#'      between the hypothesis weight and the transition matrix constraints.
#'   * assertion functions: they error when condition is not met:
#'     * `assert_graph_constraint_not_full()`: asserts there are values to
#'     optimise (neither `hyp_constraint` or `trans_constraint` are "full").
#'     * `assert_hyp_constraint()`: wraps:
#'        * `assert_hyp_constraint_values()`: asserts values are between 0
#'        and 1.
#'        * `assert_hyp_constraint_na()`: asserts there isn't a single
#'        optimisable (i.e. `NA`) value.
#'        * `assert_hyp_constraint_sum()`: asserts the sum of the vector is not
#'        greater than 1.
#'     * `assert_trans_constr()`: wraps
#'       * `assert_trans_constr_diagonal()`: asserts the diagonal is 0.
#'       * `assert_trans_constr_square()`: asserts the transition matrix is
#'       square.
#'       * `assert_trans_constr_values()`: asserts the values are between 0 and
#'       1.
#'       * `assert_trans_constr_row_na()`: asserts the rows do not have a single
#'       optimisable (i.e. `NA`) value.
#'       * `assert_trans_constr_row_sums()`: asserts the sum of the rows is not
#'       greater than 1.
#'    * `assert_hc_tc_consistency()`: asserts the transition constraint is an
#'     `n * n` matrix where `n` is the number of hypotheses (the length of the
#'     hypothesis weight constraint vector)
#'
#' @param x (list) a list, containing 2 named elements `hyp_constraint` and
#'   `trans_constraint`. A precursor of a `multigrain_graph_constraint` object.
#'   Can also be a `multigrain_graph_constraint` object.
#' @inheritParams graph_constraint diagnose
#' @param call The execution environment of a currently running function. Since
#'   this is a helper for `graph_constraint()` we don't want the errors to
#'   surface from `validate_trans_constraint()`, but rather from
#'   `graph_constraint()`. The corresponding function call is retrieved and
#'   mentioned in the error messages as the source of the error.
#'
#' @returns the input `x`
#' @noRd
#'
#' @examples
#' gc <- graph_constraint_free(4)
#'
#' multigrain:::validate_graph_constraint(gc)
validate_graph_constraint <- function(
    x,
    diagnose = FALSE,
    call = rlang::caller_env()
) {
    if (diagnose) {
        diagnose_graph_constraint_not_full(x)
        diagnose_hyp_constraint(x)
        diagnose_trans_constr(x)
        diagnose_hc_tc_consistency(x)
    }

    assert_graph_constraint_not_full(x, call = call)
    assert_hyp_constraint(x, call = call)
    assert_trans_constr(x, call = call)
    assert_hc_tc_consistency(x, call = call)

    x
}

info_diagnose <- glue::glue(
    "For a more detailed diagnosis run `graph_constraint()` with \\
    `diagnose = TRUE`."
)

assert_hc_tc_consistency <- function(x, call = rlang::caller_env()) {
    hc <- x$hyp_constraint
    tc <- x$trans_constraint

    if (length(hc) != nrow(tc) || length(hc) != ncol(tc)) {
        cli::cli_abort(
            c(
                "The hypothesis weight and transition matrix constraints must \\
                be dimensionally consistent. The transition matrix constraint \\
                has {cli::col_red(nrow(tc))} rows and \\
                {cli::col_red(ncol(tc))} columns while the hypothesis weight \\
                vector contains {cli::col_green(length(hc))} elements.",
                i = info_diagnose
            ),
            call = call
        )
    }

    invisible(x)
}

diagnose_hc_tc_consistency <- function(x) {
    hc <- x$hyp_constraint
    tc <- x$trans_constraint
    m <- length(hc)

    cli::cli_h3(
        "Hypothesis weight and transition matrix constraints consistency:"
    )

    cli::cli_alert_success(
        "The hypothesis weights vector has {.field {m}} elements."
    )

    if (m != nrow(tc) || m != ncol(tc)) {
        cli::cli_alert_danger(
            "The transition matrix must have {.field {m}} rows and \\
            {.field {m}} columns. It has {cli::col_red(nrow(tc))} rows and \\
            {cli::col_red(ncol(tc))} columns. "
        )
    } else {
        cli::cli_alert_success(
            "The transition matrix has {.field {m}} columns and {.field {m}} \\
            rows."
        )
    }

    invisible(x)
}

# nolint start: object_length_linter
assert_graph_constraint_not_full <- function(x, call = rlang::caller_env()) {
    # nolint end
    nas_in_hyp_constraint <- sum(is.na(x$hyp_constraint))
    nas_in_trans_constraint <- sum(is.na(x$trans_constraint))

    if (nas_in_hyp_constraint + nas_in_trans_constraint == 0) {
        cli::cli_abort(
            c(
                "No hypothesis or transition weights to optimise.",
                i = info_diagnose
            ),
            call = call
        )
    }

    invisible(x)
}

# nolint start: object_length_linter
diagnose_graph_constraint_not_full <- function(x) {
    # nolint end
    nas_in_hyp_constraint <- sum(is.na(x$hyp_constraint))
    nas_in_trans_constraint <- sum(is.na(x$trans_constraint))

    if (nas_in_hyp_constraint + nas_in_trans_constraint == 0) {
        cli::cli_alert_danger(
            "The {.code multigrain_graph_constraint} object is \\
            {.emph complete}. There are no hypothesis or transition weights \\
            to optimise."
        )
    }

    invisible(x)
}
