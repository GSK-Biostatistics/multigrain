# `multigrain_family_spec` class -----------------------------------------

# Within-family tests. The position in this vector (0-based) is the method
# code used by the C++ routine `adjust_pvals_families()`.
family_tests <- c("hommel", "hochberg", "holm")

## constructor ------------------------------------------------------------

new_family_spec <- function(members = list(), test = character()) {
    stopifnot(
        is.list(members),
        is.character(test),
        length(test) == length(members)
    )

    structure(
        list(
            members = members,
            test = test
        ),
        class = "multigrain_family_spec",
        m = sum(lengths(members))
    )
}

## user-facing constructor ------------------------------------------------

#' Create a _family specification_ for graphs of families of hypotheses
#'
#' A _family specification_ groups hypotheses into families. When it is passed
#' to [graph_optimise()] or [calc_power_pvals()], each node of the graph is a
#' family rather than a single hypothesis, and the graph is tested with the
#' procedure of Maurer and Bretz (2014):
#'
#' * within a family, the hypotheses are tested with `test` at the family's
#'   local significance level;
#' * a family passes its weight on only when every hypothesis in it is
#'   rejected.
#'
#' @details
#' Families are combined with a weighted Bonferroni test, as the hypotheses of
#' a standard graph are. Within each family, `test` is applied to that
#' family's p-values alone:
#'   * `"hommel"`: the Hommel (1988) procedure, i.e. the closed test with Simes
#'     intersection tests.
#'   * `"hochberg"`: the Hochberg (1988) step-up procedure.
#'   * `"holm"`: the Holm (1979) step-down procedure.
#'
#' Hommel and Hochberg control the family-wise error rate when the test
#' statistics within a family are positively dependent, for example one-sided
#' tests based on multivariate normal statistics with non-negative
#' correlations. This cannot be checked from the p-values; use `"holm"` for a
#' family where negative dependence is possible. For a family with a single
#' hypothesis every test gives the same result.
#'
#' The hypotheses are the columns of the p-value matrix. In [trial_success()],
#' `r1, r2, ...` still refer to these individual hypotheses; the success of a
#' whole family is written as a conjunction, for example `r2 && r3`.
#'
#' @param members A list of integer vectors, one per family, giving the
#'   hypotheses (columns of the p-value matrix) in each family. Together they
#'   must contain each of `1, ..., m` exactly once. If the list is named, the
#'   names are used as family names; otherwise families are named `"F1"`,
#'   `"F2"`, and so on.
#' @inheritParams rlang::args_dots_empty
#' @param test The within-family test: `"hommel"` (default), `"hochberg"` or
#'   `"holm"`. Either a single test used for every family, or one test per
#'   family.
#'
#' @returns A `multigrain_family_spec` object containing:
#'   * `members`: a named list with the (sorted) hypotheses in each family.
#'   * `test`: a named character vector with the within-family test of each
#'     family.
#'
#' @references
#' Maurer, W. and Bretz, F. (2014). A note on testing families of hypotheses
#' using graphical procedures. *Statistics in Medicine*, 33(30), 5340--5346.
#' \doi{10.1002/sim.6307}
#'
#' Hommel, G. (1988). A stagewise rejective multiple test procedure based on
#' a modified Bonferroni test. *Biometrika*, 75(2), 383--386.
#' \doi{10.1093/biomet/75.2.383}
#'
#' Hochberg, Y. (1988). A sharper Bonferroni procedure for multiple tests of
#' significance. *Biometrika*, 75(4), 800--802. \doi{10.1093/biomet/75.4.800}
#'
#' Holm, S. (1979). A simple sequentially rejective multiple test procedure.
#' *Scandinavian Journal of Statistics*, 6(2), 65--70.
#'
#' @seealso [graph_optimise()] and [calc_power_pvals()], which take the family
#'   specification through their `families` argument.
#'
#' @export
#' @examples
#' # Five hypotheses in three families: H1 alone, H2 to H4 together and H5
#' # alone. Hommel is used within every family.
#' family_spec(list(primary = 1, secondary = 2:4, other = 5))
#'
#' # One test per family
#' family_spec(
#'   list(primary = 1, secondary = 2:4, other = 5),
#'   test = c("hommel", "hochberg", "holm")
#' )
family_spec <- function(members, ..., test = "hommel") {
    rlang::check_dots_empty()
    check_family_members(members)
    check_character(test)

    num_fam <- length(members)

    if (!length(test) %in% c(1L, num_fam)) {
        cli::cli_abort(
            "{.arg test} must have length 1 or {num_fam} (one per family), \\
            not {length(test)}."
        )
    }
    test <- vapply(
        rep_len(test, num_fam),
        rlang::arg_match0,
        character(1),
        values = family_tests,
        arg_nm = "test",
        error_call = rlang::current_env(),
        USE.NAMES = FALSE
    )

    fam_names <- names(members) %||% build_family_names(num_fam)
    members <- lapply(members, function(x) sort(as.integer(x)))
    names(members) <- fam_names
    names(test) <- fam_names

    new_family_spec(members = members, test = test)
}

## helpers ----------------------------------------------------------------

#' Check that `members` is a partition of `1, ..., m`
#'
#' @param members A list of integer vectors, as passed to [family_spec()].
#'
#' @returns `NULL`, invisibly, or an error.
#' @noRd
check_family_members <- function(members, call = rlang::caller_env()) {
    if (!rlang::is_bare_list(members) || rlang::is_empty(members)) {
        cli::cli_abort(
            "{.arg members} must be a non-empty list of integer vectors, \\
            not {.obj_type_friendly {members}}.",
            call = call
        )
    }

    for (i in seq_along(members)) {
        check_integerish(
            members[[i]],
            arg = glue::glue("members[[{i}]]"),
            call = call
        )
    }

    # Indices reported in messages are kept as character, so that cli
    # pluralises on how many there are rather than on their value.
    empty <- as.character(which(lengths(members) == 0L))
    if (length(empty) > 0L) {
        cli::cli_abort(
            "Every family must contain at least one hypothesis. \\
            Famil{?y/ies} {empty} {?is/are} empty.",
            call = call
        )
    }

    idx <- unlist(members, use.names = FALSE)
    if (anyNA(idx) || any(idx < 1)) {
        cli::cli_abort(
            "The hypotheses in {.arg members} must be positive whole numbers \\
            (columns of the p-value matrix).",
            call = call
        )
    }

    duplicates <- as.character(sort(unique(idx[duplicated(idx)])))
    if (length(duplicates) > 0L) {
        cli::cli_abort(
            "Each hypothesis must belong to exactly one family. \\
            Hypothes{?is/es} {duplicates} {?appears/appear} in more than one \\
            family.",
            call = call
        )
    }

    not_covered <- as.character(setdiff(seq_len(max(idx)), idx))
    if (length(not_covered) > 0L) {
        cli::cli_abort(
            "{.arg members} must contain every hypothesis from 1 to \\
            {max(idx)}. {cli::qty(length(not_covered))}Hypothes{?is/es} \\
            {not_covered} {?is/are} not in any family.",
            call = call
        )
    }

    fam_names <- names(members)
    if (!is.null(fam_names) && (anyNA(fam_names) || any(fam_names == ""))) {
        cli::cli_abort(
            "{.arg members} must be fully named or not named at all.",
            call = call
        )
    }
    if (anyDuplicated(fam_names) > 0L) {
        cli::cli_abort(
            "Family names in {.arg members} must be unique.",
            call = call
        )
    }

    invisible(NULL)
}

is_family_spec <- function(x) {
    inherits(x, "multigrain_family_spec")
}

check_family_spec <- function(
    families,
    arg = rlang::caller_arg(families),
    call = rlang::caller_env(),
    allow_null = FALSE
) {
    if (!missing(families)) {
        if (is_family_spec(families)) {
            return(invisible(NULL))
        }

        if (allow_null && is.null(families)) {
            return(invisible(NULL))
        }
    }

    rlang::stop_input_type(
        families,
        "a multigrain family specification object",
        allow_null = allow_null,
        arg = arg,
        call = call
    )
}

#' Helper for building family names
#'
#' @param num_fam (integerish) number of families
#'
#' @returns A character vector of names `"F1"`, `"F2"`, ... .
#' @noRd
build_family_names <- function(num_fam) {
    glue::glue("F{rlang::seq2(1, num_fam)}") |>
        vctrs::vec_cast(character())
}

#' Get the number of hypotheses in a `multigrain_family_spec`
#'
#' @param families A `multigrain_family_spec` object.
#'
#' @returns An integer: the number of hypotheses (columns of the p-value
#'   matrix) covered by the families.
#' @noRd
family_spec_get_m <- function(families) {
    check_family_spec(families)

    attr(families, "m")
}

#' Get the family of each hypothesis
#'
#' @param families A `multigrain_family_spec` object.
#'
#' @returns An integer vector of length m: the (1-based) family of each
#'   hypothesis.
#' @noRd
family_spec_get_index <- function(families) {
    check_family_spec(families)

    index <- integer(family_spec_get_m(families))
    index[unlist(families$members)] <- rep(
        seq_along(families$members),
        lengths(families$members)
    )
    index
}

#' Check that p-values and graph match a `multigrain_family_spec`
#'
#' @param families A `multigrain_family_spec` object.
#' @param pvals The p-value matrix, with one column per hypothesis.
#' @param num_nodes The number of nodes of the graph (one per family).
#'
#' @returns `NULL`, invisibly, or an error.
#' @noRd
check_families_dims <- function(
    families,
    pvals,
    num_nodes,
    call = rlang::caller_env()
) {
    m <- family_spec_get_m(families)
    num_fam <- length(families$members)

    if (ncol(pvals) != m) {
        cli::cli_abort(
            "{.arg pvals} has {ncol(pvals)} column{?s} but {.arg families} \\
            covers {m} hypothes{?is/es}.",
            call = call
        )
    }

    if (num_nodes != num_fam) {
        cli::cli_abort(
            c(
                "The graph has {num_nodes} node{?s} but {.arg families} has \\
                {num_fam} famil{?y/ies}.",
                i = "With {.arg families}, the graph has one node per family."
            ),
            call = call
        )
    }

    invisible(NULL)
}

#' Prepare p-values for the shortcut kernel
#'
#' Without families the p-values are used as they are. With families, each
#' family's p-values are replaced by their within-family adjusted p-values
#' (row by row, with the family's `test`), which is what
#' `graph_shortcut_families()` expects.
#'
#' Callers pass raw p-values and call this right before the kernel (or once,
#' outside the objective function), so `pvals` always means raw p-values
#' elsewhere in the package.
#'
#' @param pvals (numeric matrix) Raw p-values (`nsim × m`).
#' @param families `NULL` or a `multigrain_family_spec` object.
#'
#' @returns A numeric matrix with the same dimensions as `pvals`.
#' @noRd
.prepare_pvals <- function(pvals, families) {
    if (is.null(families)) {
        return(pvals)
    }

    adjust_pvals_families(
        pvals,
        family = family_spec_get_index(families) - 1L,
        test = match(families$test, family_tests) - 1L
    )
}

# Print method ------------------------------------------------------------

#' @export
print.multigrain_family_spec <- function(x, ...) {
    if (is.null(x)) {
        return()
    }

    cli::cat_line(cli::format_inline("{.cls {class(x)}}"))
    cli::cat_line("Hypotheses (columns of the p-value matrix) in each family:")
    format_families(x) |>
        cli::cat_line()
    invisible(x)
}

# Summary method ---------------------------------------------------------

#' @export
summary.multigrain_family_spec <- function(object, ...) {
    if (is.null(object)) {
        return()
    }

    cli::cat_line(cli::style_underline("\nFamilies of hypotheses"), ":")
    format_families(object) |>
        cli::cat_line()
    invisible(object)
}

# One line per family: "name (test): 2 3 4".
format_families <- function(families) {
    paste0(
        format(names(families$members)),
        " (",
        families$test,
        "): ",
        vapply(families$members, paste, character(1), collapse = " ")
    )
}
