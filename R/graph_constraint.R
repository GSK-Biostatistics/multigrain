# `multigrain_graph_constraint` class -------------------------------------

## constructor ------------------------------------------------------------

new_graph_constraint <- function(
    hyp_constraint = double(),
    trans_constraint = matrix(double()),
    names = character(),
    tolerance = sqrt(.Machine$double.eps)
) {
    stopifnot(
        is.double(hyp_constraint),
        is.double(trans_constraint),
        is.matrix(trans_constraint),
        is.character(names),
        is.double(tolerance)
    )

    if (!rlang::is_empty(names)) {
        names(hyp_constraint) <- names
        dimnames(trans_constraint) <- list(names, names)
    }

    structure(
        list(
            hyp_constraint = hyp_constraint,
            trans_constraint = trans_constraint
        ),
        class = "multigrain_graph_constraint",
        m = length(hyp_constraint),
        tolerance = tolerance
    )
}

## user-facing constructor ------------------------------------------------

#' Create a _graph constraint_ object for optimisation procedures
#'
#' A _graph constraint_ object defines constraints on the hypothesis weight
#' vector and transition matrix for optimisation of graph-based multiple testing
#' procedures.
#'
#' The object contains a constraint on the weights (`hyp_constraint`) and a
#' constraint on the transition matrix (`trans_constraint`).
#'
#' @details
#' The _graph constraint_ object is used to define constraints on both the
#' hypothesis weight vector and the transition matrix in graph-based
#' optimisation procedures. The `graph_optimise()` function will read the graph
#' constraints and only optimise free parameters (specified by `NA` in
#' `hyp_constraint` and `trans_constraint`).
#'
#' Either `hyp_constraint` or `trans_constraint` must be provided (they can't
#' both be `NULL` at the same time). If only one is provided, the
#' _graph constraint_ object will allow `graph_optimise()` to optimise any
#' parameter (i.e., no constraints will be specified) in the other one.
#'
#' `hyp_constraint` and `trans_constraint` must refer to the same number of
#' hypotheses.
#'
#' @param hyp_constraint A numeric vector defining the constraints on the
#'   hypothesis weight vector. If `NULL`, all hypothesis weights are free
#'   parameters to be optimised.
#' @param trans_constraint A numeric matrix defining the constraints on the
#'   transition matrix between hypotheses. If `NULL`, all transition matrix
#'   weights are free parameters to be optimised, except for diagonal elements
#'   which remain set to 0.
#' @inheritParams rlang::args_dots_empty
#' @param names An optional character vector containing hypotheses' names. If
#'   not provided it defaults to `"auto"` meaning the hypotheses will be
#'   automatically named `"H1"`, `"H2"`, etc.
#' @param diagnose A logical value enabling detailed diagnosis. Default is
#'   `FALSE`.
#' @param tolerance numeric >= 0. Differences smaller than `tolerance` will not
#'   be reported. The default value is close to `1.5e-8` -
#'   i.e. `sqrt(.Machine$double.eps)` (the standard R definition of
#'   "practically equal", as used by [base::all.equal()]).
#'
#' @returns A multigrain _graph constraint_ object (an S3 list with class
#'   `multigrain_graph_constraint`) containing:
#'   * `hyp_constraint`: a numeric vector representing the constraints on the
#'   hypothesis weight vector.
#'   * `trans_constraint`: a numeric matrix representing the constraints on the
#'   transition matrix.
#' If an element is `NA`, it is a free parameter to be optimised by
#' `graph_optimise()`.
#'
#' @references
#' Xi, D. and Chen, Y. (2024). Optimal weighted Bonferroni tests and
#' their graphical extensions. *Statistics in Medicine*, 43(3),
#' 475--500. <https://doi.org/10.1002/sim.9958>
#'
#' @export
#' @examples
#' # Create a graph_constraint object with predefined weight constraints
#' graph_constraint(hyp_constraint = c(NA, 0.4, NA))
#'
graph_constraint <- function(
    hyp_constraint = NULL,
    trans_constraint = NULL,
    ...,
    names = "auto",
    diagnose = FALSE,
    tolerance = sqrt(.Machine$double.eps)
) {
    if (is.null(hyp_constraint) && is.null(trans_constraint)) {
        cli::cli_abort(
            "{.arg hyp_constraint} and {.arg trans_constraint} cannot both be \\
            {.code NULL} at the same time. At least one must be supplied."
        )
    }

    if (is.null(hyp_constraint)) {
        num_hyp <- nrow(trans_constraint)
        hyp_constraint <- hyp_constraint_free(num_hyp)
    }

    if (is.null(trans_constraint)) {
        num_hyp <- length(hyp_constraint)
        trans_constraint <- trans_constraint_free(num_hyp)
    }

    rlang::check_dots_empty()
    check_character(names)
    check_double(hyp_constraint, allow_na = TRUE)

    hyp_constraint <- vctrs::vec_cast(hyp_constraint, double())

    check_double_matrix(trans_constraint)
    check_logical(diagnose, allow_na = FALSE)
    rlang::check_number_decimal(tolerance, min = 0)

    num_hyp <- length(hyp_constraint)

    if (missing(names)) {
        names <- build_hyp_names(num_hyp)
    }

    if (length(names) != num_hyp) {
        cli::cli_abort(
            "{.arg names} must have {num_hyp} elements. It has {length(names)}."
        )
    }

    validate_graph_constraint(
        new_graph_constraint(
            hyp_constraint = hyp_constraint,
            trans_constraint = trans_constraint,
            names = names,
            tolerance = tolerance
        ),
        diagnose = diagnose
    )
}


## constructor for an unconstrained graph ---------------------------------

#' Create an unconstrained _graph constraint_
#'
#' Create a graph constraint that allows for all hypotheses transitions to be
#' optimised.
#'
#' @param num_hyp An integer denoting the number of hypotheses. Must be greater
#'   than or equal to 2.
#' @inheritParams graph_constraint names
#'
#' @returns An unconstrained `multigrain_graph_constraint` object. In
#'   `hyp_constraint` all values are `NA` and similarly in the
#'   `trans_constraint`, except for the diagonal which is set to `0`. If we have
#'   only 2 hypotheses, then the matrix will have 0 on the diagonal and the
#'   other 2 elements are set to 1.
#'
#' @export
#' @examples
#' # Create a graph constraint object with 3 hypotheses and no constraints
#' graph_constraint_free(3)
#'
#' # Create a graph constraint object with 2 hypotheses results in set values in
#' # the transition matrix
#' graph_constraint_free(2)
graph_constraint_free <- function(num_hyp, names = "auto") {
    rlang::check_number_whole(num_hyp, min = 2)
    check_character(names)

    hyp_constraint <- hyp_constraint_free(num_hyp)
    trans_constraint <- trans_constraint_free(num_hyp)

    if (missing(names)) {
        names <- build_hyp_names(num_hyp)
    }

    if (length(names) != num_hyp) {
        cli::cli_abort(
            "{.arg names} must have {num_hyp} elements. It has {length(names)}."
        )
    }

    graph_constraint(
        hyp_constraint = hyp_constraint,
        trans_constraint = trans_constraint,
        names = names
    )
}

## helpers ----------------------------------------------------------------

is_graph_constraint <- function(x) {
    inherits(x, "multigrain_graph_constraint")
}

check_graph_constraint <- function(
    graph_constraint,
    arg = rlang::caller_arg(graph_constraint),
    call = rlang::caller_env(),
    allow_null = FALSE
) {
    if (!missing(graph_constraint)) {
        if (is_graph_constraint(graph_constraint)) {
            return(invisible(NULL))
        }

        if (allow_null && is.null(graph_constraint)) {
            return(invisible(NULL))
        }
    }

    rlang::stop_input_type(
        graph_constraint,
        "a multigrain graph constraint object",
        allow_null = allow_null,
        arg = arg,
        call = call
    )
}


#' Helper for building hypothesis names
#'
#' @param num_hyp (integerish) number of hypotheses
#'
#' @returns A character vector of names the same length as `hyp_constraint`.
#'   Names are `"H1"`, ... .
#' @noRd
build_hyp_names <- function(num_hyp) {
    glue::glue("H{rlang::seq2(1, num_hyp)}") |>
        vctrs::vec_cast(character())
}

#' Get the `multigrain_graph_constraint` names
#'
#' @inheritParams graph_optimise graph_constraint
#'
#' @returns A character vector of names extracted from the `hyp_constraint`
#'   element.
#' @noRd
graph_constraint_get_names <- function(graph_constraint) {
    check_graph_constraint(graph_constraint)

    names(graph_constraint$hyp_constraint)
}

#' Get the `multigrain_graph_constraint` dimensions
#'
#' @inheritParams graph_optimise graph_constraint
#'
#' @returns An integer representing the number of hypotheses.
#' @noRd
graph_constraint_get_m <- function(graph_constraint) {
    check_graph_constraint(graph_constraint)

    attr(graph_constraint, "m")
}

#' Get the `multigrain_graph_constraint` tolerance
#'
#' @inheritParams graph_optimise graph_constraint
#'
#' @returns A numeric scalar; the sum-to-one tolerance stored on the constraint.
#' @noRd
graph_constraint_get_tolerance <- function(graph_constraint) {
    check_graph_constraint(graph_constraint)

    attr(graph_constraint, "tolerance")
}

#' Helper function to build `hyp_constraint` from `m`
#'
#' @param num_hyp (integerish) number of hypotheses
#'
#' @returns A numeric vector (all `NA`) the same length as the number of
#' hypotheses.
#' @noRd
hyp_constraint_free <- function(num_hyp) {
    rep(NA_real_, num_hyp)
}

#' Helper function to build `trans_constraint` from `m`
#'
#' @param num_hyp (integerish) number of hypotheses
#'
#' @returns A `m * m` numeric matrix vector (all `NA`) with 0 on the diagonal.
#' If `m` is 2, then all remaining `NA`s are set to 1.
#' @noRd
trans_constraint_free <- function(num_hyp) {
    trans_constraint <- rep(NA_real_, num_hyp * num_hyp) |>
        matrix(
            nrow = num_hyp,
            byrow = TRUE
        )

    diag(trans_constraint) <- 0

    if (num_hyp == 2) {
        trans_constraint[is.na(trans_constraint)] <- 1
    }

    trans_constraint
}

# Print method ------------------------------------------------------------

#' @export
print.multigrain_graph_constraint <- function(x, ...) {
    if (is.null(x)) {
        return()
    }

    cli::cat_line(cli::format_inline("{.cls {class(x)}}"))

    cli::cat_line("Constraints on hypothesis weights:")
    print(x$hyp_constraint)
    cli::cat_line()
    cli::cat_line("Constraints on transition matrix:")
    print(x$trans_constraint)
    invisible(x)
}


# Summary method ---------------------------------------------------------

#' @export
summary.multigrain_graph_constraint <- function(object, ...) {
    if (is.null(object)) {
        return()
    }

    cli::cat_line(cli::style_underline("\nGraph constraints"), ":")

    cli::cat_line("Constraints on hypothesis weights:")
    print(object$hyp_constraint)

    cli::cat_line("\nConstraints on transition matrix:")
    print(object$trans_constraint)
    invisible(object)
}

# Set and get methods -----------------------------------------------------

# nolint start: object_length_linter

#' @export
`[<-.multigrain_graph_constraint` <- function(x, i, ..., value) {
    x[[i, ...]] <- value

    x
}

#' @export
`[[<-.multigrain_graph_constraint` <- function(x, i, ..., value) {
    dots <- rlang::list2(...)

    if (!"tolerance" %in% names(dots)) {
        dots$tolerance <- attr(x, "tolerance")
    }

    if (!"diagnose" %in% names(dots)) {
        dots$diagnose <- FALSE
    }

    hyp_constraint_names <- names(x$hyp_constraint)

    if (is.vector(value) && rlang::is_named(value)) {
        hyp_constraint_names <- names(value)
    }

    if (is.matrix(value) && !rlang::is_empty(dimnames(value)[[1]])) {
        hyp_constraint_names <- dimnames(value)[[1]]
    }

    if (i == "hyp_constraint") {
        output <- rlang::inject(
            graph_constraint(
                hyp_constraint = value,
                trans_constraint = x$trans_constraint,
                names = hyp_constraint_names,
                !!!(dots)
            )
        )
    }

    if (i == "trans_constraint") {
        output <- rlang::inject(
            graph_constraint(
                hyp_constraint = x$hyp_constraint,
                trans_constraint = value,
                names = hyp_constraint_names,
                !!!(dots)
            )
        )
    }

    output
}

#' @export
`$<-.multigrain_graph_constraint` <- function(x, i, value) {
    x[[i, ]] <- value

    x
}

# nolint end

# Find closest constraint-obeying graph to an existing graph conditional on
# constraints

# Proportional fill of free entries with fixed entries pinned; all inputs >= 0
# Closest graph subject to constraints (minimal change):
# - No constraints  -> return input unchanged
# - Some fixed vals -> pin them; scale only the free part to make sum = 1
closest_graph_to_constraints <- function(
    gc,
    w,
    G,
    tolerance = sqrt(.Machine$double.eps)
) {
    stopifnot(is_graph_constraint(gc))
    m <- graph_constraint_get_m(gc)
    stopifnot(is.numeric(w), length(w) == m, is.matrix(G), dim(G) == c(m, m))

    # weights
    w_out <- w
    fixed_w <- which(!is.na(gc$hyp_constraint))
    if (length(fixed_w)) {
        w_out[fixed_w] <- gc$hyp_constraint[fixed_w]
        s_fixed <- sum(w_out[fixed_w])
        if (s_fixed > 1 + tolerance) {
            stop(
                sprintf(
                    "Infeasible weight constraints: fixed sum=%.6f (>1).",
                    s_fixed
                ),
                call. = FALSE
            )
        }
        free_w <- setdiff(seq_len(m), fixed_w)
        if (length(free_w)) {
            s_free_input <- sum(w[free_w])
            # if no free mass but need to fill, split evenly
            if (s_free_input <= tolerance) {
                w_out[free_w] <- (1 - s_fixed) / length(free_w)
            } else {
                w_out[free_w] <- w[free_w] * ((1 - s_fixed) / s_free_input)
            }
        } else if (abs(s_fixed - 1) > tolerance) {
            # no free entries; must already sum to 1
            stop(
                sprintf(
                    "Infeasible: no free weights and fixed sum=%.6f (!=1).",
                    s_fixed
                ),
                call. = FALSE
            )
        }
    }

    # transition matrix
    g_out <- G
    tc <- gc$trans_constraint

    for (i in seq_len(m)) {
        # fixed indices in row i
        fixed_g <- which(!is.na(tc[i, ]))
        if (length(fixed_g)) {
            g_out[i, fixed_g] <- tc[i, fixed_g]
            s_fixed <- sum(g_out[i, fixed_g])
            if (s_fixed > 1 + tolerance) {
                stop(
                    sprintf(
                        "Infeasible row %d constraints: fixed sum=%.6f (>1).",
                        i,
                        s_fixed
                    ),
                    call. = FALSE
                )
            }
            free_g <- setdiff(seq_len(m), fixed_g)
            if (length(free_g)) {
                s_free_input <- sum(G[i, free_g])
                if (s_free_input <= tolerance) {
                    g_out[i, free_g] <- (1 - s_fixed) / length(free_g)
                } else {
                    g_out[i, free_g] <- G[i, free_g] *
                        ((1 - s_fixed) / s_free_input)
                }
            } else if (abs(s_fixed - 1) > tolerance) {
                stop(
                    sprintf(
                        "Infeasible: row %d has no free entries and fixed sum=%.6f (!=1).", # nolint
                        i,
                        s_fixed
                    ),
                    call. = FALSE
                )
            }
        }
    }

    list(
        hyp_weight = w_out,
        trans_matrix = g_out
    )
}
