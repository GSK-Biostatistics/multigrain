# Function to convert vector missing diagonal back to matrix
vec_to_mat <- function(vec, n, byrow = TRUE) {
    # Create empty vector of zeros of length n^2
    x <- rep(0, n^2)
    j <- 1

    # Replace zeros with vec elements unless diagonal and return nxn matrix
    for (i in 1:n^2) {
        if (!(i %in% get_diag_idx(n))) {
            x[i] <- vec[j]
            j <- j + 1
        }
    }

    matrix(
        x,
        nrow = n,
        byrow = byrow
    )
}

# Function to get diagonal indices of matrix
get_diag_idx <- function(n) {
    k <- 1:n
    indices <- n * (k - 1) + k
    indices
}


#' Calculate non-centrality parameter
#'
#' Calculate the non-centrality parameter of the test statistic that corresponds
#' to the statistical power for a given number of hypotheses in a trial.
#'
#' @details
#' `calc_ncp()` uses the quantiles of the standard normal distribution to
#' compute the non-centrality parameter. It is given by the difference between
#' the critical value for a standard normal distribution and the quantile
#' corresponding to the specified power.
#'
#' @param power A numeric vector indicating each hypothesis' statistical power
#'   (probability of rejecting the null hypothesis if it false).
#' @inheritParams graph_optimise alpha
#'
#' @returns A numeric vector of non-centrality parameters corresponding to each
#' element in `power`.
#'
#' @export
#' @examples
#' # Basic usage
#' calc_ncp(power = c(0.8, 0.9))
#'
#' # Custom significance level
#' calc_ncp(power = c(0.8, 0.9), alpha = 0.01)
calc_ncp <- function(power, alpha = 0.025) {
    check_double(power)
    rlang::check_number_decimal(alpha, min = 0, max = 1)

    rep(stats::qnorm(1 - alpha), length(power)) - stats::qnorm(1 - power)
}


t_stat_power <- function(npc, alpha = 0.025) {
    1 - stats::pnorm(stats::qnorm(1 - alpha) - npc)
}


#' Check the validity of a graph-based MTP
#'
#' Check that a given weight vector and transition matrix satisfy several
#' conditions to determine if they form a valid graph. It verifies that the
#' matrix is square, the dimensions match the length of the weight vector,
#' the weight vector elements and matrix elements lie within the interval
#' \[0, 1\], and that the weight vector sums to 1. Additionally, it checks the
#' matrix is properly normalised, with all diagonal elements equal to 0 and each
#' row summing to 1.
#'
#' @param hyp_weight A numeric vector representing the weights of the nodes
#'   (hypotheses) in the graph. All elements must be in the interval \[0, 1\]
#'   and the vector must sum to 1.
#' @param trans_matrix A numeric matrix representing the transition matrix
#'   between nodes. All elements must be in the interval \[0, 1\], with diagonal
#'   elements equal to 0, and each row must sum to 1.
#' @param sum_to_one_constraint A logical indicating whether to allow graphs
#'   where both the transition matrix rows are not constrained to sum-to-one,
#'   for example in a fixed sequence. Defaults to `TRUE`.
#' @param tolerance numeric >= 0. The tolerance when evaluating the sum-to-one
#'   constraints of the hypothesis weights and transition matrix rows. The
#'   default value is close to `1.5e-8` - i.e. `sqrt(.Machine$double.eps)` (the
#'   standard R definition of  "practically equal", as used by
#'   [base::all.equal()]).
#'
#' @returns A logical value: `TRUE` if the graph is valid, otherwise `FALSE`.
#' The function issues warnings if any of the validity checks fail.
#'
#' @details
#' The function performs the following checks to determine graph validity:
#'   * `trans_matrix` is a square matrix.
#'   * the length of `hyp_weight` matches the number of rows and columns in
#'   `trans_matrix`.
#'   * all diagonal elements of `trans_matrix` are 0.
#'   * all elements of `hyp_weight` and `trans_matrix` lie in the interval
#'   \[0, 1\].
#'   * the sum of `hyp_weight` equals 1.
#'   * each row of `trans_matrix` sums to 1.
#'
#' If any of these conditions are not satisfied, the function returns `FALSE`
#' and issues the corresponding warning message.
#'
#' @export
#' @examples
#' hyp_weight <- c(0.4, 0.3, 0.3)
#' trans_matrix <- matrix(
#'                     c(0, 0.5, 0.5, 0.3, 0, 0.7, 0.6, 0.4, 0),
#'                     nrow = 3,
#'                     byrow = TRUE
#'                 )
#' is_graph_valid(hyp_weight, trans_matrix)
is_graph_valid <- function(
    hyp_weight,
    trans_matrix,
    sum_to_one_constraint = TRUE,
    tolerance = sqrt(.Machine$double.eps)
) {
    check_logical(sum_to_one_constraint)
    rlang::check_number_decimal(tolerance, min = 0)

    # Check if the matrix is square
    if (!is.matrix(trans_matrix)) {
        warning(
            "trans_matrix is not a matrix.",
            call. = FALSE
        )
        return(FALSE)
    }
    if (nrow(trans_matrix) != ncol(trans_matrix)) {
        warning(
            "`trans_matrix` is not a square matrix.",
            call. = FALSE
        )
        return(FALSE)
    }

    # Check if the length of the weight vector matches the matrix dimensions
    if (length(hyp_weight) != nrow(trans_matrix)) {
        warning(
            "Length of `hyp_weight` does not match the dimensions of `trans_matrix`.", # nolint
            call. = FALSE
        )
        return(FALSE)
    }

    # Check if matrix diagonals are all 0
    if (any(diag(trans_matrix) != 0)) {
        warning(
            "Diagonals of `trans_matrix` are not all zero.",
            call. = FALSE
        )
        return(FALSE)
    }

    # Check if all elements of the weight vector are in [0, 1]
    if (any(hyp_weight < 0 | hyp_weight > 1)) {
        warning(
            "`hyp_weight` contains elements outside the interval [0, 1].",
            call. = FALSE
        )
        return(FALSE)
    }

    # Check if all elements of the matrix are in [0, 1]
    if (any(trans_matrix < 0 | trans_matrix > 1)) {
        warning(
            "`trans_matrix` contains elements outside the interval [0, 1].",
            call. = FALSE
        )
        return(FALSE)
    }

    s <- sum(hyp_weight)
    if (is.na(s) || abs(s - 1) > tolerance) {
        warning(
            "`hyp_weight` does not sum to 1 within tolerance.",
            call. = FALSE
        )
        return(FALSE)
    }

    if (sum_to_one_constraint) {
        row_sums <- rowSums(trans_matrix)
        if (anyNA(row_sums) || any(abs(row_sums - 1) > tolerance)) {
            warning(
                "One or more rows of `trans_matrix` do not sum to 1 within tolerance.", # nolint: line_length_linter
                call. = FALSE
            )
            return(FALSE)
        }
    }

    TRUE
}


#' Normalise graph weights to sum to a target value
#'
#' Ensures that a numeric vector of graph weights is rescaled so its elements
#' sum to a target value (default 1). This is useful when exporting optimised
#' graphs (e.g. to gMCP), where weights must form a valid probability vector.
#'
#' @details
#' Many downstream tools (e.g. gMCP) require hypothesis weights to sum to 1.
#' Direct division by `sum(x)` may leave tiny discrepancies due to
#' floating-point rounding. This helper fixes such issues automatically,
#' ensuring exported graphs are valid.
#'
#' The algorithm first scales all free (non-fixed) elements proportionally so
#' they sum to `target - sum(x[fixed_idx])`. It then computes the largest
#' free element (the "anchor") as the exact complement of all other elements
#' (`target - sum(x[-anchor])`), absorbing any remaining rounding residual. A
#' single additive fallback pass handles the rare case where the prior scaling
#' introduced enough rounding for the complement to land 1 ULP off.
#'
#' @note
#' Two edge cases return early without enforcing `target`:
#'
#' * If all elements of `x` are zero, the zero vector is returned
#'   unchanged (proportional scaling is undefined for an all-zero input).
#' * If every index is in `fixed_idx` (no free elements), the input
#'   is returned as-is — the caller is responsible for ensuring
#'   `sum(x) == target` when no elements may be modified.
#'
#' @param x A numeric vector of weights to be normalised.
#' @inheritParams rlang::args_dots_empty
#' @param fixed_idx An integer vector representing the indices of elements that
#'   must not be modified (e.g. diagonal entries in a transition matrix, or
#'   positions locked by [graph_constraint()]). Defaults to `integer(0)` (no
#'   fixed elements).
#' @param target A number representing the desired sum for the output vector.
#'   Default is `1`.
#' @param tolerance numeric >= 0. The tolerance to be used when checking whether
#'   the sum has converged to `target`. The default value is close to `1.5e-8` -
#'   i.e. `sqrt(.Machine$double.eps)` (the standard R definition of
#'   "practically equal", as used by [base::all.equal()]).
#'
#' @returns A numeric vector of the same length as `x`, adjusted so that
#'   `sum(x)` is within `tolerance` of `target`.
#'
#' @export
#' @examples
#'
#' x <- c(0.2, 0.3, 0.500000001)
#' print(sum(x) == 1)
#'
#' x_sum_to_1 <- normalise_sum(x)
#' print(all.equal(sum(x_sum_to_1), 1))
#'
#' # With fixed elements: indices 1 and 3 are held constant
#' x2 <- c(0.25, 0.45, 0.30)
#' normalise_sum(x2, fixed_idx = c(1L, 3L))
normalise_sum <- function(
    x,
    ...,
    fixed_idx = integer(0),
    target = 1,
    tolerance = sqrt(.Machine$double.eps)
) {
    check_double(x)
    rlang::check_dots_empty()
    check_integerish(fixed_idx)
    rlang::check_number_decimal(target)
    rlang::check_number_decimal(tolerance)

    if (sum(x) == 0) {
        return(x)
    }

    free <- setdiff(seq_along(x), fixed_idx)
    if (length(free) == 0L) {
        return(x)
    }

    # Scale free elements proportionally to approximate the target
    target_free <- target - sum(x[fixed_idx])
    if (target_free < -tolerance) {
        stop(
            "Fixed elements sum to ",
            sum(x[fixed_idx]),
            " which exceeds target (",
            target,
            "). ",
            "Cannot normalise: reduce fixed element values or increase target.",
            call. = FALSE
        )
    }
    s_free <- sum(x[free])
    if (s_free > 0) {
        x[free] <- x[free] * (target_free / s_free)
    }

    anchor <- free[which.max(x[free])]
    x[anchor] <- target - sum(x[-anchor])

    # Single additive fallback: the proportional scaling step may introduce
    # enough rounding that the direct complement lands 1 ULP off.
    if (!check_equal(sum(x), target, tolerance = tolerance)) {
        x[anchor] <- x[anchor] + (target - sum(x))
    }

    # Clamp if rounding pushed the anchor negative
    if (x[anchor] < 0) {
        x[anchor] <- 0
    }

    x
}

bullets_with_header <- function(header, x) {
    if (length(x) == 0) {
        return()
    }

    cli::cat_line(cli::format_inline("{.strong {header}}"))
    nms <- names(x)
    vals <- purrr::map_chr(x, as_simple)
    # options
    glue::glue("{cli::col_blue(nms)}: {vals}") |>
        cli::cat_bullet()
}

as_simple <- function(x) {
    # nolint start: unnecessary_nesting_linter
    if (is.atomic(x) && length(x) == 1) {
        if (is.character(x)) {
            paste0('"', x, '"')
        } else {
            format(x)
        }
    } else {
        if (length(x) == 0) {
            "<unset>"
        } else {
            paste0("<", class(x)[[1L]], ">")
        }
    }
    # nolint end
}

modify_list <- function(
    .x,
    ...,
    call = rlang::caller_env()
) {
    dots <- rlang::list2(...)
    if (length(dots) == 0) {
        return(.x)
    }

    if (!rlang::is_named(dots)) {
        cli::cli_abort(
            "All components of {.arg ...} must be named.",
            call = call
        )
    }

    out <- .x[!names(.x) %in% names(dots)]

    out <- c(out, purrr::compact(dots))

    if (length(out) == 0) {
        names(out) <- NULL
    }

    out
}
