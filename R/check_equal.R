#' Check for equality with tolerance
#'
#' Account for small floating point numeric differences between numbers:
#'   * `check_equal()` and `vec_check_equal()` check for equality.
#'   * `check_not_equal()` and `vec_check_not_equal()` check for inequality.
#'   * `vec_check_equal()` and `vec_check_not_equal()` are vectorised (they
#'   return a vector of the same length as the input), while `check_equal()` and
#'   `check_not_equal()` return a single value.
#'
#' @param vec (numeric) vector to compare.
#' @param compare_to (numeric) scalar numeric to compare to.
#' @param tolerance numeric >= 0. The tolerance to use for the comparison. If
#'   unspecified the base R default (`sqrt(.Machine$double.eps)`) will be used.
#'   Passed down to `all.equal()`.
#' @param check_attributes (logical) indicates whether to include attributes in
#'   the check.
#'
#' @returns
#'   * `check_equal()` returns a logical scalar (vector of length one)
#'   indicating if all the values of the input vector (`vec`) are equal to the
#'   `compare_to` value.
#'   * `vec_check_equal()` returns a logical vector of the same length as the
#'   input vector (`vec`).
#'   * `check_not_equal()` is the opposite of `check_equal()` as it checks for
#'   inequality. Similarly, it returns a single value (logical scalar).
#'   * `vec_check_not_equal()` returns a logical vector of the same length as
#'   the input vector (`vec`).
#'
#' @noRd
#'
#' @examples
#' # check_equal()
#' multigrain:::check_equal(1 + 1e-12, 1)
#'
#' multigrain:::check_equal(1 + 1e-12, 1, tolerance = 0)
#'
#' # vectorised check_equal()
#' multigrain:::vec_check_equal(c(1 + 1e-12, 1 + 1e-11), 1)
#'
#' multigrain:::vec_check_equal(c(1 + 1e-12, 1 + 1e-11), 1, tolerance = 1e-13)
#'
#' # check_not_equal()
#' multigrain:::check_not_equal(1 + 1e-12, 1)
#'
#' multigrain:::check_not_equal(1 + 1e-12, 1, tolerance = 0)
#'
#' # vectorised check_not_equal()
#' multigrain:::vec_check_not_equal(c(1 + 1e-12, 1 + 1e-11), 1)
#'
#' multigrain:::vec_check_not_equal(
#'     c(1 + 1e-12, 1 + 1e-11),
#'     1,
#'     tolerance = 1e-13
#' )
check_equal <- function(
    vec,
    compare_to = 1,
    tolerance = sqrt(.Machine$double.eps),
    check_attributes = FALSE
) {
    output <- all.equal(
        vec,
        rep(compare_to, length(vec)),
        tolerance = tolerance,
        check.attributes = check_attributes
    )

    if (is.character(output)) {
        output <- FALSE
    }

    output
}

# a vectorised version of check_equal
vec_check_equal <- function(
    vec,
    compare_to = 1,
    tolerance = sqrt(.Machine$double.eps),
    check_attributes = FALSE
) {
    purrr::map_lgl(
        vec,
        check_equal,
        compare_to = compare_to,
        tolerance = tolerance,
        check_attributes = check_attributes
    )
}

check_not_equal <- function(
    vec,
    compare_to = 1,
    tolerance = sqrt(.Machine$double.eps),
    check_attributes = FALSE
) {
    !check_equal(
        vec = vec,
        compare_to = compare_to,
        tolerance = tolerance,
        check_attributes = check_attributes
    )
}

# a vectorised version of check_not_equal
vec_check_not_equal <- function(
    vec,
    compare_to = 1,
    tolerance = sqrt(.Machine$double.eps),
    check_attributes = FALSE
) {
    !vec_check_equal(
        vec = vec,
        compare_to = compare_to,
        tolerance = tolerance,
        check_attributes = check_attributes
    )
}
