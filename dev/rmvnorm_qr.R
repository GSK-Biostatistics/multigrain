#' Quasi-random multivariate normal draws (Sobol' + $\Phi^{-1}$)
#'
#' Generates `n` draws from $\mathcal{N}_p(\mu, \Sigma)$ using a
#' low-discrepancy Sobol' sequence with a digital shift, mapped through the
#' standard-normal quantile and then linearly transformed to impose the target
#' covariance. This is a lightweight, dependency-free replacement for
#' `gMCP::rqmvnorm()` based on randomised QMC (RQMC).
#'
#' @param n Integer. Number of draws to return.
#' @param mean Numeric vector of length $p$. Mean vector $\mu$.
#' @param sigma Symmetric, positive semi-definite $p \times p$
#'   covariance matrix $\Sigma$.
#' @param round_pow2 Logical; default `FALSE`. If `TRUE`, for
#'   quasi-random generation the function rounds `n` up to the next power
#'   of two (at least $2^{10}$), mimicking the embedded nets used by classic
#'   lattice rules. Leave `FALSE` for exact-`n` Sobol' draws.
#' @param randomize Character; one of `"digital.shift"` (default) or
#'   `"none"`. Randomisation is recommended to obtain unbiased estimators
#'   and valid error assessment under RQMC.
#'
#' @return An $n \times p$ matrix. Column names are propagated from
#'   `mean` if present.
#'
#' @details
#' Points $U \in (0,1)^p$ are generated via a Sobol' sequence
#' ([qrng::sobol()]), optionally with a digital shift;
#' then $Z = \Phi^{-1}(U)$ and finally $X = Z L^\top + \mu$, where
#' $L$ is a factor of $\Sigma$. By default a pivoted Cholesky is
#' attempted; if $\Sigma$ is numerically rank-deficient the function falls
#' back to an SVD-based square root.
#'
#' Compared with rank-1 lattice rules with Cranley-Patterson random shifts (as
#' used in {gMCP}), randomised Sobol' provides similar variance reduction
#' for smooth Gaussian problems and does not require the sample size to be a
#' power of two. Setting `round_pow2 = TRUE` retains the power-of-two
#' behaviour when needed for strict comparability or embedded-net workflows.
#'
#' @references
#'
#' - Sobol, I. M. (1967). On the distribution of points in a cube and the
#'   approximation of integrals. *USSR Comput. Math. Math. Phys.*
#'
#' - Cranley, R., & Patterson, T. N. L. (1976). Randomization of number
#'   theoretic methods for multiple integration. *SIAM J. Numer. Anal.*
#'
#' - Owen, A. B. (1995). Randomly permuted (t, m, s)-nets and (t, s)-sequences.
#'   In *Monte Carlo and Quasi-Monte Carlo Methods*.
#'
#' - Joe, S., & Kuo, F. Y. (2008). Constructing Sobol sequences with better
#'   two-dimensional projections. *ACM TOMS*.
#'
#' - Cools, R., Kuo, F. Y., & Nuyens, D. (2006). Constructing embedded lattice
#'   rules for multivariate integration. *SIAM J. Sci. Comput.*
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' mu <- c(0, 0, 0)
#' S  <- matrix(c(1, .3, .1, .3, 1, .4, .1, .4, 1), 3, 3)
#' x  <- rqmvnorm_qr(5000, mean = mu, sigma = S)  # exact-n Sobol'
#' colMeans(x)             # ≈ mu
#' cov(x)                  # ≈ S
#'
#' # Power-of-two rounding for embedded-net workflows:
#' x2 <- rqmvnorm_qr(5000, mean = mu, sigma = S, round_pow2 = TRUE)
#' }
#'
#' @noRd
rmvnorm_qr <- function(
    n,
    mean,
    sigma,
    round_pow2 = FALSE,
    randomize = c("digital.shift", "none")
) {
    randomize <- match.arg(randomize)

    if (!is.matrix(sigma)) {
        stop(
            "`sigma` must be a matrix.",
            call. = FALSE
        )
    }
    if (!isSymmetric(sigma, tolerance = sqrt(.Machine$double.eps))) {
        stop(
            "`sigma` must be symmetric.",
            call. = FALSE
        )
    }
    if (length(mean) != nrow(sigma)) {
        stop(
            "`mean` and `sigma` have non-conforming sizes.",
            call. = FALSE
        )
    }
    p <- length(mean)

    ## Choose effective n
    if (round_pow2) {
        log2n <- max(ceiling(log2(n)), 10L)
        if (log2n > 30L) {
            stop(
                "Requested n rounds above 2^30; refuse to generate that many points.",
                call. = FALSE
            )
        }
        n_eff <- 2^log2n
    } else {
        n_eff <- as.integer(n)
    }

    ## Factor Sigma: try pivoted Cholesky, else SVD-based square-root
    use_chol <- TRUE
    l <- try(chol(sigma, pivot = TRUE), silent = TRUE)
    if (inherits(l, "try-error")) {
        use_chol <- FALSE
    } else {
        rnk <- attr(l, "rank")
        if (!is.null(rnk) && rnk < p) use_chol <- FALSE
    }

    if (!use_chol) {
        sv <- svd(sigma)
        if (!all(sv$d >= -sqrt(.Machine$double.eps) * abs(sv$d[1]))) {
            stop(
                "`sigma` has negative eigenvalues (beyond numerical tolerance).",
                call. = FALSE
            )
        }
        sqrt_sigma <- t(sv$v %*% (t(sv$u) * sqrt(pmax(sv$d, 0))))
    }

    ## Sobol’ points in (0,1)^p, optional digital shift
    u <- qrng::sobol(
        n = n_eff,
        d = p,
        randomize = if (randomize == "none") NULL else "digital.shift"
    )
    z <- stats::qnorm(u)

    ## Impose covariance and add mean
    if (use_chol) {
        piv <- attr(l, "pivot")
        y <- z %*% l
        x <- y[, order(piv), drop = FALSE]
    } else {
        x <- z %*% sqrt_sigma
    }

    x <- sweep(x, 2, mean, "+")
    if (!is.null(names(mean))) {
        colnames(x) <- names(mean)
    }

    if (n_eff != n) {
        x <- x[seq_len(n), , drop = FALSE]
    }
    x
}
