#' Simulate raw p-values for a group sequential design
#'
#' `simulate_pvalues_gsd()` simulates raw p-values at several analyses
#' ("looks") of a group sequential trial, under the canonical joint normal
#' model for the test statistics. It is the group sequential counterpart of
#' [simulate_pvalues()]: `power_nominal` and `corr_matrix` describe the
#' *final-look* statistics, and the information fractions say how much of that
#' information is available at each analysis.
#'
#' @details The non-centrality parameter \eqn{\Delta_i} of each hypothesis is
#'   obtained from its nominal power with [calc_ncp()], exactly as in
#'   [simulate_pvalues()]. Writing \eqn{t_{i,k}} for the information fraction
#'   of hypothesis \eqn{i} at analysis \eqn{k}, the test statistics follow the
#'   canonical joint distribution
#'
#'   \deqn{E[Z_{i,k}] = \Delta_i \sqrt{t_{i,k}}, \qquad
#'   \mathrm{Corr}(Z_{i,k}, Z_{j,l}) = \rho_{ij}
#'   \sqrt{\min(t_{i,k}, t_{j,l}) / \max(t_{i,k}, t_{j,l})},}
#'
#'   where \eqn{\rho_{ij}} is the entry of `corr_matrix` for hypotheses
#'   \eqn{i} and \eqn{j}. Two analyses of the same hypothesis at the same
#'   information fraction are therefore perfectly correlated, and two different
#'   hypotheses observed at the same information fraction have correlation
#'   \eqn{\rho_{ij}}. P-values are \eqn{1 - \Phi(Z_{i,k})}.
#'
#'   Only the analyses that carry distinct information are simulated. Once a
#'   hypothesis reaches full information its statistic no longer changes, so
#'   the p-value at the analysis where it matured is copied forward to every
#'   later analysis; drawing those columns instead would make the covariance
#'   singular. An analysis where a hypothesis has no data (`NA` information
#'   fraction) before it matures is `NA` in the returned array;
#'   [transform_pvalues_gsd()] turns that into a repeated p-value of 1, which
#'   is "cannot reject".
#'
#' @param power_nominal A numeric vector of nominal power values for each
#'   hypothesis, at its final analysis.
#' @inheritParams rlang::args_dots_empty
#' @param alpha A number giving the overall one-sided significance level used
#'   to convert `power_nominal` into non-centrality parameters. Default is
#'   `0.025`.
#' @param corr_matrix A numeric matrix representing the correlation matrix
#'   \eqn{\Sigma} of the test statistics at the final analysis.
#' @param info_frac Information fractions. Either a numeric vector of length
#'   `K`, applied to every hypothesis, or an `m` by `K` numeric matrix with
#'   `NA` at analyses where a hypothesis has no data. They must be positive and
#'   non-decreasing over the analyses where a hypothesis has data, and every
#'   hypothesis must have data at some analysis.
#' @param nsim An integer indicating the number of simulations to run. Default
#'   is `1e5`.
#'
#' @returns A numeric array of raw p-values with dimensions `c(nsim, m, K)`:
#'   simulations by hypotheses by analyses. The `m` by `K` matrix of
#'   information fractions is attached as an `info_frac` attribute, where
#'   [transform_pvalues_gsd()] picks it up.
#'
#' @export
#' @examples
#' # two hypotheses analysed at 60% and 100% of the planned information
#' corr <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#'
#' set.seed(1)
#' pvals <- simulate_pvalues_gsd(
#'     c(0.8, 0.9),
#'     corr_matrix = corr,
#'     info_frac = c(0.6, 1),
#'     nsim = 100
#' )
#'
#' dim(pvals)
#' attr(pvals, "info_frac")
#'
#' # the second hypothesis has no data at the interim analysis, and the first
#' # matures there: its p-value is carried forward to the final analysis
#' set.seed(1)
#' pvals_na <- simulate_pvalues_gsd(
#'     c(0.8, 0.9),
#'     corr_matrix = corr,
#'     info_frac = rbind(c(1, 1), c(NA, 1)),
#'     nsim = 100
#' )
#'
#' identical(pvals_na[, 1L, 1L], pvals_na[, 1L, 2L])
#' all(is.na(pvals_na[, 2L, 1L]))
simulate_pvalues_gsd <- function(
    power_nominal,
    ...,
    alpha = 0.025,
    corr_matrix = diag(length(power_nominal)),
    info_frac, # nolint: function_argument_linter. required, named after `...`
    nsim = 1e5
) {
    check_double(power_nominal)
    rlang::check_dots_empty()
    rlang::check_required(info_frac)
    rlang::check_number_decimal(alpha, min = 0, max = 1)
    check_double_matrix(corr_matrix)
    rlang::check_number_whole(nsim, min = 1)

    m <- length(power_nominal)
    if (!identical(dim(corr_matrix), c(m, m))) {
        cli::cli_abort(
            "{.arg corr_matrix} must be a {m} by {m} matrix, one row and \\
            column per value in {.arg power_nominal}."
        )
    }

    nsim <- as.integer(nsim)
    info_frac <- .gsd_sim_info_frac(info_frac, m = m)
    n_look <- ncol(info_frac)

    # the analyses carrying distinct information, and the one at which each
    # hypothesis matures (shared with the transform, design record 4.1)
    look_info <- lapply(
        seq_len(m),
        function(i) .gsd_looks(info_frac[i, ], hyp = i)
    )
    looks <- lapply(look_info, `[[`, "looks")
    maturity <- vapply(look_info, `[[`, integer(1L), "maturity")

    # one multivariate normal over every (hypothesis, distinct analysis) pair
    hyp_index <- rep(seq_len(m), lengths(looks))
    info_vec <- unlist(
        lapply(seq_len(m), function(i) info_frac[i, looks[[i]]]),
        use.names = FALSE
    )

    ncp <- calc_ncp(power_nominal, alpha = alpha)
    mean_vec <- ncp[hyp_index] * sqrt(info_vec)
    covariance <- corr_matrix[hyp_index, hyp_index, drop = FALSE] *
        sqrt(
            outer(info_vec, info_vec, pmin) / outer(info_vec, info_vec, pmax)
        )

    pvals <- stats::pnorm(
        mvtnorm::rmvnorm(nsim, mean = mean_vec, sigma = covariance),
        lower.tail = FALSE
    )

    out <- array(NA_real_, dim = c(nsim, m, n_look))
    column <- 0L
    for (i in seq_len(m)) {
        for (k in looks[[i]]) {
            column <- column + 1L
            out[, i, k] <- pvals[, column]
        }
        # a matured statistic does not change: copy its p-value forward
        if (maturity[[i]] < n_look) {
            later <- seq.int(maturity[[i]] + 1L, n_look)
            out[, i, later] <- out[, i, maturity[[i]]]
        }
    }

    attr(out, "info_frac") <- info_frac

    out
}


# The m by K matrix of information fractions, from a length-K vector recycled
# over hypotheses or from a matrix supplied as is. The number of analyses is
# whatever the user supplied, so it is read from `info_frac` itself.
.gsd_sim_info_frac <- function(info_frac, m, call = rlang::caller_env()) {
    if (!is.numeric(info_frac) || length(info_frac) == 0L) {
        cli::cli_abort(
            "{.arg info_frac} must be a non-empty numeric vector or matrix.",
            call = call
        )
    }

    if (is.matrix(info_frac)) {
        if (nrow(info_frac) != m) {
            cli::cli_abort(
                "{.arg info_frac} must have one row per hypothesis: {m}, \\
                not {nrow(info_frac)}.",
                call = call
            )
        }
        out <- info_frac
    } else {
        out <- matrix(
            info_frac,
            nrow = m,
            ncol = length(info_frac),
            byrow = TRUE
        )
    }

    storage.mode(out) <- "double"

    out
}
