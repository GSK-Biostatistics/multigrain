#' Simulate raw p-values
#'
#' `simulate_pvalues()` simulates raw p-values from under the alternative
#' hypotheses and the assumption that the distribution of test statistics is a
#' multivariate normal distribution. The most important input parameters are
#' the nominal power of each hypothesis, the significance level, and a
#' correlation matrix.
#'
#' @details It starts by calculating the non-centrality parameter \eqn{\Delta},
#'   of the a test statistic z - where \eqn{z \sim N(\Delta, 1)} - based on the
#'   nominal power (probability of rejection under the alternative, unadjusted
#'   for multiplicity) and significance level. It then generates raw p-values by
#'   simulating the test statistics using a multivariate normal distribution
#'   with the given correlation matrix.
#'
#'   `simulate_pvalues()` assumes a point global alternative, using a vector of
#'   transformed p-values \eqn{(\Phi^{-1}(1-p_1), \ldots, \Phi^{-1}(1-p_m))}
#'   which follows a multivariate normal distribution with a correlation
#'   matrix \eqn{\Sigma}. Here, \eqn{\Phi^{-1}} signifies the inverse function
#'   of the standard normal distribution.
#'
#'   This assumption holds, for example, when \eqn{p_1, \ldots, p_m} are the raw
#'   p-values derived from one-sided z-tests for distinct hypotheses.
#'
#'   Note: applying the transformation
#'   \eqn{\Phi^{-1}(1-p_i)} to p-values from two-sided tests does not generally
#'   result in a multivariate normal distribution.
#'
#' @param power_nominal A numeric vector of nominal power values for each
#'   hypothesis.
#' @param alpha A numeric value indicating the significance level.
#' @param corr_matrix A numeric matrix representing the correlation matrix
#'   \eqn{\Sigma} of the test statistics.
#' @param nsim An integer indicating the number of simulations to run. Default
#'   is `1e5`.
#'
#' @return A matrix where each row represents a set of simulated raw p-values
#'   for each hypothesis.
#'
#' @export
#' @examples
#' # Define parameters for simulation
#' nominal_power <- c(0.8, 0.85, 0.9)
#' alpha_level <- 0.025
#' corr <- matrix(
#'   c(1, 0.5, 0.5,
#'     0.5, 1, 0.5,
#'     0.5, 0.5, 1),
#'    nrow = 3
#'  )
#' num_simulations <- 1000
#'
#' # Simulate raw p-values
#' pvals <- simulate_pvalues(nominal_power, alpha_level, corr, num_simulations)
simulate_pvalues <- function(
    power_nominal,
    alpha,
    corr_matrix = diag(length(power_nominal)),
    nsim = 1e5
) {
    check_double(power_nominal)
    rlang::check_number_decimal(alpha, min = 0, max = 1)
    check_double_matrix(corr_matrix)
    rlang::check_number_whole(nsim)

    # Calculate non-centrality parameter (ncp) using the calc_ncp function
    ncp <- calc_ncp(power_nominal, alpha)

    # Simulate raw p-values from the multivariate normal distribution
    pvals <- stats::pnorm(
        mvtnorm::rmvnorm(nsim, mean = ncp, sigma = corr_matrix),
        lower.tail = FALSE
    )

    pvals
}
