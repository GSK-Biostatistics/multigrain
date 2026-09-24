#' Create a parallelised group sequential objective function
#'
#' Group sequential twin of `create_obj_func()`. Constructs a closure that
#' captures everything needed for fitness evaluation, so that optimisers
#' (`GA::ga`, `nloptr::nloptr`) do not need to forward large objects like
#' `pvals` via `...`. The only differences from the fixed-sample factory are
#' the kernel it calls and the matrix it hands to the gain: the group
#' sequential kernel returns decision times as well as rejections, and a
#' `multigrain_trial_success_gsd` object scores the times.
#'
#' @param m Number of hypotheses.
#' @param power_criterion Function accessed from [trial_success_gsd()] used to
#'   calculate the trial success measure. It takes the kernel's integer matrix
#'   of decision times (`0` = never rejected).
#' @param hyp_constraint (numeric) Vector of hypothesis weight constraints.
#' @param trans_constraint (numeric matrix) Transition matrix constraints.
#' @param pvals (numeric matrix) The reshaped repeated p-values
#'   (`nsim` by `m * K`) returned by `.gsd_kernel_matrix()`. Column
#'   `(k - 1) * m + i` is hypothesis `i` at analysis `k`. The caller has
#'   already asserted that it holds no `NA`.
#' @param K (integer) Number of analyses.
#' @param alpha (numeric scalar) Overall one-sided significance level. Default
#'   is 0.025.
#' @inheritParams graph_optimise_gsd
#'
#' @returns A function with signature `function(x)` that evaluates the trial
#'   success measure for a given encoded parameter vector `x`. All other
#'   inputs (`alpha`, `pvals`, `K`, constraints, thread count) are captured in
#'   the closure at creation time.
#' @noRd
create_obj_func_gsd <- function(
    m,
    power_criterion,
    hyp_constraint,
    trans_constraint,
    pvals,
    K, # nolint: object_name_linter. K is the design record's symbol
    alpha = 0.025,
    num_threads = 1L
) {
    force(power_criterion)
    force(hyp_constraint)
    force(trans_constraint)
    force(alpha)
    force(pvals)
    force(K)
    force(num_threads)

    use_parallel <- num_threads >= 2L

    function(x) {
        theta <- split_theta(x, hyp_constraint)
        hyp_weight <- recover_full_weights(theta$w_pars, hyp_constraint)
        trans_matrix <- recover_full_trans_matrix(
            theta$g_pars,
            trans_constraint
        )

        if (anyNA(hyp_weight) || anyNA(trans_matrix)) {
            return(-1e6)
        }
        if (any(hyp_weight < 0)) {
            return(sum(hyp_weight[hyp_weight < 0]))
        }
        if (any(trans_matrix < 0)) {
            return(sum(trans_matrix[trans_matrix < 0]))
        }
        if (any(hyp_weight > 1)) {
            return(-sum(hyp_weight[hyp_weight > 1]))
        }
        if (any(trans_matrix > 1)) {
            return(-sum(trans_matrix[trans_matrix > 1]))
        }

        hyp_weight[hyp_weight < 1e-4] <- 0
        trans_matrix[trans_matrix < 1e-5] <- 0

        res <- if (use_parallel) {
            graph_shortcut_gsd_parallel(
                pvals = pvals,
                alpha = alpha,
                w = hyp_weight,
                G = trans_matrix,
                K = K,
                num_threads = num_threads,
                grain_size = -1L
            )
        } else {
            graph_shortcut_gsd(
                pvals = pvals,
                alpha = alpha,
                w = hyp_weight,
                G = trans_matrix,
                K = K
            )
        }

        power_criterion(res$time)
    }
}
