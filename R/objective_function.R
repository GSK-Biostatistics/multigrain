#' Create a parallelised objective function used for optimisation
#'
#' Constructs a closure that captures all data needed for fitness evaluation.
#' The returned function has signature `function(x)` — it takes only the
#' encoded parameter vector — so that optimisers (`GA::ga`, `nloptr::nloptr`)
#' do not need to forward large objects like `pvals` via `...`.
#'
#' @param m Number of hypotheses.
#' @param power_criterion Function accessed from [trial_success()] used to
#'   calculate trial success measure.
#' @param hyp_constraint (numeric) Vector of hypothesis weight constraints.
#' @param trans_constraint (numeric matrix) Transition matrix constraints.
#' @param alpha (numeric scalar) Overall one-sided significance level.
#' @param pvals (numeric matrix) Matrix of p-values (`nsim × m`). Each row
#'   is a simulated trial; each column corresponds to a hypothesis.
#' @inheritParams graph_optimise
#'
#' @return A function with signature `function(x)` that evaluates the trial
#'   success measure for a given encoded parameter vector `x`. All other
#'   inputs (`alpha`, `pvals`, constraints, thread count) are captured in the
#'   closure at creation time.
#' @noRd
create_obj_func <- function(
    m,
    power_criterion,
    hyp_constraint,
    trans_constraint,
    alpha,
    pvals,
    num_threads = 1L
) {
    force(power_criterion)
    force(hyp_constraint)
    force(trans_constraint)
    force(alpha)
    force(pvals)
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

        rej_matrix <- if (use_parallel) {
            graph_shortcut_parallel(
                pvals = pvals,
                alpha = alpha,
                w = hyp_weight,
                G = trans_matrix,
                num_threads = num_threads,
                grain_size = -1L
            )
        } else {
            graph_shortcut(
                pvals = pvals,
                alpha = alpha,
                w = hyp_weight,
                G = trans_matrix
            )
        }

        power_criterion(rej_matrix)
    }
}

#' Split encoded parameter vector into hypothesis weights, and transition
#' matrix.
#'
#' @param theta (numeric) Stacked sample-size-optimisation vector.
#' @param hyp_constraint (numeric) Vector; NA marks free weights, otherwise
#'     fixed.
#'
#' @return list(w_pars = numeric(), g_pars = numeric()).
#' @noRd
split_theta <- function(theta, hyp_constraint) {
    stopifnot(is.numeric(theta), length(theta) >= 1L)

    free_w <- sum(is.na(hyp_constraint))
    w_len <- max(free_w - 1L, 0L)

    if (w_len) {
        w_pars <- theta[seq_len(w_len)]
        g_pars <- theta[-seq_len(w_len)]
    } else {
        w_pars <- numeric(0)
        g_pars <- theta
    }

    list(w_pars = w_pars, g_pars = g_pars)
}


#' Recover full weights from starting values or solution
#'
#' @param x A vector containing the free hypothesis weight parameters being
#'   optimised by [GA::ga()] or [nloptr::nloptr()].
#' @param hyp_constraint A vector containing constraints on hypothesis weights.
#'
#' @return A numeric vector with the same length as `hyp_constraint`.
#' @noRd
recover_full_weights <- function(x, hyp_constraint) {
    hyp_constraint_full <- numeric(length = length(hyp_constraint))
    free_vars_indices <- which(is.na(hyp_constraint))
    hyp_constraint_full[!is.na(hyp_constraint)] <- hyp_constraint[
        !is.na(hyp_constraint)
    ]

    hyp_constraint_full[free_vars_indices] <- c(x, NA)

    hyp_constraint_full[free_vars_indices[length(free_vars_indices)]] <-
        1 - sum(hyp_constraint_full, na.rm = TRUE)

    hyp_constraint_full
}

#' Recover full transition matrix from starting values or solution
#'
#' @param x A vector containing the free transition matrix weight parameters
#'   being optimised by [GA::ga()] or [nloptr::nloptr()].
#' @param trans_constraint A vector containing constraints on transition matrix.
#'
#' @return A matrix with the same dimensions as `trans_constraint`.
#' @noRd
recover_full_trans_matrix <- function(x, trans_constraint) {
    stopifnot(is.matrix(trans_constraint))
    m <- nrow(trans_constraint)
    n <- ncol(trans_constraint)

    G <- matrix(0, nrow = m, ncol = n)

    # Identify free elements
    free_idx_list <- apply(
        trans_constraint,
        1,
        function(x) which(is.na(x)),
        simplify = FALSE
    )
    fixed_idx_list <- apply(
        trans_constraint,
        1,
        function(x) which(!is.na(x)),
        simplify = FALSE
    )

    # How many params do we expect (k-1 per row w/ k free)
    n_params_needed <- sum(pmax(lengths(free_idx_list) - 1L, 0L))
    if (length(x) != n_params_needed) {
        stop(
            "Length of 'x' (",
            length(x),
            ") does not match number of required G parameters (",
            n_params_needed,
            ").",
            call. = FALSE
        )
    }

    x_pos <- 1L

    for (i in seq_len(m)) {
        fixed_cols <- fixed_idx_list[[i]]
        if (length(fixed_cols)) {
            G[i, fixed_cols] <- trans_constraint[i, fixed_cols]
        }

        free_cols <- free_idx_list[[i]]
        k <- length(free_cols)

        if (k == 0L) {
            next
        } else if (k == 1L) {
            G[i, free_cols] <- 1 - sum(G[i, ])
        } else {
            take <- k - 1L
            G[i, free_cols[seq_len(take)]] <- x[x_pos:(x_pos + take - 1L)]
            x_pos <- x_pos + take

            G[i, free_cols[k]] <- 1 - sum(G[i, ])
        }
    }

    G
}
