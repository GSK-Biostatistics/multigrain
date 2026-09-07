# Truncated Cauchy perturbation centered on parent value.
# Draws from Cauchy(centre, scale) via rejection sampling.
.cauchy_perturb <- function(centre, lower, upper, scale = 1.0) {
    for (i in seq_len(100L)) {
        proposal <- centre + scale * tan(pi * (stats::runif(1) - 0.5))
        if (proposal >= lower && proposal <= upper) {
            return(proposal)
        }
    }
    # Fallback: uniform on [lower, upper]
    stats::runif(1, lower, upper)
}


# Row index of each free transition parameter in the encoding (NA for
# hypothesis-weight parameters). Mirrors the parameter order used by
# recover_full_trans_matrix(): the weight parameters come first, then each row
# of G contributes k - 1 parameters for its k free entries, rows in order.
#
# @param hyp_constraint (numeric) Vector of hypothesis weight constraints.
# @param trans_constraint (numeric matrix) Transition matrix constraints.
# @returns An integer vector as long as the encoded parameter vector.
.g_param_rows <- function(hyp_constraint, trans_constraint) {
    n_w <- max(sum(is.na(hyp_constraint)) - 1L, 0L)
    k <- rowSums(is.na(trans_constraint))
    rows <- rep.int(seq_along(k), pmax(k - 1L, 0L))
    c(rep(NA_integer_, n_w), rows)
}


# Target sum for a row's free parameters under the derived-entry zeroing move.
#
# The derived entry of row i is 1 - sum(G[i, ]), i.e. one minus the row's
# fixed entries and its free parameters. Making it 5e-6 therefore means
# scaling the parameters to sum to 1 - 5e-6 - (fixed entries of that row).
#
# @param trans_constraint (numeric matrix) Transition matrix constraints.
# @returns A numeric vector with one target per row of `trans_constraint`.
.zeroing_row_targets <- function(trans_constraint) {
    1 - 5e-6 - rowSums(trans_constraint, na.rm = TRUE)
}


# Factory: returns a multi-parameter Cauchy perturbation mutation closure.
#
# Perturbs each parameter within graph with probability p_param_mutate
#
# Designed for use with p_mutation = 1.0 (every individual mutated every
# generation), with exploitation handled by Nelder-Mead local search
# (optim = TRUE, poptim = 0.2) and elitism preserving the best solution(s).
#
# --- What the mutation does today -------------------------------------------
#
# The closure the GA calls on one parent per mutation selects each parameter
# with probability 0.1 (at least one), and replaces each selected value with a
# draw from a Cauchy distribution centred on the current value and truncated to
# [0, 1]. That is a good operator for moving weight around continuously. It is
# a poor operator for removing an edge, because "removed" means the decoded
# entry falls below 1e-5, and a Cauchy draw lands in [0, 1e-5) about once in
# 100,000 tries. Crossover does not help: GA's default real-valued crossover
# takes convex combinations of two parents, so a child entry is zero only where
# both parents are already zero. Nelder-Mead and COBYLA are continuous methods
# and never aim for the threshold. With the lexicographic score of `.lexico()`
# the GA would therefore *prefer* sparser graphs but almost never *produce* one
# to compare, and the score would be inert during the global search apart from
# the seeds.
#
# --- What is added ----------------------------------------------------------
#
# In `graph_simplify()`, each mutation call first draws a coin with probability
# `p_zero`. If it comes up, the call performs a zeroing move instead of a
# Cauchy perturbation, of one of two kinds:
#
#   1. Zero a parameter. Pick, uniformly, one free transition parameter
#      currently at or above 1e-5 and set it to exactly zero. That removes the
#      corresponding edge outright.
#   2. Zero a derived entry. Each row of the transition matrix has one entry
#      that is not a parameter: it is computed as one minus the sum of the
#      row's parameters. Setting a parameter to zero cannot remove that edge.
#      The move instead picks a row and rescales its parameters so they sum to
#      1 - 5e-6 (less any entries the constraint pins in that row). The derived
#      entry then equals 5e-6: positive, so it is not penalised as a negative
#      entry, and below 1e-5, so it does not count as an edge and is snapped to
#      exact zero on output.
#
# If the coin does not come up, or there is nothing to zero (no parameter above
# the threshold, or a row whose parameters are all zero), the call falls
# through to the ordinary Cauchy perturbation.
#
# With `p_zero = 0` -- what `graph_optimise()` asks for -- the factory returns
# the Cauchy closure itself, unwrapped, so `graph_optimise()` runs the same
# code and consumes the same random numbers as before this move existed.
#
# @param p_param_mutate Per-parameter mutation probability. Default 0.1.
# @param scale Cauchy scale parameter. Default 1.0 (matching ESCH).
# @param p_zero Probability that a call performs a zeroing move rather than a
#   Cauchy perturbation. Default 0, which disables the move entirely.
# @param param_rows Row map from `.g_param_rows()`. Required when `p_zero` > 0.
# @param row_target Per-row parameter-sum targets from
#   `.zeroing_row_targets()`. Required when `p_zero` > 0.
# @returns A function(object, parent) for the [GA::ga()] mutation slot.
.make_cauchy_mutation_multi <- function(
    p_param_mutate = 0.1,
    scale = 1.0,
    p_zero = 0,
    param_rows = NULL,
    row_target = NULL
) {
    force(p_param_mutate)
    force(scale)
    force(p_zero)
    force(param_rows)
    force(row_target)

    cauchy_mutation <- function(object, parent) {
        parent_vec <- as.numeric(object@population[parent, ])
        d <- length(parent_vec)
        lower <- object@lower
        upper <- object@upper

        # Independently decide which parameters to mutate
        mutate_mask <- stats::runif(d) < p_param_mutate

        # Guarantee at least one mutation (avoid wasting a fitness evaluation)
        if (!any(mutate_mask)) {
            mutate_mask[sample.int(d, 1L)] <- TRUE
        }

        # Perturb selected parameters around their current values
        for (j in which(mutate_mask)) {
            parent_vec[j] <- .cauchy_perturb(
                parent_vec[j],
                lower[j],
                upper[j],
                scale
            )
        }

        parent_vec
    }

    if (p_zero <= 0) {
        return(cauchy_mutation)
    }

    g_idx <- which(!is.na(param_rows))
    rows <- unique(param_rows[g_idx])

    function(object, parent) {
        if (stats::runif(1) >= p_zero) {
            return(cauchy_mutation(object, parent))
        }

        parent_vec <- as.numeric(object@population[parent, ])

        if (stats::runif(1) < 0.5) {
            # 1. Zero a parameter, removing that edge outright.
            candidates <- g_idx[parent_vec[g_idx] >= 1e-5]
            if (length(candidates) == 0L) {
                return(cauchy_mutation(object, parent))
            }
            j <- candidates[sample.int(length(candidates), 1L)]
            parent_vec[j] <- 0
        } else {
            # 2. Zero a row's derived entry by rescaling its parameters.
            i <- rows[sample.int(length(rows), 1L)]
            idx <- g_idx[param_rows[g_idx] == i]
            s <- sum(parent_vec[idx])
            if (s <= 0 || row_target[i] <= 0) {
                return(cauchy_mutation(object, parent))
            }
            parent_vec[idx] <- pmin(parent_vec[idx] * (row_target[i] / s), 1)
        }

        parent_vec
    }
}
