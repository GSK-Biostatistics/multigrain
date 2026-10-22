# Split-merge (SM) mutation for the GA.
#
# SM works on the GA's own encoding: free weight parameters first, then each
# transition row's free parameters, as split_theta(), recover_full_weights()
# and recover_full_trans_matrix() decode them. A group (the free weights, or
# the free entries of one transition row) stores one parameter per free entry
# less one; the last free entry is derived so the group sums to its free mass
# c (one minus its fixed entries). The free entries are c * s, with s on the
# unit simplex.
#
# Each call makes one move in one group: with probability p_del a merge (the
# smaller of two active entries gives its mass to the larger), otherwise a
# transfer (a reflected Cauchy step in the split of two entries' total mass).
# SM works on decoded shares, so no entry is special by where it is stored, and
# every child is a valid graph that differs from its parent. It replaced the
# per-parameter truncated Cauchy mutation (archived in
# dev/archived_code/cauchy_mutation.R) after screen 2 of multigrain-benchmark,
# where p_del = 0.45, scale = 0.2 was chosen.
#
# Designed for use with pmutation = 1 (every individual mutated every
# generation), with exploitation handled by Nelder-Mead local search
# (optim = TRUE, poptim = 0.2) and elitism preserving the best solution(s).

# The variable groups of the latent vector, in the order the decoder reads it.
#
# Every group with two or more free entries takes up k - 1 latent positions.
# Only groups with free mass c > 0 can be mutated, so only those are returned.
# Thresholds are those below which the objective treats an entry as zero.
#
# @returns A list of groups, each list(latent, c, threshold).
.sm_groups <- function(hyp_constraint, trans_constraint) {
    groups <- list()
    pos <- 0L

    k_w <- sum(is.na(hyp_constraint))
    if (k_w >= 2L) {
        groups[[length(groups) + 1L]] <- list(
            latent = pos + seq_len(k_w - 1L),
            c = 1 - sum(hyp_constraint, na.rm = TRUE),
            threshold = 1e-4
        )
        pos <- pos + k_w - 1L
    }
    for (i in seq_len(nrow(trans_constraint))) {
        k <- sum(is.na(trans_constraint[i, ]))
        if (k >= 2L) {
            groups[[length(groups) + 1L]] <- list(
                latent = pos + seq_len(k - 1L),
                c = 1 - sum(trans_constraint[i, ], na.rm = TRUE),
                threshold = 1e-5
            )
            pos <- pos + k - 1L
        }
    }

    Filter(function(g) g$c > 0, groups)
}


# Euclidean projection onto the unit simplex.
.simplex_project <- function(v) {
    u <- sort(v, decreasing = TRUE)
    rho <- max(which(u - (cumsum(u) - 1) / seq_along(u) > 0))
    pmax(v - (sum(u[seq_len(rho)]) - 1) / rho, 0)
}


# A group's free entries, c * s, derived entry last. May be negative.
.sm_free_entries <- function(x, g) {
    p <- x[g$latent]
    c(p, g$c - sum(p))
}


# Decode a group to shares s on the simplex, projecting an invalid one.
.sm_to_shares <- function(x, g) {
    v <- .sm_free_entries(x, g) / g$c
    if (any(v < 0)) {
        return(.simplex_project(v))
    }
    v
}


# Encode shares s back to the group's latent parameters.
#
# Rounding can leave the derived entry a few ulps below zero, which the
# objective would reject as an invalid graph. The decoder sums the fixed
# entries too, in its own order, so keep the parameters a relative 1e-14
# inside c and take any excess off the largest parameter.
.sm_from_shares <- function(x, g, s) {
    p <- g$c * s[-length(s)]
    for (attempt in 1:3) {
        excess <- sum(p) - g$c * (1 - 1e-14)
        if (excess <= 0) {
            break
        }
        j <- which.max(p)
        p[j] <- max(p[j] - excess, 0)
    }
    x[g$latent] <- pmin(pmax(p, 0), 1)
    x
}


# Fold any real number into [0, 1]: v mod 2, then 2 - v above 1. The identity
# on [0, 1]. NaN for a non-finite v.
.reflect01 <- function(v) {
    v <- v %% 2
    ifelse(v > 1, 2 - v, v)
}


# Entries of s whose free mass c * s is at or above the group's threshold.
.sm_active <- function(g, s) {
    which(g$c * s >= g$threshold)
}


# Merge in one group: two distinct active entries, picked uniformly; the
# smaller's mass goes to the larger (the second on a tie).
.sm_merge <- function(g, s) {
    a <- .sm_active(g, s)
    pair <- a[sample.int(length(a), 2L)]
    small <- if (s[pair[1L]] < s[pair[2L]]) pair[1L] else pair[2L]
    large <- setdiff(pair, small)
    s[large] <- s[large] + s[small]
    s[small] <- 0
    s
}


# Transfer in one group. The donor is uniform among the active entries (the
# largest entry if none is active), the recipient uniform among the other free
# entries. The donor's share u of the pair's total moves by a reflected Cauchy
# step, redrawn until it lands strictly inside (0, 1), so both entries end
# positive.
.sm_transfer <- function(g, s, scale) {
    a <- .sm_active(g, s)
    i <- if (length(a)) a[sample.int(length(a), 1L)] else which.max(s)
    others <- setdiff(seq_along(s), i)
    j <- others[sample.int(length(others), 1L)]
    total <- s[i] + s[j]
    u <- s[i] / total
    repeat {
        u_new <- .reflect01(u + scale * stats::rcauchy(1))
        if (is.finite(u_new) && u_new > 0 && u_new < 1) {
            break
        }
    }
    s[i] <- total * u_new
    s[j] <- total - s[i]
    s
}


# Factory: returns a split-merge mutation closure.
#
# @param hyp_constraint,trans_constraint The graph constraint's weight and
#   transition constraints (NA marks a free entry).
# @param p_del Probability of a merge. Default 0.45. If no group has two
#   active entries a transfer is made instead.
# @param scale Scale of the Cauchy step in a transfer. Default 0.2.
# @returns A function(object, parent) for the [GA::ga()] mutation slot. If no
#   group can be mutated the closure returns the parent unchanged.
.make_split_merge_mutation <- function(
    hyp_constraint,
    trans_constraint,
    p_del = 0.45,
    scale = 0.2
) {
    stopifnot(
        is.numeric(p_del),
        length(p_del) == 1L,
        p_del >= 0,
        p_del <= 1,
        is.numeric(scale),
        length(scale) == 1L,
        scale > 0
    )
    force(p_del)
    force(scale)
    groups <- .sm_groups(hyp_constraint, trans_constraint)

    function(object, parent) {
        x <- as.numeric(object@population[parent, ])
        if (!length(groups)) {
            return(x)
        }
        shares <- lapply(groups, function(g) .sm_to_shares(x, g))

        if (stats::runif(1) < p_del) {
            eligible <- which(vapply(
                seq_along(groups),
                function(k) length(.sm_active(groups[[k]], shares[[k]])) >= 2L,
                logical(1)
            ))
            if (length(eligible)) {
                k <- eligible[sample.int(length(eligible), 1L)]
                return(.sm_from_shares(
                    x,
                    groups[[k]],
                    .sm_merge(groups[[k]], shares[[k]])
                ))
            }
        }

        # Transfer (also the fallback when no group can merge)
        k <- sample.int(length(groups), 1L)
        .sm_from_shares(
            x,
            groups[[k]],
            .sm_transfer(groups[[k]], shares[[k]], scale)
        )
    }
}
