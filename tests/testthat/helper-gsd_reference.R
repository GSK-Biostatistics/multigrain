# A direct implementation of the group sequential graphical procedure of
# Maurer and Bretz (2013), with nominal boundaries computed on the fly at the
# hypothesis's current weight. It is the oracle for `graph_shortcut_gsd()`:
# where the kernel runs the fixed-sample cascade on repeated p-values, this
# compares each raw p-value with the boundary its current level buys.
#
# Ported from the `direct()` function and its helpers (`nom_bounds()`,
# `upd()`, `bound()`, `eff_l()`) in Appendix B of `dev/gsd_design_record.md`,
# with the semantics kept exactly as the record states them. Every hypothesis
# must reach full information at some analysis, as in the record's setup.

# Cumulative spend of a spending function, whether it returns a plain numeric
# vector or an object with a `spend` element.
gsd_reference_spend <- function(spending, level, info) {
    spent <- spending(level, info)
    if (is.list(spent) && !is.null(spent[["spend"]])) {
        spent <- spent[["spend"]]
    }

    as.numeric(spent)
}

# Forward map: nominal p-value boundaries at each analysis for a one-sided
# alpha-spending design holding level `level` at information fractions `info`.
gsd_reference_bounds <- function(level, info, spending) {
    increments <- diff(c(0, gsd_reference_spend(spending, level, info)))
    bounds <- gsDesign::gsBound1(
        theta = 0,
        I = info,
        a = rep(-20, length(info)),
        probhi = increments
    )$b

    stats::pnorm(bounds, lower.tail = FALSE)
}

# Reshape the transformed array of a `multigrain_pvals_gsd` object (or a bare
# nsim by m by K array) into the nsim by m*K matrix the kernel takes. Column
# (k - 1) * m + i is hypothesis i at analysis k, which R's column-major layout
# gives for free.
gsd_pvals_matrix <- function(x) {
    values <- if (inherits(x, "multigrain_pvals_gsd")) x$pvals else x
    dims <- dim(values)
    dim(values) <- c(dims[[1L]], dims[[2L]] * dims[[3L]])

    values
}

# The Bretz et al. (2009) graph update after rejecting hypothesis `i`. `rej`
# already has `i` set.
gsd_reference_update <- function(w, G, i, rej) {
    m <- length(w)
    w_new <- w + w[[i]] * G[i, ]
    w_new[[i]] <- 0
    G_new <- matrix(0, m, m)

    for (l in seq_len(m)) {
        for (j in seq_len(m)) {
            if (l != j && !rej[[l]] && !rej[[j]]) {
                denom <- 1 - G[l, i] * G[i, l]
                G_new[l, j] <- if (denom > 0) {
                    (G[l, j] + G[l, i] * G[i, j]) / denom
                } else {
                    0
                }
            }
        }
    }

    list(w = w_new, G = G_new)
}

# One trial of the procedure. `pk` is the m by K matrix of that trial's raw
# p-values; `bound(i, l, a)` returns the nominal boundary of hypothesis i at
# analysis l when it holds level a; `eff_l(i, k)` is the analysis whose
# evidence hypothesis i uses at analysis k (its own, or the one at which it
# matured).
gsd_reference_trial <- function(pk, w0, G0, alpha, look_back, bound, eff_l) {
    m <- length(w0)
    n_look <- ncol(pk)
    w <- w0
    G <- G0
    rej <- rep(FALSE, m)
    decision <- integer(m)

    for (k in seq_len(n_look)) {
        repeat {
            hit <- NA_integer_
            for (i in which(!rej)) {
                a <- w[[i]] * alpha
                if (a <= 0) {
                    next
                }
                look_idx <- if (look_back[[i]]) {
                    seq_len(eff_l(i, k))
                } else {
                    eff_l(i, k)
                }
                crossed <- pk[i, look_idx] <
                    vapply(look_idx, function(l) bound(i, l, a), double(1L))
                if (any(crossed)) {
                    hit <- i
                    break
                }
            }
            if (is.na(hit)) {
                break
            }
            rej[[hit]] <- TRUE
            decision[[hit]] <- k
            updated <- gsd_reference_update(w, G, hit, rej)
            w <- updated$w
            G <- updated$G
        }
    }

    list(rejected = rej, time = decision)
}

#' The Maurer and Bretz procedure with boundaries computed on the fly
#'
#' @param p_list A list of `m` raw p-value matrices, one per hypothesis, each
#'   `nsim` by `K`. A matured hypothesis must have its p-value carried forward.
#' @param tmat An `m` by `K` matrix of information fractions.
#' @param sfs A list of `m` spending functions.
#' @param w0,G0 The initial hypothesis weights and transition matrix.
#' @param alpha The overall one-sided level.
#' @param look_back A logical of length 1 or `m`.
#' @param cache An environment memoising boundary evaluations; pass the same
#'   one across configurations to avoid recomputing `gsBound1()` calls.
#'
#' @returns A list with `rejected` (`nsim` by `m` logical) and `time`
#'   (`nsim` by `m` integer, 0 = never rejected).
#' @noRd
gsd_reference_direct <- function(
    p_list,
    tmat,
    sfs,
    w0,
    G0,
    alpha = 0.025,
    look_back = FALSE,
    cache = new.env(parent = emptyenv())
) {
    m <- nrow(tmat)
    n_sim <- nrow(p_list[[1L]])
    if (length(look_back) == 1L) {
        look_back <- rep(look_back, m)
    }

    mature_k <- apply(tmat, 1L, function(t) which(t >= 1)[[1L]])
    eff_l <- function(i, k) min(k, mature_k[[i]])

    bound <- function(i, l, a) {
        key <- sprintf("%d_%d_%.15g", i, l, a)
        value <- cache[[key]]
        if (is.null(value)) {
            value <- gsd_reference_bounds(
                a,
                tmat[i, seq_len(l)],
                sfs[[i]]
            )[[l]]
            cache[[key]] <- value
        }

        value
    }

    rejected <- matrix(FALSE, nrow = n_sim, ncol = m)
    decision_time <- matrix(0L, nrow = n_sim, ncol = m)

    for (n in seq_len(n_sim)) {
        pk <- do.call(rbind, lapply(p_list, function(p) p[n, ]))
        trial <- gsd_reference_trial(
            pk,
            w0 = w0,
            G0 = G0,
            alpha = alpha,
            look_back = look_back,
            bound = bound,
            eff_l = eff_l
        )
        rejected[n, ] <- trial$rejected
        decision_time[n, ] <- trial$time
    }

    list(rejected = rejected, time = decision_time)
}
