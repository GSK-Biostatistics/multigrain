# Exercise 07: look-back (sequential p-values)
#
# Default semantics (look_back = FALSE): a hypothesis is rejected at look k
# only on the strength of its look-k evidence. That is Algorithm 1 of Maurer
# and Bretz: at each analysis, compare p_{i,k} with the boundary at the
# *current* allocation.
#
# Look-back (look_back = TRUE): if a hypothesis gains alpha at a later look
# (because something else was rejected and recycled), it may be rejected on
# the strength of its evidence at ANY look so far. In repeated-p currency
# that is the running minimum p^s_{i,k} = min(p^r_{i,1}, ..., p^r_{i,k}),
# the "sequential p-value". No new machinery: one cummin per row.

suppressMessages(devtools::load_all(quiet = TRUE))
library(gsDesign)

# --- 1. It is literally a running minimum --------------------------------

set.seed(5)
t   <- c(1/3, 2/3, 1)
raw <- simulate_pvalues_gsd(c(0.8, 0.8), info_frac = t, nsim = 6)

rep_p <- transform_pvalues_gsd(raw, spending = sfLDOF)$pvals
seq_p <- transform_pvalues_gsd(raw, spending = sfLDOF, look_back = TRUE)$pvals

cat("H1, repeated p (rows = simulations, cols = looks):\n")
print(signif(rep_p[, 1, ], 4))
cat("\nH1, sequential p (look_back = TRUE):\n")
print(signif(seq_p[, 1, ], 4))

by_hand <- t(apply(rep_p[, 1, ], 1, cummin))
cat("\nequal to cummin along each row?", isTRUE(all.equal(seq_p[, 1, ], by_hand)), "\n\n")

# >>> Question: sequential p-values are non-increasing in k. What does that
#     mean for decision times under look-back versus without?

# --- 2. A case where it changes the decision --------------------------------
#
# Build one trial by hand. Two hypotheses, alpha split 0.0125 each, with an
# edge H2 -> H1 (H2's alpha goes to H1 when H2 is rejected).
#
#   H1: strong evidence at look 2, none to speak of at look 3.
#   H2: rejected only at look 3.
#
# Without look-back, H1 at look 3 is judged on its look-3 evidence, which
# is not good enough even at the full 0.025. With look-back, H1 is judged on
# the best it has shown so far (look 2), which IS good enough at 0.025.

alpha <- 0.025
w     <- c(0.5, 0.5)

one <- array(1, dim = c(1, 2, 3))
one[1, 1, ] <- c(0.05, 0.003, 0.03)     # H1 raw p at looks 1..3
one[1, 2, ] <- c(0.20, 0.10,  0.005)    # H2 raw p at looks 1..3

pr <- transform_pvalues_gsd(one, info_frac = t, spending = sfLDOF)$pvals[1, , ]
ps <- transform_pvalues_gsd(one, info_frac = t, spending = sfLDOF,
                            look_back = TRUE)$pvals[1, , ]
dimnames(pr) <- dimnames(ps) <- list(c("H1", "H2"), paste0("look", 1:3))
cat("repeated p:\n");   print(signif(pr, 4))
cat("\nsequential p:\n"); print(signif(ps, 4))

# A tiny two-hypothesis cascade, run look by look, carrying the graph state.
run_looks <- function(P) {
    w   <- c(0.5, 0.5)
    rej <- c(FALSE, FALSE)
    when <- c(0L, 0L)
    for (k in 1:3) {
        repeat {
            j <- which(!rej & P[, k] <= w * alpha)[1]
            if (is.na(j)) break
            rej[j] <- TRUE; when[j] <- k
            other <- 3 - j
            if (!rej[other]) w[other] <- w[other] + w[j]   # single edge each way
            w[j] <- 0
        }
    }
    list(rejected = rej, at_look = when)
}

cat("\nwithout look-back:\n"); print(run_looks(pr))
cat("\nwith look-back:\n");    print(run_looks(ps))

# >>> Question: trace it. At look 2, H1's repeated p is ~0.015 > 0.0125: not
#     rejected (close, but H1 only holds half the alpha). At look 3, H2 is
#     rejected (~0.005 < 0.0125) and passes its alpha to H1, which now holds
#     0.025. Without look-back H1 is judged on its look-3 repeated p, which
#     is 1: raw 0.03 is above the look-3 boundary even at the full alpha
#     (about 0.023), so no allocation could reject it. With look-back H1 is
#     judged on its sequential p, min(1, 0.015, 1) = 0.015 < 0.025: rejected.
#
#     Same data, same graph, different answer. The protocol has to say which
#     rule the trial uses.

# --- 3. Why is it off by default? -----------------------------------------
#
# Both semantics control the FWER (design record 4.2 and 4.9). Look-back is
# uniformly more powerful, but Algorithm 1 of Maurer and Bretz (the
# established, regulator-familiar version) does not do it. The protocol has
# to say which one the trial uses; the package makes it a per-hypothesis
# switch and defaults to the conservative choice.

no_lb <- run_looks(pr); lb <- run_looks(ps)
stopifnot(
    isTRUE(all.equal(seq_p[, 1, ], by_hand)),
    all(seq_p <= rep_p),
    identical(no_lb$rejected, c(FALSE, TRUE)),
    identical(lb$rejected, c(TRUE, TRUE))
)
cat("\nOK\n")
