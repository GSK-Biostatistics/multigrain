# Exercise 08: why it works
#
# Two checks. The first is the theorem the whole design rests on; the second
# is the consequence a regulator cares about.
#
# THEOREM (Maurer and Bretz 2013, section 3.3). If the boundary alpha*_k(a)
# is non-decreasing in a (well ordered), then for every p and every a,
#
#        p <= alpha*_k(a)     <=>     p^r_k <= a
#
# where p^r_k is the repeated p-value, the smallest a' with alpha*_k(a') >= p.
# Proof in one line: alpha*_k is monotone, so {a : alpha*_k(a) >= p} is an
# upper set [p^r_k, alpha], and "a is in that set" is "a >= p^r_k".
#
# Consequence: the graphical procedure never needs a boundary. Whatever
# allocation w_i * alpha a hypothesis holds at the moment it is tested, the
# decision "p_{i,k} <= alpha*_{i,k}(w_i alpha)" is the decision
# "p^r_{i,k} <= w_i alpha", and the latter is a fixed-sample comparison.

suppressMessages(devtools::load_all(quiet = TRUE))
library(gsDesign)

alpha <- 0.025
t     <- c(1/3, 2/3, 1)

nom_bounds <- function(gamma, t, sf = sfLDOF) {
    inc <- diff(c(0, sf(gamma, t)$spend))
    b <- gsBound1(theta = 0, I = t, a = rep(-20, length(t)), probhi = inc)$b
    pnorm(b, lower.tail = FALSE)
}

# --- 1. The equivalence, checked on random (p, allocation) pairs -----------
#
# For each random p at each look and each random allocation a in (0, alpha]:
# decide both ways and count disagreements. Draw p on a log scale so the
# interesting region (tiny p) is well covered.

set.seed(6)
n_p <- 300
n_a <- 20
p_draw <- 10^runif(n_p, -7, log10(0.3))
a_draw <- 10^runif(n_a, -6, log10(alpha))

raw <- array(1, dim = c(n_p, 1, 3))
for (k in 1:3) raw[, 1, k] <- p_draw
pr <- transform_pvalues_gsd(raw, info_frac = t, spending = sfLDOF)$pvals[, 1, ]

# boundaries at each drawn allocation: n_a x K
B_at_a <- t(sapply(a_draw, nom_bounds, t = t))

disagree <- 0L; total <- 0L; close_calls <- 0L
for (ia in seq_len(n_a)) {
    for (k in 1:3) {
        direct   <- p_draw <= B_at_a[ia, k]
        via_rep  <- pr[, k] <= a_draw[ia]
        d <- direct != via_rep
        disagree <- disagree + sum(d)
        total    <- total + length(d)
        # how many decisions sit within interpolation error of the threshold?
        close_calls <- close_calls + sum(abs(pr[, k] / a_draw[ia] - 1) < 1e-5)
    }
}
cat(sprintf("decisions compared: %d\n", total))
cat(sprintf("disagreements     : %d\n", disagree))
cat(sprintf("within 1e-5 of the threshold (the only place interpolation could bite): %d\n\n",
            close_calls))

# >>> Question: if there were a disagreement, where could it come from?
#     (Only from interpolation error, and only when p^r is within ~1e-6
#     relative of the allocation. Design record 4.1 argues this affects
#     about one simulated row in a million.)

# --- 2. Type I error under the global null ---------------------------------
#
# If the transform is right, then rejecting H at look k whenever
# p^r_k <= alpha (the whole level, single hypothesis, no graph) is exactly
# the alpha-spending test, whose type I error is alpha by construction.
# Simulate under H0 and check the *any-look* rejection rate.

set.seed(7)
nsim <- 2e5
# power_nominal = alpha means "no effect": the noncentrality is zero.
null <- simulate_pvalues_gsd(alpha, info_frac = t, nsim = nsim)
prn  <- transform_pvalues_gsd(null, spending = sfLDOF)$pvals[, 1, ]

reject_any <- apply(prn <= alpha, 1, any)
rate <- mean(reject_any)
se   <- sqrt(alpha * (1 - alpha) / nsim)
cat(sprintf("nominal alpha        : %.4f\n", alpha))
cat(sprintf("simulated FWER       : %.4f  (MC s.e. %.4f)\n", rate, se))
cat(sprintf("rejections by look   : %s\n",
            paste(sprintf("%.4f", colMeans(prn <= alpha & !t(apply(cbind(FALSE, prn[, -3] <= alpha), 1, cummax)))), collapse = ", ")))
cat("(the by-look rates should sum to the FWER and match the LDOF increments)\n")
cat(sprintf("LDOF increments      : %s\n\n",
            paste(sprintf("%.4f", diff(c(0, sfLDOF(alpha, t)$spend))), collapse = ", ")))

# --- 3. And with the naive (wrong) approach ----------------------------------
#
# Compare each raw p with alpha at every look. This is what the transform
# saves you from: repeated testing at the same level.

naive <- mean(apply(null[, 1, ] <= alpha, 1, any))
cat(sprintf("naive 'test at alpha every look' FWER: %.4f  (inflated)\n", naive))

# >>> Question: about how much bigger than alpha is the naive rate with
#     three looks? Why is it not simply 3 * alpha? (The looks are correlated:
#     correlation sqrt(t_k / t_l).)

stopifnot(
    disagree == 0L,
    abs(rate - alpha) < 4 * se,
    naive > alpha + 4 * se
)
cat("\nOK\n")
