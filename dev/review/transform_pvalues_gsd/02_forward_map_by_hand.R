# Exercise 02: the forward map, by hand
#
# Everything in the transform rests on one function: given an allocated level
# gamma, a spending function phi, and information fractions t, what are the
# nominal p-value boundaries alpha*_k(gamma) at each look k?
#
# This exercise computes that for two looks without gsDesign, using only
# qnorm() and a bivariate normal, and then shows gsBound1() gives the same.
#
# The model. Under H0 the z-statistics at looks 1 and 2 are jointly normal,
# mean 0, variance 1, correlation sqrt(t1 / t2) (the "canonical joint
# distribution" of group sequential theory; for t2 = 1 this is sqrt(t1)).
#
# The spending function says how much type I error may be used by each look:
#   phi(gamma, t1)                    at look 1
#   phi(gamma, t2) - phi(gamma, t1)   at look 2 (the increment)
#
# The boundaries c1, c2 (z-scale) solve
#   P(Z1 > c1)                 = phi(gamma, t1)
#   P(Z1 <= c1, Z2 > c2)       = phi(gamma, t2) - phi(gamma, t1)
# and the p-value boundaries are alpha*_k = 1 - pnorm(c_k).

suppressMessages(devtools::load_all(quiet = TRUE))
library(gsDesign)
library(mvtnorm)

gamma <- 0.025
t     <- c(0.5, 1)
phi   <- sfLDOF                       # Lan-DeMets O'Brien-Fleming

# --- 1. How much alpha does each look get? ----------------------------------

cum   <- phi(gamma, t)$spend           # cumulative spend at each look
inc   <- diff(c(0, cum))               # per-look increments
cat("cumulative spend:", cum, "\n")
cat("increments      :", inc, "\n")
cat("sum of increments == gamma?", isTRUE(all.equal(sum(inc), gamma)), "\n\n")

# >>> Question: LDOF spends very little early. What fraction of gamma is
#     available at the halfway look? What would sfLDPocock give? Try it.

# --- 2. Look 1 boundary: a one-dimensional problem --------------------------

c1 <- qnorm(1 - inc[1])
b1 <- 1 - pnorm(c1)                    # equals inc[1] by construction
cat("look 1: c1 =", c1, " boundary =", b1, "\n")

# --- 3. Look 2 boundary: a two-dimensional root-find ------------------------
#
# Find c2 such that P(Z1 <= c1, Z2 > c2) = inc[2].

rho   <- sqrt(t[1] / t[2])
Sigma <- matrix(c(1, rho, rho, 1), 2)

p_cross_at_2 <- function(c2) {
    pmvnorm(lower = c(-Inf, c2), upper = c(c1, Inf), sigma = Sigma)[1]
}
c2 <- uniroot(function(z) p_cross_at_2(z) - inc[2], c(0, 6), tol = 1e-12)$root
b2 <- 1 - pnorm(c2)
cat("look 2: c2 =", c2, " boundary =", b2, "\n")
cat("        rho =", rho, "\n\n")

# >>> Question: b2 is below gamma (0.025) even though all the alpha has been
#     spent by look 2. Why? (Hint: some of it was spent at look 1, and the
#     two looks are correlated.) What happens to b2 as t1 -> 0?

# --- 4. Compare with gsBound1 ----------------------------------------------

gs <- gsBound1(theta = 0, I = t, a = rep(-20, 2), probhi = inc)
b_gs <- pnorm(gs$b, lower.tail = FALSE)
cat("by hand  :", format(c(b1, b2), digits = 10), "\n")
cat("gsBound1 :", format(b_gs, digits = 10), "\n")
cat("max abs diff:", max(abs(c(b1, b2) - b_gs)), "\n\n")

# `a = rep(-20, K)` is a lower boundary so far below zero that it never
# binds; gsBound1 always wants one. probhi are the crossing probabilities.

# --- 5. The boundary is monotone in gamma ----------------------------------
#
# This is the property the whole transform depends on. Recompute for a few
# levels and watch each look's boundary rise with gamma.

nom_bounds <- function(gamma, t, sf = sfLDOF) {
    inc <- diff(c(0, sf(gamma, t)$spend))
    b <- gsBound1(theta = 0, I = t, a = rep(-20, length(t)), probhi = inc)$b
    pnorm(b, lower.tail = FALSE)
}

levels <- c(0.001, 0.005, 0.0125, 0.025)
tab <- t(sapply(levels, nom_bounds, t = t))
dimnames(tab) <- list(paste("gamma =", levels), c("look 1", "look 2"))
print(signif(tab, 4))

# >>> Question: pick a raw p-value at look 2, say p = 0.02. Reading down the
#     "look 2" column, at which gamma does the boundary first exceed 0.02?
#     That gamma is (approximately) the repeated p-value. Exercise 04 makes
#     this precise.

stopifnot(
    isTRUE(all.equal(sum(inc), gamma)),
    max(abs(c(b1, b2) - b_gs)) < 1e-8,
    all(diff(tab[, 1]) > 0),
    all(diff(tab[, 2]) > 0)
)
cat("\nOK\n")
