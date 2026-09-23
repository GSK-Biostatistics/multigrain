# Exercise 05: the whole transform in 30 lines
#
# Put exercises 02-04 together into a hand-rolled transform for the simple
# case (every hypothesis has data at every look, all mature at the last
# look, no look-back), and check it against transform_pvalues_gsd() and
# against Maurer and Bretz (2013) Table 2.

suppressMessages(devtools::load_all(quiet = TRUE))
library(gsDesign)

# --- 1. The hand-rolled version --------------------------------------------

nom_bounds <- function(gamma, t, sf) {
    inc <- diff(c(0, sf(gamma, t)$spend))
    b <- gsBound1(theta = 0, I = t, a = rep(-20, length(t)), probhi = inc)$b
    pnorm(b, lower.tail = FALSE)
}

my_transform <- function(pvals, t, sf, alpha = 0.025, G = 1024) {
    grid <- exp(seq(log(1e-14), log(alpha), length.out = G))
    B    <- t(sapply(grid, nom_bounds, t = t, sf = sf))     # G x K
    out  <- pvals
    for (k in seq_along(t)) {
        col  <- B[, k]
        keep <- !duplicated(col) & col > 0
        p    <- pvals[, , k]
        r    <- exp(approx(log(col[keep]), log(grid[keep]), xout = log(p), rule = 1)$y)
        r[p >= max(col)]       <- 1
        r[p < min(col[keep])]  <- 1e-14
        out[, , k] <- r
    }
    out
}

# --- 2. Compare with the package on random data ----------------------------

set.seed(3)
t   <- c(1/3, 2/3, 1)
raw <- simulate_pvalues_gsd(c(0.8, 0.9), info_frac = t, nsim = 200)

mine <- my_transform(raw, t, sfLDOF)
pkg  <- transform_pvalues_gsd(raw, spending = sfLDOF)$pvals

cat("max abs difference, hand vs package:", max(abs(mine - pkg)), "\n")
cat("max rel difference (where < 1)     :",
    max(abs(mine - pkg) / pkg, na.rm = TRUE), "\n\n")

# >>> Question: the differences are ~1e-16, i.e. floating-point only. What
#     did the package do that my_transform() did not? (Look at the callflow
#     document: argument checking, spending check, per-hypothesis looks,
#     maturity handling, look-back. None of them bite in this simple case.)

# --- 3. Maurer and Bretz (2013) Table 2 ------------------------------------
#
# Their worked example: four hypotheses, LDOF, three looks at 1/3, 2/3, 1.
# Raw p-values at looks 1 and 2, and the repeated p-values they report
# (computed with ADDPLAN).

p_look1 <- c(0.0062, 0.017, 0.009, 0.13)
p_look2 <- c(0.0002, 0.0035, 0.002, 0.06)
paper1  <- c(0.1141, 0.1683, 0.1316, 0.382)
paper2  <- c(0.0024, 0.0172, 0.0117, 0.1285)

# One "simulation", four hypotheses, three looks (look 3 filled with 1s).
mb <- array(1, dim = c(1, 4, 3))
mb[1, , 1] <- p_look1
mb[1, , 2] <- p_look2

tr <- transform_pvalues_gsd(mb, info_frac = t, spending = sfLDOF)
print(data.frame(
    hyp = paste0("H", 1:4),
    raw1 = p_look1, rep1 = signif(tr$pvals[1, , 1], 4), paper1 = paper1,
    raw2 = p_look2, rep2 = signif(tr$pvals[1, , 2], 4), paper2 = paper2
))

# >>> Question: every look-1 value, and H4 at look 2, came back as 1, while
#     the paper reports 0.114, 0.168, ... What happened?
#
#     The table only runs up to alpha = 0.025. A repeated p-value above alpha
#     means "no allocation the graph could ever give this hypothesis would
#     reject it", and the transform reports that as 1. The paper's numbers
#     are the mathematically exact inverse over (0, 1); the package's are the
#     same thing truncated at the only threshold that matters. For the
#     decision they are equivalent: 0.114 > 0.025 and 1 > 0.025 both mean
#     "not rejected".
#
#     To reproduce the paper's numbers, widen the grid by passing a larger
#     alpha. (This is only for checking; in a real design alpha = 0.025.)

tr_wide <- transform_pvalues_gsd(mb, info_frac = t, spending = sfLDOF, alpha = 0.5)
print(data.frame(
    hyp = paste0("H", 1:4),
    rep1 = signif(tr_wide$pvals[1, , 1], 4), paper1 = paper1,
    rep2 = signif(tr_wide$pvals[1, , 2], 4), paper2 = paper2
))

# >>> Question: H1 at look 1 has raw p = 0.0062 but repeated p = 0.114. In
#     words: at a third of the information, under O'Brien-Fleming spending,
#     a p-value of 0.006 would only be significant if H1 held an allocation
#     of 0.114, far more than the 0.0125 it has. So it is not rejected at
#     look 1. By look 2, raw 0.0002 -> repeated 0.0024 < 0.0125: rejected.
#
#     The residual differences (up to ~1e-4) are ADDPLAN versus gsDesign
#     numerics, not interpolation; the design record notes the same gap.

# --- 4. What the graph then does with them ----------------------------------
#
# The point of the transform: the graphical procedure can now be run on the
# repeated p-values exactly as it would on fixed-sample p-values, comparing
# each against w_i * alpha. Appendix A of the design record does this and
# reproduces the paper's decision (H1, H2, H3 rejected at look 2; H4 not).
# Nothing about looks, spending, or boundaries appears in that cascade.

stopifnot(
    max(abs(mine - pkg)) < 1e-12,
    all(tr$pvals[1, , 1] == 1),
    max(abs(tr_wide$pvals[1, , 1] - paper1)) < 2e-4,
    max(abs(tr_wide$pvals[1, , 2] - paper2)) < 2e-4
)
cat("\nOK\n")
