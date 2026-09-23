# Exercise 01: one look at full information
#
# Claim: if a hypothesis is analysed once, with all its information, the group
# sequential boundary at allocated level a is a itself. There is nothing to
# invert: the repeated p-value equals the raw p-value.
#
# What to take away: this is the K = 1 case that must be bit-identical to the
# fixed-sample package, and it is why `.gsd_transform_hyp()` has a short-circuit
# rather than building a one-column table.

suppressMessages(devtools::load_all(quiet = TRUE))
library(gsDesign)

# --- 1. Ask gsBound1 for the boundary at a single look --------------------
#
# gsBound1(theta = 0, I = t, a = lower, probhi = increments) returns z-scale
# upper boundaries b such that P(cross at look k | no earlier crossing) equals
# probhi[k] under theta = 0. With one look and probhi = a, b = qnorm(1 - a).

a <- 0.0125
b <- gsBound1(theta = 0, I = 1, a = -20, probhi = a)$b
boundary <- pnorm(b, lower.tail = FALSE)

cat("level a          :", format(a, digits = 17), "\n")
cat("qnorm(1 - a)     :", format(qnorm(1 - a), digits = 17), "\n")
cat("gsBound1 b       :", format(b, digits = 17), "\n")
cat("boundary pnorm(b):", format(boundary, digits = 17), "\n")
cat("identical(boundary, a)?", identical(boundary, a), "\n")
cat("difference       :", boundary - a, "\n\n")

# >>> Question: the boundary is a "to 15 digits" but not identical. If the
#     transform built a table here and inverted it, would a raw p-value of
#     exactly 0.0125 come back as exactly 0.0125? What would that do to a
#     test that compares the K = 1 path against the fixed-sample kernel?

# --- 2. What the package does ------------------------------------------
#
# Build a 1-look array, transform it, and compare to the raw values.

set.seed(1)
raw <- array(runif(6), dim = c(3, 2, 1))          # nsim = 3, m = 2, K = 1
attr(raw, "info_frac") <- matrix(1, 2, 1)

tr <- transform_pvalues_gsd(raw, spending = sfLDOF)

cat("raw p-values:\n"); print(raw[, , 1])
cat("repeated p-values:\n"); print(tr$pvals[, , 1])
cat("identical?", identical(raw[, , 1], tr$pvals[, , 1]), "\n")
cat("table built for H1?", !is.null(tr$tables[[1]]$bounds), "\n\n")

# --- 3. The same short-circuit fires mid-array -----------------------------
#
# A hypothesis with data only at look 2, at full information, is also a
# single matured look. Its look-1 repeated p-value is 1 (no data); its look-2
# value is the raw value; and a table is still not built.

raw2 <- array(runif(4), dim = c(2, 1, 2))          # nsim = 2, m = 1, K = 2
info <- matrix(c(NA, 1), nrow = 1)                 # 1 hypothesis x 2 looks: no data
                                                   # at look 1, all of it at look 2
tr2 <- transform_pvalues_gsd(raw2, info_frac = info, spending = sfLDOF)

# each printed pair below is (simulation 1, simulation 2) for the one hypothesis
cat("info_frac:\n"); print(info)
cat("raw look 2      :", raw2[, 1, 2], "\n")
cat("repeated look 1 :", tr2$pvals[, 1, 1], "\n")
cat("repeated look 2 :", tr2$pvals[, 1, 2], "\n")
cat("table built?", !is.null(tr2$tables[[1]]$bounds), "\n")

stopifnot(
    !identical(boundary, a),
    abs(boundary - a) < 1e-15,
    identical(raw[, , 1], tr$pvals[, , 1]),
    is.null(tr$tables[[1]]$bounds),
    all(tr2$pvals[, 1, 1] == 1),
    identical(tr2$pvals[, 1, 2], raw2[, 1, 2])
)
cat("\nOK\n")
