# Exercise 04: inverting the table
#
# The repeated p-value of a raw p-value p at look k is the level gamma at
# which the look-k boundary equals p. The table from exercise 03 gives the
# boundary at each grid level; `.gsd_invert()` interpolates between grid
# points to find the gamma for an arbitrary p.
#
# The interpolation is linear in (log boundary, log level), because both
# axes span many orders of magnitude and the curve is close to a straight
# line on log-log axes (exercise 03's plot).

suppressMessages(devtools::load_all(quiet = TRUE))
library(gsDesign)

alpha <- 0.025
t     <- c(0.5, 1)
tab   <- multigrain:::.gsd_boundary_table(
    t_look = t, spending = sfLDOF, alpha = alpha, grid_size = 1024, hyp = 1
)

# --- 1. Invert a few p-values at look 2 ------------------------------------

p <- c(0.02, 0.01, 0.001, 1e-6)
pr <- multigrain:::.gsd_invert(tab$bounds[, 2], grid = tab$grid, p = p)
print(data.frame(raw_p = p, repeated_p = signif(pr, 6)))

# >>> Question: repeated_p is slightly *larger* than raw_p at look 2. Why?
#     (At look 2 with LDOF, the boundary at level gamma is a bit below gamma,
#     so you need a bit more than p of allocation to reach a boundary of p.)

# --- 2. Check against a root-finder ----------------------------------------
#
# The "true" answer is the gamma solving boundary_k(gamma) = p exactly.

nom_bounds <- function(gamma, t, sf = sfLDOF) {
    inc <- diff(c(0, sf(gamma, t)$spend))
    b <- gsBound1(theta = 0, I = t, a = rep(-20, length(t)), probhi = inc)$b
    pnorm(b, lower.tail = FALSE)
}
exact <- function(p, k) {
    uniroot(function(g) nom_bounds(g, t)[k] - p, c(1e-13, alpha), tol = 1e-12)$root
}
pr_exact <- sapply(p, exact, k = 2)
print(data.frame(
    raw_p = p,
    interpolated = signif(pr, 8),
    uniroot = signif(pr_exact, 8),
    rel_error = signif(abs(pr - pr_exact) / pr_exact, 2)
))

# --- 3. Error against grid size -------------------------------------------
#
# The design record quotes ~1e-6 relative error at G = 1024. Reproduce the
# trend with a handful of random p-values.

set.seed(2)
p_test <- 10^runif(40, -8, log10(0.02))
truth  <- sapply(p_test, exact, k = 2)
for (G in c(128, 256, 512, 1024)) {
    tb <- multigrain:::.gsd_boundary_table(
        t_look = t, spending = sfLDOF, alpha = alpha, grid_size = G, hyp = 1
    )
    est <- multigrain:::.gsd_invert(tb$bounds[, 2], grid = tb$grid, p = p_test)
    cat(sprintf("G = %4d  max relative error = %.2e\n", G,
                max(abs(est - truth) / truth)))
}

# >>> Question: the error roughly quarters each time G doubles. What does
#     that tell you about the interpolation scheme?

# --- 4. The two clamps -----------------------------------------------------
#
# Above the table: p is bigger than the boundary at the largest level alpha,
# so no allocation up to alpha would reject it. Repeated p = 1.
#
# Below the table: p is smaller than the smallest positive boundary. The
# "true" repeated p is below 1e-14. It is floored at gsd_grid_min = 1e-14,
# NOT 0.

top    <- max(tab$bounds[, 2])
bottom <- min(tab$bounds[tab$bounds[, 2] > 0, 2])
cat("\nlook-2 table runs from", bottom, "to", top, "\n")

edge <- c(top * 1.01, 0.5, bottom * 0.5, 1e-100, NA)
pr_edge <- multigrain:::.gsd_invert(tab$bounds[, 2], grid = tab$grid, p = edge)
print(data.frame(raw_p = edge, repeated_p = pr_edge))

# Why not floor at 0? The optimiser's snaps can leave a hypothesis an
# allocation as tiny as ~2.5e-11. With a 0 floor, ANY p below the table would
# be rejected at that allocation, even if its true repeated p is, say, 5e-11
# (which is above 2.5e-11 and should NOT reject). With the 1e-14 floor the
# comparison 1e-14 <= 2.5e-11 still rejects, which is correct whenever the
# true value is below 1e-14, and only wrong (conservatively) in the sliver
# between. Design record 4.1, "The floor is gamma_min rather than 0".

cat("\nfloor value:", multigrain:::gsd_grid_min, "\n")

# --- 5. The underflow filter -----------------------------------------------
#
# Look 1's column has a run of 2.75e-89 at the bottom (exercise 03). approx()
# needs strictly increasing x, so .gsd_invert() drops duplicated and
# non-positive entries before interpolating. Show it still inverts look 1.

p1  <- c(1e-3, 1e-5, 1e-8)
pr1 <- multigrain:::.gsd_invert(tab$bounds[, 1], grid = tab$grid, p = p1)
ex1 <- sapply(p1, exact, k = 1)
print(data.frame(raw_p = p1, interpolated = signif(pr1, 6), uniroot = signif(ex1, 6)))

# >>> Question: a look-1 p of 1e-3 has a repeated p around 0.02: the raw value
#     looks impressive but is barely enough at the full 0.025. Explain that
#     to a clinician in one sentence.

stopifnot(
    all(abs(pr - pr_exact) / pr_exact < 1e-5),
    pr_edge[1] == 1, pr_edge[2] == 1,
    pr_edge[3] == multigrain:::gsd_grid_min,
    pr_edge[4] == multigrain:::gsd_grid_min,
    pr_edge[5] == 1,
    all(abs(pr1 - ex1) / ex1 < 1e-5)
)
cat("\nOK\n")
