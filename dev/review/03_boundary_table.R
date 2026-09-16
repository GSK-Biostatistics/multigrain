# Exercise 03: the boundary table
#
# `.gsd_boundary_table()` calls the forward map of exercise 02 once per grid
# level and stacks the results: a matrix with one row per level and one column
# per look. This exercise builds the same table, looks at it, and then breaks
# the well-ordering check on purpose.
#
# The grid is log-spaced from gsd_grid_min = 1e-14 to alpha. Log spacing
# matters because boundaries span many orders of magnitude (the LDOF look-1
# boundary at gamma = 1e-6 is around 1e-20) and we want uniform *relative*
# accuracy when we interpolate later.

suppressMessages(devtools::load_all(quiet = TRUE))
library(gsDesign)

alpha <- 0.025
t     <- c(0.5, 1)

# --- 1. Build the table with the package's own helper -----------------------

tab <- multigrain:::.gsd_boundary_table(
    t_look = t, spending = sfLDOF, alpha = alpha, grid_size = 128, hyp = 1
)
str(tab)

cat("\nfirst 5 rows (levels near 1e-14):\n")
print(cbind(level = tab$grid[1:5], tab$bounds[1:5, ]))
cat("\nlast 5 rows (levels near alpha):\n")
n <- length(tab$grid)
print(cbind(level = tab$grid[(n - 4):n], tab$bounds[(n - 4):n, ]))

# >>> Question: look at the look-1 column in the first rows. It is a constant
#     2.75e-89. That is not a real boundary; it is gsBound1 underflowing when
#     asked for a crossing probability of ~1e-30. Why does this not matter
#     for the transform? (Exercise 04 shows what .gsd_invert() does with it.)

# --- 2. Plot it -------------------------------------------------------------
#
# On log-log axes each column should be a smooth, increasing curve. The
# look-2 curve sits close to the diagonal (boundary ~ level) because LDOF
# spends almost everything at the end; the look-1 curve is far below it.

keep1 <- tab$bounds[, 1] > 1e-80
if (interactive()) {
    plot(tab$grid, tab$bounds[, 2], log = "xy", type = "l", lwd = 2,
         xlab = "allocated level gamma", ylab = "nominal boundary",
         main = "Boundary table, LDOF, t = (0.5, 1)")
    lines(tab$grid[keep1], tab$bounds[keep1, 1], lwd = 2, col = "steelblue")
    abline(0, 1, lty = 3)
    legend("topleft", c("look 2", "look 1", "boundary = level"),
           col = c("black", "steelblue", "black"), lty = c(1, 1, 3), lwd = c(2, 2, 1))
} else {
    cat("\n(run interactively to see the plot)\n")
}

# --- 3. Monotone in every column? -------------------------------------------

cat("\nlook 1 non-decreasing:", !is.unsorted(tab$bounds[, 1]), "\n")
cat("look 2 non-decreasing:", !is.unsorted(tab$bounds[, 2]), "\n")

# --- 4. Break the well-ordering check ---------------------------------------
#
# A "well ordered" spending function is one whose boundaries rise with the
# level at every look (Maurer and Bretz condition 2). Here is one that is
# not: it spends *less* at look 1 as the level goes up. Such a function is
# not a valid alpha-spending design for the graphical procedure, and the
# transform must refuse it rather than silently produce nonsense.

bad_sf <- function(a, t) {
    early <- pmin(a, 0.002 * (1 - a / 0.025))   # decreasing in a
    ifelse(t >= 1, a, early)
}
cat("\nbad_sf look-1 spend at a = 0.005, 0.015, 0.025:",
    bad_sf(c(0.005, 0.015, 0.025), 0.5), "\n")

res <- tryCatch(
    multigrain:::.gsd_boundary_table(
        t_look = t, spending = bad_sf, alpha = alpha, grid_size = 64, hyp = 1
    ),
    error = function(e) e
)
cat("\nwell-ordering check said:\n")
cat(conditionMessage(res), "\n")

# >>> Question: which analysis did the error name, and why that one and not
#     the other?

stopifnot(
    !is.unsorted(tab$bounds[, 1]),
    !is.unsorted(tab$bounds[, 2]),
    inherits(res, "error"),
    grepl("well ordered", conditionMessage(res))
)
cat("\nOK\n")
