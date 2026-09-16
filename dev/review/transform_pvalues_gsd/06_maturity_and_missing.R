# Exercise 06: endpoints that start late or finish early
#
# Real designs are ragged: PFS may be fully mature at the first analysis
# while OS is still accruing; a secondary endpoint may have no data at all
# until look 2. The transform handles this with two conventions (design
# record, "Two conventions the implementer must know"):
#
#   * no data yet          -> repeated p = 1 ("cannot reject")
#   * matured at look k    -> repeated p at looks > k is a copy of look k
#
# `.gsd_looks()` decides which looks carry distinct information for each
# hypothesis; `.gsd_transform_hyp()` applies the conventions.

suppressMessages(devtools::load_all(quiet = TRUE))
library(gsDesign)

# --- 1. What .gsd_looks() returns ------------------------------------------

show_looks <- function(t_row) {
    r <- multigrain:::.gsd_looks(t_row, hyp = 1)
    cat(sprintf("info_frac = (%s)  ->  looks = {%s}, maturity = %d\n",
                paste(format(t_row), collapse = ", "),
                paste(r$looks, collapse = ", "), r$maturity))
}
show_looks(c(0.5, 1))          # ordinary
show_looks(c(1, 1))            # matured at look 1; look 2 carries nothing new
show_looks(c(NA, 1))           # no data until look 2
show_looks(c(0.4, 0.7))        # never reaches 1: last look is "maturity"
show_looks(c(NA, 0.6, 1))      # starts late, matures at 3
show_looks(c(0.5, 1, 1, 1))    # matured at 2; three trailing copies

# >>> Question: for (0.5, 1, 1, 1), why must the boundary table be built for
#     looks {1, 2} only and not {1, 2, 3, 4}? (Design record: "Never ask
#     gsDesign for a boundary with a repeated information fraction of 1".)
#     Try it directly:

inc <- diff(c(0, sfLDOF(0.025, c(0.5, 1, 1))$spend))
cat("\nspend increments for t = (0.5, 1, 1):", inc, "\n")
b <- gsBound1(theta = 0, I = c(0.5, 1, 1), a = rep(-20, 3), probhi = inc)$b
cat("gsBound1 z-boundaries:", b, "\n")
cat("as p-value boundaries:", pnorm(b, lower.tail = FALSE), "\n")
cat("(a zero increment gives an unrejectable boundary at look 3)\n\n")

# --- 2. The Example 5 shape: PFS matures early, OS runs long ----------------

set.seed(4)
nsim <- 5
K    <- 3
info <- rbind(
    PFS = c(1,   1,   1),      # complete at the first analysis
    OS  = c(0.4, 0.7, 1)       # accrues through all three
)
raw <- array(runif(nsim * 2 * K), dim = c(nsim, 2, K))
# make PFS's raw p identical across looks, as a matured endpoint's would be
raw[, 1, 2] <- raw[, 1, 1]
raw[, 1, 3] <- raw[, 1, 1]

tr <- transform_pvalues_gsd(raw, info_frac = info, spending = sfLDOF)
summary(tr)

cat("\nPFS raw (look 1):     ", signif(raw[, 1, 1], 4), "\n")
cat("PFS repeated, look 1: ", signif(tr$pvals[, 1, 1], 4), "\n")
cat("PFS repeated, look 2: ", signif(tr$pvals[, 1, 2], 4), "\n")
cat("PFS repeated, look 3: ", signif(tr$pvals[, 1, 3], 4), "\n")

# >>> Question: PFS's repeated p equals its raw p at every look. Why (two
#     reasons: single matured look -> short-circuit; then copied forward)?
#     And why is it useful that it is *not* set to 1 after look 1? (Hint:
#     what if OS is rejected at look 3 and recycles alpha to PFS?)

# --- 3. A hypothesis with nothing at look 1 --------------------------------

info2 <- rbind(
    A = c(0.5, 1),
    B = c(NA,  1)              # B has no data until the final analysis
)
raw2 <- array(runif(nsim * 2 * 2), dim = c(nsim, 2, 2))
tr2  <- transform_pvalues_gsd(raw2, info_frac = info2, spending = sfLDOF)

cat("\nB repeated, look 1:", tr2$pvals[, 2, 1], "\n")
cat("B raw,      look 2:", signif(raw2[, 2, 2], 4), "\n")
cat("B repeated, look 2:", signif(tr2$pvals[, 2, 2], 4), "\n")
cat("(look 1 is all 1s: no data means cannot reject; look 2 is the raw value,\n",
    " because a single look at full information needs no table)\n")

# --- 4. What if you supply different raw p-values after maturity? ----------
#
# The transform ignores them (the maturity value is copied forward) and
# warns, because the caller has supplied data the procedure will not use.

raw3 <- raw
raw3[, 1, 3] <- 0.5           # PFS "changes" at look 3 -- not possible if matured
res <- withCallingHandlers(
    transform_pvalues_gsd(raw3, info_frac = info, spending = sfLDOF),
    warning = function(w) { cat("\nWARNING:", conditionMessage(w), "\n"); invokeRestart("muffleWarning") }
)
cat("PFS repeated, look 3, still:", signif(res$pvals[, 1, 3], 4), "\n")

stopifnot(
    identical(tr$pvals[, 1, 1], raw[, 1, 1]),
    identical(tr$pvals[, 1, 3], raw[, 1, 1]),
    all(tr2$pvals[, 2, 1] == 1),
    identical(tr2$pvals[, 2, 2], raw2[, 2, 2]),
    identical(res$pvals[, 1, 3], raw[, 1, 1]),
    is.null(tr$tables[[1]]$bounds),      # PFS: no table
    !is.null(tr$tables[[2]]$bounds)      # OS: table over looks 1..3
)
cat("\nOK\n")
