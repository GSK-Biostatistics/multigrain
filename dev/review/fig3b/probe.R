suppressPackageStartupMessages({
    library(gsDesign)
    library(mvtnorm)
    library(Rcpp)
})
S <- "dev/review/fig3b" # run from the repo root; originally the session scratchpad
sourceCpp(file.path(S, "quick_pfs_os.cpp"))
ref <- readRDS(file.path(S, "Ex5-results.rds"))
cat("== Reference w* (N = 5e6, OF boundary) ==\n")
print(ref$w_star, row.names = FALSE)

ALPHA <- 0.025
t_os <- 0.7
rho <- 0.5
timing <- c(t_os, 1)
calc_ncp <- function(power, alpha) qnorm(1 - alpha) - qnorm(1 - power)
Delta_P <- calc_ncp(0.98, ALPHA)
Delta_O <- calc_ncp(0.93, ALPHA)

of_bounds <- function(a, timing) {
    x <- gsDesign(k = 2, test.type = 1, alpha = a, sfu = "OF", timing = timing)
    1 - pnorm(x$upper$bound)
}
ldof_bounds <- function(a, timing) {
    inc <- diff(c(0, sfLDOF(a, timing)$spend))
    b <- gsBound1(theta = 0, I = timing, a = rep(-20, 2), probhi = inc)$b
    pnorm(b, lower.tail = FALSE)
}
cat("\n== OF vs LDOF nominal bounds at timing (0.7, 1) ==\n")
for (a in c(0.025, 0.0125, 0.005, 1e-3, 1e-4, 1e-6)) {
    cat(sprintf(
        "alpha=%-8g OF: %s   LDOF: %s\n", a,
        paste(format(of_bounds(a, timing), digits = 4), collapse = " "),
        paste(format(ldof_bounds(a, timing), digits = 4), collapse = " ")
    ))
}
# implied spending of the OF boundary (cumulative), as a spending function
of_spend <- function(a, t) {
    x <- gsDesign(k = length(t), test.type = 1, alpha = a, sfu = "OF", timing = t)
    cumsum(x$upper$spend)
}
cat("\nOF implied cumulative spend at alpha=0.025:", format(of_spend(0.025, timing), digits = 10), "\n")
cat("sfLDOF cumulative spend at alpha=0.025:    ", format(sfLDOF(0.025, timing)$spend, digits = 10), "\n")
cat("OF spend(a,1) - a, relative:", (of_spend(0.025, timing)[2] - 0.025) / 0.025, "\n")
tm <- system.time(for (i in 1:20) of_spend(0.01, timing))
cat("time per OF gsDesign() call (ms):", 1000 * tm[["elapsed"]] / 20, "\n")

# ---- sweep: paper's procedure, with either boundary family, at given N/seed ----
sweep <- function(N, seed, bounds_fn, w1_grid) {
    set.seed(seed)
    mu <- c(Delta_P, Delta_O * sqrt(t_os), Delta_O)
    Sigma <- matrix(c(1, sqrt(t_os) * rho, rho, sqrt(t_os) * rho, 1, sqrt(t_os), rho, sqrt(t_os), 1), 3, byrow = TRUE)
    Z <- rmvnorm(N, mean = mu, sigma = Sigma)
    pv <- cbind(1 - pnorm(Z[, 1]), 1 - pnorm(Z[, 2]), 1 - pnorm(Z[, 3]))
    full <- bounds_fn(ALPHA, timing)
    probs <- t(vapply(w1_grid, function(w1) {
        w2 <- 1 - w1
        init <- if (w2 <= 0) c(0, 0) else bounds_fn(ALPHA * w2, timing)
        out <- graph_2h_time(c(w1, w2), pv, os_bounds = rbind(init, full), alpha = ALPHA)
        c(mean(out$time[, 1] == 1L), mean(out$time[, 1] == 2L), mean(out$time[, 2] == 1L), mean(out$time[, 2] == 2L))
    }, numeric(4)))
    res <- ref$w_star
    res$w1_hat <- NA_real_
    res$gain_gap <- NA_real_
    for (j in seq_len(nrow(res))) {
        r <- res$r[j]; d <- res$delta[j]
        vp <- 1 / (1 + r); vo <- r / (1 + r)
        gain <- vp * (probs[, 1] + d * probs[, 2]) + vo * (probs[, 3] + d * probs[, 4])
        res$w1_hat[j] <- w1_grid[which.max(gain)]
        # gain at the reference w* (nearest grid point) vs the max on the grid
        at_ref <- gain[which.min(abs(w1_grid - res$w1_star[j]))]
        res$gain_gap[j] <- max(gain) - at_ref
    }
    res
}
w1_grid <- seq(0, 1, by = 0.005)

cat("\n== OF boundary, N = 1e5, seeds 1..3: |w1_hat - w1_star| ==\n")
for (s in 1:3) {
    r <- sweep(1e5, s, of_bounds, w1_grid)
    cat(sprintf("seed %d: max abs dev %.3f; mean %.3f; max gain gap %.2e; devs: %s\n", s,
        max(abs(r$w1_hat - r$w1_star)), mean(abs(r$w1_hat - r$w1_star)), max(r$gain_gap),
        paste(format(round(r$w1_hat - r$w1_star, 3)), collapse = " ")))
}
cat("\n== LDOF spending, N = 1e5, seed 1 ==\n")
r <- sweep(1e5, 1, ldof_bounds, w1_grid)
cat(sprintf("max abs dev %.3f; mean %.3f; max gain gap %.2e\n", max(abs(r$w1_hat - r$w1_star)), mean(abs(r$w1_hat - r$w1_star)), max(r$gain_gap)))
print(r[, c("r", "delta", "w1_star", "w1_hat", "gain_gap")], row.names = FALSE)
cat("\n== LDOF spending, N = 1e6, seed 1 ==\n")
r <- sweep(1e6, 1, ldof_bounds, w1_grid)
cat(sprintf("max abs dev %.3f; mean %.3f; max gain gap %.2e\n", max(abs(r$w1_hat - r$w1_star)), mean(abs(r$w1_hat - r$w1_star)), max(r$gain_gap)))
print(r[, c("r", "delta", "w1_star", "w1_hat", "gain_gap")], row.names = FALSE)
