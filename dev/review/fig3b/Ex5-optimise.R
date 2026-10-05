# Ex5-optimise.R  --  Example 5: Group sequential design, PFS/OS (Section 3.5)

source("00-shared-settings.R")
library(mvtnorm)
library(gsDesign)
library(Rcpp)

sourceCpp("./src/quick_pfs_os.cpp")
set.seed(202605)

t_os <- 0.7
rho <- 0.5
timing <- c(t_os, 1)
N_SIM <- 5e6

Delta_P <- calc_ncp(0.98, alpha = ALPHA) # ~ 4.01
Delta_O <- calc_ncp(0.93, alpha = ALPHA) # ~ 3.44

get_of_bounds <- function(alpha_one_sided, timing) {
  x <- gsDesign(
    k = 2,
    test.type = 1,
    alpha = alpha_one_sided,
    sfu = "OF",
    timing = timing
  )
  1 - pnorm(x$upper$bound)
}

os_bounds_of_from_w2 <- function(w2, alpha, timing) {
  full <- get_of_bounds(alpha, timing)
  init <- if (w2 <= 0) c(0, 0) else get_of_bounds(alpha * w2, timing)
  rbind(init, full)
}

make_cov_3 <- function(t_os, rho) {
  matrix(
    c(
      1,
      sqrt(t_os) * rho,
      rho,
      sqrt(t_os) * rho,
      1,
      sqrt(t_os),
      rho,
      sqrt(t_os),
      1
    ),
    nrow = 3,
    byrow = TRUE
  )
}

# Simulate joint (Z_PFS, Z_OS_interim, Z_OS_final) once; reuse across the grid.
mu <- c(Delta_P, Delta_O * sqrt(t_os), Delta_O)
Sigma <- make_cov_3(t_os, rho)
Z <- rmvnorm(N_SIM, mean = mu, sigma = Sigma)
pvals <- cbind(
  pP = 1 - pnorm(Z[, 1]),
  pO1 = 1 - pnorm(Z[, 2]),
  pO2 = 1 - pnorm(Z[, 3])
)

# Rejection times depend on w1 alone, so the grid is swept once and the
# probability of rejection at each analysis recorded.
rejection_probs <- function(w1) {
  osb <- os_bounds_of_from_w2(1 - w1, ALPHA, timing)
  out <- graph_2h_time(c(w1, 1 - w1), pvals, os_bounds = osb, alpha = ALPHA)
  tP <- out$time[, 1]
  tO <- out$time[, 2]
  c(
    pfs_interim = mean(tP == 1L),
    pfs_final = mean(tP == 2L),
    os_interim = mean(tO == 1L),
    os_final = mean(tO == 2L)
  )
}

w1_grid <- seq(0, 1, by = 0.0005)
probs <- t(vapply(w1_grid, rejection_probs, numeric(4)))

# Value 1 for rejection at the interim, delta at the final analysis, 0 if never
# rejected. Endpoint values are normalised so that V_P + V_O = 1.
Egain <- function(r, delta) {
  V_P <- 1 / (1 + r)
  V_O <- r / (1 + r)
  V_P *
    (probs[, "pfs_interim"] + delta * probs[, "pfs_final"]) +
    V_O * (probs[, "os_interim"] + delta * probs[, "os_final"])
}

find_w1_star <- function(r, delta) {
  gain <- Egain(r, delta)
  best <- which.max(gain)
  data.frame(
    r = r,
    delta = delta,
    w1_star = w1_grid[best],
    Egain_star = gain[best]
  )
}

r_grid <- c(0.25, 0.5, 1, 2, 4, 8, 16)
delta_grid <- c(1.0, 0.75, 0.5)

grid <- expand.grid(r = r_grid, delta = delta_grid, KEEP.OUT.ATTRS = FALSE)
res <- do.call(rbind, Map(find_w1_star, grid$r, grid$delta))
res <- res[order(res$delta, res$r), ]

results <- list(
  example = 5L,
  inputs = list(
    alpha = ALPHA,
    t_os = t_os,
    rho = rho,
    timing = timing,
    Delta_P = Delta_P,
    Delta_O = Delta_O,
    n_sim = N_SIM,
    s1 = 37,
    s2 = 47
  ),
  r_grid = r_grid,
  delta_grid = delta_grid,
  w_star = res
)
saveRDS(results, file.path(RESULTS_DIR, "Ex5-results.rds"))

cat("\nExample 5 -- NCPs (Table 6):\n")
cat(sprintf(
  "  Delta_PFS = %.2f (target 4.01),  Delta_OS = %.2f (target 3.44)\n",
  Delta_P,
  Delta_O
))
cat(sprintf("  Interim OS mean = %.2f (target 2.87)\n", Delta_O * sqrt(t_os)))
cat("\nOptimal w*_PFS by (r, delta):\n")
print(res, row.names = FALSE)
