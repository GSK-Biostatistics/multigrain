# check1_tables.R -- Does inverting gsDesign boundaries reproduce Maurer & Bretz
# (2013) Tables 1 and 2, and does the fixed-sample cascade on repeated p-values
# reproduce their case-study rejections?
suppressPackageStartupMessages(library(gsDesign))
options(digits = 6)

# Forward map: nominal p-value boundaries alpha*_k(gamma) for a one-sided
# alpha-spending design at information fractions t, spending function sf.
nom_bounds <- function(gamma, t, sf = sfLDOF) {
  K <- length(t)
  inc <- diff(c(0, sf(gamma, t)$spend))                     # per-analysis spend
  b <- gsBound1(theta = 0, I = t, a = rep(-20, K), probhi = inc)$b
  pnorm(b, lower.tail = FALSE)
}

t3 <- c(1/3, 2/3, 1)
cat("== Maurer-Bretz Table 1 (LDOF, t = 1/3, 2/3, 1) ==\n")
cat("gamma=0.0125  :", format(nom_bounds(0.0125, t3), digits = 4),  "  paper: 0.00002 0.0022 -\n")
cat("gamma=0.01875 :", format(nom_bounds(0.01875, t3), digits = 4), "  paper k=2: 0.004\n")
cat("gamma=0.00625 :", format(nom_bounds(0.00625, t3), digits = 4), "  paper k=2: 0.0008\n")
cat("gamma=0.025   :", format(nom_bounds(0.025, t3), digits = 4),   "  paper k=2: 0.006, k=3: 0.02313\n")
d <- gsDesign(k = 3, test.type = 1, alpha = 0.0125, timing = t3, sfu = sfLDOF)
cat("cross-check gsDesign():", format(pnorm(d$upper$bound, lower.tail = FALSE), digits = 4), "\n")

# Inverse: repeated p-value = the gamma solving alpha*_k(gamma) = p
rep_p <- function(p, k, t, sf = sfLDOF) {
  uniroot(function(g) nom_bounds(g, t, sf)[k] - p, c(1e-12, 0.999), tol = 1e-10)$root
}
cat("\n== Maurer-Bretz Table 2 (repeated p-values, ADDPLAN) ==\n")
p1 <- c(0.0062, 0.017, 0.009, 0.13); p2 <- c(0.0002, 0.0035, 0.002, 0.06)
r1 <- sapply(p1, rep_p, k = 1, t = t3); r2 <- sapply(p2, rep_p, k = 2, t = t3)
cat("k=1:", format(r1, digits = 4), "  paper: 0.1141 0.1683 0.1316 0.382\n")
cat("k=2:", format(r2, digits = 4), "  paper: 0.0024 0.0172 0.0117 0.1285\n")

cat("\n== Fixed-sample cascade on the k=2 repeated p-values (Figure 1 graph) ==\n")
w <- c(0.5, 0.5, 0, 0)
G <- rbind(c(0, .5, .5, 0), c(.5, 0, 0, .5), c(0, 1, 0, 0), c(1, 0, 0, 0))
alpha <- 0.025; p <- r2; rej <- rep(FALSE, 4); order <- integer(0)
repeat {
  j <- which(!rej & p <= w * alpha)[1]; if (is.na(j)) break
  rej[j] <- TRUE; order <- c(order, j)
  wn <- w + w[j] * G[j, ]; wn[j] <- 0
  Gn <- G
  for (l in 1:4) for (k in 1:4)
    Gn[l, k] <- if (l == k || rej[l] || rej[k]) 0 else (G[l, k] + G[l, j] * G[j, k]) / (1 - G[l, j] * G[j, l])
  w <- wn; G <- Gn
}
cat("rejected:", rej, " order:", order, "  paper: H1, H2, H3 rejected; H4 retained\n")
