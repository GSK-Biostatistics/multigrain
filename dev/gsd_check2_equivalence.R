# check2_equivalence.R -- Is "transform once + K passes of the fixed-sample
# cascade" identical to a direct boundary-on-the-fly implementation of the
# Maurer-Bretz procedure, for mixed per-hypothesis look-back and matured
# hypotheses?  Also: gsBound1 edge cases, grid accuracy, clamp safety, cost.
suppressPackageStartupMessages(library(gsDesign))
options(digits = 6, warn = -1)
alpha <- 0.025
g_min <- 1e-14
nom_bounds <- function(gamma, t, sf = sfLDOF) {
    K <- length(t)
    inc <- diff(c(0, sf(gamma, t)$spend))
    b <- gsBound1(theta = 0, I = t, a = rep(-20, K), probhi = inc)$b
    pnorm(b, lower.tail = FALSE)
}

cat("== gsBound1 edge cases ==\n")
cat(
    "duplicate t=c(0.7,1,1):",
    format(nom_bounds(0.025, c(0.7, 1, 1)), digits = 4),
    "  <- third look unrejectable: never pass repeated t=1\n"
)
cat(
    "t > 1, t=c(0.5,1.2):   ",
    format(nom_bounds(0.025, c(0.5, 1.2)), digits = 4),
    "\n"
)
cat(
    "K=1, t=1, gamma=0.0125:",
    format(nom_bounds(0.0125, 1), digits = 17),
    " (identical to gamma?",
    identical(nom_bounds(0.0125, 1), 0.0125),
    ")\n"
)
for (g in c(1e-4, 1e-8, 1e-14)) {
    cat(sprintf(
        "gamma=%g, t=(1/3,2/3,1): %s\n",
        g,
        paste(
            format(nom_bounds(g, c(1 / 3, 2 / 3, 1)), digits = 3),
            collapse = " "
        )
    ))
}

cat(
    "\n== Grid size vs inverse accuracy (t=(0.7,1), levels in [1e-14, alpha], log-log linear) ==\n"
)
t2 <- c(0.7, 1)
rep_p <- function(p, k, t) {
    uniroot(
        function(g) nom_bounds(g, t)[k] - p,
        c(1e-15, 0.999),
        tol = 1e-13
    )$root
}
set.seed(2)
ptest <- 10^runif(200, -9, log10(0.03))
ex <- lapply(1:2, function(k) sapply(ptest, rep_p, k = k, t = t2))
for (Gn in c(256, 512, 1024)) {
    grid <- exp(seq(log(g_min), log(alpha), length.out = Gn))
    fwd <- do.call(rbind, lapply(grid, nom_bounds, t = t2))
    for (k in 1:2) {
        f <- fwd[, k]
        ok <- !duplicated(f) & f > 0
        ap <- exp(
            approx(log(f[ok]), log(grid[ok]), xout = log(ptest), rule = 1)$y
        )
        sel <- ex[[k]] <= alpha & !is.na(ap)
        cat(sprintf(
            "G=%4d k=%d: n=%3d  max rel err %.2e\n",
            Gn,
            k,
            sum(sel),
            max(abs(ex[[k]] - ap)[sel] / ex[[k]][sel])
        ))
    }
}

cat(
    "\n== End-to-end equivalence: direct vs transform-once + cascade (m=3, K=3, N=400) ==\n"
)
set.seed(11)
N <- 400
m <- 3
K <- 3
tmat <- rbind(c(1 / 3, 2 / 3, 1), c(0.5, 1, 1), c(1, 1, 1)) # H2 matures at k=2, H3 at k=1
sfs <- list(sfLDOF, function(a, t) sfHSD(a, t, param = -2), sfLDOF)
w0 <- c(0.5, 0.3, 0.2)
G0 <- rbind(c(0, .5, .5), c(.5, 0, .5), c(.5, .5, 0))
Delta <- c(2.6, 2.4, 2.9)
mature_k <- apply(tmat, 1, function(t) which(t >= 1)[1])
eff_l <- function(i, k) min(k, mature_k[i])
cache <- new.env()
bound <- function(i, l, a) {
    key <- sprintf("%d_%d_%.15g", i, l, a)
    v <- cache[[key]]
    if (is.null(v)) {
        v <- nom_bounds(a, tmat[i, 1:l], sfs[[i]])[l]
        assign(key, v, envir = cache)
    }
    v
}
P <- lapply(1:m, function(i) {
    kk <- mature_k[i]
    t <- tmat[i, 1:kk]
    S <- outer(t, t, function(a, b) sqrt(pmin(a, b) / pmax(a, b)))
    Z <- mvtnorm::rmvnorm(N, mean = Delta[i] * sqrt(t), sigma = S)
    p <- pnorm(Z, lower.tail = FALSE)
    cbind(p, matrix(p[, kk], N, K - kk))
}) # matured p carried forward
upd <- function(w, G, i, rej) {
    wn <- w + w[i] * G[i, ]
    wn[i] <- 0
    Gn <- matrix(0, m, m)
    for (l in 1:m) {
        for (j in 1:m) {
            if (l != j && !rej[l] && !rej[j]) {
                d <- 1 - G[l, i] * G[i, l]
                Gn[l, j] <- if (d > 0) (G[l, j] + G[l, i] * G[i, j]) / d else 0
            }
        }
    }
    list(w = wn, G = Gn)
}
# Direct Maurer-Bretz procedure: boundaries computed on the fly at the current weight
direct <- function(pk, lb) {
    w <- w0
    G <- G0
    rej <- rep(FALSE, m)
    tm <- integer(m)
    for (k in 1:K) {
        repeat {
            hit <- NA
            for (i in which(!rej)) {
                a <- w[i] * alpha
                if (a <= 0) {
                    next
                }
                ls <- if (lb[i]) 1:eff_l(i, k) else eff_l(i, k)
                if (any(pk[i, ls] < sapply(ls, function(l) bound(i, l, a)))) {
                    hit <- i
                    break
                }
            }
            if (is.na(hit)) {
                break
            }
            rej[hit] <- TRUE
            tm[hit] <- k
            u <- upd(w, G, hit, rej)
            w <- u$w
            G <- u$G
        }
    }
    c(rej, tm)
}
# Transform: grid tables per hypothesis, inverse with floor clamp and first-look short-circuit
grids <- lapply(1:m, function(i) {
    g <- exp(seq(log(g_min), log(alpha), length.out = 1024))
    kk <- mature_k[i]
    list(
        g = g,
        fwd = do.call(
            rbind,
            lapply(g, function(gm) nom_bounds(gm, tmat[i, 1:kk], sfs[[i]]))
        )
    )
})
inv <- function(i, l, p) {
    if (l == 1 && tmat[i, 1] >= 1) {
        return(p)
    } # short-circuit: boundary is the level itself
    g <- grids[[i]]$g
    f <- grids[[i]]$fwd[, l]
    ok <- !duplicated(f) & f > 0
    out <- exp(approx(log(f[ok]), log(g[ok]), xout = log(p), rule = 1)$y)
    out[p >= max(f)] <- 1
    out[p < min(f[ok])] <- g_min
    out
}
transform <- function(lb) {
    A <- array(NA_real_, c(N, m, K))
    for (i in 1:m) {
        for (k in 1:K) {
            l <- eff_l(i, k)
            A[, i, k] <- inv(i, l, P[[i]][, l])
        }
    }
    for (i in which(lb)) {
        A[, i, ] <- t(apply(A[, i, ], 1, cummin))
    }
    A
}
kern <- function(pr) {
    w <- w0
    G <- G0
    rej <- rep(FALSE, m)
    tm <- integer(m)
    for (k in 1:K) {
        repeat {
            i <- which(!rej & pr[, k] < w * alpha)[1]
            if (is.na(i)) {
                break
            }
            rej[i] <- TRUE
            tm[i] <- k
            u <- upd(w, G, i, rej)
            w <- u$w
            G <- u$G
        }
    }
    c(rej, tm)
}
for (lb in list(c(F, F, F), c(T, T, T), c(T, F, T), c(F, T, F))) {
    A <- transform(lb)
    D <- t(sapply(1:N, function(n) {
        direct(rbind(P[[1]][n, ], P[[2]][n, ], P[[3]][n, ]), lb)
    }))
    Tr <- t(sapply(1:N, function(n) kern(A[n, , ])))
    cat(sprintf(
        "look_back=%s: rejection mismatches %d/%d, time mismatches %d/%d; local power %s; mean decision time %s\n",
        paste(substr(as.character(lb), 1, 1), collapse = ""),
        sum(D[, 1:m] != Tr[, 1:m]),
        N * m,
        sum(D[, -(1:m)] != Tr[, -(1:m)]),
        N * m,
        paste(format(colMeans(Tr[, 1:m]), digits = 3), collapse = " "),
        paste(
            format(
                colMeans(
                    replace(Tr[, -(1:m)], Tr[, -(1:m)] == 0, NA),
                    na.rm = TRUE
                ),
                digits = 3
            ),
            collapse = " "
        )
    ))
}
cat(
    "matured-at-first-look hypothesis: transformed column identical() to raw p:",
    identical(transform(c(F, F, F))[, 3, 1], P[[3]][, 1]),
    "\n"
)

cat(
    "\n== Clamp safety (red-team check): tiny allocation 1e-4 * 1e-5 * alpha after the kernel's snaps ==\n"
)
a_tiny <- 1e-4 * 1e-5 * alpha
p_row <- 5e-11 # true p^r at a t=1 look ~ p = 5e-11 > a_tiny
old_floor <- 1e-10
g_old <- exp(seq(log(old_floor), log(alpha), length.out = 1024))
f_old <- sapply(g_old, function(gm) nom_bounds(gm, t2)[2])
pr_old <- if (p_row < min(f_old)) 0 else NA # old rule: below table -> 0
f_new <- grids[[1]]$fwd[, 3]
ok <- !duplicated(f_new) & f_new > 0
pr_new <- exp(
    approx(log(f_new[ok]), log(grids[[1]]$g[ok]), xout = log(p_row), rule = 1)$y
)
cat(sprintf(
    "allocation %.2e, p = %.0e: old rule p^r = %g -> rejects (wrong: %s); new rule p^r = %.2e -> rejects %s\n",
    a_tiny,
    p_row,
    pr_old,
    pr_old < a_tiny,
    pr_new,
    pr_new < a_tiny
))
set.seed(3)
cat(sprintf(
    "fraction of p < 1e-10 at NCP 4: %.4f\n",
    mean(pnorm(rnorm(1e6, 4), lower.tail = FALSE) < 1e-10)
))

cat("\n== Cost ==\n")
tm <- system.time(
    fw <- do.call(
        rbind,
        lapply(
            exp(seq(log(g_min), log(alpha), length.out = 1024)),
            nom_bounds,
            t = c(1 / 3, 2 / 3, 1)
        )
    )
)
cat(sprintf(
    "boundary table, K=3, 1024 levels: %.2f s (%.2f ms per gsBound1 call)\n",
    tm[["elapsed"]],
    1000 * tm[["elapsed"]] / 1024
))
g <- grids[[1]]
p <- runif(1e6)
ok <- !duplicated(g$fwd[, 3])
tm <- system.time(
    for (j in 1:12) {
        exp(approx(log(g$fwd[ok, 3]), log(g$g[ok]), xout = log(p), rule = 2)$y)
    }
)
cat(sprintf(
    "inverse interpolation, 12 slices x 1e6 values: %.2f s\n",
    tm[["elapsed"]]
))
