# Helpers for the efficacy stopping-rule tests of `graph_shortcut_gsd()`.
#
# `gsd_apply_stopping()` is the `apply_stopping()` post-step of the script
# 04_stopping_rules.R in dev/review/gsd_user_testing, generalised to K
# analyses and one stopping condition per interim: an efficacy stopping rule
# is a deterministic edit of the always-continue decision-time matrix (README
# of that directory, section 7.3).
#
# The `gsd_ref_*()` functions are verbatim copies of `ref_sf_ldof()`,
# `ref_simulate()`, `ref_bounds()` and `ref_test()` from
# `dev/review/gsd_user_testing/00_reference.R`, renamed with a `gsd_` prefix.
# They are copied rather than sourced because `dev/` is build-ignored, so the
# installed tests cannot reach it. They form an independent K = 2 reference
# (no gsDesign, no transform, no C++ kernel) that really stops analysing a
# trial when `stop_fun` fires.

# `stop_when` is a list of K - 1 functions `function(tau, k)`, one per interim,
# each returning a logical per trial: TRUE stops the trial after look k given
# the decision times so far. Each must look only at `tau > 0 & tau <= k`.
# Returns the edited decision times and the look S at which each trial stops
# under the rule (K if it never fires). The all-rejected stop of the kernel is
# not applied here.
gsd_apply_stopping <- function(
    tau,
    stop_when,
    K # nolint: object_name_linter.
) {
    stop_look <- rep(K, nrow(tau))
    for (k in rev(seq_len(K - 1L))) {
        stop_look[stop_when[[k]](tau, k)] <- k
    }
    tau[tau > stop_look] <- 0L # nothing can be rejected after the trial stops
    list(time = tau, stop_look = stop_look)
}


# ---- verbatim copies from dev/review/gsd_user_testing/00_reference.R -------
# nolint start: line_length_linter, object_name_linter, commas_linter, cyclocomp_linter.

gsd_ref_sf_ldof <- function(a, t) {
    2 *
        stats::pnorm(
            stats::qnorm(1 - a / 2) / sqrt(pmin(t, 1)),
            lower.tail = FALSE
        )
}

# Z_{i,1} = sqrt(t1) * Delta_i + e_{i,1},
# Z_{i,2} = sqrt(t1) * Z_{i,1} + sqrt(1 - t1) * (sqrt(1 - t1) * Delta_i + e_{i,2}),
# with e_{.,1}, e_{.,2} independent N(0, corr). This gives the canonical joint
# distribution for a common interim information fraction t1 and final 1.
gsd_ref_simulate <- function(power_nom, corr, t1, nsim, alpha = 0.025) {
    m <- length(power_nom)
    delta <- stats::qnorm(1 - alpha) + stats::qnorm(power_nom)
    L <- chol(corr)
    e1 <- matrix(stats::rnorm(nsim * m), nsim, m) %*% L
    e2 <- matrix(stats::rnorm(nsim * m), nsim, m) %*% L
    mu <- matrix(delta, nsim, m, byrow = TRUE)
    z1 <- sqrt(t1) * mu + e1
    z2 <- sqrt(t1) * z1 + sqrt(1 - t1) * (sqrt(1 - t1) * mu + e2)
    out <- array(NA_real_, dim = c(nsim, m, 2L))
    out[,, 1] <- stats::pnorm(z1, lower.tail = FALSE)
    out[,, 2] <- stats::pnorm(z2, lower.tail = FALSE)
    out
}

# Look 1: P(Z1 > c1) = spend(t1). Look 2: P(Z1 <= c1, Z2 > c2) = spend(1) -
# spend(t1), with Corr(Z1, Z2) = sqrt(t1). The bivariate probability is the
# one-dimensional integral of phi(z) * P(Z2 > c2 | Z1 = z) over z <= c1.
#
# Two special cases need no integration. An endpoint with no data at the
# interim (t1 = NA) cannot be rejected there and is tested once, at its full
# level, at the final analysis. An endpoint that is complete at the interim
# (t1 >= 1) is tested at its full level at both analyses on the same p-value.
gsd_ref_bounds <- function(a, t1, sf) {
    if (a <= 0) {
        return(c(0, 0))
    }
    if (is.na(t1)) {
        return(c(0, a))
    }
    if (t1 >= 1) {
        return(c(a, a))
    }
    cum <- sf(a, c(t1, 1))
    c1 <- stats::qnorm(cum[1], lower.tail = FALSE)
    rho <- sqrt(t1)
    cross2 <- function(c2) {
        stats::integrate(
            function(z) {
                stats::dnorm(z) *
                    stats::pnorm(
                        (c2 - rho * z) / sqrt(1 - rho^2),
                        lower.tail = FALSE
                    )
            },
            lower = -Inf,
            upper = c1,
            rel.tol = 1e-12,
            abs.tol = 0,
            subdivisions = 500L
        )$value
    }
    target <- cum[2] - cum[1]
    c2 <- stats::uniroot(
        function(x) cross2(x) - target,
        c(0, 10),
        tol = 1e-13
    )$root
    stats::pnorm(c(c1, c2), lower.tail = FALSE)
}

# raw: nsim x m x 2 array. t1: the interim information fraction, one value or
# one per hypothesis (NA = no data at the interim). sf: one spending function
# or a list of m.
# Returns the nsim x m integer matrix of decision times (0 = never).
# Boundaries are cached by (hypothesis, level): for a fixed graph the level a
# hypothesis can hold takes only a handful of values.
#
# stop_fun, if given, is a function of the decision times so far (a length-m
# integer vector, 0 = not yet rejected) that returns TRUE when the trial stops
# after the current analysis. A stopped trial is not analysed again, so nothing
# can be rejected at a later analysis.
#
# stop_after_1, if given, is a logical vector with one entry per trial: TRUE
# stops that trial after the interim whatever was rejected (a futility rule
# evaluated on the raw data beforehand).
gsd_ref_test <- function(
    raw,
    w,
    G,
    t1,
    sf,
    alpha = 0.025,
    strict = TRUE,
    stop_fun = NULL,
    stop_after_1 = NULL
) {
    nsim <- dim(raw)[1]
    m <- dim(raw)[2]
    if (is.function(sf)) {
        sf <- rep(list(sf), m)
    }
    t1 <- rep_len(t1, m)
    cache <- vector("list", m)
    for (i in seq_len(m)) {
        cache[[i]] <- new.env(parent = emptyenv())
    }
    bound <- function(i, a) {
        key <- format(a, digits = 17)
        b <- cache[[i]][[key]]
        if (is.null(b)) {
            b <- gsd_ref_bounds(a, t1[i], sf[[i]])
            assign(key, b, envir = cache[[i]])
        }
        b
    }
    below <- if (strict) `<` else `<=`
    tau <- matrix(0L, nsim, m)
    for (s in seq_len(nsim)) {
        a <- w * alpha
        g <- G
        live <- rep(TRUE, m)
        for (k in 1:2) {
            repeat {
                j <- 0L
                for (i in which(live)) {
                    if (
                        a[i] > 0 &&
                            !is.na(raw[s, i, k]) &&
                            below(raw[s, i, k], bound(i, a[i])[k])
                    ) {
                        j <- i
                        break
                    }
                }
                if (j == 0L) {
                    break
                }
                tau[s, j] <- k
                live[j] <- FALSE
                a_new <- a + a[j] * g[j, ]
                a_new[j] <- 0
                g_new <- matrix(0, m, m)
                for (l in which(live)) {
                    for (q in which(live)) {
                        if (l != q) {
                            den <- 1 - g[l, j] * g[j, l]
                            if (den > 0) {
                                g_new[l, q] <- (g[l, q] + g[l, j] * g[j, q]) /
                                    den
                            }
                        }
                    }
                }
                a <- a_new
                g <- g_new
            }
            if (!is.null(stop_fun) && stop_fun(tau[s, ])) {
                break
            }
            if (k == 1L && !is.null(stop_after_1) && stop_after_1[s]) break
        }
    }
    tau
}
# nolint end
