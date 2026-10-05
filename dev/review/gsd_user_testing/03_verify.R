# 03_verify.R
# Independent checks on the two optimised designs of 01 and 02.
#
#   A. transform + kernel against an independent reference (00_reference.R):
#      Maurer and Bretz Algorithm 1 on RAW p-values, boundaries by quadrature.
#   B. the compiled gain against the same gain written in plain R.
#   C. the expected gain on fresh draws from an independent simulator.
#   D. graphicalMCP::graph_test_shortcut_gsd() on a subsample.
#   E. familywise error rate of the optimised graphs under null configurations.
#   F. internal consistency of the reported power fields.
#   G. the optimised graph against the simple graphs on fresh draws, as paired
#      differences with standard errors (is the advantage real out of sample?).
#
# Run from the repository root after 01 and 02:
#   Rscript dev/review/gsd_user_testing/03_verify.R

suppressPackageStartupMessages({
    library(multigrain)
    library(gsDesign)
})

here <- "dev/review/gsd_user_testing"
if (!dir.exists(here)) here <- "."
source(file.path(here, "00_reference.R"))

kernel_times <- function(pv, w, G, alpha = pv$alpha) {
    x <- pv$pvals
    dim(x) <- c(pv$nsim, pv$m * pv$K)
    multigrain:::graph_shortcut_gsd(x, alpha, w, G, pv$K)$time
}
ldof_plain <- function(a, t) gsDesign::sfLDOF(a, t)$spend

G2 <- rbind(c(0, 1), c(1, 0))
G3_holm <- matrix(0.5, 3, 3); diag(G3_holm) <- 0
cases <- list(
    two = list(
        file = "01_two_hyp_two_stage_results.rds",
        per_trial = function(tau, value, d) {
            dd <- c(0, d)
            value[1] * dd[tau[, 1] + 1L] + value[2] * dd[tau[, 2] + 1L]
        },
        comparators = list(
            `Holm` = list(w = c(0.5, 0.5), G = G2),
            `fixed sequence H1 -> H2` = list(w = c(1, 0), G = G2),
            `fixed sequence H2 -> H1` = list(w = c(0, 1), G = G2)
        ),
        gain_r = function(tau, value, d) ref_gain_additive(tau, value, d),
        gain_pkg = function(value, d) {
            trial_success_gsd(!!value[1] * d(t1) + !!value[2] * d(t2), d = d,
                              verbose = "silent")
        }
    ),
    three = list(
        file = "02_three_hyp_two_stage_results.rds",
        per_trial = function(tau, value, d) {
            dd <- c(0, d)
            value[1] * dd[tau[, 1] + 1L] +
                (tau[, 1] > 0) * (value[2] * dd[tau[, 2] + 1L] + value[3] * dd[tau[, 3] + 1L])
        },
        comparators = list(
            `Holm` = list(w = rep(1, 3) / 3, G = G3_holm),
            `fixed sequence H1 -> H2 -> H3` =
                list(w = c(1, 0, 0), G = rbind(c(0, 1, 0), c(0, 0, 1), c(1, 0, 0))),
            `H1 gate, then Holm on H2, H3` =
                list(w = c(1, 0, 0), G = rbind(c(0, 0.5, 0.5), c(0, 0, 1), c(0, 1, 0)))
        ),
        gain_r = function(tau, value, d) {
            dd <- c(0, d)
            mean(value[1] * dd[tau[, 1] + 1L] +
                     (tau[, 1] > 0) * (value[2] * dd[tau[, 2] + 1L] +
                                           value[3] * dd[tau[, 3] + 1L]))
        },
        gain_pkg = function(value, d) {
            trial_success_gsd(
                !!value[1] * d(t1) + r1 * (!!value[2] * d(t2) + !!value[3] * d(t3)),
                d = d, verbose = "silent"
            )
        }
    )
)

results <- list()
for (nm in names(cases)) {
    cs  <- cases[[nm]]
    sav <- readRDS(file.path(here, cs$file))
    ds  <- sav$design
    opt <- sav$optimised$discounted
    w <- unname(opt$hyp_weight); G <- unname(opt$trans_matrix)
    m <- length(w); d <- c(1, ds$delta)
    cat("\n==================== ", nm, " hypotheses ====================\n", sep = "")

    # the same draws as the optimisation script
    set.seed(ds$seed)
    raw <- simulate_pvalues_gsd(ds$power_nom, alpha = ds$alpha, corr_matrix = ds$corr,
                                info_frac = ds$info_frac, nsim = ds$nsim)
    pv  <- transform_pvalues_gsd(raw, spending = sfLDOF, alpha = ds$alpha)
    tau_pkg <- kernel_times(pv, w, G)

    # ---- A. independent reference on the raw p-values -----------------------
    t_ref <- system.time(
        tau_ref <- ref_test(unclass(raw), w, G, t1 = ds$info_frac[1],
                            sf = ref_sf_ldof, alpha = ds$alpha)
    )[["elapsed"]]
    n_diff <- sum(tau_pkg != tau_ref)
    cat(sprintf("A. reference vs transform + kernel: %d of %d decision times differ (%d trials; reference took %.0f s)\n",
                n_diff, length(tau_ref), ds$nsim, t_ref))
    if (n_diff > 0) {
        bad <- which(rowSums(tau_pkg != tau_ref) > 0)
        cat("   differing trials:", utils::head(bad, 10), "\n")
        for (b in utils::head(bad, 3)) {
            cat("   trial", b, " raw p:", signif(raw[b, , ], 8),
                " pkg:", tau_pkg[b, ], " ref:", tau_ref[b, ], "\n")
        }
    }

    # ---- B. compiled gain against plain R ------------------------------------
    gain    <- cs$gain_pkg(ds$value, d)
    g_cpp   <- gain$func(tau_pkg)
    g_r     <- cs$gain_r(tau_pkg, ds$value, d)
    cat(sprintf("B. compiled gain %.10f, plain R %.10f, difference %.1e; reported by optimiser %.10f\n",
                g_cpp, g_r, g_cpp - g_r, opt$power$trial_success))

    # ---- C. fresh draws from the independent simulator ----------------------
    set.seed(ds$seed + 1000)
    n_new   <- 2e5
    raw_new <- ref_simulate(ds$power_nom, ds$corr, ds$info_frac[1], n_new, ds$alpha)
    tau_new <- ref_test(raw_new, w, G, t1 = ds$info_frac[1], sf = ref_sf_ldof,
                        alpha = ds$alpha)
    per_trial <- if (nm == "two") {
        rowSums(sapply(1:2, function(i) ds$value[i] * c(0, d)[tau_new[, i] + 1L]))
    } else {
        dd <- c(0, d)
        ds$value[1] * dd[tau_new[, 1] + 1L] +
            (tau_new[, 1] > 0) * (ds$value[2] * dd[tau_new[, 2] + 1L] +
                                      ds$value[3] * dd[tau_new[, 3] + 1L])
    }
    se_new <- stats::sd(per_trial) / sqrt(n_new)
    se_in  <- stats::sd(per_trial) / sqrt(ds$nsim)
    z <- (g_cpp - mean(per_trial)) / sqrt(se_new^2 + se_in^2)
    cat(sprintf("C. independent simulator + reference, %d fresh trials: gain %.5f (s.e. %.5f); package in-sample %.5f; z = %.2f\n",
                n_new, mean(per_trial), se_new, g_cpp, z))
    # the same fresh draws through the package (transform + kernel)
    attr(raw_new, "info_frac") <- matrix(ds$info_frac, m, 2, byrow = TRUE)
    pv_new   <- transform_pvalues_gsd(raw_new, spending = sfLDOF, alpha = ds$alpha)
    tau_new2 <- kernel_times(pv_new, w, G)
    cat(sprintf("   same fresh draws through transform + kernel: %d of %d decision times differ\n",
                sum(tau_new2 != tau_new), length(tau_new)))

    # ---- D. graphicalMCP oracle on a subsample ------------------------------
    if (requireNamespace("graphicalMCP", quietly = TRUE)) {
        set.seed(7)
        idx   <- sample.int(ds$nsim, 300)
        graph <- graphicalMCP::graph_create(w, G)
        n_rej_diff <- 0L; n_time_diff <- 0L; n_rej <- 0L
        for (s in idx) {
            p_s <- matrix(raw[s, , ], m, 2, dimnames = list(names(graph$hypotheses), NULL))
            o <- suppressWarnings(graphicalMCP::graph_test_shortcut_gsd(
                graph = graph, p = p_s, alpha = ds$alpha,
                info_frac = ds$info_frac, spending_fn = ldof_plain,
                look_back = FALSE
            ))
            rej <- unname(o$outputs$rejected)
            n_rej_diff  <- n_rej_diff + sum(rej != (tau_pkg[s, ] > 0))
            n_time_diff <- n_time_diff +
                sum(as.integer(unname(o$outputs$decision_at))[rej] != tau_pkg[s, rej])
            n_rej <- n_rej + sum(rej)
        }
        cat(sprintf("D. graphicalMCP on 300 trials: %d rejection disagreements of %d; %d decision-time disagreements of %d rejections\n",
                    n_rej_diff, 300L * m, n_time_diff, n_rej))
    } else {
        cat("D. graphicalMCP not installed; skipped\n")
    }

    # ---- E. familywise error rate under null configurations -----------------
    # A hypothesis is made true by giving it nominal power alpha (zero
    # non-centrality). Every non-empty set of true nulls is tried.
    n_null <- 1e6
    fwer <- list()
    for (code in seq_len(2^m - 1)) {
        is_null <- as.logical(bitwAnd(code, 2^(seq_len(m) - 1)))
        pw <- ifelse(is_null, ds$alpha, ds$power_nom)
        set.seed(500 + code)
        raw0 <- simulate_pvalues_gsd(pw, alpha = ds$alpha, corr_matrix = ds$corr,
                                     info_frac = ds$info_frac, nsim = n_null)
        pv0  <- transform_pvalues_gsd(raw0, spending = sfLDOF, alpha = ds$alpha)
        t0   <- kernel_times(pv0, w, G)
        fwer[[paste(which(is_null), collapse = ",")]] <-
            mean(rowSums(t0[, is_null, drop = FALSE] > 0) > 0)
    }
    se0 <- sqrt(ds$alpha * (1 - ds$alpha) / n_null)
    cat(sprintf("E. FWER of the optimised graph, %d trials per configuration (s.e. about %.5f):\n",
                n_null, se0))
    for (k in names(fwer)) cat(sprintf("   true nulls {%s}: %.5f\n", k, fwer[[k]]))

    # ---- F. internal consistency of the power fields -------------------------
    pw <- calc_power_pvals_gsd(pv, hyp_weight = w, trans_matrix = G,
                               custom_power = list(gain = gain))
    td <- pw$time_distribution
    checks <- c(
        `time_distribution rows sum to 1` =
            isTRUE(all.equal(unname(rowSums(td)), rep(1, m))),
        `local_power equals last column of local_power_by_analysis` =
            isTRUE(all.equal(unname(pw$local_power),
                             unname(pw$local_power_by_analysis[, 2]))),
        `time_distribution matches kernel times` =
            isTRUE(all.equal(unname(td[, 2]), colMeans(tau_pkg == 1))) &&
            isTRUE(all.equal(unname(td[, 3]), colMeans(tau_pkg == 2))),
        `mean_decision_look matches kernel times` =
            isTRUE(all.equal(unname(pw$mean_decision_look),
                             sapply(seq_len(m), function(i) mean(tau_pkg[tau_pkg[, i] > 0, i])))),
        `custom_power gain equals optimiser gain` =
            isTRUE(all.equal(pw$gain, opt$power$trial_success))
    )
    cat("F. consistency checks:\n")
    for (k in names(checks)) cat(sprintf("   %-60s %s\n", k, if (checks[[k]]) "ok" else "FAIL"))

    # ---- G. advantage over simple graphs on fresh draws ----------------------
    set.seed(ds$seed + 2000)
    n_g   <- 1e6
    raw_g <- simulate_pvalues_gsd(ds$power_nom, alpha = ds$alpha, corr_matrix = ds$corr,
                                  info_frac = ds$info_frac, nsim = n_g)
    pv_g  <- transform_pvalues_gsd(raw_g, spending = sfLDOF, alpha = ds$alpha)
    g_opt <- cs$per_trial(kernel_times(pv_g, w, G), ds$value, d)
    cat(sprintf("G. fresh package draws, %d trials: optimised graph gain %.5f (s.e. %.5f); in-sample %.5f\n",
                n_g, mean(g_opt), stats::sd(g_opt) / sqrt(n_g), g_cpp))
    paired <- do.call(rbind, lapply(names(cs$comparators), function(cn) {
        cg <- cs$comparators[[cn]]
        g_c <- cs$per_trial(kernel_times(pv_g, cg$w, cg$G), ds$value, d)
        data.frame(comparator = cn, gain = mean(g_c),
                   optimised_minus_comparator = mean(g_opt - g_c),
                   se = stats::sd(g_opt - g_c) / sqrt(n_g))
    }))
    paired$z <- paired$optimised_minus_comparator / paired$se
    print(paired, digits = 4, row.names = FALSE)

    results[[nm]] <- list(
        paired = paired, fresh_package_gain = mean(g_opt),
        n_time_diff_reference = n_diff, gain_cpp = g_cpp, gain_r = g_r,
        fresh_gain = mean(per_trial), fresh_se = se_new, z = z,
        fwer = fwer, checks = checks
    )
}

saveRDS(results, file.path(here, "03_verify_results.rds"))
