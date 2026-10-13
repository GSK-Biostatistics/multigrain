# `calc_power_pvals_gsd()` summarises both outputs of the group sequential
# kernel: the rejection indicators and the decision times (design record,
# section 4.7).

cp_pvals_gsd <- withr::with_seed(17, {
    raw <- simulate_pvalues_gsd(
        power_nominal = c(0.9, 0.85, 0.8),
        corr_matrix = diag(3),
        info_frac = c(0.5, 1),
        nsim = 2000
    )
    transform_pvalues_gsd(
        raw,
        spending = gsDesign::sfLDOF,
        grid_size = 128L
    )
})

cp_pvals_k1 <- withr::with_seed(19, {
    raw <- simulate_pvalues_gsd(
        power_nominal = c(0.9, 0.85, 0.8),
        corr_matrix = diag(3),
        info_frac = 1,
        nsim = 2000
    )
    transform_pvalues_gsd(
        raw,
        spending = gsDesign::sfLDOF,
        grid_size = 128L
    )
})

cp_w <- c(0.4, 0.3, 0.3)
cp_g <- matrix(
    c(
        0, 0.5, 0.5,
        0.5, 0, 0.5,
        0.5, 0.5, 0
    ),
    nrow = 3,
    byrow = TRUE
)

cp_gain <- trial_success_gsd(
    d(t1) + d(t2) + d(t3),
    d = c(1, 0.75),
    verbose = "silent"
)

cp_time <- function(
    pvals,
    w = cp_w,
    trans = cp_g,
    alpha = pvals$alpha,
    stop_rule = NULL
) {
    graph_shortcut_gsd(
        pvals = .gsd_kernel_matrix(pvals),
        alpha = alpha,
        w = w,
        G = trans,
        K = pvals$K,
        stop_rule = stop_rule
    )
}

# The user's rule: stop at the interim once H1 and one of H2 and H3 are
# rejected. `cp_stop_when` is the same condition written directly on the
# decision times, for `gsd_apply_stopping()`, independently of
# `.gsd_stop_table()`.
cp_rule <- stopping_rule_gsd(r1 && (r2 || r3))
cp_table <- .gsd_stop_table(cp_rule$conditions, 3L, 2L)
cp_stop_when <- list(
    function(tau, k) {
        rej <- tau > 0L & tau <= k
        rej[, 1L] & (rej[, 2L] | rej[, 3L])
    }
)


test_that("every element has the documented shape and dimnames", {
    out <- calc_power_pvals_gsd(cp_pvals_gsd, cp_w, cp_g)

    expect_length(out$local_power, 3L)
    expect_null(names(out$local_power))

    expect_identical(dim(out$local_power_by_analysis), c(3L, 2L))
    expect_identical(
        dimnames(out$local_power_by_analysis),
        list(c("H1", "H2", "H3"), c("analysis 1", "analysis 2"))
    )

    expect_length(out$exp_rejections, 1L)
    expect_length(out$disj_power, 1L)
    expect_length(out$conj_power, 1L)

    expect_length(out$mean_decision_look, 3L)

    expect_identical(dim(out$time_distribution), c(3L, 3L))
    expect_identical(
        dimnames(out$time_distribution),
        list(c("H1", "H2", "H3"), c("never", "analysis 1", "analysis 2"))
    )

    expect_type(out$stop_distribution, "double")
    expect_named(out$stop_distribution, c("analysis 1", "analysis 2"))
    expect_length(out$mean_stop_look, 1L)
})


test_that("the stopping fields sum to one and agree with the kernel", {
    no_rule <- calc_power_pvals_gsd(cp_pvals_gsd, cp_w, cp_g)
    res <- cp_time(cp_pvals_gsd)

    # sums of proportions: equality holds to floating point, not bitwise
    # nolint start: expect_identical_linter
    expect_equal(sum(no_rule$stop_distribution), 1)
    expect_equal(
        no_rule$mean_stop_look,
        sum(seq_len(2L) * no_rule$stop_distribution)
    )
    # nolint end

    # without a rule a trial stops early only once it has rejected every
    # hypothesis
    all_rejected <- rowSums(res$time > 0L) == 3L
    expected_look <- ifelse(all_rejected, apply(res$time, 1L, max), 2L)
    expect_identical(res$stop_look, as.integer(expected_look))
    expect_identical(
        no_rule$stop_distribution,
        c(
            "analysis 1" = mean(res$stop_look == 1L),
            "analysis 2" = mean(res$stop_look == 2L)
        )
    )
    expect_identical(no_rule$mean_stop_look, mean(res$stop_look))
})


test_that("with a stopping rule every field is the stopped trials'", {
    no_rule <- cp_time(cp_pvals_gsd)
    edit <- gsd_apply_stopping(no_rule$time, cp_stop_when, 2L)
    stopped <- cp_time(cp_pvals_gsd, stop_rule = cp_table)

    out <- calc_power_pvals_gsd(
        cp_pvals_gsd,
        cp_w,
        cp_g,
        stopping_rule = cp_rule
    )

    # not vacuous: the rule fires and removes final-analysis rejections
    expect_gt(out$stop_distribution[["analysis 1"]], 0)
    expect_true(any(no_rule$time != edit$time))

    expect_identical(out$local_power, colMeans(edit$time > 0L))
    expect_identical(
        out$stop_distribution,
        c(
            "analysis 1" = mean(stopped$stop_look == 1L),
            "analysis 2" = mean(stopped$stop_look == 2L)
        )
    )
    expect_identical(out$mean_stop_look, mean(stopped$stop_look))

    by_hand_time <- matrix(NA_real_, nrow = 3L, ncol = 3L)
    by_hand_power <- matrix(NA_real_, nrow = 3L, ncol = 2L)
    for (i in seq_len(3L)) {
        tau <- edit$time[, i]
        for (k in 0:2) {
            by_hand_time[i, k + 1L] <- mean(tau == k)
        }
        for (k in 1:2) {
            by_hand_power[i, k] <- mean(tau >= 1L & tau <= k)
        }
    }
    expect_identical(unname(out$time_distribution), by_hand_time)
    expect_identical(unname(out$local_power_by_analysis), by_hand_power)
})


test_that("every kind of custom measure sees the stopped trials", {
    stopped <- cp_time(cp_pvals_gsd, stop_rule = cp_table)

    gain_s <- trial_success_gsd(
        d(t1) + d(t2) + d(t3) + 0.1 * (s == 1),
        d = c(1, 0.8),
        verbose = "silent"
    )
    # a fixed-sample gain must be defined over all three hypotheses
    fixed_gain <- trial_success(r1 + r2 + r3, verbose = "silent")
    early_stop <- function(x, s) as.numeric(s == 1L)

    out <- calc_power_pvals_gsd(
        cp_pvals_gsd,
        cp_w,
        cp_g,
        custom_power = list(
            gain_s = gain_s,
            fixed = fixed_gain,
            early_stop = early_stop
        ),
        stopping_rule = cp_rule
    )

    disc <- c(1, 0.8)
    tau <- stopped$time
    per_trial <- rowSums(ifelse(tau > 0L, disc[pmax(tau, 1L)], 0)) +
        0.1 * (stopped$stop_look == 1L)
    expect_equal(out$gain_s, mean(per_trial), tolerance = 1e-12)
    expect_identical(out$gain_s, gain_s$func(tau, stopped$stop_look))

    expect_identical(out$fixed, fixed_gain$func(stopped$rejected))

    expect_identical(
        out$early_stop,
        out$stop_distribution[["analysis 1"]]
    )
})


test_that("d(s) is the discount of the stopping look", {
    # `d(s)` alone refers to no hypothesis; the zero term fixes m = 3
    gain_ds <- trial_success_gsd(
        d(s) + 0 * (r1 + r2 + r3),
        d = c(1, 0.75),
        verbose = "silent"
    )
    disc <- c(1, 0.75)

    no_rule <- calc_power_pvals_gsd(
        cp_pvals_gsd,
        cp_w,
        cp_g,
        custom_power = list(gain = gain_ds)
    )
    ruled <- calc_power_pvals_gsd(
        cp_pvals_gsd,
        cp_w,
        cp_g,
        custom_power = list(gain = gain_ds),
        stopping_rule = cp_rule
    )

    # a mean against a weighted sum of proportions: equal to floating point
    # nolint start: expect_identical_linter
    expect_equal(no_rule$gain, sum(disc * no_rule$stop_distribution))
    expect_equal(ruled$gain, sum(disc * ruled$stop_distribution))
    # nolint end
    expect_gt(ruled$gain, no_rule$gain)
})


test_that("a graph over the wrong number of hypotheses aborts in R", {
    # before this check the mismatch surfaced only as the kernel's C++ error
    # "G must be an m x m matrix"
    expect_error(
        calc_power_pvals_gsd(cp_pvals_gsd, c(0.5, 0.5), cp_g),
        "`hyp_weight` describes m = 2 hypotheses, but the p-values hold m = 3",
        fixed = TRUE
    )
    expect_error(
        calc_power_pvals_gsd(
            cp_pvals_gsd,
            cp_w,
            matrix(c(0, 1, 1, 0), nrow = 2)
        ),
        paste(
            "`trans_matrix` describes m = 2 hypotheses, but the p-values",
            "hold m = 3"
        ),
        fixed = TRUE
    )
    expect_error(
        calc_power_pvals_gsd(cp_pvals_gsd, cp_w, cp_g[, 1:2]),
        "`trans_matrix` describes m = 2 hypotheses",
        fixed = TRUE
    )
})


test_that("the stopping rule is validated against the p-values", {
    expect_error(
        calc_power_pvals_gsd(
            cp_pvals_gsd,
            cp_w,
            cp_g,
            stopping_rule = quote(r1)
        ),
        "multigrain_stopping_rule_gsd"
    )
    expect_error(
        calc_power_pvals_gsd(
            cp_pvals_gsd,
            cp_w,
            cp_g,
            stopping_rule = stopping_rule_gsd(r1 || r4)
        ),
        "m = 3"
    )
})


test_that("the summaries are internally consistent", {
    out <- calc_power_pvals_gsd(cp_pvals_gsd, cp_w, cp_g)

    # sums of proportions: equality holds to floating point, not bitwise
    # nolint start: expect_identical_linter
    expect_equal(
        out$local_power,
        unname(out$local_power_by_analysis[, 2L])
    )
    expect_equal(rowSums(out$time_distribution), c(H1 = 1, H2 = 1, H3 = 1))
    expect_equal(
        unname(1 - out$time_distribution[, "never"]),
        out$local_power
    )
    expect_equal(out$exp_rejections, sum(out$local_power))
    # nolint end
})


test_that("mean_decision_look matches a hand computation", {
    res <- cp_time(cp_pvals_gsd)
    out <- calc_power_pvals_gsd(cp_pvals_gsd, cp_w, cp_g)

    by_hand <- vapply(
        seq_len(3L),
        function(i) {
            tau <- res$time[, i]
            mean(tau[tau > 0L])
        },
        numeric(1L)
    )

    expect_identical(out$mean_decision_look, by_hand)
})


test_that("mean_decision_look is NA for an unreachable hypothesis", {
    w <- c(0.5, 0.5, 0)
    trans <- rbind(
        c(0, 1, 0),
        c(1, 0, 0),
        c(1, 0, 0)
    )

    out <- calc_power_pvals_gsd(cp_pvals_gsd, w, trans)

    expect_identical(out$local_power[[3L]], 0)
    expect_identical(out$mean_decision_look[[3L]], NA_real_)
    expect_identical(out$time_distribution[3L, "never"], 1)
})


test_that("at K = 1 the shared fields equal calc_power_pvals()", {
    gsd <- calc_power_pvals_gsd(cp_pvals_k1, cp_w, cp_g)
    fixed <- calc_power_pvals(cp_pvals_k1$pvals[, , 1L], cp_w, cp_g)

    expect_identical(gsd$local_power, fixed$local_power)
    expect_identical(gsd$exp_rejections, fixed$exp_rejections)
    expect_identical(gsd$disj_power, fixed$disj_power)
    expect_identical(gsd$conj_power, fixed$conj_power)
})


test_that("a GSD gain is scored on the decision times", {
    res <- cp_time(cp_pvals_gsd)

    out <- calc_power_pvals_gsd(
        cp_pvals_gsd,
        cp_w,
        cp_g,
        custom_power = list(gain = cp_gain)
    )

    expect_identical(out$gain, cp_gain$func(res$time, res$stop_look))
})


test_that("a fixed-sample gain is scored on the rejection indicators", {
    fixed_gain <- trial_success(r1 + r2 + r3, verbose = "silent")

    # at K = 1 it must agree with calc_power_pvals() on the same matrix
    gsd_k1 <- calc_power_pvals_gsd(
        cp_pvals_k1,
        cp_w,
        cp_g,
        custom_power = list(total = fixed_gain)
    )
    fixed_k1 <- calc_power_pvals(
        cp_pvals_k1$pvals[, , 1L],
        cp_w,
        cp_g,
        custom_power = list(total = fixed_gain)
    )
    expect_identical(gsd_k1$total, fixed_k1$total)

    # at K = 2 it is the compiled function on the kernel's `rejected` matrix
    res <- cp_time(cp_pvals_gsd)
    gsd_k2 <- calc_power_pvals_gsd(
        cp_pvals_gsd,
        cp_w,
        cp_g,
        custom_power = list(total = fixed_gain)
    )
    expect_identical(gsd_k2$total, fixed_gain$func(res$rejected))
})


test_that("a plain function receives the decision-time row and the look", {
    res <- cp_time(cp_pvals_gsd)
    early <- function(x, s) as.numeric(x[1L] == 1L)

    out <- calc_power_pvals_gsd(
        cp_pvals_gsd,
        cp_w,
        cp_g,
        custom_power = list(early = early)
    )

    expect_identical(out$early, mean(res$time[, 1L] == 1L))

    # an unnamed measure is reported as funcN
    unnamed <- calc_power_pvals_gsd(
        cp_pvals_gsd,
        cp_w,
        cp_g,
        custom_power = list(early)
    )
    expect_named(
        unnamed,
        c(
            "local_power", "local_power_by_analysis", "exp_rejections",
            "disj_power", "conj_power", "mean_decision_look",
            "time_distribution", "stop_distribution", "mean_stop_look",
            "func1"
        )
    )

    # a one-argument function is refused: every plain function receives the
    # stopping look as its second argument
    expect_error(
        calc_power_pvals_gsd(
            cp_pvals_gsd,
            cp_w,
            cp_g,
            custom_power = list(old = function(x) x[1L])
        ),
        "unused argument"
    )
})


test_that("a gain with the wrong dimensions aborts", {
    wrong_m <- trial_success_gsd(
        d(t1) + d(t2),
        d = c(1, 0.75),
        verbose = "silent"
    )
    wrong_k <- trial_success_gsd(
        d(t1) + d(t2) + d(t3),
        d = c(1, 0.9, 0.75),
        verbose = "silent"
    )

    expect_error(
        calc_power_pvals_gsd(
            cp_pvals_gsd,
            cp_w,
            cp_g,
            custom_power = list(gain = wrong_m)
        ),
        "number of hypotheses"
    )
    expect_error(
        calc_power_pvals_gsd(
            cp_pvals_gsd,
            cp_w,
            cp_g,
            custom_power = list(gain = wrong_k)
        ),
        "number of analyses"
    )
})


test_that("a fixed-sample gain with the wrong m aborts", {
    # its compiled function indexes the `rejected` matrix by hypothesis
    # without a bounds check, so a gain over more hypotheses than the graph
    # has would read past the row
    wide <- trial_success(r1 + r2 + r3, verbose = "silent")
    pvals_2 <- withr::with_seed(29, {
        raw <- simulate_pvalues_gsd(
            power_nominal = c(0.9, 0.8),
            corr_matrix = diag(2),
            info_frac = c(0.5, 1),
            nsim = 500
        )
        transform_pvalues_gsd(
            raw,
            spending = gsDesign::sfLDOF,
            grid_size = 64L
        )
    })
    w2 <- c(0.5, 0.5)
    g2 <- matrix(c(0, 1, 1, 0), nrow = 2)

    expect_error(
        calc_power_pvals_gsd(
            pvals_2,
            w2,
            g2,
            custom_power = list(total = wide)
        ),
        "custom_power$total",
        fixed = TRUE
    )
    expect_error(
        calc_power_pvals_gsd(
            pvals_2,
            w2,
            g2,
            custom_power = list(total = wide)
        ),
        "number of hypotheses"
    )

    # a matching gain is unaffected
    narrow <- trial_success(r1 + r2, verbose = "silent")
    out <- calc_power_pvals_gsd(
        pvals_2,
        w2,
        g2,
        custom_power = list(total = narrow)
    )
    expect_type(out$total, "double")
    expect_length(out$total, 1L)
})


test_that("the same graph twice gives identical output", {
    a <- calc_power_pvals_gsd(
        cp_pvals_gsd,
        cp_w,
        cp_g,
        custom_power = list(gain = cp_gain)
    )
    b <- calc_power_pvals_gsd(
        cp_pvals_gsd,
        cp_w,
        cp_g,
        custom_power = list(gain = cp_gain)
    )

    expect_identical(a, b)
})


test_that("alpha is bounded by the level of the p-value object", {
    lower <- calc_power_pvals_gsd(cp_pvals_gsd, cp_w, cp_g, alpha = 0.01)
    full <- calc_power_pvals_gsd(cp_pvals_gsd, cp_w, cp_g)

    expect_true(all(lower$local_power <= full$local_power))

    expect_error(
        calc_power_pvals_gsd(cp_pvals_gsd, cp_w, cp_g, alpha = 0.05),
        "0.025"
    )
})


test_that("a plain matrix is not accepted as pvals", {
    expect_error(
        calc_power_pvals_gsd(cp_pvals_gsd$pvals[, , 1L], cp_w, cp_g),
        "multigrain_pvals_gsd"
    )
})
