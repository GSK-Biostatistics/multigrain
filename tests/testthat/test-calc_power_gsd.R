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

cp_time <- function(pvals, w = cp_w, trans = cp_g, alpha = pvals$alpha) {
    graph_shortcut_gsd(
        pvals = .gsd_kernel_matrix(pvals),
        alpha = alpha,
        w = w,
        G = trans,
        K = pvals$K
    )
}


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

    expect_identical(out$gain, cp_gain$func(res$time))
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


test_that("a plain function receives the decision-time row", {
    res <- cp_time(cp_pvals_gsd)
    early <- function(x) as.numeric(x[1L] == 1L)

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
            "time_distribution", "func1"
        )
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
