# Direct tests of the shared validation helpers in `R/check_gsd.R`. Each runs
# once on the public path of the group sequential consumers, before the
# kernel is called.

chk_pvals_gsd <- withr::with_seed(61, {
    raw <- simulate_pvalues_gsd(
        power_nominal = c(0.9, 0.85, 0.8),
        corr_matrix = diag(3),
        info_frac = c(0.5, 1),
        nsim = 200
    )
    transform_pvalues_gsd(
        raw,
        spending = gsDesign::sfLDOF,
        grid_size = 64L
    )
})

chk_gain_gsd <- trial_success_gsd(
    d(t1) + d(t2) + d(t3),
    d = c(1, 0.75),
    verbose = "silent"
)


test_that("is_pvals_gsd() recognises the class only", {
    expect_true(is_pvals_gsd(chk_pvals_gsd))
    expect_false(is_pvals_gsd(chk_pvals_gsd$pvals))
    expect_false(is_pvals_gsd(unclass(chk_pvals_gsd)))
    expect_false(is_pvals_gsd(NULL))
})


test_that("check_pvals_gsd() passes the class and refuses anything else", {
    expect_null(check_pvals_gsd(chk_pvals_gsd))
    expect_invisible(check_pvals_gsd(chk_pvals_gsd))

    expect_null(check_pvals_gsd(NULL, allow_null = TRUE))
    expect_error(check_pvals_gsd(NULL), "multigrain_pvals_gsd")
    expect_error(check_pvals_gsd(NULL), "not `NULL`", fixed = TRUE)

    plain <- chk_pvals_gsd$pvals
    expect_error(check_pvals_gsd(plain), "`plain`", fixed = TRUE)
    expect_error(check_pvals_gsd(plain), "transform_pvalues_gsd()")
    expect_error(check_pvals_gsd(list()), "multigrain_pvals_gsd")

    # a missing argument of the calling function is reported as absent
    consumer <- function(pvals) check_pvals_gsd(pvals)
    expect_error(consumer(), "multigrain_pvals_gsd")
    expect_error(consumer(), "absent")
})


test_that("check_trial_success_gsd() refuses a fixed-sample gain", {
    expect_null(check_trial_success_gsd(chk_gain_gsd))

    fixed_gain <- trial_success(r1 + r2 + r3, verbose = "silent")
    expect_error(
        check_trial_success_gsd(fixed_gain),
        "scores rejection indicators only"
    )
    expect_error(
        check_trial_success_gsd(fixed_gain),
        "`fixed_gain`",
        fixed = TRUE
    )

    expect_error(
        check_trial_success_gsd(function(x) 1),
        "multigrain_trial_success_gsd"
    )
    expect_error(check_trial_success_gsd(NULL), "multigrain_trial_success_gsd")
    consumer <- function(trial_success) {
        check_trial_success_gsd(trial_success)
    }
    expect_error(consumer(), "multigrain_trial_success_gsd")
    expect_error(consumer(), "absent")
})


test_that(".gsd_check_gain_dims() checks m and, when fixed, K", {
    expect_null(.gsd_check_gain_dims(chk_gain_gsd, chk_pvals_gsd))

    wrong_m <- trial_success_gsd(
        d(t1) + d(t2),
        d = c(1, 0.75),
        verbose = "silent"
    )
    expect_error(
        .gsd_check_gain_dims(wrong_m, chk_pvals_gsd),
        "defined over m = 2, but the p-values hold m = 3"
    )

    wrong_k <- trial_success_gsd(
        d(t1) + d(t2) + d(t3),
        d = c(1, 0.9, 0.75),
        verbose = "silent"
    )
    expect_error(
        .gsd_check_gain_dims(wrong_k, chk_pvals_gsd),
        "defined over K = 3, but the p-values hold K = 2"
    )
    expect_error(
        .gsd_check_gain_dims(wrong_k, chk_pvals_gsd),
        "Its discount table (`d`) has length 3.",
        fixed = TRUE
    )

    # a gain with no discount table fixes no K, so that branch is skipped
    no_k <- trial_success_gsd(r1 + r2 + r3, verbose = "silent")
    expect_null(no_k$K)
    expect_null(.gsd_check_gain_dims(no_k, chk_pvals_gsd))
})


test_that(".gsd_check_alpha() resolves and bounds the level", {
    expect_identical(.gsd_check_alpha(NULL, chk_pvals_gsd), 0.025)
    expect_identical(.gsd_check_alpha(0.025, chk_pvals_gsd), 0.025)
    expect_identical(.gsd_check_alpha(0.01, chk_pvals_gsd), 0.01)

    expect_error(.gsd_check_alpha(0.05, chk_pvals_gsd), "0.025")
    expect_error(.gsd_check_alpha(0, chk_pvals_gsd), "greater than 0")
    expect_error(.gsd_check_alpha("0.01", chk_pvals_gsd), "number")
})


test_that(".gsd_kernel_matrix() reshapes analysis by analysis", {
    # value 100 * k + 10 * i + n: trial n, hypothesis i, analysis k
    values <- array(NA_real_, dim = c(4L, 2L, 2L))
    for (n in 1:4) {
        for (i in 1:2) {
            for (k in 1:2) {
                values[n, i, k] <- (100 * k + 10 * i + n) / 1000
            }
        }
    }
    pvals <- list(pvals = values, nsim = 4L, m = 2L, K = 2L)

    out <- .gsd_kernel_matrix(pvals)
    expect_identical(dim(out), c(4L, 4L))
    # column (k - 1) * m + i is hypothesis i at analysis k
    expect_identical(out[, 1L], (110 + 1:4) / 1000)
    expect_identical(out[, 2L], (120 + 1:4) / 1000)
    expect_identical(out[, 3L], (210 + 1:4) / 1000)
    expect_identical(out[, 4L], (220 + 1:4) / 1000)

    with_na <- pvals
    with_na$pvals[2L, 1L, 2L] <- NA_real_
    expect_error(.gsd_kernel_matrix(with_na), "never rejectable")

    out_of_range <- pvals
    out_of_range$pvals[3L, 2L, 1L] <- 1.5
    expect_error(
        .gsd_kernel_matrix(out_of_range),
        "outside [0, 1]",
        fixed = TRUE
    )
})


test_that(".gsd_check_graph_dims() names both numbers", {
    expect_null(.gsd_check_graph_dims(3L, chk_pvals_gsd, "{.arg hyp_weight}"))

    expect_error(
        .gsd_check_graph_dims(2L, chk_pvals_gsd, "{.arg hyp_weight}"),
        "`hyp_weight` and `pvals` disagree on the number of hypotheses",
        fixed = TRUE
    )
    expect_error(
        .gsd_check_graph_dims(4L, chk_pvals_gsd, "{.arg graph_constraint}"),
        "describes m = 4 hypotheses, but the p-values hold m = 3",
        fixed = TRUE
    )
})


test_that(".gsd_stop_table() tabulates a rule on every rejection set", {
    # row mask + 1, hypothesis i at bit i - 1: the rule holds on {H1, H2},
    # {H1, H3} and {H1, H2, H3}, rows 4, 6 and 8
    one <- .gsd_stop_table(list(quote(r1 && (r2 || r3))), 3L, 2L)
    expect_identical(
        one,
        matrix(c(FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE))
    )

    two <- .gsd_stop_table(
        list(quote(r1), quote(r1 + r2 + r3 >= 2)),
        3L,
        3L
    )
    expect_identical(dim(two), c(8L, 2L))
    expect_identical(two[, 1L], rep(c(FALSE, TRUE), 4L))
    expect_identical(
        two[, 2L],
        c(FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE)
    )

    # K = 1: no interim, no column
    empty <- .gsd_stop_table(list(), 3L, 1L)
    expect_identical(dim(empty), c(8L, 0L))
    expect_type(empty, "logical")
})


test_that(".gsd_stop_table() refuses a bad condition or size", {
    expect_error(
        .gsd_stop_table(list(quote(r1 && NA)), 2L, 2L),
        "must return a single"
    )
    expect_error(
        .gsd_stop_table(list(quote(r1 && NA)), 2L, 2L),
        "returned NA"
    )
    expect_error(
        .gsd_stop_table(list(quote(c(r1, r2))), 2L, 2L),
        "must return a single"
    )
    expect_error(
        .gsd_stop_table(list(quote(r1)), 17L, 2L),
        "at most 16 hypotheses"
    )
    expect_error(
        .gsd_stop_table(list(quote(r1)), 2L, 3L),
        "one condition per interim analysis"
    )
})


test_that("a graph over more hypotheses than the p-values aborts in R", {
    # a one-hypothesis design: before the graph-dimension check, a
    # two-hypothesis graph reached the kernel and failed there
    pvals_1 <- withr::with_seed(71, {
        raw <- simulate_pvalues_gsd(0.9, info_frac = c(0.5, 1), nsim = 200)
        transform_pvalues_gsd(
            raw,
            spending = gsDesign::sfLDOF,
            grid_size = 64L
        )
    })
    gain_1 <- trial_success_gsd(d(t1), d = c(1, 0.75), verbose = "silent")

    expect_error(
        graph_optimise_gsd(
            pvals = pvals_1,
            graph_constraint = graph_constraint_free(2),
            trial_success = gain_1,
            global_search = FALSE,
            verbose = "silent"
        ),
        "describes m = 2 hypotheses, but the p-values hold m = 1",
        fixed = TRUE
    )
    expect_error(
        calc_power_pvals_gsd(
            pvals_1,
            hyp_weight = c(0.5, 0.5),
            trans_matrix = matrix(c(0, 1, 1, 0), nrow = 2)
        ),
        "describes m = 2 hypotheses, but the p-values hold m = 1",
        fixed = TRUE
    )
})
