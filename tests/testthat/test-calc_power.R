# Trial success function
custom_objective <- trial_success((r1 && r2) || r3, verbose = "silent")

# Test calc_power_pvals()
test_that("calc_power_pvals throws errors for incorrect graph structure", {
    pvals <- matrix(runif(100), nrow = 10)
    hyp_weight <- c(0.5, 0.3, 0.2)
    trans_matrix <- matrix(
        c(
            0, 0.5, 0.5,
            0.5, 0, 0.5,
            0.5, 0.5, 0
        ),
        ncol = 3
    )

    expect_error(
        calc_power_pvals(
            pvals,
            hyp_weight,
            trans_matrix,
            custom_power = list(
                custom_obj = custom_objective
            )
        ),
        "G must be an m x m matrix matching ncol(pvals)",
        fixed = TRUE
    )

    expect_error(
        calc_power_pvals(pvals),
        "`hyp_weight` must be a double, not absent."
    )
})

test_that("calc_power_pvals calculates power metrics", {
    pvals <- matrix(runif(30), nrow = 10)
    hyp_weight <- c(0.5, 0.3, 0.2)
    trans_matrix <- matrix(
        c(
            0, 0.5, 0.5,
            0.5, 0, 0.5,
            0.5, 0.5, 0
        ),
        ncol = 3
    )

    result <- calc_power_pvals(
        pvals,
        hyp_weight,
        trans_matrix,
        alpha = 0.05,
        custom_power = list(
            custom_obj = custom_objective
        )
    )

    expect_type(result, "list")
    expect_named(
        result,
        c(
            "local_power",
            "exp_rejections",
            "disj_power",
            "conj_power",
            "custom_obj"
        )
    )
    expect_length(result$local_power, 3)
    expect_type(result$exp_rejections, "double")
})


test_that("calc_power_pvals stops suboptimal graphs as default", {
    pvals <- matrix(runif(9000), nrow = 1000, ncol = 9)
    fs <- gMCPLite::fixedSequence(9)

    expect_snapshot(error = TRUE, {
        calc_power_pvals(
            pvals = pvals,
            hyp_weight = fs@weights,
            trans_matrix = fs@m,
            custom_power = list(
                custom_obj = custom_objective
            ),
            sum_to_one_constraint = TRUE
        )
    })
})

test_that("calc_power_pvals allows fixed sequence sum_to_one_constraint=TRUE", {
    pvals <- matrix(
        rbeta(90000, shape1 = 0.24, shape2 = 0.65),
        nrow = 10000,
        ncol = 9
    )
    fs <- gMCPLite::fixedSequence(9)

    expect_no_error(
        result <- calc_power_pvals(
            pvals = pvals,
            hyp_weight = fs@weights,
            trans_matrix = fs@m,
            custom_power = list(
                custom_obj = custom_objective
            ),
            sum_to_one_constraint = FALSE
        )
    )

    expect_type(result, "list")
    expect_named(
        result,
        c(
            "local_power",
            "exp_rejections",
            "disj_power",
            "conj_power",
            "custom_obj"
        )
    )
})


test_that("calc_power_pvals allows list of anonymous functions", {
    pvals <- matrix(
        rbeta(50000, shape1 = 0.24, shape2 = 0.65),
        nrow = 10000,
        ncol = 5
    )
    bh <- gMCPLite::BonferroniHolm(5)

    expect_no_error(
        result <- calc_power_pvals(
            pvals = pvals,
            hyp_weight = bh@weights,
            trans_matrix = bh@m,
            custom_power = list(
                custom_obj1 = function(x) {
                    0.25 * (x[1] + x[2] + x[3] + as.double(x[4] || x[5]))
                },
                custom_obj2 = function(x) x[1] && x[2] && x[3] && x[4] && x[5]
            ),
            sum_to_one_constraint = FALSE
        )
    )

    expect_type(result, "list")
    expect_named(
        result,
        c(
            "local_power",
            "exp_rejections",
            "disj_power",
            "conj_power",
            "custom_obj1",
            "custom_obj2"
        )
    )

    # Check custom_obj2 is the same as conj_power
    expect_equal(
        result$conj_power,
        result$custom_obj2,
        tolerance = sqrt(.Machine$double.eps)
    )
})


test_that("calc_power_pvals with list of anonymous funcs and trial_success", {
    pvals <- matrix(
        rbeta(30000, shape1 = 0.24, shape2 = 0.65),
        nrow = 10000,
        ncol = 3
    )
    bh <- gMCPLite::BonferroniHolm(3)

    expect_no_error(
        result <- calc_power_pvals(
            pvals = pvals,
            hyp_weight = bh@weights,
            trans_matrix = bh@m,
            custom_power = list(
                anon_obj = function(x) {
                    (x[1] && x[2]) || x[3]
                },
                custom_objective = custom_objective
            ),
            sum_to_one_constraint = FALSE
        )
    )

    expect_type(result, "list")
    expect_named(
        result,
        c(
            "local_power",
            "exp_rejections",
            "disj_power",
            "conj_power",
            "anon_obj",
            "custom_objective"
        )
    )

    # Check anon_obj is the same as custom_objective
    expect_equal(
        result$anon_obj,
        result$custom_objective,
        tolerance = sqrt(.Machine$double.eps)
    )
})


test_that("calc_power_pvals allows single anonymous function", {
    pvals <- matrix(
        rbeta(40000, shape1 = 0.24, shape2 = 0.65),
        nrow = 10000,
        ncol = 4
    )
    bh <- gMCPLite::BonferroniHolm(4)

    expect_no_error(
        result <- calc_power_pvals(
            pvals = pvals,
            hyp_weight = bh@weights,
            trans_matrix = bh@m,
            custom_power = function(x) {
                0.5 * (as.double((x[1] && x[2]) || x[3]) + x[4])
            },
            sum_to_one_constraint = FALSE
        )
    )

    expect_type(result, "list")
    expect_named(
        result,
        c(
            "local_power",
            "exp_rejections",
            "disj_power",
            "conj_power",
            "custom_power"
        )
    )
})


test_that("calc_power_pvals allows single trial_success function", {
    pvals <- matrix(
        rbeta(40000, shape1 = 0.24, shape2 = 0.65),
        nrow = 10000,
        ncol = 4
    )
    bh <- gMCPLite::BonferroniHolm(4)

    expect_no_error(
        result <- calc_power_pvals(
            pvals = pvals,
            hyp_weight = bh@weights,
            trans_matrix = bh@m,
            custom_power = custom_objective,
            sum_to_one_constraint = FALSE
        )
    )

    expect_type(result, "list")
    expect_named(
        result,
        c(
            "local_power",
            "exp_rejections",
            "disj_power",
            "conj_power",
            "custom_power"
        )
    )
})

test_that("calc_power_pvals validates input types", {
    pvals <- matrix(runif(30), nrow = 10)
    hyp_weight <- c(1 / 3, 1 / 3, 1 / 3)
    trans_matrix <- matrix(0, nrow = 3, ncol = 3)

    expect_error(
        calc_power_pvals("not_a_matrix", hyp_weight, trans_matrix),
        '`pvals` must be a double matrix, not the string "not_a_matrix"'
    )

    expect_error(
        calc_power_pvals(pvals, hyp_weight = "foo"),
        '`hyp_weight` must be a double, not the string "foo"'
    )

    expect_error(
        calc_power_pvals(pvals, hyp_weight, trans_matrix = "bar"),
        '`trans_matrix` must be a double matrix, not the string "bar"'
    )

    expect_error(
        calc_power_pvals(pvals, hyp_weight, trans_matrix, alpha = "bad"),
        '`alpha` must be a number, not the string "bad"'
    )

    expect_error(
        calc_power_pvals(
            pvals,
            hyp_weight,
            trans_matrix,
            sum_to_one_constraint = "yes"
        ),
        '`sum_to_one_constraint` must be a logical vector, not the string "yes"'
    )
})


test_that("calc_power_pvals complains", {
    # when custom_power is not a list or trial_success
    pvals <- matrix(
        rbeta(40000, shape1 = 0.24, shape2 = 0.65),
        nrow = 10000,
        ncol = 4
    )
    bh <- gMCPLite::BonferroniHolm(4)

    expect_error(
        calc_power_pvals(
            pvals = pvals,
            hyp_weight = bh@weights,
            trans_matrix = bh@m,
            custom_power = "foo",
            sum_to_one_constraint = FALSE
        ),
        "must be a function"
    )
})

test_that("calc_power_pvals complains when users passes anything via dots", {
    pvals <- matrix(
        rbeta(40000, shape1 = 0.24, shape2 = 0.65),
        nrow = 10000,
        ncol = 4
    )
    alpha <- 0.024
    bh <- gMCPLite::BonferroniHolm(4)

    expect_snapshot(error = TRUE, {
        calc_power_pvals(
            pvals = pvals,
            hyp_weight = bh@weights,
            trans_matrix = bh@m,
            alpha
        )
    })
})

# --- .auto_name_custom_power() unit tests ---

test_that(".auto_name_custom_power returns empty list for NULL input", {
    out <- .auto_name_custom_power(NULL)
    expect_type(out, "list")
    expect_length(out, 0L)
})

test_that(".auto_name_custom_power wraps a single function in a named list", {
    fn <- function(x) x[1]
    out <- .auto_name_custom_power(fn)
    expect_type(out, "list")
    expect_named(out, "custom_power")
    expect_identical(out$custom_power, fn)
})

test_that(".auto_name_custom_power wraps a trial_success in a named list", {
    out <- .auto_name_custom_power(custom_objective)
    expect_type(out, "list")
    expect_named(out, "custom_power")
    expect_true(is_trial_success(out$custom_power))
})

test_that(".auto_name_custom_power errors on non-function, non-list input", {
    expect_error(.auto_name_custom_power("bad"), "must be a function")
    expect_error(.auto_name_custom_power(42L), "must be a function")
})

test_that(".auto_name_custom_power errors on invalid list elements", {
    fn <- function(x) x[1]
    expect_error(
        .auto_name_custom_power(list(fn, "bad")),
        "must be a function"
    )
    expect_error(
        .auto_name_custom_power(list(a = fn, b = 42L)),
        "must be a function"
    )
})

test_that(".auto_name_custom_power assigns funcN names to fully unnamed list", {
    fn1 <- function(x) x[1]
    fn2 <- function(x) x[2]
    out <- .auto_name_custom_power(list(fn1, fn2))
    expect_named(out, c("func1", "func2"))
})

test_that(".auto_name_custom_power fills only blank positions", {
    # while preserving existing names
    fn1 <- function(x) x[1]
    fn2 <- function(x) x[2]
    fn3 <- function(x) x[3]
    out <- .auto_name_custom_power(list(alpha = fn1, fn2, gamma = fn3))
    expect_named(out, c("alpha", "func2", "gamma"))
})

test_that(".auto_name_custom_power does not change a fully named list", {
    fn1 <- function(x) x[1]
    fn2 <- function(x) x[2]
    inp <- list(a = fn1, b = fn2)
    out <- .auto_name_custom_power(inp)
    expect_identical(out, inp)
})


# --- .eval_custom_power() unit tests ---

test_that(".eval_custom_power with anonymous functions", {
    # evaluates anonymous functions row-wise and returns means
    rej <- matrix(
        c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE),
        nrow = 3, ncol = 2
    )
    fns <- list(h1 = function(x) as.numeric(x[1]))
    out <- .eval_custom_power(fns, rej)
    expect_named(out, "h1")
    expect_equal(out$h1, mean(rej[, 1]), tolerance = sqrt(.Machine$double.eps))
})

test_that(".eval_custom_power with trial_success objects", {
    # evaluates trial_success objects against full matrix
    rej <- matrix(
        c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE),
        nrow = 3, ncol = 3
    )
    out <- .eval_custom_power(list(ts = custom_objective), rej)
    expect_named(out, "ts")
    expect_type(out$ts, "double")
    expect_gte(out$ts, 0)
    expect_lte(out$ts, 1)
})

test_that(".eval_custom_power complains", {
    # errors when a list element is neither function nor trial_success
    rej <- matrix(c(TRUE, FALSE, TRUE, TRUE), nrow = 2, ncol = 2)
    expect_error(
        .eval_custom_power(list(bad = "not a function"), rej),
        "must be a function"
    )
})
