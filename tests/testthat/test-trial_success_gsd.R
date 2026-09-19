# Mean with plain double accumulation in row order, as the compiled loop does
# (R's mean() accumulates in long double and can differ in the last bit)
mean_double <- function(x) {
    Reduce(`+`, x) / length(x)
}

# Decision-time matrix covering every (t1, t2) pair in 0..2 (0 = never)
gsd_time_grid <- function() {
    tm <- as.matrix(expand.grid(t1 = 0:2, t2 = 0:2))
    dimnames(tm) <- NULL
    storage.mode(tm) <- "integer"
    tm
}

# --- the four gain functions of design record section 4.6 -------------------

test_that("Example 5 gain compiles and matches R on a hand-built time matrix", {
    tm <- gsd_time_grid()
    v_pfs <- 0.4
    v_os <- 1
    d <- c(1, 0.75)
    dd <- c(0, d)

    ts <- trial_success_gsd(
        !!v_pfs * d(t1) + !!v_os * d(t2),
        d = d,
        verbose = "silent"
    )

    expect_s3_class(ts, "multigrain_trial_success_gsd")
    expect_s3_class(ts, "multigrain_trial_success")
    expect_identical(ts$m, 2L)
    expect_identical(ts$K, 2L)
    expect_identical(ts$tables, list(d = d))
    expect_identical(ts$objective, "0.4 * d(t1) + 1 * d(t2)")

    expected <- mean_double(v_pfs * dd[tm[, 1] + 1] + v_os * dd[tm[, 2] + 1])
    expect_identical(ts$func(tm), expected)
})

test_that("supplement dual-primary gain with same-look bonus matches R", {
    tm <- gsd_time_grid()
    nu1 <- 0.5
    nu2 <- 0.3
    nu12 <- 0.2
    delta <- 0.8
    delta12 <- 0.6
    dd <- c(0, 1, delta)

    ts <- trial_success_gsd(
        !!nu1 * d(t1) + !!nu2 * d(t2) +
            !!nu12 * ((t1 == 1 && t2 == 1) + !!delta12 * (t1 == 2 && t2 == 2)),
        d = c(1, delta),
        verbose = "silent"
    )

    t1 <- tm[, 1]
    t2 <- tm[, 2]
    expected <- mean_double(
        nu1 * dd[t1 + 1] + nu2 * dd[t2 + 1] +
            nu12 * ((t1 == 1 & t2 == 1) + delta12 * (t1 == 2 & t2 == 2))
    )
    expect_identical(ts$func(tm), expected)
})

test_that("supplement co-primary gain at the same look matches R", {
    tm <- gsd_time_grid()
    nu12 <- 0.7
    delta12 <- 0.5

    ts <- trial_success_gsd(
        !!nu12 * ((t1 == 1 && t2 == 1) + !!delta12 * (t1 == 2 && t2 == 2)),
        verbose = "silent"
    )
    expect_null(ts$K)

    t1 <- tm[, 1]
    t2 <- tm[, 2]
    expected <- mean_double(
        nu12 * ((t1 == 1 & t2 == 1) + delta12 * (t1 == 2 & t2 == 2))
    )
    expect_identical(ts$func(tm), expected)
})

test_that("supplement PFS and OS gain matches R", {
    tm <- gsd_time_grid()
    a1 <- 1
    b1 <- 0.5
    a2 <- 0.8
    b2 <- 0.4
    c2 <- 0.3
    b3 <- 0.2

    ts <- trial_success_gsd(
        (t2 == 1) * (!!a1 + !!b1 * (t1 == 1)) +
            (t2 == 2) * (!!a2 + !!b2 * (t1 == 1) + !!c2 * (t1 == 2)) +
            (t2 == 0) * (t1 == 1) * !!b3,
        K = 2,
        verbose = "silent"
    )
    expect_identical(ts$K, 2L)

    t1 <- tm[, 1]
    t2 <- tm[, 2]
    expected <- mean_double(
        (t2 == 1) * (a1 + b1 * (t1 == 1)) +
            (t2 == 2) * (a2 + b2 * (t1 == 1) + c2 * (t1 == 2)) +
            (t2 == 0) * (t1 == 1) * b3
    )
    expect_identical(ts$func(tm), expected)
})

# --- discount tables --------------------------------------------------------

test_that("d(0) is 0 and table values round-trip exactly", {
    # an increasing table on purpose, to check the round trip of the values;
    # `.gsd_gain_tables()` warns about the ordering, which is not the point
    # of this test
    ts <- suppressWarnings(
        trial_success_gsd(d(t1), d = c(1 / 3, 2 / 3), verbose = "silent")
    )

    never <- matrix(0L, nrow = 3L, ncol = 1L)
    expect_identical(ts$func(never), 0)

    expect_identical(ts$func(matrix(1L, 2L, 1L)), 1 / 3)
    expect_identical(ts$func(matrix(2L, 2L, 1L)), 2 / 3)

    mixed <- matrix(c(0L, 1L, 2L, 0L), ncol = 1L)
    expect_identical(ts$func(mixed), (0 + 1 / 3 + 2 / 3 + 0) / 4)
})

test_that("two tables of the same length are both usable", {
    # `e` is outside [0, 1] on purpose: a value table is allowed, it only
    # warns
    ts <- suppressWarnings(
        trial_success_gsd(
            d(t1) + e(t2),
            d = c(1, 0.5),
            e = c(2, 1),
            verbose = "silent"
        )
    )
    expect_identical(ts$K, 2L)
    expect_named(ts$tables, c("d", "e"))

    tm <- gsd_time_grid()
    expected <- mean(c(0, 1, 0.5)[tm[, 1] + 1] + c(0, 2, 1)[tm[, 2] + 1])
    expect_identical(ts$func(tm), expected)
})

# --- K handling -------------------------------------------------------------

test_that("K is inferred from tables, checked against an explicit K", {
    tables_only <- trial_success_gsd(d(t1), d = c(1, 1, 1), verbose = "silent")
    expect_identical(tables_only$K, 3L)

    agree <- trial_success_gsd(d(t1), d = c(1, 1, 1), K = 3, verbose = "silent")
    expect_identical(agree$K, 3L)

    expect_error(
        trial_success_gsd(d(t1), d = c(1, 0.75), K = 3, verbose = "silent"),
        "`K` is 3 but the discount table `d` has length 2"
    )
    expect_error(
        trial_success_gsd(
            d(t1) + e(t1),
            d = c(1, 0.75),
            e = c(1, 1, 1),
            verbose = "silent"
        ),
        "All discount tables must have the same length"
    )
})

test_that("K without tables is stored or left NULL", {
    expect_null(trial_success_gsd(r1 + r2, verbose = "silent")$K)
    expect_identical(trial_success_gsd(r1, K = 3L, verbose = "silent")$K, 3L)
    expect_identical(trial_success_gsd(r1, K = 3, verbose = "silent")$K, 3L)
    expect_error(
        trial_success_gsd(r1, K = 0, verbose = "silent"),
        "`K` must be a whole number larger than or equal to 1"
    )
    expect_error(
        trial_success_gsd(r1, K = 1.5, verbose = "silent"),
        "`K` must be a whole number"
    )
})

# --- validation errors ------------------------------------------------------

test_that("un-injected symbols error with the !! hint", {
    w <- 2
    expect_error(
        trial_success_gsd(w * r1, verbose = "silent"),
        "not a rejection indicator"
    )
    expect_error(
        trial_success_gsd(w * r1, verbose = "silent"),
        "!!w",
        fixed = TRUE
    )
    expect_error(
        trial_success_gsd(r1 * (w + t2), verbose = "silent"),
        "not a rejection indicator"
    )
})

test_that("unsupported functions and operators error", {
    expect_error(
        trial_success_gsd(sqrt(t1), verbose = "silent"),
        "Unsupported operator or function `sqrt`"
    )
    expect_error(
        trial_success_gsd(min(t1, t2), verbose = "silent"),
        "Unsupported operator or function `min`"
    )
    expect_error(
        trial_success_gsd(!r1, verbose = "silent"),
        "Unsupported operator or function `!`"
    )
    expect_error(
        trial_success_gsd("sqrt(t1)", verbose = "silent"),
        "Unsupported operator or function `sqrt`"
    )
})

test_that("discount table misuse errors", {
    expect_error(
        trial_success_gsd(r1, c(1, 0.75), verbose = "silent"),
        "must be named"
    )
    expect_error(
        trial_success_gsd(d(t1), d = "a", verbose = "silent"),
        "must be a non-empty numeric vector"
    )
    expect_error(
        trial_success_gsd(d(t1), d = c(1, NA), verbose = "silent"),
        "must be a non-empty numeric vector"
    )
    expect_error(
        trial_success_gsd(d(t1), d = numeric(), verbose = "silent"),
        "must be a non-empty numeric vector"
    )
    expect_error(
        trial_success_gsd(t1(t1), t1 = c(1, 1), verbose = "silent"),
        "not allowed"
    )
    expect_error(
        trial_success_gsd(and(t1), and = c(1, 1), verbose = "silent"),
        "not allowed"
    )
    expect_error(
        trial_success_gsd(d.1(t1), d.1 = c(1, 1), verbose = "silent"),
        "not allowed"
    )
    expect_error(
        trial_success_gsd(d(r1), d = c(1, 0.75), verbose = "silent"),
        "applied to a single decision time symbol"
    )
    expect_error(
        trial_success_gsd(d(1), d = c(1, 0.75), verbose = "silent"),
        "applied to a single decision time symbol"
    )
    expect_error(
        trial_success_gsd(d(t1 + 1), d = c(1, 0.75), verbose = "silent"),
        "applied to a single decision time symbol"
    )
    expect_error(
        trial_success_gsd(d(t1, t2), d = c(1, 0.75), verbose = "silent"),
        "applied to a single decision time symbol"
    )
    expect_error(
        trial_success_gsd(d * r1, d = c(1, 0.75), verbose = "silent"),
        "cannot stand alone"
    )
    expect_error(
        trial_success_gsd(d(t1), verbose = "silent"),
        "Unsupported operator or function `d`"
    )
})

test_that("injected vectors, empty expressions and bad types error", {
    d <- c(1, 0.75)
    expect_error(
        trial_success_gsd(!!d * r1, verbose = "silent"),
        "must be a single number"
    )
    expect_error(
        trial_success_gsd("4", verbose = "silent"),
        "must reference at least one rejection indicator"
    )
    expect_error(
        trial_success_gsd(42, verbose = "silent"),
        "must be an expression or a character string"
    )
    two_strings <- c("r1", "r2")
    expect_error(
        trial_success_gsd(!!two_strings, verbose = "silent"),
        "must be a single string"
    )
    expect_error(
        resolve_expr_gsd(list(a = 1)),
        "must be an expression or a character string"
    )
    expect_error(
        resolve_expr_gsd(NULL),
        "must be an expression or a character string"
    )
})

test_that("the bool rule is unchanged: && and || need boolean operands", {
    expect_error(
        trial_success_gsd(t1 && r2, verbose = "silent"),
        "`&&` (AND) only allowed between booleans",
        fixed = TRUE
    )
    expect_error(
        trial_success_gsd(5 || r1, verbose = "silent"),
        "`||` (OR) only allowed between booleans",
        fixed = TRUE
    )
    expect_error(
        trial_success_gsd(d(t1) || r1, d = c(1, 1), verbose = "silent"),
        "`||` (OR) only allowed between booleans",
        fixed = TRUE
    )
    # R precedence: (r1 + r2) && r3 is real && bool
    expect_error(
        trial_success_gsd(r1 + r2 && r3, verbose = "silent"),
        "`&&` (AND) only allowed between booleans",
        fixed = TRUE
    )
})

# --- precedence and the r/t relationship ------------------------------------

test_that("comparisons combine with && at R precedence (record 8.10)", {
    tm <- gsd_time_grid()
    t1 <- tm[, 1]
    t2 <- tm[, 2]

    both_first <- trial_success_gsd(t1 == 1 && t2 == 1, verbose = "silent")
    expect_identical(both_first$func(tm), mean(t1 == 1 & t2 == 1))

    with_r <- trial_success_gsd((t1 == 1) && r2, verbose = "silent")
    expect_identical(with_r$func(tm), mean(t1 == 1 & t2 > 0))

    # word operators in a string
    words <- trial_success_gsd("t1 == 1 and t2 == 1", verbose = "silent")
    expect_identical(words$func(tm), mean(t1 == 1 & t2 == 1))
    either <- trial_success_gsd("t1 == 1 OR t2 == 2", verbose = "silent")
    expect_identical(either$func(tm), mean(t1 == 1 | t2 == 2))

    # every comparison operator, against a literal or an expression
    ops <- trial_success_gsd(
        (t1 != 0) + (t1 < 2) + (t1 <= 1) + (t2 > t1) + (t2 >= t1 + 1),
        verbose = "silent"
    )
    expect_identical(
        ops$func(tm),
        mean((t1 != 0) + (t1 < 2) + (t1 <= 1) + (t2 > t1) + (t2 >= t1 + 1))
    )
})

test_that("r<i> is t<i> > 0, and m spans both symbol kinds", {
    tm <- gsd_time_grid()
    by_r <- trial_success_gsd(r1 + r2, verbose = "silent")
    by_t <- trial_success_gsd((t1 > 0) + (t2 > 0), verbose = "silent")
    expect_identical(by_r$func(tm), by_t$func(tm))
    expect_identical(by_r$func(tm), mean(rowSums(tm > 0)))

    expect_identical(trial_success_gsd(t2 * r1, verbose = "silent")$m, 2L)
    expect_warning(
        ts <- trial_success_gsd(t3 + r1, verbose = "silent"),
        "missing 2"
    )
    expect_identical(ts$m, 3L)
})

test_that("|| and && follow R precedence (review finding)", {
    tm <- as.matrix(expand.grid(t1 = 0:1, t2 = 0:1, t3 = 0:1))
    dimnames(tm) <- NULL
    storage.mode(tm) <- "integer"
    r1 <- tm[, 1] > 0
    r2 <- tm[, 2] > 0
    r3 <- tm[, 3] > 0

    ts <- trial_success_gsd(r1 || r2 && r3, verbose = "silent")
    expect_identical(ts$func(tm), mean(r1 | (r2 & r3)))
    expect_false(ts$func(tm) == mean((r1 | r2) & r3))

    words <- trial_success_gsd("r1 or r2 and r3", verbose = "silent")
    expect_identical(words$func(tm), mean(r1 | (r2 & r3)))
})

test_that("the compiled function refuses a matrix with too few columns", {
    ts <- trial_success_gsd(r1 + r2, verbose = "silent")
    expect_error(
        ts$func(matrix(1L, 3L, 1L)),
        "has 1 column(s) but the trial success function refers to 2",
        fixed = TRUE
    )
    # extra columns are ignored
    expect_identical(ts$func(matrix(1L, 3L, 5L)), 2)
})

test_that("unary minus compiles and evaluates", {
    tm <- gsd_time_grid()
    ts <- trial_success_gsd(-r1 + 1, verbose = "silent")
    expect_identical(ts$func(tm), mean(1 - (tm[, 1] > 0)))
    ts2 <- trial_success_gsd(-(t1 - t2), verbose = "silent")
    expect_identical(ts2$func(tm), mean(-(tm[, 1] - tm[, 2])))
})

test_that("numeric literals of every spelling become valid C++", {
    tm <- gsd_time_grid()
    tiny <- 1e-05
    big <- 1e+06
    ts <- trial_success_gsd(
        !!tiny * r1 + !!big * r2 + 3 * t1,
        verbose = "silent"
    )
    expect_identical(
        ts$func(tm),
        mean_double(tiny * (tm[, 1] > 0) + big * (tm[, 2] > 0) + 3 * tm[, 1])
    )
    expect_identical(
        replace_indices_gsd("r1 + 5"),
        "double(t(i, 0) > 0) + 5.0"
    )

    # injected constants round-trip exactly, like table values
    third <- 1 / 3
    exact <- trial_success_gsd(!!third * r1, verbose = "silent")
    expect_identical(exact$func(matrix(1L, 4L, 1L)), third)
    expect_identical(.gsd_gain_cpp_number(0.75), "0.75")
    expect_identical(.gsd_gain_cpp_number(1 / 3), "0.3333333333333333")
    expect_identical(.gsd_gain_cpp_number(2), "2.0")
    expect_identical(.gsd_gain_cpp_number(1e-05), "1e-05")
    expect_identical(.gsd_gain_cpp_number(1e+06), "1000000.0")
    expect_identical(.gsd_gain_cpp_number(0.1 + 0.2), "0.30000000000000004")
    # the display string keeps the fixed-sample deparse, the code is exact
    expect_identical(exact$objective, "0.333333333333333 * r1")
    expect_match(exact$cpp_code, "0.3333333333333333 * double", fixed = TRUE)
    expect_identical(
        replace_indices_gsd("TRUE * r1"),
        "1.0 * double(t(i, 0) > 0)"
    )
})

# --- generated code ---------------------------------------------------------

test_that("replace_indices_gsd emits the record's C++ forms", {
    expect_identical(replace_indices_gsd("r1"), "double(t(i, 0) > 0)")
    expect_identical(replace_indices_gsd("t1"), "double(t(i, 0))")
    expect_identical(replace_indices_gsd("t10"), "double(t(i, 9))")
    expect_identical(
        replace_indices_gsd("d(t2)", table_names = "d"),
        "d_tab[t(i, 1)]"
    )
    expect_identical(
        replace_indices_gsd("t1 == 1"),
        "double(double(t(i, 0)) == 1.0)"
    )
    expect_identical(
        replace_indices_gsd("t1 == 1 && t2 == 1"),
        "double(double(t(i, 0)) == 1.0) * double(double(t(i, 1)) == 1.0)"
    )
    expect_identical(
        replace_indices_gsd("r1 || r2"),
        "std_min(double(1), double(t(i, 0) > 0) + double(t(i, 1) > 0))"
    )
    expect_identical(
        replace_indices_gsd("1 / 3 * (r1 + r2)"),
        "1.0 / 3.0 * (double(t(i, 0) > 0) + double(t(i, 1) > 0))"
    )
})

test_that("cpp_code snapshots", {
    ts <- trial_success_gsd(
        0.4 * d(t1) + 1 * d(t2),
        d = c(1, 0.75),
        verbose = "silent"
    )
    expect_snapshot(cat(ts$cpp_code))

    ts2 <- trial_success_gsd(
        (t1 == 1 && t2 == 1) + 0.5 * (t1 == 2 && t2 == 2),
        verbose = "silent"
    )
    expect_snapshot(cat(ts2$cpp_code))
})

test_that("compiling twice gives identical results", {
    tm <- gsd_time_grid()
    dt <- c(1, 0.75)
    a <- trial_success_gsd(0.4 * d(t1) + d(t2), d = dt, verbose = "silent")
    b <- trial_success_gsd(0.4 * d(t1) + d(t2), d = dt, verbose = "silent")
    expect_identical(a$cpp_code, b$cpp_code)
    expect_identical(a$func(tm), b$func(tm))
    expect_identical(a$func(tm), a$func(tm))
})

# --- class, checks, methods -------------------------------------------------

test_that("class helpers accept the GSD object", {
    ts <- trial_success_gsd(r1 + r2, verbose = "silent")
    expect_true(is_trial_success_gsd(ts))
    expect_true(is_trial_success(ts))
    expect_false(is_trial_success_gsd(trial_success(r1, verbose = "silent")))
    expect_false(is_trial_success_gsd("foo"))
    expect_no_error(check_trial_success(ts))
})

test_that("print and summary methods", {
    ts <- trial_success_gsd(
        0.4 * d(t1) + d(t2),
        d = c(1, 0.75),
        verbose = "silent"
    )
    expect_snapshot(print(ts))
    expect_snapshot(summary(ts))
    expect_null(print.multigrain_trial_success_gsd(NULL))
    expect_null(summary.multigrain_trial_success_gsd(NULL))

    no_k <- trial_success_gsd(r1, verbose = "silent")
    expect_snapshot(print(no_k))
})

test_that("verbose handling matches trial_success()", {
    expect_snapshot(trial_success_gsd(r1 + r2, verbose = "info"))
    expect_snapshot(trial_success_gsd(r1 + r2, verbose = "silent"))
    expect_snapshot(trial_success_gsd(r1 + r2, verbose = TRUE))
    expect_snapshot(trial_success_gsd(r1 + r2, verbose = FALSE))
    expect_snapshot(error = TRUE, {
        trial_success_gsd(r1 + r2, verbose = 2)
    })
})

# --- discount table warnings (design record, section 6 P3 [Rev 2026-09-18]) --

test_that("an increasing discount table warns naming the analyses", {
    expect_warning(
        trial_success_gsd(d(t1), d = c(0.5, 1), verbose = "silent"),
        "increases from analysis 1 (0.5) to analysis 2 (1)",
        fixed = TRUE
    )
})

test_that("a discount table outside [0, 1] warns naming both values", {
    w <- tryCatch(
        trial_success_gsd(d(t1), d = c(1.2, -0.1), verbose = "silent"),
        warning = function(cnd) conditionMessage(cnd)
    )
    expect_true(grepl("outside [0, 1]", w, fixed = TRUE))
    expect_true(grepl("1.2", w, fixed = TRUE))
    expect_true(grepl("-0.1", w, fixed = TRUE))
})

test_that("a leading zero gets the shifted-by-one hint", {
    w <- tryCatch(
        trial_success_gsd(d(t1), d = c(0, 1, 0.75), verbose = "silent"),
        warning = function(cnd) conditionMessage(cnd)
    )
    expect_true(grepl("increases from analysis 1 (0) to analysis 2 (1)", w,
        fixed = TRUE
    ))
    expect_true(grepl("shifts every value by one analysis", w, fixed = TRUE))
})

test_that("a well-behaved discount table does not warn", {
    expect_no_warning(
        trial_success_gsd(d(t1), d = c(1, 0.75), verbose = "silent")
    )
    expect_no_warning(
        trial_success_gsd(d(t1), d = c(1, 1), verbose = "silent")
    )
})

test_that("both checks tolerate rounding", {
    # a table computed rather than typed must not be flagged
    expect_no_warning(
        trial_success_gsd(d(t1), d = c(1, 1 + 1e-12), verbose = "silent")
    )
    expect_no_warning(
        trial_success_gsd(d(t1), d = c(0.5, -1e-12), verbose = "silent")
    )

    # `c(-1e-12, 0.5)` is inside the range to tolerance, but it does increase,
    # so it warns about the ordering only
    w <- tryCatch(
        trial_success_gsd(d(t1), d = c(-1e-12, 0.5), verbose = "silent"),
        warning = function(cnd) conditionMessage(cnd)
    )
    expect_false(grepl("outside [0, 1]", w, fixed = TRUE))
    expect_true(grepl("increases from analysis 1", w, fixed = TRUE))

    expect_warning(
        trial_success_gsd(d(t1), d = c(1, 1 + 1e-6), verbose = "silent"),
        "outside [0, 1]",
        fixed = TRUE
    )
})

test_that("one warning carries both messages", {
    w <- tryCatch(
        trial_success_gsd(d(t1), d = c(-0.5, 1.5), verbose = "silent"),
        warning = function(cnd) conditionMessage(cnd)
    )
    expect_true(grepl("outside [0, 1]", w, fixed = TRUE))
    expect_true(grepl("increases from analysis 1", w, fixed = TRUE))
})

test_that("only the offending table is named", {
    w <- tryCatch(
        trial_success_gsd(
            d(t1) + e(t2),
            d = c(1, 0.75),
            e = c(0.5, 1),
            verbose = "silent"
        ),
        warning = function(cnd) conditionMessage(cnd)
    )
    expect_true(grepl("`e`", w, fixed = TRUE))
    expect_false(grepl("`d`", w, fixed = TRUE))
})

test_that("the compiled function still works after a warning", {
    ts <- suppressWarnings(
        trial_success_gsd(d(t1), d = c(0.5, 1), verbose = "silent")
    )
    tm <- gsd_time_grid()
    dd <- c(0, 0.5, 1)
    expect_identical(
        ts$func(tm),
        mean_double(lapply(seq_len(nrow(tm)), \(i) dd[[tm[i, 1L] + 1L]]))
    )
})
