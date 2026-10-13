# A hand-built `multigrain_pvals_gsd` carrying only the fields that
# `.gsd_check_stopping_rule_dims()` reads
gsd_rule_pvals <- function(m = 3L, K = 3L) { # nolint: object_name_linter.
    structure(list(m = m, K = K), class = "multigrain_pvals_gsd")
}

# --- the grammar accepts ----------------------------------------------------

test_that("expressions in the grammar are accepted", {
    one <- stopping_rule_gsd(r1)
    expect_s3_class(one, "multigrain_stopping_rule_gsd")
    expect_identical(one$conditions, list(quote(r1)))
    expect_identical(one$labels, "r1")
    expect_identical(one$m_min, 1L)
    expect_null(one$n_interim)

    nested <- stopping_rule_gsd(r1 && (r2 || r3))
    expect_identical(nested$conditions, list(quote(r1 && (r2 || r3))))
    expect_identical(nested$labels, "r1 && (r2 || r3)")
    expect_identical(nested$m_min, 3L)
    expect_null(nested$n_interim)

    negated <- stopping_rule_gsd(!r1)
    expect_identical(negated$labels, "!r1")
    expect_identical(negated$m_min, 1L)

    counting <- stopping_rule_gsd(r1 + r2 + r3 >= 2)
    expect_identical(counting$labels, "r1 + r2 + r3 >= 2")
    expect_identical(counting$m_min, 3L)

    # the grammar admits a logical literal; the rule needs an r<i> as well
    expect_no_error(.gsd_rule_validate(TRUE))
    with_true <- stopping_rule_gsd(r2 && TRUE)
    expect_identical(with_true$labels, "r2 && TRUE")
    expect_identical(with_true$m_min, 2L)

    k <- 1
    injected <- stopping_rule_gsd(!!k <= r1 + r2)
    expect_identical(injected$labels, "1 <= r1 + r2")
    expect_identical(injected$m_min, 2L)
    expect_null(injected$n_interim)
})

test_that("strings are parsed with the word operators", {
    words <- stopping_rule_gsd("r1 and (r2 or r3)")
    expect_identical(words$conditions, list(quote(r1 && (r2 || r3))))
    expect_identical(words$labels, "r1 && (r2 || r3)")
    expect_identical(words$m_min, 3L)
    expect_null(words$n_interim)

    negated <- stopping_rule_gsd("not r1")
    expect_identical(negated$conditions, list(quote(!r1)))
    expect_identical(negated$labels, "!r1")
    expect_identical(negated$m_min, 1L)

    upper <- stopping_rule_gsd("r1 AND NOT r2 Or r3")
    expect_identical(upper$conditions, list(quote(r1 && !r2 || r3)))
})

test_that("a list gives one condition per interim analysis", {
    mixed <- stopping_rule_gsd(list(r1, "r2"))
    expect_identical(mixed$conditions, list(quote(r1), quote(r2)))
    expect_identical(mixed$labels, c("r1", "r2"))
    expect_identical(mixed$m_min, 2L)
    expect_identical(mixed$n_interim, 2L)

    two <- stopping_rule_gsd(list(r1, r1 && (r2 || r3)))
    expect_identical(two$labels, c("r1", "r1 && (r2 || r3)"))
    expect_identical(two$m_min, 3L)
    expect_identical(two$n_interim, 2L)

    # a list of one is one interim, not a recycled condition
    single <- stopping_rule_gsd(list(r4))
    expect_identical(single$m_min, 4L)
    expect_identical(single$n_interim, 1L)
})

test_that("a character vector gives one condition per interim analysis", {
    conds <- c("r1", "r1 and r2")
    rule <- stopping_rule_gsd(!!conds)
    expect_identical(rule$conditions, list(quote(r1), quote(r1 && r2)))
    expect_identical(rule$labels, c("r1", "r1 && r2"))
    expect_identical(rule$m_min, 2L)
    expect_identical(rule$n_interim, 2L)
})

test_that("a literal c() of strings is read like list()", {
    rule <- stopping_rule_gsd(c("r1", "r2 && r3"))
    expect_identical(rule$conditions, list(quote(r1), quote(r2 && r3)))
    expect_identical(rule$labels, c("r1", "r2 && r3"))
    expect_identical(rule$m_min, 3L)
    expect_identical(rule$n_interim, 2L)

    # nolint next: unnecessary_concatenation_linter. the c() call is the input
    expect_identical(stopping_rule_gsd(c("r2"))$n_interim, 1L)

    # anything but single strings in c() is an unsupported call
    expect_error(
        stopping_rule_gsd(c(r1, "r2")),
        "Unsupported operator or function `c`"
    )
    expect_error(
        stopping_rule_gsd(c()), # nolint: unnecessary_concatenation_linter.
        "Unsupported operator or function `c`"
    )
})

test_that("an injected list is read like list()", {
    conds <- list(quote(r1), "r1 and r2")
    rule <- stopping_rule_gsd(!!conds)
    expect_identical(rule$conditions, list(quote(r1), quote(r1 && r2)))
    expect_identical(rule$labels, c("r1", "r1 && r2"))
    expect_identical(rule$m_min, 2L)
    expect_identical(rule$n_interim, 2L)

    bad <- list(quote(r1), mean)
    expect_error(
        stopping_rule_gsd(!!bad),
        "Element 2 of `condition` must be an expression or a single string",
        fixed = TRUE
    )
    empty <- list()
    expect_error(stopping_rule_gsd(!!empty), "at least one condition")
})

# --- the grammar rejects ----------------------------------------------------

test_that("decision times, the stopping look and tables are refused", {
    expect_error(
        stopping_rule_gsd(r1 && t1 == 1),
        "decision times are not available to a stopping rule"
    )
    expect_error(
        stopping_rule_gsd(s == 1 && r1),
        "the stopping look is what the rule determines"
    )
    expect_error(
        stopping_rule_gsd(d(r1)),
        "Unsupported operator or function `d`"
    )
    expect_error(
        stopping_rule_gsd(d(r1)),
        "no discount tables in a stopping rule"
    )
})

test_that("r0 and zero-padded indices are not rejection indicators", {
    expect_error(stopping_rule_gsd(r0), "Symbol `r0` in a stopping condition")
    expect_error(stopping_rule_gsd(r0), "!!r0", fixed = TRUE)
    expect_error(
        stopping_rule_gsd(r1 && r01),
        "Symbol `r01` in a stopping condition"
    )
})

test_that("un-injected symbols error with the !! hint", {
    n_stop <- 2
    expect_error(
        stopping_rule_gsd(r1 + r2 >= n_stop),
        "not a rejection indicator"
    )
    expect_error(
        stopping_rule_gsd(r1 + r2 >= n_stop),
        "!!n_stop",
        fixed = TRUE
    )
    expect_error(
        stopping_rule_gsd(list(r1, "r2 && x")),
        "Symbol `x` in a stopping condition"
    )
})

test_that("vector literals, constant conditions and & are refused", {
    expect_error(
        stopping_rule_gsd(!!c(1, 2) * r1),
        "not a vector of length 2"
    )
    expect_error(
        stopping_rule_gsd(TRUE),
        "Stopping condition `TRUE` must refer to at least one rejection"
    )
    expect_error(
        stopping_rule_gsd(list(r1, "TRUE")),
        "Stopping condition `TRUE` must refer to at least one rejection"
    )
    expect_error(
        stopping_rule_gsd(r1 & r2),
        "Unsupported operator or function `&`"
    )
    expect_error(
        stopping_rule_gsd("r1 | r2"),
        "Unsupported operator or function `|`"
    )
})

test_that("inputs that are not conditions are refused", {
    f <- function(x) x
    expect_error(
        stopping_rule_gsd(!!f),
        "not a function",
        fixed = TRUE
    )
    expect_error(
        stopping_rule_gsd(list(r1, !!f)),
        "Element 2 of `condition` must be an expression or a single string",
        fixed = TRUE
    )
    expect_error(stopping_rule_gsd(list()), "at least one condition")
    expect_error(stopping_rule_gsd(!!character()), "at least one condition")
    expect_error(stopping_rule_gsd(), "`condition` is absent")
    expect_error(
        stopping_rule_gsd("r1 &&"),
        "Could not parse the stopping condition \"r1 &&\"",
        fixed = TRUE
    )
    expect_error(
        stopping_rule_gsd(list(r1, "r2 or")),
        "Could not parse the stopping condition \"r2 or\"",
        fixed = TRUE
    )
})

# --- settling against the p-values ------------------------------------------

test_that("a single condition is recycled to every interim analysis", {
    rule <- stopping_rule_gsd(r1 && r2)
    settled <- .gsd_check_stopping_rule_dims(rule, gsd_rule_pvals())
    expect_identical(settled, list(quote(r1 && r2), quote(r1 && r2)))

    k4 <- .gsd_check_stopping_rule_dims(rule, gsd_rule_pvals(K = 4L))
    expect_length(k4, 3L)

    # K = 1: no interim, so the rule is accepted and has no effect
    k1 <- .gsd_check_stopping_rule_dims(rule, gsd_rule_pvals(K = 1L))
    expect_identical(k1, list())
})

test_that("a list of K - 1 conditions passes, a list of K warns", {
    pvals <- gsd_rule_pvals()
    two <- stopping_rule_gsd(list(r1, r1 && r2))
    expect_identical(
        .gsd_check_stopping_rule_dims(two, pvals),
        list(quote(r1), quote(r1 && r2))
    )

    three <- stopping_rule_gsd(list(r1, r1 && r2, r3))
    expect_warning(
        settled <- .gsd_check_stopping_rule_dims(three, pvals),
        "The final analysis always ends the trial"
    )
    expect_identical(settled, list(quote(r1), quote(r1 && r2)))
    expect_warning(
        .gsd_check_stopping_rule_dims(three, pvals),
        "the last condition of `three`, `r3`, is ignored",
        fixed = TRUE
    )
})

test_that("any other number of conditions aborts naming L, K and K - 1", {
    pvals <- gsd_rule_pvals()
    short <- stopping_rule_gsd(list(r1))
    expect_error(
        .gsd_check_stopping_rule_dims(short, pvals),
        "must give one condition per interim analysis"
    )
    expect_error(
        .gsd_check_stopping_rule_dims(short, pvals),
        "It gives 1 condition."
    )

    long <- stopping_rule_gsd(list(r1, r2, r3, r1))
    cnd <- rlang::catch_cnd(.gsd_check_stopping_rule_dims(long, pvals))
    expect_s3_class(cnd, "rlang_error")
    msg <- gsub("\\s+", " ", conditionMessage(cnd))
    expect_match(msg, "It gives 4 conditions.", fixed = TRUE)
    expect_match(msg, "holds K = 3 analyses.", fixed = TRUE)
    expect_match(msg, "The rule needs 2 conditions", fixed = TRUE)
})

test_that("a rule beyond the hypotheses of the p-values aborts", {
    rule <- stopping_rule_gsd(r1 && r4)
    expect_error(
        .gsd_check_stopping_rule_dims(rule, gsd_rule_pvals()),
        "`rule` refers to hypothesis `r4` but `pvals` holds m = 3 hypotheses",
        fixed = TRUE
    )
    expect_error(
        .gsd_check_stopping_rule_dims(rule, gsd_rule_pvals(m = 1L)),
        "holds m = 1 hypothesis.",
        fixed = TRUE
    )
    expect_length(
        .gsd_check_stopping_rule_dims(rule, gsd_rule_pvals(m = 4L)),
        2L
    )
})

test_that("the settled conditions give the stopping truth table", {
    rule <- stopping_rule_gsd(r1 && (r2 || r3))
    settled <- .gsd_check_stopping_rule_dims(
        rule,
        gsd_rule_pvals(m = 3L, K = 2L)
    )
    expect_identical(
        .gsd_stop_table(settled, 3L, 2L),
        matrix(c(FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE), ncol = 1)
    )
})

# --- class, checks, methods -------------------------------------------------

test_that("class helpers", {
    rule <- stopping_rule_gsd(r1)
    expect_true(is_stopping_rule_gsd(rule))
    expect_false(is_stopping_rule_gsd(quote(r1)))
    expect_false(is_stopping_rule_gsd(NULL))

    expect_no_error(check_stopping_rule_gsd(rule))
    expect_error(
        check_stopping_rule_gsd("r1"),
        "must be a multigrain_stopping_rule_gsd object"
    )
    expect_error(
        check_stopping_rule_gsd(NULL),
        "must be a multigrain_stopping_rule_gsd object"
    )
    expect_no_error(check_stopping_rule_gsd(NULL, allow_null = TRUE))
    expect_no_error(check_stopping_rule_gsd(rule, allow_null = TRUE))
    expect_error(
        check_stopping_rule_gsd(1, allow_null = TRUE),
        "or `NULL`"
    )
})

test_that("print and summary methods", {
    recycled <- stopping_rule_gsd(r1 && (r2 || r3))
    expect_snapshot(print(recycled))
    expect_snapshot(summary(recycled))

    two <- stopping_rule_gsd(list(r1, "r1 and (r2 or r3)"))
    expect_snapshot(print(two))
    expect_snapshot(summary(two))

    expect_null(print.multigrain_stopping_rule_gsd(NULL))
    expect_null(summary.multigrain_stopping_rule_gsd(NULL))
})
