# Test validate_expr_symbols()
# should accept
test_that("validate_expr_symbols accepts pure r<digit> expressions", {
    expect_invisible(validate_expr_symbols(quote(r1 + r2 + r3)))
    expect_invisible(validate_expr_symbols(quote(r1 * r2)))
    expect_invisible(validate_expr_symbols(quote(r99)))
})

test_that("validate_expr_symbols accepts numeric literals alongside r<digit>", {
    expect_invisible(validate_expr_symbols(quote(2 * r1 + 0.5 * r2)))
    expect_invisible(validate_expr_symbols(quote(r1 + 3)))
})

test_that("validate_expr_symbols accepts expressions with parentheses", {
    expect_invisible(validate_expr_symbols(quote((r1 + r2) * r3)))
    expect_invisible(validate_expr_symbols(quote(((r1)))))
})

test_that("validate_expr_symbols errors on bare non-r symbol", {
    expect_error(
        validate_expr_symbols(quote(w * r1)),
        "in the trial success expression is not a rejection indicator"
    )
})

test_that("validate_expr_symbols error message suggests !! syntax", {
    expect_error(
        validate_expr_symbols(quote(alpha * r1 + r2)),
        "!!alpha",
        fixed = TRUE
    )
})

test_that("validate_expr_symbols errors on bare 'r' (no digits)", {
    expect_error(
        validate_expr_symbols(quote(r * r1)),
        "not a rejection indicator"
    )
})

test_that("validate_expr_symbols catches symbol nested in parentheses", {
    expect_error(
        validate_expr_symbols(quote(r1 * (w + r2))),
        "not a rejection indicator"
    )
})

test_that("validate_expr_symbols catches multiple bad symbols (first hit)", {
    expect_error(
        validate_expr_symbols(quote(a * r1 + b * r2)),
        "not a rejection indicator"
    )
})

test_that("validate_expr_symbols accepts bare numeric and logical literals", {
    expect_invisible(validate_expr_symbols(42))
    expect_invisible(validate_expr_symbols(3.14))
    expect_invisible(validate_expr_symbols(TRUE))
})

test_that("validate_expr_symbols errors on unexpected node types", {
    expect_error(
        validate_expr_symbols("a string"),
        "Unexpected element"
    )
})

test_that("validate_expr_symbols rejects unsupported functions", {
    expect_error(
        validate_expr_symbols(quote(sqrt(r1))),
        "Unsupported operator or function"
    )
    expect_error(
        validate_expr_symbols(quote(min(r1, r2))),
        "Unsupported operator or function"
    )
    expect_error(
        validate_expr_symbols(quote(log(r1 + r2))),
        "Unsupported operator or function"
    )
})


# --- resolve_expr: literal expressions ---
test_that("resolve_expr deparses a simple expression", {
    expect_identical(resolve_expr(quote(r1 + r2 + r3)), "r1 + r2 + r3")
})

test_that("resolve_expr deparses expression with numeric literals", {
    expect_identical(resolve_expr(quote(2 * r1 + r2)), "2 * r1 + r2")
})


# --- resolve_expr: post-!! expressions (simulating what enexpr produces) ---
test_that("resolve_expr handles scalar already inlined by !!", {
    # After enexpr(!!w * r1 + r2) with w = 3, the language object is 3 * r1 + r2
    expect_identical(resolve_expr(quote(3 * r1 + r2)), "3 * r1 + r2")
})

test_that("resolve_expr handles inlined value inside parentheses", {
    expect_identical(resolve_expr(quote(r1 * (0.5 + r2))), "r1 * (0.5 + r2)")
})

test_that("resolve_expr handles multiple inlined values", {
    expect_identical(
        resolve_expr(quote(2 * r1 + 0.5 * r2)),
        "2 * r1 + 0.5 * r2"
    )
})


test_that("resolve_expr handles fractional inlined value", {
    alpha <- 1 / 3
    # Simulate what enexpr(!!alpha * (r1 + r2)) produces
    expr <- bquote(.(alpha) * (r1 + r2))
    result <- resolve_expr(expr)
    expect_type(result, "character")
    expect_length(result, 1L)
    # Round-trip: the result should parse back to a valid expression
    expect_invisible(validate_expr_symbols(str2lang(result)))
})

# --- resolve_expr: string input ---
test_that("resolve_expr passes through string input", {
    expect_identical(resolve_expr("r1 + r2"), "r1 + r2")
})


# --- resolve_expr: error paths ---
test_that("resolve_expr errors on un-resolved symbol", {
    expect_error(
        resolve_expr(quote(w * r1)),
        "not a rejection indicator"
    )
})

test_that("resolve_expr error suggests !! syntax", {
    expect_error(
        resolve_expr(quote(w * r1)),
        "!!w",
        fixed = TRUE
    )
})

test_that("resolve_expr errors on non-language, non-string input", {
    expect_error(
        resolve_expr(42),
        "must be an expression or a character string"
    )
    expect_error(
        resolve_expr(TRUE),
        "must be an expression or a character string"
    )
    expect_error(
        resolve_expr(list(a = 1)),
        "must be an expression or a character string"
    )
    expect_error(
        resolve_expr(NULL),
        "must be an expression or a character string"
    )
})


# Test replace_r_indices
test_that("replace_r_indices correctly transforms indices", {
    # Test with a simple single index
    expect_identical(
        replace_r_indices("r1 + r2"),
        "double(x(i, 0)) + double(x(i, 1))"
    )

    # Test with multiple indices and operations
    # nolint start: line_length_linter
    expect_identical(
        replace_r_indices("r1*r10 - (r2 + r20) * 0.5"),
        "double(x(i, 0)) * double(x(i, 9)) - (double(x(i, 1)) + double(x(i, 19))) * 0.5"
    )
    # nolint end

    # Test with a complex expression including multiple occurrences of the
    # same index
    expect_identical(
        replace_r_indices("r1 + r1 * r3"),
        "double(x(i, 0)) + double(x(i, 0)) * double(x(i, 2))"
    )

    # Test input with no r indices
    expect_identical(
        replace_r_indices("x + y"),
        "x + y"
    )

    # Test with indices formatted with additional spaces (should not occur but
    # testing robustness)
    expect_identical(
        replace_r_indices("r1 +       r2 * 2"),
        "double(x(i, 0)) + double(x(i, 1)) * 2.0"
    )

    expect_identical(
        replace_r_indices("5      *       r2 || r3"),
        "5.0 * std_min(double(1), double(x(i, 1)) + double(x(i, 2)))"
    )

    # Test edge case of higher numbered indices
    expect_identical(
        replace_r_indices("r99 + r100"),
        "double(x(i, 98)) + double(x(i, 99))"
    )
})

# Stop forbidden Boolean operations.
test_that("Allowed boolean arithmetic transforms correctly", {
    # 1) bool * bool => bool
    expect_identical(
        replace_r_indices("r1 * r2"),
        "double(x(i, 0)) * double(x(i, 1))"
    )

    # 2) bool + bool => real
    expect_identical(
        replace_r_indices("r1 + r2"),
        "double(x(i, 0)) + double(x(i, 1))"
    )

    # 3) bool - bool => real
    expect_identical(
        replace_r_indices("r10 - r20"),
        "double(x(i, 9)) - double(x(i, 19))"
    )

    # 4) (bool * bool) - (bool * bool) => real
    expect_identical(
        replace_r_indices("r1*r10 - r2*r20"),
        "double(x(i, 0)) * double(x(i, 9)) - double(x(i, 1)) * double(x(i, 19))"
    )

    # 5) Mixing bool & real => real
    # r1 is bool, 5.0 is real => allowed
    expect_identical(
        replace_r_indices("r1 + 5"),
        "double(x(i, 0)) + 5.0"
    )

    # 6) Boolean operators with booleans: OK
    # nolint start: line_length_linter
    expect_identical(
        replace_r_indices("(r1 || r2) && (r3 || r4)"),
        "std_min(double(1), double(x(i, 0)) + double(x(i, 1))) * std_min(double(1), double(x(i, 2)) + double(x(i, 3)))"
    )
    # nolint end
})

test_that("Forbidden boolean arithmetic catches errors", {
    # 1) Boolean logic with a real => error
    # e.g. r1 is bool, 5.0 is real => (bool && real) => not allowed
    expect_error(
        replace_r_indices("r1 && 5"),
        "`&&` (AND) only allowed between booleans",
        fixed = TRUE
    )

    # 2) (real || bool) => error
    expect_error(
        replace_r_indices("5 || r1"),
        "`||` (OR) only allowed between booleans",
        fixed = TRUE
    )
})

test_that("replace_r_indices matrix index translation is correct", {
    expect_identical(replace_r_indices("r1"), "double(x(i, 0))")
    expect_identical(replace_r_indices("r2"), "double(x(i, 1))")
    expect_identical(replace_r_indices("r10"), "double(x(i, 9))")
})

test_that("replace_r_indices Logical OR is correctly transformed", {
    expect_identical(
        replace_r_indices("r1 || r2"),
        "std_min(double(1), double(x(i, 0)) + double(x(i, 1)))"
    )

    # nolint start: line_length_linter
    expect_identical(
        replace_r_indices("r1 || r2 + r3"),
        "std_min(double(1), double(x(i, 0)) + double(x(i, 1))) + double(x(i, 2))"
    )

    expect_identical(
        replace_r_indices("r1 + r2 || r3"),
        "double(x(i, 0)) + std_min(double(1), double(x(i, 1)) + double(x(i, 2)))"
    )
    # nolint end
})

test_that("replace_r_indices Logical AND is correctly transformed", {
    expect_identical(
        replace_r_indices("r1 && r2"),
        "double(x(i, 0)) * double(x(i, 1))"
    )
    expect_identical(
        replace_r_indices("r1 && r2 + r3"),
        "double(x(i, 0)) * double(x(i, 1)) + double(x(i, 2))"
    )
    expect_identical(
        replace_r_indices("r1 + r2 && r3"),
        "double(x(i, 0)) + double(x(i, 1)) * double(x(i, 2))"
    )
})

# nolint start: line_length_linter
test_that("replace_r_indices check for successful integer -> decimal transform", {
    expect_identical(
        replace_r_indices("1 / 3 * (r1 + r2)"),
        "1.0 / 3.0 * (double(x(i, 0)) + double(x(i, 1)))"
    )

    expect_identical(
        replace_r_indices("0.25 * (r1 || r2) || (r3 && r4)"),
        "0.25 * std_min(double(1), std_min(double(1), double(x(i, 0)) + double(x(i, 1))) + double(x(i, 2)) * double(x(i, 3)))"
    )
    # nolint end
})

test_that("replace_r_indices handles complex logical expressions correctly", {
    input_expr <- "r1 || r2 && r3 + r4"
    # nolint start: line_length_linter
    expected_output <- "std_min(double(1), double(x(i, 0)) + double(x(i, 1))) * double(x(i, 2)) + double(x(i, 3))"
    expect_identical(replace_r_indices(input_expr), expected_output)

    input_expr <- "(r1 || r2) * 5/2 && r3 + r4"
    expected_output <- "`&&` (AND) only allowed between booleans"
    expect_error(replace_r_indices(input_expr), expected_output, fixed = TRUE)

    input_expr <- "r1 || r2 || r3 + r4"
    expected_output <- "std_min(double(1), std_min(double(1), double(x(i, 0)) + double(x(i, 1))) + double(x(i, 2))) + double(x(i, 3))"
    expect_identical(replace_r_indices(input_expr), expected_output)

    input_expr <- "(r1 && r2) || (r3 || r4) + r5"
    expected_output <- "std_min(double(1), double(x(i, 0)) * double(x(i, 1)) + std_min(double(1), double(x(i, 2)) + double(x(i, 3)))) + double(x(i, 4))"
    expect_identical(replace_r_indices(input_expr), expected_output)

    input_expr <- "2 * (r1 || r2) && (r3 || r4) * 3 / 4 + r5"
    expected_output <- "2.0 * (std_min(double(1), double(x(i, 0)) + double(x(i, 1))) * std_min(double(1), double(x(i, 2)) + double(x(i, 3)))) * 3.0 / 4.0 + double(x(i, 4))"
    # nolint end

    expect_identical(replace_r_indices(input_expr), expected_output)
})


test_that("count_unique_indices identifies and handles indices correctly", {
    # Test with a straightforward case
    expect_identical(count_unique_indices("r1 + r2 + r3"), 3L)

    # Test with non-sequential indices
    expect_identical(count_unique_indices("r1 + r5 + r10"), 10L) |>
        suppressWarnings()

    # Test with missing indices and check for warnings
    expect_warning(
        count_unique_indices("r1 / r3 * r4"),
        "Missing indices in the sequence. Expected every index from 1 to 4, but missing 2" # nolint
    )

    # Test with repeated indices
    suppressWarnings(
        expect_identical(
            count_unique_indices(
                "r2 + r2 + r3"
            ),
            3L
        )
    )

    # check white space
    expect_identical(
        count_unique_indices(
            "r1 /r3 + r6 * r2 * r4 * r5"
        ),
        6L
    )
})


test_that("new_trial_success creates valid objects", {
    # Note: the compiled powerFunc now accepts LogicalMatrix, but Rcpp will
    # coerce a NumericMatrix argument. These tests pass numeric 0/1 matrices
    # and check behavioural equivalence -- the coercion is lossless because
    # all entries are in {0, 1}.
    test_matrix <- matrix(
        c(
            1, 0, 0,
            0, 1, 1,
            1, 1, 1,
            0, 0, 0,
            0, 0, 1,
            1, 0, 1
        ),
        ncol = 3,
        byrow = TRUE
    )

    test_expr <- "1/3 * (r1 + r2 + r3)"
    test_obj <- new_trial_success(test_expr, verbose = "silent")

    expect_s3_class(test_obj, "multigrain_trial_success")
    expect_identical(test_obj$m, 3L)
    expect_identical(test_obj$objective, test_expr)
    expect_identical(test_obj$func(test_matrix), 0.5)

    # test white space and unquoted
    # fmt: skip
    test_obj2 <- trial_success(r1 * r2 +        r3, verbose = "silent")
    expect_identical(test_obj2$func(test_matrix), 5 / 6)

    # test when expression is quoted
    test_obj3 <- trial_success("r1 + r2 + r3", verbose = "silent")
    expect_equal(test_obj3$func(test_matrix), 1.5)

    # test when expression is quoted and uses boolean logic
    test_obj4 <- trial_success("(r1 || r2) + r3", verbose = "silent")
    expect_identical(
        test_obj4$func(test_matrix),
        sum(apply(test_matrix, 1, function(x) (x[1] || x[2]) + x[3])) / 6
    )

    # test when expression when not quoted and uses boolean logic
    test_obj_4_5 <- trial_success(r1 || (r2 && r3), verbose = "silent")
    expect_identical(
        test_obj_4_5$func(test_matrix),
        sum(apply(test_matrix, 1, function(x) x[1] || (x[2] && x[3]))) / 6
    )

    # test when indices are missing
    expect_warning(test_obj5 <- trial_success("r2", verbose = "silent"))

    expect_identical(test_obj5$func(test_matrix), 2 / 6)
    expect_identical(test_obj5$m, 2L)
})


test_that("powerFunc matches numeric and logical matrix inputs", {
    # Regression guard for the LogicalMatrix template change:
    # numeric 0/1 input (coerced by Rcpp) and native logical input must agree.
    m_num <- matrix(
        c(
            1, 0, 0,
            0, 1, 1,
            1, 1, 1,
            0, 0, 0,
            0, 0, 1,
            1, 0, 1
        ),
        ncol = 3,
        byrow = TRUE
    )
    m_log <- matrix(
        as.logical(m_num),
        nrow = nrow(m_num),
        ncol = ncol(m_num)
    )

    # Cover each code path in the expression transformer
    exprs <- c(
        "r1 + r2 + r3",
        "r1 * r2 + r3",
        "r1 && r2 && r3",
        "r1 || r2 || r3",
        "(r1 && r2) || r3",
        "2 * r1 + 0.5 * r2 + r3"
    )

    for (e in exprs) {
        expect_snapshot({
            ts <- new_trial_success(e)
        })
        expect_identical(
            ts$func(m_num),
            ts$func(m_log),
            info = paste("expression:", e)
        )
    }
})

test_that("is_trial_success", {
    expect_true(
        is_trial_success(
            trial_success(r1 + r2 + r3 + r4, verbose = "silent")
        )
    )

    expect_false(
        is_trial_success("foo")
    )
})

test_that("check_trial_success", {
    expect_no_error(
        check_trial_success(
            trial_success(r1 + r2 + r3 + r4, verbose = "silent")
        )
    )

    expect_no_error(
        check_trial_success(
            NULL,
            allow_null = TRUE
        )
    )

    expect_error(
        check_trial_success(
            NULL
        ),
        "`NULL` must be a multigrain trial success object, not `NULL`."
    )

    expect_error(
        check_trial_success("foo"),
        '`"foo"` must be a multigrain trial success object'
    )

    expect_error(
        check_trial_success(2),
        "`2` must be a multigrain trial success object, not the number 2."
    )
})

test_that("trial_success print and summary methods", {
    ts <- trial_success(r1 + r2 + r3 + r4, verbose = "silent")

    expect_snapshot(print(ts))
    expect_null(print.multigrain_trial_success(NULL))

    expect_snapshot(summary(ts))
    expect_null(summary.multigrain_trial_success(NULL))
})

test_that("trial_success chatty", {
    # the new approach is with verbose as character
    expect_snapshot(
        trial_success(r1 + r2 + r3, verbose = "info")
    )

    expect_snapshot(
        trial_success(r1 + r2 + r3, verbose = "detail")
    )

    expect_snapshot(
        trial_success(r1 + r2 + r3, verbose = "silent")
    )

    # TRUE and FALSE still work
    expect_snapshot(
        trial_success(r1 + r2 + r3, verbose = TRUE)
    )

    expect_snapshot(
        trial_success(r1 + r2 + r3, verbose = FALSE)
    )

    expect_snapshot(error = TRUE, {
        trial_success(r1 + r2 + r3 + r4, verbose = 2)
    })
})

## E2E tests
## End-to-end: !! injection
test_that("trial_success end-to-end: !! injection", {
    w <- 2
    ts <- trial_success(!!w * r1 + r2, verbose = "silent")
    m <- rbind(
        c(1, 0),
        c(1, 1),
        c(1, 0)
    )
    expect_identical(ts$func(m), (2 + 2 + 1 + 2) / 3)
    expect_identical(ts$m, 2L)
})

test_that("trial_success errors on un-unquoted symbol", {
    w <- 2
    expect_error(
        trial_success(w * r1, verbose = "silent"),
        "in the trial success expression is not a rejection indicator"
    )
})

test_that("trial_success error suggests !! fix", {
    w <- 2
    expect_error(
        trial_success(w * r1, verbose = "silent"),
        "!!w",
        fixed = TRUE
    )
})

test_that("trial_success errors on expression with no r<digit> tokens", {
    expect_error(
        trial_success("4", verbose = "silent"),
        "must reference at least one hypothesis indicator"
    )
    expect_error(
        trial_success("x + y", verbose = "silent"),
        "must reference at least one hypothesis indicator"
    )
})

## r1 bug issue #131
test_that("r1 in global env no longer breaks trial_success(r1)", {
    had_r1 <- exists("r1", envir = .GlobalEnv, inherits = FALSE)
    old <- if (had_r1) get("r1", envir = .GlobalEnv, inherits = FALSE) else NULL
    on.exit({
        if (had_r1) {
            assign("r1", old, envir = .GlobalEnv)
        } else {
            rm("r1", envir = .GlobalEnv, inherits = FALSE)
        }
    })

    assign("r1", 999, envir = .GlobalEnv)

    ts <- trial_success(r1, verbose = "silent")
    m <- matrix(c(1, 0, 0, 1, 1, 1), ncol = 1, byrow = TRUE)
    expect_identical(ts$func(m), 2 / 3)
})
