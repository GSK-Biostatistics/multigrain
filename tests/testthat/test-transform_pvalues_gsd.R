# Reference values are Maurer and Bretz (2013) Tables 1 and 2, reproduced from
# `gsDesign` in Appendix A of `dev/gsd_design_record.md`.

test_that("Table 1 boundaries are reproduced to four significant figures", {
    info <- c(1 / 3, 2 / 3, 1)

    # the last row of a boundary table is the boundary at the full level, so
    # building one table per level of Table 1 reproduces the published rows
    bounds_at <- function(level) {
        bound_table <- .gsd_boundary_table(
            info,
            spending = gsDesign::sfLDOF,
            alpha = level,
            grid_size = 16L,
            hyp = 1L
        )
        signif(bound_table$bounds[nrow(bound_table$bounds), ], 4)
    }

    expect_equal(bounds_at(0.0125), c(1.517e-05, 2.215e-03, 1.180e-02))
    expect_equal(bounds_at(0.01875), c(4.679e-05, 3.976e-03, 1.750e-02))
    expect_equal(bounds_at(0.00625), c(2.179e-06, 8.105e-04, 5.986e-03))
    expect_equal(bounds_at(0.025), c(1.035e-04, 6.012e-03, 2.313e-02))
})

test_that("Table 2 repeated p-values are reproduced by the transform", {
    info <- c(1 / 3, 2 / 3, 1)
    raw <- array(0.5, dim = c(1L, 4L, 3L))
    raw[1L, , 1L] <- c(0.0062, 0.017, 0.009, 0.13)
    raw[1L, , 2L] <- c(0.0002, 0.0035, 0.002, 0.06)

    # the published repeated p-values run up to 0.382, so the tables have to be
    # built over levels reaching that far: anything above `alpha` maps to 1
    out <- transform_pvalues_gsd(
        raw,
        info_frac = info,
        spending = gsDesign::sfLDOF,
        alpha = 0.5,
        grid_size = 1024L
    )

    expect_equal(
        signif(out$pvals[1L, , 1L], 4),
        c(0.1141, 0.1682, 0.1315, 0.3820)
    )
    # Appendix A reports 0.002393 0.017161 0.011649 0.128460 at analysis 2
    expect_equal(
        signif(out$pvals[1L, , 2L], 4),
        c(0.002393, 0.01716, 0.01165, 0.1285)
    )
})

test_that("a spending function that is not well ordered aborts", {
    # switching family with the level makes the boundary at analysis 1 fall as
    # the allocated level rises, which is condition (2) of Maurer and Bretz
    not_well_ordered <- function(a, t) {
        if (a < 0.005) {
            gsDesign::sfLDPocock(a, t)$spend
        } else {
            gsDesign::sfLDOF(a, t)$spend
        }
    }

    expect_error(
        transform_pvalues_gsd(
            array(0.01, dim = c(1L, 2L, 3L)),
            info_frac = c(1 / 3, 2 / 3, 1),
            spending = list(gsDesign::sfLDOF, not_well_ordered),
            grid_size = 64L
        ),
        "hypothesis 2"
    )
})

test_that("the grid inverse matches a uniroot inverse", {
    info <- c(0.7, 1)
    alpha <- 0.025

    nominal_bounds <- function(level) {
        increments <- diff(c(0, gsDesign::sfLDOF(level, info)$spend))
        bounds <- gsDesign::gsBound1(
            theta = 0,
            I = info,
            a = rep(-20, length(info)),
            probhi = increments
        )$b
        stats::pnorm(bounds, lower.tail = FALSE)
    }
    exact_inverse <- function(p, k) {
        stats::uniroot(
            function(level) nominal_bounds(level)[[k]] - p,
            c(1e-15, 0.999),
            tol = 1e-13
        )$root
    }

    set.seed(2)
    p_test <- 10^stats::runif(200, -9, log10(0.03))
    bound_table <- .gsd_boundary_table(
        info,
        spending = gsDesign::sfLDOF,
        alpha = alpha,
        grid_size = 1024L,
        hyp = 1L
    )

    for (k in seq_along(info)) {
        exact <- vapply(p_test, exact_inverse, double(1L), k = k)
        approximate <- .gsd_invert(
            bound_table$bounds[, k],
            bound_table$grid,
            p_test
        )
        # levels above alpha are outside the table and map to 1 by design
        keep <- exact < alpha
        expect_lt(
            max(abs(exact[keep] - approximate[keep]) / exact[keep]),
            1e-5
        )
    }
})

test_that("a single full-information analysis returns the raw p-values", {
    set.seed(11)
    raw <- array(stats::runif(20), dim = c(10L, 2L, 1L))

    out <- transform_pvalues_gsd(
        raw,
        info_frac = 1,
        spending = gsDesign::sfLDOF,
        grid_size = 8L
    )

    expect_identical(out$pvals[, 1L, 1L], raw[, 1L, 1L])
    expect_identical(out$pvals[, 2L, 1L], raw[, 2L, 1L])
    expect_null(out$tables[[1L]]$bounds)
})

test_that("matured p-values are copied forward, with a warning", {
    raw <- array(0.001, dim = c(2L, 1L, 3L))
    raw[, 1L, 2L] <- 0.02
    raw[, 1L, 3L] <- 0.03

    expect_warning(
        out <- transform_pvalues_gsd(
            raw,
            info_frac = c(0.5, 1, 1),
            spending = gsDesign::sfLDOF,
            grid_size = 256L
        ),
        "hypothesis 1"
    )

    # the hypothesis matures at analysis 2, so analysis 3 reuses its value and
    # the raw p-value of 0.03 is ignored
    expect_identical(out$tables[[1L]]$maturity, 2L)
    expect_identical(out$pvals[, 1L, 3L], out$pvals[, 1L, 2L])

    quiet <- array(0.001, dim = c(2L, 1L, 2L))
    expect_no_warning(transform_pvalues_gsd(
        quiet,
        info_frac = c(1, 1),
        spending = gsDesign::sfLDOF,
        grid_size = 8L
    ))
})

test_that("analyses without data give 1, and look-back carries across them", {
    raw <- array(c(3e-4, 8e-4), dim = c(2L, 1L, 3L))

    # no data at analysis 1: the hypothesis cannot be rejected there, and with
    # look_back = FALSE analysis 3 uses analysis 3 evidence only
    late <- transform_pvalues_gsd(
        raw,
        info_frac = matrix(c(NA, 0.5, 1), nrow = 1L),
        spending = gsDesign::sfLDOF,
        grid_size = 256L
    )
    expect_identical(late$pvals[, 1L, 1L], c(1, 1))
    expect_true(all(late$pvals[, 1L, 3L] < late$pvals[, 1L, 2L]))

    # no data at analysis 2: with look_back = FALSE it cannot be rejected
    # there, and with look_back = TRUE the analysis 1 sequential p-value is
    # carried into analysis 2
    gap_rep <- transform_pvalues_gsd(
        raw,
        info_frac = matrix(c(0.5, NA, 1), nrow = 1L),
        spending = gsDesign::sfLDOF,
        grid_size = 256L
    )
    expect_identical(gap_rep$pvals[, 1L, 2L], c(1, 1))

    gap_seq <- transform_pvalues_gsd(
        raw,
        info_frac = matrix(c(0.5, NA, 1), nrow = 1L),
        spending = gsDesign::sfLDOF,
        look_back = TRUE,
        grid_size = 256L
    )
    expect_identical(gap_seq$pvals[, 1L, 2L], gap_rep$pvals[, 1L, 1L])
    expect_identical(
        gap_seq$pvals[, 1L, 3L],
        pmin(gap_rep$pvals[, 1L, 1L], gap_rep$pvals[, 1L, 3L])
    )
})

test_that("the spending function must spend the whole level by t = 1", {
    stops_short <- function(a, t) 0.5 * a * t

    expect_error(
        transform_pvalues_gsd(
            array(0.01, dim = c(1L, 1L, 2L)),
            info_frac = c(0.5, 1),
            spending = stops_short,
            grid_size = 8L
        ),
        "hypothesis 1"
    )
})

test_that("plain-numeric spending and information fractions above 1 work", {
    # cumulative spend returned as a plain numeric vector, with no `$spend`
    plain_of <- function(a, t) {
        pmin(2 * (1 - stats::pnorm(stats::qnorm(1 - a / 2) / sqrt(t))), a)
    }

    plain <- transform_pvalues_gsd(
        array(0.001, dim = c(1L, 1L, 2L)),
        info_frac = c(0.5, 1),
        spending = plain_of,
        grid_size = 128L
    )
    expect_true(all(plain$pvals > 0 & plain$pvals <= 1))

    over_one <- transform_pvalues_gsd(
        array(0.001, dim = c(1L, 1L, 2L)),
        info_frac = c(0.5, 1.2),
        spending = gsDesign::sfLDOF,
        grid_size = 128L
    )
    expect_identical(over_one$tables[[1L]]$maturity, 2L)
    expect_true(all(over_one$pvals > 0 & over_one$pvals <= 1))
})

test_that("a single interim analysis inflates the p-value", {
    raw <- array(c(5e-4, 1e-3), dim = c(2L, 1L, 1L))

    out <- transform_pvalues_gsd(
        raw,
        info_frac = 0.5,
        spending = gsDesign::sfLDOF,
        grid_size = 512L
    )

    # the boundary at an interim analysis is the partial spend, not the level
    expect_true(all(out$pvals[, 1L, 1L] > raw[, 1L, 1L]))
})

test_that("repeated p-values are floored at 1e-14, never at zero", {
    raw <- array(NA_real_, dim = c(3L, 1L, 3L))
    raw[, 1L, 1L] <- c(5e-11, 1e-100, 0)
    raw[, 1L, 2L] <- c(5e-11, 1e-100, 0)
    raw[, 1L, 3L] <- c(5e-11, 1e-100, 0)

    out <- transform_pvalues_gsd(
        raw,
        info_frac = c(1 / 3, 2 / 3, 1),
        spending = gsDesign::sfLDOF,
        grid_size = 1024L
    )

    expect_true(all(out$pvals > 0))
    expect_true(all(out$pvals >= 1e-14))
    # p = 5e-11 at an interim analysis is inside the table, not floored
    expect_gt(out$pvals[1L, 1L, 1L], 1e-14)
    # p-values below the whole table take the floor rather than zero
    expect_identical(out$pvals[2L, 1L, 1L], 1e-14)
    expect_identical(out$pvals[3L, 1L, 1L], 1e-14)

    # a p-value at a single full-information analysis maps to itself
    single <- transform_pvalues_gsd(
        array(5e-11, dim = c(1L, 1L, 1L)),
        info_frac = 1,
        spending = gsDesign::sfLDOF,
        grid_size = 8L
    )
    expect_identical(single$pvals[1L, 1L, 1L], 5e-11)
})

test_that("the transform agrees with graphicalMCP", {
    skip_if_not_installed("graphicalMCP")

    # Convention exercised: every hypothesis has the same information
    # fractions, both analyses carry data, and none matures before the final
    # analysis, so graphicalMCP's NA padding and its matured-hypothesis
    # convention (section 4.9 of the design record) never come into play.
    # `look_back = FALSE` is `repeated_p()`; `look_back = TRUE` is
    # `sequential_p()`. `spending_of()` returns a plain numeric vector.
    info <- c(0.5, 1)
    spending <- graphicalMCP::spending_of

    set.seed(21)
    raw <- array(NA_real_, dim = c(5L, 2L, 2L))
    raw[, , 1L] <- 10^stats::runif(10, -5, -3.2)
    raw[, , 2L] <- 10^stats::runif(10, -4, -2)

    repeated <- transform_pvalues_gsd(
        raw,
        info_frac = info,
        spending = spending,
        grid_size = 1024L
    )
    sequential <- transform_pvalues_gsd(
        raw,
        info_frac = info,
        spending = spending,
        look_back = TRUE,
        grid_size = 1024L
    )

    for (n in seq_len(dim(raw)[[1L]])) {
        for (i in seq_len(dim(raw)[[2L]])) {
            oracle_first <- graphicalMCP::repeated_p(
                raw[n, i, 1L],
                info[[1L]],
                spending
            )
            oracle_last <- graphicalMCP::repeated_p(
                raw[n, i, ],
                info,
                spending
            )
            oracle_seq <- graphicalMCP::sequential_p(
                raw[n, i, ],
                info,
                spending
            )

            expect_lt(abs(repeated$pvals[n, i, 1L] - oracle_first), 1e-5)
            expect_lt(abs(repeated$pvals[n, i, 2L] - oracle_last), 1e-5)
            expect_lt(abs(sequential$pvals[n, i, 2L] - oracle_seq), 1e-5)
        }
    }
})

test_that("print() and summary() work with NA padding", {
    set.seed(4)
    raw <- array(stats::runif(3 * 6 * 3), dim = c(3L, 6L, 3L))
    raw[, 2L, 3L] <- raw[, 2L, 2L]
    raw[, 3L, 2L] <- raw[, 3L, 1L]
    raw[, 3L, 3L] <- raw[, 3L, 1L]
    info <- rbind(
        c(1 / 3, 2 / 3, 1),
        c(0.5, 1, 1),
        c(1, 1, 1),
        c(NA, 0.5, 1),
        c(0.5, NA, 1),
        c(NA, NA, 1)
    )

    out <- transform_pvalues_gsd(
        raw,
        info_frac = info,
        spending = gsDesign::sfLDOF,
        look_back = c(TRUE, FALSE, FALSE, TRUE, TRUE, FALSE),
        grid_size = 32L
    )

    expect_output(print(out), "multigrain_pvals_gsd")
    expect_output(print(out), "H6")
    expect_output(summary(out), "Nominal boundaries")
    expect_output(summary(out), "matured_at")
})
