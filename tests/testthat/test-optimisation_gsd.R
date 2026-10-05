# `graph_optimise_gsd()` is `graph_optimise()` with the group sequential
# kernel and gain substituted in (design record, sections 4.5 and 4.7).

opt_pvals_gsd <- withr::with_seed(101, {
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

opt_gain_gsd <- trial_success_gsd(
    d(t1) + d(t2) + d(t3),
    d = c(1, 0.75),
    verbose = "silent"
)

opt_control <- multigrain_control() |>
    control_global(maxiter = 20, run = 5, popSize = 30) |>
    control_local(maxeval = 200) |>
    control_nsim_global(1000)

opt_result <- withr::with_seed(202, {
    graph_optimise_gsd(
        pvals = opt_pvals_gsd,
        graph_constraint = graph_constraint_free(3),
        trial_success = opt_gain_gsd,
        control = opt_control,
        verbose = "silent"
    )
})


test_that("the optimiser returns a populated multigrain_graph_optimal", {
    expect_s3_class(opt_result, "multigrain_graph_optimal")
    expect_true(is_graph_valid(
        opt_result$hyp_weight,
        opt_result$trans_matrix
    ))

    expect_type(opt_result$power$trial_success, "double")
    expect_length(opt_result$power$trial_success, 1L)

    expect_identical(dim(opt_result$power$local_power_by_analysis), c(3L, 2L))
    # a sum of proportions: equality holds to floating point, not bitwise
    # nolint start: expect_identical_linter
    expect_equal(
        rowSums(opt_result$power$time_distribution),
        c(H1 = 1, H2 = 1, H3 = 1)
    )
    # nolint end

    expect_true(is_trial_success_gsd(opt_result$trial_success))
    expect_identical(opt_result$trial_success, opt_gain_gsd)

    expect_true(opt_result$solution$opt_source %in% c("local", "global"))
    expect_true(opt_result$global_search)
})


test_that("global_search = FALSE runs the local stage only", {
    result <- withr::with_seed(303, {
        graph_optimise_gsd(
            pvals = opt_pvals_gsd,
            graph_constraint = graph_constraint_free(3),
            trial_success = opt_gain_gsd,
            global_search = FALSE,
            control = opt_control,
            verbose = "silent"
        )
    })

    expect_s3_class(result, "multigrain_graph_optimal")
    expect_false(result$global_search)
    expect_null(result$global_output)
    expect_identical(result$solution$opt_source, "local")
    expect_true(is_graph_valid(result$hyp_weight, result$trans_matrix))
})


test_that("a user start_graph is accepted", {
    start <- list(
        list(
            hyp_weight = c(0.6, 0.2, 0.2),
            trans_matrix = rbind(
                c(0, 0.5, 0.5),
                c(0.5, 0, 0.5),
                c(0.5, 0.5, 0)
            )
        )
    )

    result <- withr::with_seed(404, {
        graph_optimise_gsd(
            pvals = opt_pvals_gsd,
            graph_constraint = graph_constraint_free(3),
            trial_success = opt_gain_gsd,
            start_graph = start,
            global_search = FALSE,
            control = opt_control,
            verbose = "silent"
        )
    })

    expect_s3_class(result, "multigrain_graph_optimal")
    expect_identical(result$start_graph, start)
    expect_true(is_graph_valid(result$hyp_weight, result$trans_matrix))
})


test_that("print() and summary() of the result succeed", {
    expect_output(print(opt_result))
    expect_output(summary(opt_result))
})


test_that("print() and summary() survive m = 4, K = 3 with NA padding", {
    # design record, section 8 item 16
    pvals_4 <- withr::with_seed(505, {
        raw <- simulate_pvalues_gsd(
            power_nominal = c(0.9, 0.85, 0.8, 0.75),
            corr_matrix = diag(4),
            info_frac = rbind(
                c(1 / 3, 2 / 3, 1),
                c(0.5, 1, 1),
                c(NA, 0.6, 1),
                c(1, 1, 1)
            ),
            nsim = 500
        )
        transform_pvalues_gsd(
            raw,
            spending = gsDesign::sfLDOF,
            grid_size = 64L
        )
    })

    gain_4 <- trial_success_gsd(
        d(t1) + d(t2) + d(t3) + d(t4),
        d = c(1, 0.9, 0.75),
        verbose = "silent"
    )

    gc <- graph_constraint_free(4)
    w <- rep(0.25, 4L)
    trans <- matrix(1 / 3, nrow = 4L, ncol = 4L)
    diag(trans) <- 0

    power <- calc_power_pvals_gsd(
        pvals_4,
        hyp_weight = w,
        trans_matrix = trans,
        custom_power = list(trial_success = gain_4)
    )

    names(w) <- graph_constraint_get_names(gc)
    dimnames(trans) <- list(names(w), names(w))

    obj <- graph_optimal(
        hyp_weight = w,
        trans_matrix = trans,
        constraints = gc,
        trial_success = gain_4,
        power = power,
        solution = list(
            opt_source = "local",
            graph_valid = c(local = TRUE, global = NA)
        ),
        global_search = FALSE
    )

    expect_output(print(obj))
    expect_output(summary(obj))
    expect_output(print(pvals_4))
    expect_output(summary(pvals_4))
})


test_that(".sample_pvals_gsd() keeps the array and the class", {
    # design record, section 8 item 12
    sampled <- withr::with_seed(606, .sample_pvals_gsd(opt_pvals_gsd, 500L))

    expect_s3_class(sampled, "multigrain_pvals_gsd")
    expect_identical(dim(sampled$pvals), c(500L, 3L, 2L))
    expect_identical(sampled$nsim, 500L)
    expect_identical(sampled$alpha, opt_pvals_gsd$alpha)
    expect_identical(sampled$K, opt_pvals_gsd$K)
    expect_identical(sampled$tables, opt_pvals_gsd$tables)

    # rows come from the whole range, not just the first `nsim`
    idx <- withr::with_seed(606, sample.int(opt_pvals_gsd$nsim, 500L))
    expect_identical(
        sampled$pvals,
        opt_pvals_gsd$pvals[idx, , , drop = FALSE]
    )

    # K = 1 must not collapse to a matrix
    pvals_k1 <- withr::with_seed(707, {
        raw <- simulate_pvalues_gsd(
            power_nominal = c(0.9, 0.8),
            corr_matrix = diag(2),
            info_frac = 1,
            nsim = 200
        )
        transform_pvalues_gsd(raw, spending = gsDesign::sfLDOF)
    })
    small_k1 <- .sample_pvals_gsd(pvals_k1, 50L)
    expect_identical(dim(small_k1$pvals), c(50L, 2L, 1L))

    # m = 1 likewise
    pvals_m1 <- withr::with_seed(808, {
        raw <- simulate_pvalues_gsd(
            power_nominal = 0.9,
            corr_matrix = diag(1),
            info_frac = c(0.5, 1),
            nsim = 200
        )
        transform_pvalues_gsd(
            raw,
            spending = gsDesign::sfLDOF,
            grid_size = 64L
        )
    })
    small_m1 <- .sample_pvals_gsd(pvals_m1, 50L)
    expect_identical(dim(small_m1$pvals), c(50L, 1L, 2L))
})


test_that("control_prepare_dims() matches control_prepare()", {
    from_matrix <- control_prepare(
        multigrain_control(),
        matrix(0.5, nrow = 2000L, ncol = 3L)
    )
    from_dims <- control_prepare_dims(
        multigrain_control(),
        nsim = 2000L,
        m = 3L
    )
    expect_identical(from_dims, from_matrix)

    user_matrix <- control_prepare(
        opt_control,
        matrix(0.5, nrow = 2000L, ncol = 3L)
    )
    user_dims <- control_prepare_dims(opt_control, nsim = 2000L, m = 3L)
    expect_identical(user_dims, user_matrix)
})


test_that("mismatched dimensions and gain classes abort", {
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
    fixed_gain <- trial_success(r1 + r2 + r3, verbose = "silent")

    expect_error(
        graph_optimise_gsd(
            opt_pvals_gsd,
            graph_constraint_free(3),
            wrong_m,
            verbose = "silent"
        ),
        "number of hypotheses"
    )
    expect_error(
        graph_optimise_gsd(
            opt_pvals_gsd,
            graph_constraint_free(3),
            wrong_k,
            verbose = "silent"
        ),
        "number of analyses"
    )
    expect_error(
        graph_optimise_gsd(
            opt_pvals_gsd,
            graph_constraint_free(3),
            fixed_gain,
            verbose = "silent"
        ),
        "trial_success_gsd"
    )
    expect_error(
        graph_optimise_gsd(
            opt_pvals_gsd,
            graph_constraint_free(3),
            opt_gain_gsd,
            alpha = 0.05,
            verbose = "silent"
        ),
        "0.025"
    )
    expect_error(
        graph_optimise_gsd(
            opt_pvals_gsd$pvals[, , 1L],
            graph_constraint_free(3),
            opt_gain_gsd,
            verbose = "silent"
        ),
        "multigrain_pvals_gsd"
    )
})


# --- Figure 3b of the manuscript (design record, section 6 P4 gate) ---------
#
# Example 5: PFS matured at the interim, OS at information fractions
# (0.7, 1) with LDOF spending, rho = 0.5, nominal powers 0.98 and 0.93,
# N = 1e5 and one fixed seed. The gate is stated on the gain, not the argmax:
# the expected gain is flat near its maximum, so the argmax at N = 1e5 is
# noisy by up to 0.13 while the gain at the reference optimum is within
# 2.2e-4 of the grid maximum (record, Appendix B "Figure 3b argmax noise").

test_that("Figure 3b: the gain at the paper's optimum is a grid maximum", {
    skip_on_cran()

    ref <- readRDS(test_path("data", "gsd_example5_reference.rds"))
    inputs <- ref$inputs

    fig3b_pvals <- withr::with_seed(20260918, {
        raw <- simulate_pvalues_gsd(
            power_nominal = c(0.98, 0.93),
            corr_matrix = matrix(
                c(1, inputs$rho, inputs$rho, 1),
                nrow = 2
            ),
            info_frac = rbind(c(1, 1), c(inputs$t_os, 1)),
            alpha = inputs$alpha,
            nsim = 1e5
        )
        transform_pvalues_gsd(
            raw,
            spending = gsDesign::sfLDOF,
            alpha = inputs$alpha,
            grid_size = 1024L
        )
    })

    # the paper's non-centrality parameters
    expect_equal(
        calc_ncp(c(0.98, 0.93), alpha = inputs$alpha),
        c(inputs$Delta_P, inputs$Delta_O),
        tolerance = 1e-2
    )

    full_recycle <- matrix(c(0, 1, 1, 0), nrow = 2)
    w1_grid <- seq(0, 1, by = 0.005)

    # One kernel run per grid point; the 21 gains are then formed in R from
    # the decision-time distribution, so no gain is compiled per cell.
    time_dist <- lapply(w1_grid, function(w1) {
        calc_power_pvals_gsd(
            fig3b_pvals,
            hyp_weight = c(w1, 1 - w1),
            trans_matrix = full_recycle
        )$time_distribution
    })

    pick <- function(i, k) {
        vapply(time_dist, function(x) x[i, k + 1L], numeric(1L))
    }
    p_pfs <- list(pick(1L, 1L), pick(1L, 2L))
    p_os <- list(pick(2L, 1L), pick(2L, 2L))

    gain_curve <- function(r, delta) {
        v_pfs <- 1 / (1 + r)
        v_os <- r / (1 + r)
        v_pfs * (p_pfs[[1L]] + delta * p_pfs[[2L]]) +
            v_os * (p_os[[1L]] + delta * p_os[[2L]])
    }

    cells <- ref$w_star
    cells$gap <- NA_real_
    cells$w1_hat <- NA_real_
    for (j in seq_len(nrow(cells))) {
        curve <- gain_curve(cells$r[[j]], cells$delta[[j]])
        at_ref <- curve[which.min(abs(w1_grid - cells$w1_star[[j]]))]
        cells$gap[[j]] <- max(curve) - at_ref
        cells$w1_hat[[j]] <- w1_grid[[which.max(curve)]]
    }

    cat("\nFigure 3b: gain at the reference optimum vs the grid maximum\n")
    print(
        cells[, c("r", "delta", "w1_star", "w1_hat", "gap")],
        row.names = FALSE
    )

    for (j in seq_len(nrow(cells))) {
        expect_lte(
            cells$gap[[j]],
            5e-4,
            label = sprintf(
                "gain gap at r = %g, delta = %g",
                cells$r[[j]],
                cells$delta[[j]]
            )
        )
    }

    # reported, not asserted (record, section 6 P4 item 3)
    for (delta in ref$delta_grid) {
        sub <- cells[cells$delta == delta, ]
        sub <- sub[order(sub$r), ]
        cat(sprintf(
            "delta = %.2f: argmax non-increasing in r: %s\n",
            delta,
            all(diff(sub$w1_hat) <= 0)
        ))
    }

    # Three cells run through the compiled gain and the optimiser.
    three <- data.frame(r = c(1, 4, 8), delta = c(1, 0.75, 0.5))
    fig3b_control <- multigrain_control() |>
        control_global(popSize = 30, run = 20, maxiter = 100) |>
        control_nsim_global(2e4)
    gc_two <- graph_constraint(c(NA, NA), full_recycle)

    for (j in seq_len(nrow(three))) {
        r <- three$r[[j]]
        delta <- three$delta[[j]]
        v_pfs <- 1 / (1 + r)
        v_os <- r / (1 + r)

        gain <- trial_success_gsd(
            !!v_pfs * d(t1) + !!v_os * d(t2),
            d = c(1, delta),
            verbose = "silent"
        )
        curve <- gain_curve(r, delta)

        for (idx in c(21L, 141L)) {
            compiled <- calc_power_pvals_gsd(
                fig3b_pvals,
                hyp_weight = c(w1_grid[[idx]], 1 - w1_grid[[idx]]),
                trans_matrix = full_recycle,
                custom_power = list(g = gain)
            )$g
            expect_equal(compiled, curve[[idx]], tolerance = 1e-12)
        }

        optimised <- withr::with_seed(909 + j, {
            graph_optimise_gsd(
                pvals = fig3b_pvals,
                graph_constraint = gc_two,
                trial_success = gain,
                global_search = TRUE,
                control = fig3b_control,
                verbose = "silent"
            )
        })

        reference_w1 <- ref$w_star$w1_star[
            ref$w_star$r == r & ref$w_star$delta == delta
        ]
        cat(sprintf(
            "r = %g, delta = %.2f: optimised w_PFS = %.4f (reference %.4f)\n",
            r,
            delta,
            optimised$hyp_weight[[1L]],
            reference_w1
        ))

        expect_lte(
            abs(max(curve) - optimised$power$trial_success),
            5e-4,
            label = sprintf(
                "optimised gain at r = %g, delta = %g",
                r,
                delta
            )
        )
    }
})
