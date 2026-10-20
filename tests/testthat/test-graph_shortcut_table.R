# Table kernel for graph_shortcut() / graph_shortcut_parallel(). The table
# kernel (precomputed critical values for every intersection) must give
# exactly the same rejections as the sequential kernel, which is the oracle.

shortcut_kernels <- c("sequential", "table", "auto")

# Every kernel through every backend: serial, num_threads = 1 and 2.
shortcut_all_kernels <- function(pvals, alpha, w, G) {
    out <- list()
    for (k in shortcut_kernels) {
        withr::with_options(list(multigrain.shortcut_kernel = k), {
            out[[paste0(k, "_serial")]] <- graph_shortcut(pvals, alpha, w, G)
            out[[paste0(k, "_t1")]] <- graph_shortcut_parallel(
                pvals,
                alpha,
                w,
                G,
                num_threads = 1
            )
            out[[paste0(k, "_t2")]] <- graph_shortcut_parallel(
                pvals,
                alpha,
                w,
                G,
                num_threads = 2
            )
        })
    }
    out
}

expect_all_identical <- function(outs) {
    for (nm in names(outs)[-1]) {
        testthat::expect_identical(outs[[nm]], outs[[1]], label = nm)
    }
}

pvals_high_power <- withr::with_seed(11, {
    z <- matrix(stats::rnorm(2000 * 6, mean = 3), ncol = 6)
    stats::pnorm(z, lower.tail = FALSE)
})

### AGREEMENT ###

test_that("forced kernels agree on random graphs (m = 2..6)", {
    for (m in 2:6) {
        for (seed in 1:4) {
            gr <- make_test_graph(m, seed = 100 * m + seed)
            expect_all_identical(shortcut_all_kernels(
                pvals_fixture[1:2000, 1:m],
                0.025,
                gr$w,
                gr$G
            ))
            expect_all_identical(
                shortcut_all_kernels(pvals_high_power[, 1:m], 0.025, gr$w, gr$G)
            )
        }
    }
})

test_that("forced kernels agree for m = 1", {
    outs <- shortcut_all_kernels(
        pvals_fixture[1:2000, 1, drop = FALSE],
        0.025,
        1,
        matrix(0, 1, 1)
    )
    expect_all_identical(outs)
    expect_identical(
        as.vector(outs[[1]]),
        pvals_fixture[1:2000, 1] < 0.025
    )
})

test_that("auto agrees with forced kernels either side of table_max_m", {
    withr::local_options(list(multigrain.shortcut_table_max_m = 3))
    for (m in 3:4) {
        gr <- make_test_graph(m, seed = m)
        expect_all_identical(
            shortcut_all_kernels(pvals_high_power[, 1:m], 0.025, gr$w, gr$G)
        )
    }
})

### EDGE CASES ###

test_that("p = 0 with level 0 is never rejected and the walk terminates", {
    w <- c(0, 0.5, 0.5)
    G <- rbind(c(0, 0.5, 0.5), c(0, 0, 1), c(0, 1, 0))
    pvals <- pvals_high_power[, 1:3]
    pvals[, 1] <- 0

    outs <- shortcut_all_kernels(pvals, 0.025, w, G)
    expect_all_identical(outs)
    expect_false(any(outs[[1]][, 1]))
    expect_true(any(outs[[1]][, 2]))
})

test_that("p equal to alpha * w at the full set is not rejected", {
    alpha <- 0.025
    w <- c(0.5, 0.25, 0.25)
    G <- rbind(c(0, 0.5, 0.5), c(0.5, 0, 0.5), c(0.5, 0.5, 0))
    pvals <- rbind(
        c(w[1] * alpha, 1, 1),
        c(1e-6, w[2] * alpha, 1),
        c(w[1] * alpha, w[2] * alpha, w[3] * alpha)
    )

    outs <- shortcut_all_kernels(pvals, alpha, w, G)
    expect_all_identical(outs)
    expect_identical(
        as_int_mat(outs[[1]]),
        rbind(c(0L, 0L, 0L), c(1L, 1L, 0L), c(0L, 0L, 0L))
    )
})

test_that("degenerate denominators (g_ij = g_ji = 1) agree", {
    G2 <- matrix(c(0, 1, 1, 0), 2, 2)
    expect_all_identical(
        shortcut_all_kernels(pvals_high_power[, 1:2], 0.025, c(0.5, 0.5), G2)
    )

    G4 <- rbind(
        c(0, 1, 0, 0),
        c(1, 0, 0, 0),
        c(0, 0, 0, 1),
        c(0.5, 0.5, 0, 0)
    )
    w4 <- c(0.25, 0.25, 0.25, 0.25)
    expect_all_identical(
        shortcut_all_kernels(pvals_high_power[, 1:4], 0.025, w4, G4)
    )
    expect_all_identical(
        shortcut_all_kernels(pvals_fixture[1:2000, 1:4], 0.025, w4, G4)
    )
})

test_that("zero weights and isolated nodes are never rejected", {
    w <- c(0.5, 0.5, 0, 0)
    G <- rbind(
        c(0, 0.5, 0.5, 0),
        c(1, 0, 0, 0),
        c(1, 0, 0, 0),
        c(0, 0, 0, 0)
    )
    pvals <- pvals_high_power[, 1:4]
    pvals[, 4] <- 1e-12

    outs <- shortcut_all_kernels(pvals, 0.025, w, G)
    expect_all_identical(outs)
    expect_false(any(outs[[1]][, 4]))
    expect_true(any(outs[[1]][, 3]))
})

test_that("rows of G summing to less than 1 agree", {
    for (m in 3:6) {
        gr <- make_test_graph(m, seed = 200 + m)
        G <- gr$G * 0.6
        expect_all_identical(
            shortcut_all_kernels(pvals_high_power[, 1:m], 0.025, gr$w, G)
        )
    }
})

test_that("weights near the objective's zeroing thresholds agree", {
    power_fn <- trial_success(r1 || r2 || r3, verbose = "silent")$func
    gc <- graph_constraint_free(3)
    pvals <- pvals_fixture[1:5000, 1:3]

    w_cases <- list(
        c(0.5, 0.5 - 0.99e-4, 0.99e-4),
        c(0.5, 0.5 - 1.01e-4, 1.01e-4)
    )
    G_cases <- list(
        rbind(c(0, 1 - 0.99e-5, 0.99e-5), c(0.5, 0, 0.5), c(0.5, 0.5, 0)),
        rbind(c(0, 1 - 1.01e-5, 1.01e-5), c(0.5, 0, 0.5), c(0.5, 0.5, 0))
    )

    for (w0 in w_cases) {
        for (G0 in G_cases) {
            x0 <- create_start_params(gc, w0 = w0, G0 = G0)
            vals <- list()
            for (k in shortcut_kernels) {
                for (nt in 1:2) {
                    obj <- create_obj_func(
                        3,
                        power_criterion = power_fn,
                        hyp_constraint = gc$hyp_constraint,
                        trans_constraint = gc$trans_constraint,
                        pvals = pvals,
                        num_threads = nt
                    )
                    vals[[paste0(k, "_", nt)]] <- withr::with_options(
                        list(multigrain.shortcut_kernel = k),
                        obj(x0)
                    )
                }
            }
            expect_all_identical(vals)
        }
    }
})

### OPTIONS AND KERNEL CHOICE ###

test_that("invalid multigrain.shortcut_kernel values give a clear error", {
    gr <- make_test_graph(3, seed = 1)
    pvals <- pvals_fixture[1:10, 1:3]
    for (bad in list("fast", NA_character_, c("auto", "table"), 1, TRUE)) {
        withr::with_options(list(multigrain.shortcut_kernel = bad), {
            expect_error(
                graph_shortcut(pvals, 0.025, gr$w, gr$G),
                "multigrain.shortcut_kernel"
            )
            expect_error(
                graph_shortcut_parallel(pvals, 0.025, gr$w, gr$G),
                "multigrain.shortcut_kernel"
            )
        })
    }
})

test_that("invalid multigrain.shortcut_table_max_m values give a clear error", {
    gr <- make_test_graph(3, seed = 1)
    pvals <- pvals_fixture[1:10, 1:3]
    for (bad in list(-1, 21, NA, NA_real_, "10", 2.5, TRUE, c(3, 4))) {
        withr::with_options(list(multigrain.shortcut_table_max_m = bad), {
            expect_error(
                graph_shortcut(pvals, 0.025, gr$w, gr$G),
                "multigrain.shortcut_table_max_m"
            )
        })
    }
})

test_that("kernel choice follows the options and the size rule", {
    # Defaults: table for small m at large N, sequential beyond the cap
    expect_identical(graph_shortcut_kernel(1e6, 3), "table")
    expect_identical(
        graph_shortcut_kernel(.Machine$integer.max, 21),
        "sequential"
    )
    expect_identical(graph_shortcut_kernel(1, 10), "sequential")
    expect_identical(graph_shortcut_kernel(1e6, 0), "sequential")

    withr::with_options(list(multigrain.shortcut_table_max_m = 0), {
        expect_identical(graph_shortcut_kernel(1e6, 1), "sequential")
        expect_identical(graph_shortcut_kernel(1e6, 3), "sequential")
    })
    withr::with_options(list(multigrain.shortcut_table_max_m = 5L), {
        expect_identical(graph_shortcut_kernel(1e6, 4), "table")
        expect_identical(graph_shortcut_kernel(1e6, 5), "table")
        expect_identical(graph_shortcut_kernel(1e6, 6), "sequential")
    })
    withr::with_options(list(multigrain.shortcut_kernel = "sequential"), {
        expect_identical(graph_shortcut_kernel(1e6, 3), "sequential")
    })
    withr::with_options(list(multigrain.shortcut_kernel = "table"), {
        expect_identical(graph_shortcut_kernel(1, 10), "table")
        # table_max_m applies to "auto" only
        withr::with_options(list(multigrain.shortcut_table_max_m = 0), {
            expect_identical(graph_shortcut_kernel(1, 3), "table")
        })
    })
})

test_that("auto uses the table kernel exactly when 2^m <= 2 * N", {
    # TABLE_KAPPA = 2, see the comment on it in src/graph_shortcut.cpp
    for (m in c(8, 14, 17, 20)) {
        expect_identical(graph_shortcut_kernel(2^(m - 1), m), "table")
        expect_identical(graph_shortcut_kernel(2^(m - 1) - 1, m), "sequential")
    }
})

test_that("forcing the table kernel above m = 20 errors before allocating", {
    withr::local_options(list(multigrain.shortcut_kernel = "table"))
    expect_error(graph_shortcut_kernel(10, 21), "at most 20")
    pvals <- matrix(0.5, nrow = 2, ncol = 21)
    G <- matrix(1 / 20, 21, 21)
    diag(G) <- 0
    expect_error(
        graph_shortcut(pvals, 0.025, rep(1 / 21, 21), G),
        "at most 20"
    )
    expect_error(
        graph_shortcut_parallel(pvals, 0.025, rep(1 / 21, 21), G),
        "at most 20"
    )
})

### TABLE BUILDER ###

test_that("graph_shortcut_table() matches calc_local_weights() * alpha", {
    alpha <- 0.025
    for (m in 2:5) {
        gr <- make_test_graph(m, seed = 300 + m)
        tab <- graph_shortcut_table(alpha, gr$w, gr$G)
        expect_identical(dim(tab), as.integer(c(2^m, m)))
        expect_identical(tab[1, ], rep(0, m))

        lw <- calc_local_weights(gr$w, gr$G)
        member <- lw[, seq_len(m), drop = FALSE]
        mask <- as.vector(member %*% 2^(seq_len(m) - 1))
        expect_equal(
            unname(tab[mask + 1, , drop = FALSE]),
            unname(lw[, m + seq_len(m), drop = FALSE]) * alpha,
            tolerance = 1e-12
        )
    }
})

test_that("graph_shortcut_table() full-set row is exactly w * alpha", {
    gr <- make_test_graph(4, seed = 7)
    tab <- graph_shortcut_table(0.025, gr$w, gr$G)
    expect_identical(tab[2^4, ], unname(gr$w) * 0.025)
})
