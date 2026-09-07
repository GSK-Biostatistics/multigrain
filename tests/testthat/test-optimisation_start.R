test_that(".build_start_matrix produces a valid output - no constraints", {
    gc <- graph_constraint_free(4)
    start_graph <- list(
        list(
            hyp_weight = c(0.1, 0.5, 0.4, 0),
            trans_matrix = gMCPLite::parallelGatekeeping()@m
        ),
        list(
            hyp_weight = c(0.6, 0.2, 0.2, 0.2),
            trans_matrix = gMCPLite::parallelGatekeeping()@m
        ),
        list(
            hyp_weight = c(0.6, 0.2, 0.2, 0.2),
            trans_matrix = gMCPLite::generalSuccessive(
                gamma = 0.4,
                delta = 0.2
            )@m
        )
    )

    s <- .build_start_matrix(gc, start_graph)
    expect_true(is.matrix(s))
    expect_gt(nrow(s), 1)
    expect_identical(ncol(s), length(create_start_params(gc)))
    expect_true(all(is.finite(s)))
    # no duplicate rows
    expect_identical(nrow(unique(s)), nrow(s))
    expect_identical(nrow(s), 5L)
    expect_identical(ncol(s), 11L)
})

test_that(".build_start_matrix produces a valid output - WITH constraints", {
    gc <- graph_constraint(
        hyp_constraint = c(NA, 0.2, 0.2, NA),
        trans_constraint = rbind(
            c(0, 0.8, NA, NA),
            c(0.8, 0, NA, NA),
            c(0.8, NA, 0, NA),
            c(0.8, NA, NA, 0)
        )
    )
    start_graph <- list(
        list(
            hyp_weight = c(0.1, 0.5, 0.4, 0),
            trans_matrix = gMCPLite::parallelGatekeeping()@m
        ),
        list(
            hyp_weight = c(0.6, 0.2, 0.2, 0.2),
            trans_matrix = gMCPLite::parallelGatekeeping()@m
        ),
        list(
            hyp_weight = c(0.6, 0.2, 0.2, 0.2),
            trans_matrix = gMCPLite::generalSuccessive(
                gamma = 0.4,
                delta = 0.2
            )@m
        )
    )

    s <- .build_start_matrix(gc, start_graph)
    expect_true(is.matrix(s))
    expect_gt(nrow(s), 1)
    expect_identical(ncol(s), length(create_start_params(gc)))
    expect_true(all(is.finite(s)))
    # no duplicate rows
    expect_identical(nrow(unique(s)), nrow(s))
    expect_identical(nrow(s), 5L)
    expect_identical(ncol(s), 5L)
})

test_that(".build_start_matrix with default start graph", {
    gc <- graph_constraint_free(5)

    expect_snapshot({
        .build_start_matrix(
            gc,
            start_graph = NULL
        )
    })

    expect_snapshot({
        .build_start_matrix(
            gc,
            start_graph = list(
                list(
                    hyp_weight = NULL,
                    trans_matrix = NULL
                )
            )
        )
    })
})

test_that(".is_default_start_graph", {
    expect_true(
        .is_default_start_graph(NULL)
    )

    expect_true(
        .is_default_start_graph(
            list(
                list(
                    hyp_weight = NULL,
                    trans_matrix = NULL
                )
            )
        )
    )

    expect_false(
        .is_default_start_graph("foo")
    )
})

test_that(".validate_start_graphs with default start graphs", {
    expect_no_error(
        .validate_start_graphs(NULL)
    )

    expect_no_error(
        .validate_start_graphs(
            list(
                list(
                    hyp_weight = NULL,
                    trans_matrix = NULL
                )
            )
        )
    )
})

test_that(".validate_start_graphs with other start graphs", {
    hyp_w <- c(0.1, 0.2, NA, NA, NA)
    trans_m <- matrix(
        rep_len(
            c(0.1, NA),
            length.out = 25
        ),
        nrow = 5,
        ncol = 5
    )
    diag(trans_m) <- 0

    expect_no_error(
        .validate_start_graphs(
            list(
                list(
                    hyp_weight = hyp_w,
                    trans_matrix = trans_m
                )
            ),
            m = 5
        )
    )

    expect_snapshot(error = TRUE, {
        .validate_start_graphs(
            list(
                list(
                    hyp_weight = c(0.1, 0.2, NA, NA),
                    trans_matrix = trans_m
                )
            ),
            m = 5
        )
    })

    expect_snapshot(error = TRUE, {
        .validate_start_graphs(
            list(
                list(
                    hyp_weight = hyp_w,
                    trans_matrix = hyp_w
                )
            ),
            m = 5
        )
    })

    trans_m_4 <- matrix(
        rep_len(
            c(0.1, NA),
            length.out = 16
        ),
        nrow = 4,
        ncol = 4
    )
    diag(trans_m) <- 0

    expect_snapshot(error = TRUE, {
        .validate_start_graphs(
            list(
                list(
                    hyp_weight = hyp_w,
                    trans_matrix = trans_m_4
                )
            ),
            m = 5
        )
    })
})


# ---- .encode_graph ----

# The seven-edge reference graph of the design record's appendix check 9.
simplify_ref_graph <- function() {
    list(
        hyp_weight = c(0.9951, 0, 0, 0.0049),
        trans_matrix = rbind(
            c(0, 0.6567, 0.3433, 0),
            c(0.2834, 0, 0.3354, 0.3812),
            c(1, 0, 0, 0),
            c(0, 1, 0, 0)
        )
    )
}

test_that(".encode_graph guards derived entries that should be zero", {
    gc <- graph_constraint_free(4)
    ref <- simplify_ref_graph()
    expect_true(is_graph_valid(ref$hyp_weight, ref$trans_matrix))

    x <- .encode_graph(gc, ref$hyp_weight, ref$trans_matrix)
    expect_type(x, "double")
    expect_null(attributes(x))
    expect_length(x, length(create_start_params(gc)))

    theta <- split_theta(x, gc$hyp_constraint)
    w_dec <- recover_full_weights(theta$w_pars, gc$hyp_constraint)
    G_dec <- recover_full_trans_matrix(theta$g_pars, gc$trans_constraint)

    # Nothing decodes negative, which is the failure mode the guard removes.
    expect_true(all(G_dec >= 0))
    expect_true(all(w_dec >= 0))

    # Derived entry of each row: rows 1-3 end at column 4, row 4 at column 3.
    derived <- c(G_dec[1, 4], G_dec[2, 4], G_dec[3, 4], G_dec[4, 3])
    zero_in_ref <- c(TRUE, FALSE, TRUE, TRUE)
    expect_true(all(abs(derived[zero_in_ref] - 5e-6) < 1e-12))
    # The one derived entry that is genuinely non-zero is left alone.
    expect_equal(derived[2], 0.3812, tolerance = 1e-9)

    # The derived weight is above 1e-4, so the weight guard must not fire.
    expect_equal(w_dec[4], 0.0049, tolerance = 1e-9)

    # Same edge count before and after: the guarded entries sit below the
    # 1e-5 threshold, so they are not counted and are snapped to exact zero.
    G_thresh <- G_dec
    G_thresh[G_thresh < 1e-5] <- 0
    expect_identical(
        sum(G_thresh[is.na(gc$trans_constraint)] != 0),
        sum(ref$trans_matrix[is.na(gc$trans_constraint)] != 0)
    )
})

test_that(".encode_graph leaves the objective value unchanged", {
    gc <- graph_constraint_free(4)
    ref <- simplify_ref_graph()
    ts <- trial_success(
        0.25 * (2 * (r1 && r2) + r1 * r3 + r2 * r4),
        verbose = "silent"
    )
    pv <- withr::with_seed(2, {
        corr <- matrix(0.2, 4, 4)
        diag(corr) <- 1
        sims <- mvtnorm::rmvnorm(
            1e4,
            mean = calc_ncp(c(0.93, 0.91, 0.90, 0.85)),
            sigma = corr
        )
        stats::pnorm(sims, lower.tail = FALSE)
    })

    obj <- create_obj_func(
        4,
        power_criterion = ts$func,
        hyp_constraint = gc$hyp_constraint,
        trans_constraint = gc$trans_constraint,
        pvals = pv
    )

    plain <- as.numeric(create_start_params(
        gc,
        w0 = ref$hyp_weight,
        G0 = ref$trans_matrix,
        sum_to_one_constraint = FALSE
    ))
    guarded <- .encode_graph(gc, ref$hyp_weight, ref$trans_matrix)

    expect_identical(obj(guarded), obj(plain))
    # ... and both equal the trial success of the reference itself.
    u_ref <- ts$func(
        graph_shortcut(pv, 0.025, ref$hyp_weight, ref$trans_matrix)
    )
    expect_identical(obj(guarded), u_ref)
})

test_that(".encode_graph guards a derived hypothesis weight of zero", {
    gc <- graph_constraint_free(3)
    # w[3] is the derived weight and is zero in the graph being encoded.
    w <- c(0.6, 0.4, 0)
    G <- rbind(c(0, 0.5, 0.5), c(0.5, 0, 0.5), c(0.5, 0.5, 0))
    expect_true(is_graph_valid(w, G))

    x <- .encode_graph(gc, w, G)
    theta <- split_theta(x, gc$hyp_constraint)
    w_dec <- recover_full_weights(theta$w_pars, gc$hyp_constraint)

    expect_true(all(w_dec >= 0))
    expect_true(abs(w_dec[3] - 5e-5) < 1e-12)
    # Below the 1e-4 weight threshold, so it is still zeroed by the objective.
    expect_lt(w_dec[3], 1e-4)
})

test_that(".encode_graph respects a row's pinned non-zero entries", {
    gc <- graph_constraint(
        hyp_constraint = c(NA, NA, NA, NA),
        trans_constraint = rbind(
            c(0, 0.5, NA, NA),
            c(NA, 0, NA, NA),
            c(NA, NA, 0, NA),
            c(NA, NA, NA, 0)
        )
    )
    # Row 1's derived entry (column 4) is zero and must be guarded to 5e-6
    # while G[1, 2] stays pinned at 0.5.
    w <- c(0.25, 0.25, 0.25, 0.25)
    G <- rbind(
        c(0, 0.5, 0.5, 0),
        c(1 / 3, 0, 1 / 3, 1 / 3),
        c(1 / 3, 1 / 3, 0, 1 / 3),
        c(1 / 3, 1 / 3, 1 / 3, 0)
    )
    expect_true(is_graph_valid(w, G))

    x <- .encode_graph(gc, w, G)
    theta <- split_theta(x, gc$hyp_constraint)
    G_dec <- recover_full_trans_matrix(theta$g_pars, gc$trans_constraint)

    expect_true(all(G_dec >= 0))
    expect_identical(G_dec[1, 2], 0.5)
    expect_true(abs(G_dec[1, 4] - 5e-6) < 1e-12)
})


# ---- .build_simplify_seeds ----

test_that(".build_simplify_seeds puts the reference first and dedupes", {
    gc <- graph_constraint_free(4)
    ref <- simplify_ref_graph()
    pop_size <- 200L

    seeds <- .build_simplify_seeds(gc, ref, pop_size = pop_size)

    expect_true(is.matrix(seeds))
    expect_identical(typeof(seeds), "double")
    expect_identical(ncol(seeds), length(create_start_params(gc)))
    expect_identical(
        seeds[1, ],
        .encode_graph(gc, ref$hyp_weight, ref$trans_matrix)
    )
    expect_lte(nrow(seeds), pop_size)
    expect_identical(nrow(unique(seeds)), nrow(seeds))
    expect_true(all(is.finite(seeds)))

    # Every neighbour row decodes to a graph with fewer edges than the
    # reference, or is one of the .build_start_matrix() seeds.
    n_ref_edges <- sum(ref$trans_matrix != 0)
    neighbours <- seq_len(nrow(seeds) - 2L)[-1L]
    for (r in neighbours) {
        theta <- split_theta(seeds[r, ], gc$hyp_constraint)
        G_dec <- recover_full_trans_matrix(theta$g_pars, gc$trans_constraint)
        G_dec[G_dec < 1e-5] <- 0
        expect_lt(sum(G_dec != 0), n_ref_edges)
    }
})

test_that(".build_simplify_seeds truncates to pop_size", {
    gc <- graph_constraint_free(4)
    ref <- list(
        hyp_weight = c(0.25, 0.25, 0.25, 0.25),
        trans_matrix = matrix(1 / 3, 4, 4) - diag(1 / 3, 4)
    )
    # A3: a dense reference with many neighbours and a deliberately small
    # popSize must still produce a legal suggestions matrix.
    seeds <- .build_simplify_seeds(gc, ref, pop_size = 4L)
    expect_identical(nrow(seeds), 4L)
    expect_identical(
        seeds[1, ],
        .encode_graph(gc, ref$hyp_weight, ref$trans_matrix)
    )
})

test_that(".build_simplify_seeds appends a population of matching width", {
    gc <- graph_constraint_free(4)
    ref <- simplify_ref_graph()
    n_par <- length(create_start_params(gc))

    without <- .build_simplify_seeds(gc, ref, pop_size = 200L)
    pop <- withr::with_seed(3, matrix(stats::runif(20 * n_par), ncol = n_par))
    with_pop <- .build_simplify_seeds(
        gc,
        ref,
        pop_size = 200L,
        population = pop
    )

    expect_identical(nrow(with_pop), nrow(without) + 20L)
    expect_identical(with_pop[seq_len(nrow(without)), ], without)
})

test_that(".build_simplify_seeds ignores a population of the wrong width", {
    gc <- graph_constraint_free(4)
    ref <- simplify_ref_graph()

    # A8: the stored population came from a different constraint.
    bad_pop <- withr::with_seed(3, matrix(stats::runif(20 * 3), ncol = 3))
    expect_no_error(
        seeds <- .build_simplify_seeds(
            gc,
            ref,
            pop_size = 200L,
            population = bad_pop
        )
    )
    expect_identical(seeds, .build_simplify_seeds(gc, ref, pop_size = 200L))
})

test_that(".build_simplify_seeds skips candidates that keep the edge count", {
    # Row 3 is a single edge of weight 1: dropping it leaves no free
    # recipient, so .redistribute_mass() returns a row summing to zero and the
    # candidate must be skipped rather than seeded.
    gc <- graph_constraint_free(3)
    ref <- list(
        hyp_weight = c(0.5, 0.5, 0),
        trans_matrix = rbind(c(0, 0.5, 0.5), c(0.5, 0, 0.5), c(1, 0, 0))
    )
    expect_true(is_graph_valid(ref$hyp_weight, ref$trans_matrix))

    seeds <- .build_simplify_seeds(gc, ref, pop_size = 200L)

    # Five free non-zero entries, but dropping G[3, 1] leaves row 3 with no
    # free recipient, so only four neighbours are seeded alongside the
    # reference and the two .build_start_matrix() rows.
    expect_identical(nrow(seeds), 1L + 4L + nrow(.build_start_matrix(gc, NULL)))

    n_ref_edges <- sum(ref$trans_matrix != 0)
    for (r in seq_len(nrow(seeds))) {
        theta <- split_theta(seeds[r, ], gc$hyp_constraint)
        G_dec <- recover_full_trans_matrix(theta$g_pars, gc$trans_constraint)
        G_dec[G_dec < 1e-5] <- 0
        # No row collapsed to all zeros. Guarded rows sum to 1 - 5e-6 rather
        # than exactly 1; param_to_solution() renormalises them later.
        expect_true(all(rowSums(G_dec) > 0.99))
        if (r > 1L && r <= 5L) {
            expect_lt(sum(G_dec != 0), n_ref_edges)
        }
    }
})
