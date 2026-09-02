# --------------------------------------------------------------------------
# Integration tests for graph_optimise_n() — these construct real
# trial_success() objects, which triggers Rcpp compilation, so they are kept
# separate from the unit tests. All runs use nsim <= 1e3 with tolerances
# sized accordingly.
# --------------------------------------------------------------------------

skip_on_cran()

disjunctive_ts_3m <- trial_success(r1 || r2 || r3, verbose = FALSE)
conjunctive_ts_3m <- trial_success(r1 && r2 && r3, verbose = FALSE)

# small search budget: integration tests exercise the wiring, not convergence
small_control <- multigrain_control() |>
    control_global(
        popSize = 20L,
        generations_per_block = 5L,
        max_blocks = 3L,
        run = 15L
    )

test_that("graph_optimise_n: end-to-end disjunctive target", {
    result <- graph_optimise_n(
        graph_constraint = graph_constraint_free(3),
        effect_size = rep(0.2, 3),
        trial_success = disjunctive_ts_3m,
        target = 0.8,
        n_range = c(50L, 500L),
        nsim = 1000L,
        num_threads = 1L,
        control = small_control,
        seed = 42,
        verbose = FALSE
    )

    expect_s3_class(result, "multigrain_graph_optimal")
    expect_identical(result$N, result$n_final)
    expect_gte(result$n_final, 50L)
    expect_lte(result$n_final, 500L)
    expect_lte(result$n_final, result$n_init)
    expect_true(all(diff(result$phase2$history$n) <= 0))

    # the returned graph is valid and feasible at n_final (exact on the
    # test's own CRN draw; the MC tolerance covers re-simulation only)
    expect_equal(sum(result$hyp_weight), 1, tolerance = 1e-12)
    expect_equal(
        unname(rowSums(result$trans_matrix)),
        rep(1, 3),
        tolerance = 1e-12
    )
    expect_true(is_graph_valid(result$hyp_weight, result$trans_matrix))
    expect_true(result$solution$graph_valid[["sample_size"]])
    expect_gte(result$power$trial_success, 0.8 - 0.05)

    # the returned object must be slim: graph_optimise_n() runs a
    # Cauchy-mutation search (no GA::ga / nloptr objects), so it records the
    # phase-1/phase-2 search history rather than global_output / local_output,
    # and must never embed the simulated p-value matrix. A round-trip
    # serialisation of a small run stays well under 5 MB (the legacy bloated
    # object embedded the CRN p-values and was ~1 GB).
    expect_null(result$GA_output)
    expect_null(result$nloptr_output)
    expect_false(is.null(result$phase2))
    expect_lt(length(serialize(result, NULL)), 5e6)
})

test_that("graph_optimise_n: marginal-only run, print and summary", {
    floors <- c(0.7, NA, 0.6)

    result <- graph_optimise_n(
        graph_constraint = graph_constraint_free(3),
        effect_size = rep(0.2, 3),
        local_power_target = floors,
        n_range = c(50L, 500L),
        nsim = 1000L,
        control = small_control,
        seed = 43,
        verbose = FALSE
    )

    expect_s3_class(result, "multigrain_graph_optimal")
    expect_null(result$trial_success)
    expect_null(result$power$trial_success)
    expect_true(result$solution$graph_valid[["sample_size"]])
    # every non-NA floor is met on the optimisation's own CRN draw
    idx <- which(!is.na(floors))
    expect_true(all(result$power$local_power[idx] >= floors[idx]))

    # print/summary must handle a NULL trial_success
    expect_no_error(invisible(capture.output(print(result))))
    expect_no_error(invisible(capture.output(summary(result))))
    printed <- paste(capture.output(print(result)), collapse = "\n")
    expect_match(printed, "Selected sample size")
    expect_match(printed, "sample-size optimisation")
    expect_match(printed, "Marginal power floors")
})

test_that("graph_optimise_n: known answer for a conjunctive target", {
    # For equal effects under independence, no valid graph can beat the
    # full-recycling chain for conjunctive power, so the minimal n solves
    # pnorm(0.15 * sqrt(n) - qnorm(0.975))^3 = 0.8. The fixed-sequence seed
    # in the default bank attains this bound, making it a known answer.
    pow <- function(n) pnorm(0.15 * sqrt(n) - qnorm(0.975))^3
    n_ref <- ceiling(
        uniroot(function(n) pow(n) - 0.8, c(50, 2000), tol = 1e-6)$root
    )

    result <- graph_optimise_n(
        graph_constraint = graph_constraint_free(3),
        effect_size = rep(0.15, 3),
        trial_success = conjunctive_ts_3m,
        target = 0.8,
        n_range = c(450L, 600L),
        nsim = 1000L,
        control = small_control,
        seed = 44,
        verbose = FALSE
    )

    # QMC at nsim = 1e3: allow +-10 of the analytic boundary
    expect_lte(abs(result$n_final - n_ref), 10L)
    expect_gte(result$power$trial_success, 0.8 - 0.05)
})

test_that("graph_optimise_n: seeded runs reproduce and match across threads", {
    run_once <- function(num_threads) {
        graph_optimise_n(
            graph_constraint = graph_constraint_free(3),
            effect_size = rep(0.2, 3),
            trial_success = disjunctive_ts_3m,
            target = 0.8,
            n_range = c(50L, 500L),
            nsim = 1000L,
            num_threads = num_threads,
            control = small_control,
            seed = 99,
            verbose = FALSE
        )
    }

    set.seed(123) # ensure an RNG state exists to capture
    rng_before <- .Random.seed
    serial_1 <- run_once(1L)
    expect_identical(rng_before, .Random.seed) # caller RNG state restored

    serial_2 <- run_once(1L)
    expect_identical(serial_1$n_final, serial_2$n_final)
    expect_identical(serial_1$hyp_weight, serial_2$hyp_weight)
    expect_identical(serial_1$trans_matrix, serial_2$trans_matrix)

    parallel_run <- run_once(2L)
    expect_identical(serial_1$n_final, parallel_run$n_final)
    expect_identical(serial_1$hyp_weight, parallel_run$hyp_weight)
    expect_identical(serial_1$trans_matrix, parallel_run$trans_matrix)
})
