# Trial success functions to be used for tests
disjunctive_3m_power <- trial_success(r1 || r2 || r3, verbose = FALSE)
conjunctive_4m_power <- trial_success(r1 && r2 && r3 && r4, verbose = FALSE)
avg_6m_power <- trial_success(r1 + r2 + r3 + r4 + r5 + r6, verbose = FALSE)

# simulate p-values for power calculations for rest of tests
power_vec <- c(0.9, 0.8, 0.6, 0.85, 0.85, 0.5)
ncp_vec <- rep(stats::qnorm(1 - 0.025), 6) - stats::qnorm(1 - power_vec)
corr_test <- matrix(
    0.2,
    nrow = 6,
    ncol = 6,
    byrow = TRUE
)
diag(corr_test) <- 1

pvals <- withr::with_seed(5, {
    sims <- mvtnorm::rmvnorm(2^20, mean = ncp_vec, sigma = corr_test)
    pvals <- stats::pnorm(sims, lower.tail = FALSE)
})

expect_ga_result <- function(x, m) {
    testthat::expect_type(x, "list")
    testthat::expect_length(x$ga_hyp_weight, m)

    testthat::expect_identical(dim(x$ga_trans_matrix), c(m, m))

    testthat::expect_type(x$is_graph_valid, "logical")
    testthat::expect_s4_class(x$ga_output, "ga")
}

expect_local_result <- function(x, m) {
    testthat::expect_type(x, "list")
    testthat::expect_length(x$local_hyp_weight, m)

    testthat::expect_identical(dim(x$local_trans_matrix), c(m, m))

    testthat::expect_type(x$is_graph_valid, "logical")
    testthat::expect_type(x$local_output, "list")
}


### OPTIMISATION
# Create setup for optimisation #
m <- 4
alpha <- 0.025
power_vector <- c(0.93, 0.91, 0.90, 0.85)

corr_mat <- matrix(
    0.2,
    nrow = m,
    ncol = m,
    byrow = TRUE
)
diag(corr_mat) <- 1

sims <- mvtnorm::rmvnorm(
    5e4,
    sigma = as.matrix(corr_mat),
    mean = calc_ncp(power_vector)
)
pvals_4m <- stats::pnorm(sims, lower.tail = FALSE)
pvals_rand <- pvals_4m[sample.int(1e4), ]

no_constr <- graph_constraint_free(m)
constraint_4m <- graph_constraint(hyp_constraint = c(1, 0, 0, 0))


## GA::ga optimisation
test_that(".graph_optimise_ga returns valid object", {
    ctrl <- multigrain_control() |>
        control_global(
            maxiter = 10,
            popSize = 20,
            run = 5
        ) |>
        control_prepare(
            pvals = pvals_rand
        )

    set.seed(1)

    expect_snapshot(
        {
            ga_res <- .graph_optimise_ga(
                pvals = pvals_rand,
                graph_constraint = no_constr,
                trial_success = conjunctive_4m_power,
                alpha = alpha,
                nsim = ctrl$nsim_global,
                global_opts = ctrl$global_opt
            )
        },
        transform = fix_duration
    )

    expect_ga_result(ga_res, m = 4L)

    expect_gt(ga_res$ga_output@iter, 0) # at l

    if (isTRUE(ga_res$is_graph_valid)) {
        expect_gt(ga_res$ga_trial_success, 0)
    }
})

## nloptr::nloptr optimisation
test_that(".graph_optimise_local returns valid object", {
    ctrl <- multigrain_control() |>
        control_local(maxeval = 200, print_level = 0) |>
        control_prepare(
            pvals = pvals_rand
        )

    set.seed(3)

    expect_snapshot(
        {
            loc_res <- .graph_optimise_local(
                pvals = pvals_rand,
                graph_constraint = no_constr,
                trial_success = conjunctive_4m_power,
                alpha = alpha,
                local_opts = ctrl$local_opt
            )
        },
        transform = fix_duration
    )

    expect_local_result(loc_res, m = 4L)

    expect_gt(nchar(loc_res$local_output$message), 0)

    if (isTRUE(loc_res$is_graph_valid)) {
        expect_gt(loc_res$local_trial_success, 0)
    }
})


test_that("x0 passed successfully from .graph_optimise_ga to local", {
    ctrl <- multigrain_control() |>
        control_global(
            maxiter = 20,
            popSize = 10,
            run = 5
        ) |>
        control_prepare(
            pvals = pvals_rand
        )

    set.seed(3)
    expect_snapshot(
        {
            ga_x0 <- .graph_optimise_ga(
                pvals = pvals_rand,
                graph_constraint = no_constr,
                trial_success = conjunctive_4m_power,
                alpha = alpha,
                global_opts = ctrl$global_opt,
                nsim = ctrl$nsim_global
            )$ga_output@solution[1, ] |>
                as.vector()
        },
        transform = fix_duration
    )

    ctrl <- multigrain_control() |>
        control_local(maxeval = 100, print_level = 0) |>
        control_prepare(
            pvals = pvals_rand
        )

    expect_snapshot(
        {
            loc_res <- .graph_optimise_local(
                pvals = pvals_rand,
                graph_constraint = no_constr,
                trial_success = conjunctive_4m_power,
                alpha = alpha,
                local_opts = ctrl$local_opt,
                x0 = ga_x0
            )
        },
        transform = fix_duration
    )

    # Should not change x0 length / dims
    expect_length(loc_res$local_output$solution, length(ga_x0))
})

test_that("fix_duration", {
    expect_identical(
        fix_duration("v Evaluating trial success of pruned graph [34ms]"),
        "v Evaluating trial success of pruned graph [10ms]"
    )

    expect_identical(
        fix_duration(
            "v Evaluating trial success of locally optimised graph [1m 18.9s]"
        ),
        "v Evaluating trial success of locally optimised graph [10ms]"
    )

    expect_identical(
        fix_duration(
            "
            i Running local optimization
            v Running local optimization [1.8s]

            i Evaluating trial success of locally optimised graph
            v Evaluating trial success of locally optimised graph [389ms]

            i Evaluating trial success of pruned graph
            v Evaluating trial success of pruned graph [34ms]
            "
        ),
        "
            i Running local optimization
            v Running local optimization [10ms]

            i Evaluating trial success of locally optimised graph
            v Evaluating trial success of locally optimised graph [10ms]

            i Evaluating trial success of pruned graph
            v Evaluating trial success of pruned graph [10ms]
            "
    )
})

test_that("Optimise 4m conjunctive power: local search only", {
    ctrl <- multigrain_control() |>
        control_global(run = 7) |>
        control_nsim_local(2e4)

    set.seed(10)
    expect_snapshot(
        {
            graph_optimise_4m_result <- graph_optimise(
                pvals = pvals_4m,
                alpha = alpha,
                graph_constraint = no_constr,
                trial_success = conjunctive_4m_power,
                control = ctrl,
                global_search = FALSE
            )
        },
        transform = fix_duration
    )

    expect_s3_class(graph_optimise_4m_result, "multigrain_graph_optimal")
    expect_false(graph_optimise_4m_result$global_search)
    expect_length(
        graph_optimise_4m_result$hyp_weight,
        m
    )
})


test_that("Optimise 4m conjunctive power: include global search", {
    ctrl <- multigrain_control() |>
        control_global(run = 7) |>
        control_nsim_global(2e4) |>
        control_nsim_local(2e4)

    set.seed(20)
    expect_snapshot(
        {
            graph_optimise_4m_result <- graph_optimise(
                pvals = pvals_4m,
                alpha = alpha,
                graph_constraint = no_constr,
                trial_success = conjunctive_4m_power,
                control = ctrl
            )
        },
        transform = fix_duration
    )

    expect_s3_class(graph_optimise_4m_result, "multigrain_graph_optimal")
    expect_true(graph_optimise_4m_result$global_search)
    expect_length(
        graph_optimise_4m_result$hyp_weight,
        m
    )
})

test_that("Optimise 4m conjunctive power: include global search & verbose", {
    ctrl <- multigrain_control() |>
        control_global(run = 7) |>
        control_nsim_global(2e4) |>
        control_nsim_local(2e4)

    set.seed(20)
    graph_optimise_4m_result <- graph_optimise(
        pvals = pvals_4m,
        alpha = alpha,
        graph_constraint = no_constr,
        trial_success = conjunctive_4m_power,
        control = ctrl,
        verbose = FALSE
    )

    expect_s3_class(graph_optimise_4m_result, "multigrain_graph_optimal")
    expect_true(graph_optimise_4m_result$global_search)
    expect_length(
        graph_optimise_4m_result$hyp_weight,
        m
    )

    # verbose must be logical
    expect_snapshot(error = TRUE, {
        graph_optimise(
            pvals = pvals_4m,
            alpha = alpha,
            graph_constraint = no_constr,
            trial_success = conjunctive_4m_power,
            control = ctrl,
            verbose = "foo"
        )
    })
})

test_that("Optimise 2m: E2E check", {
    set.seed(1)
    ts <- trial_success(r1 || r2, verbose = FALSE)

    pwr_vector <- c(0.975, 0.50)
    corr_matrix <- matrix(c(1, 0.2, 0.2, 1), nrow = 2)

    pvals_2 <- simulate_pvalues(
        power_nominal = pwr_vector,
        alpha = 0.025,
        corr_matrix = corr_matrix,
        nsim = 1e+04
    )

    gc <- graph_constraint(c(NA, NA), matrix(c(0, 1, 1, 0), nrow = 2))

    ctrl <- multigrain_control() |>
        control_global(run = 20) |>
        control_nsim_global(1e3)

    set.seed(123)

    expect_snapshot(
        {
            graph_2m <- graph_optimize(
                pvals_2,
                graph_constraint = gc,
                trial_success = ts,
                alpha = 0.025,
                control = ctrl
            )
        },
        transform = fix_duration
    )

    expect_length(graph_2m$hyp_weight, 2L)
    expect_s3_class(graph_2m, "multigrain_graph_optimal")
    expect_true(graph_2m$global_search)
    expect_identical(unname(graph_2m$trans_matrix), unname(gc$trans_constraint))
})


test_that("Test graph_optimise argument checks", {
    ctrl <- multigrain_control() |>
        control_global(run = 7) |>
        control_nsim_global(2e4) |>
        control_nsim_local(2e4)

    set.seed(20)

    args <- list(
        pvals = pvals_rand,
        alpha = alpha,
        graph_constraint = no_constr,
        trial_success = conjunctive_4m_power,
        control = ctrl
    )

    alpha_args <- args
    alpha_args$alpha <- c(0.24, 0.66)
    expect_error(
        do.call(graph_optimise, alpha_args),
        "`alpha` must be a number, not a double vector."
    )

    alpha_args$alpha <- function(x) x + 2
    expect_error(
        do.call(graph_optimise, alpha_args),
        "`alpha` must be a number, not a function."
    )
})


## Add tests to check warnings for matching dimensions of pvals/graph_constraint
test_that("Prevent pvals & trial_success dimension mismatch", {
    set.seed(19)

    alpha <- 0.025
    no_constraint_3m <- graph_constraint_free(3)

    args <- list(
        pvals = pvals[, 1:3],
        alpha = alpha,
        graph_constraint = no_constraint_3m,
        trial_success = avg_6m_power
    ) # custom power function with 6 hypotheses

    expect_error(
        do.call(graph_optimise, args),
        "`multigrain_trial_success` dimensions (m) must match the number of columns in `pvals`", # nolint
        fixed = TRUE
    )
})

test_that(".graph_optimise_ga returns valid object", {
    set.seed(1)

    m <- 4
    pvals_4m <- pvals[1:1e4, 1:4]
    alpha <- 0.025
    no_constr <- graph_constraint_free(m)

    ctrl <- multigrain_control() |>
        control_global(
            maxiter = 10,
            popSize = 20,
            run = 5
        ) |>
        control_prepare(
            pvals = pvals_4m
        )

    expect_snapshot(
        {
            ga_res <- .graph_optimise_ga(
                pvals = pvals_4m,
                graph_constraint = no_constr,
                trial_success = conjunctive_4m_power,
                alpha = alpha,
                global_opts = ctrl$global_opt,
                num_threads = cran_cores(),
                nsim = ctrl$nsim_global
            )
        },
        transform = fix_duration
    )

    expect_ga_result(ga_res, m = 4L)
    expect_gt(ga_res$ga_output@iter, 0)

    if (isTRUE(ga_res$is_graph_valid)) {
        expect_gt(ga_res$ga_trial_success, 0)
    }
})

test_that(".graph_optimise_local returns valid object", {
    set.seed(3)

    m <- 4
    pvals_4m <- pvals[1:1e4, 1:4]
    alpha <- 0.025
    no_constr <- graph_constraint_free(m)

    ctrl <- multigrain_control() |>
        control_local(
            max_eval = 200,
            print_level = 0
        ) |>
        control_prepare(
            pvals = pvals_4m
        )

    expect_snapshot(
        {
            loc_res <- .graph_optimise_local(
                pvals = pvals_4m,
                graph_constraint = no_constr,
                trial_success = conjunctive_4m_power,
                alpha = alpha,
                local_opts = ctrl$local_opt,
                num_threads = cran_cores()
            )
        },
        transform = fix_duration
    )

    expect_local_result(loc_res, m = 4L)
    expect_gt(nchar(loc_res$local_output$message), 0)

    if (isTRUE(loc_res$is_graph_valid)) {
        expect_gt(loc_res$local_trial_success, 0)
    }
})

test_that("x0 passed successfully from .graph_optimise_ga to local", {
    set.seed(3)

    m <- 4
    pvals_4m <- pvals[1:1e4, 1:4]
    alpha <- 0.025
    no_constr <- graph_constraint_free(m)

    ctrl <- multigrain_control() |>
        control_global(
            maxiter = 10,
            popSize = 20,
            run = 5
        ) |>
        control_local(
            max_eval = 100,
            print_level = 0
        ) |>
        control_prepare(
            pvals = pvals_4m
        )

    expect_snapshot(
        {
            ga_x0 <- .graph_optimise_ga(
                pvals = pvals_4m,
                graph_constraint = no_constr,
                trial_success = conjunctive_4m_power,
                alpha = alpha,
                global_opts = ctrl$global_opt,
                num_threads = 1L,
                nsim = ctrl$nsim_global
            )$ga_output@solution[1, ] |>
                as.vector()
        },
        transform = fix_duration
    )

    expect_snapshot(
        {
            loc_res <- .graph_optimise_local(
                pvals = pvals_4m,
                graph_constraint = no_constr,
                trial_success = conjunctive_4m_power,
                alpha = alpha,
                local_opts = ctrl$local_opt,
                num_threads = 1L,
                x0 = ga_x0
            )
        },
        transform = fix_duration
    )

    expect_length(loc_res$local_output$solution, length(ga_x0))
})

test_that("graph_optimise with num_threads: local search only", {
    set.seed(10)

    m <- 4
    alpha <- 0.025
    no_constr <- graph_constraint_free(m)

    ctrl <- multigrain_control() |>
        control_local(print_level = 0) |>
        control_nsim_local(2e4)

    expect_snapshot(
        {
            graph_optimise_4m_result <- graph_optimise(
                pvals = pvals[, 1:4],
                alpha = alpha,
                graph_constraint = no_constr,
                trial_success = conjunctive_4m_power,
                num_threads = cran_cores(),
                control = ctrl,
                global_search = FALSE
            )
        },
        transform = fix_duration
    )

    expect_s3_class(graph_optimise_4m_result, "multigrain_graph_optimal")
    expect_length(graph_optimise_4m_result$hyp_weight, m)
})

test_that("graph_optimise with num_threads: include global search", {
    set.seed(20)

    m <- 4
    alpha <- 0.025
    no_constr <- graph_constraint_free(m)

    ctrl <- multigrain_control() |>
        control_local(print_level = 0) |>
        control_global(run = 7) |>
        control_nsim_global(2e4) |>
        control_nsim_local(2e4)

    expect_snapshot(
        {
            graph_optimise_4m_result <- graph_optimise(
                pvals = pvals[, 1:4],
                alpha = alpha,
                graph_constraint = no_constr,
                trial_success = conjunctive_4m_power,
                num_threads = cran_cores(),
                control = ctrl
            )
        },
        transform = fix_duration
    )

    expect_s3_class(graph_optimise_4m_result, "multigrain_graph_optimal")
    expect_length(graph_optimise_4m_result$hyp_weight, m)
})


test_that("graph_optimise parallel vs serial identical results (1 thread)", {
    m <- 3
    alpha <- 0.025
    pvals_3m <- pvals[1:5e3, 1:3]
    no_constr <- graph_constraint_free(m)

    # Serial optimisation
    ctrl <- multigrain_control() |>
        control_local(max_eval = 100, print_level = 0) |>
        control_nsim_local(nrow(pvals_3m))

    set.seed(42)
    expect_snapshot(
        {
            result_serial <- graph_optimise(
                pvals = pvals_3m,
                alpha = alpha,
                graph_constraint = no_constr,
                trial_success = disjunctive_3m_power,
                control = ctrl,
                global_search = FALSE
            )
        },
        transform = fix_duration
    )

    # Parallel optimization with 1 thread (should be identical)
    set.seed(42) # Reset seed

    expect_snapshot(
        {
            result_parallel_1t <- graph_optimise(
                pvals = pvals_3m,
                alpha = alpha,
                graph_constraint = no_constr,
                trial_success = disjunctive_3m_power,
                num_threads = 1L,
                control = ctrl,
                global_search = FALSE
            )
        },
        transform = fix_duration
    )

    # Results should be identical
    expect_identical(
        result_serial$hyp_weight,
        result_parallel_1t$hyp_weight
    )
    expect_identical(
        result_serial$trans_matrix,
        result_parallel_1t$trans_matrix
    )
    expect_identical(
        result_serial$power$trial_success,
        result_parallel_1t$power$trial_success
    )
})

test_that("graph_optimise: 2m E2E check with parallelization", {
    set.seed(1)

    ts <- trial_success(r1 || r2, verbose = FALSE)
    pwr_vector <- c(0.975, 0.50)
    corr_matrix <- matrix(c(1, 0.2, 0.2, 1), nrow = 2)

    pvals_2m <- simulate_pvalues(
        power_nominal = pwr_vector,
        alpha = 0.025,
        corr_matrix = corr_matrix,
        nsim = 1e+04
    )

    gc <- graph_constraint(
        c(NA, NA),
        matrix(
            c(0, 1, 1, 0),
            nrow = 2
        )
    )

    ctrl <- multigrain_control() |>
        control_global(run = 20) |>
        control_nsim_global(1e3)

    set.seed(123)

    expect_snapshot(
        {
            graph_2m <- graph_optimise(
                pvals_2m,
                graph_constraint = gc,
                trial_success = ts,
                alpha = 0.025,
                num_threads = 2L,
                control = ctrl
            )
        },
        transform = fix_duration
    )

    expect_length(graph_2m$hyp_weight, 2L)
    expect_s3_class(graph_2m, "multigrain_graph_optimal")
    expect_identical(unname(graph_2m$trans_matrix), unname(gc$trans_constraint))
})


## ------ NAN guards in create_obj_func ------ ##
test_that("create_obj_func returns penalty when x contains NaN", {
    pvals_4m <- pvals[, 1:4]
    alpha <- 0.025

    no_constraint_4m <- graph_constraint_free(4)

    obj_func <- create_obj_func(
        4,
        power_criterion = conjunctive_4m_power$func,
        hyp_constraint = no_constraint_4m$hyp_constraint,
        trans_constraint = no_constraint_4m$trans_constraint,
        alpha = alpha,
        pvals = pvals_4m,
        num_threads = 1L
    )

    x0 <- create_start_params(no_constraint_4m)
    x_nan <- x0
    x_nan[1] <- NaN

    result <- obj_func(
        x = x_nan
    )

    expect_identical(result, -1e6)
})


test_that("graph_optimise stores settings correctly", {
    set.seed(50)

    m <- 3
    alpha <- 0.025
    pvals_3m <- pvals[1:5e3, 1:3]
    no_constr <- graph_constraint_free(m)
    ctrl <- multigrain_control() |>
        control_local(
            max_eval = 50,
            print_level = 0
        ) |>
        control_nsim_local(2e3)

    expect_snapshot(
        {
            result <- graph_optimise(
                pvals = pvals_3m,
                alpha = alpha,
                graph_constraint = no_constr,
                trial_success = disjunctive_3m_power,
                num_threads = cran_cores(),
                control = ctrl,
                global_search = FALSE
            )
        },
        transform = fix_duration
    )

    expect_s3_class(result, "multigrain_graph_optimal")
})


# Test for pre-pruning bug #
test_that("graph_optimise $power reflects pruned graph, not pre-pruned", {
    ctrl <- multigrain_control() |>
        control_global(
            maxiter = 30,
            run = 20
        ) |>
        control_nsim_global(1e4) |>
        control_nsim_local(2e4)

    set.seed(20)

    expect_snapshot(
        {
            result <- graph_optimise(
                pvals = pvals[, 1:6],
                alpha = 0.025,
                graph_constraint = graph_constraint_free(6),
                trial_success = avg_6m_power,
                control = ctrl
            )
        },
        transform = fix_duration
    )

    # Re-evaluate on the returned (pruned) weights
    check_power <- calc_power_pvals(
        pvals = pvals[, 1:6],
        alpha = 0.025,
        hyp_weight = result$hyp_weight,
        trans_matrix = result$trans_matrix,
        custom_power = list(
            trial_success = avg_6m_power
        )
    )

    expect_identical(
        result$power$trial_success,
        check_power$trial_success
    )
})


#### Issue with pruning graph

power_vec <- c(0.9, 0.8, 0.6, 0.85, 0.85, 0.5)
ncp_vec <- rep(stats::qnorm(1 - 0.025), 6) - stats::qnorm(1 - power_vec)
corr_test <- matrix(
    0.2,
    nrow = 6,
    ncol = 6,
    byrow = TRUE
)
diag(corr_test) <- 1

pvals_6m <- withr::with_seed(5, {
    sims <- mvtnorm::rmvnorm(2^20, mean = ncp_vec, sigma = corr_test)
    pvals_6m <- stats::pnorm(sims, lower.tail = FALSE)
})

test_that("deprecation message for power_nsim_local", {
    expect_snapshot(error = TRUE, {
        graph_optimise_4m_result <- graph_optimise(
            pvals = pvals,
            alpha = alpha,
            graph_constraint = no_constr,
            trial_success = conjunctive_4m_power,
            power_nsim_local = 200
        )
    })
})

test_that("deprecation message for power_nsim_global", {
    expect_snapshot(error = TRUE, {
        graph_optimise_4m_result <- graph_optimise(
            pvals = pvals,
            alpha = alpha,
            graph_constraint = no_constr,
            trial_success = conjunctive_4m_power,
            power_nsim_global = 200
        )
    })
})
