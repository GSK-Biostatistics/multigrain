# The fixed-sample consumers must refuse a gain built with
# `trial_success_gsd()`: it inherits from `multigrain_trial_success`, so it
# passes `check_trial_success()`, but its compiled function expects the group
# sequential kernel's decision-time matrix (design record, section 4.10).

gsd_gain_guard <- trial_success_gsd(
    d(t1) + d(t2) + d(t3),
    d = c(1, 0.75),
    verbose = "silent"
)

fixed_gain_guard <- trial_success(r1 + r2 + r3, verbose = "silent")

pvals_guard <- withr::with_seed(42, {
    simulate_pvalues(
        power_nominal = c(0.9, 0.8, 0.7),
        corr_matrix = diag(3),
        nsim = 1000
    )
})

hyp_weight_guard <- c(0.4, 0.3, 0.3)
trans_matrix_guard <- matrix(
    c(
        0, 0.5, 0.5,
        0.5, 0, 0.5,
        0.5, 0.5, 0
    ),
    nrow = 3,
    byrow = TRUE
)


test_that("graph_optimise() refuses a GSD trial success object", {
    expect_error(
        graph_optimise(
            pvals = pvals_guard,
            graph_constraint = graph_constraint_free(3),
            trial_success = gsd_gain_guard,
            verbose = "silent"
        ),
        "was created with `trial_success_gsd()`",
        fixed = TRUE
    )

    expect_error(
        graph_optimise(
            pvals = pvals_guard,
            graph_constraint = graph_constraint_free(3),
            trial_success = gsd_gain_guard,
            verbose = "silent"
        ),
        "graph_optimise_gsd"
    )
})


test_that("calc_power_pvals() refuses a GSD object passed on its own", {
    expect_error(
        calc_power_pvals(
            pvals_guard,
            hyp_weight = hyp_weight_guard,
            trans_matrix = trans_matrix_guard,
            custom_power = gsd_gain_guard
        ),
        "was created with `trial_success_gsd()`",
        fixed = TRUE
    )
})


test_that("calc_power_pvals() refuses a GSD object inside a list", {
    expect_error(
        calc_power_pvals(
            pvals_guard,
            hyp_weight = hyp_weight_guard,
            trans_matrix = trans_matrix_guard,
            custom_power = list(
                average = function(x) sum(x),
                gain = gsd_gain_guard
            )
        ),
        "gain"
    )

    expect_error(
        calc_power_pvals(
            pvals_guard,
            hyp_weight = hyp_weight_guard,
            trans_matrix = trans_matrix_guard,
            custom_power = list(
                average = function(x) sum(x),
                gain = gsd_gain_guard
            )
        ),
        "calc_power_pvals_gsd"
    )

    # an unnamed list element is named by position
    expect_error(
        calc_power_pvals(
            pvals_guard,
            hyp_weight = hyp_weight_guard,
            trans_matrix = trans_matrix_guard,
            custom_power = list(gsd_gain_guard)
        ),
        "Element 1"
    )
})


test_that("a fixed-sample gain still works in both consumers", {
    guard_control <- control_local(multigrain_control(), maxeval = 50)

    power <- calc_power_pvals(
        pvals_guard,
        hyp_weight = hyp_weight_guard,
        trans_matrix = trans_matrix_guard,
        custom_power = list(trial_success = fixed_gain_guard)
    )

    expect_type(power$trial_success, "double")
    expect_length(power$trial_success, 1L)
    expect_length(power$local_power, 3L)

    result <- graph_optimise(
        pvals = pvals_guard,
        graph_constraint = graph_constraint_free(3),
        trial_success = fixed_gain_guard,
        global_search = FALSE,
        control = guard_control,
        verbose = "silent"
    )

    expect_s3_class(result, "multigrain_graph_optimal")
    expect_true(is_graph_valid(result$hyp_weight, result$trans_matrix))
})
