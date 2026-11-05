# to ensure consistent cli output (prevents local vs GHA whitespace differences)
cli::start_app()
on.exit(cli::stop_app(), add = TRUE)

test_that("new_graph_optimal", {
    go <- new_graph_optimal()

    expect_s3_class(go, "multigrain_graph_optimal")

    expect_named(
        go,
        c(
            "hyp_weight",
            "trans_matrix",
            "constraints",
            "trial_success",
            "power",
            "solution",
            "global_search",
            "control",
            "global_output",
            "local_output",
            "start_graph"
        )
    )
})

test_that("graph_optimal", {
    graph_opt <- graph_optimal(
        hyp_weight = c(0.1, 0.2, 0.3, 0.4),
        trans_matrix = matrix(
            c(
                0, 0.1, 0.2, 0.7,
                0.2, 0, 0.4, 0.4,
                0.5, 0.1, 0, 0.4,
                0.8, 0.1, 0.1, 0
            ),
            nrow = 4
        )
    )

    expect_s3_class(graph_opt, "multigrain_graph_optimal")
})

test_that("graph_optimal with solution", {
    expect_no_error(
        new_graph_optimal(
            solution = list(
                opt_source = "local",
                graph_valid = c(
                    local = TRUE,
                    global = FALSE
                )
            )
        )
    )
})

test_that("graph_optimal with global_search", {
    expect_no_error(
        new_graph_optimal(
            global_search = TRUE
        )
    )
})

test_that("graph_optimal_get_control", {
    graph_custom_power <- readRDS(test_path("data", "graph_custom_power.rds"))

    graph_opt_ctrl <- graph_optimal_get_control(graph_custom_power)

    expect_s3_class(graph_opt_ctrl, "multigrain_control")

    expect_identical(
        graph_opt_ctrl$nsim_local,
        100000L
    )

    expect_identical(
        graph_opt_ctrl$nsim_global,
        50000L
    )

    expect_identical(
        graph_opt_ctrl$local_opt,
        list(
            algorithm = "NLOPT_LN_COBYLA",
            xtol_rel = 5e-08,
            xtol_abs = 5e-09,
            maxeval = 5000,
            print_level = 1L
        )
    )

    expect_identical(
        graph_opt_ctrl$global_opt,
        list(
            pcrossover = 0.2,
            pmutation = 0.8,
            maxiter = 1e+05,
            popSize = 200L,
            run = 200,
            monitor = TRUE,
            optimArgs = list(
                method = "Nelder-Mead",
                poptim = 0.2,
                pressel = 0.6
            )
        )
    )
})

test_that("summarise helpers", {
    expect_snapshot(
        summarise_solution_source(NULL)
    )

    expect_snapshot(
        summarise_power_object(NULL)
    )

    expect_snapshot({
        summarise_solution_source(
            list(
                opt_source = "GA"
            )
        )
    })

    expect_snapshot({
        summarise_solution_source(
            list(
                opt_source = "GA_minN"
            )
        )
    })

    expect_snapshot({
        summarise_solution_source(
            list(
                opt_source = "foo"
            )
        )
    })
})

test_that("graph_optimal print & summary methods", {
    disjunctive_3m_power <- trial_success(r1 || r2 || r3, verbose = "silent")

    obj <- graph_optimal(
        hyp_weight = c(H1 = 0.5, H2 = 0.3, H3 = 0.2),
        trans_matrix = diag(3),
        trial_success = disjunctive_3m_power,
        power = list(trial_success = 0.85),
        solution = list(
            opt_source = "local",
            graph_valid = c(
                local = TRUE,
                global = TRUE
            )
        ),
        global_search = FALSE
    )
    expect_snapshot(print(obj))
    expect_null(print.multigrain_graph_optimal(NULL))

    expect_snapshot(summary(obj))
    expect_null(summary.multigrain_graph_optimal(NULL))

    graph_custom_power <- readRDS(test_path("data", "graph_custom_power.rds"))

    expect_snapshot(print(graph_custom_power))
    expect_snapshot(summary(graph_custom_power))

    expect_snapshot(print(graph_optimal_example))
    expect_snapshot(summary(graph_optimal_example))
})

test_that("graph_optimal print & summary methods, group sequential", {
    pvals <- withr::with_seed(31, {
        simulate_pvalues_gsd(c(0.9, 0.8), info_frac = c(0.5, 1), nsim = 500)
    })
    pvals <- transform_pvalues_gsd(
        pvals,
        spending = gsDesign::sfLDOF,
        grid_size = 64L
    )
    gain <- trial_success_gsd(
        d(t1) + d(t2),
        d = c(1, 0.75),
        verbose = "silent"
    )
    rule <- stopping_rule_gsd(r1)
    hyp_weight <- c(0.5, 0.5)
    trans_matrix <- rbind(c(0, 1), c(1, 0))

    power <- calc_power_pvals_gsd(
        pvals,
        hyp_weight,
        trans_matrix,
        custom_power = list(trial_success = gain),
        stopping_rule = rule
    )
    obj <- graph_optimal(
        hyp_weight = hyp_weight,
        trans_matrix = trans_matrix,
        trial_success = gain,
        power = power,
        stopping_rule = rule
    )
    expect_identical(obj$stopping_rule, rule)
    expect_snapshot(print(obj))
    expect_snapshot(summary(obj))

    # without a rule every trial stops at the last analysis unless it
    # rejected every hypothesis at the interim
    power_no_rule <- calc_power_pvals_gsd(
        pvals,
        hyp_weight,
        trans_matrix,
        custom_power = list(trial_success = gain)
    )
    obj_no_rule <- graph_optimal(
        hyp_weight = hyp_weight,
        trans_matrix = trans_matrix,
        trial_success = gain,
        power = power_no_rule
    )
    expect_null(obj_no_rule$stopping_rule)
    expect_snapshot(summary(obj_no_rule))
})

test_that("summarise_power_object shows the group sequential fields", {
    look_names <- c("analysis 1", "analysis 2")
    power_object <- list(
        local_power_by_analysis = matrix(
            c(0.12345, 0.2, 0.8, 0.65432),
            nrow = 2,
            dimnames = list(c("H1", "H2"), look_names)
        ),
        time_distribution = matrix(
            c(0.2, 0.34568, 0.12345, 0.2, 0.67655, 0.45432),
            nrow = 2,
            dimnames = list(c("H1", "H2"), c("never", look_names))
        ),
        stop_distribution = c("analysis 1" = 0.25, "analysis 2" = 0.75),
        mean_stop_look = 1.75
    )
    expect_snapshot(summarise_power_object(power_object))
})

test_that("is_graph_optimal", {
    graph_custom_power <- readRDS(test_path("data", "graph_custom_power.rds"))

    expect_true(is_graph_optimal(graph_custom_power))
    expect_false(is_graph_optimal("foo"))
})

test_that("check_graph_optimal", {
    graph_custom_power <- readRDS(test_path("data", "graph_custom_power.rds"))

    expect_no_error(
        check_graph_optimal(
            graph_custom_power
        )
    )

    expect_error(
        check_graph_optimal(2),
        "`2` must be a multigrain graph optimal object, not the number 2."
    )
})
