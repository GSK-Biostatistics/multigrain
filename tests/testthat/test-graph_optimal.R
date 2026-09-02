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

test_that("graph_optimal complains", {
    expect_error(
        graph_optimal(hyp_weight = "foo"),
        "`hyp_weight` must be a double"
    )

    hyp_w <- c(0.1, 0.2, NA, NA, NA)

    expect_error(
        graph_optimal(
            hyp_weight = hyp_w,
            trans_matrix = "foo"
        ),
        "`trans_matrix` must be a double matrix"
    )

    trans_m <- matrix(rep_len(c(0.1, NA), length.out = 25), nrow = 5, ncol = 5)
    diag(trans_m) <- 0

    expect_error(
        graph_optimal(
            hyp_weight = hyp_w,
            trans_matrix = trans_m,
            constraints = "foo"
        ),
        "`constraints` must be a multigrain graph constraint object"
    )

    expect_error(
        graph_optimal(
            hyp_weight = hyp_w,
            trans_matrix = trans_m,
            trial_success = "foo"
        ),
        "`trial_success` must be a multigrain trial success object"
    )

    expect_error(
        graph_optimal(
            hyp_weight = hyp_w,
            trans_matrix = trans_m,
            control = "foo"
        ),
        "`control` must be a multigrain control object"
    )

    expect_error(
        graph_optimal(
            hyp_weight = hyp_w,
            trans_matrix = trans_m,
            global_output = "foo"
        ),
        "`global_output` must be a GA object"
    )

    expect_error(
        graph_optimal(
            hyp_weight = hyp_w,
            trans_matrix = trans_m,
            local_output = "foo"
        ),
        "`local_output` must be a nloptr object"
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

    expect_snapshot({
        graph_optimal_get_control(graph_custom_power)
    })

    expect_snapshot({
        graph_optimal_get_control(graph_optimal_example)
    })
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
    disjunctive_3m_power <- trial_success(r1 || r2 || r3, verbose = FALSE)

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
