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
            "start_graph",
            "alpha",
            "sparsity"
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


# ---- alpha and sparsity on the object ----

test_that("graph_optimise() stores alpha and a NULL sparsity", {
    pv <- withr::with_seed(1, {
        sims <- mvtnorm::rmvnorm(
            1000,
            mean = calc_ncp(c(0.9, 0.8, 0.7)),
            sigma = diag(3)
        )
        stats::pnorm(sims, lower.tail = FALSE)
    })
    ts <- trial_success(r1 + r2 + r3, verbose = "silent")
    ctrl <- multigrain_control() |> control_local(maxeval = 50)

    res <- graph_optimise(
        pvals = pv,
        graph_constraint = graph_constraint_free(3),
        trial_success = ts,
        alpha = 0.01,
        global_search = FALSE,
        control = ctrl,
        verbose = "silent"
    )

    # The previous element names, plus the two new ones at the end.
    expect_named(
        res,
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
            "start_graph",
            "alpha",
            "sparsity"
        )
    )
    expect_identical(res$alpha, 0.01)
    expect_null(res$sparsity)
})

test_that("graph_optimal() validates alpha and defaults both new elements", {
    obj <- graph_optimal(
        hyp_weight = c(0.5, 0.5),
        trans_matrix = matrix(c(0, 1, 1, 0), nrow = 2)
    )
    expect_null(obj$alpha)
    expect_null(obj$sparsity)

    expect_error(
        graph_optimal(
            hyp_weight = c(0.5, 0.5),
            trans_matrix = matrix(c(0, 1, 1, 0), nrow = 2),
            alpha = 1.5
        ),
        class = "rlang_error"
    )
})


# ---- reporting the sparsity element ----

sparsity_example_object <- function() {
    constraints <- graph_constraint_free(3)
    graph_optimal(
        hyp_weight = c(1, 0, 0),
        trans_matrix = rbind(c(0, 1, 0), c(0, 0, 1), c(1, 0, 0)),
        constraints = constraints,
        trial_success = trial_success(r1 + r2 + r3, verbose = "silent"),
        power = list(
            local_power = c(0.9, 0.8, 0.7),
            exp_rejections = 2.4,
            disj_power = 0.95,
            conj_power = 0.6,
            trial_success = 0.8041
        ),
        solution = list(
            opt_source = "simplify:local",
            graph_valid = c(local = TRUE, global = NA)
        ),
        global_search = FALSE,
        alpha = 0.025,
        sparsity = list(
            gain_tolerance = 1e-3,
            reference = list(
                hyp_weight = c(1, 0, 0),
                trans_matrix = rbind(c(0, 0.5, 0.5), c(0.5, 0, 0.5), c(1, 0, 0))
            ),
            gain_reference = 0.8048,
            gain = 0.8041,
            gain_loss = 0.0007,
            gain_loss_fraction = 0.0007 / 0.8048,
            budget = 1e-3 * 0.8048,
            n_edges_reference = 12L,
            n_edges = 7L,
            n_edges_free_reference = 12L,
            n_edges_free = 7L,
            prune_loss = 0.0007,
            source = "local"
        )
    )
}

test_that("summarise_sparsity is silent when there is no sparsity element", {
    expect_snapshot(summarise_sparsity(NULL))
})

test_that("print reports the simplification when sparsity is present", {
    expect_snapshot(print(sparsity_example_object()))
})

test_that("summary reports the simplification when sparsity is present", {
    expect_snapshot(summary(sparsity_example_object()))
})

test_that("print is unchanged when sparsity is NULL", {
    obj <- sparsity_example_object()
    obj$sparsity <- NULL
    out <- utils::capture.output(print(obj))
    expect_false(any(grepl("Simplified from", out, fixed = TRUE)))
})
