# nolint start: object_usage_linter. Helpers live in helper-hsd_opt.R.
# Phase P3 of dev/hsd_opt_build_plan.md: the driver of the joint graph and
# spending-parameter search (Stage G, Stage L, the profile and the dispatch).
# Fixtures are in helper-hsd_opt.R; every run uses the tiny optimiser
# settings of test-optimisation_gsd.R.

opt_hsd_gain_3 <- trial_success_gsd(
    d(t1) + d(t2) + d(t3),
    d = c(1, 0.75),
    verbose = "silent"
)

opt_hsd_gain_4 <- trial_success_gsd(
    d(t1) + d(t2) + d(t3) + d(t4),
    d = c(1, 0.9),
    verbose = "silent"
)

opt_hsd_rule <- stopping_rule_gsd(r1)

opt_hsd_run <- function(
    spec = hsd_opt_spec(),
    gain = opt_hsd_gain_3,
    seed = 202,
    ...
) {
    withr::with_seed(seed, {
        graph_optimise_gsd(
            pvals = spec,
            graph_constraint = graph_constraint_free(spec$m),
            trial_success = gain,
            control = hsd_opt_control(),
            verbose = "silent",
            ...
        )
    })
}

opt_hsd_result <- function() {
    hsd_opt_cached("opt_result_d1", opt_hsd_run())
}

opt_hsd_result_d2 <- function() {
    hsd_opt_cached(
        "opt_result_d2",
        opt_hsd_run(
            spec = hsd_opt_spec_d2(),
            gain = opt_hsd_gain_4,
            seed = 303,
            stopping_rule = opt_hsd_rule
        )
    )
}

# The prepared local options of the tiny control at the d = 1 fixture.
opt_hsd_local_opts <- function(spec = hsd_opt_spec()) {
    control_prepare_dims(
        hsd_opt_control(),
        nsim = spec$nsim,
        m = spec$m
    )$local_opt
}


test_that("(a) step 1 is the fixed local driver's nloptr call", {
    spec <- hsd_opt_spec()
    gc <- graph_constraint_free(3)
    local_opts <- opt_hsd_local_opts(spec)
    x0 <- .hsd_encode_graph(
        gc,
        c(0.5, 0.3, 0.2),
        rbind(c(0, 0.5, 0.5), c(0.7, 0, 0.3), c(0.4, 0.6, 0))
    )

    for (node in c(-4, -3, -2)) {
        refined <- withr::with_seed(11, {
            idx <- sample.int(spec$nsim, size = 1500L, replace = FALSE)
            ctx <- .hsd_local_context(
                spec,
                graph_constraint = gc,
                trial_success = opt_hsd_gain_3,
                alpha = spec$alpha,
                local_opts = local_opts,
                idx = idx
            )
            kmat <- .gsd_kernel_matrix(.hsd_node_array(ctx, node))
            obj <- .hsd_objective(ctx, kmat)
            list(
                ours = .hsd_refine_node(x0, obj = obj, kmat = kmat, ctx = ctx),
                by_hand = nloptr::nloptr(
                    x0 = x0,
                    eval_f = function(x) -obj(x),
                    lb = rep(0, length(x0)),
                    ub = rep(1, length(x0)),
                    opts = local_opts
                )
            )
        })

        # the by-hand call on the identical closure, start and options
        expect_identical(
            refined$ours$nloptr_output$solution,
            refined$by_hand$solution
        )
        expect_identical(
            refined$ours$nloptr_output$objective,
            refined$by_hand$objective
        )
        expect_gt(refined$ours$n_eval, 0L)

        # and the fixed local driver on the resolved array, same seed: it
        # draws the same rows, so it is the same optimisation
        fixed <- withr::with_seed(11, {
            .graph_optimise_local_gsd(
                pvals = resolve_pvalues_gsd(spec, node),
                graph_constraint = gc,
                trial_success = opt_hsd_gain_3,
                local_opts = local_opts,
                alpha = spec$alpha,
                nsim = 1500L,
                x0 = x0,
                verbose = "silent"
            )
        })
        expect_identical(
            refined$ours$nloptr_output$solution,
            fixed$local_output$solution
        )
        expect_identical(
            refined$ours$nloptr_output$objective,
            fixed$local_output$objective
        )
        expect_identical(refined$ours$hyp_weight, fixed$local_hyp_weight)
        expect_identical(refined$ours$trans_matrix, fixed$local_trans_matrix)
    }
})

test_that("(a) a warm-start encoding decodes to its graph", {
    gc <- graph_constraint_free(3)
    graphs <- list(
        list(
            w = c(1 / 3, 1 / 3, 1 / 3),
            G = rbind(c(0, 0.5, 0.5), c(0.5, 0, 0.5), c(0.5, 0.5, 0))
        ),
        list(
            w = c(0.6, 0.4, 0),
            G = rbind(c(0, 1, 0), c(0, 0, 1), c(1, 0, 0))
        ),
        list(
            w = c(0.2, 0.25, 0.55),
            G = rbind(c(0, 0.3, 0.7), c(0.9, 0, 0.1), c(0.45, 0.55, 0))
        )
    )

    for (g in graphs) {
        x <- .hsd_encode_graph(gc, g$w, g$G)
        expect_type(x, "double")
        expect_null(attributes(x))
        sol <- param_to_solution(x, gc)
        # a decode is arithmetic: equal to floating point, not bitwise
        # nolint start: expect_identical_linter
        expect_equal(unname(sol$hyp_weight), g$w)
        expect_equal(unname(sol$trans_matrix), g$G)
        # nolint end
    }
})

test_that("(b) the reported power is calc_power_pvals_gsd() at the optimum", {
    res <- opt_hsd_result()

    expect_identical(
        res$power,
        calc_power_pvals_gsd(
            resolve_pvalues_gsd(hsd_opt_spec(), res$spending_param),
            res$hyp_weight,
            res$trans_matrix,
            custom_power = list(trial_success = opt_hsd_gain_3)
        )
    )

    res_d2 <- opt_hsd_result_d2()
    expect_identical(
        res_d2$power,
        calc_power_pvals_gsd(
            resolve_pvalues_gsd(hsd_opt_spec_d2(), res_d2$spending_param),
            res_d2$hyp_weight,
            res_d2$trans_matrix,
            custom_power = list(trial_success = opt_hsd_gain_4),
            stopping_rule = opt_hsd_rule
        )
    )
})

test_that("(c) the same seed gives the same result", {
    first <- opt_hsd_result()
    second <- opt_hsd_run()

    keep <- setdiff(names(first), "trial_success")
    expect_named(second, names(first))
    expect_identical(second[keep], first[keep])
})

test_that("(d) the result fields at m = 4, d = 2", {
    spec <- hsd_opt_spec_d2()
    res <- opt_hsd_result_d2()
    lattice <- spec$param$lattice

    expect_s3_class(res, "multigrain_graph_optimal")
    expect_length(res$spending_param, 2L)
    expect_named(res$spending_param, c("gamma1", "b"))
    for (r in 1:2) {
        expect_true(res$spending_param[[r]] %in% lattice[[r]])
        expect_gte(res$spending_param[[r]], spec$param$bounds[r, 1L])
        expect_lte(res$spending_param[[r]], spec$param$bounds[r, 2L])
    }

    # H1 fixed and H4 short-circuited
    expect_named(res$hyp_spending_param, paste0("H", 1:4))
    expect_identical(
        unname(res$hyp_spending_param),
        c(NA, unname(res$spending_param), NA)
    )
    expect_identical(
        res$spending,
        resolve_pvalues_gsd(spec, res$spending_param)$spending
    )
    expect_identical(res$spending_bounds, spec$param$bounds)

    profile <- res$spending_profile
    expect_s3_class(profile, "data.frame")
    expect_named(
        profile,
        c("param", "value", "gain", "gain_fixed_graph", "diff", "se", "n_eval")
    )
    expect_identical(nrow(profile), sum(lengths(lattice)))
    expect_identical(profile$param, rep(c("gamma1", "b"), lengths(lattice)))
    expect_identical(profile$value, unname(unlist(lattice)))
    expect_true(all(profile$gain >= profile$gain_fixed_graph))
    expect_type(profile$n_eval, "integer")
    expect_length(attr(profile, "graphs"), nrow(profile))

    at_star <- profile$value == res$spending_param[profile$param]
    expect_identical(sum(at_star), 2L)
    expect_identical(profile$diff[at_star], c(0, 0))
    expect_identical(profile$se[at_star], c(0, 0))
    expect_type(res$profile_converged, "logical")
    expect_length(res$profile_converged, 1L)
    if (res$profile_converged) {
        expect_true(all(profile$diff >= 0))
    }

    search <- res$spending_search
    expect_named(
        search,
        c("param", "value", "visits", "best_fitness_subsample")
    )
    expect_identical(nrow(search), sum(lengths(lattice)))
    expect_type(search$visits, "integer")
    expect_gt(sum(search$visits[search$param == "b"]), 0L)
})

test_that("(d) a d = 1 profile passes through the optimum", {
    res <- opt_hsd_result()
    profile <- res$spending_profile
    at_star <- profile$value == res$spending_param[["sec"]]

    expect_identical(nrow(profile), 5L)
    expect_identical(profile$value, seq(-4, -2, by = 0.5))
    expect_identical(profile$diff[at_star], 0)
    expect_identical(profile$se[at_star], 0)
    expect_true(all(profile$gain >= profile$gain_fixed_graph))
    expect_true(res$profile_converged)
    # the reported graph re-scores at the optimum's own profile row
    expect_identical(
        res$spending_param,
        c(sec = res$spending_param[["sec"]])
    )
})

test_that("(e) global_search = FALSE starts from the default graph", {
    spec <- hsd_opt_spec()
    gc <- graph_constraint_free(3)
    res <- opt_hsd_run(seed = 404, global_search = FALSE)

    expect_s3_class(res, "multigrain_graph_optimal")
    expect_false(res$global_search)
    expect_null(res$global_output)
    expect_identical(res$solution$opt_source, "local")
    expect_true(is_graph_valid(res$hyp_weight, res$trans_matrix))
    expect_identical(
        res$local_output$x0,
        as.numeric(create_start_params(gc))
    )
    expect_true(all(is.na(res$spending_search$visits)))
    expect_true(all(is.na(res$spending_search$best_fitness_subsample)))
    # Stage L starts at the spec's start node and profiles every node
    expect_identical(nrow(res$spending_profile), 5L)
    expect_true(res$spending_param[["sec"]] %in% spec$param$lattice$sec)
})

test_that("(e) global_search = FALSE honours a user start_graph", {
    gc <- graph_constraint_free(3)
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
    res <- opt_hsd_run(seed = 404, global_search = FALSE, start_graph = start)

    expect_identical(res$start_graph, start)
    expect_identical(
        res$local_output$x0,
        .hsd_encode_graph(gc, start[[1]]$hyp_weight, start[[1]]$trans_matrix)
    )
    # not the Holm seed the fixed path would use
    expect_false(identical(
        res$local_output$x0,
        as.numeric(.build_start_matrix(gc, start)[1, ])
    ))
    expect_true(is_graph_valid(res$hyp_weight, res$trans_matrix))
})

test_that("(f) a stopping rule runs through the whole pipeline", {
    res <- opt_hsd_result_d2()

    expect_s3_class(res$stopping_rule, "multigrain_stopping_rule_gsd")
    expect_named(res$power$stop_distribution, c("analysis 1", "analysis 2"))
    # a sum of proportions: equal to floating point, not bitwise
    # nolint next: expect_identical_linter.
    expect_equal(sum(res$power$stop_distribution), 1)
    expect_gt(res$power$stop_distribution[[1L]], 0)
})

test_that("the default popSize counts each spending parameter", {
    spec <- hsd_opt_spec_d2()
    control <- multigrain_control() |>
        control_global(maxiter = 2, run = 2) |>
        control_local(maxeval = 50) |>
        control_nsim_global(300)
    res <- withr::with_seed(1, {
        graph_optimise_gsd(
            pvals = spec,
            graph_constraint = graph_constraint_free(4),
            trial_success = opt_hsd_gain_4,
            control = control,
            verbose = "silent"
        )
    })

    # m = 4 and d = 2: 40 * 6 = 240, where the fixed rule gives 200
    expect_identical(res$control$global_opt$popSize, 240L)
    expect_identical(res$global_output@popSize, 240L)
})

test_that("a user popSize survives", {
    control <- hsd_opt_control() |> control_global(popSize = 31)
    res <- withr::with_seed(2, {
        graph_optimise_gsd(
            pvals = hsd_opt_spec(),
            graph_constraint = graph_constraint_free(3),
            trial_success = opt_hsd_gain_3,
            control = control,
            verbose = "silent"
        )
    })

    expect_identical(res$control$global_opt$popSize, 31)
    expect_identical(res$global_output@popSize, 31)
})

test_that("the dispatch checks dots and validates on the spec", {
    spec <- hsd_opt_spec()

    expect_error(
        graph_optimise_gsd(
            spec,
            graph_constraint_free(3),
            opt_hsd_gain_3,
            unused = 1
        ),
        class = "rlib_error_dots_nonempty"
    )
    expect_error(
        graph_optimise_gsd(
            spec,
            graph_constraint_free(2),
            opt_hsd_gain_3,
            verbose = "silent"
        ),
        "disagree on the number of"
    )
    expect_error(
        graph_optimise_gsd(
            spec,
            graph_constraint_free(3),
            opt_hsd_gain_3,
            alpha = 0.05,
            verbose = "silent"
        ),
        "between 0 and 0.025"
    )
})

test_that("batch-means differences are paired block means", {
    time_a <- matrix(rep_len(c(1L, 2L), 400L), nrow = 200L)
    time_b <- time_a
    time_b[1:20, 1L] <- 2L
    func <- function(time) mean(time[, 1L] == 1L)

    same <- .hsd_block_se(time_a, time_a, func)
    expect_identical(same, list(diff = 0, se = 0))

    out <- .hsd_block_se(time_a, time_b, func, B = 10)
    blocks <- vapply(
        split(seq_len(200L), rep(1:10, each = 20L)),
        function(i) {
            func(time_a[i, , drop = FALSE]) - func(time_b[i, , drop = FALSE])
        },
        double(1L)
    )
    expect_identical(out$diff, mean(blocks))
    expect_identical(out$se, stats::sd(blocks) / sqrt(10))
    # never more blocks than rows
    expect_length(.hsd_block_gains(time_a[1:5, ], func = func, B = 100), 5L)
})

test_that("the plateau is the connected run around the optimum", {
    value <- seq(-4, 0, by = 0.5)
    diff <- c(1e-3, 4e-4, 0, 1e-4, 6e-4, 2e-4, 1e-4, 3e-3, 0)

    expect_identical(.hsd_plateau(value, diff, at = -3), c(-3.5, -2.5))
    expect_identical(.hsd_plateau(value, diff, tol = 1e-3, at = -3), c(-4, -1))
    expect_identical(
        .hsd_plateau(rev(value), rev(diff), at = -3),
        c(-3.5, -2.5)
    )
})

test_that("graph_optimal() attaches a named spending list only", {
    base <- list(
        hyp_weight = c(0.5, 0.5),
        trans_matrix = matrix(c(0, 1, 1, 0), nrow = 2L)
    )
    plain <- do.call(graph_optimal, base)
    expect_named(plain, names(new_graph_optimal()))

    with_spending <- do.call(
        graph_optimal,
        c(base, list(spending = list(spending_param = c(gamma = -2))))
    )
    expect_identical(with_spending$spending_param, c(gamma = -2))

    expect_error(
        do.call(graph_optimal, c(base, list(spending = list(-2)))),
        "named list"
    )
})

test_that("print() and summary() of a real result succeed", {
    res <- opt_hsd_result()

    expect_output(print(res), "Spending functions:", fixed = TRUE)
    expect_output(
        summary(res),
        "Spending-parameter profile (sec)",
        fixed = TRUE
    )
})

test_that("the P3 files use no superassignment", {
    files <- c(
        test_path("..", "..", "R", "optimisation_gsd_opt.R"),
        test_path("..", "..", "R", "objective_function_gsd_opt.R"),
        test_path("..", "..", "R", "transform_pvalues_gsd_opt.R")
    )
    files <- files[file.exists(files)]
    skip_if(length(files) == 0L, "package sources not available")

    for (f in files) {
        code <- readLines(f, warn = FALSE)
        code <- sub("#.*$", "", code)
        expect_false(any(grepl("<<-", code, fixed = TRUE)), label = f)
    }
})
# nolint end
