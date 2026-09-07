# Shared fixture: m = 4, nsim = 1e4, seeded so every test sees the same trials.
simplify_pvals <- withr::with_seed(2, {
    corr <- matrix(0.2, nrow = 4, ncol = 4)
    diag(corr) <- 1
    sims <- mvtnorm::rmvnorm(
        1e4,
        mean = calc_ncp(c(0.93, 0.91, 0.90, 0.85)),
        sigma = corr
    )
    stats::pnorm(sims, lower.tail = FALSE)
})

simplify_ts <- trial_success(
    0.25 * (2 * (r1 && r2) + r1 * r3 + r2 * r4),
    verbose = "silent"
)

# Deliberately small budgets: these tests check invariants, not search quality.
simplify_control <- function() {
    multigrain_control() |>
        control_global(run = 4, maxiter = 8, popSize = 20) |>
        control_local(maxeval = 400)
}

simplify_reference <- function(seed = 1, global_search = TRUE) {
    withr::with_seed(seed, {
        graph_optimise(
            pvals = simplify_pvals,
            graph_constraint = graph_constraint_free(4),
            trial_success = simplify_ts,
            control = simplify_control(),
            global_search = global_search,
            verbose = "silent"
        )
    })
}


# ---- the cap, end to end ----

test_that("graph_simplify() holds the cap and never adds edges", {
    ref <- simplify_reference()

    for (lambda in c(0, 1e-3, 1e-2, 1)) {
        res <- withr::with_seed(7, {
            graph_simplify(
                ref,
                simplify_pvals,
                gain_tolerance = lambda,
                control = simplify_control(),
                verbose = "silent"
            )
        })
        sp <- res$sparsity

        expect_s3_class(res, "multigrain_graph_optimal")
        expect_true(is_graph_valid(
            unname(res$hyp_weight),
            unname(res$trans_matrix)
        ))

        # The cap is an exact inequality, checked with no tolerance at all.
        expect_gte(sp$gain, (1 - lambda) * sp$gain_reference)
        expect_lte(sp$n_edges, sp$n_edges_reference)
        expect_lte(sp$n_edges_free, sp$n_edges_free_reference)

        # $power describes the graph returned, not the one the search saw.
        fresh <- calc_power_pvals(
            simplify_pvals,
            hyp_weight = unname(res$hyp_weight),
            trans_matrix = unname(res$trans_matrix),
            alpha = res$alpha,
            custom_power = list(trial_success = simplify_ts)
        )
        expect_identical(res$power$trial_success, fresh$trial_success)
        expect_identical(sp$gain, res$power$trial_success)

        # Internal consistency of the reported numbers.
        expect_identical(sp$gain_loss, sp$gain_reference - sp$gain)
        expect_identical(
            sp$gain_loss_fraction,
            sp$gain_loss / sp$gain_reference
        )
        expect_identical(sp$budget, lambda * sp$gain_reference)
        expect_identical(sp$gain_tolerance, lambda)
    }
})

test_that("gain_tolerance = 1 leaves exactly one free edge per row", {
    ref <- simplify_reference()
    res <- withr::with_seed(7, {
        graph_simplify(
            ref,
            simplify_pvals,
            gain_tolerance = 1,
            control = simplify_control(),
            verbose = "silent"
        )
    })

    expect_identical(res$sparsity$n_edges_free, 4L)
    expect_identical(
        rowSums(unname(res$trans_matrix) != 0),
        rep(1, 4)
    )
    expect_false(identical(res$sparsity$source, "reference"))
})

test_that("the cap holds across seeds and budgets", {
    ref <- simplify_reference()

    for (seed in 1:5) {
        for (lambda in c(0, 1e-3, 5e-3, 1e-2, 1)) {
            res <- withr::with_seed(seed, {
                graph_simplify(
                    ref,
                    simplify_pvals,
                    gain_tolerance = lambda,
                    global_search = FALSE,
                    control = simplify_control(),
                    verbose = "silent"
                )
            })
            expect_gte(
                res$sparsity$gain,
                (1 - lambda) * res$sparsity$gain_reference
            )
            expect_lte(res$sparsity$n_edges, res$sparsity$n_edges_reference)
        }
    }
})

test_that("a returned reference is exactly the graph that went in", {
    ref <- simplify_reference()
    res <- withr::with_seed(7, {
        graph_simplify(
            ref,
            simplify_pvals,
            gain_tolerance = 0,
            global_search = FALSE,
            control = simplify_control(),
            verbose = "silent"
        )
    })

    if (identical(res$sparsity$source, "reference")) {
        expect_identical(unname(res$hyp_weight), unname(ref$hyp_weight))
        expect_identical(unname(res$trans_matrix), unname(ref$trans_matrix))
        expect_identical(res$solution$opt_source, "reference")
        expect_identical(res$sparsity$gain, res$sparsity$gain_reference)
    } else {
        # Otherwise the result must be strictly better under the same rule.
        expect_true(
            res$sparsity$n_edges < res$sparsity$n_edges_reference ||
                res$sparsity$gain >= res$sparsity$gain_reference
        )
    }
})


# ---- reporting ----

test_that("the sparsity element carries the reference graph and its names", {
    ref <- simplify_reference()
    res <- withr::with_seed(7, {
        graph_simplify(
            ref,
            simplify_pvals,
            gain_tolerance = 1e-2,
            global_search = FALSE,
            control = simplify_control(),
            verbose = "silent"
        )
    })

    gc_names <- graph_constraint_get_names(ref$constraints)
    expect_named(res$sparsity$reference$hyp_weight, gc_names)
    expect_identical(
        dimnames(res$sparsity$reference$trans_matrix),
        list(gc_names, gc_names)
    )
    expect_identical(
        unname(res$sparsity$reference$trans_matrix),
        unname(ref$trans_matrix)
    )
    expect_named(res$hyp_weight, gc_names)

    expect_setequal(
        names(res$sparsity),
        c(
            "gain_tolerance", "reference", "gain_reference", "gain",
            "gain_loss", "gain_loss_fraction", "budget",
            "n_edges_reference", "n_edges", "n_edges_free_reference",
            "n_edges_free", "prune_loss", "source"
        )
    )
    expect_true(res$sparsity$source %in% c("global", "local", "reference"))
})

test_that("the result prints its simplification", {
    ref <- simplify_reference()
    res <- withr::with_seed(7, {
        graph_simplify(
            ref,
            simplify_pvals,
            gain_tolerance = 1e-2,
            global_search = FALSE,
            control = simplify_control(),
            verbose = "silent"
        )
    })
    out <- utils::capture.output(print(res))
    expect_true(any(grepl("Simplified from", out, fixed = TRUE)))
    expect_true(any(grepl("cap 1.0%", out, fixed = TRUE)))
})


# ---- inputs ----

test_that("alpha comes from the object and can be overridden", {
    ref <- simplify_reference()
    expect_identical(ref$alpha, 0.025)

    res <- withr::with_seed(7, {
        graph_simplify(
            ref,
            simplify_pvals,
            gain_tolerance = 1e-3,
            alpha = 0.01,
            global_search = FALSE,
            control = simplify_control(),
            verbose = "silent"
        )
    })
    expect_identical(res$alpha, 0.01)

    # The reference gain must be measured at the alpha actually used.
    expect_identical(
        res$sparsity$gain_reference,
        calc_power_pvals(
            simplify_pvals,
            hyp_weight = unname(ref$hyp_weight),
            trans_matrix = unname(ref$trans_matrix),
            alpha = 0.01,
            custom_power = list(trial_success = simplify_ts)
        )$trial_success
    )
})

test_that("an object without alpha errors unless alpha is supplied", {
    ref <- simplify_reference()
    ref$alpha <- NULL

    expect_error(
        graph_simplify(ref, simplify_pvals, verbose = "silent"),
        "does not store"
    )
    expect_no_error(
        withr::with_seed(7, {
            graph_simplify(
                ref,
                simplify_pvals,
                alpha = 0.025,
                global_search = FALSE,
                control = simplify_control(),
                verbose = "silent"
            )
        })
    )
})

test_that("pvals with the wrong number of columns errors", {
    ref <- simplify_reference()
    expect_error(
        graph_simplify(ref, simplify_pvals[, 1:3], verbose = "silent"),
        "one column per hypothesis"
    )
})

test_that("a non-positive reference gain errors before any search", {
    zero_ts <- trial_success(0 * r1 + 0 * r2 + 0 * r3 + 0 * r4,
        verbose = "silent"
    )
    obj <- graph_optimal(
        hyp_weight = c(0.25, 0.25, 0.25, 0.25),
        trans_matrix = matrix(1 / 3, 4, 4) - diag(1 / 3, 4),
        constraints = graph_constraint_free(4),
        trial_success = zero_ts,
        alpha = 0.025
    )
    expect_error(
        graph_simplify(obj, simplify_pvals, verbose = "silent"),
        "non-positive trial success"
    )
})

test_that("an object without its constraint or trial success errors", {
    ref <- simplify_reference()
    stripped <- ref
    stripped$constraints <- NULL
    expect_error(
        graph_simplify(stripped, simplify_pvals, verbose = "silent"),
        "constraint and trial-success"
    )
})

test_that("num_threads does not change the result for the same seed", {
    ref <- simplify_reference()
    args <- list(
        gain_tolerance = 1e-2,
        global_search = FALSE,
        control = simplify_control(),
        verbose = "silent"
    )
    one <- withr::with_seed(7, {
        do.call(graph_simplify, c(list(ref, simplify_pvals), args,
            list(num_threads = 1L)
        ))
    })
    two <- withr::with_seed(7, {
        do.call(graph_simplify, c(list(ref, simplify_pvals), args,
            list(num_threads = 2L)
        ))
    })

    expect_identical(one$hyp_weight, two$hyp_weight)
    expect_identical(one$trans_matrix, two$trans_matrix)
    expect_identical(one$sparsity$gain, two$sparsity$gain)
})


# ---- degenerate cases ----

test_that("nothing removable warns and returns the reference", {
    # m = 2: each row has a single free entry, so no edge can go.
    pv2 <- withr::with_seed(3, {
        sims <- mvtnorm::rmvnorm(
            2000,
            mean = calc_ncp(c(0.9, 0.8)),
            sigma = diag(2)
        )
        stats::pnorm(sims, lower.tail = FALSE)
    })
    ts2 <- trial_success(r1 + r2, verbose = "silent")
    ref2 <- withr::with_seed(1, {
        graph_optimise(
            pvals = pv2,
            graph_constraint = graph_constraint_free(2),
            trial_success = ts2,
            control = multigrain_control() |> control_local(maxeval = 50),
            global_search = FALSE,
            verbose = "silent"
        )
    })

    expect_warning(
        res <- graph_simplify(ref2, pv2, verbose = "silent"),
        "No removable edges"
    )
    expect_identical(unname(res$hyp_weight), unname(ref2$hyp_weight))
    expect_identical(unname(res$trans_matrix), unname(ref2$trans_matrix))
    expect_identical(res$sparsity$source, "reference")
    expect_identical(res$sparsity$gain_loss, 0)
    expect_identical(res$sparsity$prune_loss, 0)
    expect_identical(res$solution$opt_source, "reference")
})

test_that("a different pvals recomputes the reference gain on it", {
    ref <- simplify_reference()
    other <- withr::with_seed(99, {
        corr <- matrix(0.2, nrow = 4, ncol = 4)
        diag(corr) <- 1
        sims <- mvtnorm::rmvnorm(
            1e4,
            mean = calc_ncp(c(0.93, 0.91, 0.90, 0.85)),
            sigma = corr
        )
        stats::pnorm(sims, lower.tail = FALSE)
    })

    res <- withr::with_seed(7, {
        graph_simplify(
            ref,
            other,
            gain_tolerance = 1e-2,
            global_search = FALSE,
            control = simplify_control(),
            verbose = "silent"
        )
    })

    expect_identical(
        res$sparsity$gain_reference,
        calc_power_pvals(
            other,
            hyp_weight = unname(ref$hyp_weight),
            trans_matrix = unname(ref$trans_matrix),
            alpha = 0.025,
            custom_power = list(trial_success = simplify_ts)
        )$trial_success
    )
    # ... and the cap holds on the sample it was measured on.
    expect_gte(res$sparsity$gain, 0.99 * res$sparsity$gain_reference)
})

test_that("simplifying twice compounds the cap", {
    ref <- simplify_reference()
    once <- withr::with_seed(7, {
        graph_simplify(
            ref,
            simplify_pvals,
            gain_tolerance = 1e-2,
            global_search = FALSE,
            control = simplify_control(),
            verbose = "silent"
        )
    })
    twice <- withr::with_seed(7, {
        graph_simplify(
            once,
            simplify_pvals,
            gain_tolerance = 1e-2,
            global_search = FALSE,
            control = simplify_control(),
            verbose = "silent"
        )
    })

    # The second call's reference is the first call's result.
    expect_identical(twice$sparsity$gain_reference, once$sparsity$gain)
    expect_gte(twice$sparsity$gain, 0.99 * twice$sparsity$gain_reference)
    expect_lte(twice$sparsity$n_edges, once$sparsity$n_edges)
})

test_that("a stored population of the wrong width is ignored", {
    # A8: the object was optimised under a different constraint.
    ref <- simplify_reference()
    other <- simplify_reference(seed = 2)
    ref$global_output <- other$global_output
    ref$global_output@population <- other$global_output@population[, 1:3]

    expect_no_error(
        withr::with_seed(7, {
            graph_simplify(
                ref,
                simplify_pvals,
                gain_tolerance = 1e-3,
                control = simplify_control(),
                verbose = "silent"
            )
        })
    )
})

test_that("a stored run of 1 still halves to a legal value", {
    # A10: halving floors at 1 rather than reaching 0.
    ref <- simplify_reference()
    ref$control$global_opt$run <- 1L

    expect_no_error(
        res <- withr::with_seed(7, {
            graph_simplify(
                ref,
                simplify_pvals,
                gain_tolerance = 1e-3,
                verbose = "silent"
            )
        })
    )
    expect_gte(res$control$global_opt$run, 1L)
})
