# Phase P2 of dev/hsd_opt_build_plan.md: the objective closure of the joint
# graph and spending-parameter search, and the slice bank it reads. The gate
# is the node oracle: at every node the closure is the fixed closure on the
# resolved p-values. Fixtures are in helper-hsd_opt.R.

opt_obj_start <- list(list(hyp_weight = NULL, trans_matrix = NULL))

# Everything a test needs for one spec: rows, bank, gain, constraint, graph
# gene count, a closure factory and the fixed reference at a node.
opt_obj_setup <- function(spec, gain) {
    idx <- withr::with_seed(3, sample.int(spec$nsim, 1000))
    constraint <- graph_constraint_free(spec$m)
    bank <- .hsd_bank(spec, idx)
    p_graph <- ncol(.build_start_matrix(constraint, opt_obj_start))

    make <- function(
        num_threads = 1L,
        stop_rule = NULL,
        start_param = spec$param$start
    ) {
        create_obj_func_gsd_opt(
            spec,
            bank = bank,
            power_criterion = gain$func,
            hyp_constraint = constraint$hyp_constraint,
            trans_constraint = constraint$trans_constraint,
            alpha = spec$alpha,
            num_threads = num_threads,
            stop_rule = stop_rule,
            p_graph = p_graph,
            start_param = start_param
        )
    }

    ref <- function(node, num_threads = 1L, stop_rule = NULL) {
        resolved <- resolve_pvalues_gsd(spec, node)
        resolved$pvals <- resolved$pvals[idx, , , drop = FALSE]
        resolved$nsim <- length(idx)
        create_obj_func_gsd(
            spec$m,
            gain$func,
            constraint$hyp_constraint,
            constraint$trans_constraint,
            pvals = .gsd_kernel_matrix(resolved),
            K = spec$K,
            alpha = spec$alpha,
            num_threads = num_threads,
            stop_rule = stop_rule
        )
    }

    list(
        spec = spec,
        idx = idx,
        constraint = constraint,
        bank = bank,
        p_graph = p_graph,
        make = make,
        ref = ref
    )
}

opt_obj_gain_3 <- trial_success_gsd(r1 + r2 + r3, K = 2, verbose = "silent")

opt_obj_gain_4 <- trial_success_gsd(
    d(t1) + d(t2) + d(t3) + d(t4),
    d = c(1, 0.9),
    verbose = "silent"
)

# TRUE when graph genes `x` decode to a graph the closure scores, i.e. when
# the evaluation reaches the kernel (the penalty checks of the closure).
opt_obj_feasible <- function(x, constraint) {
    theta <- split_theta(x, constraint$hyp_constraint)
    w <- recover_full_weights(theta$w_pars, constraint$hyp_constraint)
    G <- recover_full_trans_matrix(theta$g_pars, constraint$trans_constraint)

    !anyNA(w) && !anyNA(G) && all(w >= 0 & w <= 1) && all(G >= 0 & G <= 1)
}

opt_obj_d1 <- function() {
    hsd_opt_cached(
        # nolint: object_usage_linter.
        "obj_setup_d1",
        opt_obj_setup(hsd_opt_spec(), opt_obj_gain_3) # nolint: object_usage_linter.
    )
}

opt_obj_d2 <- function() {
    hsd_opt_cached(
        # nolint: object_usage_linter.
        "obj_setup_d2",
        opt_obj_setup(hsd_opt_spec_d2(), opt_obj_gain_4) # nolint: object_usage_linter.
    )
}

# `n` random graph gene vectors, one per row. Uniform genes on [0, 1] mostly
# decode to infeasible graphs (the free weights must sum to at most 1), which
# return a penalty before the kernel; so three rows in four are drawn on
# [0, 1 / (m - 1)], which always decodes to a feasible graph under
# `graph_constraint_free(m)`, and every fourth on [0, 1].
opt_obj_genes <- function(p_graph, m, n = 20L, seed = 11) {
    withr::with_seed(seed, {
        genes <- matrix(stats::runif(n * p_graph), nrow = n)
        narrow <- seq_len(n) %% 4L != 0L
        genes[narrow, ] <- genes[narrow, ] / (m - 1)
        genes
    })
}

# Every node combination of a spec: a matrix with one row per combination,
# columns the node indices per parameter.
opt_obj_node_grid <- function(spec) {
    as.matrix(expand.grid(lapply(spec$param$lattice, seq_along)))
}

# Run the node oracle on `setup` and return the number of comparisons made.
opt_obj_oracle <- function(setup, num_threads = 1L, stop_rule = NULL) {
    spec <- setup$spec
    lattice <- spec$param$lattice
    n_nodes <- lengths(lattice)
    fn <- setup$make(num_threads = num_threads, stop_rule = stop_rule)
    genes <- opt_obj_genes(setup$p_graph, spec$m)
    nodes <- opt_obj_node_grid(spec)
    n_cmp <- 0L
    reached <- apply(genes, 1L, opt_obj_feasible, setup$constraint)
    # not vacuous: most comparisons reach the kernel
    expect_gte(sum(reached), 15L) # nolint: object_usage_linter.

    for (k in seq_len(nrow(nodes))) {
        at <- nodes[k, ]
        node <- vapply(seq_along(at), function(r) lattice[[r]][[at[[r]]]], 1)
        ref <- setup$ref(node, num_threads = 1L, stop_rule = stop_rule)
        u <- (at - 0.5) / n_nodes
        for (j in seq_len(nrow(genes))) {
            x <- genes[j, ]
            expect_identical(fn(c(x, u)), ref(x)) # nolint: object_usage_linter.
            n_cmp <- n_cmp + 1L
        }
    }

    n_cmp
}


test_that(".hsd_node_index() gives equal mass per node and clamps", {
    expect_identical(
        .hsd_node_index(c(0, 0.19, 0.2, 0.99, 1), 5),
        c(1L, 1L, 2L, 5L, 5L)
    )
    expect_identical(.hsd_node_index(c(-0.3, 1.7), 5), c(1L, 5L))
    expect_identical(.hsd_node_index((seq_len(21) - 0.5) / 21, 21), 1:21)
    expect_identical(.hsd_node_index(c(0.5, 0.5), c(3, 5)), c(2L, 3L))
})

test_that("node oracle: d = 1, every node equals the fixed closure", {
    expect_identical(opt_obj_oracle(opt_obj_d1()), 5L * 20L)
})

test_that("node oracle: d = 2, every node pair equals the fixed closure", {
    setup <- opt_obj_d2()
    # H1 fixed, H2 and H3 each on their own parameter, H4 short-circuited
    expect_identical(setup$bank$hyp, 2:3)
    expect_identical(opt_obj_oracle(setup), 9L * 20L)
})

test_that("node oracle holds with a stopping rule and with two threads", {
    spec <- hsd_opt_spec()
    rule <- .gsd_stop_table(stopping_rule_gsd(r1)$conditions, spec$m, spec$K)

    expect_identical(opt_obj_oracle(opt_obj_d1(), stop_rule = rule), 100L)
    expect_identical(opt_obj_oracle(opt_obj_d1(), num_threads = 2L), 100L)
    expect_identical(
        opt_obj_oracle(opt_obj_d2(), num_threads = 2L, stop_rule = NULL),
        180L
    )
})

test_that("parameter genes outside [0, 1] map to the end nodes", {
    setup <- opt_obj_d1()
    fn <- setup$make()
    lattice <- setup$spec$param$lattice[[1L]]
    ref_1 <- setup$ref(lattice[[1L]])
    ref_m <- setup$ref(lattice[[5L]])
    genes <- opt_obj_genes(setup$p_graph, setup$spec$m, n = 5L)

    for (j in seq_len(nrow(genes))) {
        x <- genes[j, ]
        expect_identical(fn(c(x, -0.3)), ref_1(x))
        expect_identical(fn(c(x, 1.7)), ref_m(x))
        expect_identical(fn(c(x, 0)), ref_1(x))
        expect_identical(fn(c(x, 1)), ref_m(x))
    }

    # the end nodes really differ, so the checks are not vacuous
    values_1 <- apply(genes, 1L, ref_1)
    values_m <- apply(genes, 1L, ref_m)
    expect_false(identical(values_1, values_m))
})

test_that("an infeasible graph returns the penalty, no kernel call", {
    setup <- opt_obj_d1()
    fn <- setup$make()
    ref <- setup$ref(setup$spec$param$start)
    arch <- attr(fn, "archive")

    # free constraint, m = 3: genes 1-2 are the weights; a negative gene gives
    # a negative weight
    x <- c(-0.4, 0.3, 0.5, 0.5, 0.5)
    expect_length(x, setup$p_graph)
    penalty <- fn(c(x, 0.1))

    expect_identical(penalty, ref(x))
    expect_lt(penalty, 0)
    expect_identical(arch$visits[[1L]], integer(5L))
    expect_identical(arch$best[[1L]], rep(-Inf, 5L))
})

test_that("the archive matches a by-hand replay of 30 evaluations", {
    for (setup in list(opt_obj_d1(), opt_obj_d2())) {
        spec <- setup$spec
        n_nodes <- lengths(spec$param$lattice, use.names = FALSE)
        d <- length(n_nodes)
        fn <- setup$make()
        arch <- attr(fn, "archive")

        expect_identical(attr(fn, "d"), d)
        expect_identical(attr(fn, "p_graph"), setup$p_graph)

        best <- lapply(n_nodes, function(n) rep(-Inf, n))
        visits <- lapply(n_nodes, integer)
        genes <- lapply(n_nodes, function(n) {
            matrix(NA_real_, nrow = n, ncol = setup$p_graph)
        })

        gene_draws <- opt_obj_genes(setup$p_graph, spec$m, n = 30L, seed = 5)
        # every fifth evaluation is infeasible and must not be recorded
        gene_draws[seq(5L, 30L, by = 5L), 1L] <- -0.2
        reached <- apply(gene_draws, 1L, opt_obj_feasible, setup$constraint)
        expect_false(any(reached[seq(5L, 30L, by = 5L)]))
        expect_gte(sum(reached), 18L)
        u_draws <- withr::with_seed(6, {
            matrix(stats::runif(30L * d, -0.2, 1.2), nrow = 30L)
        })

        for (j in seq_len(30L)) {
            xg <- gene_draws[j, ]
            fit <- fn(c(xg, u_draws[j, ]))
            if (!reached[[j]]) {
                next
            }
            at <- .hsd_node_index(u_draws[j, ], n_nodes)
            for (r in seq_len(d)) {
                visits[[r]][[at[[r]]]] <- visits[[r]][[at[[r]]]] + 1L
                if (fit > best[[r]][[at[[r]]]]) {
                    best[[r]][[at[[r]]]] <- fit
                    genes[[r]][at[[r]], ] <- xg
                }
            }
        }

        names(best) <- names(visits) <- names(genes) <- spec$param$names
        expect_identical(arch$best, best)
        expect_identical(arch$visits, visits)
        expect_identical(arch$genes, genes)
        expect_identical(sum(arch$visits[[1L]]), sum(reached))
    }
})

test_that("skipping unchanged nodes does not change results", {
    setup <- opt_obj_d2()
    spec <- setup$spec
    lattice <- spec$param$lattice
    n_nodes <- lengths(lattice, use.names = FALSE)
    fn <- setup$make(start_param = c(-2.5, -3))
    refs <- list()
    genes <- opt_obj_genes(setup$p_graph, spec$m, n = 50L, seed = 9)
    # runs of repeated node pairs, alternating between them
    q_seq <- withr::with_seed(10, {
        cbind(
            rep(sample.int(3L, 10L, replace = TRUE), each = 5L),
            rep(c(1L, 1L, 3L, 2L, 2L), 10L)
        )
    })

    for (j in seq_len(50L)) {
        at <- q_seq[j, ]
        key <- paste(at, collapse = ",")
        if (is.null(refs[[key]])) {
            refs[[key]] <- setup$ref(c(
                lattice[[1L]][[at[[1L]]]],
                lattice[[2L]][[at[[2L]]]]
            ))
        }
        x <- genes[j, ]
        expect_identical(fn(c(x, (at - 0.5) / n_nodes)), refs[[key]](x))
    }
})

test_that("the working matrix is not copied per evaluation", {
    skip_if_not(capabilities("profmem"))

    for (setup in list(opt_obj_d1(), opt_obj_d2())) {
        spec <- setup$spec
        n_nodes <- lengths(spec$param$lattice, use.names = FALSE)
        fn <- setup$make()
        genes <- opt_obj_genes(setup$p_graph, spec$m, n = 10L, seed = 4)
        # every parameter moves at every evaluation
        nodes <- (seq_len(10L) %% 2L) * (n_nodes[[1L]] - 1L) + 1L
        state <- attr(fn, "archive")
        expect_identical(state, environment(fn)$state)

        out <- utils::capture.output({
            tracemem(state$mat)
            for (j in seq_len(10L)) {
                at <- pmin(nodes[[j]], n_nodes)
                fn(c(genes[j, ], (at - 0.5) / n_nodes))
            }
            untracemem(state$mat)
        })

        expect_identical(out, character())
        # not vacuous: the evaluations reached the kernel and moved nodes
        reached <- apply(genes, 1L, opt_obj_feasible, setup$constraint)
        expect_gte(sum(reached), 8L)
        expect_identical(sum(state$visits[[1L]]), sum(reached))
        expect_identical(state$last_q, rep(1L, length(n_nodes)))
    }
})

test_that("the closure factory uses no superassignment", {
    path <- test_path("..", "..", "R", "objective_function_gsd_opt.R")
    skip_if_not(file.exists(path))

    expect_false(any(grepl("<<-", readLines(path), fixed = TRUE)))
})

test_that("bank slices, byte count and the memory guard", {
    for (setup in list(opt_obj_d1(), opt_obj_d2())) {
        spec <- setup$spec
        bank <- setup$bank
        n_nodes <- lengths(spec$param$lattice, use.names = FALSE)
        hyp <- which(!is.na(spec$group))

        expect_identical(bank$idx, setup$idx)
        expect_identical(bank$hyp, hyp)
        expect_length(bank$slices, length(hyp))
        for (a in seq_along(hyp)) {
            slices <- bank$slices[[a]]
            expect_length(slices, n_nodes[[spec$group[[hyp[[a]]]]]])
            for (s in slices) {
                expect_identical(dim(s), c(1000L, spec$K))
            }
        }
        expect_identical(
            bank$bytes,
            8 * spec$K * 1000 * sum(n_nodes[spec$group[hyp]])
        )
    }
    # d = 1: two hypotheses on 5 nodes; d = 2: two hypotheses on 3 nodes
    expect_identical(opt_obj_d1()$bank$bytes, 8 * 2 * 1000 * 10)
    expect_identical(opt_obj_d2()$bank$bytes, 8 * 2 * 1000 * 6)

    expect_no_condition(.check_bank_memory(1.9e9))
    expect_snapshot(.check_bank_memory(2.5e9))
    expect_snapshot(error = TRUE, .check_bank_memory(4.5e9))
    expect_snapshot(error = TRUE, .check_bank_memory(1.5e9, max_bank_gb = 1))
})
