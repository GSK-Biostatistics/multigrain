# Phase P1 of dev/hsd_opt_build_plan.md: the spending spec built by
# `transform_pvalues_gsd()` when a hypothesis is tagged with `sf_hsd_opt()`,
# and `resolve_pvalues_gsd()`. Fixtures are in helper-hsd_opt.R.

test_that("sf_hsd_opt() returns a tagged sfHSD spending function", {
    f <- sf_hsd_opt(id = 3, bounds = c(-4L, 1L), step = 0.5, start = -1)

    expect_s3_class(f, "multigrain_sf_hsd_opt")
    expect_type(f, "closure")
    expect_true(is_sf_hsd_opt(f))
    expect_false(is_sf_hsd_opt(gsDesign::sfLDOF))
    expect_identical(attr(f, "id"), "3")
    expect_identical(attr(f, "bounds"), c(-4, 1))
    expect_identical(attr(f, "step"), 0.5)
    expect_identical(attr(f, "start"), -1)
    expect_identical(
        f(0.025, c(0.5, 1), param = -2),
        gsDesign::sfHSD(0.025, c(0.5, 1), -2)$spend
    )

    g <- sf_hsd_opt()
    expect_null(attr(g, "id"))
    expect_null(attr(g, "start"))
    expect_identical(attr(g, "bounds"), c(-4, 1))
    expect_identical(attr(g, "step"), 0.25)

    expect_true(.has_hsd_opt(g))
    expect_true(.has_hsd_opt(list(gsDesign::sfLDOF, g)))
    expect_false(.has_hsd_opt(list(gsDesign::sfLDOF, gsDesign::sfLDPocock)))
    expect_false(.has_hsd_opt(gsDesign::sfLDOF))
})

test_that("sf_hsd_opt() validation messages name the limits", {
    expect_snapshot(error = TRUE, sf_hsd_opt(bounds = c(-11, 1)))
    expect_snapshot(error = TRUE, sf_hsd_opt(bounds = c(-4, 2.5)))
    expect_snapshot(error = TRUE, sf_hsd_opt(bounds = c(-1, -2)))
    expect_snapshot(error = TRUE, sf_hsd_opt(bounds = c(-1, -1)))
    expect_snapshot(error = TRUE, sf_hsd_opt(bounds = c(-4, Inf)))
    expect_snapshot(error = TRUE, sf_hsd_opt(bounds = -4))
    expect_snapshot(error = TRUE, sf_hsd_opt(step = 0.1))
    expect_snapshot(error = TRUE, sf_hsd_opt(step = c(0.5, 0.25)))
    expect_snapshot(error = TRUE, sf_hsd_opt(start = 1.5))
    expect_snapshot(error = TRUE, sf_hsd_opt(start = NA_real_))
    expect_snapshot(error = TRUE, sf_hsd_opt(id = c("a", "b")))
    expect_snapshot(error = TRUE, sf_hsd_opt(id = NA))

    # the corners of the allowed range are accepted
    expect_no_error(sf_hsd_opt(bounds = c(-10, 2), step = 0.125))
})

test_that("group members must agree on bounds, step and start", {
    raw <- hsd_opt_raw()

    expect_snapshot(
        error = TRUE,
        transform_pvalues_gsd(
            raw,
            spending = list(
                gsDesign::sfLDOF,
                sf_hsd_opt(id = "sec", bounds = c(-4, -2), step = 0.5),
                sf_hsd_opt(id = "sec", bounds = c(-4, -1), step = 0.5)
            ),
            grid_size = 16L
        )
    )
    expect_snapshot(
        error = TRUE,
        transform_pvalues_gsd(
            raw,
            spending = list(
                gsDesign::sfLDOF,
                sf_hsd_opt(id = "sec", step = 0.5),
                sf_hsd_opt(id = "sec", step = 0.25)
            ),
            grid_size = 16L
        )
    )
    expect_snapshot(
        error = TRUE,
        transform_pvalues_gsd(
            raw,
            spending = list(
                gsDesign::sfLDOF,
                sf_hsd_opt(id = "gamma2", bounds = c(-4, -2), step = 0.5),
                sf_hsd_opt(bounds = c(-4, -2), step = 0.5)
            ),
            grid_size = 16L
        )
    )
})

test_that("a single recycled tag is one shared parameter", {
    raw <- hsd_opt_raw()
    tag <- sf_hsd_opt(bounds = c(-4, -2), step = 0.5)

    spec <- transform_pvalues_gsd(raw, spending = tag, grid_size = 16L)
    expect_s3_class(spec, "multigrain_pvals_gsd_spec")
    expect_identical(spec$group, rep(1L, 3L))
    expect_identical(spec$param$names, "gamma")
    expect_length(spec$param$lattice, 1L)

    # the same object listed twice without an id is two parameters
    spec2 <- transform_pvalues_gsd(
        raw,
        spending = list(gsDesign::sfLDOF, tag, tag),
        grid_size = 16L
    )
    expect_identical(spec2$group, c(NA, 1L, 2L))
    expect_identical(spec2$param$names, c("gamma1", "gamma2"))

    # ids are compared as strings
    groups <- .hsd_groups(
        list(sf_hsd_opt(id = 1), gsDesign::sfLDOF, sf_hsd_opt(id = "1")),
        m = 3L
    )
    expect_identical(groups$group, c(1L, NA, 1L))
    expect_identical(groups$names, "1")
})

test_that("lattices are exact multiples of step plus both bounds", {
    expect_identical(.hsd_lattice(c(-4, 1), 0.25), (-16:4) * 0.25)
    expect_identical(.hsd_lattice(c(-4, -2), 0.5), c(-4, -3.5, -3, -2.5, -2))
    expect_identical(
        .hsd_lattice(c(-3.1, -2.2), 0.5),
        c(-3.1, -3, -2.5, -2.2)
    )
    expect_identical(.hsd_lattice(c(0.1, 0.2), 0.5), c(0.1, 0.2))
    expect_length(.hsd_lattice(c(-10, 2), 0.125), 97L)

    lattice <- c(-3.5, -3, -2.5)
    expect_identical(.hsd_snap(-3.25, lattice), -3.5)
    expect_identical(.hsd_snap(-3.2, lattice), -3)
    expect_identical(.hsd_snap(-10, lattice), -3.5)

    # start = NULL is the node nearest the midpoint; a given start is snapped
    spec <- transform_pvalues_gsd(
        hsd_opt_raw(),
        spending = list(
            sf_hsd_opt(id = "a"),
            sf_hsd_opt(id = "b", start = 0.1),
            sf_hsd_opt(id = "c", bounds = c(-4, -1), step = 0.5)
        ),
        grid_size = 8L
    )
    expect_identical(unname(spec$param$start), c(-1.5, 0, -2.5))
})

test_that("an untagged call is untouched", {
    raw <- hsd_opt_raw()
    out <- transform_pvalues_gsd(
        raw,
        spending = gsDesign::sfLDOF,
        grid_size = 64L
    )
    expect_s3_class(out, "multigrain_pvals_gsd", exact = TRUE)

    info <- matrix(c(0.5, 1), nrow = 3L, ncol = 2L, byrow = TRUE)
    by_hand <- array(1, dim = dim(raw), dimnames = dimnames(raw))
    tables <- vector("list", 3L)
    for (i in 1:3) {
        transformed <- .gsd_transform_hyp(
            matrix(raw[, i, ], nrow = 2000L, ncol = 2L),
            t_row = info[i, ],
            spending = gsDesign::sfLDOF,
            alpha = 0.025,
            grid_size = 64L,
            look_back = FALSE,
            hyp = i
        )
        by_hand[, i, ] <- transformed$values
        tables[[i]] <- transformed$table
    }
    expected <- new_pvals_gsd(
        pvals = by_hand,
        nsim = 2000L,
        m = 3L,
        n_look = 2L,
        alpha = 0.025,
        info_frac = info,
        look_back = rep(FALSE, 3L),
        spending = rep("gsDesign::sfLDOF", 3L),
        tables = tables
    )

    expect_identical(out, expected)
})

test_that("an all-short-circuited spec is the fixed transform", {
    raw <- hsd_opt_raw()[, 1:2, , drop = FALSE]
    info <- matrix(c(0.5, 1, NA, 1), nrow = 2L, byrow = TRUE)

    # the message says why there is nothing to optimise
    expect_snapshot(
        out <- transform_pvalues_gsd(
            raw,
            info_frac = info,
            spending = list(
                gsDesign::sfLDOF,
                sf_hsd_opt(bounds = c(-4, -2), step = 0.5)
            ),
            grid_size = 64L
        )
    )
    fixed <- transform_pvalues_gsd(
        raw,
        info_frac = info,
        spending = list(gsDesign::sfLDOF, hsd_fixed_fn(-3)),
        grid_size = 64L
    )

    expect_s3_class(out, "multigrain_pvals_gsd", exact = TRUE)
    expect_identical(out$spending, c("gsDesign::sfLDOF", "sfHSD(param = -3)"))
    out$spending <- NULL
    fixed$spending <- NULL
    expect_identical(out, fixed)
})

test_that("node tables are shared across hypotheses on one schedule", {
    spec <- hsd_opt_spec()
    lattice <- spec$param$lattice$sec

    expect_identical(lattice, c(-4, -3.5, -3, -2.5, -2))
    expect_identical(spec$group, c(NA, 1L, 1L))
    expect_identical(spec$schedule_key[[2L]], spec$schedule_key[[3L]])
    expect_true(is.na(spec$schedule_key[[1L]]))
    expect_length(spec$node_tables, length(lattice))

    for (node in lattice) {
        key <- .hsd_table_key(spec$schedule_key[[2L]], node)
        for (i in 2:3) {
            expect_identical(
                spec$node_tables[[key]],
                .gsd_boundary_table(
                    c(0.5, 1),
                    spending = hsd_fixed_fn(node),
                    alpha = 0.025,
                    grid_size = 64L,
                    hyp = i
                )
            )
        }
    }
})

test_that("the spec holds the raw draws under `raw`, not `pvals`", {
    spec <- hsd_opt_spec()

    expect_s3_class(spec, "multigrain_pvals_gsd_spec", exact = TRUE)
    expect_null(spec$pvals)
    expect_identical(spec$raw, hsd_opt_raw())
    expect_identical(spec$nsim, 2000L)
    expect_identical(spec$m, 3L)
    expect_identical(spec$K, 2L)
    expect_identical(spec$grid_size, 64L)
    expect_identical(spec$short_circuit, rep(FALSE, 3L))
    expect_null(spec$fixed[[2L]])
    expect_null(spec$fixed_tables[[3L]])
    expect_identical(dim(spec$fixed[[1L]]), c(2000L, 2L))
    expect_identical(
        spec$spending,
        c(
            "gsDesign::sfLDOF",
            "sfHSD(param = sec in [-4, -2])",
            "sfHSD(param = sec in [-4, -2])"
        )
    )
})

test_that("building a spec draws no random numbers", {
    raw <- hsd_opt_raw()
    withr::with_seed(1, {
        before <- .Random.seed
        transform_pvalues_gsd(
            raw,
            spending = list(
                gsDesign::sfLDOF,
                sf_hsd_opt(bounds = c(-4, -2), step = 0.5),
                sf_hsd_opt(bounds = c(-4, -2), step = 0.5)
            ),
            grid_size = 16L
        )
        expect_identical(.Random.seed, before)
    })
})

test_that("resolve_pvalues_gsd() is the fixed transform at every node", {
    raw <- hsd_opt_raw()
    spec <- hsd_opt_spec()

    for (node in spec$param$lattice$sec) {
        resolved <- resolve_pvalues_gsd(spec, node)
        fixed <- transform_pvalues_gsd(
            raw,
            spending = list(
                gsDesign::sfLDOF,
                hsd_fixed_fn(node),
                hsd_fixed_fn(node)
            ),
            grid_size = 64L
        )

        expect_s3_class(resolved, "multigrain_pvals_gsd")
        for (field in c("pvals", "tables", "info_frac", "look_back", "alpha")) {
            expect_identical(resolved[[field]], fixed[[field]])
        }
        for (field in c("nsim", "m", "K")) {
            expect_identical(resolved[[field]], fixed[[field]])
        }
        label <- sprintf("sfHSD(param = %s)", format(node))
        expect_identical(
            resolved$spending,
            c("gsDesign::sfLDOF", label, label)
        )
    }

    expect_identical(
        resolve_pvalues_gsd(spec, -3.5)$spending[[2L]],
        "sfHSD(param = -3.5)"
    )
    # a name is accepted, and a value within 1e-12 of a node is that node
    expect_identical(
        resolve_pvalues_gsd(spec, c(sec = -3 + 1e-13)),
        resolve_pvalues_gsd(spec, -3)
    )
})

test_that("resolve_pvalues_gsd() matches the fixed path for two parameters", {
    raw <- hsd_opt_raw_d2()
    info <- hsd_opt_info_frac_d2()
    spec <- hsd_opt_spec_d2()

    expect_identical(spec$group, c(NA, 1L, 2L, NA))
    expect_identical(spec$short_circuit, c(FALSE, FALSE, FALSE, TRUE))
    expect_identical(spec$param$names, c("gamma1", "b"))
    expect_identical(spec$spending[[4L]], "sfHSD(param = -2.5)")
    expect_identical(
        spec$param$lattice,
        list(gamma1 = c(-3, -2.5, -2), b = c(-3, -2.5, -2))
    )
    # H2 and H3 share a schedule, so the two parameters share tables
    expect_length(spec$node_tables, 3L)

    for (g1 in spec$param$lattice$gamma1) {
        for (g2 in spec$param$lattice$b) {
            resolved <- resolve_pvalues_gsd(spec, c(b = g2, gamma1 = g1))
            fixed <- transform_pvalues_gsd(
                raw,
                info_frac = info,
                spending = list(
                    gsDesign::sfLDOF,
                    hsd_fixed_fn(g1),
                    hsd_fixed_fn(g2),
                    hsd_fixed_fn(-2.5)
                ),
                grid_size = 64L
            )
            resolved$spending <- NULL
            fixed$spending <- NULL
            expect_identical(resolved, fixed)
        }
    }

    expect_identical(
        resolve_pvalues_gsd(spec, c(-3, -2))$spending,
        c(
            "gsDesign::sfLDOF",
            "sfHSD(param = -3)",
            "sfHSD(param = -2)",
            "sfHSD(param = -2.5)"
        )
    )
})

test_that("fixed hypotheses are labelled by their list argument or name", {
    raw <- hsd_opt_raw()
    tag <- sf_hsd_opt(bounds = c(-4, -2), step = 0.5)

    named <- transform_pvalues_gsd(
        raw,
        spending = list(primary = gsDesign::sfLDOF, a = tag, b = tag),
        grid_size = 8L
    )
    expect_identical(named$spending[[1L]], "primary")

    spending <- list(gsDesign::sfLDPocock, tag, tag)
    by_value <- transform_pvalues_gsd(raw, spending = spending, grid_size = 8L)
    expect_identical(by_value$spending[[1L]], "spending[[1]]")
})

test_that("off-lattice and malformed parameter values are refused", {
    spec <- hsd_opt_spec()

    expect_snapshot(error = TRUE, resolve_pvalues_gsd(spec, -3.2))
    expect_snapshot(error = TRUE, resolve_pvalues_gsd(spec, -5))
    expect_snapshot(error = TRUE, resolve_pvalues_gsd(spec, c(-3, -2)))
    expect_snapshot(error = TRUE, resolve_pvalues_gsd(spec, c(prim = -3)))
    expect_snapshot(error = TRUE, resolve_pvalues_gsd(spec, NA_real_))
    expect_snapshot(error = TRUE, resolve_pvalues_gsd(hsd_opt_raw(), -3))
})

test_that("a node that is not well ordered aborts naming it", {
    # a test double: below level 0.005 it spends like Pocock and above like
    # HSD at node -2, which makes the first boundary fall as the level rises
    not_well_ordered <- structure(
        function(alpha, t, param) {
            if (param > -2.5 && alpha < 0.005) {
                gsDesign::sfLDPocock(alpha, t)$spend
            } else {
                gsDesign::sfHSD(alpha, t, param)$spend
            }
        },
        class = c("multigrain_sf_hsd_opt", "function"),
        id = "sec",
        bounds = c(-4, -2),
        step = 0.5
    )

    expect_snapshot(
        error = TRUE,
        transform_pvalues_gsd(
            hsd_opt_raw(),
            spending = list(
                gsDesign::sfLDOF,
                gsDesign::sfLDOF,
                not_well_ordered
            ),
            grid_size = 64L
        )
    )
})

test_that("a table without two distinct positive boundaries is refused", {
    tab <- list(
        grid = c(1e-3, 1e-2, 2.5e-2),
        bounds = matrix(c(1e-4, 1e-3, 1e-2, 0, 0, 1e-2), nrow = 3L)
    )

    expect_snapshot(
        error = TRUE,
        .hsd_check_invertible(tab, hyp = 2L, node = -3, name = "sec")
    )
    expect_true(.hsd_check_invertible(
        list(grid = tab$grid, bounds = tab$bounds[, c(1L, 1L)]),
        hyp = 2L,
        node = -3
    ))
})

test_that("a tagged function cannot reach the fixed transform", {
    expect_error(
        .gsd_spend(sf_hsd_opt(), alpha = 0.025, t_look = 1, hyp = 1L),
        "bug in multigrain"
    )
    node_fn <- .hsd_node_fn(sf_hsd_opt(), -2)
    expect_identical(
        .gsd_spend(node_fn, alpha = 0.025, t_look = c(0.5, 1), hyp = 1L),
        gsDesign::sfHSD(0.025, c(0.5, 1), -2)$spend
    )
})

test_that("check_pvals_gsd() refuses a spec with a hint", {
    spec <- hsd_opt_spec()

    expect_snapshot(error = TRUE, check_pvals_gsd(spec))
})

test_that("the spec prints", {
    expect_snapshot(print(hsd_opt_spec()))
    expect_snapshot(print(hsd_opt_spec_d2()))
})
