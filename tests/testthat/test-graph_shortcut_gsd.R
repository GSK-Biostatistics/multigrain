# Gate for phase P2 of `dev/gsd_design_record.md`. `graph_shortcut_gsd()` is
# the fixed-sample cascade of `graph_shortcut()` wrapped in a loop over
# analyses, run on the nsim by m*K matrix of transformed p-values that
# `transform_pvalues_gsd()` produces. The reference implementation it is
# checked against lives in `helper-gsd_reference.R`.

# LDOF cumulative spend as a plain numeric vector: the form both our transform
# and `graphicalMCP` accept, so the oracle tests use one spending function.
gsd_ldof <- function(a, t) gsDesign::sfLDOF(a, t)$spend

# The Appendix B simulation, reproduced exactly: three hypotheses, three
# analyses, H2 maturing at analysis 2 and H3 at analysis 1.
gsd_appendix_b <- function() {
    tmat <- rbind(c(1 / 3, 2 / 3, 1), c(0.5, 1, 1), c(1, 1, 1))
    sfs <- list(
        gsDesign::sfLDOF,
        function(a, t) gsDesign::sfHSD(a, t, param = -2),
        gsDesign::sfLDOF
    )
    n_sim <- 400L
    m <- 3L
    n_look <- 3L
    delta <- c(2.6, 2.4, 2.9)
    mature <- apply(tmat, 1L, function(t) which(t >= 1)[[1L]])

    set.seed(11)
    p_list <- lapply(seq_len(m), function(i) {
        kk <- mature[[i]]
        t_i <- tmat[i, seq_len(kk)]
        corr <- outer(
            t_i,
            t_i,
            function(a, b) sqrt(pmin(a, b) / pmax(a, b))
        )
        z <- mvtnorm::rmvnorm(
            n_sim,
            mean = delta[[i]] * sqrt(t_i),
            sigma = corr
        )
        p <- stats::pnorm(z, lower.tail = FALSE)
        if (kk < n_look) {
            p <- cbind(p, matrix(p[, kk], nrow = n_sim, ncol = n_look - kk))
        }
        p
    })

    raw_pvals <- array(NA_real_, dim = c(n_sim, m, n_look))
    for (i in seq_len(m)) {
        raw_pvals[, i, ] <- p_list[[i]]
    }

    list(
        p_list = p_list,
        raw_pvals = raw_pvals,
        tmat = tmat,
        sfs = sfs,
        w0 = c(0.5, 0.3, 0.2),
        G0 = rbind(c(0, 0.5, 0.5), c(0.5, 0, 0.5), c(0.5, 0.5, 0)),
        n_sim = n_sim,
        m = m,
        n_look = n_look
    )
}

# A p-value matrix with the ordinary fixture rows, rows well below 1e-12, and
# a row of exact zeros.
gsd_extreme_pvals <- function() {
    tiny <- withr::with_seed(
        31,
        matrix(10^stats::runif(40L * 4L, -20, -12), nrow = 40L, ncol = 4L)
    )

    # nolint start: object_usage_linter. `pvals_fixture` is a test helper.
    rbind(
        pvals_fixture[seq_len(2000L), seq_len(4L)],
        tiny,
        matrix(0, nrow = 1L, ncol = 4L)
    )
    # nolint end
}


### K = 1 IDENTITY ###

test_that("K = 1 is identical to graph_shortcut(), including p below 1e-12", {
    pvals <- gsd_extreme_pvals()
    alpha <- 0.025

    for (seed in 1:4) {
        gr <- make_test_graph(4L, seed = seed)
        reference <- graph_shortcut(pvals, alpha, gr$w, gr$G)
        out <- graph_shortcut_gsd(pvals, alpha, gr$w, gr$G, K = 1L)

        expect_identical(out$rejected, reference)
        expect_identical(
            out$time,
            matrix(
                as.integer(reference),
                nrow = nrow(pvals),
                ncol = ncol(pvals)
            )
        )
    }
})

test_that("the transform then the kernel is identical at K = 1", {
    pvals <- gsd_extreme_pvals()
    alpha <- 0.025

    transformed <- transform_pvalues_gsd(
        array(pvals, dim = c(nrow(pvals), 4L, 1L)),
        info_frac = 1,
        spending = gsDesign::sfLDOF,
        alpha = alpha,
        grid_size = 8L
    )
    reshaped <- gsd_pvals_matrix(transformed)

    # a single full-information analysis short-circuits, so the transform is
    # the identity and the reshape is a no-op
    expect_identical(reshaped, pvals)

    for (seed in 1:4) {
        gr <- make_test_graph(4L, seed = seed)
        expect_identical(
            graph_shortcut_gsd(reshaped, alpha, gr$w, gr$G, K = 1L)$rejected,
            graph_shortcut(pvals, alpha, gr$w, gr$G)
        )
    }
})


### R REFERENCE EQUIVALENCE ###

test_that("the kernel matches the boundary-on-the-fly reference", {
    setup <- gsd_appendix_b()
    alpha <- 0.025
    cache <- new.env(parent = emptyenv())

    configs <- list(
        c(FALSE, FALSE, FALSE),
        c(TRUE, TRUE, TRUE),
        c(TRUE, FALSE, TRUE),
        c(FALSE, TRUE, FALSE)
    )

    for (look_back in configs) {
        label <- paste(substr(as.character(look_back), 1L, 1L), collapse = "")

        transformed <- transform_pvalues_gsd(
            setup$raw_pvals,
            info_frac = setup$tmat,
            spending = setup$sfs,
            alpha = alpha,
            look_back = look_back,
            grid_size = 1024L
        )
        out <- graph_shortcut_gsd(
            gsd_pvals_matrix(transformed),
            alpha,
            setup$w0,
            setup$G0,
            K = setup$n_look
        )
        reference <- gsd_reference_direct(
            setup$p_list,
            setup$tmat,
            setup$sfs,
            setup$w0,
            setup$G0,
            alpha = alpha,
            look_back = look_back,
            cache = cache
        )

        expect_identical(out$rejected, reference$rejected, info = label)
        expect_identical(out$time, reference$time, info = label)
    }
})


### PARALLEL EQUALS SERIAL ###

test_that("the parallel kernel is identical to the serial one", {
    setup <- gsd_appendix_b()
    alpha <- 0.025

    transformed <- transform_pvalues_gsd(
        setup$raw_pvals,
        info_frac = setup$tmat,
        spending = setup$sfs,
        alpha = alpha,
        look_back = c(TRUE, FALSE, TRUE),
        grid_size = 1024L
    )
    reshaped <- gsd_pvals_matrix(transformed)
    serial <- graph_shortcut_gsd(
        reshaped,
        alpha,
        setup$w0,
        setup$G0,
        K = setup$n_look
    )

    for (threads in c(1L, 2L, 4L, 8L)) {
        expect_identical(
            graph_shortcut_gsd_parallel(
                reshaped,
                alpha,
                setup$w0,
                setup$G0,
                K = setup$n_look,
                num_threads = threads
            ),
            serial,
            info = paste("num_threads =", threads)
        )
    }

    expect_identical(
        graph_shortcut_gsd_parallel(
            reshaped,
            alpha,
            setup$w0,
            setup$G0,
            K = setup$n_look,
            num_threads = 4L,
            grain_size = 7L
        ),
        serial
    )
})


### ZERO WEIGHT ###

test_that("a zero-weight hypothesis never rejects without recycled alpha", {
    # H3 carries no weight and no edge points into it, so its allocation is 0
    # at every analysis even though its transformed p-value sits on the floor
    w <- c(0.5, 0.5, 0)
    G <- rbind(c(0, 1, 0), c(1, 0, 0), c(0, 0, 0))
    alpha <- 0.025

    pvals <- matrix(1, nrow = 1L, ncol = 9L)
    pvals[1L, c(1L, 5L)] <- 1e-6 # H1 at analysis 1, H2 at analysis 2
    pvals[1L, c(3L, 6L, 9L)] <- 1e-14 # H3 on the floor at every analysis

    out <- graph_shortcut_gsd(pvals, alpha, w, G, K = 3L)

    expect_identical(out$rejected[1L, ], c(TRUE, TRUE, FALSE))
    expect_identical(out$time[1L, ], c(1L, 2L, 0L))
})

test_that("a zero-weight hypothesis rejects exactly when alpha reaches it", {
    # alpha reaches H3 only once H2 is rejected, which happens at analysis 2
    w <- c(0.5, 0.5, 0)
    G <- rbind(c(0, 1, 0), c(0, 0, 1), c(0, 0, 0))
    alpha <- 0.025

    pvals <- matrix(1, nrow = 1L, ncol = 9L)
    pvals[1L, 1L] <- 1e-6 # H1 rejects at analysis 1
    pvals[1L, 2L] <- 0.03 # H2 not yet
    pvals[1L, 5L] <- 1e-6 # H2 rejects at analysis 2
    pvals[1L, c(3L, 6L, 9L)] <- 1e-14 # H3 on the floor throughout

    out <- graph_shortcut_gsd(pvals, alpha, w, G, K = 3L)

    expect_identical(out$rejected[1L, ], c(TRUE, TRUE, TRUE))
    expect_identical(out$time[1L, ], c(1L, 2L, 2L))
})


### CLAMP SAFETY ###

test_that("p = 5e-11 does not reject on a propagated allocation of 2.5e-11", {
    # H1 has weight 1e-4 and matures at analysis 1, so its local level is
    # 1e-4 * alpha = 2.5e-6 and its transformed p-value is its raw one. The
    # single edge out of H1 carries 1e-5 of that to H2, which therefore holds
    # 1e-4 * 1e-5 * alpha = 2.5e-11 once H1 is rejected.
    alpha <- 0.025
    w <- c(1e-4, 0)
    G <- rbind(c(0, 1e-5), c(0, 0))

    raw_pvals <- array(1, dim = c(3L, 2L, 3L))
    raw_pvals[, 1L, ] <- 1e-12 # H1 rejects at analysis 1
    raw_pvals[1L, 2L, 3L] <- 5e-11 # above the allocation: must not reject
    raw_pvals[2L, 2L, 3L] <- 1e-11 # below it: must reject
    raw_pvals[3L, 2L, 3L] <- 1e-100 # floored at 1e-14, still below it

    transformed <- transform_pvalues_gsd(
        raw_pvals,
        info_frac = rbind(c(1, 1, 1), c(1 / 3, 2 / 3, 1)),
        spending = gsDesign::sfLDOF,
        alpha = alpha,
        grid_size = 1024L
    )
    repeated <- transformed$pvals[, 2L, 3L]
    out <- graph_shortcut_gsd(
        gsd_pvals_matrix(transformed),
        alpha,
        w,
        G,
        K = 3L
    )

    # the transformed value for the red-team row is the raw one to within
    # interpolation error, and it exceeds the allocation
    expect_lt(abs(repeated[[1L]] - 5e-11) / 5e-11, 1e-5)
    expect_gt(repeated[[1L]], 2.5e-11)
    expect_lt(repeated[[2L]], 2.5e-11)
    expect_identical(repeated[[3L]], 1e-14)

    # H1 is rejected in every row, so the allocation really is the propagated
    # one; rows 2 and 3 bracket it from below and row 1 from above
    expect_true(all(out$rejected[, 1L]))
    expect_identical(out$rejected[, 2L], c(FALSE, TRUE, TRUE))
    expect_identical(out$time[, 2L], c(0L, 3L, 3L))
})


### ORACLE: MAURER AND BRETZ CASE STUDY ###

test_that("the kernel matches graphicalMCP on the Maurer-Bretz case study", {
    skip_if_not_installed("graphicalMCP")

    # Figure 1 graph and the Table 2 raw p-values at analyses 1 and 2. The
    # published case study has two analyses; analysis 3 is padded with a raw
    # p-value of 1 on both sides so that the oracle and the kernel see the
    # same K = 3 problem. A p-value of 1 lies above every boundary, so its
    # repeated p-value is 1 and it cannot reject.
    #
    # Argument mapping for `graphicalMCP::graph_test_shortcut_gsd()`:
    #   `graph` is the graph_create() object built from w and G;
    #   `p` is the m by K matrix of RAW p-values, with rows named as the
    #     graph's hypotheses and no column names, because its validator
    #     compares `is.na(p)` with `is.na(info_frac)` using `identical()`,
    #     which also compares dimnames;
    #   `info_frac` is a length-K vector, recycled over hypotheses;
    #   `spending_fn` is one function returning a PLAIN numeric vector of
    #     cumulative spend, since `gsDesign::sfLDOF` itself returns an object
    #     and silently yields NA boundaries;
    #   `look_back` is FALSE, matching our transform's `look_back = FALSE`.
    # The oracle compares with `<=` where the kernel uses `<`; the case-study
    # p-values are far from every boundary, so the two agree.
    w <- c(0.5, 0.5, 0, 0)
    G <- rbind(
        c(0, 0.5, 0.5, 0),
        c(0.5, 0, 0, 0.5),
        c(0, 1, 0, 0),
        c(1, 0, 0, 0)
    )
    alpha <- 0.025
    info <- c(1 / 3, 2 / 3, 1)
    p1 <- c(0.0062, 0.017, 0.009, 0.13)
    p2 <- c(0.0002, 0.0035, 0.002, 0.06)

    raw_pvals <- array(1, dim = c(1L, 4L, 3L))
    raw_pvals[1L, , 1L] <- p1
    raw_pvals[1L, , 2L] <- p2

    transformed <- transform_pvalues_gsd(
        raw_pvals,
        info_frac = info,
        spending = gsd_ldof,
        alpha = alpha,
        grid_size = 1024L
    )
    out <- graph_shortcut_gsd(
        gsd_pvals_matrix(transformed),
        alpha,
        w,
        G,
        K = 3L
    )

    graph <- graphicalMCP::graph_create(w, G)
    oracle_p <- cbind(p1, p2, rep(1, 4L))
    dimnames(oracle_p) <- list(names(graph$hypotheses), NULL)
    oracle <- graphicalMCP::graph_test_shortcut_gsd(
        graph = graph,
        p = oracle_p,
        alpha = alpha,
        info_frac = info,
        spending_fn = gsd_ldof,
        look_back = FALSE
    )

    # nothing crosses at analysis 1 (the Table 2 values all exceed 0.0125);
    # H1, H2 and H3 are rejected at analysis 2 and H4 is retained
    expect_identical(out$rejected[1L, ], c(TRUE, TRUE, TRUE, FALSE))
    expect_identical(out$time[1L, ], c(2L, 2L, 2L, 0L))

    rejected <- unname(oracle$outputs$rejected)
    expect_identical(out$rejected[1L, ], rejected)
    # `decision_at` is also set for hypotheses that were merely tested, so it
    # is compared only where the oracle rejected
    expect_identical(
        out$time[1L, rejected],
        as.integer(unname(oracle$outputs$decision_at)[rejected])
    )

    # our repeated p-values agree with the oracle's, except that our boundary
    # tables stop at alpha: a repeated p-value the oracle puts above alpha is
    # "not rejectable at any allocation" and comes back as 1 (section 4.1)
    oracle_repeated <- unname(oracle$outputs$repeated_p[, 1:2])
    ours <- transformed$pvals[1L, , 1:2]
    above <- oracle_repeated > alpha

    expect_identical(ours[above], rep(1, sum(above)))
    expect_lt(
        max(
            abs(ours[!above] - oracle_repeated[!above]) /
                oracle_repeated[!above]
        ),
        1e-4
    )
})


### ORACLE: EXAMPLE 5 ###

test_that("a matured hypothesis is rejected when alpha is recycled to it", {
    skip_if_not_installed("graphicalMCP")

    # Example 5 of the manuscript: PFS (H1) matures at analysis 1 and holds no
    # weight; OS (H2) is analysed at half information and then at full
    # information, and recycles its whole level to PFS when it is rejected.
    #
    # `graphicalMCP` pads a matured hypothesis's later analyses with NA and,
    # with `look_back = FALSE`, never re-tests it there; only
    # `look_back = c(TRUE, FALSE)` carries its sequential p-value forward and
    # so reproduces the manuscript's convention (design record, section 4.9).
    # Our transform copies the matured repeated p-value forward instead, so
    # the same configuration is `look_back = FALSE` on our side.
    alpha <- 0.025
    w <- c(0, 1)
    G <- rbind(c(0, 1), c(1, 0))
    info <- rbind(c(1, NA), c(0.5, 1))
    p_pfs <- 0.02
    p_os <- c(0.5, 0.005)

    raw_pvals <- array(NA_real_, dim = c(1L, 2L, 2L))
    raw_pvals[1L, 1L, 1L] <- p_pfs
    raw_pvals[1L, 2L, ] <- p_os

    transformed <- transform_pvalues_gsd(
        raw_pvals,
        info_frac = info,
        spending = gsd_ldof,
        alpha = alpha,
        grid_size = 1024L
    )
    out <- graph_shortcut_gsd(
        gsd_pvals_matrix(transformed),
        alpha,
        w,
        G,
        K = 2L
    )

    expect_identical(out$rejected[1L, ], c(TRUE, TRUE))
    expect_identical(out$time[1L, ], c(2L, 2L))

    graph <- graphicalMCP::graph_create(w, G)
    oracle_p <- rbind(c(p_pfs, NA), p_os)
    dimnames(oracle_p) <- list(names(graph$hypotheses), NULL)
    oracle_info <- info
    dimnames(oracle_info) <- list(names(graph$hypotheses), NULL)

    matched <- graphicalMCP::graph_test_shortcut_gsd(
        graph = graph,
        p = oracle_p,
        alpha = alpha,
        info_frac = oracle_info,
        spending_fn = gsd_ldof,
        look_back = c(TRUE, FALSE)
    )
    expect_identical(unname(matched$outputs$rejected), c(TRUE, TRUE))
    expect_identical(
        as.integer(unname(matched$outputs$decision_at)),
        c(2L, 2L)
    )

    # the convention matters: with look_back = FALSE the oracle treats the
    # matured hypothesis as having no data at analysis 2 and never rejects it
    unmatched <- graphicalMCP::graph_test_shortcut_gsd(
        graph = graph,
        p = oracle_p,
        alpha = alpha,
        info_frac = oracle_info,
        spending_fn = gsd_ldof,
        look_back = c(FALSE, FALSE)
    )
    expect_identical(unname(unmatched$outputs$rejected), c(FALSE, TRUE))
})

test_that("a matured hypothesis rejected outright keeps decision time 1", {
    skip_if_not_installed("graphicalMCP")

    alpha <- 0.025
    w <- c(0.5, 0.5)
    G <- rbind(c(0, 1), c(1, 0))
    info <- rbind(c(1, NA), c(0.5, 1))
    p_pfs <- 0.001
    p_os <- c(0.5, 0.005)

    raw_pvals <- array(NA_real_, dim = c(1L, 2L, 2L))
    raw_pvals[1L, 1L, 1L] <- p_pfs
    raw_pvals[1L, 2L, ] <- p_os

    transformed <- transform_pvalues_gsd(
        raw_pvals,
        info_frac = info,
        spending = gsd_ldof,
        alpha = alpha,
        grid_size = 1024L
    )
    out <- graph_shortcut_gsd(
        gsd_pvals_matrix(transformed),
        alpha,
        w,
        G,
        K = 2L
    )

    expect_identical(out$rejected[1L, ], c(TRUE, TRUE))
    expect_identical(out$time[1L, ], c(1L, 2L))

    graph <- graphicalMCP::graph_create(w, G)
    oracle_p <- rbind(c(p_pfs, NA), p_os)
    dimnames(oracle_p) <- list(names(graph$hypotheses), NULL)
    oracle_info <- info
    dimnames(oracle_info) <- list(names(graph$hypotheses), NULL)

    oracle <- graphicalMCP::graph_test_shortcut_gsd(
        graph = graph,
        p = oracle_p,
        alpha = alpha,
        info_frac = oracle_info,
        spending_fn = gsd_ldof,
        look_back = c(TRUE, FALSE)
    )
    expect_identical(unname(oracle$outputs$rejected), c(TRUE, TRUE))
    expect_identical(as.integer(unname(oracle$outputs$decision_at)), c(1L, 2L))
})


### DETERMINISM ###

test_that("repeated calls give identical results", {
    setup <- gsd_appendix_b()
    alpha <- 0.025
    transformed <- transform_pvalues_gsd(
        setup$raw_pvals,
        info_frac = setup$tmat,
        spending = setup$sfs,
        alpha = alpha,
        grid_size = 256L
    )
    reshaped <- gsd_pvals_matrix(transformed)

    expect_identical(
        graph_shortcut_gsd(reshaped, alpha, setup$w0, setup$G0, K = 3L),
        graph_shortcut_gsd(reshaped, alpha, setup$w0, setup$G0, K = 3L)
    )
    expect_identical(
        graph_shortcut_gsd_parallel(
            reshaped,
            alpha,
            setup$w0,
            setup$G0,
            K = 3L,
            num_threads = 4L
        ),
        graph_shortcut_gsd_parallel(
            reshaped,
            alpha,
            setup$w0,
            setup$G0,
            K = 3L,
            num_threads = 4L
        )
    )
})


### INPUT VALIDATION ###

test_that("the kernel validates its inputs", {
    pvals <- matrix(stats::runif(60L), nrow = 10L, ncol = 6L)
    w <- c(0.5, 0.5, 0)
    G <- matrix(0, nrow = 3L, ncol = 3L)

    expect_error(
        graph_shortcut_gsd(pvals, 0.025, w, G, K = 4L),
        "multiple of K"
    )
    expect_error(graph_shortcut_gsd(pvals, 0.025, w, G, K = 0L), "K must be")
    expect_error(
        graph_shortcut_gsd(pvals, 0.025, w, matrix(0, 4L, 4L), K = 2L),
        "m x m"
    )
    expect_error(
        graph_shortcut_gsd(pvals, 0.025, c(0.5, 0.5), G, K = 2L),
        "length"
    )

    expect_error(
        graph_shortcut_gsd_parallel(pvals, 0.025, w, G, K = 4L),
        "multiple of K"
    )
    expect_error(
        graph_shortcut_gsd_parallel(pvals, 0.025, w, G, K = 0L),
        "K must be"
    )
    expect_error(
        graph_shortcut_gsd_parallel(pvals, 0.025, w, matrix(0, 4L, 4L), K = 2L),
        "m x m"
    )
    expect_error(
        graph_shortcut_gsd_parallel(pvals, 0.025, c(0.5, 0.5), G, K = 2L),
        "length"
    )
    expect_error(
        graph_shortcut_gsd_parallel(
            pvals,
            0.025,
            w,
            G,
            K = 2L,
            num_threads = 0L
        ),
        "num_threads must be >= 1"
    )
    expect_error(
        graph_shortcut_gsd_parallel(
            pvals,
            0.025,
            w,
            G,
            K = 2L,
            num_threads = 1L,
            grain_size = 0L
        ),
        "grain_size must be >= 1"
    )
})
