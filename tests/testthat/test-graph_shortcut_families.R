# Kernels for families of hypotheses (Maurer and Bretz 2014): the within-family
# adjustment `adjust_pvals_families()` and the shortcut
# `graph_shortcut_families()` / `graph_shortcut_families_parallel()`.
# Helpers (oracle, random partitions and graphs) are in helper-families.R.

### INPUTS ###

test_that("adjust_pvals_families validates family and test codes", {
    pvals <- matrix(stats::runif(30), nrow = 10, ncol = 3)

    expect_error(
        adjust_pvals_families(pvals, family = c(0L, 1L), test = c(0L, 0L)),
        "family must have length"
    )
    expect_error(
        adjust_pvals_families(pvals, family = c(0L, 1L, 2L), test = c(0L, 0L)),
        "family must take values"
    )
    expect_error(
        adjust_pvals_families(pvals, family = c(0L, 0L, 1L), test = c(0L, 3L)),
        "test must take values"
    )
})

test_that("graph_shortcut_families validates its inputs", {
    pvals <- matrix(stats::runif(30), nrow = 10, ncol = 3)
    w <- c(0.5, 0.5)
    G <- matrix(c(0, 1, 1, 0), nrow = 2)

    expect_error(
        graph_shortcut_families(pvals, c(0L, 0L, 1L), 0.025, w, diag(3)),
        "M x M"
    )
    expect_error(
        graph_shortcut_families(pvals, c(0L, 1L), 0.025, w, G),
        "family must have length"
    )
    expect_error(
        graph_shortcut_families(pvals, c(0L, 1L, 2L), 0.025, w, G),
        "family must take values"
    )
    expect_error(
        graph_shortcut_families(pvals, c(0L, 0L, 0L), 0.025, w, G),
        "at least one hypothesis"
    )
    expect_error(
        graph_shortcut_families_parallel(
            pvals,
            c(0L, 0L, 1L),
            0.025,
            w,
            G,
            num_threads = 0
        ),
        "num_threads must be >= 1"
    )
    expect_error(
        graph_shortcut_families_parallel(
            pvals,
            c(0L, 0L, 1L),
            0.025,
            w,
            G,
            num_threads = 1,
            grain_size = 0
        ),
        "grain_size must be >= 1"
    )
})

### WITHIN-FAMILY ADJUSTMENT ###

test_that("adjust_pvals_families equals stats::p.adjust (k = 1..4)", {
    withr::local_seed(1)

    for (test in c("hommel", "hochberg", "holm")) {
        for (k in 1:4) {
            pvals <- matrix(stats::runif(2000 * k)^2, ncol = k)
            expected <- t(apply(pvals, 1, stats::p.adjust, method = test))
            if (k == 1) {
                expected <- t(expected)
            }

            expect_identical(
                adjust_families_r(pvals, list(seq_len(k)), test),
                unname(expected),
                info = paste(test, "k =", k)
            )
        }
    }
})

test_that("adjust_pvals_families adjusts each family on its own", {
    withr::local_seed(2)
    pvals <- matrix(stats::runif(500 * 4)^2, ncol = 4)
    members <- list(c(1, 3), 4, 2)
    test <- c("holm", "hommel", "hochberg")

    expected <- pvals
    for (f in seq_along(members)) {
        idx <- members[[f]]
        expected[, idx] <- t(apply(
            pvals[, idx, drop = FALSE],
            1,
            stats::p.adjust,
            method = test[f]
        ))
    }

    expect_identical(adjust_families_r(pvals, members, test), expected)
})

### OUTPUTS ###

test_that("graph_shortcut_families returns a LogicalMatrix with correct dims", {
    out <- graph_shortcut_families(
        pvals_fixture[1:10, 1:4],
        family = c(0L, 1L, 1L, 2L),
        alpha = 0.025,
        w = c(1, 0, 0),
        G = matrix(c(0, 0, 1, 1, 0, 0, 0, 1, 0), nrow = 3)
    )
    expect_true(is.matrix(out))
    expect_type(out, "logical")
    expect_identical(dim(out), c(10L, 4L))
})

### EQUIVALENCE ###

test_that("all-singleton families are identical to graph_shortcut", {
    for (m in 2:4) {
        gr <- make_test_graph(m, seed = 300 + m)
        pvals <- pvals_fixture[1:5000, seq_len(m)]

        expect_identical(
            graph_shortcut_families(pvals, 0:(m - 1), 0.025, gr$w, gr$G),
            graph_shortcut(pvals, 0.025, gr$w, gr$G),
            info = paste("m =", m)
        )
    }
})

test_that("graph_shortcut_families matches literal Algorithm 5.7", {
    withr::local_seed(57)

    for (test in c("hommel", "hochberg", "holm")) {
        for (case in 1:10) {
            m <- sample(2:4, 1)
            num_fam <- sample.int(m, 1)
            members <- random_partition(m, num_fam)
            gr <- random_family_graph(num_fam, full_rows = case %% 2 == 0)
            pvals <- matrix(stats::runif(300 * m)^4, ncol = m)

            out <- shortcut_families_r(pvals, members, test, gr$w, gr$G)
            for (order in c("lowest", "random")) {
                oracle <- alg57_oracle(
                    pvals,
                    members,
                    test,
                    gr$w,
                    gr$G,
                    order = order
                )
                expect_identical(
                    unname(out),
                    oracle,
                    info = paste(test, "case", case, order)
                )
            }
        }
    }
})

test_that("worked example (families design, section 1.4)", {
    # Peanut primary; peanut secondaries; cashew, milk, egg primaries; egg,
    # milk and cashew secondaries. Hommel in every family.
    members <- list(1, 2:5, 6:8, 9:12, 13:16, 17:20)
    pvals <- matrix(
        c(
            0.004,
            0.003, 0.006, 0.015, 0.20,
            0.002, 0.009, 0.016,
            0.001, 0.002, 0.003, 0.004,
            0.004, 0.03, 0.5, 0.6,
            0.002, 0.003, 0.010, 0.012
        ),
        nrow = 1
    )
    w <- c(1, 0, 0, 0, 0, 0)
    G <- matrix(0, 6, 6)
    G[1, 2:3] <- c(0.25, 0.75)
    G[2, 3] <- 1
    G[3, 4:6] <- 1 / 3
    G[4, c(2, 5, 6)] <- c(0.001, 0.4995, 0.4995)
    G[5, c(4, 6)] <- 0.5
    G[6, c(4, 5)] <- 0.5

    out <- shortcut_families_r(pvals, members, "hommel", w, G)

    expect_identical(which(out), c(1L, 6:12, 17L, 18L))
    expect_identical(
        unname(out),
        alg57_oracle(pvals, members, "hommel", w, G, order = "random")
    )
})

### CASE CHECK ###

test_that("a family passes its weight on only when all members are rejected", {
    # Family 2 = {H2, H3} receives alpha from H1; H4 only from family 2.
    members <- list(1, 2:3, 4)
    w <- c(1, 0, 0)
    G <- matrix(
        c(
            0, 1, 0,
            0, 0, 1,
            1, 0, 0
        ),
        nrow = 3,
        byrow = TRUE
    )
    pvals <- rbind(
        c(0.001, 0.001, 0.010, 0.001), # family 2 emptied -> H4 rejected
        c(0.001, 0.001, 0.030, 0.001), # H2 only -> nothing reaches H4
        c(0.030, 0.001, 0.001, 0.001) # H1 fails -> nothing downstream
    )

    out <- shortcut_families_r(pvals, members, "hommel", w, G)

    expected <- rbind(
        c(TRUE, TRUE, TRUE, TRUE),
        c(TRUE, TRUE, FALSE, FALSE),
        c(FALSE, FALSE, FALSE, FALSE)
    )
    expect_identical(unname(out), expected)
})

test_that("a family that is never reached rejects nothing", {
    members <- list(1:2, 3, 4)
    w <- c(1, 0, 0)
    G <- matrix(
        c(
            0, 1, 0,
            1, 0, 0,
            1, 0, 0
        ),
        nrow = 3,
        byrow = TRUE
    )
    pvals <- matrix(1e-10, nrow = 5, ncol = 4)

    out <- shortcut_families_r(pvals, members, "hommel", w, G)

    expect_true(all(out[, 1:3]))
    expect_false(any(out[, 4]))
})

test_that("rows of G summing below 1 match the oracle", {
    withr::local_seed(11)
    members <- list(c(1, 4), 2, 3)
    w <- c(0.6, 0.4, 0)
    G <- matrix(
        c(
            0, 0.3, 0.2,
            0.1, 0, 0.5,
            0.4, 0.4, 0
        ),
        nrow = 3,
        byrow = TRUE
    )
    pvals <- matrix(stats::runif(2000 * 4)^4, ncol = 4)

    expect_identical(
        unname(shortcut_families_r(pvals, members, "holm", w, G)),
        alg57_oracle(pvals, members, "holm", w, G)
    )
})

test_that("graph_shortcut_families uses strict `<` at the family level", {
    # Family {H1, H2} with Holm: q = 2 * 0.0125 = 0.025 = w * alpha.
    pvals <- matrix(c(0.0125, 0.0125), nrow = 1)
    padj <- adjust_families_r(pvals, list(1:2), "holm")
    expect_identical(padj, matrix(c(0.025, 0.025), nrow = 1))

    out <- graph_shortcut_families(padj, c(0L, 0L), 0.025, 1, matrix(0))
    expect_false(any(out))
})

test_that("read-only inputs are not mutated", {
    pvals <- pvals_fixture[1:500, 1:4]
    fam_index <- c(0L, 1L, 1L, 2L)
    w <- c(0.5, 0.5, 0)
    G <- matrix(
        c(
            0, 0.5, 0.5,
            0.5, 0, 0.5,
            0.5, 0.5, 0
        ),
        nrow = 3,
        byrow = TRUE
    )

    p_before <- pvals
    family_before <- fam_index
    w_before <- w
    G_before <- G

    padj <- adjust_pvals_families(pvals, fam_index, test = c(0L, 0L, 0L))
    padj_before <- padj
    invisible(graph_shortcut_families(padj, fam_index, 0.025, w, G))
    invisible(
        graph_shortcut_families_parallel(
            padj,
            fam_index,
            0.025,
            w,
            G,
            num_threads = 2
        )
    )

    expect_identical(pvals, p_before)
    expect_identical(padj, padj_before)
    expect_identical(fam_index, family_before)
    expect_identical(w, w_before)
    expect_identical(G, G_before)
})

### SERIAL equivalent to PARALLEL ###

test_that("graph_shortcut_families_parallel equals the serial kernel", {
    withr::local_seed(21)
    members <- list(c(1, 3), 2, 4)
    gr <- random_family_graph(3)
    padj <- adjust_families_r(pvals_fixture[1:1e4, 1:4], members, "hommel")
    fam_index <- family_index0(members)

    out_ref <- graph_shortcut_families(padj, fam_index, 0.025, gr$w, gr$G)

    for (nt in c(1L, 2L, 4L)) {
        expect_identical(
            graph_shortcut_families_parallel(
                padj,
                fam_index,
                0.025,
                gr$w,
                gr$G,
                num_threads = nt
            ),
            out_ref,
            info = paste("num_threads =", nt)
        )
    }
    for (gs in c(1L, 500L, 10000L)) {
        expect_identical(
            graph_shortcut_families_parallel(
                padj,
                fam_index,
                0.025,
                gr$w,
                gr$G,
                num_threads = 2,
                grain_size = gs
            ),
            out_ref,
            info = paste("grain_size =", gs)
        )
    }
})
