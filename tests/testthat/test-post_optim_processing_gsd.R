# The greedy pruning of `prune_graph_gsd()` accepts a candidate only when the
# group sequential gain does not fall, so the gain of the pruned graph can
# never be below the gain of the graph it started from (design record,
# section 6 P4).

pp_pvals_gsd <- withr::with_seed(23, {
    raw <- simulate_pvalues_gsd(
        power_nominal = c(0.9, 0.85, 0.8),
        corr_matrix = diag(3),
        info_frac = c(0.5, 1),
        nsim = 1000
    )
    transform_pvalues_gsd(
        raw,
        spending = gsDesign::sfLDOF,
        grid_size = 128L
    )
})

pp_gain <- trial_success_gsd(
    d(t1) + d(t2) + d(t3),
    d = c(1, 0.75),
    verbose = "silent"
)

pp_value <- function(w, trans) {
    calc_power_pvals_gsd(
        pp_pvals_gsd,
        hyp_weight = w,
        trans_matrix = trans,
        custom_power = pp_gain
    )$custom_power
}

# Graphs with epsilon weights and epsilon edges planted, which greedy pruning
# should be able to remove.
pp_graphs <- list(
    list(
        w = c(0.5, 0.4999, 0.0001),
        trans = rbind(
            c(0, 0.999, 0.001),
            c(0.999, 0, 0.001),
            c(0.5, 0.5, 0)
        )
    ),
    list(
        w = c(0.001, 0.4995, 0.4995),
        trans = rbind(
            c(0, 0.5, 0.5),
            c(0.001, 0, 0.999),
            c(0.001, 0.999, 0)
        )
    ),
    list(
        w = c(0.3333, 0.3333, 0.3334),
        trans = rbind(
            c(0, 0.5, 0.5),
            c(0.5, 0, 0.5),
            c(0.5, 0.5, 0)
        )
    )
)


test_that("prune_graph_gsd() never lowers the gain and stays valid", {
    gc <- graph_constraint_free(3)

    for (i in seq_along(pp_graphs)) {
        g <- pp_graphs[[i]]
        before <- pp_value(g$w, g$trans)

        pruned <- prune_graph_gsd(
            pvals = pp_pvals_gsd,
            hyp_weight = g$w,
            trans_matrix = g$trans,
            trial_success = pp_gain,
            graph_constraint = gc,
            gamma = 1,
            verbose = "silent"
        )

        after <- pp_value(pruned$hyp_weight, pruned$trans_matrix)

        expect_gte(after, before)
        expect_true(
            is_graph_valid(pruned$hyp_weight, pruned$trans_matrix),
            label = paste("graph", i)
        )
    }
})


test_that("prune_graph_gsd() leaves fixed weights and edges untouched", {
    gc <- graph_constraint(
        hyp_constraint = c(NA, 0.2, NA),
        trans_constraint = rbind(
            c(0, NA, NA),
            c(1, 0, 0),
            c(NA, NA, 0)
        )
    )

    w <- c(0.7999, 0.2, 0.0001)
    trans <- rbind(
        c(0, 0.999, 0.001),
        c(1, 0, 0),
        c(0.001, 0.999, 0)
    )

    pruned <- prune_graph_gsd(
        pvals = pp_pvals_gsd,
        hyp_weight = w,
        trans_matrix = trans,
        trial_success = pp_gain,
        graph_constraint = gc,
        gamma = 1,
        verbose = "silent"
    )

    expect_identical(pruned$hyp_weight[[2L]], 0.2)
    expect_identical(pruned$trans_matrix[2L, ], c(1, 0, 0))
    expect_true(is_graph_valid(pruned$hyp_weight, pruned$trans_matrix))
})


test_that("prune_graph_gsd() respects a marginal power constraint", {
    gc <- graph_constraint_free(3)

    w <- c(0.5, 0.4999, 0.0001)
    trans <- rbind(
        c(0, 0.999, 0.001),
        c(0.999, 0, 0.001),
        c(0.5, 0.5, 0)
    )

    # a constraint the starting graph already satisfies on H3
    start_power <- calc_power_pvals_gsd(
        pp_pvals_gsd,
        hyp_weight = w,
        trans_matrix = trans
    )$local_power

    constraint <- c(NA, NA, start_power[[3L]])

    pruned <- prune_graph_gsd(
        pvals = pp_pvals_gsd,
        hyp_weight = w,
        trans_matrix = trans,
        trial_success = pp_gain,
        graph_constraint = gc,
        gamma = 1,
        power_constraint = constraint,
        verbose = "silent"
    )

    final_power <- calc_power_pvals_gsd(
        pp_pvals_gsd,
        hyp_weight = pruned$hyp_weight,
        trans_matrix = pruned$trans_matrix
    )$local_power

    expect_gte(final_power[[3L]], constraint[[3L]])
    expect_true(is_graph_valid(pruned$hyp_weight, pruned$trans_matrix))
})


test_that(".try_prune_gsd() reports acceptance on the gain", {
    w <- c(0.4, 0.3, 0.3)
    trans <- rbind(
        c(0, 0.5, 0.5),
        c(0.5, 0, 0.5),
        c(0.5, 0.5, 0)
    )
    value <- pp_value(w, trans)

    accepted <- .try_prune_gsd(
        pvals = pp_pvals_gsd,
        hyp_weight = w,
        trans_matrix = trans,
        trial_success = pp_gain,
        power_best = value,
        constrained_idx = integer(0),
        power_constraint = NULL
    )
    expect_true(accepted$accepted)
    expect_identical(accepted$power_best, value)

    rejected <- .try_prune_gsd(
        pvals = pp_pvals_gsd,
        hyp_weight = w,
        trans_matrix = trans,
        trial_success = pp_gain,
        power_best = value + 1,
        constrained_idx = integer(0),
        power_constraint = NULL
    )
    expect_false(rejected$accepted)
    expect_identical(rejected$power_best, value + 1)
})
