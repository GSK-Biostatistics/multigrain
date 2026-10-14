# Phase P4 of dev/hsd_opt_build_plan.md: print() and summary() of a result
# of the optimised-spending path. Every object is hand-built with
# graph_optimal(); no optimiser runs here.

hsd_print_gain_3 <- trial_success_gsd(
    d(t1) + d(t2) + d(t3),
    d = c(1, 0.75),
    verbose = "silent"
)

hsd_print_gain_4 <- trial_success_gsd(
    d(t1) + d(t2) + d(t3) + d(t4),
    d = c(1, 0.9),
    verbose = "silent"
)

# The profile rows of one parameter: re-optimised differences `diff` and
# fixed-graph differences `fixed_diff` from the optimum's gain `g_star`.
hsd_print_profile_rows <- function(name, value, diff, fixed_diff, g_star) {
    data.frame(
        param = rep(name, length(value)),
        value = value,
        gain = g_star - diff,
        gain_fixed_graph = g_star - fixed_diff,
        diff = diff,
        se = ifelse(diff == 0, 0, 1e-4 + abs(diff) / 7),
        n_eval = as.integer(40L + seq_along(value))
    )
}

hsd_print_graphs <- function(n, m) {
    rep(
        list(list(
            hyp_weight = rep(1 / m, m),
            trans_matrix = (1 - diag(m)) / (m - 1)
        )),
        n
    )
}

hsd_print_obj_d1 <- function(stopping_rule = NULL) {
    value <- seq(-4, 1, by = 0.25)
    # plateau of the re-optimised profile: [-3.25, -1.25]; of the fixed
    # graph: [-2.5, -2]
    prof <- hsd_print_profile_rows(
        "sec",
        value = value,
        diff = 4.5e-4 * (value + 2.25)^2,
        fixed_diff = 4e-3 * (value + 2.25)^2,
        g_star = 0.71234
    )
    attr(prof, "graphs") <- hsd_print_graphs(nrow(prof), 3L)

    coverage <- data.frame(
        param = rep("sec", length(value)),
        value = value,
        visits = c(0L, 0L, 0L, seq_len(18L)),
        best_fitness_subsample = c(NA, NA, NA, rep(0.7, 18L)),
        stringsAsFactors = FALSE
    )

    names_3 <- paste0("H", 1:3)
    graph_optimal(
        hyp_weight = c(H1 = 0.6, H2 = 0.4, H3 = 0),
        trans_matrix = matrix(
            c(0, 0.5, 0.5, 1, 0, 0, 0.5, 0.5, 0),
            nrow = 3L,
            byrow = TRUE,
            dimnames = list(names_3, names_3)
        ),
        constraints = graph_constraint_free(3),
        trial_success = hsd_print_gain_3,
        power = list(trial_success = 0.71234),
        stopping_rule = stopping_rule,
        spending = list(
            spending_param = c(sec = -2.25),
            hyp_spending_param = c(H1 = NA, H2 = -2.25, H3 = -2.25),
            spending = c(
                "sfLDOF",
                "sfHSD(param = -2.25)",
                "sfHSD(param = -2.25)"
            ),
            spending_bounds = matrix(
                c(-4, 1),
                nrow = 1L,
                dimnames = list("sec", NULL)
            ),
            spending_profile = prof,
            spending_search = coverage,
            profile_converged = TRUE
        )
    )
}

hsd_print_obj_d2 <- function() {
    value_1 <- seq(-3, -1, by = 0.5)
    value_2 <- seq(-3, -2, by = 0.5)
    rows_1 <- hsd_print_profile_rows(
        "gamma1",
        value = value_1,
        diff = c(2e-3, 3e-4, 0, 1e-4, 9e-4),
        fixed_diff = c(6e-3, 2e-3, 0, 4e-4, 3e-3),
        g_star = 0.6543
    )
    # a negative difference after the final recompute
    rows_2 <- hsd_print_profile_rows(
        "b",
        value = value_2,
        diff = c(7e-4, 0, -2e-4),
        fixed_diff = c(1e-3, 0, 1e-4),
        g_star = 0.6543
    )
    prof <- rbind(rows_1, rows_2)
    attr(prof, "graphs") <- hsd_print_graphs(nrow(prof), 4L)

    coverage <- data.frame(
        param = rep(c("gamma1", "b"), c(5L, 3L)),
        value = c(value_1, value_2),
        visits = NA_integer_,
        best_fitness_subsample = NA_real_,
        stringsAsFactors = FALSE
    )

    names_4 <- paste0("H", 1:4)
    graph_optimal(
        hyp_weight = c(H1 = 0.5, H2 = 0.2, H3 = 0.2, H4 = 0.1),
        trans_matrix = matrix(
            1 / 3,
            nrow = 4L,
            ncol = 4L,
            dimnames = list(names_4, names_4)
        ) -
            diag(1 / 3, 4L),
        constraints = graph_constraint_free(4),
        trial_success = hsd_print_gain_4,
        power = list(trial_success = 0.6543),
        spending = list(
            spending_param = c(gamma1 = -2, b = -2.5),
            hyp_spending_param = c(H1 = NA, H2 = -2, H3 = -2.5, H4 = NA),
            spending = c(
                "sfLDOF",
                "sfHSD(param = -2)",
                "sfHSD(param = -2.5)",
                "sfHSD(param = -2.5)"
            ),
            spending_bounds = matrix(
                c(-3, -3, -1, -2),
                nrow = 2L,
                dimnames = list(c("gamma1", "b"), NULL)
            ),
            spending_profile = prof,
            spending_search = coverage,
            profile_converged = FALSE
        )
    )
}


test_that("print() and summary() of a d = 1 result", {
    obj <- hsd_print_obj_d1()

    expect_snapshot(print(obj))
    expect_snapshot(summary(obj))
})

test_that("print() and summary() of a d = 2 result", {
    obj <- hsd_print_obj_d2()

    expect_snapshot(print(obj))
    expect_snapshot(summary(obj))
})

test_that("the plateau stops at a gap in the profile", {
    value <- seq(-3, 0, by = 0.5)
    # the run around -2 ends at -1.5: -1 is above the tolerance, although
    # -0.5 and 0 are back below it
    diff <- c(8e-4, 1e-4, 0, 2e-4, 9e-4, 1e-4, 0)

    expect_identical(.hsd_plateau(value, diff, at = -2), c(-2.5, -1.5))
    expect_identical(.hsd_plateau(value, diff, at = 0), c(-0.5, 0))
})

test_that("the spending block comes before the stopping-rule block", {
    obj <- hsd_print_obj_d1(stopping_rule = stopping_rule_gsd(r1))
    out <- paste(utils::capture.output(print(obj)), collapse = "\n")

    at_gain <- regexpr("Trial success function:", out, fixed = TRUE)
    at_spending <- regexpr("Spending functions:", out, fixed = TRUE)
    at_rule <- regexpr("Stopping rule:", out, fixed = TRUE)
    expect_gt(at_gain, 0)
    expect_lt(at_gain, at_spending)
    expect_lt(at_spending, at_rule)

    out <- paste(utils::capture.output(summary(obj)), collapse = "\n")
    at_profile <- regexpr("Spending-parameter profile (sec)", out, fixed = TRUE)
    at_rule <- regexpr("Stopping rule (group sequential)", out, fixed = TRUE)
    expect_gt(at_profile, 0)
    expect_lt(at_profile, at_rule)
})

test_that("a hypothesis map the values do not settle is shown as such", {
    # H1 is the first member of gamma1; either of H2 and H3 may open b
    groups <- .hsd_infer_groups(c(-2, -2, -2), c(-2, -2))
    expect_identical(groups, list(1L, 1:2, 1:2))
    # with two members the map is forced
    groups <- .hsd_infer_groups(c(-2, NA, -2), c(-2, -2))
    expect_identical(groups, list(1L, integer(), 2L))

    # distinct values settle every hypothesis
    groups <- .hsd_infer_groups(c(NA, -2, -2.5, -2), c(-2, -2.5))
    expect_identical(groups, list(integer(), 1L, 2L, 1L))

    obj <- hsd_print_obj_d1()
    obj$spending_param <- c(sec = -2.25, ter = -2.25)
    obj$spending_bounds <- rbind(obj$spending_bounds, ter = c(-4, 1))
    obj$hyp_spending_param[["H1"]] <- -2.25
    obj$spending[[1L]] <- "sfHSD(param = -2.25)"
    out <- utils::capture.output(.hsd_cat_print(obj))
    expect_match(out, "[optimised: sec or ter]", fixed = TRUE, all = FALSE)
    expect_match(out, "  ter: bounds -4 to 1", fixed = TRUE, all = FALSE)
})

test_that("a fixed-path result prints no spending block", {
    obj <- hsd_print_obj_d1()
    fields <- c(
        "spending_param",
        "hyp_spending_param",
        "spending",
        "spending_bounds",
        "spending_profile",
        "spending_search",
        "profile_converged"
    )
    obj[fields] <- NULL

    expect_no_match(
        utils::capture.output(print(obj)),
        "Spending",
        fixed = TRUE
    )
    expect_no_match(
        utils::capture.output(summary(obj)),
        "Spending",
        fixed = TRUE
    )
})

test_that("the P4 files use no superassignment", {
    files <- c(
        test_path("..", "..", "R", "optimisation_gsd_opt.R"),
        test_path("..", "..", "R", "graph_optimal.R")
    )
    files <- files[file.exists(files)]
    skip_if(length(files) == 0L, "package sources not available")

    for (f in files) {
        code <- readLines(f, warn = FALSE)
        code <- sub("#.*$", "", code)
        expect_false(any(grepl("<<-", code, fixed = TRUE)), label = f)
    }
})
