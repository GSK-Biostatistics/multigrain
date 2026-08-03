test_that("choose_graph returns local when both valid, local >= GA", {
    ga_result <- list(
        ga_hyp_weight = c(0.5, 0.5),
        ga_trans_matrix = matrix(c(0, 1, 1, 0), 2),
        ga_trial_success = 0.7,
        is_graph_valid = TRUE
    )
    local_result <- list(
        local_hyp_weight = c(0.6, 0.4),
        local_trans_matrix = matrix(c(0, 1, 1, 0), 2),
        local_trial_success = 0.8,
        is_graph_valid = TRUE
    )
    sel <- choose_graph(ga_result, local_result)
    expect_identical(sel$source, "local")
    expect_identical(sel$hyp_weight, c(0.6, 0.4))
})

test_that("choose_graph returns GA when both valid and GA > local", {
    ga_result <- list(
        ga_hyp_weight = c(0.5, 0.5),
        ga_trans_matrix = matrix(c(0, 1, 1, 0), 2),
        ga_trial_success = 0.9,
        is_graph_valid = TRUE
    )
    local_result <- list(
        local_hyp_weight = c(0.6, 0.4),
        local_trans_matrix = matrix(c(0, 1, 1, 0), 2),
        local_trial_success = 0.8,
        is_graph_valid = TRUE
    )
    sel <- choose_graph(ga_result, local_result)
    expect_identical(sel$source, "global")
    expect_identical(sel$hyp_weight, c(0.5, 0.5))
})

test_that("choose_graph returns local when both valid, equal trial_success", {
    ga_result <- list(
        ga_hyp_weight = c(0.5, 0.5),
        ga_trans_matrix = matrix(c(0, 1, 1, 0), 2),
        ga_trial_success = 0.8,
        is_graph_valid = TRUE
    )
    local_result <- list(
        local_hyp_weight = c(0.6, 0.4),
        local_trans_matrix = matrix(c(0, 1, 1, 0), 2),
        local_trial_success = 0.8,
        is_graph_valid = TRUE
    )
    sel <- choose_graph(ga_result, local_result)
    expect_identical(sel$source, "local")
    expect_identical(sel$hyp_weight, c(0.6, 0.4))
})

test_that("choose_graph returns GA when local is invalid", {
    ga_result <- list(
        ga_hyp_weight = c(0.5, 0.5),
        ga_trans_matrix = matrix(c(0, 1, 1, 0), 2),
        ga_trial_success = 0.7,
        is_graph_valid = TRUE
    )
    local_result <- list(
        local_hyp_weight = c(0.6, 0.4),
        local_trans_matrix = matrix(c(0, 1, 1, 0), 2),
        local_trial_success = NULL,
        is_graph_valid = FALSE
    )
    sel <- choose_graph(ga_result, local_result)
    expect_identical(sel$source, "global")
})

test_that("choose_graph warns and returns local when both invalid", {
    ga_result <- list(
        ga_hyp_weight = c(0.5, 0.5),
        ga_trans_matrix = matrix(c(0, 1, 1, 0), 2),
        ga_trial_success = NULL,
        is_graph_valid = FALSE
    )
    local_result <- list(
        local_hyp_weight = c(0.6, 0.4),
        local_trans_matrix = matrix(c(0, 1, 1, 0), 2),
        local_trial_success = NULL,
        is_graph_valid = FALSE
    )
    expect_warning(
        sel <- choose_graph(ga_result, local_result),
        "Neither GA nor local"
    )
    expect_identical(sel$source, "local")
})

test_that("choose_graph returns local when ga_result is NULL", {
    local_result <- list(
        local_hyp_weight = c(0.6, 0.4),
        local_trans_matrix = matrix(c(0, 1, 1, 0), 2),
        local_trial_success = 0.8,
        is_graph_valid = TRUE
    )
    sel <- choose_graph(NULL, local_result)
    expect_identical(sel$source, "local")
    expect_identical(sel$hyp_weight, c(0.6, 0.4))
})
