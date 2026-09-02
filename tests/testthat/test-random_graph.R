test_that("random_graph", {
    set.seed(1)
    first_random_graph <- graph_random(5)

    expect_named(
        first_random_graph$hyp_weight,
        c("H1", "H2", "H3", "H4", "H5")
    )

    expect_identical(
        dimnames(first_random_graph$trans_matrix),
        list(
            c("H1", "H2", "H3", "H4", "H5"),
            c("H1", "H2", "H3", "H4", "H5")
        )
    )

    set.seed(1)
    second_random_graph <- graph_random(5)

    expect_identical(first_random_graph, second_random_graph)

    set.seed(1)
    expect_snapshot({
        graph_random(6)
    })
})

test_that("random_graph with user supplied names", {
    graph_names <- c("A", "B", "C", "D", "E")

    random_graph <- graph_random(
        5,
        names = graph_names
    )

    expect_named(random_graph$hyp_weight, graph_names)

    expect_identical(
        dimnames(random_graph$trans_matrix),
        list(graph_names, graph_names)
    )
})

test_that("random_graph complains with names length mismatch", {
    expect_error(
        graph_random(
            5,
            names = c("A", "B", "C", "D")
        ),
        "`names` must have 5 elements. It has 4."
    )
})

test_that("random_graph with graph_constraint", {
    gc <- graph_constraint(
        hyp_constraint = c(0.5, NA, NA),
        trans_constraint = matrix(
            c(0, NA, NA, NA, 0, NA, NA, NA, 0),
            3, 3
        )
    )

    set.seed(42)
    rg <- graph_random(graph_constraint = gc)

    expect_length(rg$hyp_weight, 3)
    expect_equal(rg$hyp_weight[[1]], 0.5)
    expect_equal(sum(rg$hyp_weight), 1, tolerance = sqrt(.Machine$double.eps))
    expect_equal(
        dim(rg$trans_matrix),
        c(3, 3),
        tolerance = sqrt(.Machine$double.eps)
    )
    expect_true(all(diag(rg$trans_matrix) == 0))
    expect_equal(
        unname(rowSums(rg$trans_matrix)),
        c(1, 1, 1),
        tolerance = 1e-10
    )
})

test_that("random_graph infers m from graph_constraint", {
    gc <- graph_constraint(hyp_constraint = c(NA, NA, NA, NA))
    rg <- graph_random(graph_constraint = gc)
    expect_length(rg$hyp_weight, 4)
})

test_that("random_graph errors when m and graph_constraint disagree", {
    gc <- graph_constraint(hyp_constraint = c(NA, NA, NA))
    expect_error(
        graph_random(m = 5, graph_constraint = gc),
        "5.*3"
    )
})

test_that("random_graph errors when neither m nor graph_constraint", {
    expect_error(
        graph_random(),
        "must be supplied"
    )
})

test_that("random_graph uses names from graph_constraint", {
    gc <- graph_constraint(
        hyp_constraint = c(NA, NA, NA),
        names = c("A", "B", "C")
    )
    rg <- graph_random(graph_constraint = gc)
    expect_named(
        rg$hyp_weight,
        c("A", "B", "C")
    )
    expect_identical(
        dimnames(rg$trans_matrix),
        list(
            c("A", "B", "C"),
            c("A", "B", "C")
        )
    )
})

test_that("graph_optimal_random()", {
    a <- graph_optimal_random(5)

    expect_s3_class(a, "multigrain_graph_optimal")
})

test_that("random_weights with NULL hyp_constraint", {
    set.seed(1)
    expect_equal(
        random_weights(5),
        c(
            0.114424861968674,
            0.160372265625078,
            0.246879579197571,
            0.39140549984016,
            0.0869177933685177
        )
    )
})

test_that("random_weights complains", {
    # when sum of fixed constraints is greater than 1
    expect_error(
        random_weights(5, c(0.1, 0.2, 0.3, 0.4, 0.5)),
        "Sum of fixed constraints exceeds 1"
    )

    # when sum of constraints is less than 1 and there are no unconstrained
    # elements
    expect_error(
        random_weights(5, c(0.1, 0.2, 0.3, 0.2, 0.1)),
        "Constraints sum to less than 1, but no unconstrained elements"
    )
})

test_that("random_weights normalises", {
    expect_equal(
        random_weights(3, c(0.2, 0.3, 0.5 - 0.000000000000001)),
        c(0.2, 0.3, 0.5)
    )
})

test_that("random_transitions with NULL trans_constraint", {
    set.seed(1)
    expect_equal(
        random_transitions(3),
        matrix(
            c(
                0, 0.386785759694967, 0.18333527390685,
                0.41639759109067, 0, 0.81666472609315,
                0.58360240890933, 0.613214240305033, 0
            ),
            nrow = 3
        ),
        tolerance = 1e-8
    )
})

test_that("random_transitions complains", {
    # when sum of fixed constraints is greater than 1
    expect_error(
        random_transitions(
            3,
            matrix(
                c(
                    0, 0.8, 0.3,
                    0.4, 0, 0.6,
                    0.1, 0.9, 0
                ),
                nrow = 3,
                byrow = TRUE
            )
        ),
        "Sum of fixed constraints in row 1 exceeds 1"
    )

    # when sum of constraints is less than 1 and there are no unconstrained
    # elements
    expect_error(
        random_transitions(
            3,
            matrix(
                c(
                    0, 0.8, 0.2,
                    0.4, 0, 0.5,
                    0.1, 0.9, 0
                ),
                nrow = 3,
                byrow = TRUE
            )
        ),
        "Constraints in row 2 sum to less than 1, but no unconstrained elements"
    )
})
