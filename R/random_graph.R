#' Generate random graph
#'
#' Generate a random graph, consisting of a vector of hypothesis weights and a
#' transition matrix. The weights and transition matrix are randomly generated
#' respecting any constraints provided, ensuring that the sum of the weights
#' equals 1 and each row of the transition matrix sums to 1.
#'
#' @param m An integer representing the number of hypotheses. It defines both
#'   the length of the weight vector and the dimensions (`m x m`) of the
#'   transition matrix. Optional when `graph_constraint` is supplied (inferred
#'   from constraint dimensions). If both `m` and `graph_constraint` are
#'   supplied, they must agree.
#' @param graph_constraint An optional graph constraint object created by
#'   [graph_constraint()]. When supplied, fixed elements are honoured and only
#'   free (`NA`) positions are randomised.
#' @inheritParams graph_constraint names
#'
#' @return A list containing:
#'
#' * `hyp_weight`: A numeric vector of length `m` representing the generated
#'   hypothesis weights.
#' * `trans_matrix`: A numeric matrix of dimension `m x m` representing the
#'   generated transition matrix.
#'
#' @export
#' @examples
#' # Generate a random graph for 5 hypotheses
#' random_graph <- graph_random(5)
#' print(random_graph$hyp_weight)  # Prints the weight vector
#' print(random_graph$trans_matrix)  # Prints the transition matrix
#'
#' # Generate a random graph respecting constraints
#' gc <- graph_constraint(
#'     hyp_constraint = c(0.5, NA, NA),
#'     trans_constraint = matrix(c(0, NA, NA, NA, 0, NA, NA, NA, 0), 3, 3)
#' )
#' random_graph <- graph_random(graph_constraint = gc)
graph_random <- function(
    m = NULL,
    graph_constraint = NULL,
    names = "auto"
) {
    rlang::check_number_whole(m, allow_null = TRUE)
    check_graph_constraint(graph_constraint, allow_null = TRUE)
    check_character(names)

    if (is.null(m) && is.null(graph_constraint)) {
        cli::cli_abort(
            "Either {.arg m} or {.arg graph_constraint} must be supplied. \\
            They cannot both be `NULL` at the same time."
        )
    }

    graph_constraint <- graph_constraint %||% graph_constraint_free(m)
    gc_m <- graph_constraint_get_m(graph_constraint)
    gc_names <- graph_constraint_get_names(graph_constraint)
    m <- m %||% gc_m

    if (m != gc_m) {
        cli::cli_abort(
            "{.arg m} is {m} but {.arg graph_constraint} has {gc_m} hypotheses."
        )
    }

    if (missing(names)) {
        names <- gc_names
    }

    if (length(names) != m) {
        cli::cli_abort(
            "{.arg names} must have {m} elements. It has {length(names)}."
        )
    }

    hyp_constraint <- graph_constraint$hyp_constraint
    trans_constraint <- graph_constraint$trans_constraint

    hyp_weight <- random_weights(m, hyp_constraint)
    names(hyp_weight) <- names

    trans_matrix <- random_transitions(m, trans_constraint)
    dimnames(trans_matrix) <- list(names, names)

    list(
        hyp_weight = hyp_weight,
        trans_matrix = trans_matrix
    )
}

#' Generate a random graph optimal
#'
#' A helper for testing and examples. It is not recommended to use it in routine
#' practice.
#'
#' @param ... Arguments passed down to [graph_random()].
#'
#' @returns a minimal `multigrain_graph_optimal`. To be used only for examples
#'   and tests.
#'
#' @export
#' @examples
#' graph_optimal_random(5)
graph_optimal_random <- function(...) {
    random_graph <- graph_random(...)

    graph_optimal(
        hyp_weight = random_graph$hyp_weight,
        trans_matrix = random_graph$trans_matrix
    )
}

# Internal helpers --------------------------------------------------------

#' Generate random weight vector with constraints
#'
#' @param m Integer. The number of hypothesis weights.
#' @param hyp_constraint A numeric vector of length `m` with constraints on the
#'   weights. Use `NA` for unconstrained elements. If `NULL`, all weights are
#'   free.
#' @return A numeric vector of length `m` summing to 1.
#' @noRd
random_weights <- function(m, hyp_constraint = NULL) {
    if (is.null(hyp_constraint)) {
        hyp_constraint <- hyp_constraint_free(m)
    }

    hyp_weights <- numeric(m)
    free_idx <- which(is.na(hyp_constraint))
    fixed_idx <- which(!is.na(hyp_constraint))

    hyp_weights[fixed_idx] <- hyp_constraint[fixed_idx]
    total_fixed <- sum(hyp_weights[fixed_idx])

    if (total_fixed > 1) {
        cli::cli_abort("Sum of fixed constraints exceeds 1.")
    }

    remaining <- 1 - total_fixed

    if (length(free_idx) == 0L) {
        if (abs(total_fixed - 1) > 1e-8) {
            cli::cli_abort(
                "Constraints sum to less than 1, but no unconstrained elements."
            )
        }
        return(normalise_sum(hyp_weights))
    }

    rand <- stats::runif(length(free_idx))
    rand <- rand * remaining / sum(rand)
    # Anchor last free element for exact sum
    rand[length(rand)] <- remaining - sum(rand[-length(rand)])
    hyp_weights[free_idx] <- rand

    normalise_sum(hyp_weights, fixed_idx = fixed_idx)
}


#' Generate random transition matrix with constraints
#'
#' @param m Integer. The dimension ($m \times m$) of the transition matrix.
#' @param trans_constraint A numeric matrix of dimension ($m \times m$) with
#'   constraints. Use `NA` for unconstrained elements. If `NULL`, diagonal
#'   elements default to 0 and all off-diagonal elements are free.
#' @return A numeric matrix of dimension ($m \times m$) with rows summing to 1
#'   and zero diagonal.
#' @noRd
random_transitions <- function(m, trans_constraint = NULL) {
    if (is.null(trans_constraint)) {
        trans_constraint <- trans_constraint_free(m)
    }

    rand_G <- matrix(0, m, m)

    for (i in seq_len(m)) {
        row_constraints <- trans_constraint[i, ]
        # Diagonal is always 0
        row_constraints[i] <- 0

        free_idx <- which(is.na(row_constraints))
        fixed_idx <- which(!is.na(row_constraints))

        rand_G[i, fixed_idx] <- row_constraints[fixed_idx]
        total_fixed <- sum(rand_G[i, fixed_idx])

        if (total_fixed > 1) {
            cli::cli_abort("Sum of fixed constraints in row {i} exceeds 1.")
        }

        remaining <- 1 - total_fixed

        if (length(free_idx) == 0L) {
            if (abs(total_fixed - 1) > 1e-8) {
                cli::cli_abort(
                    "Constraints in row {i} sum to less than 1, \\
                    but no unconstrained elements."
                )
            }
            next
        }

        rand <- stats::runif(length(free_idx))
        rand <- rand * remaining / sum(rand)
        rand_G[i, free_idx] <- rand

        # normalise_sum with fixed positions (including diagonal)
        rand_G[i, ] <- normalise_sum(
            rand_G[i, ],
            fixed_idx = fixed_idx
        )
    }

    rand_G
}
