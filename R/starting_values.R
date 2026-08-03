# Function to finding starting weights
# Split alpha equally if no user input
default_start_weights <- function(hyp_constraint) {
    # cast as replace_na does not coerce (NA is logical, so it will fail when
    # hyp_constraint contains only NAs)
    hyp_constraint <- as.numeric(hyp_constraint)

    number_nas <- sum(is.na(hyp_constraint))
    remaining_alpha <- 1 - sum(hyp_constraint, na.rm = TRUE)
    w_split <- remaining_alpha / number_nas
    output <- tidyr::replace_na(hyp_constraint, w_split)

    output
}


# Function to finding starting transition weights
# Split edge weights equally across all
# edges leaving a hypothesis
default_start_edges <- function(trans_matrix) {
    output <- trans_matrix
    nas_per_row <- rowSums(is.na(output))
    remaining_alpha <- 1 - rowSums(output, na.rm = TRUE)

    for (i in seq_len(nrow(output))) {
        output[i, ][is.na(output[i, ])] <- remaining_alpha[i] / nas_per_row[i]
    }

    output
}


#' Create reduced dimension starting values for optimisation procedures
#'
#' Initialize starting values for optimisation by setting the hypothesis weight
#' vector and transition matrix according to the constraints defined by the
#' `multigrain_graph_constraint` object. If optional initial values are provided
#' for the weights (`w0`) and transition matrix (`G0`), they are validated and
#' used. If not, default starting values are generated according to the
#' specified constraints.
#'
#' @param graph_constraint A `multigrain_graph_constraint` object defining
#'   constraints on hypothesis weights and the transition matrix.
#' @param w0 (optional) Hypothesis weight vector used as starting values for
#' optimisation. If not provided, the remaining \eqn{\alpha} will be uniformly
#' distributed across unconstrained hypotheses.
#' @param G0 (optional) Transition matrix used as starting values for
#' optimisation. If not provided, for each row \eqn{i}, the remaining transition
#' weights will be uniformly distributed across unconstrained elements of the
#' matrix row \eqn{g_{i1}, ..., g_{im}}.
#'
#' @details
#' The function creates an object of class `optParams` that contains the
#' reduced-dimension starting values for optimisation based on the constraints
#' provided by the `multigrain_graph_constraint` object. Specifically, it
#' generates starting values for both the hypothesis weight vector (`w0_opt`)
#' and the transition matrix (`G0_opt`). These starting values are either
#' provided explicitly by the user or are computed by uniformly distributing the
#' remaining \eqn{\alpha} across unconstrained elements.
#'
#' If the provided weight vector `w0` or transition matrix `G0` does not match
#' the fixed elements in the `multigrain_graph_constraint`, or if they do not
#' sum to 1 where required, an error will be thrown. Default starting values
#' will be computed if `w0` or `G0` are not supplied.
#'
#' @return
#' An object of class `optParams`, which contains the following components:
#'   * `w0_opt`: A numeric vector representing the unconstrained elements of
#'   the hypothesis weight vector.
#'   * `G0_opt`: A numeric vector representing the unconstrained elements of
#'   the transition matrix.
#'   * `graph_constraint`: The `multigrain_graph_constraint` object used to
#'   define the optimisation constraints. This object can be passed to
#'   optimisation procedures to define the starting values.
#'
#' @noRd
create_start_params <- function(
    graph_constraint,
    w0 = NULL,
    G0 = NULL,
    sum_to_one_constraint = TRUE
) {
    hyp_constraint <- graph_constraint$hyp_constraint
    trans_constraint <- graph_constraint$trans_constraint

    tolerance <- attr(graph_constraint, "tolerance")

    # nolint start: line_length_linter
    # Set starting hypothesis weights
    if (!is.null(w0)) {
        stopifnot(
            "Fixed weights in weight constraint vector must match elements in w0" = hyp_constraint[
                !is.na(hyp_constraint)
            ] ==
                w0[!is.na(hyp_constraint)]
        )

        # Sum-to-one check with tolerance (if enabled)
        if (sum_to_one_constraint && abs(sum(w0) - 1) > tolerance) {
            stop(
                "Weights must sum to 1 (within tolerance).",
                call. = FALSE
            )
        }

        w0_opt <- w0[is.na(hyp_constraint)][
            seq_along(w0[is.na(hyp_constraint)]) - 1
        ]
    } else {
        w0_opt <- default_start_weights(hyp_constraint)[is.na(hyp_constraint)]
        w0_opt <- w0_opt[seq_along(w0_opt) - 1]
    }

    # Set starting transition weights
    if (!is.null(G0)) {
        for (i in seq_len(nrow(G0))) {
            stopifnot(
                "Fixed elements in transition constraint matrix must match elements in G0" = G0[
                    i,
                ][!is.na(trans_constraint[i, ])] ==
                    trans_constraint[i, ][!is.na(trans_constraint[i, ])]
            )
        }
    }
    # nolint end

    # Row sums must be 1 within tolerance (if enabled)
    if (!is.null(G0) && sum_to_one_constraint) {
        rs <- rowSums(G0)
        if (any(abs(rs - 1) > tolerance)) {
            stop(
                "Rows of transition matrix must sum to 1 (within tolerance).",
                call. = FALSE
            )
        }
    }

    if (is.null(G0)) {
        G0 <- default_start_edges(trans_constraint)
    }

    free_g <- apply(trans_constraint, 1, function(x) sum(is.na(x)))
    G0_opt <- list()

    for (i in seq_len(nrow(G0))) {
        take <- max(free_g[i] - 1L, 0L)
        G0_opt[[i]] <- G0[i, ][is.na(trans_constraint[i, ])][seq_len(take)]
    }

    G0_opt <- unlist(G0_opt, use.names = FALSE)

    x0 <- structure(
        c(w0_opt, G0_opt),
        w0_opt = w0_opt,
        G0_opt = G0_opt,
        graph_constraint = graph_constraint,
        class = "optParams"
    )

    x0
}
