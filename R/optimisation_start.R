# Helper function to create start graph matrix for GA::ga()
.build_start_matrix <- function(graph_constraint, start_graph) {
    stopifnot(is_graph_constraint(graph_constraint))
    m <- graph_constraint_get_m(graph_constraint)

    seeds <- list()

    # 1) Unconstrained baseline seed (respects constraint)
    seeds[[length(seeds) + 1L]] <- create_start_params(graph_constraint)

    # 2) Fixed-sequence seed
    w_fs <- c(1, rep(0, m - 1))
    g_fs <- matrix(0, m, m)
    diag(g_fs) <- 0
    if (m > 1) {
        for (i in seq_len(m - 1)) {
            g_fs[i, i + 1] <- 1
        }
    }
    proj_grph <- closest_graph_to_constraints(graph_constraint, w_fs, g_fs)
    seeds[[length(seeds) + 1L]] <- create_start_params(
        graph_constraint,
        w0 = proj_grph$hyp_weight,
        G0 = proj_grph$trans_matrix,
        sum_to_one_constraint = FALSE
    )

    # 3) User-provided start_graphs
    if (!is.null(start_graph) && length(start_graph)) {
        user_rows <- lapply(start_graph, function(g) {
            if (is.null(g$hyp_weight) && is.null(g$trans_matrix)) {
                return(NULL)
            }
            w0 <- if (!is.null(g$hyp_weight)) {
                as.numeric(g$hyp_weight)
            } else {
                rep(1, m)
            }
            if (length(w0) != m) {
                stop(
                    "start_graph$hyp_weight has wrong length.",
                    call. = FALSE
                )
            }
            G0 <- if (!is.null(g$trans_matrix)) {
                as.matrix(g$trans_matrix)
            } else {
                start_mat <- matrix(0, m, m)
                diag(start_mat) <- 0
                if (m > 1) {
                    start_mat[row(start_mat) != col(start_mat)] <- 1 / (m - 1)
                }
                start_mat
            }

            if (!all(dim(G0) == c(m, m))) {
                stop(
                    "start_graph$trans_matrix has wrong dim.",
                    call. = FALSE
                )
            }
            proj_grph <- closest_graph_to_constraints(graph_constraint, w0, G0)
            create_start_params(
                graph_constraint,
                w0 = proj_grph$hyp_weight,
                G0 = proj_grph$trans_matrix,
                sum_to_one_constraint = FALSE
            )
        })
        user_rows <- Filter(Negate(is.null), user_rows)
        seeds <- c(seeds, user_rows)
    }

    start_mat <- do.call(rbind, seeds)
    start_mat <- unique(start_mat, MARGIN = 1)
    storage.mode(start_mat) <- "double"
    start_mat
}

# Default placeholder detector
.is_default_start_graph <- function(x) {
    is.null(x) ||
        identical(x, list(list(hyp_weight = NULL, trans_matrix = NULL)))
}


# Validate user-supplied start_graphs for dimension compatibility
.validate_start_graphs <- function(
    start_graph,
    m,
    call = rlang::caller_env()
) {
    if (.is_default_start_graph(start_graph)) {
        return(invisible(NULL))
    }

    for (i in seq_along(start_graph)) {
        g <- start_graph[[i]]
        w <- g$hyp_weight
        G <- g$trans_matrix

        arg_name_w <- glue::glue("start_graph[[{i}]]$hyp_weight")

        check_double(
            w,
            allow_null = TRUE,
            arg = arg_name_w,
            call = call
        )

        check_length(
            w,
            m = m,
            allow_null = TRUE,
            arg = arg_name_w,
            call = call
        )

        arg_name_g <- glue::glue("start_graph[[{i}]]$trans_matrix")

        check_double_matrix(
            G,
            allow_null = TRUE,
            arg = arg_name_g,
            call = call
        )

        check_dim(
            G,
            m = m,
            allow_null = TRUE,
            arg = arg_name_g,
            call = call
        )
    }

    invisible(NULL)
}


#' Encode a full graph for the optimisers, guarding derived entries
#'
#' A derived entry (the last free entry of a row, or the last free weight) is
#' one minus the rest. When it should be zero, the rest is rescaled so the
#' derived entry decodes to a small positive value below the zeroing threshold
#' (5e-6 for edges, 5e-5 for weights) instead of an exact zero that floating
#' point could turn into a negative and so into a penalty.
#'
#' @param graph_constraint A `multigrain_graph_constraint` object.
#' @param hyp_weight (numeric) Hypothesis weights of the graph to encode.
#' @param trans_matrix (numeric matrix) Its transition matrix.
#'
#' @returns A plain numeric vector of encoded parameters.
#' @noRd
.encode_graph <- function(graph_constraint, hyp_weight, trans_matrix) {
    hc <- graph_constraint$hyp_constraint
    tc <- graph_constraint$trans_constraint

    x <- as.numeric(create_start_params(
        graph_constraint,
        w0 = hyp_weight,
        G0 = trans_matrix,
        sum_to_one_constraint = FALSE
    ))

    n_w <- max(sum(is.na(hc)) - 1L, 0L)
    if (n_w > 0L) {
        last_w <- max(which(is.na(hc)))
        target <- 1 - 5e-5 - sum(hc, na.rm = TRUE)
        w_sum <- sum(x[seq_len(n_w)])
        if (hyp_weight[last_w] < 1e-4 && w_sum > 0 && target > 0) {
            x[seq_len(n_w)] <- x[seq_len(n_w)] * (target / w_sum)
        }
    }

    param_rows <- .g_param_rows(hc, tc)
    row_target <- .zeroing_row_targets(tc)
    g_idx <- which(!is.na(param_rows))

    for (i in unique(param_rows[g_idx])) {
        idx <- g_idx[param_rows[g_idx] == i]
        last_g <- max(which(is.na(tc[i, ])))
        g_sum <- sum(x[idx])
        if (trans_matrix[i, last_g] < 1e-5 && g_sum > 0 && row_target[i] > 0) {
            x[idx] <- x[idx] * (row_target[i] / g_sum)
        }
    }

    x
}


#' Seed matrix for the simplification GA
#'
#' Rows, most valuable first: the encoded reference graph; its single-edge
#' removal neighbours; the usual seeds from `.build_start_matrix()`; and the
#' rows of the reference object's stored GA population when it is present and
#' of matching width. Duplicates are dropped and the result is truncated to
#' `pop_size`, because `GA::ga()` errors when `suggestions` has more rows than
#' `popSize`.
#'
#' @param graph_constraint A `multigrain_graph_constraint` object.
#' @param ref_graph (list) The reference graph, with `hyp_weight` and
#'   `trans_matrix`.
#' @param pop_size (integer) The GA's population size.
#' @param start_graph Optional list of user-supplied start graphs, as stored on
#'   the reference object.
#' @param population Optional matrix, the stage-1 GA's final population.
#'
#' @returns A numeric matrix with at most `pop_size` rows.
#' @noRd
.build_simplify_seeds <- function(
    graph_constraint,
    ref_graph,
    pop_size,
    start_graph = NULL,
    population = NULL
) {
    tc <- graph_constraint$trans_constraint
    fixed_edge <- !is.na(tc)
    G_ref <- ref_graph$trans_matrix
    w_ref <- ref_graph$hyp_weight
    tolerance <- sqrt(.Machine$double.eps)

    seeds <- list(.encode_graph(graph_constraint, w_ref, G_ref))

    cand <- which(G_ref != 0 & !fixed_edge, arr.ind = TRUE)
    for (k in seq_len(nrow(cand))) {
        i <- cand[k, 1L]
        j <- cand[k, 2L]
        G_try <- G_ref
        G_try[i, ] <- .redistribute_mass(
            G_ref[i, ],
            drop_idx = j,
            fixed_idx = which(fixed_edge[i, ])
        )
        # Skip candidates that do not reduce the edge count, and those that
        # leave the row unable to sum to one (no free recipient).
        if (sum(G_try != 0) >= sum(G_ref != 0)) {
            next
        }
        if (abs(sum(G_try[i, ]) - 1) > tolerance) {
            next
        }
        seeds[[length(seeds) + 1L]] <- .encode_graph(
            graph_constraint,
            w_ref,
            G_try
        )
    }

    mat <- do.call(rbind, seeds)
    mat <- rbind(mat, .build_start_matrix(graph_constraint, start_graph))

    if (!is.null(population) && ncol(population) == ncol(mat)) {
        mat <- rbind(mat, population)
    }

    mat <- unique(mat, MARGIN = 1)
    storage.mode(mat) <- "double"
    mat[seq_len(min(nrow(mat), pop_size)), , drop = FALSE]
}
