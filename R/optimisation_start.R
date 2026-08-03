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
