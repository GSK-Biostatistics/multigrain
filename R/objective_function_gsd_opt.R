# The objective closure of the joint (graph, spending parameter) search and
# the slice bank it reads (dev/hsd_opt_design_revised.md, section 1.3 "Stage
# G" and section 5.2; build plan, phase P2). Nothing here draws random
# numbers or touches the fixed-spending path.

# Map clamped parameter genes to lattice nodes. Each of the `M` nodes owns an
# interval of width 1 / M of the gene, so every node has the same prior mass
# under a uniform gene; `u = 1` would be node `M + 1` and is mapped to node
# `M`. Genes outside [0, 1] are clamped first. Vectorised over `u` and `M`.
.hsd_node_index <- function(u, M) {
    # nolint: object_name_linter.
    u <- pmin(pmax(u, 0), 1)

    as.integer(pmin(M, 1L + floor(u * M)))
}

# Bytes held by a slice bank: one `n_idx` by `K` double matrix per active
# hypothesis and node of its lattice.
.hsd_bank_bytes <- function(spec, n_idx) {
    n_nodes <- lengths(spec$param$lattice)
    active <- spec$group[!is.na(spec$group)]

    8 * spec$K * n_idx * sum(n_nodes[active])
}

# Refuse a bank that would not fit comfortably in memory. Warns above 2 GB
# and aborts above `max_bank_gb` GB (1 GB = 1e9 bytes).
.check_bank_memory <- function(
    bytes,
    max_bank_gb = 4,
    call = rlang::caller_env()
) {
    remedies <- c(
        i = "Use a coarser {.arg step} or narrower {.arg bounds} in \\
        {.fn sf_hsd_opt}, or a smaller {.arg nsim_global} \\
        (see {.fn control_nsim_global})."
    )
    # nolint start: object_usage_linter. used inside the cli messages
    size <- format(bytes / 1e9, digits = 3)
    # nolint end

    if (bytes > max_bank_gb * 1e9) {
        cli::cli_abort(
            c(
                "The spending-parameter slice bank would need {size} GB, \\
                above the limit of {max_bank_gb} GB.",
                remedies
            ),
            call = call
        )
    }
    if (bytes > 2e9) {
        cli::cli_warn(
            c(
                "The spending-parameter slice bank needs {size} GB.",
                remedies
            ),
            call = call
        )
    }

    invisible(bytes)
}

# The slice bank of Stage G: for every active hypothesis `i` and every node
# `q` of its parameter's lattice, the `length(idx)` by `K` matrix of repeated
# p-values on rows `idx` of the raw draws (in `idx` order), from the stored
# node table. `slices[[a]][[q]]` belongs to hypothesis `hyp[[a]]`.
.hsd_bank <- function(spec, idx, call = rlang::caller_env()) {
    n_idx <- length(idx)
    bytes <- .hsd_bank_bytes(spec, n_idx)
    .check_bank_memory(bytes, call = call)

    hyp <- which(!is.na(spec$group))
    slices <- vector("list", length(hyp))

    for (a in seq_along(hyp)) {
        i <- hyp[[a]]
        raw_i <- spec$raw[idx, i, , drop = FALSE]
        raw_i <- matrix(raw_i, nrow = n_idx, ncol = spec$K)
        lattice <- spec$param$lattice[[spec$group[[i]]]]

        slices[[a]] <- lapply(lattice, function(node) {
            tab <- spec$node_tables[[.hsd_table_key(
                spec$schedule_key[[i]],
                node
            )]]
            .gsd_invert_hyp(
                raw_i,
                t_row = spec$info_frac[i, ],
                tab = tab,
                look_back = spec$look_back[[i]]
            )$values
        })
    }

    list(slices = slices, hyp = hyp, idx = idx, bytes = bytes)
}

# The kernel matrix of `resolve_pvalues_gsd(spec, start_param)` on rows
# `idx`. A helper, so that neither the full resolved array nor the subset
# outlives the call.
.hsd_working_matrix <- function(spec, start_param, idx) {
    resolved <- resolve_pvalues_gsd(spec, start_param)
    resolved$pvals <- resolved$pvals[idx, , , drop = FALSE]
    resolved$nsim <- length(idx)

    .gsd_kernel_matrix(resolved)
}

# Overwrite columns `cols[[j]]` of the working matrix `state$mat` with
# `values[[j]]`, in place. The matrix is taken out of the environment before
# it is modified, so the local name is its only reference and `[<-` does not
# duplicate it; it is put back before returning. Inside a function,
# `state$mat[, j] <- v` would copy the whole matrix on every call, even
# though the same statement at top level does not (d03 probed it only at top
# level).
.hsd_set_columns <- function(state, cols, values) {
    mat <- state$mat
    state$mat <- NULL
    for (j in seq_along(cols)) {
        mat[, cols[[j]]] <- values[[j]]
    }
    state$mat <- mat

    invisible(state)
}

# Fill the archive of the objective closure into the environment `arch`: per
# parameter, the best fitness seen at each node, the number of evaluations
# that reached the kernel there, and the graph genes of the best one.
.hsd_archive_init <- function(arch, param_names, n_nodes, p_graph) {
    arch$best <- stats::setNames(
        lapply(n_nodes, function(n) rep(-Inf, n)),
        param_names
    )
    arch$visits <- stats::setNames(lapply(n_nodes, integer), param_names)
    arch$genes <- stats::setNames(
        lapply(n_nodes, function(n) {
            matrix(NA_real_, nrow = n, ncol = p_graph)
        }),
        param_names
    )

    arch
}

# Record one evaluation that reached the kernel at nodes `node_q`.
.hsd_archive_update <- function(arch, node_q, fit, genes) {
    for (r in seq_along(node_q)) {
        node <- node_q[[r]]
        arch$visits[[r]][[node]] <- arch$visits[[r]][[node]] + 1L
        if (fit > arch$best[[r]][[node]]) {
            arch$best[[r]][[node]] <- fit
            arch$genes[[r]][node, ] <- genes
        }
    }

    invisible(arch)
}

#' Create the objective of the joint graph and spending-parameter search
#'
#' Twin of `create_obj_func_gsd()` for a `multigrain_pvals_gsd_spec`. The
#' encoded vector is the graph genes followed by one gene per spending
#' parameter, `x = c(graph genes [p_graph], u_1, ..., u_d)` (build plan,
#' reconciliation 2). Each parameter gene is clamped to the unit interval
#' and then mapped to a lattice node by `.hsd_node_index()`. The active
#' hypotheses' columns of the working kernel matrix are overwritten from
#' the bank slices of those nodes, and the graph genes are decoded and
#' penalised exactly as
#' `create_obj_func_gsd()` does (the block is copied, not shared).
#'
#' All mutable state lives in one environment, `state`, created in the
#' factory body: the working matrix `state$mat` (built once here, at
#' `start_param`, and never bound to a second name), the node of each
#' parameter it currently holds, `state$last_q`, and the archive. It is
#' updated with ordinary `<-` through the environment; there is no
#' superassignment. Columns of the matrix are overwritten by
#' `.hsd_set_columns()`, which takes the matrix out of `state` first so that
#' it is modified in place, never copied.
#'
#' @param spec A `multigrain_pvals_gsd_spec`.
#' @param bank The slice bank from `.hsd_bank(spec, idx)`; its rows are the
#'   rows the closure scores.
#' @param power_criterion,hyp_constraint,trans_constraint As for
#'   `create_obj_func_gsd()`.
#' @param alpha,num_threads,stop_rule As for `create_obj_func_gsd()`.
#' @param p_graph (integer) Number of graph genes, `ncol(.build_start_matrix(
#'   graph_constraint, start_graph))`.
#' @param start_param Lattice values, one per parameter, at which the working
#'   matrix is first filled. Defaults to the spec's start nodes.
#'
#' @returns A function `function(x)` returning the trial success measure (or
#'   the penalty of an infeasible graph), with attributes `archive` (the
#'   `state` environment: per parameter, lists `best` fitness, `visits` and
#'   best graph `genes` at each node, over evaluations that reached the
#'   kernel; it also holds `mat` and `last_q`), `d` and `p_graph`.
#' @noRd
create_obj_func_gsd_opt <- function(
    spec,
    bank,
    power_criterion,
    hyp_constraint,
    trans_constraint,
    alpha,
    p_graph,
    num_threads = 1L,
    stop_rule = NULL,
    start_param = spec$param$start
) {
    force(spec)
    force(bank)
    force(power_criterion)
    force(hyp_constraint)
    force(trans_constraint)
    force(alpha)
    force(p_graph)
    force(num_threads)
    force(stop_rule)
    force(start_param)

    K <- spec$K # nolint: object_name_linter. K is the design record's symbol
    m <- spec$m
    use_parallel <- num_threads >= 2L

    param <- spec$param
    d <- length(param$names)
    n_nodes <- lengths(param$lattice, use.names = FALSE)
    gene_pos <- p_graph + seq_len(d)

    hyp <- bank$hyp
    slices <- bank$slices
    hyp_group <- spec$group[hyp]
    cols <- lapply(hyp, function(i) (seq_len(K) - 1L) * m + i)
    members <- lapply(seq_len(d), function(r) which(hyp_group == r))

    start_nodes <- .hsd_resolve_nodes(param, start_param)

    state <- new.env(parent = emptyenv())
    .hsd_archive_init(state, param$names, n_nodes, p_graph)
    state$last_q <- vapply(
        seq_len(d),
        function(r) match(start_nodes[[r]], param$lattice[[r]]),
        integer(1L)
    )
    state$mat <- .hsd_working_matrix(spec, start_param, bank$idx)

    fn <- function(x) {
        xg <- x[seq_len(p_graph)]
        node_q <- .hsd_node_index(x[gene_pos], n_nodes)

        theta <- split_theta(xg, hyp_constraint)
        hyp_weight <- recover_full_weights(theta$w_pars, hyp_constraint)
        trans_matrix <- recover_full_trans_matrix(
            theta$g_pars,
            trans_constraint
        )

        if (anyNA(hyp_weight) || anyNA(trans_matrix)) {
            return(-1e6)
        }
        if (any(hyp_weight < 0)) {
            return(sum(hyp_weight[hyp_weight < 0]))
        }
        if (any(trans_matrix < 0)) {
            return(sum(trans_matrix[trans_matrix < 0]))
        }
        if (any(hyp_weight > 1)) {
            return(-sum(hyp_weight[hyp_weight > 1]))
        }
        if (any(trans_matrix > 1)) {
            return(-sum(trans_matrix[trans_matrix > 1]))
        }

        hyp_weight[hyp_weight < 1e-4] <- 0
        trans_matrix[trans_matrix < 1e-5] <- 0

        # only the hypotheses whose parameter moved since the last call
        moved <- unlist(members[node_q != state$last_q])
        if (length(moved)) {
            .hsd_set_columns(
                state,
                cols = cols[moved],
                values = lapply(moved, function(a) {
                    slices[[a]][[node_q[[hyp_group[[a]]]]]]
                })
            )
            state$last_q <- node_q
        }

        res <- if (use_parallel) {
            graph_shortcut_gsd_parallel(
                pvals = state$mat,
                alpha = alpha,
                w = hyp_weight,
                G = trans_matrix,
                K = K,
                num_threads = num_threads,
                grain_size = -1L,
                stop_rule = stop_rule
            )
        } else {
            graph_shortcut_gsd(
                pvals = state$mat,
                alpha = alpha,
                w = hyp_weight,
                G = trans_matrix,
                K = K,
                stop_rule = stop_rule
            )
        }

        fit <- power_criterion(res$time)

        .hsd_archive_update(state, node_q, fit, xg)

        fit
    }

    attr(fn, "archive") <- state
    attr(fn, "d") <- d
    attr(fn, "p_graph") <- p_graph

    fn
}
