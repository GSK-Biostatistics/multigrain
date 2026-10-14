# Optimised Hwang-Shih-DeCani spending parameters: the tagged spending
# function, the spec object that `transform_pvalues_gsd()` returns when any
# hypothesis is tagged, and `resolve_pvalues_gsd()`, which turns a spec and a
# choice of parameter values back into a `multigrain_pvals_gsd` object
# (dev/hsd_opt_design_revised.md, sections 1.1 and 1.2).

# The lattice steps allowed by `sf_hsd_opt()`: binary fractions, so that
# `k * step` is exact and nodes coincide across parameters.
hsd_opt_steps <- c(0.5, 0.25, 0.125)

# The range of parameter values a user may ask for.
hsd_opt_limits <- c(-10, 2)


#' Hwang-Shih-DeCani spending with an optimised parameter
#'
#' `sf_hsd_opt()` tags a hypothesis whose Hwang-Shih-DeCani spending parameter
#' should be chosen by [graph_optimise_gsd()] jointly with the graph, rather
#' than fixed in advance. Pass it, alone or inside the list of spending
#' functions, to [transform_pvalues_gsd()], which then returns a
#' `multigrain_pvals_gsd_spec` object instead of transformed p-values.
#'
#' @details The parameter is searched over a lattice: the multiples of `step`
#'   inside `bounds`, together with both bounds. Every node's boundary table is
#'   built, and checked for well ordering, when the spec is created, so every
#'   value the search can try is a valid group sequential transform.
#'
#'   Tagged hypotheses with the same `id` share one parameter. Tagged
#'   hypotheses without an `id` each get their own, except that a single
#'   `sf_hsd_opt()` passed as the whole `spending` argument (and so recycled
#'   across hypotheses) is one parameter shared by all of them. Hypotheses
#'   sharing a parameter must use identical `bounds`, `step` and `start`.
#'
#'   A tagged hypothesis with a single analysis at full information is not
#'   affected by its spending function and takes no part in the search; if
#'   that leaves no parameter to optimise, [transform_pvalues_gsd()] returns an
#'   ordinary `multigrain_pvals_gsd` object with a message.
#'
#' @param id `NULL`, or a single string or number naming the parameter.
#'   Hypotheses with the same `id` share one parameter, and the `id` is the
#'   name under which its value is reported.
#' @param bounds A numeric vector of length 2: the lowest and highest value of
#'   the parameter to search. The lower bound must be at least -10, the upper
#'   at most 2, and the lower below the upper. Default is `c(-4, 1)`.
#' @param step The spacing of the lattice of values searched: one of `0.5`,
#'   `0.25` or `0.125`. Default is `0.25`.
#' @param start `NULL` or a number inside `bounds`: the value the search
#'   starts from, snapped to the nearest lattice node. `NULL` starts from the
#'   node nearest the midpoint of `bounds`.
#'
#' @returns A spending function `function(alpha, t, param)` of class
#'   `multigrain_sf_hsd_opt`, returning
#'   `gsDesign::sfHSD(alpha, t, param)$spend`, with attributes `id`, `bounds`,
#'   `step` and `start`.
#'
#' @seealso [gsDesign::sfHSD()], [transform_pvalues_gsd()],
#'   [resolve_pvalues_gsd()]
#' @export
#' @examples
#' # one parameter shared by hypotheses 2 and 3, searched over -4, -3.5, ..., -2
#' spending <- list(
#'     gsDesign::sfLDOF,
#'     sf_hsd_opt(id = "secondary", bounds = c(-4, -2), step = 0.5),
#'     sf_hsd_opt(id = "secondary", bounds = c(-4, -2), step = 0.5)
#' )
#'
#' set.seed(1)
#' raw <- array(stats::runif(100 * 3 * 2), dim = c(100L, 3L, 2L))
#' spec <- transform_pvalues_gsd(
#'     raw,
#'     info_frac = c(0.5, 1),
#'     spending = spending,
#'     grid_size = 64L
#' )
#' spec
sf_hsd_opt <- function(
    id = NULL,
    bounds = c(-4, 1),
    step = 0.25,
    start = NULL
) {
    id <- .hsd_check_id(id)
    bounds <- .hsd_check_bounds(bounds)
    step <- .hsd_check_step(step)
    start <- .hsd_check_start(start, bounds = bounds)

    fn <- function(alpha, t, param) {
        gsDesign::sfHSD(alpha, t, param)$spend
    }
    # the function's environment would otherwise hold the arguments, so that
    # two equal tags would not be `identical()`
    environment(fn) <- asNamespace("multigrain")
    attr(fn, "id") <- id
    attr(fn, "bounds") <- bounds
    attr(fn, "step") <- step
    attr(fn, "start") <- start
    class(fn) <- c("multigrain_sf_hsd_opt", "function")

    fn
}

# Validation of the `sf_hsd_opt()` arguments; each returns the normalised
# value.
.hsd_check_id <- function(id, call = rlang::caller_env()) {
    if (is.null(id)) {
        return(NULL)
    }
    valid <- (is.character(id) || is.numeric(id)) &&
        length(id) == 1L &&
        !is.na(id)
    if (!valid) {
        cli::cli_abort(
            "{.arg id} must be {.code NULL} or a single non-missing string \\
            or number.",
            call = call
        )
    }

    as.character(id)
}

.hsd_check_bounds <- function(bounds, call = rlang::caller_env()) {
    valid <- is.numeric(bounds) &&
        length(bounds) == 2L &&
        all(is.finite(bounds))
    if (!valid) {
        cli::cli_abort(
            "{.arg bounds} must be a finite numeric vector of length 2.",
            call = call
        )
    }
    bounds <- as.numeric(bounds)
    if (bounds[[1L]] < hsd_opt_limits[[1L]]) {
        cli::cli_abort(
            "The lower bound in {.arg bounds} must be at least \\
            {hsd_opt_limits[[1L]]}, not {bounds[[1L]]}.",
            call = call
        )
    }
    if (bounds[[2L]] > hsd_opt_limits[[2L]]) {
        cli::cli_abort(
            "The upper bound in {.arg bounds} must be at most \\
            {hsd_opt_limits[[2L]]}, not {bounds[[2L]]}.",
            call = call
        )
    }
    if (bounds[[1L]] >= bounds[[2L]]) {
        cli::cli_abort(
            "The lower bound in {.arg bounds} must be below the upper \\
            bound, not {bounds[[1L]]} and {bounds[[2L]]}.",
            call = call
        )
    }

    bounds
}

.hsd_check_step <- function(step, call = rlang::caller_env()) {
    valid <- is.numeric(step) &&
        length(step) == 1L &&
        !is.na(step) &&
        as.numeric(step) %in% hsd_opt_steps
    if (!valid) {
        cli::cli_abort(
            c(
                "{.arg step} must be one of {.or {hsd_opt_steps}}.",
                i = "Binary fractions keep every lattice node exact."
            ),
            call = call
        )
    }

    as.numeric(step)
}

.hsd_check_start <- function(start, bounds, call = rlang::caller_env()) {
    if (is.null(start)) {
        return(NULL)
    }
    valid <- is.numeric(start) && length(start) == 1L && is.finite(start)
    if (!valid) {
        cli::cli_abort(
            "{.arg start} must be {.code NULL} or a single finite number.",
            call = call
        )
    }
    start <- as.numeric(start)
    if (start < bounds[[1L]] || start > bounds[[2L]]) {
        cli::cli_abort(
            "{.arg start} must lie inside {.arg bounds} \\
            [{bounds[[1L]]}, {bounds[[2L]]}], not {start}.",
            call = call
        )
    }

    start
}

# TRUE for a function tagged by `sf_hsd_opt()`.
is_sf_hsd_opt <- function(x) {
    inherits(x, "multigrain_sf_hsd_opt")
}

# TRUE if `spending` is a tagged function or a list holding at least one.
.has_hsd_opt <- function(spending) {
    if (is_sf_hsd_opt(spending)) {
        return(TRUE)
    }

    is.list(spending) && any(vapply(spending, is_sf_hsd_opt, logical(1L)))
}

# The two-argument spending function a fixed helper receives for one node
# (build plan, reconciliation 1).
.hsd_node_fn <- function(f, node) {
    force(f)
    force(node)
    function(a, t) f(a, t, param = node)
}


# Groups and lattices ----------------------------------------------------

# The group map of a length-m spending list: `group[i]` is NA for an untagged
# entry and 1..d for a tagged one. Entries sharing a non-NULL `id` share a
# group, `id = NULL` entries each get their own, and a single tagged function
# recycled across hypotheses (`recycled = TRUE`) is one group. Group members
# must agree on `bounds`, `step` and `start`.
.hsd_groups <- function(
    spending,
    m,
    recycled = FALSE,
    call = rlang::caller_env()
) {
    tagged <- vapply(spending, is_sf_hsd_opt, logical(1L))
    keys <- rep(NA_character_, m)

    for (i in which(tagged)) {
        id <- attr(spending[[i]], "id")
        keys[[i]] <- if (!is.null(id)) {
            paste0("id:", id)
        } else if (recycled) {
            "recycled"
        } else {
            paste0("hyp:", i)
        }
    }

    group <- match(keys, unique(keys[tagged]))
    d <- if (any(tagged)) max(group, na.rm = TRUE) else 0L
    first <- match(seq_len(d), group)

    for (r in seq_len(d)) {
        members <- which(group == r)
        ref <- spending[[first[[r]]]]
        for (field in c("bounds", "step", "start")) {
            same <- vapply(
                spending[members],
                function(f) identical(attr(f, field), attr(ref, field)),
                logical(1L)
            )
            if (!all(same)) {
                cli::cli_abort(
                    c(
                        "Hypotheses {members} share a spending parameter but \\
                        differ in {.arg {field}} of {.fn sf_hsd_opt}.",
                        i = "Hypotheses with the same {.arg id} must use \\
                        identical {.arg bounds}, {.arg step} and \\
                        {.arg start}."
                    ),
                    call = call
                )
            }
        }
    }

    refs <- spending[first]
    ids <- lapply(refs, attr, which = "id")
    param_names <- .hsd_param_names(ids, call = call)
    bounds <- matrix(
        unlist(lapply(refs, attr, which = "bounds")),
        ncol = 2L,
        byrow = TRUE,
        dimnames = list(param_names, c("lower", "upper"))
    )
    start_values <- vapply(
        refs,
        function(f) attr(f, "start") %||% NA_real_,
        double(1L)
    )

    list(
        group = group,
        names = param_names,
        ids = ids,
        bounds = bounds,
        step = vapply(refs, attr, double(1L), which = "step"),
        start = start_values
    )
}

# Parameter names: the `id` where one was given, else "gamma" when there is
# one parameter and "gamma1", "gamma2", ... when there are several.
.hsd_param_names <- function(ids, call = rlang::caller_env()) {
    d <- length(ids)
    out <- if (d == 1L) "gamma" else paste0("gamma", seq_len(d))
    named <- !vapply(ids, is.null, logical(1L))
    out[named] <- unlist(ids[named])

    if (anyDuplicated(out)) {
        cli::cli_abort(
            c(
                "Spending parameter names must be unique.",
                x = "{.val {unique(out[duplicated(out)])}} names more than \\
                one parameter.",
                i = "Give each {.fn sf_hsd_opt} group a distinct {.arg id}."
            ),
            call = call
        )
    }

    out
}

# The lattice of one parameter: the multiples of `step` inside `bounds`, plus
# both bounds. `k * step` is exact because `step` is a binary fraction.
.hsd_lattice <- function(bounds, step) {
    lower <- bounds[[1L]]
    upper <- bounds[[2L]]
    k_lo <- ceiling(lower / step)
    k_hi <- floor(upper / step)
    inner <- if (k_lo <= k_hi) step * seq(k_lo, k_hi) else double()

    sort(union(c(lower, upper), inner))
}

# The node of `lattice` (sorted) nearest to `value`, ties to the lower node.
.hsd_snap <- function(value, lattice) {
    lattice[[which.min(abs(lattice - value))]]
}

# Key of a node table: the information schedule and the node value.
.hsd_table_key <- function(schedule_key, node) {
    paste0(schedule_key, "|", format(node, digits = 15))
}


# Node tables ------------------------------------------------------------

# Build-time guarantee that `.gsd_invert()` cannot abort inside the search:
# every column must hold at least two distinct positive boundaries.
.hsd_check_invertible <- function(
    tab,
    hyp,
    node,
    name = "gamma",
    looks = seq_len(ncol(tab$bounds)),
    call = rlang::caller_env()
) {
    for (l in seq_len(ncol(tab$bounds))) {
        column <- tab$bounds[, l]
        if (sum(!duplicated(column) & column > 0) < 2L) {
            cli::cli_abort(
                c(
                    "The boundary table of hypothesis {hyp} cannot be \\
                    inverted at {name} = {node}.",
                    x = "Fewer than two distinct positive boundaries were \\
                    found at analysis {looks[[l]]}.",
                    i = "Narrow {.arg bounds} of {.fn sf_hsd_opt}."
                ),
                call = call
            )
        }
    }

    invisible(TRUE)
}

# One boundary table per distinct (information schedule, node) pair among the
# active tagged hypotheses, each checked for well ordering (inside
# `.gsd_boundary_table()`) and invertibility. `alpha` and `grid_size` are
# constant within a spec, so the schedule alone keys a table.
.hsd_node_tables <- function(
    spending,
    info_frac,
    group,
    lattice,
    names,
    alpha,
    grid_size,
    call = rlang::caller_env()
) {
    m <- length(group)
    schedule_key <- rep(NA_character_, m)
    tables <- list()

    for (i in which(!is.na(group))) {
        t_row <- info_frac[i, ]
        looks <- .gsd_looks(t_row, hyp = i, call = call)$looks
        schedule_key[[i]] <- paste(
            format(t_row[looks], digits = 15),
            collapse = ","
        )
        f <- spending[[i]]
        name <- names[[group[[i]]]]

        for (node in lattice[[group[[i]]]]) {
            node_fn <- .hsd_node_fn(f, node)
            .gsd_check_spending(node_fn, alpha = alpha, hyp = i, call = call)

            key <- .hsd_table_key(schedule_key[[i]], node)
            if (!is.null(tables[[key]])) {
                next
            }
            tab <- rlang::try_fetch(
                .gsd_boundary_table(
                    t_row[looks],
                    spending = node_fn,
                    alpha = alpha,
                    grid_size = grid_size,
                    hyp = i,
                    call = call
                ),
                error = function(cnd) {
                    cli::cli_abort(
                        c(
                            "Hypothesis {i} is not well ordered at \\
                            {name} = {node}.",
                            i = "Narrow {.arg bounds} of {.fn sf_hsd_opt}."
                        ),
                        parent = cnd,
                        call = call
                    )
                }
            )
            .hsd_check_invertible(
                tab,
                hyp = i,
                node = node,
                name = name,
                looks = looks,
                call = call
            )
            tables[[key]] <- tab
        }
    }

    list(tables = tables, schedule_key = schedule_key)
}

# Steps 4 to 7 of `.gsd_transform_hyp()` on a prebuilt boundary table:
# inversion per look, matured copy-forward, NA looks = 1, running minimum
# under look-back. Duplicated from `.gsd_transform_hyp()`, not shared, so the
# fixed path is untouched (house rule).
.gsd_invert_hyp <- function(raw, t_row, tab, look_back) {
    n_look <- ncol(raw)
    look_info <- .gsd_looks(t_row)
    looks <- look_info$looks
    maturity <- look_info$maturity

    out <- matrix(1, nrow = nrow(raw), ncol = n_look)

    for (l in seq_along(looks)) {
        out[, looks[[l]]] <- .gsd_invert(
            tab$bounds[, l],
            grid = tab$grid,
            p = raw[, looks[[l]]]
        )
    }
    if (maturity < n_look) {
        out[, seq.int(maturity + 1L, n_look)] <- out[, maturity]
    }

    if (look_back && n_look > 1L) {
        for (k in seq.int(2L, n_look)) {
            out[, k] <- pmin(out[, k], out[, k - 1L])
        }
    }

    list(
        values = out,
        table = c(list(looks = looks, maturity = maturity), tab)
    )
}


# The spec builder -------------------------------------------------------

# Called by `transform_pvalues_gsd()` once `spending` is a length-m list that
# holds at least one `sf_hsd_opt()`. Validates and transforms the fixed
# hypotheses as the fixed path does, short-circuits tagged hypotheses that
# spending cannot affect, and builds every node table up front. Draws no
# random numbers.
.transform_pvalues_gsd_spec <- function(
    pvals,
    info_frac,
    look_back,
    spending,
    alpha,
    grid_size,
    spending_label,
    recycled = FALSE,
    call = rlang::caller_env()
) {
    dims <- dim(pvals)
    n_sim <- dims[[1L]]
    m <- dims[[2L]]
    n_look <- dims[[3L]]

    groups <- .hsd_groups(spending, m = m, recycled = recycled, call = call)
    group <- groups$group
    tagged <- !is.na(group)
    lattice <- Map(
        .hsd_lattice,
        lapply(seq_len(nrow(groups$bounds)), function(r) groups$bounds[r, ]),
        groups$step
    )
    start_nodes <- vapply(
        seq_along(lattice),
        function(r) {
            value <- groups$start[[r]]
            if (is.na(value)) {
                value <- mean(groups$bounds[r, ])
            }
            .hsd_snap(value, lattice[[r]])
        },
        double(1L)
    )

    hyp_labels <- .hsd_fixed_labels(
        spending,
        label = spending_label,
        expr = .hsd_spending_expr(call)
    )

    for (i in which(!tagged)) {
        .gsd_check_spending(spending[[i]], alpha = alpha, hyp = i, call = call)
    }

    fixed <- vector("list", m)
    fixed_tables <- vector("list", m)
    short_circuit <- rep(FALSE, m)

    for (i in seq_len(m)) {
        raw_i <- matrix(pvals[, i, ], nrow = n_sim, ncol = n_look)
        t_row <- info_frac[i, ]
        hyp_spending <- spending[[i]]

        if (tagged[[i]]) {
            look_info <- .gsd_looks(t_row, hyp = i, call = call)
            looks <- look_info$looks
            if (!(length(looks) == 1L && t_row[[looks]] >= 1)) {
                .gsd_warn_matured(
                    raw_i,
                    hyp = i,
                    maturity = look_info$maturity,
                    n_look = n_look
                )
                next
            }
            # one look at full information: the boundary is the level for
            # any spending function, so the start node stands in
            short_circuit[[i]] <- TRUE
            node <- start_nodes[[group[[i]]]]
            hyp_spending <- .hsd_node_fn(hyp_spending, node)
            .gsd_check_spending(
                hyp_spending,
                alpha = alpha,
                hyp = i,
                call = call
            )
            hyp_labels[[i]] <- sprintf(
                "sfHSD(param = %s)",
                format(node, digits = 15)
            )
        }

        transformed <- .gsd_transform_hyp(
            raw_i,
            t_row = t_row,
            spending = hyp_spending,
            alpha = alpha,
            grid_size = grid_size,
            look_back = look_back[[i]],
            hyp = i
        )
        fixed[[i]] <- transformed$values
        fixed_tables[[i]] <- transformed$table
    }

    # drop short-circuited hypotheses from the search; a group left empty
    # disappears and the rest are renumbered 1..d
    group[short_circuit] <- NA_integer_
    kept <- sort(unique(group[!is.na(group)]))
    group <- match(group, kept)
    d <- length(kept)

    if (d == 0L) {
        cli::cli_inform(c(
            i = "Every hypothesis tagged with {.fn sf_hsd_opt} has a single \\
            analysis at full information, where spending has no effect.",
            " " = "Returning a {.cls multigrain_pvals_gsd} object; there is \\
            no spending parameter to optimise."
        ))
        out <- array(1, dim = dims, dimnames = dimnames(pvals))
        for (i in seq_len(m)) {
            out[, i, ] <- fixed[[i]]
        }
        return(new_pvals_gsd(
            pvals = out,
            nsim = n_sim,
            m = m,
            n_look = n_look,
            alpha = alpha,
            info_frac = info_frac,
            look_back = look_back,
            spending = hyp_labels,
            tables = fixed_tables
        ))
    }

    param_names <- .hsd_param_names(groups$ids[kept], call = call)
    bounds <- groups$bounds[kept, , drop = FALSE]
    rownames(bounds) <- param_names
    lattice <- stats::setNames(lattice[kept], param_names)

    node_tables <- .hsd_node_tables(
        spending,
        info_frac = info_frac,
        group = group,
        lattice = lattice,
        names = param_names,
        alpha = alpha,
        grid_size = grid_size,
        call = call
    )

    new_pvals_gsd_spec(
        raw = pvals,
        nsim = n_sim,
        m = m,
        n_look = n_look,
        alpha = alpha,
        info_frac = info_frac,
        look_back = look_back,
        grid_size = grid_size,
        fixed = fixed,
        fixed_tables = fixed_tables,
        param = list(
            names = param_names,
            bounds = bounds,
            step = stats::setNames(groups$step[kept], param_names),
            start = stats::setNames(start_nodes[kept], param_names),
            lattice = lattice
        ),
        group = group,
        schedule_key = node_tables$schedule_key,
        node_tables = node_tables$tables,
        spending = .hsd_labels(
            hyp_labels,
            group = group,
            names = param_names,
            bounds = bounds
        ),
        short_circuit = short_circuit,
        spending_label = spending_label
    )
}

# Low-level constructor for a `multigrain_pvals_gsd_spec` object. The raw
# draws are stored under `raw`, never `pvals`, and the class does not inherit
# `multigrain_pvals_gsd`, so a spec cannot reach the kernel untransformed
# (build plan, reconciliation 3). Does no validation.
new_pvals_gsd_spec <- function(
    raw = array(double(), dim = c(0L, 0L, 0L)),
    nsim = integer(),
    m = integer(),
    n_look = integer(),
    alpha = double(),
    info_frac = matrix(double()),
    look_back = logical(),
    grid_size = integer(),
    fixed = list(),
    fixed_tables = list(),
    param = list(),
    group = integer(),
    schedule_key = character(),
    node_tables = list(),
    spending = character(),
    short_circuit = logical(),
    spending_label = character()
) {
    structure(
        list(
            raw = raw,
            nsim = nsim,
            m = m,
            K = n_look,
            alpha = alpha,
            info_frac = info_frac,
            look_back = look_back,
            grid_size = grid_size,
            fixed = fixed,
            fixed_tables = fixed_tables,
            param = param,
            group = group,
            schedule_key = schedule_key,
            node_tables = node_tables,
            spending = spending,
            short_circuit = short_circuit,
            spending_label = spending_label
        ),
        class = "multigrain_pvals_gsd_spec"
    )
}


# Labels -----------------------------------------------------------------

# The `spending` expression of the `transform_pvalues_gsd()` call whose frame
# is `env`, or NULL if it cannot be recovered.
.hsd_spending_expr <- function(env) {
    frame_call <- tryCatch(rlang::frame_call(env), error = function(cnd) NULL)
    if (!is.call(frame_call)) {
        return(NULL)
    }

    frame_call[["spending"]]
}

# Labels for the hypotheses with fixed spending (build plan, reconciliation
# 7): the list names if all are set, else the i-th argument of a literal
# `list(...)` call as written, else the fixed path's labels.
.hsd_fixed_labels <- function(spending, label, expr = NULL) {
    m <- length(spending)
    if (!is.null(names(spending)) && all(nzchar(names(spending)))) {
        return(names(spending))
    }

    is_list_call <- is.call(expr) &&
        identical(expr[[1L]], quote(list)) &&
        length(expr) == m + 1L
    if (is_list_call) {
        return(vapply(
            as.list(expr)[-1L],
            deparse1,
            character(1L),
            collapse = " ",
            width.cutoff = 500L
        ))
    }

    .gsd_spending_labels(spending, label = label)
}

# Per-hypothesis spending labels. Fixed and short-circuited hypotheses keep
# `fixed_labels`; an active hypothesis reads `sfHSD(param = <value>)` when
# `values` (one per parameter) is given, and
# `sfHSD(param = <name> in [lower, upper])` on the spec.
.hsd_labels <- function(fixed_labels, group, names, bounds, values = NULL) {
    out <- fixed_labels

    for (i in which(!is.na(group))) {
        r <- group[[i]]
        out[[i]] <- if (is.null(values)) {
            sprintf(
                "sfHSD(param = %s in [%s, %s])",
                names[[r]],
                format(bounds[r, 1L], digits = 15),
                format(bounds[r, 2L], digits = 15)
            )
        } else {
            sprintf("sfHSD(param = %s)", format(values[[r]], digits = 15))
        }
    }

    out
}


# Resolve ----------------------------------------------------------------

#' Resolve a spending spec at chosen parameter values
#'
#' When any spending function passed to [transform_pvalues_gsd()] is tagged
#' with [sf_hsd_opt()], the transform returns a `multigrain_pvals_gsd_spec`:
#' the raw p-values, the hypotheses with fixed spending already transformed,
#' and the boundary table of every lattice node of every tagged parameter.
#' `resolve_pvalues_gsd()` picks one node per parameter and returns the
#' `multigrain_pvals_gsd` object that [transform_pvalues_gsd()] would have
#' returned had the tagged hypotheses used
#' `function(a, t) gsDesign::sfHSD(a, t, value)$spend`, using the stored
#' tables, so it does not rebuild any boundary.
#'
#' @param spec A `multigrain_pvals_gsd_spec` object from
#'   [transform_pvalues_gsd()].
#' @param spending_param A numeric vector with one value per spending
#'   parameter of `spec`. If named, the names must be the parameter names (in
#'   any order). Every value must be a node of its parameter's lattice.
#' @param x A `multigrain_pvals_gsd_spec` object.
#' @param ... Not used.
#'
#' @returns A `multigrain_pvals_gsd` object; the `spending` label of each
#'   tagged hypothesis reads, for example, `sfHSD(param = -2.5)`.
#'
#' @seealso [sf_hsd_opt()], [transform_pvalues_gsd()]
#' @export
#' @examples
#' spending <- list(
#'     gsDesign::sfLDOF,
#'     sf_hsd_opt(id = "secondary", bounds = c(-4, -2), step = 0.5)
#' )
#'
#' set.seed(1)
#' raw <- array(stats::runif(100 * 2 * 2), dim = c(100L, 2L, 2L))
#' spec <- transform_pvalues_gsd(
#'     raw,
#'     info_frac = c(0.5, 1),
#'     spending = spending,
#'     grid_size = 64L
#' )
#'
#' resolve_pvalues_gsd(spec, c(secondary = -3.5))
resolve_pvalues_gsd <- function(spec, spending_param) {
    if (!inherits(spec, "multigrain_pvals_gsd_spec")) {
        rlang::stop_input_type(
            spec,
            "a multigrain_pvals_gsd_spec object (see sf_hsd_opt())"
        )
    }
    nodes <- .hsd_resolve_nodes(spec$param, spending_param)

    raw_draws <- spec$raw
    out <- array(1, dim = dim(raw_draws), dimnames = dimnames(raw_draws))
    tables <- spec$fixed_tables

    for (i in seq_len(spec$m)) {
        r <- spec$group[[i]]
        if (is.na(r)) {
            out[, i, ] <- spec$fixed[[i]]
            next
        }
        tab <- spec$node_tables[[.hsd_table_key(
            spec$schedule_key[[i]],
            nodes[[r]]
        )]]
        transformed <- .gsd_invert_hyp(
            matrix(raw_draws[, i, ], nrow = spec$nsim, ncol = spec$K),
            t_row = spec$info_frac[i, ],
            tab = tab,
            look_back = spec$look_back[[i]]
        )
        out[, i, ] <- transformed$values
        tables[[i]] <- transformed$table
    }

    new_pvals_gsd(
        pvals = out,
        nsim = spec$nsim,
        m = spec$m,
        n_look = spec$K,
        alpha = spec$alpha,
        info_frac = spec$info_frac,
        look_back = spec$look_back,
        spending = .hsd_labels(
            spec$spending,
            group = spec$group,
            names = spec$param$names,
            bounds = spec$param$bounds,
            values = nodes
        ),
        tables = tables
    )
}

# Validate `spending_param` against the parameters of a spec and return the
# lattice node of each, in the spec's parameter order.
.hsd_resolve_nodes <- function(
    param,
    spending_param,
    call = rlang::caller_env()
) {
    param_names <- param$names
    d <- length(param_names)

    valid <- is.numeric(spending_param) &&
        length(spending_param) == d &&
        all(is.finite(spending_param))
    if (!valid) {
        cli::cli_abort(
            "{.arg spending_param} must be {d} finite number{?s}, one per \\
            spending parameter ({.val {param_names}}).",
            call = call
        )
    }

    given <- names(spending_param)
    if (!is.null(given)) {
        if (!setequal(given, param_names) || anyDuplicated(given)) {
            cli::cli_abort(
                c(
                    "The names of {.arg spending_param} must be the \\
                    spending parameter names.",
                    x = "Expected {.val {param_names}}, got {.val {given}}."
                ),
                call = call
            )
        }
        spending_param <- spending_param[param_names]
    }

    nodes <- double(d)
    for (r in seq_len(d)) {
        value <- spending_param[[r]]
        lattice <- param$lattice[[r]]
        hit <- which(abs(value - lattice) <= 1e-12)
        if (length(hit) == 0L) {
            .hsd_abort_off_lattice(
                value,
                name = param_names[[r]],
                lattice = lattice,
                step = param$step[[r]],
                call = call
            )
        }
        nodes[[r]] <- lattice[[hit[[1L]]]]
    }

    nodes
}

.hsd_abort_off_lattice <- function(value, name, lattice, step, call) {
    below <- lattice[lattice < value]
    above <- lattice[lattice > value]
    nearest <- c(below[length(below)], above[1L])
    nearest_text <- if (length(nearest) == 1L) {
        paste("The nearest node is", nearest)
    } else {
        paste("The nearest nodes are", nearest[[1L]], "and", nearest[[2L]])
    }

    cli::cli_abort(
        c(
            "{.arg spending_param} value {value} for {name} is not a node \\
            of its lattice.",
            i = paste0(
                nearest_text,
                "; the lattice runs from {lattice[[1L]]} to \\
                {lattice[[length(lattice)]]} in steps of {step}."
            )
        ),
        call = call
    )
}


# Print ------------------------------------------------------------------

#' @rdname resolve_pvalues_gsd
#' @export
print.multigrain_pvals_gsd_spec <- function(x, ...) {
    if (is.null(x)) {
        return()
    }

    cli::cat_line(cli::format_inline("{.cls {class(x)}}"))
    cli::cat_line(sprintf(
        "%d simulation(s), %d hypothes(es), %d analys(es); alpha = %s",
        x$nsim,
        x$m,
        x$K,
        format(x$alpha)
    ))
    cli::cat_line("Information fractions (NA = no data at that analysis):")
    print(.gsd_hyp_matrix(x$info_frac, m = x$m, n_look = x$K))
    cli::cat_line(
        "Look-back: ",
        toString(paste0(.gsd_hyp_labels(x$m), " = ", x$look_back))
    )

    cli::cat_line("Spending:")
    hyp_lines <- paste0(.gsd_hyp_labels(x$m), ": ", x$spending)
    hyp_lines[x$short_circuit] <- paste(
        hyp_lines[x$short_circuit],
        "[short-circuited: one look at full information]"
    )
    cli::cat_line(hyp_lines)

    cli::cat_line("Spending parameters:")
    param <- x$param
    for (r in seq_along(param$names)) {
        cli::cat_line(sprintf(
            "%s: bounds [%s, %s], step %s, %d nodes, start %s",
            param$names[[r]],
            format(param$bounds[r, 1L]),
            format(param$bounds[r, 2L]),
            format(param$step[[r]]),
            length(param$lattice[[r]]),
            format(param$start[[r]])
        ))
    }
    cli::cat_line("Boundary tables built: ", length(x$node_tables))

    invisible(x)
}
