# Smallest allocated level represented in a boundary table. Repeated p-values
# that fall below the table are floored here rather than at zero, so that a
# tiny recycled allocation cannot reject them (design record, section 4.1).
gsd_grid_min <- 1e-14

#' Transform group sequential p-values into repeated p-values
#'
#' `transform_pvalues_gsd()` converts an array of nominal (raw) p-values from a
#' group sequential design into *repeated* p-values: for each hypothesis and
#' each analysis, the smallest level at which that p-value would cross its
#' alpha-spending boundary. Because the boundary is monotone in the allocated
#' level, the group sequential rejection rule "\eqn{p_{i,k}} below the boundary
#' for level \eqn{w_i \alpha}" becomes the fixed-sample rule
#' "\eqn{p^r_{i,k} \le w_i \alpha}", so a graphical procedure can be run on the
#' transformed values without ever evaluating a boundary.
#'
#' @details Boundaries are computed with [gsDesign::gsBound1()] on a log-spaced
#'   grid of allocated levels running from `1e-14` up to `alpha`, and the
#'   resulting monotone table is inverted by log-log linear interpolation.
#'
#'   Per hypothesis, writing \eqn{D} for the analyses with a non-missing
#'   information fraction and \eqn{k^{mat}} for the first of those that reaches
#'   full information (or the last analysis in \eqn{D} if none does):
#'
#'   * an analysis without data (`NA` information fraction) before
#'     \eqn{k^{mat}} is given a repeated p-value of 1, i.e. "cannot reject";
#'   * analyses after \eqn{k^{mat}} copy the repeated p-value of \eqn{k^{mat}}
#'     forward, because the statistic of a matured endpoint does not change (it
#'     may still be rejected later if recycling gives it more alpha). A
#'     non-missing raw p-value supplied after maturity is ignored, with a
#'     warning if it differs from the maturity value;
#'   * a hypothesis with a single analysis at full information is returned
#'     unchanged, since there the boundary *is* the level;
#'   * if `look_back` is `TRUE` for a hypothesis its repeated p-values are
#'     replaced by their running minimum over analyses (the *sequential*
#'     p-value), which allows a rejection on the strength of the evidence at
#'     any earlier analysis.
#'
#'   If a boundary table is not non-decreasing in the allocated level the
#'   spending function is not "well ordered", the group sequential graphical
#'   procedure is not valid for it, and the transform aborts.
#'
#' @param pvals A numeric array of raw p-values with dimensions
#'   `c(nsim, m, K)`: simulations by hypotheses by analyses.
#' @inheritParams rlang::args_dots_empty
#' @param info_frac Information fractions. Either a numeric vector of length
#'   `K`, applied to every hypothesis, or an `m` by `K` numeric matrix with
#'   `NA` at analyses where a hypothesis has no data. Defaults to the
#'   `info_frac` attribute of `pvals`.
#' @param spending An alpha-spending function, or a list of `m` of them (one
#'   per hypothesis). A spending function takes a level and a vector of
#'   information fractions and returns the cumulative level spent, either as a
#'   plain numeric vector or as an object with a `spend` element (as the
#'   `gsDesign` spending functions do). It must spend the whole level by full
#'   information.
#' @param alpha A number giving the overall one-sided significance level.
#'   Default is `0.025`.
#' @param look_back A logical of length 1 or `m`. `TRUE` switches a hypothesis
#'   to sequential p-values (see Details). Default is `FALSE`.
#' @param grid_size A whole number giving the number of allocated levels in
#'   each boundary table. Default is `1024`.
#'
#' @returns A `multigrain_pvals_gsd` object made up of:
#'   * `pvals`: the transformed array, of the same dimensions as the input.
#'   * `nsim`, `m`, `K`: simulations, hypotheses and analyses.
#'   * `alpha`: the level the boundary tables were built up to.
#'   * `info_frac`: the `m` by `K` matrix of information fractions.
#'   * `look_back`: the per-hypothesis semantics actually used.
#'   * `spending`: a label per spending function.
#'   * `tables`: per hypothesis, the analyses that carry distinct information,
#'     the analysis at which it matured, and the grid of allocated levels with
#'     the nominal boundaries computed at each of them.
#'
#' @export
#' @examples
#' # two hypotheses analysed at 60% and 100% of the planned information
#' set.seed(1)
#' raw <- array(stats::runif(100 * 2 * 2), dim = c(100L, 2L, 2L))
#'
#' pvals_gsd <- transform_pvalues_gsd(
#'     raw,
#'     info_frac = c(0.6, 1),
#'     spending = gsDesign::sfLDOF,
#'     grid_size = 128L
#' )
#'
#' pvals_gsd
transform_pvalues_gsd <- function(
    pvals,
    ...,
    info_frac = NULL,
    spending, # nolint: function_argument_linter. required, named after `...`
    alpha = 0.025,
    look_back = FALSE,
    grid_size = 1024L
) {
    rlang::check_dots_empty()
    rlang::check_required(spending)
    spending_label <- rlang::as_label(rlang::enquo(spending))

    check_double(pvals)
    rlang::check_number_decimal(alpha, min = 0, max = 1)
    rlang::check_number_whole(grid_size, min = 2)
    grid_size <- as.integer(grid_size)

    dims <- dim(pvals)
    if (length(dims) != 3L) {
        cli::cli_abort(
            "{.arg pvals} must be a three-dimensional array, \\
            {.code c(nsim, m, K)}, not one with {length(dims)} dimension{?s}."
        )
    }
    if (any(pvals < 0 | pvals > 1, na.rm = TRUE)) {
        cli::cli_abort("{.arg pvals} must contain p-values between 0 and 1.")
    }

    n_sim <- dims[[1L]]
    m <- dims[[2L]]
    n_look <- dims[[3L]]

    if (is.null(info_frac)) {
        info_frac <- attr(pvals, "info_frac")
    }
    info_frac <- .gsd_info_frac(info_frac, m = m, n_look = n_look)
    look_back <- .gsd_look_back(look_back, m = m)
    spending <- .gsd_spending_list(spending, m = m)

    for (i in seq_len(m)) {
        .gsd_check_spending(spending[[i]], alpha = alpha, hyp = i)
    }

    out <- array(1, dim = dims, dimnames = dimnames(pvals))
    tables <- vector("list", m)

    for (i in seq_len(m)) {
        transformed <- .gsd_transform_hyp(
            matrix(pvals[, i, ], nrow = n_sim, ncol = n_look),
            t_row = info_frac[i, ],
            spending = spending[[i]],
            alpha = alpha,
            grid_size = grid_size,
            look_back = look_back[[i]],
            hyp = i
        )
        out[, i, ] <- transformed$values
        tables[[i]] <- transformed$table
    }

    new_pvals_gsd(
        pvals = out,
        nsim = n_sim,
        m = m,
        n_look = n_look,
        alpha = alpha,
        info_frac = info_frac,
        look_back = look_back,
        spending = .gsd_spending_labels(spending, label = spending_label),
        tables = tables
    )
}

# Steps 1 and 4 to 7 of section 4.1, for one hypothesis. `raw` is the
# simulations by analyses matrix of that hypothesis's raw p-values.
.gsd_transform_hyp <- function(
    raw,
    t_row,
    spending,
    alpha,
    grid_size,
    look_back,
    hyp
) {
    n_look <- ncol(raw)
    look_info <- .gsd_looks(t_row, hyp = hyp)
    looks <- look_info$looks
    maturity <- look_info$maturity
    .gsd_warn_matured(raw, hyp = hyp, maturity = maturity, n_look = n_look)

    out <- matrix(1, nrow = nrow(raw), ncol = n_look)

    if (length(looks) == 1L && t_row[[looks]] >= 1) {
        # the boundary at level a is a itself: nothing to invert
        tab <- NULL
        out[, seq.int(looks, n_look)] <- raw[, looks]
    } else {
        tab <- .gsd_boundary_table(
            t_row[looks],
            spending = spending,
            alpha = alpha,
            grid_size = grid_size,
            hyp = hyp
        )
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


# Constructor ------------------------------------------------------------

new_pvals_gsd <- function(
    pvals = array(double(), dim = c(0L, 0L, 0L)),
    nsim = integer(),
    m = integer(),
    n_look = integer(),
    alpha = double(),
    info_frac = matrix(double()),
    look_back = logical(),
    spending = character(),
    tables = list()
) {
    structure(
        list(
            pvals = pvals,
            nsim = nsim,
            m = m,
            K = n_look,
            alpha = alpha,
            info_frac = info_frac,
            look_back = look_back,
            spending = spending,
            tables = tables
        ),
        class = "multigrain_pvals_gsd"
    )
}


#' @export
print.multigrain_pvals_gsd <- function(x, ...) {
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

    invisible(x)
}

#' @export
summary.multigrain_pvals_gsd <- function(object, ...) {
    if (is.null(object)) {
        return()
    }

    print(object)

    cli::cat_line(cli::style_underline("\nPer-hypothesis detail"), ":")
    detail <- data.frame(
        hypothesis = .gsd_hyp_labels(object$m),
        analyses = vapply(
            object$tables,
            function(tab) toString(tab$looks),
            character(1L)
        ),
        matured_at = vapply(
            object$tables,
            function(tab) tab$maturity,
            integer(1L)
        ),
        look_back = object$look_back,
        spending = object$spending
    )
    print(detail, row.names = FALSE)

    cli::cat_line(cli::format_inline(
        "\nNominal boundaries at the full level ({format(object$alpha)}):"
    ))
    print(.gsd_hyp_matrix(
        .gsd_full_bounds(object),
        m = object$m,
        n_look = object$K
    ))

    invisible(object)
}


# Boundary tables --------------------------------------------------------

# Steps 2 and 3 of section 4.1: a monotone table of nominal boundary against
# allocated level, one row per grid level and one column per analysis.
.gsd_boundary_table <- function(
    t_look,
    spending,
    alpha,
    grid_size,
    hyp = NULL,
    call = rlang::caller_env()
) {
    grid_levels <- exp(seq(
        log(gsd_grid_min),
        log(alpha),
        length.out = grid_size
    ))
    lower <- rep(-20, length(t_look))
    bounds <- matrix(NA_real_, nrow = grid_size, ncol = length(t_look))

    for (g in seq_len(grid_size)) {
        cum_spend <- .gsd_spend(
            spending,
            alpha = grid_levels[[g]],
            t_look = t_look,
            hyp = hyp,
            call = call
        )
        b <- gsDesign::gsBound1(
            theta = 0,
            I = t_look,
            a = lower,
            probhi = diff(c(0, cum_spend))
        )$b
        bounds[g, ] <- stats::pnorm(b, lower.tail = FALSE)
    }

    .gsd_check_well_ordered(bounds, hyp = hyp, call = call)

    list(grid = grid_levels, bounds = bounds)
}

# Step 4 of section 4.1: invert one column of a boundary table.
.gsd_invert <- function(
    bounds,
    grid,
    p,
    floor = gsd_grid_min,
    call = rlang::caller_env()
) {
    keep <- !duplicated(bounds) & bounds > 0
    if (sum(keep) < 2L) {
        cli::cli_abort(
            c(
                "The boundary table cannot be inverted.",
                x = "Fewer than two distinct positive boundaries were found.",
                i = "Try a larger {.arg grid_size} or a later information \\
                fraction."
            ),
            call = call
        )
    }

    out <- exp(stats::approx(
        log(bounds[keep]),
        log(grid[keep]),
        xout = log(p),
        rule = 1
    )$y)
    # above the table nothing rejects at any allocation up to alpha; below it
    # the level is floored rather than set to zero
    out[which(p >= max(bounds))] <- 1
    out[which(p < min(bounds[keep]))] <- floor
    out[which(is.na(p))] <- 1

    out
}

# The phi(a, 1) = a check of section 4.1.
.gsd_check_spending <- function(
    spending,
    alpha,
    hyp = NULL,
    call = rlang::caller_env()
) {
    spent <- .gsd_spend(
        spending,
        alpha = alpha,
        t_look = 1,
        hyp = hyp,
        call = call
    )

    if (abs(spent - alpha) > 1e-6 * alpha) {
        cli::cli_abort(
            c(
                "Spending function for hypothesis {hyp} does not spend \\
                the whole level by full information.",
                x = "It spends {spent} of {alpha} at an information \\
                fraction of 1.",
                i = "A spending function that stops short would silently \\
                waste alpha."
            ),
            call = call
        )
    }

    invisible(TRUE)
}

.gsd_check_well_ordered <- function(
    bounds,
    hyp = NULL,
    call = rlang::caller_env()
) {
    for (l in seq_len(ncol(bounds))) {
        column <- bounds[, l]
        tol <- sqrt(.Machine$double.eps) * max(column)
        if (any(diff(column) < -tol)) {
            cli::cli_abort(
                c(
                    "Spending function for hypothesis {hyp} is not well \\
                    ordered.",
                    x = "Its nominal boundary at analysis {l} decreases \\
                    as the allocated level increases.",
                    i = "The group sequential graphical procedure is not \\
                    valid for such a spending function."
                ),
                call = call
            )
        }
    }

    invisible(TRUE)
}


# Helpers ----------------------------------------------------------------

.gsd_spend <- function(
    spending,
    alpha,
    t_look,
    hyp = NULL,
    call = rlang::caller_env()
) {
    spent <- spending(alpha, t_look)
    if (is.list(spent) && !is.null(spent[["spend"]])) {
        spent <- spent[["spend"]]
    }
    spent <- as.numeric(spent)

    if (length(spent) != length(t_look) || anyNA(spent)) {
        cli::cli_abort(
            c(
                "The spending function for hypothesis {hyp} is not usable.",
                x = "It must return {length(t_look)} non-missing cumulative \\
                spend value{?s}, not {length(spent)}."
            ),
            call = call
        )
    }

    spent
}

# The looks that carry distinct information for one hypothesis, and the
# analysis at which it matures (section 4.1).
.gsd_looks <- function(t_row, hyp = NULL, call = rlang::caller_env()) {
    have <- which(!is.na(t_row))
    if (length(have) == 0L) {
        cli::cli_abort(
            "Hypothesis {hyp} has no analysis with an information fraction.",
            call = call
        )
    }

    t_have <- t_row[have]
    if (any(t_have <= 0)) {
        cli::cli_abort(
            "Information fractions for hypothesis {hyp} must be positive.",
            call = call
        )
    }
    if (is.unsorted(t_have)) {
        cli::cli_abort(
            "Information fractions for hypothesis {hyp} must be \\
            non-decreasing.",
            call = call
        )
    }

    full <- which(t_have >= 1)[1L]
    maturity <- if (is.na(full)) have[[length(have)]] else have[[full]]

    list(looks = have[have <= maturity], maturity = maturity)
}

.gsd_warn_matured <- function(raw, hyp, maturity, n_look) {
    if (maturity >= n_look) {
        return(invisible(FALSE))
    }

    reference <- raw[, maturity]
    for (k in seq.int(maturity + 1L, n_look)) {
        later <- raw[, k]
        differs <- !is.na(later) & (is.na(reference) | later != reference)
        if (any(differs)) {
            cli::cli_warn(c(
                "Raw p-values for hypothesis {hyp} after analysis \\
                {maturity} differ from the value at which it matured.",
                i = "The value at analysis {maturity} is used at every \\
                later analysis."
            ))
            return(invisible(TRUE))
        }
    }

    invisible(FALSE)
}

.gsd_info_frac <- function(
    info_frac,
    m,
    n_look,
    call = rlang::caller_env()
) {
    if (is.null(info_frac)) {
        cli::cli_abort(
            c(
                "{.arg info_frac} must be supplied.",
                i = "Either pass it directly or attach it to {.arg pvals} \\
                as an {.field info_frac} attribute."
            ),
            call = call
        )
    }
    if (!is.numeric(info_frac)) {
        cli::cli_abort("{.arg info_frac} must be numeric.", call = call)
    }

    if (is.matrix(info_frac)) {
        if (!identical(dim(info_frac), c(m, n_look))) {
            cli::cli_abort(
                "{.arg info_frac} must be a {m} by {n_look} matrix.",
                call = call
            )
        }
        out <- info_frac
    } else {
        if (length(info_frac) != n_look) {
            cli::cli_abort(
                "{.arg info_frac} must have {n_look} value{?s}, not \\
                {length(info_frac)}.",
                call = call
            )
        }
        out <- matrix(info_frac, nrow = m, ncol = n_look, byrow = TRUE)
    }

    storage.mode(out) <- "double"
    out
}

.gsd_look_back <- function(look_back, m, call = rlang::caller_env()) {
    check_logical(look_back, call = call)
    if (anyNA(look_back)) {
        cli::cli_abort("{.arg look_back} must not be missing.", call = call)
    }
    if (length(look_back) == 1L) {
        look_back <- rep(look_back, m)
    }
    if (length(look_back) != m) {
        cli::cli_abort(
            "{.arg look_back} must have 1 or {m} values, not \\
            {length(look_back)}.",
            call = call
        )
    }

    look_back
}

.gsd_spending_list <- function(spending, m, call = rlang::caller_env()) {
    if (is.function(spending)) {
        return(rep(list(spending), m))
    }
    all_functions <- is.list(spending) &&
        all(vapply(spending, is.function, logical(1L)))
    if (!all_functions) {
        cli::cli_abort(
            "{.arg spending} must be a function or a list of functions.",
            call = call
        )
    }
    if (length(spending) != m) {
        cli::cli_abort(
            "{.arg spending} must hold 1 or {m} functions, not \\
            {length(spending)}.",
            call = call
        )
    }

    spending
}

.gsd_spending_labels <- function(spending, label) {
    m <- length(spending)
    if (!is.null(names(spending)) && all(nzchar(names(spending)))) {
        return(names(spending))
    }
    if (m > 1L && length(unique(spending)) > 1L) {
        return(sprintf("%s[[%d]]", label, seq_len(m)))
    }

    rep(label, m)
}

.gsd_hyp_labels <- function(m) {
    paste0("H", seq_len(m))
}

.gsd_hyp_matrix <- function(x, m, n_look) {
    dimnames(x) <- list(
        .gsd_hyp_labels(m),
        paste0("analysis ", seq_len(n_look))
    )
    x
}

# The nominal boundary of each hypothesis at each analysis when it holds the
# full level: the last row of its boundary table, copied forward past maturity.
.gsd_full_bounds <- function(object) {
    out <- matrix(NA_real_, nrow = object$m, ncol = object$K)

    for (i in seq_len(object$m)) {
        tab <- object$tables[[i]]
        full <- if (is.null(tab$bounds)) {
            object$alpha
        } else {
            tab$bounds[nrow(tab$bounds), ]
        }
        out[i, tab$looks] <- full
        if (tab$maturity < object$K) {
            out[i, seq.int(tab$maturity + 1L, object$K)] <- out[i, tab$maturity]
        }
    }

    out
}
