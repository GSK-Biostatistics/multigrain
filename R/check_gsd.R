# Shared validation for the group sequential consumers of a
# `multigrain_pvals_gsd` object (design record, sections 4.7 and 6 P4). Every
# check here runs once, on the public path, before the kernel is called: the
# objective closure repeats none of them.

is_pvals_gsd <- function(x) {
    inherits(x, "multigrain_pvals_gsd")
}

check_pvals_gsd <- function(
    pvals,
    arg = rlang::caller_arg(pvals),
    call = rlang::caller_env(),
    allow_null = FALSE
) {
    if (!missing(pvals)) {
        if (is_pvals_gsd(pvals)) {
            return(invisible(NULL))
        }

        if (allow_null && is.null(pvals)) {
            return(invisible(NULL))
        }

        if (inherits(pvals, "multigrain_pvals_gsd_spec")) {
            cli::cli_abort(
                c(
                    "{.arg {arg}} is a spending spec, not a transformed \\
                    p-value object.",
                    i = "Pass it to {.fn graph_optimise_gsd} to optimise the \\
                    spending parameters, or call {.fn resolve_pvalues_gsd} \\
                    with chosen values."
                ),
                call = call
            )
        }
    }

    rlang::stop_input_type(
        pvals,
        "a multigrain_pvals_gsd object (see transform_pvalues_gsd())",
        allow_null = allow_null,
        arg = arg,
        call = call
    )
}

# The group sequential gain. A fixed-sample `multigrain_trial_success` is
# refused with its own message: its compiled function expects a logical
# rejection matrix, and the GSD consumers only ever have decision times.
check_trial_success_gsd <- function(
    trial_success,
    arg = rlang::caller_arg(trial_success),
    call = rlang::caller_env()
) {
    if (!missing(trial_success)) {
        if (is_trial_success_gsd(trial_success)) {
            return(invisible(NULL))
        }

        if (is_trial_success(trial_success)) {
            cli::cli_abort(
                c(
                    "{.arg {arg}} was created with {.fn trial_success} and \\
                    scores rejection indicators only.",
                    x = "A group sequential gain is evaluated on the analysis \\
                    at which each hypothesis was rejected.",
                    i = "Use {.fn trial_success_gsd}; a rejection-only gain \\
                    such as {.code r1 + r2} can be written there too."
                ),
                call = call
            )
        }
    }

    rlang::stop_input_type(
        trial_success,
        "a multigrain_trial_success_gsd object (see trial_success_gsd())",
        allow_null = FALSE,
        arg = arg,
        call = call
    )
}

# The gain must be defined over the same hypotheses, and over the same number
# of analyses whenever it fixes one. The second check is not cosmetic: a gain
# carrying discount tables indexes a C array with the decision time, so a time
# above its `K` reads out of bounds (design record, section 10 item 6).
.gsd_check_gain_dims <- function(
    trial_success,
    pvals,
    arg = rlang::caller_arg(trial_success),
    call = rlang::caller_env()
) {
    if (trial_success$m != pvals$m) {
        cli::cli_abort(
            c(
                "{.arg {arg}} and {.arg pvals} disagree on the number of \\
                hypotheses.",
                x = "The trial success function is defined over \\
                m = {trial_success$m}, but the p-values hold \\
                m = {pvals$m}."
            ),
            call = call
        )
    }

    if (!is.null(trial_success$K) && trial_success$K != pvals$K) {
        n_tab <- length(trial_success$tables)
        cli::cli_abort(
            c(
                "{.arg {arg}} and {.arg pvals} disagree on the number of \\
                analyses.",
                x = "The trial success function is defined over \\
                K = {trial_success$K}, but the p-values hold \\
                K = {pvals$K}.",
                i = if (n_tab > 0L) {
                    sprintf(
                        "Its discount %s (%s) %s length %d.",
                        if (n_tab == 1L) "table" else "tables",
                        toString(sprintf("`%s`", names(trial_success$tables))),
                        if (n_tab == 1L) "has" else "have",
                        trial_success$K
                    )
                }
            ),
            call = call
        )
    }

    invisible(NULL)
}

# Resolve the level to test at. `NULL` means the level the boundary tables
# were built up to; anything larger could never reject, because repeated
# p-values are capped at 1 above the top of the table (section 4.1, step 4).
.gsd_check_alpha <- function(alpha, pvals, call = rlang::caller_env()) {
    if (is.null(alpha)) {
        return(pvals$alpha)
    }

    rlang::check_number_decimal(alpha, min = 0, max = pvals$alpha, call = call)

    if (alpha <= 0 || alpha > pvals$alpha) {
        cli::cli_abort(
            c(
                "{.arg alpha} must be greater than 0 and at most \\
                {pvals$alpha}, not {alpha}.",
                i = "The repeated p-values were built up to \\
                {pvals$alpha} and are capped at 1 above it, so no allocation \\
                beyond that level could ever reject."
            ),
            call = call
        )
    }

    alpha
}

# The graph must be defined over the hypotheses of the p-values. Before this
# check a mismatch surfaced only as the kernel's C++ error "G must be an m x m
# matrix" (review finding), which names neither argument. `what` is a
# preformatted string such as "{.arg graph_constraint}" or "{.arg hyp_weight}".
.gsd_check_graph_dims <- function(
    m_graph,
    pvals,
    what,
    call = rlang::caller_env()
) {
    if (m_graph != pvals$m) {
        cli::cli_abort(
            c(
                paste(
                    what,
                    "and {.arg pvals} disagree on the number of \\
                hypotheses."
                ),
                x = paste(
                    what,
                    "describes {cli::qty(m_graph)}m = {m_graph} \\
                hypothes{?is/es}, but the p-values hold \\
                {cli::qty(pvals$m)}m = {pvals$m}."
                )
            ),
            call = call
        )
    }

    invisible(NULL)
}

# The `dim<-` reshape of section 4.1: column (k - 1) * m + i is hypothesis i
# at analysis k, which R's column-major layout gives without a rearrangement.
# This is also the single `anyNA()` assertion of section 10 item 11: the
# kernel compares with `<`, so an `NA` would silently mean "never reject".
# The range check goes with it (P4 decision, beyond the record): the array
# walk is already paid for, a negative repeated p-value is silently a
# rejection at any positive allocation, and one above 1 is silently "never".
.gsd_kernel_matrix <- function(pvals, call = rlang::caller_env()) {
    values <- pvals$pvals
    dim(values) <- c(pvals$nsim, pvals$m * pvals$K)

    if (anyNA(values)) {
        cli::cli_abort(
            c(
                "The transformed p-value array contains {.val {NA}}.",
                x = "The kernel would treat these as never rejectable.",
                i = "{.fn transform_pvalues_gsd} never emits {.val {NA}}; \\
                check how {.arg pvals} was built."
            ),
            call = call
        )
    }

    if (any(values < 0 | values > 1)) {
        cli::cli_abort(
            c(
                "The transformed p-value array contains values outside \\
                [0, 1].",
                x = "The kernel would read a negative value as a rejection \\
                at any positive allocation, and a value above 1 as never \\
                rejectable.",
                i = "{.fn transform_pvalues_gsd} emits repeated p-values in \\
                [1e-14, 1]; check how {.arg pvals} was built."
            ),
            call = call
        )
    }

    values
}

# The truth table of an efficacy stopping rule (one condition per interim
# analysis, a language object over `r1..rm`), evaluated once here so that the
# kernel only indexes it. Row index convention: row `mask + 1` holds the
# rejection set whose bitmask is `mask`, with hypothesis i (1-based) at bit
# i - 1; column k holds interim analysis k. The kernel
# (`src/graph_shortcut_gsd.cpp`) uses the same convention.
.gsd_stop_table <- function(
    conditions,
    m,
    K, # nolint: object_name_linter.
    call = rlang::caller_env()
) {
    n_cond <- length(conditions)

    if (m > 16L && n_cond > 0L) {
        cli::cli_abort(
            c(
                "A stopping rule is evaluated on every rejection set, 2^m \\
                rows per interim; m = {m} gives {2^m} rows.",
                i = "Stopping rules are limited to at most 16 hypotheses."
            ),
            call = call
        )
    }

    if (n_cond != K - 1L) {
        cli::cli_abort(
            c(
                "A stopping rule needs one condition per interim analysis.",
                x = "K = {K} analyses need {K - 1} condition{?s}, not \\
                {n_cond}."
            ),
            call = call
        )
    }

    n_row <- 2^m
    if (n_cond == 0L) {
        return(matrix(logical(), nrow = n_row, ncol = 0L))
    }

    syms <- paste0("r", seq_len(m))
    stop_table <- matrix(NA, nrow = n_row, ncol = n_cond)

    for (mask in seq_len(n_row) - 1L) {
        rejected <- bitwAnd(mask, bitwShiftL(1L, seq_len(m) - 1L)) > 0L
        env <- new.env(parent = baseenv())
        for (i in seq_len(m)) {
            assign(syms[[i]], rejected[[i]], envir = env)
        }

        for (k in seq_len(n_cond)) {
            result <- eval(conditions[[k]], env)
            if (!is.logical(result) || length(result) != 1L || is.na(result)) {
                .gsd_stop_condition_abort(
                    conditions[[k]],
                    result,
                    rejected,
                    k,
                    call = call
                )
            }
            stop_table[mask + 1L, k] <- result
        }
    }

    stop_table
}

.gsd_stop_condition_abort <- function(condition, result, rejected, k, call) {
    # nolint start: object_usage_linter. used inside the cli message
    cond_text <- deparse1(condition)
    set_text <- if (any(rejected)) {
        paste0("{", paste0("H", which(rejected), collapse = ", "), "}")
    } else {
        "none"
    }
    result_text <- if (is.logical(result) && length(result) == 1L) {
        "NA"
    } else {
        obj_type_friendly(result)
    }
    # nolint end

    cli::cli_abort(
        c(
            "Stopping condition {.code {cond_text}} must return a single \\
            {.code TRUE} or {.code FALSE}.",
            x = "At interim {k} with rejection set {set_text} it returned \\
            {result_text}."
        ),
        call = call
    )
}
