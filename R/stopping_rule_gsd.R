#' Create an efficacy stopping rule for a group sequential design
#'
#' Create a rule that stops a simulated group sequential trial for efficacy
#' at an interim analysis, given the hypotheses rejected so far. The rule is
#' one Boolean condition per interim analysis: when the condition for interim
#' analysis `k` is true, the trial stops after analysis `k` and no hypothesis
#' can be rejected at a later analysis. The final analysis always ends the
#' trial, so it has no condition.
#'
#' A condition is written using:
#'
#' * `r1, r2, ..., rm`: the indicator that hypothesis \eqn{H_i} has been
#'   rejected at the current interim analysis or at an earlier one;
#' * the logical operators `&&`, `||` and `!` (or the words `and`, `or` and
#'   `not` in a string);
#' * comparisons (`== != < <= > >=`) and arithmetic (`+ - * /`), so that a
#'   count of rejections such as `r1 + r2 + r3 >= 2` can be tested;
#' * numeric or logical scalars, and parentheses.
#'
#' Every condition must refer to at least one `r<i>`. Decision times (`t1`,
#' ...) and discount tables, which a group sequential gain can use (see
#' [trial_success_gsd()]), are not available to a stopping rule.
#'
#' @param condition The stopping condition, or one condition per interim
#'   analysis. To inject values from your R session, use rlang's injection
#'   operator [`!!`][rlang::injection-operator]. One of:
#'   * an expression, such as `r1 && (r2 || r3)`, or a single string, such as
#'     `"r1 and (r2 or r3)"`: one condition, applied at every interim
#'     analysis;
#'   * a call to `list()` with one expression or string per interim analysis,
#'     such as `list(r1, r1 && (r2 || r3))`, or a list of them injected with
#'     `!!`;
#'   * a character vector with one string per interim analysis, such as
#'     `c("r1", "r1 and (r2 or r3)")`, or one injected with `!!`.
#'
#' @returns A `multigrain_stopping_rule_gsd` object made up of:
#'   * `conditions`: a list of the conditions, as language objects.
#'   * `labels`: the conditions deparsed to strings, for display.
#'   * `m_min`: the largest hypothesis index referred to by an `r<i>` symbol,
#'     the fewest hypotheses the rule can be used with.
#'   * `n_interim`: `NULL` for a single condition applied at every interim
#'     analysis, otherwise the number of conditions supplied.
#'
#' @details
#'
#'   ## When the rule is applied
#'
#'   The rule does not store the number of hypotheses `m` or of analyses `K`;
#'   both are settled when it meets a `multigrain_pvals_gsd` object (see
#'   [transform_pvalues_gsd()]). The rule may then refer to no hypothesis
#'   beyond `m`. A single condition is recycled to each of the `K - 1` interim
#'   analyses; otherwise one condition per interim analysis is needed. A list
#'   of `K` conditions is accepted with a warning: the final analysis always
#'   ends the trial, so its condition is dropped. With `K = 1` there is no
#'   interim analysis, and a single condition is accepted and has no effect.
#'   The rule is evaluated once on every possible set of rejected hypotheses,
#'   `2^m` sets per interim analysis, so it can be used with at most 16
#'   hypotheses.
#'
#'   Stopping only removes rejections, so error control is unaffected: a
#'   stopped trial keeps the rejections it has already declared and declares
#'   none at a later analysis.
#'
#'   A trial that has rejected every hypothesis stops at that analysis whether
#'   or not a rule is supplied, as nothing is left to test.
#'
#'   ## Operator precedence
#'
#'   Operator precedence follows R's own rules. Parenthesise to make intent
#'   explicit: write `r1 && (r2 || r3)` rather than relying on `&&` binding
#'   tighter than `||`.
#'
#' @seealso [trial_success_gsd()] for the group sequential gain,
#'   [transform_pvalues_gsd()] for the p-values the rule is settled against.
#'
#' @export
#' @examples
#' # Stop at any interim analysis once H1 and at least one of H2 and H3 are
#' # rejected
#' stopping_rule_gsd(r1 && (r2 || r3))
#'
#' # One condition per interim analysis (K = 3): stop at the first interim
#' # once H1 is rejected, at the second only once H2 or H3 is rejected too
#' stopping_rule_gsd(list(r1, r1 && (r2 || r3)))
#'
#' # The first rule, written as a string
#' stopping_rule_gsd("r1 and (r2 or r3)")
#'
#' # Stop once at least two of three hypotheses are rejected
#' stopping_rule_gsd(r1 + r2 + r3 >= 2)
stopping_rule_gsd <- function(condition) {
    # Call flow: the captured argument is split into its conditions (parsing
    # strings), each condition is checked against the grammar and for at least
    # one `r<i>`, and the object records the largest hypothesis index. The
    # number of analyses is not known here: `.gsd_check_stopping_rule_dims()`
    # settles the conditions against the p-values.
    rlang::check_required(condition)
    captured <- .gsd_rule_capture(rlang::enexpr(condition))
    conditions <- captured$conditions

    cond_labels <- vapply(
        conditions,
        deparse1,
        character(1L),
        width.cutoff = 500L
    )

    m_min <- 0L
    for (j in seq_along(conditions)) {
        .gsd_rule_validate(conditions[[j]])
        m_min <- max(
            m_min,
            .gsd_rule_indices(conditions[[j]], cond_labels[[j]])
        )
    }

    structure(
        list(
            conditions = conditions,
            labels = cond_labels,
            m_min = m_min,
            n_interim = captured$n_interim
        ),
        class = "multigrain_stopping_rule_gsd"
    )
}


#' @export
print.multigrain_stopping_rule_gsd <- function(x, ...) {
    if (is.null(x)) {
        return()
    }

    cli::cat_line(cli::format_inline("{.cls {class(x)}}"))
    .gsd_rule_cat_detail(x)
    invisible(x)
}

#' @export
summary.multigrain_stopping_rule_gsd <- function(object, ...) {
    if (is.null(object)) {
        return()
    }

    cli::cat_line(
        cli::style_underline("\nStopping rule (group sequential)"),
        ":"
    )
    .gsd_rule_cat_detail(object)
    invisible(object)
}

# Shared body of print() and summary(): the conditions, then the hypotheses
# they reach.
.gsd_rule_cat_detail <- function(x) {
    if (is.null(x$n_interim)) {
        cli::cat_line("Every interim: ", x$labels[[1L]])
    } else {
        for (k in seq_along(x$labels)) {
            cli::cat_line("Interim ", k, ": ", x$labels[[k]])
        }
    }
    cli::cat_line("Hypotheses referenced: up to r", x$m_min)
    invisible(NULL)
}

is_stopping_rule_gsd <- function(x) {
    inherits(x, "multigrain_stopping_rule_gsd")
}

check_stopping_rule_gsd <- function(
    stopping_rule,
    arg = rlang::caller_arg(stopping_rule),
    call = rlang::caller_env(),
    allow_null = FALSE
) {
    if (!missing(stopping_rule)) {
        if (is_stopping_rule_gsd(stopping_rule)) {
            return(invisible(NULL))
        }

        if (allow_null && is.null(stopping_rule)) {
            return(invisible(NULL))
        }
    }

    rlang::stop_input_type(
        stopping_rule,
        "a multigrain_stopping_rule_gsd object (see stopping_rule_gsd())",
        allow_null = allow_null,
        arg = arg,
        call = call
    )
}

# Settle a stopping rule against the p-values: the rule may refer to no
# hypothesis beyond `m`, and the result holds one condition per interim
# analysis (`K - 1`, an empty list when `K = 1`), ready for
# `.gsd_stop_table(conditions, pvals$m, pvals$K)`.
.gsd_check_stopping_rule_dims <- function(
    stopping_rule,
    pvals,
    arg = rlang::caller_arg(stopping_rule),
    call = rlang::caller_env()
) {
    m_min <- stopping_rule$m_min
    if (m_min > pvals$m) {
        cli::cli_abort(
            "{.arg {arg}} refers to hypothesis {.code r{m_min}} but \\
            {.arg pvals} holds m = {pvals$m} hypothes{?is/es}.",
            call = call
        )
    }

    conditions <- stopping_rule$conditions
    n_interim <- stopping_rule$n_interim
    n_look <- pvals$K
    n_int <- n_look - 1L

    if (is.null(n_interim)) {
        return(rep(conditions, n_int))
    }

    if (n_interim == n_int) {
        return(conditions)
    }

    if (n_interim == n_look) {
        cli::cli_warn(
            c(
                "The final analysis always ends the trial; the last \\
                condition of {.arg {arg}}, \\
                {.code {stopping_rule$labels[[n_interim]]}}, is ignored.",
                i = "{.arg pvals} holds K = {n_look} analys{?is/es}. A rule \\
                gives one condition per interim analysis."
            ),
            call = call
        )
        return(conditions[-n_interim])
    }

    cli::cli_abort(
        c(
            "{.arg {arg}} must give one condition per interim analysis.",
            x = "It gives {n_interim} condition{?s}.",
            i = "{.arg pvals} holds K = {n_look} analys{?is/es}. The rule \\
            needs {n_int} condition{?s}, or a single condition (not in a \\
            list) applied at every interim analysis."
        ),
        call = call
    )
}


# ---------------------------------------------------------------------------
# Capture, parsing and validation of the conditions.
# ---------------------------------------------------------------------------

# Operators a stopping condition may use, besides parentheses. They are
# interpolated as values in the cli messages: a literal `<` inside cli markup
# would be read as an internal delimiter.
gsd_rule_logic_ops <- c("&&", "||", "!")
gsd_rule_compare_ops <- c("==", "!=", "<", "<=", ">", ">=")
gsd_rule_arith_ops <- c("+", "-", "*", "/")

# Split the captured argument into its conditions. Returns a list of the
# conditions (language objects or scalar literals, not yet validated) and
# `n_interim`: `NULL` for a single condition to be recycled, otherwise the
# number of conditions supplied.
.gsd_rule_capture <- function(cond_lang, call = rlang::caller_env()) {
    items <- .gsd_rule_list_items(cond_lang)
    if (!is.null(items)) {
        if (length(items) == 0L) {
            cli::cli_abort(
                "{.arg condition} must hold at least one condition, not an \\
                empty {.fn list}.",
                call = call
            )
        }
        conditions <- lapply(
            seq_along(items),
            function(j) .gsd_rule_list_item(items[[j]], j, call = call)
        )
        return(list(conditions = conditions, n_interim = length(items)))
    }

    if (is.character(cond_lang)) {
        if (length(cond_lang) == 0L) {
            cli::cli_abort(
                "{.arg condition} must hold at least one condition, not an \\
                empty character vector.",
                call = call
            )
        }
        conditions <- lapply(cond_lang, .gsd_rule_parse_string, call = call)
        n_interim <- if (length(cond_lang) == 1L) NULL else length(cond_lang)
        return(list(conditions = conditions, n_interim = n_interim))
    }

    if (.gsd_rule_is_expr(cond_lang)) {
        return(list(conditions = list(cond_lang), n_interim = NULL))
    }

    cli::cli_abort(
        "{.arg condition} must be an expression, a string, a call to \\
        {.fn list} or a character vector, not \\
        {.obj_type_friendly {cond_lang}}.",
        call = call
    )
}

# The per-interim elements when the captured argument is a list: a call to
# `list(...)`, a call to `c(...)` whose arguments are all single string
# literals, or an injected list value. `NULL` for any other argument (a `c()`
# call with anything else then fails the grammar, as an unsupported call).
.gsd_rule_list_items <- function(cond_lang) {
    if (is.list(cond_lang)) {
        return(cond_lang)
    }
    if (!is.call(cond_lang)) {
        return(NULL)
    }

    items <- as.list(cond_lang)[-1L]
    if (identical(cond_lang[[1L]], quote(list))) {
        return(items)
    }
    is_string <- function(x) is.character(x) && length(x) == 1L
    if (
        identical(cond_lang[[1L]], quote(c)) &&
            length(items) > 0L &&
            all(vapply(items, is_string, logical(1L)))
    ) {
        return(items)
    }

    NULL
}

# An element of a `list(...)` call: an expression or a single string.
.gsd_rule_list_item <- function(item, j, call) {
    if (is.character(item) && length(item) == 1L) {
        return(.gsd_rule_parse_string(item, call = call))
    }
    if (.gsd_rule_is_expr(item)) {
        return(item)
    }

    cli::cli_abort(
        "Element {j} of {.arg condition} must be an expression or a single \\
        string, not {.obj_type_friendly {item}}.",
        call = call
    )
}

# Something the grammar validator can judge: a call, a symbol, or a numeric
# or logical literal (a vector literal is refused there, with a hint).
.gsd_rule_is_expr <- function(x) {
    is.call(x) || is.symbol(x) || is.numeric(x) || is.logical(x)
}

# Parse one string after mapping the word operators to R's. Precedence is
# then R's own, as for `trial_success_gsd()`.
.gsd_rule_parse_string <- function(text, call) {
    if (is.na(text)) {
        cli::cli_abort(
            "A stopping condition must be a string, not {.val {NA}}.",
            call = call
        )
    }

    fixed <- gsub("\\b[Aa][Nn][Dd]\\b", "&&", text)
    fixed <- gsub("\\b[Oo][Rr]\\b", "||", fixed)
    fixed <- gsub("\\b[Nn][Oo][Tt]\\b", "!", fixed)

    tryCatch(
        str2lang(fixed),
        error = function(cnd) {
            cli::cli_abort(
                "Could not parse the stopping condition {.str {text}}.",
                parent = cnd,
                call = call
            )
        }
    )
}

#' Validate that a stopping condition uses only the grammar
#'
#' Stopping-rule twin of `validate_expr_symbols_gsd()`. After `!!` unquoting,
#' a condition may contain only `r<digit>` symbols, numeric or logical
#' scalars, the logical (`&&`, `||`, `!`), comparison and arithmetic
#' operators, and parentheses.
#'
#' @param expr A language object (post-unquoting) or a scalar literal.
#' @returns Invisible `NULL` on success; errors with a helpful message
#' otherwise.
#' @noRd
.gsd_rule_validate <- function(expr, call = rlang::caller_env()) {
    if (is.call(expr)) {
        .gsd_rule_validate_call(expr, call)
    } else if (is.symbol(expr)) {
        .gsd_rule_validate_symbol(expr, call)
    } else if (is.numeric(expr) || is.logical(expr)) {
        if (length(expr) != 1L) {
            cli::cli_abort(
                c(
                    "A literal in a stopping condition must be a single \\
                    number or logical value, not a vector of length \\
                    {length(expr)}.",
                    i = "Give one condition per interim analysis in a \\
                    {.fn list} instead of injecting a vector of values."
                ),
                call = call
            )
        }
    } else {
        cli::cli_abort(
            "Unexpected element of type {.code {typeof(expr)}} in a stopping \\
            condition.",
            call = call
        )
    }
    invisible(NULL)
}

# A call must be one of the allowed operators, applied to sub-expressions
# that are validated in turn.
.gsd_rule_validate_call <- function(expr, call) {
    fn <- expr[[1L]]
    fn_txt <- if (is.symbol(fn)) as.character(fn) else deparse1(fn)
    allowed <- c(
        gsd_rule_logic_ops,
        gsd_rule_compare_ops,
        gsd_rule_arith_ops,
        "("
    )

    if (!is.symbol(fn) || !(fn_txt %in% allowed)) {
        cli::cli_abort(
            c(
                "Unsupported operator or function {.code {fn_txt}} in a \\
                stopping condition.",
                i = "Only logical ({.code {gsd_rule_logic_ops}}), comparison \\
                ({.code {gsd_rule_compare_ops}}) and arithmetic \\
                ({.code {gsd_rule_arith_ops}}) operators and parentheses are \\
                allowed.",
                "*" = "There are no discount tables in a stopping rule."
            ),
            call = call
        )
    }

    for (arg in as.list(expr)[-1L]) {
        .gsd_rule_validate(arg, call = call)
    }
    invisible(NULL)
}

# A bare symbol must be `r<i>`. The decision times of the gain grammar get
# their own message; anything else is an un-injected variable.
.gsd_rule_validate_symbol <- function(expr, call) {
    txt <- as.character(expr)
    if (grepl("^r[1-9][0-9]*$", txt)) {
        return(invisible(NULL))
    }
    if (grepl("^t[0-9]+$", txt)) {
        cli::cli_abort(
            c(
                "Symbol {.code {txt}} cannot be used in a stopping condition: \\
                decision times are not available to a stopping rule.",
                i = "A condition is evaluated at each interim analysis on the \\
                rejection indicators {.code r1}, {.code r2}, ..."
            ),
            call = call
        )
    }
    cli::cli_abort(
        c(
            "Symbol {.arg {txt}} in a stopping condition is not a rejection \\
            indicator ({.code r1}, {.code r2}, ...).",
            i = "Use {.code !!} to inject values from the calling \\
            environment, e.g. {.code stopping_rule_gsd(r1 + r2 >= !!{txt})}."
        ),
        call = call
    )
}

# The hypothesis indices a validated condition refers to, aborting when it
# refers to none (a condition that cannot depend on the trial).
.gsd_rule_indices <- function(expr, label, call = rlang::caller_env()) {
    syms <- grep("^r[1-9][0-9]*$", all.names(expr), value = TRUE)
    if (length(syms) == 0L) {
        cli::cli_abort(
            c(
                "Stopping condition {.code {label}} must refer to at least \\
                one rejection indicator ({.code r1}, {.code r2}, ...).",
                i = "A condition that refers to none cannot depend on the \\
                trial."
            ),
            call = call
        )
    }
    as.integer(sub("^r", "", syms))
}
