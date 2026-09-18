#' Create a trial success function for a group sequential design
#'
#' Create a user-defined **trial-success utility** \eqn{\psi} for a group
#' sequential design. It assigns value to each simulated trial from the
#' hypotheses that were rejected *and the analysis at which each rejection was
#' declared*, so that an early rejection can be worth more than a late one. The
#' function compiles \eqn{\psi} to C++ for simulation and optimisation, exactly
#' as [trial_success()] does for a fixed-sample design.
#'
#' In code, \eqn{\psi} is written using:
#'
#' * `r1, r2, ..., rm`: the binary indicator that hypothesis \eqn{H_i} is
#'   rejected at any analysis;
#' * `t1, t2, ..., tm`: the analysis ("look") at which \eqn{H_i} was rejected,
#'   an integer from `1` to `K`, or `0` if it was never rejected. `ri` is the
#'   same as `ti > 0`;
#' * arithmetic (`+ - * /`), logical operators (`&&`, `||`, or the words `and`,
#'   `or` in a string) between logical terms, comparisons (`== != < <= > >=`),
#'   and parentheses;
#' * discount tables: any numeric vector of length `K` passed by name through
#'   `...` can be applied to a decision time as `name(ti)`. The result is the
#'   `ti`-th element of the table, or `0` when `ti` is `0` (never rejected).
#'
#' @param objective An expression or string encoding the trial-success utility
#'   \eqn{\psi} over `r1, r2, ...` and `t1, t2, ...`. To inject values from
#'   your R session, use rlang's injection operator -
#'   [`!!`][rlang::injection-operator] - (see Examples). Use `!!` only inside
#'   `objective`; it has a different meaning in ordinary arguments (see
#'   **Common mistakes**).
#' @param ... Discount tables, each a named numeric vector of length `K` (one
#'   value per analysis). Supply only the values for analyses 1 through `K`;
#'   the value for "never rejected" (`ti = 0`) is always `0` and is supplied
#'   automatically — do not include it yourself (see **Common mistakes**).
#'   Every table must be named; the names are used as functions in `objective`.
#' @param K An optional whole number giving the number of analyses. When
#'   discount tables are supplied `K` is inferred from their common length and,
#'   if also given, must agree with it. When neither is given `K` is left
#'   `NULL` and checked against the p-values at first use.
#' @inheritParams trial_success
#'
#' @returns A `multigrain_trial_success_gsd` object (which also inherits from
#'   `multigrain_trial_success`) made up of:
#'   * `func`: compiled function that evaluates \eqn{\psi} row-wise on an
#'     integer matrix of decision times (simulations by hypotheses, `0` for
#'     never rejected) and returns the mean utility.
#'   * `m`: number of hypotheses implied by the largest index among the
#'     `r<i>` and `t<i>` symbols.
#'   * `K`: number of analyses, or `NULL` if unknown.
#'   * `objective`: the original utility expression, as a display string
#'     (numbers shown to 15 significant digits; the compiled function uses
#'     the exact values).
#'   * `cpp_code`: the generated C++ source.
#'   * `tables`: the discount tables, as supplied.
#'
#' @details
#'
#'   ## Decision times
#'
#'   Each `ti` records the analysis at which the graphical procedure *declared*
#'   hypothesis \eqn{H_i} rejected. This may be later than the analysis at
#'   which the data for \eqn{H_i} matured or first crossed a boundary. For
#'   example, suppose PFS data are final at analysis 1 but PFS carries zero
#'   initial weight. PFS cannot reject until OS rejects at analysis 2 and
#'   recycles its alpha to PFS. The result is `t_PFS = 2`, even though the PFS
#'   p-value was available at analysis 1.
#'
#'   The rejection indicator `ri` is derived from `ti > 0` inside the compiled
#'   function, so `ri` and `ti` can never disagree: if `ti = 0` then
#'   `ri = FALSE`, and if `ti >= 1` then `ri = TRUE`.
#'
#'   ## Discount tables
#'
#'   A discount table maps each decision time to a value. You supply one value
#'   per analysis, for analyses 1 through `K`. The constructor automatically
#'   prepends `0.0` for `ti = 0` (never rejected), so you should **not**
#'   include a leading zero yourself. For a two-analysis design where a
#'   rejection at the final analysis is worth 75% of one at the interim, write
#'   `d = c(1, 0.75)`. Element names on the vector are ignored; only the
#'   argument name matters (the name you use to apply it in `objective`, e.g.
#'   `d(t1)`).
#'
#'   ## Operator precedence
#'
#'   Operator precedence follows R's own rules (unlike [trial_success()],
#'   whose expression-input parser treats `&&` and `||` at a different
#'   precedence level). Always parenthesise to make intent explicit: write
#'   `(t1 == 1 && t2 == 1)` rather than relying on `&&` binding tighter
#'   than `+`.
#'
#' @section Common mistakes:
#'
#'   **Including a leading zero in a discount table.** Writing
#'   `d = c(0, 1, 0.75)` does not give a two-analysis table — it gives a
#'   *three*-analysis table with `K = 3`, where analysis 1 has value 0 and
#'   analysis 2 has value 1. Every claim is mispriced and there is no warning.
#'   For `K = 2`, write `d = c(1, 0.75)`.
#'
#'   **Using `!!` inside a discount table argument.** The injection operator
#'   [`!!`][rlang::injection-operator] is only meaningful inside `objective`,
#'   where `rlang::enexpr()` captures the expression. In an ordinary argument
#'   like `d = c(1, !!delta)`, `!!` is R's double negation: `!!0.75` evaluates
#'   to `TRUE` (i.e. `1`), so the table silently becomes `c(1, 1)` and the
#'   discount is lost. Write `d = c(1, delta)` instead.
#'
#'   **Passing a GSD object to the fixed-sample optimiser.** Because
#'   `multigrain_trial_success_gsd` inherits from `multigrain_trial_success`,
#'   `graph_optimise(trial_success = <gsd_object>)` passes validation. However,
#'   the compiled function then receives a logical rejection matrix (not a
#'   decision-time matrix). Rcpp silently coerces `TRUE`/`FALSE` to `1`/`0`,
#'   so every rejection is scored as though it occurred at analysis 1 — a
#'   wrong answer with no error.
#'
#' @section Interpreting the gain:
#'
#'   The gain returned by `$func()` is the mean of \eqn{\psi} over simulated
#'   trials and is **not** a probability. Its scale depends entirely on how you
#'   define `objective`. For instance, \eqn{\psi = r_1 + r_2} at 60% power
#'   per hypothesis returns roughly 1.2, not 0.6. Only relative values matter
#'   for optimisation (the argmax is invariant to positive affine
#'   transformations of \eqn{\psi}), but normalise for reporting if you want
#'   interpretable numbers.
#'
#'   Non-monotone discount tables (e.g. `d = c(0.5, 1)`, which values a late
#'   decision more than an early one) compile without complaint. This is by
#'   design — the package does not enforce a particular preference ordering —
#'   but check that the table reflects your intent.
#'
#' @note The public GSD consumer (`calc_power_pvals_gsd()`) does not exist
#'   yet. The gain function can currently only be evaluated via the internal
#'   GSD kernel.
#'
#' @seealso [trial_success()] for the fixed-sample version.
#'
#' @export
#' @examples
#' # Manuscript Example 5: PFS (H1) and OS (H2), two analyses; a rejection at
#' # the final analysis is worth 75% of one at the interim.
#' # d = c(1, 0.75): full value at analysis 1, 75% at analysis 2.
#' # The zero for "never rejected" is supplied automatically.
#' v_pfs <- 0.4
#' v_os  <- 0.6
#' gain <- trial_success_gsd(
#'     !!v_pfs * d(t1) + !!v_os * d(t2),
#'     d = c(1, 0.75)
#' )
#' gain
#'
#' # The gain on every possible (t1, t2) pair for K = 2:
#' #
#' #   t1 \ t2 |    0      1      2
#' #   --------+---------------------
#' #      0    | 0.000  0.600  0.450
#' #      1    | 0.400  1.000  0.850
#' #      2    | 0.300  0.900  0.750
#'
#' \donttest{
#' # Mixed r and t: OS value plus a bonus only if PFS was declared at
#' # the interim
#' trial_success_gsd(r2 + 0.5 * (r1 && t1 == 1), K = 2)
#'
#' # Hurdle gain with a time-discounted base value: both claims required,
#' # and the package is worth less if the confirmatory claim is late
#' trial_success_gsd((r1 && r2) * d(t2), d = c(1, 0.75))
#'
#' # Co-primary endpoints rejected at the same analysis
#' trial_success_gsd((t1 == 1 && t2 == 1) + 0.5 * (t1 == 2 && t2 == 2))
#'
#' # Rejection indicators work as in trial_success()
#' trial_success_gsd(r1 + r2, K = 2)
#' }
trial_success_gsd <- function(
    objective,
    ...,
    K = NULL, # nolint: object_name_linter. K is the record's symbol
    verbose = multigrain_verbosity()
) {
    # Call flow: this is the only exported entry point. It normalises
    # `verbose`, checks the discount tables passed through `...`, settles `K`,
    # turns the captured expression into a string with `resolve_expr_gsd()`,
    # and hands everything to `new_trial_success_gsd()`, which generates and
    # compiles the C++.
    if (isTRUE(verbose)) {
        verbose <- "info"
    } else if (isFALSE(verbose)) {
        verbose <- "silent"
    } else {
        rlang::check_string(verbose)
        verbose <- rlang::arg_match(verbose, values = verbosity_levels)
    }

    tables <- .gsd_gain_tables(rlang::list2(...))
    K <- .gsd_gain_K(tables, K) # nolint: object_name_linter.

    expr_lang <- rlang::enexpr(objective)
    expr_string <- resolve_expr_gsd(expr_lang, table_names = names(tables))
    new_trial_success_gsd(
        expr_string,
        # The captured expression (when not a string) is compiled directly so
        # that injected constants keep full precision; the deparsed string is
        # what the object displays.
        expr_lang = if (is.language(expr_lang)) expr_lang else NULL,
        tables = tables,
        K = K,
        verbose = verbose
    )
}


#' @export
print.multigrain_trial_success_gsd <- function(x, ...) {
    if (is.null(x)) {
        return()
    }

    cli::cat_line(cli::format_inline("{.cls {class(x)}}"))
    cli::cat_line(x$objective)
    .gsd_gain_cat_detail(x)
    invisible(x)
}

#' @export
summary.multigrain_trial_success_gsd <- function(object, ...) {
    if (is.null(object)) {
        return()
    }

    cli::cat_line(
        cli::style_underline("\nTrial success function (group sequential)"),
        ":"
    )
    cli::cat_line(object$objective)
    .gsd_gain_cat_detail(object)
    invisible(object)
}

# Shared tail of print() and summary(): the number of analyses and the tables.
.gsd_gain_cat_detail <- function(x) {
    cli::cat_line(
        "Analyses (K): ",
        if (is.null(x$K)) "not set" else format(x$K)
    )
    for (nm in names(x$tables)) {
        cli::cat_line(
            "Discount table ", nm, "(t): ",
            toString(format(x$tables[[nm]]))
        )
    }
    invisible(NULL)
}

is_trial_success_gsd <- function(x) {
    inherits(x, "multigrain_trial_success_gsd")
}


# ---------------------------------------------------------------------------
# Argument helpers, called by trial_success_gsd() before the expression is
# looked at.
# ---------------------------------------------------------------------------

# Validate the discount tables collected from `...`. Returns a named list of
# finite numeric vectors that all share one length (the number of analyses).
# Names must be usable as C++ identifiers (they become `<name>_tab`) and must
# not collide with the symbol grammar (`r<i>`, `t<i>`) or the word operators
# (`and`, `or`), because those are recognised by pattern before table names.
.gsd_gain_tables <- function(tables, call = rlang::caller_env()) {
    if (length(tables) == 0L) {
        return(rlang::set_names(list(), character()))
    }

    nms <- names(tables)
    if (is.null(nms) || !all(nzchar(nms))) {
        cli::cli_abort(
            c(
                "Every discount table passed through {.arg ...} must be named.",
                i = "The name is how the table is applied in \\
                {.arg objective}, e.g. {.code d = c(1, 0.75)} for \\
                {.code d(t1)}."
            ),
            call = call
        )
    }
    if (anyDuplicated(nms) > 0L) {
        cli::cli_abort(
            "Discount table names must be unique; \\
            {.val {nms[duplicated(nms)]}} appear{?s/} more than once.",
            call = call
        )
    }

    # nolint start: nonportable_path_linter
    bad_name <- !grepl("^[A-Za-z][A-Za-z0-9_]*$", nms) |
        grepl("^[rt][0-9]+$", nms) |
        tolower(nms) %in% c("and", "or")
    # nolint end
    if (any(bad_name)) {
        cli::cli_abort(
            c(
                "Discount table name{?s} {.val {nms[bad_name]}} {?is/are} not \\
                allowed.",
                i = "A name must start with a letter, contain only letters, \\
                digits and underscores, and must not be a rejection indicator \\
                ({.code r1}, ...), a decision time ({.code t1}, ...), \\
                {.code and} or {.code or}."
            ),
            call = call
        )
    }

    for (nm in nms) {
        tab <- tables[[nm]]
        if (
            !is.numeric(tab) ||
                length(tab) == 0L ||
                !all(is.finite(tab))
        ) {
            cli::cli_abort(
                "Discount table {.arg {nm}} must be a non-empty numeric \\
                vector of finite values, not {.obj_type_friendly {tab}}.",
                call = call
            )
        }
    }

    lens <- lengths(tables)
    if (length(unique(lens)) != 1L) {
        cli::cli_abort(
            c(
                "All discount tables must have the same length (one value per \\
                analysis).",
                x = "Lengths supplied: \\
                {paste0(nms, ' = ', lens, collapse = ', ')}."
            ),
            call = call
        )
    }

    lapply(tables, as.double)
}

# Settle the number of analyses `K`: the common table length when tables are
# given (checked against an explicit `K`), otherwise the explicit `K`, or
# NULL when neither is available (design record, section 10 item 6).
.gsd_gain_K <- function( # nolint: object_name_linter.
    tables,
    K, # nolint: object_name_linter.
    call = rlang::caller_env()
) {
    if (!is.null(K)) {
        rlang::check_number_whole(K, min = 1, call = call)
        K <- as.integer(K) # nolint: object_name_linter.
    }
    if (length(tables) == 0L) {
        return(K)
    }

    k_tab <- length(tables[[1L]])
    if (!is.null(K) && K != k_tab) {
        cli::cli_abort(
            "{.arg K} is {K} but the {cli::qty(length(tables))}discount \\
            table{?s} {.arg {names(tables)}} {?has/have} length {k_tab}; \\
            they must agree.",
            call = call
        )
    }
    as.integer(k_tab)
}


# ---------------------------------------------------------------------------
# Expression capture and validation.
# ---------------------------------------------------------------------------

#' Resolve a captured GSD trial success expression to a string
#'
#' Group sequential twin of `resolve_expr()`. Takes a language object (from
#' `rlang::enexpr()`) or a string. A language object is validated with
#' `validate_expr_symbols_gsd()` and deparsed; a string passes through and is
#' validated later, after parsing, in `replace_indices_gsd()`.
#'
#' @param expr_lang A language object or a character string.
#' @param table_names Names of the discount tables supplied by the user.
#' @returns A length-1 character string of the resolved expression.
#'
#' @noRd
resolve_expr_gsd <- function(
    expr_lang,
    table_names = character(),
    call = rlang::caller_env()
) {
    if (is.character(expr_lang)) {
        if (length(expr_lang) != 1L) {
            cli::cli_abort(
                "{.arg objective} must be a single string, not a character \\
                vector of length {length(expr_lang)}.",
                call = call
            )
        }
        expr_string <- expr_lang
    } else if (is.language(expr_lang)) {
        validate_expr_symbols_gsd(expr_lang, table_names = table_names)
        expr_string <- deparse1(expr_lang, width.cutoff = 500)
    } else {
        cli::cli_abort(
            "{.arg objective} must be an expression or a character string, \\
            not {.obj_type_friendly {expr_lang}}.",
            call = call
        )
    }

    expr_string
}


# Operators the GSD grammar accepts, besides the discount-table calls.
gsd_gain_arith_ops <- c("+", "-", "*", "/")
gsd_gain_logic_ops <- c("&&", "||")
gsd_gain_compare_ops <- c("==", "!=", "<", "<=", ">", ">=")

#' Validate that a GSD trial success expression uses only the grammar
#'
#' After `!!` unquoting, the expression may contain only `r<digit>` and
#' `t<digit>` symbols, numeric or logical scalars, the arithmetic, logical and
#' comparison operators, parentheses, and calls of a discount table to a single
#' `t<digit>` symbol. Any other symbol means the user forgot to unquote.
#'
#' @param expr A language object (post-unquoting).
#' @param table_names Names of the discount tables supplied by the user.
#' @returns Invisible `NULL` on success; errors with a helpful message
#' otherwise.
#' @noRd
validate_expr_symbols_gsd <- function(expr, table_names = character()) {
    if (is.call(expr)) {
        .gsd_gain_validate_call(expr, table_names)
    } else if (is.symbol(expr)) {
        .gsd_gain_validate_symbol(expr, table_names)
    } else if (is.numeric(expr) || is.logical(expr)) {
        if (length(expr) != 1L) {
            cli::cli_abort(
                c(
                    "A literal in the trial success expression must be a \\
                    single number, not a vector of length {length(expr)}.",
                    i = "Pass a vector of per-analysis values as a named \\
                    discount table through {.arg ...} and apply it as \\
                    {.code name(t1)}."
                )
            )
        }
    } else {
        cli::cli_abort(
            "Unexpected element of type {.code {typeof(expr)}} in trial \\
            success expression."
        )
    }
    invisible(NULL)
}

# A call is either a discount table applied to one `t<i>` symbol, or one of
# the allowed operators applied to sub-expressions that are validated in turn.
.gsd_gain_validate_call <- function(expr, table_names) {
    fn <- expr[[1L]]
    if (!is.symbol(fn)) {
        cli::cli_abort(
            "Unsupported call {.code {deparse1(expr)}} in the trial success \\
            expression."
        )
    }
    fn_txt <- as.character(fn)

    if (fn_txt %in% table_names) {
        tab_args <- as.list(expr)[-1L]
        # nolint start: nonportable_path_linter
        ok <- length(tab_args) == 1L &&
            is.symbol(tab_args[[1L]]) &&
            grepl("^t\\d+$", as.character(tab_args[[1L]]))
        # nolint end
        if (!ok) {
            cli::cli_abort(
                c(
                    "A discount table must be applied to a single decision \\
                    time symbol ({.code t1}, {.code t2}, ...), \\
                    not as {.code {deparse1(expr)}}.",
                    i = "For example {.code {fn_txt}(t1)}."
                )
            )
        }
        return(invisible(NULL))
    }

    allowed <- c(
        gsd_gain_arith_ops,
        "(",
        gsd_gain_logic_ops,
        gsd_gain_compare_ops
    )
    if (!(fn_txt %in% allowed)) {
        # The operator lists are interpolated as values: a literal `<` inside
        # cli markup would be read as an internal delimiter.
        # nolint start: object_usage_linter. used inside the cli message
        arith <- gsd_gain_arith_ops
        logic <- gsd_gain_logic_ops
        compare <- gsd_gain_compare_ops
        # nolint end
        cli::cli_abort(
            c(
                "Unsupported operator or function {.code {fn_txt}} in the \\
                trial success expression.",
                i = "Only arithmetic ({.code {arith}}), logical \\
                ({.code {logic}}), comparison ({.code {compare}}), \\
                parentheses and the discount tables passed through \\
                {.arg ...} are allowed."
            )
        )
    }
    lapply(
        as.list(expr)[-1L],
        validate_expr_symbols_gsd,
        table_names = table_names
    )
    invisible(NULL)
}

# A bare symbol must be `r<i>` or `t<i>`. A table name on its own (without
# `(t1)`) gets its own hint; anything else is an un-injected variable.
.gsd_gain_validate_symbol <- function(expr, table_names) {
    txt <- as.character(expr)
    # nolint start: nonportable_path_linter, line_length_linter
    if (grepl("^[rt]\\d+$", txt)) {
        return(invisible(NULL))
    }
    if (txt %in% table_names) {
        cli::cli_abort(
            c(
                "Discount table {.arg {txt}} must be applied to a decision \\
                time, e.g. {.code {txt}(t1)}; it cannot stand alone.",
                i = "Use {.code !!} if you meant to inject a value instead."
            )
        )
    }
    cli::cli_abort(
        c(
            "Symbol {.arg {txt}} in the trial success expression is not \\
            a rejection indicator ({.code r1}, {.code r2}, ...), a decision \\
            time ({.code t1}, {.code t2}, ...), or a discount table.",
            i = "Use {.code !!} to inject values from the calling \\
            environment, e.g. {.code trial_success_gsd(!!{txt} * r1 + r2)}."
        )
    )
    # nolint end
}


# ---------------------------------------------------------------------------
# From expression string to compiled C++.
# ---------------------------------------------------------------------------

# Constructor. Called by trial_success_gsd() with a validated expression
# string, validated tables and a settled K. Counts the hypotheses, translates
# the expression to a C++ body, wraps it in a row loop over the decision-time
# matrix, compiles it with Rcpp::sourceCpp() and returns the S3 object.
new_trial_success_gsd <- function(
    expr_string,
    expr_lang = NULL,
    tables = rlang::set_names(list(), character()),
    K = NULL, # nolint: object_name_linter.
    verbose = c("info", "detail", "silent")
) {
    verbose <- rlang::arg_match(verbose)
    local_env <- new.env()

    # Number of hypotheses: the largest index over r<i> and t<i>
    m <- count_unique_indices_gsd(expr_string)

    # Expression body in C++ over the decision-time matrix `t`, from the
    # captured expression when there is one (exact constants), else from the
    # string
    cpp_body <- replace_indices_gsd(
        expr_lang %||% expr_string,
        table_names = names(tables)
    )

    # One static array per discount table. Index 0 is "never rejected" and
    # maps to 0; index k (1..K) is the table value at analysis k. Values are
    # written with enough digits to round-trip exactly.
    table_lines <- vapply(
        names(tables),
        function(nm) {
            sprintf(
                "static const double %s_tab[] = {%s};",
                nm,
                toString(c(
                    "0.0",
                    vapply(tables[[nm]], .gsd_gain_cpp_number, character(1L))
                ))
            )
        },
        character(1L)
    )
    table_block <- if (length(table_lines) == 0L) {
        ""
    } else {
        paste0(paste(table_lines, collapse = "\n"), "\n")
    }

    # nolint start: quotes_linter
    cpp_code <- sprintf(
        '
#include <Rcpp.h>
using namespace Rcpp;

#define std_min std::min

%s
// [[Rcpp::export]]
double %s(IntegerMatrix t) {
    if (t.ncol() < %d) {
        stop("the decision-time matrix has %%d column(s) but the trial success function refers to %d hypotheses", t.ncol());
    }
    int n = t.nrow();
    double total = 0.0;

    for (int i = 0; i < n; i++) {
        total += (%s); // gain of simulated trial i from its decision times
    }

    return total / n; // Return the mean
}',
        table_block,
        "powerFunc",
        m,
        m,
        cpp_body
    )
    # nolint end

    sourceCpp(code = cpp_code, env = local_env)

    if (verbose != "silent") {
        cli::cli_alert_success(
            "Trial success function compiled and sourced successfully."
        )
    }

    structure(
        list(
            func = local_env$powerFunc,
            m = m,
            K = K,
            objective = expr_string,
            cpp_code = cpp_code,
            tables = tables
        ),
        class = c("multigrain_trial_success_gsd", "multigrain_trial_success")
    )
}


# Largest hypothesis index referenced by an r<i> or t<i> symbol, with a
# warning when the sequence 1..max has gaps (as count_unique_indices()).
count_unique_indices_gsd <- function(expr_string) {
    # nolint start: nonportable_path_linter
    matches <- regmatches(
        expr_string,
        gregexpr("\\b[rt][0-9]+\\b", expr_string)
    )
    # nolint end
    all_matches <- unlist(matches)

    if (length(all_matches) == 0L) {
        cli::cli_abort(
            "{.arg objective} must reference at least one rejection \\
            indicator ({.code r1}, {.code r2}, ...) or decision time \\
            ({.code t1}, {.code t2}, ...). \\
            Expression {.code {expr_string}} contains none."
        )
    }

    numeric_indices <- as.integer(sub("^[rt]", "", all_matches))
    unique_indices <- unique(numeric_indices)
    max_index_found <- max(unique_indices, na.rm = TRUE)
    missing_indices <- setdiff(seq_len(max_index_found), unique_indices)

    if (length(missing_indices) > 0L) {
        cli::cli_warn(
            "Missing indices in the sequence. Expected every index from \\
            1 to {max_index_found}, but missing {toString(missing_indices)}."
        )
    }

    max_index_found
}


#' Translate a GSD trial success expression to a C++ body
#'
#' Group sequential twin of `replace_r_indices()`. The word operators `and` and
#' `or` (usable in string input only) are mapped to `&&` and `||`, and the
#' string is then parsed by R itself, so operator precedence is R's. (The
#' fixed-sample parser substitutes `%AND%`/`%OR%` placeholders, which R cannot
#' parse next to a comparison such as `t1 == 1`.) The parsed tree is validated
#' and rewritten by `parse_and_transform_gsd()` and deparsed to C++.
#'
#' @param expr A character string with the utility expression, or the
#'   captured language object itself (in which numeric constants are exact,
#'   where a deparsed string carries 15 significant digits).
#' @param table_names Names of the discount tables supplied by the user.
#' @returns A character string of C++ code over the row index `i` and the
#'   `IntegerMatrix t` of decision times: `r1` becomes `double(t(i, 0) > 0)`,
#'   `t1` becomes `double(t(i, 0))`, `d(t1)` becomes `d_tab[t(i, 0)]`.
#'
#' @noRd
replace_indices_gsd <- function(expr, table_names = character()) {
    if (is.character(expr)) {
        fixed_expr <- expr
        fixed_expr <- gsub("\\b[Aa][Nn][Dd]\\b", "&&", fixed_expr)
        fixed_expr <- gsub("\\b[Oo][Rr]\\b", "||", fixed_expr)
        ast <- str2lang(fixed_expr)
    } else {
        ast <- expr
    }
    validate_expr_symbols_gsd(ast, table_names = table_names)

    transformed <- parse_and_transform_gsd(ast, table_names = table_names)

    out_str <- deparse1(transformed$expr)
    out_str <- gsub("`", "", out_str, fixed = TRUE)

    # Insert space after/before slash
    out_str <- gsub("([[:alnum:]_\\)])/", "\\1 /", out_str)
    out_str <- gsub("/([[:alnum:]_\\(])", "/ \\1", out_str)
    out_str
}


# nolint start: return_linter

# Recursive rewrite of one node of the parsed expression into C++ (as an R
# call object that deparses to valid C++), tracking its type: "bool" for a
# 0/1 quantity, "real" otherwise. The type rule is that of the fixed-sample
# parser: && and || need bool operands; bool * bool is bool; every other
# arithmetic result is real; a comparison is bool; a table lookup is real.
parse_and_transform_gsd <- function(node, table_names = character()) {
    if (is.call(node)) {
        return(.gsd_gain_transform_call(node, table_names))
    } else if (is.symbol(node)) {
        return(.gsd_gain_transform_symbol(node))
    } else if (is.numeric(node)) {
        # Numeric literal => real, written as a C++ double literal
        return(list(expr = as.symbol(.gsd_gain_cpp_number(node)), type = "real"))
    } else if (is.logical(node)) {
        # TRUE/FALSE => bool 1.0/0.0
        return(list(
            expr = as.symbol(if (isTRUE(node)) "1.0" else "0.0"),
            type = "bool"
        ))
    } else {
        return(list(expr = node, type = "real"))
    }
}

.gsd_gain_transform_call <- function(node, table_names) {
    parts <- as.list(node)
    fn <- parts[[1L]]
    op_text <- as.character(fn)

    # --- Parentheses: transform the inside and return it directly ---
    if (op_text == "(") {
        return(parse_and_transform_gsd(parts[[2L]], table_names))
    }

    # --- Discount table applied to a decision time: name(t<i>) ---
    # Becomes name_tab[t(i, idx)]; index 0 (never rejected) holds 0.0.
    if (op_text %in% table_names) {
        idx <- as.integer(sub("^t", "", as.character(parts[[2L]]))) - 1
        new_call <- call(
            "[",
            as.symbol(paste0(op_text, "_tab")),
            call("t", quote(i), idx)
        )
        return(list(expr = new_call, type = "real"))
    }

    # --- Logical operators: bool operands only ---
    if (op_text %in% gsd_gain_logic_ops) {
        left <- parse_and_transform_gsd(parts[[2L]], table_names)
        right <- parse_and_transform_gsd(parts[[3L]], table_names)
        if (left$type != "bool" || right$type != "bool") {
            what <- if (op_text == "&&") "`&&` (AND)" else "`||` (OR)"
            stop(what, " only allowed between booleans.", call. = FALSE)
        }
        new_call <- if (op_text == "&&") {
            substitute(A * B, list(A = left$expr, B = right$expr))
        } else {
            bquote(std_min(double(1), .(left$expr) + .(right$expr)))
        }
        return(list(expr = new_call, type = "bool"))
    }

    # --- Comparisons: any operands, result is bool ---
    # Emitted as double(A op B) so the 0/1 result is a double like every
    # other term.
    if (op_text %in% gsd_gain_compare_ops) {
        left <- parse_and_transform_gsd(parts[[2L]], table_names)
        right <- parse_and_transform_gsd(parts[[3L]], table_names)
        new_call <- call("double", as.call(list(fn, left$expr, right$expr)))
        return(list(expr = new_call, type = "bool"))
    }

    # --- Arithmetic: unary (length 2) or binary ---
    transformed_args <- lapply(
        parts[-1L],
        parse_and_transform_gsd,
        table_names = table_names
    )
    new_call <- as.call(c(
        list(fn),
        lapply(transformed_args, `[[`, "expr")
    ))
    if (length(transformed_args) == 1L) {
        # Unary minus or plus => real
        return(list(expr = new_call, type = "real"))
    }
    combo <- combine_arithmetic(
        op_text,
        transformed_args[[1L]]$type,
        transformed_args[[2L]]$type
    )
    return(list(expr = new_call, type = combo$type))
}

.gsd_gain_transform_symbol <- function(node) {
    txt <- as.character(node)
    if (grepl("^r\\d+$", txt)) {
        # r<i> => double(t(i, idx) > 0): rejected at some analysis, bool
        idx <- as.integer(sub("^r", "", txt)) - 1
        return(list(
            expr = call("double", call(">", call("t", quote(i), idx), 0)),
            type = "bool"
        ))
    } else if (grepl("^t\\d+$", txt)) {
        # t<i> => double(t(i, idx)): the decision time itself, real
        idx <- as.integer(sub("^t", "", txt)) - 1
        return(list(
            expr = call("double", call("t", quote(i), idx)),
            type = "real"
        ))
    } else {
        # Cannot happen after validation; kept for symmetry with the
        # fixed-sample parser
        return(list(expr = node, type = "real"))
    }
}

# A numeric scalar as a C++ double literal that reads back to the same
# double: the fewest significant digits (15, 16 or 17) that round-trip, so
# 0.75 stays "0.75" while 1/3 becomes "0.33333333333333331". Integers gain
# ".0"; values already carrying a decimal point or an exponent (e.g. 1e-05)
# are kept as they are.
.gsd_gain_cpp_number <- function(x) {
    txt <- sprintf("%.15g", x)
    for (digits in 16:17) {
        if (as.numeric(txt) == x) {
            break
        }
        txt <- sprintf("%.*g", digits, x)
    }
    if (!grepl("[.eE]", txt)) {
        txt <- paste0(txt, ".0")
    }
    txt
}
# nolint end
