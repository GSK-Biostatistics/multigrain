#' Create a trial success function
#'
#' Create a user-defined **trial-success utility** \eqn{\psi} that assigns
#' value to each rejection pattern from a graphical multiple testing procedure.
#' The function compiles \eqn{\psi} to fast C++ for simulation/optimisation.
#'
#' In code, \eqn{\psi} is written using symbols `r1, r2, ..., rm`, where each
#' `ri` is the binary indicator that hypothesis \eqn{H_i} is rejected by the
#' chosen graph. You can combine these indicators with arithmetic (`+ - * /`)
#' and logical operators (`&&`, `||`, or the words `and`, `or`) to reflect
#' your design priorities (e.g., “any success”, “all successes”, weighted
#' composites).
#'
#' @param objective An expression or string encoding the trial-success utility
#'   \eqn{\psi}. The symbols `r1, r2, ...` refer to rejection indicators for
#'   the corresponding hypotheses. To inject values from your R session, use
#'   rlang's unquote operator `!!` (see Examples). Arithmetic and logical
#'   operators are allowed.
#'
#' @param verbose A logical controlling how much information to print. Defaults
#'   to `TRUE`.
#'
#' @return An object of class `multigrain_trial_success` with components:
#'   * `func`: compiled function that evaluates \eqn{\psi} row-wise on a matrix
#'     of rejection indicators and returns the mean utility (i.e., expected
#'     trial success under the simulated scenario).
#'   * `m`: number of hypotheses implied by `r1, ..., rm`.
#'   * `objective`: the original utility expression (as a string).
#'
#' @details The utility \eqn{\psi} is evaluated on each simulated rejection
#'   pattern and averaged, yielding the **expected trial success** for a given
#'   graph and data-generating scenario. This lets you optimise graphs against
#'   the utility that captures your clinical/regulatory goals, rather than a
#'   single power summary.
#'
#' @export
#' @examples
#'
#' # Expected number of rejections (equals m × average power)
#' exp_rejs <- trial_success(r1 + r2 + r3 + r4)
#'
#' \donttest{
#' # Disjunctive success: any rejection among four
#' disj_power <- trial_success(r1 || r2 || r3 || r4)
#'
#' # Conjunctive success: all four must be rejected
#' conj_power <- trial_success(r1 && r2 && r3 && r4)
#'
#' # Composite utility: require H1 AND H2 rejection to get H3 and H4 value
#' composite <- trial_success((r1 && r2) * (r3 + r4))
#'
#' # Weighted priorities (e.g., H1 gets weight 2 times H2 or H3)
#' weighted <- trial_success(2*r1 + r2 + r3)
#'
#' # Inject values from your environment with !!
#' w <- 2
#' trial_success(!!w * r1 + r2 + r3)
#'
#' # Programmatic string input
#' expr_str <- sprintf("%s * r1 + r2", w)
#' trial_success(!!expr_str)
#'
#' }
trial_success <- function(objective, verbose = TRUE) {
    check_logical(verbose, allow_na = FALSE)
    expr_string <- resolve_expr(rlang::enexpr(objective))
    new_trial_success(expr_string, verbose = verbose)
}


#' @export
print.multigrain_trial_success <- function(x, ...) {
    if (is.null(x)) {
        return()
    }

    cli::cat_line(cli::format_inline("{.cls {class(x)}}"))

    cli::cat_line(x$objective)
    invisible(x)
}

#' @export
summary.multigrain_trial_success <- function(object, ...) {
    if (is.null(object)) {
        return()
    }

    cli::cat_line(cli::style_underline("\nTrial success function"), ":")

    cli::cat_line(object$objective)
    invisible(object)
}

#' Resolve a captured trial success expression to a string
#'
#' Takes a language object (from `rlang::enexpr()`) or a string. If a language
#' object, validates that it contains only `r<digit>` symbols, numeric
#' literals, and operators, then deparses it. If a string, passes it through.
#'
#' @param expr_lang A language object or a character string.
#' @return A length-1 character string of the resolved expression.
#'
#' @noRd
resolve_expr <- function(expr_lang, call = rlang::caller_env()) {
    if (is.character(expr_lang)) {
        expr_string <- expr_lang
    } else if (is.language(expr_lang)) {
        validate_expr_symbols(expr_lang)
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


#' Validate that a trial success expression contains only r<digit> symbols
#'
#' After rlang's `!!` unquoting has resolved user variables, the expression
#' should contain only `r<digit>` symbols, numeric literals, and operators.
#' Any other symbol means the user forgot to unquote.
#'
#' @param expr A language object (post-unquoting).
#' @return Invisible `NULL` on success; errors with a helpful message otherwise.
#' @noRd
validate_expr_symbols <- function(expr) {
    if (is.call(expr)) {
        fn <- expr[[1]]
        allowed <- c("+", "-", "*", "/", "(", "&&", "||")
        if (
            is.symbol(fn) &&
                !(as.character(fn) %in% allowed)
        ) {
            # nolint start: line_length_linter
            cli::cli_abort(
                c(
                    "Unsupported operator or function {.code {as.character(fn)}} \\
                    in the trial success expression.",
                    i = "Only arithmetic ({.code +}, {.code -}, {.code *}, \\
                    {.code /}), logical ({.code &&}, {.code ||}), and \\
                    parentheses are allowed."
                )
            )
            # nolint end
        }
        lapply(expr[-1], validate_expr_symbols)
    } else if (is.symbol(expr)) {
        txt <- as.character(expr)
        # nolint start: nonportable_path_linter, line_length_linter
        if (!grepl("^r\\d+$", txt)) {
            cli::cli_abort(
                c(
                    "Symbol {.arg {txt}} in the trial success expression is not \\
                    a rejection indicator ({.code r1}, {.code r2}, ...).",
                    i = "Use {.code !!} to inject values from the calling \\
                    environment, e.g. {.code trial_success(!!{txt} * r1 + r2)}."
                )
            )
        }
        # nolint end
    } else if (is.numeric(expr) || is.logical(expr)) {} else {
        cli::cli_abort(
            "Unexpected element of type {.code {typeof(expr)}} in trial \\
            success expression."
        )
    }
    invisible(NULL)
}

# Function to create new trial success objective function
new_trial_success <- function(expr_string, verbose = TRUE) {
    local_env <- new.env()

    # Check `objective` and find maximum index (check number of hypotheses)
    m <- count_unique_indices(expr_string)

    # Convert R indexing in the string to C++ indexing handling up to r99
    cpp_body <- replace_r_indices(expr_string)

    # Create the C++ function string
    # nolint start: quotes_linter
    cpp_code <- sprintf(
        '
#include <Rcpp.h>
using namespace Rcpp;

#define std_min std::min

// [[Rcpp::export]]
double %s(LogicalMatrix x) {
    int n = x.nrow();
    double total = 0.0;

    for (int i = 0; i < n; i++) {
        total += (%s); // Apply user-defined function on each row
    }

    return total / n; // Return the mean
}',
        "powerFunc",
        cpp_body
    )
    # nolint end

    Rcpp::sourceCpp(code = cpp_code, env = local_env)

    if (verbose) {
        cli::cli_alert_success(
            "Trial success function compiled and sourced successfully."
        )
    }

    structure(
        list(
            func = local_env$powerFunc,
            m = m,
            objective = expr_string,
            cpp_code = cpp_code
        ),
        class = "multigrain_trial_success"
    )
}


#' Replace R-style indexing with C++ indexing and handle boolean operators
#'
#' Transform R-style indexing for the rejection matrix (i.e., `r1, r2, ..., rm`)
#' into C++-style indexing for use in Rcpp functions. `replace_r_indices()` also
#' processes logical operators and converts both lower and upper case "and"/"or"
#' operators to appropriate C++ logical expressions.
#'
#' @param expr_string A character string containing an expression with R-style
#'   indexing and logical operators. The function will convert the indexing to
#'   C++-style and handle logical operators: "AND", "and", "OR", and "or".
#'
#' @return A modified character string with C++-style indexing and C++ logical
#'   operators. R-style indices like `r1, r2, ..., rm` will be converted to
#'   `x(i, 0), x(i, 1), ..., x(i, n-1)`.
#'
#' @details The function looks for patterns `r<number>` in the expression and
#' converts them to C++-style matrix indexing, where the number is decremented
#' by one (to match 0-based indexing in C++). Additionally, logical operators
#' like "AND" (or "and") and "OR" (or "or") are converted to their C++
#' equivalents (`&&` and `||`, respectively).
#'
#' @examples
#' \donttest{
#'   replace_r_indices("r1 && (r2 || r3)")
#'   # Output: "x(i, 0) * (std_min(double(1), x(i, 1) + x(i, 2)))"
#'
#'   replace_r_indices("(r1 AND r2) OR r3")
#'   # Output: "std_min(double(1), (x(i, 0) * x(i, 1)) + x(i, 2))"
#'
#'   replace_r_indices("(r1 AND r2) + r3 + r4")
#'   # Output: "(x(i, 0) \* x(i, 1)) + 5 \* x(i, 2) + x(i, 3)"
#' }
#'
#' @noRd
replace_r_indices <- function(expr_string) {
    # Replace logical ops (AND, OR in any case, &&, ||) with placeholders
    fixed_expr <- expr_string
    fixed_expr <- gsub("\\b[Aa][Nn][Dd]\\b", "%AND%", fixed_expr)
    fixed_expr <- gsub("\\b[Oo][Rr]\\b", "%OR%", fixed_expr)
    fixed_expr <- gsub(
        "||",
        "%OR%",
        gsub("&&", "%AND%", fixed_expr, fixed = TRUE),
        fixed = TRUE
    )

    # Parse modified string
    ast <- str2lang(fixed_expr)

    transformed <- parse_and_transform(ast)

    out_str <- deparse1(transformed$expr)
    out_str <- gsub("`", "", out_str, fixed = TRUE)

    # Insert space after/before slash
    out_str <- gsub("([[:alnum:]_\\)])/", "\\1 /", out_str)
    out_str <- gsub("/([[:alnum:]_\\(])", "/ \\1", out_str)
    out_str
}


count_unique_indices <- function(expr_string) {
    matches <- regmatches(expr_string, gregexpr("r[0-9]+", expr_string))
    all_matches <- unlist(matches)

    if (length(all_matches) == 0L) {
        cli::cli_abort(
            "{.arg objective} must reference at least one hypothesis \\
            indicator ({.code r1}, {.code r2}, ...). \\
            Expression {.code {expr_string}} contains none."
        )
    }

    indices <- gsub("r([0-9]+)", "\\1", all_matches)
    numeric_indices <- as.integer(indices)

    # Find unique indices and determine the maximum index based on
    # sequence completeness
    unique_indices <- unique(numeric_indices)
    max_index_found <- max(unique_indices, na.rm = TRUE)
    expected_indices <- 1:max_index_found

    # Check if all indices in the sequence are present
    missing_indices <- setdiff(expected_indices, unique_indices)

    # Check for any unexpected indices that are outside the maximum expected
    # range
    # nolint start: commented_code_linter
    # unexpected_indices <- unique_indices[unique_indices > max_index_found]
    # nolint end

    if (length(missing_indices) > 0) {
        cli::cli_warn(
            "Missing indices in the sequence. Expected every index from \\
            1 to {max_index_found}, but missing {toString(missing_indices)}."
        )
    }

    # Return the maximum index found if all checks pass
    max_index_found
}


# nolint start: return_linter

# Decide result type for arithmetic operations on left_type & right_type
combine_arithmetic <- function(op, left_type, right_type) {
    # We only handle +, -, *, / in a special manner
    if (!op %in% c("+", "-", "*", "/")) {
        return(list(type = "real", error_msg = NULL))
    }
    out <- list(type = "real", error_msg = NULL)

    # If both sides are bool, decide outcome
    if (left_type == "bool" && right_type == "bool") {
        if (op == "*") {
            # bool * bool => bool
            out$type <- "bool"
        } else {
            # +, -, / => real
            out$type <- "real"
        }
    }
    out
}

# Recursive function that:
# - transforms each node to valid C++ code
# - tracks its type: "bool" or "real"
parse_and_transform <- function(node) {
    if (is.call(node)) {
        parts <- as.list(node)
        fn <- parts[[1]]

        # --- Skip parentheses ---
        # If the function is '(' it's a parenthesized expression.
        # We transform whatever is inside and return that directly.
        if (identical(fn, quote(`(`))) {
            # e.g. node is something like ( r1 %OR% r2 )
            return(parse_and_transform(parts[[2]]))
        }

        # --- Handle placeholders for OR / AND ---
        if (identical(fn, quote(`%OR%`))) {
            # (A || B)
            left <- parse_and_transform(parts[[2]])
            right <- parse_and_transform(parts[[3]])
            if (left$type != "bool" || right$type != "bool") {
                stop("`||` (OR) only allowed between booleans.", call. = FALSE)
            }
            new_call <- bquote(std_min(double(1), .(left$expr) + .(right$expr)))
            return(list(expr = new_call, type = "bool"))
        } else if (identical(fn, quote(`%AND%`))) {
            # (A && B)
            left <- parse_and_transform(parts[[2]])
            right <- parse_and_transform(parts[[3]])
            if (left$type != "bool" || right$type != "bool") {
                stop("`&&` (AND) only allowed between booleans.", call. = FALSE)
            }
            new_call <- substitute(A * B, list(A = left$expr, B = right$expr))
            return(list(expr = new_call, type = "bool"))
        } else {
            # Possibly +, -, *, /, or some other function
            op_text <- as.character(fn)
            # Transform each sub-argument
            transformed_args <- lapply(parts[-1], parse_and_transform)
            # Rebuild the call
            new_call <- as.call(c(
                list(fn),
                lapply(transformed_args, `[[`, "expr")
            ))

            # If it's arithmetic, combine types
            if (op_text %in% c("+", "-", "*", "/")) {
                left_type <- transformed_args[[1]]$type
                right_type <- transformed_args[[2]]$type
                combo <- combine_arithmetic(op_text, left_type, right_type)
                if (!is.null(combo$error_msg)) {
                    stop(combo$error_msg, call. = FALSE)
                }
                return(list(expr = new_call, type = combo$type))
            } else {
                # Some other function => treat as real
                return(list(expr = new_call, type = "real"))
            }
        }
    } else if (is.symbol(node)) {
        # 1, r2, or other symbol
        txt <- as.character(node)
        if (grepl("^r\\d+$", txt)) {
            # r<number> => double(x(i, <number>-1)), treat as bool
            #
            # The double() wrap ensures every matrix access resolves to a `double`
            # at the C++ level, which keeps the generated expression type-
            # homogeneous. Required when the compiled function takes a
            # `LogicalMatrix` (element access returns int, which breaks mixed
            # std::min / arithmetic with double literals). A no-op for
            # NumericMatrix callers.
            idx <- as.integer(sub("^r", "", txt)) - 1
            return(list(
                expr = call("double", call("x", quote(i), idx)),
                type = "bool"
            ))
        } else {
            # Some other symbol => real
            return(list(expr = node, type = "real"))
        }
    } else if (is.numeric(node)) {
        # Numeric literal => real
        txt <- as.character(node)
        if (!grepl("\\.", txt)) {
            txt <- paste0(txt, ".0")
        }
        return(list(expr = as.symbol(txt), type = "real"))
    } else {
        # Strings, etc => real
        return(list(expr = node, type = "real"))
    }
}
# nolint end

is_trial_success <- function(x) {
    inherits(x, "multigrain_trial_success")
}

check_trial_success <- function(
    trial_success,
    arg = rlang::caller_arg(trial_success),
    call = rlang::caller_env(),
    allow_null = FALSE
) {
    if (!missing(trial_success)) {
        if (is_trial_success(trial_success)) {
            return(invisible(NULL))
        }

        if (allow_null && is.null(trial_success)) {
            return(invisible(NULL))
        }
    }

    rlang::stop_input_type(
        trial_success,
        "a multigrain trial success object",
        allow_null = allow_null,
        arg = arg,
        call = call
    )
}
