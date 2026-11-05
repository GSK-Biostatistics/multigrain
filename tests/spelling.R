# nolint start: missing_package_linter, namespace_linter.
if (requireNamespace("spelling", quietly = TRUE)) {
    spelling::spell_check_test(
        vignettes = TRUE,
        error = FALSE,
        skip_on_cran = TRUE
    )
}
# nolint end
