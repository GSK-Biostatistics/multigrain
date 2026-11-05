# Shared fixtures for the optimised HSD spending tests (build plan,
# section 4, "Shared test fixture"). Every fixture is built lazily on first
# use and cached, so test files that do not use them do not pay for them.

hsd_opt_cache <- new.env(parent = emptyenv())

hsd_opt_cached <- function(key, value) {
    if (!exists(key, envir = hsd_opt_cache, inherits = FALSE)) {
        assign(key, value, envir = hsd_opt_cache)
    }
    get(key, envir = hsd_opt_cache, inherits = FALSE)
}

# The fixed-spending stand-in for a tagged hypothesis at parameter `g`.
hsd_fixed_fn <- function(g) {
    force(g)
    function(a, t) gsDesign::sfHSD(a, t, g)$spend
}

# d = 1 fixture: m = 3, K = 2, one parameter "sec" shared by H2 and H3 on the
# lattice -4, -3.5, ..., -2 (5 nodes).
hsd_opt_raw <- function() {
    hsd_opt_cached(
        "raw",
        withr::with_seed(7, {
            simulate_pvalues_gsd(
                c(0.9, 0.85, 0.8),
                corr_matrix = diag(3),
                info_frac = c(0.5, 1),
                nsim = 2000
            )
        })
    )
}

hsd_opt_spending <- function() {
    list(
        gsDesign::sfLDOF,
        sf_hsd_opt(id = "sec", bounds = c(-4, -2), step = 0.5),
        sf_hsd_opt(id = "sec", bounds = c(-4, -2), step = 0.5)
    )
}

hsd_opt_spec <- function() {
    hsd_opt_cached(
        "spec",
        transform_pvalues_gsd(
            hsd_opt_raw(),
            # written out rather than `hsd_opt_spending()` so that the label
            # of H1 is `gsDesign::sfLDOF`
            spending = list(
                gsDesign::sfLDOF,
                sf_hsd_opt(id = "sec", bounds = c(-4, -2), step = 0.5),
                sf_hsd_opt(id = "sec", bounds = c(-4, -2), step = 0.5)
            ),
            grid_size = 64L
        )
    )
}

# d = 2 fixture: m = 4, K = 2. H2 has its own unnamed parameter, H3 the
# parameter "b"; H4 (parameter "c") has one analysis at full information and
# is short-circuited, so its group disappears. Both remaining lattices are
# -3, -2.5, -2.
hsd_opt_raw_d2 <- function() {
    hsd_opt_cached(
        "raw_d2",
        withr::with_seed(7, {
            simulate_pvalues_gsd(
                c(0.9, 0.85, 0.8, 0.8),
                corr_matrix = diag(4),
                info_frac = c(0.5, 1),
                nsim = 1000
            )
        })
    )
}

hsd_opt_info_frac_d2 <- function() {
    matrix(
        c(0.5, 1, 0.5, 1, 0.5, 1, NA, 1),
        nrow = 4L,
        byrow = TRUE
    )
}

hsd_opt_spending_d2 <- function() {
    list(
        gsDesign::sfLDOF,
        sf_hsd_opt(bounds = c(-3, -2), step = 0.5),
        sf_hsd_opt(id = "b", bounds = c(-3, -2), step = 0.5),
        sf_hsd_opt(id = "c", bounds = c(-3, -2), step = 0.5)
    )
}

hsd_opt_spec_d2 <- function() {
    hsd_opt_cached(
        "spec_d2",
        transform_pvalues_gsd(
            hsd_opt_raw_d2(),
            info_frac = hsd_opt_info_frac_d2(),
            spending = list(
                gsDesign::sfLDOF,
                sf_hsd_opt(bounds = c(-3, -2), step = 0.5),
                sf_hsd_opt(id = "b", bounds = c(-3, -2), step = 0.5),
                sf_hsd_opt(id = "c", bounds = c(-3, -2), step = 0.5)
            ),
            grid_size = 64L
        )
    )
}

# The tiny optimiser settings of test-optimisation_gsd.R, for later phases.
hsd_opt_control <- function() {
    multigrain_control() |>
        control_global(maxiter = 20, run = 5, popSize = 30) |>
        control_local(maxeval = 200) |>
        control_nsim_global(1000)
}
