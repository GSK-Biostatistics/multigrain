# `create_obj_func_gsd()` must be the fixed-sample factory with one
# substitution: the group sequential kernel, and the decision-time matrix
# handed to the gain (design record, section 6 P4).

obj_pvals_gsd <- withr::with_seed(31, {
    raw <- simulate_pvalues_gsd(
        power_nominal = c(0.9, 0.85, 0.8),
        corr_matrix = diag(3),
        info_frac = c(0.5, 1),
        nsim = 1000
    )
    transform_pvalues_gsd(
        raw,
        spending = gsDesign::sfLDOF,
        grid_size = 128L
    )
})

obj_gain_gsd <- trial_success_gsd(
    d(t1) + d(t2) + d(t3),
    d = c(1, 0.75),
    verbose = "silent"
)

obj_free_3 <- graph_constraint_free(3)

# Decode an encoded parameter vector exactly as the closure does, including
# the two snaps, so the expected value can be built by hand.
obj_decode <- function(x, hyp_constraint, trans_constraint) {
    theta <- split_theta(x, hyp_constraint)
    w <- recover_full_weights(theta$w_pars, hyp_constraint)
    G <- recover_full_trans_matrix(theta$g_pars, trans_constraint)
    w[w < 1e-4] <- 0
    G[G < 1e-5] <- 0
    list(w = w, G = G)
}

obj_make_gsd <- function(pvals, num_threads = 1L, alpha = NULL) {
    create_obj_func_gsd(
        m = obj_gain_gsd$m,
        power_criterion = obj_gain_gsd$func,
        hyp_constraint = obj_free_3$hyp_constraint,
        trans_constraint = obj_free_3$trans_constraint,
        pvals = .gsd_kernel_matrix(pvals),
        K = pvals$K,
        alpha = alpha %||% pvals$alpha,
        num_threads = num_threads
    )
}

obj_x <- c(0.45, 0.3, 0.6, 0.2, 0.7)


test_that("the closure equals the kernel and gain called by hand", {
    fn <- obj_make_gsd(obj_pvals_gsd)

    decoded <- obj_decode(
        obj_x,
        obj_free_3$hyp_constraint,
        obj_free_3$trans_constraint
    )

    res <- graph_shortcut_gsd(
        pvals = .gsd_kernel_matrix(obj_pvals_gsd),
        alpha = obj_pvals_gsd$alpha,
        w = decoded$w,
        G = decoded$G,
        K = obj_pvals_gsd$K
    )

    expect_identical(fn(obj_x), obj_gain_gsd$func(res$time))
})


test_that("the same encoded vector gives an identical value twice", {
    # design record, section 8 item 11 and brief claim 5
    fn <- obj_make_gsd(obj_pvals_gsd)
    expect_identical(fn(obj_x), fn(obj_x))
})


test_that("the parallel closure equals the serial closure", {
    serial <- obj_make_gsd(obj_pvals_gsd, num_threads = 1L)
    parallel <- obj_make_gsd(obj_pvals_gsd, num_threads = 2L)

    expect_identical(parallel(obj_x), serial(obj_x))
})


test_that(".gsd_kernel_matrix() aborts on a planted NA", {
    broken <- obj_pvals_gsd
    broken$pvals[1L, 1L, 1L] <- NA_real_

    expect_error(
        .gsd_kernel_matrix(broken),
        "contains"
    )
    expect_error(
        .gsd_kernel_matrix(broken),
        "never rejectable"
    )
})


test_that(".gsd_kernel_matrix() aborts on a value outside [0, 1]", {
    # a negative repeated p-value is silently a rejection at any positive
    # allocation; one above 1 is silently "never"
    negative <- obj_pvals_gsd
    negative$pvals[1L, 1L, 1L] <- -1

    above <- obj_pvals_gsd
    above$pvals[3L, 2L, 2L] <- 2

    expect_error(
        .gsd_kernel_matrix(negative),
        "outside [0, 1]",
        fixed = TRUE
    )
    expect_error(
        .gsd_kernel_matrix(above),
        "outside [0, 1]",
        fixed = TRUE
    )

    # the boundaries themselves are fine
    edges <- obj_pvals_gsd
    edges$pvals[1L, 1L, 1L] <- 0
    edges$pvals[2L, 1L, 1L] <- 1
    expect_no_error(.gsd_kernel_matrix(edges))
})


test_that("every penalty branch matches create_obj_func()", {
    # The penalty branches sit before the kernel call and are copied verbatim
    # from `create_obj_func()`, so both closures must return the same number
    # for the same encoded vector, whatever the constraints.
    fixed_matrix <- matrix(0.5, nrow = 10L, ncol = 3L)

    free_hyp <- rep(NA_real_, 3L)
    free_trans <- matrix(NA_real_, nrow = 3L, ncol = 3L)
    diag(free_trans) <- 0

    part_trans <- rbind(
        c(0, NA, NA),
        c(0.5, 0, 0.5),
        c(0.5, 0.5, 0)
    )

    fixed_trans_big <- rbind(
        c(0, 1.5, 0),
        c(0.5, 0, 0.5),
        c(0.5, 0.5, 0)
    )

    cases <- list(
        na_weight = list(
            hyp = free_hyp,
            trans = free_trans,
            x = c(NA_real_, 0.3, 0.5, 0.5, 0.5)
        ),
        negative_weight = list(
            hyp = free_hyp,
            trans = free_trans,
            x = c(1.4, 0.2, 0.5, 0.5, 0.5)
        ),
        weight_above_one = list(
            hyp = c(1.5, 0.3, 0.2),
            trans = part_trans,
            x = 0.4
        ),
        negative_edge = list(
            hyp = free_hyp,
            trans = part_trans,
            x = c(0.4, 0.3, 1.4)
        ),
        edge_above_one = list(
            hyp = free_hyp,
            trans = fixed_trans_big,
            x = c(0.4, 0.3)
        )
    )

    for (nm in names(cases)) {
        case <- cases[[nm]]

        gsd_fn <- create_obj_func_gsd(
            m = obj_gain_gsd$m,
            power_criterion = obj_gain_gsd$func,
            hyp_constraint = case$hyp,
            trans_constraint = case$trans,
            pvals = .gsd_kernel_matrix(obj_pvals_gsd),
            K = obj_pvals_gsd$K,
            alpha = obj_pvals_gsd$alpha
        )
        fixed_fn <- create_obj_func(
            m = 3L,
            power_criterion = function(x) 0,
            hyp_constraint = case$hyp,
            trans_constraint = case$trans,
            pvals = fixed_matrix,
            alpha = 0.025
        )

        expect_identical(
            gsd_fn(case$x),
            fixed_fn(case$x),
            label = nm
        )
        # a penalty, not a gain: the kernel was never reached
        expect_lte(gsd_fn(case$x), 0)
    }
})
