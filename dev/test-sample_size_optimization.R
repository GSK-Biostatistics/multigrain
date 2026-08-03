conj_disj_ts <- trial_success((r1 && r3) || (r2 && r4))

test_that("split_theta_minN basic all-free 3x3 graph works", {
    .expected_encoded_lengths <- function(hyp_constraint, trans_constraint) {
        free_w <- sum(is.na(hyp_constraint))
        w_len <- max(free_w - 1L, 0L)
        row_free <- rowSums(is.na(trans_constraint))
        g_len <- sum(pmax(row_free - 1L, 0L))
        list(free_w = free_w, w_len = w_len, g_len = g_len)
    }

    cw <- c(NA_real_, NA_real_, NA_real_) # 3 free weights -> w_len = 2
    cG <- matrix(NA_real_, 3, 3) # each row 2 free -> 1 each -> g_len = 3
    diag(cG) <- 0

    info <- .expected_encoded_lengths(cw, cG)
    expect_identical(info$w_len, 2L, tolerance = sqrt(.Machine$double.eps))
    expect_identical(info$g_len, 3L, tolerance = sqrt(.Machine$double.eps))

    theta <- c(200, 0.3, 0.2, 0.6, 0.1, 1) # N + 2 + 3

    parts <- split_theta_minN(theta, cw, cG)
    expect_identical(parts$N, 200)
    expect_length(parts$w_pars, info$w_len)
    expect_length(parts$g_pars, info$g_len)

    # reconstruct and check simplex
    w_full <- recover_full_weights(parts$w_pars, cw)
    expect_identical(sum(w_full), 1, tolerance = sqrt(.Machine$double.eps))

    G_full <- recover_full_trans_matrix(parts$g_pars, cG)
    expect_identical(
        rowSums(G_full),
        rep(1, 3),
        tolerance = sqrt(.Machine$double.eps)
    )
})

test_that("split_theta_minN handles fixed + free weight mix", {
    set.seed(234)

    .expected_encoded_lengths <- function(hyp_constraint, trans_constraint) {
        free_w <- sum(is.na(hyp_constraint))
        w_len <- max(free_w - 1L, 0L)
        row_free <- rowSums(is.na(trans_constraint))
        g_len <- sum(pmax(row_free - 1L, 0L))
        list(free_w = free_w, w_len = w_len, g_len = g_len)
    }

    cw <- c(0.2, NA_real_, 0.2, NA_real_) # free_w = 2 -> w_len = 1
    cG <- matrix(NA_real_, 4, 4)
    cG[3, ] <- c(0, 1, 0, 0) # no free row
    diag(cG) <- 0 # all free rows (6 NA each) -> (3-1)*3 = 6 g params

    info <- .expected_encoded_lengths(cw, cG)
    expect_identical(info$w_len, 1L, tolerance = sqrt(.Machine$double.eps))
    expect_identical(info$g_len, 6L, tolerance = sqrt(.Machine$double.eps))

    rand_w <- runif(info$w_len)
    rand_w <- rand_w / sum(rand_w) / 5 # equal to 0.2
    rand_G <- runif(info$g_len)
    rand_G <- c(
        rand_G[1:2] / sum(rand_G[1:2]) / 1.25,
        rand_G[3:4] / sum(rand_G[3:4]) / 2,
        rand_G[5:6] / sum(rand_G[5:6]) / 3
    )
    theta <- c(300, rand_w, rand_G)
    parts <- split_theta_minN(theta, cw, cG)

    expect_identical(parts$N, 300)
    expect_length(parts$w_pars, info$w_len)
    expect_length(parts$g_pars, info$g_len)

    # round trip
    w_full <- recover_full_weights(parts$w_pars, cw)
    expect_identical(sum(w_full), 1, tolerance = 1e-12)
    expect_identical(w_full[!is.na(cw)], cw[!is.na(cw)])

    G_full <- recover_full_trans_matrix(parts$g_pars, cG)
    expect_identical(
        rowSums(G_full),
        rep(1, 4),
        tolerance = sqrt(.Machine$double.eps)
    )
    expect_identical(G_full[1, 4], 0.2, tolerance = 1e-12)
    expect_identical(G_full[2, 4], 0.5, tolerance = 1e-12)
    expect_identical(
        G_full[3, ],
        c(0, 1, 0, 0),
        tolerance = sqrt(.Machine$double.eps)
    )
    expect_identical(G_full[4, 3], 2 / 3, tolerance = sqrt(.Machine$double.eps))
})


test_that("split_theta_minN free_w == 1 encodes zero weight params", {
    set.seed(101)

    .expected_encoded_lengths <- function(hyp_constraint, trans_constraint) {
        free_w <- sum(is.na(hyp_constraint))
        w_len <- max(free_w - 1L, 0L)
        row_free <- rowSums(is.na(trans_constraint))
        g_len <- sum(pmax(row_free - 1L, 0L))
        list(free_w = free_w, w_len = w_len, g_len = g_len)
    }

    cw <- c(0.6, NA_real_, 0.2, 0.2) # free_w = 1 -> w_len = 0
    cG <- matrix(NA_real_, 4, 4)
    cG[3, ] <- c(0, 1, 0, 0) # no free row
    diag(cG) <- 0 # all free rows (6 NA each) -> (3-1)*3 = 6 g params

    info <- .expected_encoded_lengths(cw, cG)
    expect_identical(info$w_len, 0L, tolerance = sqrt(.Machine$double.eps))
    expect_identical(info$g_len, 6L, tolerance = sqrt(.Machine$double.eps))

    rand_G <- runif(info$g_len)
    rand_G <- c(
        rand_G[1:2] / sum(rand_G[1:2]) / 1.25,
        rand_G[3:4] / sum(rand_G[3:4]) / 2,
        rand_G[5:6] / sum(rand_G[5:6]) / 3
    )
    theta <- c(300, rand_G)
    parts <- split_theta_minN(theta, cw, cG)

    expect_identical(parts$N, 300, tolerance = sqrt(.Machine$double.eps))
    expect_length(parts$w_pars, 0L)
    expect_length(parts$g_pars, info$g_len)

    w_full <- recover_full_weights(parts$w_pars, cw)
    expect_identical(w_full, c(0.6, 0, 0.2, 0.2), tolerance = 1e-12)
    expect_identical(
        w_full[!is.na(cw)],
        cw[!is.na(cw)],
        tolerance = sqrt(.Machine$double.eps)
    )
})


test_that("split_theta_minN encodes zero transition weigts when G is fixed", {
    set.seed(234)

    .expected_encoded_lengths <- function(hyp_constraint, trans_constraint) {
        free_w <- sum(is.na(hyp_constraint))
        w_len <- max(free_w - 1L, 0L)
        row_free <- rowSums(is.na(trans_constraint))
        g_len <- sum(pmax(row_free - 1L, 0L))
        list(free_w = free_w, w_len = w_len, g_len = g_len)
    }

    cw <- c(0.2, NA_real_, 0.2, NA_real_, NA_real_) # free_w = 3 -> w_len = 2
    cG <- random_transitions(5)

    info <- .expected_encoded_lengths(cw, cG)
    expect_identical(info$w_len, 2L, tolerance = sqrt(.Machine$double.eps))
    expect_identical(info$g_len, 0L, tolerance = sqrt(.Machine$double.eps))

    rand_w <- runif(info$w_len)
    rand_w <- rand_w / sum(rand_w) / 5 # equal to 0.2
    theta <- c(300, rand_w)
    parts <- split_theta_minN(theta, cw, cG)

    expect_identical(parts$N, 300)
    expect_length(parts$w_pars, 2L)
    expect_length(parts$g_pars, 0L)

    # round trip
    w_full <- recover_full_weights(parts$w_pars, cw)
    expect_identical(sum(w_full), 1, tolerance = 1e-12)
    expect_identical(w_full[!is.na(cw)], cw[!is.na(cw)])

    G_full <- recover_full_trans_matrix(parts$g_pars, cG)
    expect_identical(G_full, cG, tolerance = 1e-12)
})


test_that("split_theta_minN respects row-wise G constraints", {
    # 4x4 example with row heterogeneity

    .expected_encoded_lengths <- function(hyp_constraint, trans_constraint) {
        free_w <- sum(is.na(hyp_constraint))
        w_len <- max(free_w - 1L, 0L)
        row_free <- rowSums(is.na(trans_constraint))
        g_len <- sum(pmax(row_free - 1L, 0L))
        list(free_w = free_w, w_len = w_len, g_len = g_len)
    }

    cG <- matrix(c(
    0, NA, NA, NA,
    0.5, 0, 0.5, 0,
    NA, NA,  0, 0.2,
    NA, 0, NA, 0
  ), nrow = 4, byrow = TRUE)

    cw <- c(NA, NA, NA, NA) # 4 free weights -> w_len = 3

    info <- .expected_encoded_lengths(cw, cG)
    expect_identical(
        info$w_len,
        3L,
        tolerance = sqrt(.Machine$double.eps)
    ) # (4-1)
    expect_identical(
        info$g_len,
        4L,
        tolerance = sqrt(.Machine$double.eps)
    )

    theta <- c(250, c(0.1, 0.2, 0.5), c(0.1, 0.2, 0.5, 0.7))
    parts <- split_theta_minN(theta, cw, cG)

    expect_identical(parts$N, 250)
    expect_length(parts$w_pars, info$w_len)
    expect_length(parts$g_pars, info$g_len)

    # round trip
    w_full <- recover_full_weights(parts$w_pars, cw)
    expect_identical(sum(w_full), 1, tolerance = 1e-12)

    G_full <- recover_full_trans_matrix(parts$g_pars, cG)
    expect_identical(rowSums(G_full), rep(1, 4), tolerance = 1e-12)
    # row2 fixed row -> unchanged
    expect_identical(G_full[2, ], c(0.5, 0, 0.5, 0), tolerance = 1e-12)
})


test_that("split_theta_minN free_w == 1 encodes zero weight params", {
    .expected_encoded_lengths <- function(hyp_constraint, trans_constraint) {
        free_w <- sum(is.na(hyp_constraint))
        w_len <- max(free_w - 1L, 0L)
        row_free <- rowSums(is.na(trans_constraint))
        g_len <- sum(pmax(row_free - 1L, 0L))
        list(free_w = free_w, w_len = w_len, g_len = g_len)
    }

    cw <- c(0.6, NA_real_, 0.2, 0.2) # free_w = 1 -> w_len = 0
    cG <- matrix(NA_real_, 4, 4) # all free -> g_len = (4-1)*4 = 12

    info <- .expected_encoded_lengths(cw, cG)
    expect_identical(
        info$w_len,
        0L,
        tolerance = sqrt(.Machine$double.eps)
    )
    expect_identical(
        info$g_len,
        12L,
        tolerance = sqrt(.Machine$double.eps)
    )

    theta <- c(150, runif(info$g_len)) # N + 0 + 12
    parts <- split_theta_minN(theta, cw, cG)

    expect_identical(parts$N, 150)
    expect_length(parts$w_pars, 0L)
    expect_length(parts$g_pars, info$g_len)

    w_full <- recover_full_weights(parts$w_pars, cw)
    expect_identical(sum(w_full), 1, tolerance = 1e-12)
    expect_identical(w_full[!is.na(cw)], cw[!is.na(cw)])
})


test_that(".optimise_N_ga returns a valid graph_optimal object with N in bounds", {
    set.seed(1001)

    alpha <- 0.025
    m <- 4
    power_vector <- c(0.95, 0.86, 0.81, 0.79)
    ncp0_vec <- rep(qnorm(1 - alpha), 4) - qnorm(1 - power_vector)
    mu_treat <- ncp0_vec * sqrt(1 / 600)
    Sigma <- diag(m)
    N_max <- 500

    gc <- graph_constraint(
        hyp_constraint = c(1, 0, 0, 0),
        trans_constraint = rbind(
            c(0, NA, NA, 0),
            c(0, 0, NA, NA),
            c(0, NA, 0, NA),
            c(0, 0, 1, 0)
        )
    )

    # small GA for speed
    ga_opts <- list(
        popSize = 150,
        maxiter = 50,
        max_stagnation = 10,
        pmutation = 0.9,
        pcrossover = 0.1
    )

    res <- .optimise_N_ga(
        ncp_vector = mu_treat,
        sigma = Sigma,
        graph_constraint = gc,
        trial_success = conj_disj_ts, # see top of script
        N_max = N_max,
        alpha = alpha,
        trial_success_target = 0.6,
        nsim_global = 1e4,
        ga_options = ga_opts,
        start_graph = list(
            list(
                hyp_weight = NULL,
                trans_matrix = NULL
            )
        ),
        verbose = TRUE
    )

    expect_s3_class(res, "graph_optimal")
    expect_type(res$N, "double")
    expect_length(res$N, 1L)
    expect_gte(res$N, 2)
    expect_lte(res$N, N_max)

    # validity of returned graph
    expect_true(is_graph_valid(res$hyp_weight, res$trans_matrix))
    expect_equal(sum(res$hyp_weight), 1, tolerance = 1e-10)
    expect_equal(rowSums(res$trans_matrix), rep(1, m), tolerance = 1e-10)

    # basic structure
    expect_true(is.list(res$power) || is.null(res$power)) # nolint
    expect_type(res$global_output, "list")
    expect_type(res$opt_settings, "list")
    expect_type(res$solution, "list")
    expect_true(all("global" %in% names(res$solution$graph_valid)))
})


test_that("optimise_N returns valid structure (GA + local) with N in bounds", {
    set.seed(1001)

    alpha <- 0.025
    m <- 4
    power_vector <- c(0.95, 0.86, 0.81, 0.79)
    ncp0_vec <- rep(qnorm(1 - alpha), 4) - qnorm(1 - power_vector)
    mu_treat <- ncp0_vec * sqrt(1 / 600)
    Sigma <- diag(m)
    N_max <- 2000

    gc <- graph_constraint(
        hyp_constraint = c(1, 0, 0, 0),
        trans_constraint = rbind(
            c(0, NA, NA, 0),
            c(0, 0, NA, NA),
            c(0, NA, 0, NA),
            c(0, 0, 1, 0)
        )
    )

    res <- optimise_N(
        ncp_vector = mu_treat,
        sigma = Sigma,
        graph_constraint = gc,
        trial_success = conj_disj_ts, # see top of script
        N_max = N_max,
        alpha = alpha,
        trial_success_target = 0.35,
        nsim_global = 2000,
        nsim_local = 2000,
        ga_options = ga_N_options(
            pop_size = 200,
            max_iter = 300,
            max_stagnation = 10,
            p_mutation = 0.9,
            p_crossover = 0.1
        ),
        start_graph = list(
            list(
                hyp_weight = NULL,
                trans_matrix = NULL
            )
        ),
        verbose = FALSE
    )

    expect_s3_class(res, "graph_optimal")
    expect_type(res$N, "double")
    expect_gte(res$N, 2)
    expect_lte(res$N, N_max)

    # local + global power summaries present
    expect_type(res$global_opt_power, "list")
    expect_type(res$local_opt_power, "list")
    expect_type(res$power, "list")

    # validity
    expect_true(is_graph_valid(res$hyp_weight, res$trans_matrix))
    expect_identical(sum(res$hyp_weight), 1, tolerance = 1e-10)
    expect_identical(
        as.vector(rowSums(res$trans_matrix)),
        rep(1, m),
        tolerance = 1e-10
    )

    # names propagated
    expect_named(res$hyp_weight, names(gc$hyp_constraint))
    expect_identical(rownames(res$trans_matrix), names(gc$hyp_constraint))
    expect_identical(colnames(res$trans_matrix), names(gc$hyp_constraint))
})


test_that("Optimise_N() for 2m: E2E check", {
    set.seed(1)
    ts <- trial_success(r1 || r2)

    ncp <- c(0.12, 0.1)
    corr_matrix <- matrix(c(1, 0.2, 0.2, 1), nrow = 2)

    gc <- graph_constraint(c(NA, NA), matrix(c(0, 1, 1, 0), nrow = 2))

    set.seed(123)
    graph_2m <- optimise_N(
        ncp,
        sigma = corr_matrix,
        graph_constraint = gc,
        trial_success = ts,
        trial_success_target = 0.5,
        N_max = 400,
        alpha = 0.025,
        ga_options = ga_N_options(
            max_stagnation = 20,
            max_iter = 100
        ),
        nsim_global = 1e3
    )

    expect_length(graph_2m$hyp_weight, 2L)
    expect_s3_class(graph_2m, "graph_optimal")
    expect_true(graph_2m$global_search)
    expect_identical(unname(graph_2m$trans_matrix), unname(gc$trans_constraint))
    expect_lte(graph_2m$N, 400)
})


# test_that("optimise_N returns expected warning when N_star > N_max", {
#   set.seed(1001)
#
#   alpha <- 0.025
#   m <- 4
#   power_vector <- c(0.95, 0.86, 0.81, 0.79)
#   ncp0_vec     <- rep(qnorm(1 - alpha), 4) - qnorm(1 - power_vector)
#   mu_treat     <- ncp0_vec * sqrt(1 / 600)
#   Sigma <- diag(m)
#   N_max <- 50
#
#   gc <- graph_constraint(
#     hyp_constraint = c(1, 0, 0, 0),
#     trans_constraint = rbind(
#       c(0, NA, NA, 0),
#       c(0, 0,  NA, NA),
#       c(0, NA, 0,  NA),
#       c(0, 0,  1,  0)
#     )
#   )
#   expect_warning(
#     res <- optimise_N( # nolint: implicit_assignment_linter.
#       ncp_vector           = mu_treat,
#       sigma                = Sigma,
#       graph_constraint     = gc,
#       trial_success        = conj_disj_ts, # see top of script
#       N_max                = N_max,
#       alpha                = alpha,
#       trial_success_target = 0.8,
#       nsim_global          = 2000,
#       nsim_local           = 2000,
#       ga_options           = list(popSize = 200, maxiter = 300, max_stagnation = 10,
#                                   pmutation = 0.9, pcrossover = 0.1),
#       start_graph          = list(list(hyp_weight = NULL, trans_matrix = NULL)),
#       verbose              = FALSE
#     ),
#     regexp = paste("Could not find graph with sample size less than N =", N_max)
#   )
#
#   expect_s3_class(res, "graph_optimal")
# })

## param_to_solution tests

# 2 hypothesis test
test_that("Correct param_to_solution_minN() behaviour for 2 hypotheses", {
    gc <- graph_constraint(c(NA, NA), matrix(c(0, 1, 1, 0), nrow = 2))
    optim_param <- c(200, 0.25)
    sol <- param_to_solution_minN(optim_param, gc)
    expect_identical(sol$N, 200)
    expect_identical(sol$hyp_weight, c(0.25, 0.75))
    expect_identical(sol$trans_matrix, unname(gc$trans_constraint))
})
