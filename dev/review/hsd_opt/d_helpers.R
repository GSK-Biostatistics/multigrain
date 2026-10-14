# d_helpers.R
# Shared fixtures and prototype helpers for the Phase D batteries d02, d04 and
# d05 (dev/review/hsd_opt/briefs/D-b_d02_d04_d05.md). Sourced by those scripts;
# not part of the package. No package code is changed: the installed library
# is used through `library(multigrain)` and `multigrain:::`.

suppressPackageStartupMessages({
    library(multigrain)
    library(gsDesign)
})

here <- "dev/review/hsd_opt"
if (!dir.exists(here)) {
    here <- "."
}

n_threads <- 2L
alpha <- 0.025
lattice <- seq(-4, 1, by = 0.25)
outer_nodes <- c(-4, -3, -2, -1, 0, 1)
plateau_tol <- 5e-4

hsd_fn <- function(g) {
    force(g)
    function(a, t) gsDesign::sfHSD(a, t, g)$spend
}

fmt_node <- function(g) formatC(g, format = "f", digits = 3)

# ---- designs -----------------------------------------------------------------

# M0 / M1 (README 7.4) at nsim = 1e5 and the README section-2 design of the
# brief. `spending(g)` returns the spending list with the optimised axis at g.
make_design <- function(name) {
    if (name %in% c("M0", "M1")) {
        set.seed(20261007)
        raw <- simulate_pvalues_gsd(
            c(H1 = 0.90, H2 = 0.80),
            alpha = alpha,
            corr_matrix = matrix(c(1, 0.5, 0.5, 1), 2),
            info_frac = c(0.5, 1),
            nsim = 1e5
        )
        gain <- trial_success_gsd(
            0.6 * d(t1) + 0.4 * d(t2),
            d = c(1, 0.8),
            verbose = "silent"
        )
        rule <- if (name == "M1") stopping_rule_gsd(r1) else NULL
        spending <- function(g) list(gsDesign::sfLDOF, hsd_fn(g))
        list(
            name = name,
            raw = raw,
            gc = graph_constraint_free(2),
            gain = gain,
            rule = rule,
            spending = spending,
            m = 2L,
            opt_cols_hyp = 2L
        )
    } else if (name == "S2") {
        set.seed(1)
        raw <- simulate_pvalues_gsd(
            c(0.9, 0.85, 0.8),
            corr_matrix = diag(3),
            info_frac = c(0.5, 1),
            nsim = 1e5
        )
        gain <- trial_success_gsd(r1 + r2 + r3, verbose = "silent")
        spending <- function(g) list(gsDesign::sfLDOF, hsd_fn(g), hsd_fn(g))
        list(
            name = name,
            raw = raw,
            gc = graph_constraint_free(3),
            gain = gain,
            rule = NULL,
            spending = spending,
            m = 3L,
            opt_cols_hyp = c(2L, 3L)
        )
    } else {
        stop("unknown design ", name)
    }
}

stop_table <- function(rule, pv) {
    if (is.null(rule)) {
        return(NULL)
    }
    multigrain:::.gsd_stop_table(
        multigrain:::.gsd_check_stopping_rule_dims(rule, pv),
        pv$m,
        pv$K
    )
}

transform_node <- function(des, g) {
    t0 <- proc.time()[["elapsed"]]
    pv <- transform_pvalues_gsd(des$raw, spending = des$spending(g), alpha = alpha)
    list(pv = pv, secs = proc.time()[["elapsed"]] - t0)
}

# ---- counting gain -----------------------------------------------------------

# A copy of the gain whose `$func` counts its calls; still a valid
# `multigrain_trial_success_gsd` (only `$func` is replaced).
counting_gain <- function(gain) {
    counter <- new.env()
    counter$n <- 0
    f <- gain$func
    out <- gain
    out$func <- function(t) {
        counter$n <- counter$n + 1
        f(t)
    }
    attr(out, "counter") <- counter
    out
}
count_of <- function(gain_c) attr(gain_c, "counter")$n
reset_count <- function(gain_c) attr(gain_c, "counter")$n <- 0

# ---- evaluation on a node's array --------------------------------------------

# Gain of a decoded graph on the full kernel matrix (no thresholding beyond
# what the kernel does; this is what calc_power_pvals_gsd() scores).
eval_graph <- function(kmat, K, gain, w, G, stop_rule) {
    res <- multigrain:::graph_shortcut_gsd_parallel(
        pvals = kmat,
        alpha = alpha,
        w = w,
        G = G,
        K = K,
        num_threads = n_threads,
        grain_size = -1L,
        stop_rule = stop_rule
    )
    gain$func(res$time)
}

encode_graph <- function(gc, w, G) {
    proj <- multigrain:::closest_graph_to_constraints(gc, w, G)
    as.numeric(multigrain:::create_start_params(
        gc,
        w0 = proj$hyp_weight,
        G0 = proj$trans_matrix,
        sum_to_one_constraint = FALSE
    ))
}

local_opts_default <- function(nsim, m) {
    multigrain:::control_prepare_dims(
        multigrain_control(),
        nsim = nsim,
        m = m,
        verbose = "silent"
    )$local_opt
}

# Warm-started COBYLA refinement of the graph on one node's full array, from
# the better of the candidate graphs. Returns gains, the graph and the number
# of objective evaluations (counted in a wrapper around the closure).
refine_node <- function(des, pv, candidates, stop_rule, local_opts) {
    kmat <- multigrain:::.gsd_kernel_matrix(pv)
    cand_vals <- vapply(
        candidates,
        function(cg) eval_graph(kmat, pv$K, des$gain, cg$w, cg$G, stop_rule),
        numeric(1)
    )
    best <- which.max(cand_vals)
    x0 <- encode_graph(des$gc, candidates[[best]]$w, candidates[[best]]$G)
    obj <- multigrain:::create_obj_func_gsd(
        m = des$gain$m,
        power_criterion = des$gain$func,
        hyp_constraint = des$gc$hyp_constraint,
        trans_constraint = des$gc$trans_constraint,
        alpha = alpha,
        pvals = kmat,
        K = pv$K,
        num_threads = n_threads,
        stop_rule = stop_rule
    )
    n_eval <- 0
    nlopt_fn <- function(x) {
        n_eval <<- n_eval + 1
        -obj(x)
    }
    t0 <- proc.time()[["elapsed"]]
    nl <- nloptr::nloptr(
        x0 = x0,
        eval_f = nlopt_fn,
        lb = rep(0, length(x0)),
        ub = rep(1, length(x0)),
        opts = local_opts
    )
    secs <- proc.time()[["elapsed"]] - t0
    sol <- multigrain:::param_to_solution(nl$solution, des$gc, process = TRUE)
    sol <- multigrain:::repair_graph(sol$hyp_weight, sol$trans_matrix, des$gc)
    refined <- eval_graph(kmat, pv$K, des$gain, sol$hyp_weight, sol$trans_matrix, stop_rule)
    all_vals <- c(cand_vals, refined)
    graphs <- c(candidates, list(list(w = sol$hyp_weight, G = sol$trans_matrix)))
    pick <- which.max(all_vals)
    list(
        cand_vals = cand_vals,
        start_from = best,
        refined = refined,
        gain = all_vals[[pick]],
        graph = graphs[[pick]],
        refined_graph = graphs[[length(graphs)]],
        n_eval = n_eval,
        cobyla_secs = secs,
        nlopt_status = nl$status,
        nlopt_iter = nl$iterations
    )
}

# Connected run of nodes within tol of the maximum, containing the argmax.
plateau_of <- function(nodes, vals, tol = plateau_tol) {
    i <- which.max(vals)
    ok <- vals >= max(vals) - tol
    lo <- i
    while (lo > 1 && ok[lo - 1]) lo <- lo - 1
    hi <- i
    while (hi < length(vals) && ok[hi + 1]) hi <- hi + 1
    c(nodes[lo], nodes[hi])
}

print_df <- function(df, digits = 6) {
    print(format(df, digits = digits), row.names = FALSE)
}
