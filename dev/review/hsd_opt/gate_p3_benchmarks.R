# gate_p3_benchmarks.R
# Dev benchmark of phase P3 (dev/hsd_opt_build_plan.md, section 4, "P3 ...
# Dev benchmark"; design section 5.3 P3 (g)). Run by the orchestrator from
# the repository root, not part of the test suite:
#
#   Rscript dev/review/hsd_opt/gate_p3_benchmarks.R
#
# Environment variables:
#   GATE_P3_SMOKE=1      nsim 5e3 and tiny optimiser settings, to prove the
#                        script runs (inequalities then meaningless: the
#                        draws differ from d02's); output goes to tempdir()
#   GATE_P3_PARTS=M1,S2,D2   which parts to run (default all three)
#   GATE_P3_FULL_GRID=1  d = 2 reference over all 25 pairs instead of the 9
#                        pairs nearest the start node
#   GATE_P3_THREADS=2    num_threads for every optimiser call
#
# Checks, each with tolerance 5e-4 on the trial success measure:
#   (1) d = 1 (M1 and the README section-2 design S2, nsim 1e5, default
#       bounds -4 to 1, step 0.25, 21 nodes): joint gain >= the best
#       outer-loop gain of d02 (full default graph_optimise_gsd() at fixed
#       gamma; also shown: the best d02 warm-started profile value) - 5e-4.
#   (2) d = 2 coarse (M1 with both hypotheses tagged in separate groups,
#       bounds -2 to 0, step 0.5: 5 x 5 nodes, nsim 2e4): joint gain >= the
#       product-grid maximum of default graph_optimise_gsd() runs at fixed
#       pairs - 5e-4.
#   (3) every part: joint gain >= the fixed-spending gain (default
#       graph_optimise_gsd() on resolve_pvalues_gsd(spec, start)) - 5e-4.
#   Wall time of the S2 run (m = 3, K = 2, nsim 1e5, d = 1) against the
#   design section 4.2 projection of 3.5 to 6.5 min (pause if > 2x).

Sys.setenv(NOT_CRAN = "true")
suppressMessages(devtools::load_all(quiet = TRUE))

here <- "dev/review/hsd_opt"
smoke <- identical(Sys.getenv("GATE_P3_SMOKE"), "1")
full_grid <- identical(Sys.getenv("GATE_P3_FULL_GRID"), "1")
parts <- strsplit(Sys.getenv("GATE_P3_PARTS", "M1,S2,D2"), ",")[[1]]
n_threads <- as.integer(Sys.getenv("GATE_P3_THREADS", "2"))
tol <- 5e-4
alpha <- 0.025

nsim_d1 <- if (smoke) 5e3 else 1e5
nsim_d2 <- if (smoke) 5e3 else 2e4
out_dir <- if (smoke) tempdir() else here
out_stem <- file.path(
    out_dir,
    if (smoke) "gate_p3_benchmarks_smoke" else "gate_p3_benchmarks"
)

control <- if (smoke) {
    multigrain_control() |>
        control_global(maxiter = 5, run = 3, popSize = 30) |>
        control_local(maxeval = 200) |>
        control_nsim_global(1000)
} else {
    multigrain_control()
}

options(multigrain_verbosity = "silent")

if (!exists("%||%", mode = "function")) {
    `%||%` <- function(x, y) if (is.null(x)) y else x
}

log_file <- paste0(out_stem, "_output.txt")
log_con <- file(log_file, open = "wt")
say <- function(...) {
    line <- sprintf(...)
    cat(line, "\n", sep = "")
    cat(line, "\n", sep = "", file = log_con)
}
elapsed <- function(expr) {
    t0 <- proc.time()[["elapsed"]]
    value <- force(expr)
    list(value = value, secs = proc.time()[["elapsed"]] - t0)
}

say(
    "gate_p3_benchmarks: %s; R %s; multigrain %s; num_threads %d",
    if (smoke) "SMOKE (nsim 5e3, tiny control)" else "full",
    getRversion(),
    utils::packageVersion("multigrain"),
    n_threads
)
say("parts: %s; d = 2 reference: %s", toString(parts),
    if (full_grid) "all 25 pairs" else "9 pairs nearest the start node")

d02 <- readRDS(file.path(here, "d02_warm_start_profile_results.rds"))


# ---- designs (copied from d_helpers.R so the draws are d02's) ---------------

simulate_m1 <- function(nsim) {
    set.seed(20261007)
    simulate_pvalues_gsd(
        c(H1 = 0.90, H2 = 0.80),
        alpha = alpha,
        corr_matrix = matrix(c(1, 0.5, 0.5, 1), 2),
        info_frac = c(0.5, 1),
        nsim = nsim
    )
}

gain_m1 <- trial_success_gsd(
    0.6 * d(t1) + 0.4 * d(t2),
    d = c(1, 0.8),
    verbose = "silent"
)
rule_m1 <- stopping_rule_gsd(r1)

simulate_s2 <- function(nsim) {
    set.seed(1)
    simulate_pvalues_gsd(
        c(0.9, 0.85, 0.8),
        corr_matrix = diag(3),
        info_frac = c(0.5, 1),
        nsim = nsim
    )
}

gain_s2 <- trial_success_gsd(r1 + r2 + r3, verbose = "silent")


# ---- one joint run plus the fixed-spending reference at the start node ------

run_joint <- function(name, spec, gc, gain, rule, seed = 1) {
    joint <- elapsed(withr::with_seed(seed, {
        graph_optimise_gsd(
            spec,
            graph_constraint = gc,
            trial_success = gain,
            stopping_rule = rule,
            num_threads = n_threads,
            control = control
        )
    }))
    res <- joint$value
    say(
        "[%s] joint: gain %.6f at %s (source %s, converged %s); %.1f s",
        name,
        res$power$trial_success,
        toString(sprintf("%s = %g", names(res$spending_param),
                         res$spending_param)),
        res$solution$opt_source,
        res$profile_converged,
        joint$secs
    )
    say("[%s] profile:", name)
    prof <- res$spending_profile
    utils::capture.output(print(prof, digits = 6, row.names = FALSE)) |>
        vapply(function(l) {
            say("  %s", l)
            ""
        }, character(1)) |>
        invisible()
    say(
        "[%s] GA visited %d of %d nodes",
        name,
        sum(res$spending_search$visits > 0, na.rm = TRUE),
        nrow(res$spending_search)
    )

    fixed <- elapsed(withr::with_seed(seed, {
        graph_optimise_gsd(
            resolve_pvalues_gsd(spec, spec$param$start),
            graph_constraint = gc,
            trial_success = gain,
            stopping_rule = rule,
            num_threads = n_threads,
            control = control
        )
    }))
    fixed_gain <- fixed$value$power$trial_success
    check3 <- res$power$trial_success >= fixed_gain - tol
    say(
        "[%s] (3) fixed spending at the start node %s: gain %.6f (%.1f s); joint - fixed = %+.2e -> %s",
        name,
        toString(spec$param$start),
        fixed_gain,
        fixed$secs,
        res$power$trial_success - fixed_gain,
        if (check3) "PASS" else "FAIL"
    )

    list(
        result = res,
        joint_secs = joint$secs,
        fixed_gain = fixed_gain,
        fixed_secs = fixed$secs,
        check3 = check3
    )
}


# ---- d = 1 parts: M1 and S2 against d02 --------------------------------------

run_d1 <- function(name) {
    raw <- if (name == "M1") simulate_m1(nsim_d1) else simulate_s2(nsim_d1)
    spending <- if (name == "M1") {
        list(gsDesign::sfLDOF, sf_hsd_opt(id = "gamma"))
    } else {
        list(gsDesign::sfLDOF, sf_hsd_opt(id = "gamma"), sf_hsd_opt(id = "gamma"))
    }
    built <- elapsed(transform_pvalues_gsd(raw, spending = spending))
    spec <- built$value
    say(
        "[%s] spec: %d nodes, %d tables, start %g; built in %.1f s",
        name,
        length(spec$param$lattice[[1]]),
        length(spec$node_tables),
        spec$param$start[[1]],
        built$secs
    )

    gc <- graph_constraint_free(spec$m)
    gain <- if (name == "M1") gain_m1 else gain_s2
    rule <- if (name == "M1") rule_m1 else NULL
    run <- run_joint(name, spec, gc, gain, rule)

    ref <- d02$designs[[name]]
    outer_best <- max(ref$outer_df$gain)
    outer_at <- ref$outer_df$node[which.max(ref$outer_df$gain)]
    profile_best <- max(ref$profile_df$P)
    profile_at <- ref$profile_df$node[which.max(ref$profile_df$P)]
    reference <- max(outer_best, profile_best)
    joint_gain <- run$result$power$trial_success
    check1 <- joint_gain >= reference - tol
    say(
        "[%s] (1) d02 best outer-loop gain %.6f (gamma %g); best d02 profile %.6f (gamma %g); joint - max = %+.2e -> %s%s",
        name,
        outer_best,
        outer_at,
        profile_best,
        profile_at,
        joint_gain - reference,
        if (check1) "PASS" else "FAIL",
        if (smoke) " (smoke: different draws, not meaningful)" else ""
    )

    total <- built$secs + run$joint_secs
    say(
        "[%s] wall time: spec %.1f s + optimiser %.1f s = %.1f min",
        name,
        built$secs,
        run$joint_secs,
        total / 60
    )
    if (name == "S2") {
        say(
            "[S2] design 4.2 projection 3.5 to 6.5 min at m = 3, K = 2, nsim 1e5, d = 1: %s",
            if (total / 60 > 2 * 6.5) "> 2x the projection: PAUSE" else "within 2x"
        )
    }

    c(run, list(
        spec_secs = built$secs,
        reference = reference,
        outer_best = outer_best,
        profile_best = profile_best,
        check1 = check1
    ))
}


# ---- d = 2 coarse part: M1 with both hypotheses tagged -----------------------

run_d2 <- function() {
    name <- "D2"
    raw <- simulate_m1(nsim_d2)
    spending <- list(
        sf_hsd_opt(id = "g1", bounds = c(-2, 0), step = 0.5),
        sf_hsd_opt(id = "g2", bounds = c(-2, 0), step = 0.5)
    )
    built <- elapsed(transform_pvalues_gsd(raw, spending = spending))
    spec <- built$value
    say(
        "[%s] spec: %s nodes, %d tables, start (%s); built in %.1f s",
        name,
        paste(lengths(spec$param$lattice), collapse = " x "),
        length(spec$node_tables),
        toString(spec$param$start),
        built$secs
    )

    gc <- graph_constraint_free(2)
    run <- run_joint(name, spec, gc, gain_m1, rule_m1)

    # the product grid of fixed transforms (resolve() is identical to the
    # fixed transform at a node, P1 gate) and a default optimiser run per pair
    lattice <- spec$param$lattice
    pairs <- expand.grid(g1 = lattice[[1]], g2 = lattice[[2]])
    if (!full_grid) {
        start_q <- mapply(match, spec$param$start, lattice)
        q1 <- match(pairs$g1, lattice[[1]])
        q2 <- match(pairs$g2, lattice[[2]])
        near <- abs(q1 - start_q[[1]]) <= 1 & abs(q2 - start_q[[2]]) <= 1
        pairs <- pairs[near, , drop = FALSE]
    }
    pairs$gain <- NA_real_
    pairs$secs <- NA_real_
    for (p in seq_len(nrow(pairs))) {
        theta <- c(g1 = pairs$g1[[p]], g2 = pairs$g2[[p]])
        one <- elapsed(withr::with_seed(1, {
            graph_optimise_gsd(
                resolve_pvalues_gsd(spec, theta),
                graph_constraint = gc,
                trial_success = gain_m1,
                stopping_rule = rule_m1,
                num_threads = n_threads,
                control = control
            )
        }))
        pairs$gain[[p]] <- one$value$power$trial_success
        pairs$secs[[p]] <- one$secs
        say(
            "[%s] grid (%g, %g): gain %.6f; %.1f s",
            name,
            theta[[1]],
            theta[[2]],
            pairs$gain[[p]],
            one$secs
        )
    }
    best <- which.max(pairs$gain)
    joint_gain <- run$result$power$trial_success
    check2 <- joint_gain >= pairs$gain[[best]] - tol
    say(
        "[%s] (2) product-grid max %.6f at (%g, %g) over %d pairs; joint - max = %+.2e -> %s",
        name,
        pairs$gain[[best]],
        pairs$g1[[best]],
        pairs$g2[[best]],
        nrow(pairs),
        joint_gain - pairs$gain[[best]],
        if (check2) "PASS" else "FAIL"
    )

    c(run, list(spec_secs = built$secs, grid = pairs, check2 = check2))
}


# ---- run ---------------------------------------------------------------------

results <- list()
for (part in parts) {
    say("")
    say("== part %s ==", part)
    results[[part]] <- switch(
        part,
        M1 = run_d1("M1"),
        S2 = run_d1("S2"),
        D2 = run_d2(),
        stop("unknown part ", part)
    )
}

say("")
say("== summary ==")
for (part in names(results)) {
    r <- results[[part]]
    checks <- c(
        "(1)" = r$check1 %||% NA,
        "(2)" = r$check2 %||% NA,
        "(3)" = r$check3
    )
    checks <- checks[!is.na(checks)]
    say(
        "%s: %s",
        part,
        toString(paste(names(checks), ifelse(checks, "PASS", "FAIL")))
    )
}

# results without the (large) optimiser outputs
slim <- lapply(results, function(r) {
    r$result <- r$result[c(
        "hyp_weight", "trans_matrix", "power", "solution", "spending_param",
        "hyp_spending_param", "spending", "spending_profile",
        "spending_search", "profile_converged"
    )]
    r
})
saveRDS(slim, paste0(out_stem, "_results.rds"))
say("written: %s, %s", log_file, paste0(out_stem, "_results.rds"))
close(log_con)
