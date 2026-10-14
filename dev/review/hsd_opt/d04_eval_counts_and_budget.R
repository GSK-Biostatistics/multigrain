# d04_eval_counts_and_budget.R
# Phase D, item 4 (design section 4.3 item 4; brief D-b).
#   counts  - Table D6: on the section-2 design (sfLDOF on all three), the
#             internal drivers of graph_optimise_gsd() in its order with a
#             counter-wrapped gain, seeds 1 and 2 and global_search = FALSE,
#             checked against the one-call result; Table D7: warm-started
#             COBYLA counts E_C at gamma in {-2.25, -1.75, -3, -1, -4, 1} from
#             the gamma = -2 optimum.
#   bank_M0, bank_M1
#           - Table D8: prototype bank closure (graph genes + one parameter
#             gene u, node q = min(21, 1 + floor(u * 21)), H2's two columns
#             copied from slice q into a 5e4 x 4 working matrix) under GA::ga()
#             exactly as .graph_optimise_ga_gsd() calls it, against the plain
#             GA driver at gamma = -2; 5 seeds each (set.seed(100 + s)); one
#             extra configuration with the gene at popSize = 300.
#   report  - merges the parts into d04_eval_counts_and_budget_results.rds and
#             prints the tables (redirect to *_output.txt).
# The tracemem()/copy numbers are not repeated: d03 is in
# (d03_memory_and_copy_output.txt).
#
# Run from the repository root:
#   Rscript dev/review/hsd_opt/d04_eval_counts_and_budget.R    # all + report

source(file.path(
    if (dir.exists("dev/review/hsd_opt")) "dev/review/hsd_opt" else ".",
    "d_helpers.R"
))

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) {
    args <- c("counts", "bank_M0", "bank_M1", "report")
}
part_file <- function(nm) {
    file.path(here, sprintf("d04_eval_counts_and_budget_results_%s.rds", nm))
}
stamp <- function(...) {
    cat(format(Sys.time(), "%H:%M:%S"), " ", sprintf(...), "\n", sep = "")
    flush.console()
}
elapsed <- function() proc.time()[["elapsed"]]

ec_nodes <- c(-2.25, -1.75, -3, -1, -4, 1)
n_full <- 1e5
n_ga <- 5e4
bank_seeds <- 1:5
bank_seeds_pop300 <- 1:3

# ---- counts: D6 and D7 -------------------------------------------------------

staged_run <- function(pv, gc, gain, seed, global_search = TRUE) {
    gain_c <- counting_gain(gain)
    control <- multigrain:::control_prepare_dims(
        multigrain_control(),
        nsim = pv$nsim,
        m = pv$m,
        verbose = "silent"
    )
    rows <- list()
    ga_result <- NULL
    x0_for_local <- NULL
    set.seed(seed)
    if (global_search) {
        t0 <- elapsed()
        ga_result <- multigrain:::.graph_optimise_ga_gsd(
            pvals = pv,
            graph_constraint = gc,
            trial_success = gain_c,
            alpha = alpha,
            nsim = control$nsim_global,
            start_graph = list(list(hyp_weight = NULL, trans_matrix = NULL)),
            global_opts = control$global_opt,
            num_threads = n_threads,
            stop_rule = NULL,
            verbose = "silent"
        )
        rows$GA <- c(
            iter = ga_result$ga_output@iter,
            evals = count_of(gain_c),
            secs = elapsed() - t0
        )
        x0_for_local <- pmin(pmax(ga_result$ga_output@solution[1, ], 0), 1)
    }
    reset_count(gain_c)
    t0 <- elapsed()
    local_result <- multigrain:::.graph_optimise_local_gsd(
        pvals = pv,
        graph_constraint = gc,
        trial_success = gain_c,
        alpha = alpha,
        nsim = control$nsim_local,
        local_opts = control$local_opt,
        num_threads = n_threads,
        x0 = x0_for_local,
        stop_rule = NULL,
        verbose = "silent"
    )
    rows$local <- c(
        iter = local_result$local_output$iterations,
        evals = count_of(gain_c),
        secs = elapsed() - t0
    )
    best <- multigrain:::choose_graph(ga_result, local_result)
    reset_count(gain_c)
    t0 <- elapsed()
    clean <- multigrain:::prune_graph_gsd(
        pvals = pv,
        hyp_weight = best$hyp_weight,
        trans_matrix = best$trans_matrix,
        alpha = alpha,
        graph_constraint = gc,
        trial_success = gain_c,
        gamma = 1,
        verbose = "silent",
        stop_rule = NULL
    )
    rows$prune <- c(iter = NA, evals = count_of(gain_c), secs = elapsed() - t0)
    reset_count(gain_c)
    t0 <- elapsed()
    final <- calc_power_pvals_gsd(
        pvals = pv,
        hyp_weight = clean$hyp_weight,
        trans_matrix = clean$trans_matrix,
        alpha = alpha,
        custom_power = list(trial_success = gain_c),
        stopping_rule = NULL
    )
    rows$final <- c(iter = NA, evals = count_of(gain_c), secs = elapsed() - t0)
    tab <- do.call(rbind, rows)
    list(
        table = data.frame(stage = rownames(tab), tab, row.names = NULL),
        hyp_weight = clean$hyp_weight,
        trans_matrix = clean$trans_matrix,
        gain = final$trial_success,
        source = best$source,
        local_ga = c(
            ga = if (is.null(ga_result)) NA else ga_result$ga_trial_success,
            local = local_result$local_trial_success
        )
    )
}

run_counts <- function() {
    des <- make_design("S2")
    t0 <- elapsed()
    pv <- transform_pvalues_gsd(
        des$raw,
        spending = gsDesign::sfLDOF,
        alpha = alpha
    )
    tr_secs <- elapsed() - t0
    stamp("[counts] transform sfLDOF: %.2f s", tr_secs)

    d6 <- list()
    for (cfg in list(
        list(seed = 1, gs = TRUE, one_call = TRUE),
        list(seed = 2, gs = TRUE, one_call = FALSE),
        list(seed = 1, gs = FALSE, one_call = TRUE)
    )) {
        lab <- sprintf("seed %d, global_search = %s", cfg$seed, cfg$gs)
        t0 <- elapsed()
        st <- staged_run(pv, des$gc, des$gain, cfg$seed, cfg$gs)
        staged_secs <- elapsed() - t0
        if (!cfg$one_call) {
            d6[[lab]] <- list(
                st = st,
                staged_secs = staged_secs,
                one_secs = NA,
                one_gain = NA,
                identical_hyp_weight = NA,
                identical_graph = NA
            )
            stamp(
                "[counts] %s: staged %.1f s (one-call check skipped for wall time)",
                lab,
                staged_secs
            )
            saveRDS(list(d6 = d6), part_file("counts_partial"))
            next
        }
        set.seed(cfg$seed)
        t0 <- elapsed()
        one <- graph_optimise_gsd(
            pv,
            des$gc,
            des$gain,
            global_search = cfg$gs,
            num_threads = n_threads,
            control = multigrain_control(),
            verbose = "silent"
        )
        one_secs <- elapsed() - t0
        same <- identical(unname(one$hyp_weight), unname(st$hyp_weight)) &&
            identical(unname(one$trans_matrix), unname(st$trans_matrix))
        d6[[lab]] <- list(
            st = st,
            staged_secs = staged_secs,
            one_secs = one_secs,
            one_gain = one$power$trial_success,
            identical_hyp_weight = identical(
                unname(one$hyp_weight),
                unname(st$hyp_weight)
            ),
            identical_graph = same
        )
        stamp(
            "[counts] %s: staged %.1f s, one-call %.1f s, identical hyp_weight %s, gain %.6f",
            lab,
            staged_secs,
            one_secs,
            d6[[lab]]$identical_hyp_weight,
            one$power$trial_success
        )
        saveRDS(list(d6 = d6), part_file("counts_partial"))
    }

    # ---- D7: E_C from the gamma = -2 optimum --------------------------------
    d02_file <- file.path(here, "d02_warm_start_profile_results_S2.rds")
    if (file.exists(d02_file)) {
        d02 <- readRDS(d02_file)
        g2 <- d02$graphs[[match(-2, d02$nodes)]]
        g2_source <- "d02 S2 profile optimum at gamma -2"
    } else {
        tn <- transform_node(des, -2)
        set.seed(1)
        r2 <- graph_optimise_gsd(
            tn$pv,
            des$gc,
            des$gain,
            num_threads = n_threads,
            verbose = "silent"
        )
        g2 <- list(w = unname(r2$hyp_weight), G = unname(r2$trans_matrix))
        g2_source <- "default run at gamma -2, seed 1"
    }
    stamp("[counts] E_C start graph: %s", g2_source)
    local_opts <- local_opts_default(n_full, 3L)
    d7 <- list()
    for (g in c(-2, ec_nodes)) {
        tn <- transform_node(des, g)
        t0 <- elapsed()
        r <- refine_node(des, tn$pv, list(g2), NULL, local_opts)
        d7[[length(d7) + 1L]] <- data.frame(
            node = g,
            evals = r$n_eval,
            cobyla_s = r$cobyla_secs,
            node_s = elapsed() - t0,
            transform_s = tn$secs,
            gain_start = r$cand_vals[[1]],
            gain_refined = r$refined,
            gain = r$gain,
            ms_per_eval = 1000 * r$cobyla_secs / r$n_eval
        )
        stamp(
            "[counts] E_C node %.2f: %d evaluations, %.1f s",
            g,
            as.integer(r$n_eval),
            r$cobyla_secs
        )
    }
    d7 <- do.call(rbind, d7)
    saveRDS(
        list(
            d6 = d6,
            d7 = d7,
            g2 = g2,
            g2_source = g2_source,
            transform_ldof_s = tr_secs
        ),
        part_file("counts")
    )
    unlink(part_file("counts_partial"))
}

# ---- bank: D8 ------------------------------------------------------------------

# Prototype bank closure. The decode/penalty block is copied from
# create_obj_func_gsd(); the working matrix is built inside the factory (d03:
# then a column assignment with <<- does not copy the matrix).
make_bank_obj <- function(
    gc,
    gain,
    fixed_kmat,
    bank,
    opt_cols,
    p_graph,
    K,
    stop_rule
) {
    hyp_constraint <- gc$hyp_constraint
    trans_constraint <- gc$trans_constraint
    power_criterion <- gain$func
    M <- length(bank)
    work <- fixed_kmat + 0 # fresh copy owned by this closure
    archive <- new.env()
    archive$best <- rep(-Inf, M)
    archive$visits <- integer(M)
    archive$n_eval <- 0
    archive$n_kernel <- 0
    fn <- function(x) {
        archive$n_eval <- archive$n_eval + 1
        xg <- x[seq_len(p_graph)]
        u <- min(max(x[[p_graph + 1L]], 0), 1)
        q <- min(M, 1L + floor(u * M))

        theta <- multigrain:::split_theta(xg, hyp_constraint)
        hyp_weight <- multigrain:::recover_full_weights(
            theta$w_pars,
            hyp_constraint
        )
        trans_matrix <- multigrain:::recover_full_trans_matrix(
            theta$g_pars,
            trans_constraint
        )

        if (anyNA(hyp_weight) || anyNA(trans_matrix)) {
            return(-1e6)
        }
        if (any(hyp_weight < 0)) {
            return(sum(hyp_weight[hyp_weight < 0]))
        }
        if (any(trans_matrix < 0)) {
            return(sum(trans_matrix[trans_matrix < 0]))
        }
        if (any(hyp_weight > 1)) {
            return(-sum(hyp_weight[hyp_weight > 1]))
        }
        if (any(trans_matrix > 1)) {
            return(-sum(trans_matrix[trans_matrix > 1]))
        }

        hyp_weight[hyp_weight < 1e-4] <- 0
        trans_matrix[trans_matrix < 1e-5] <- 0

        work[, opt_cols] <<- bank[[q]]
        res <- multigrain:::graph_shortcut_gsd_parallel(
            pvals = work,
            alpha = alpha,
            w = hyp_weight,
            G = trans_matrix,
            K = K,
            num_threads = n_threads,
            grain_size = -1L,
            stop_rule = stop_rule
        )
        val <- power_criterion(res$time)
        archive$n_kernel <- archive$n_kernel + 1
        archive$visits[[q]] <- archive$visits[[q]] + 1L
        if (val > archive$best[[q]]) {
            archive$best[[q]] <- val
        }
        val
    }
    attr(fn, "archive") <- archive
    fn
}

run_bank <- function(name) {
    des <- make_design(name)
    nodes <- lattice
    M <- length(nodes)
    t0 <- elapsed()
    pvs <- lapply(nodes, function(g) transform_node(des, g)$pv)
    tr_secs <- elapsed() - t0
    stamp("[bank %s] %d transforms %.1f s", name, M, tr_secs)
    kfull <- lapply(pvs, multigrain:::.gsd_kernel_matrix)
    K <- pvs[[1]]$K
    m <- pvs[[1]]$m
    opt_cols <- as.vector(outer(des$opt_cols_hyp, (seq_len(K) - 1L) * m, `+`))
    opt_cols <- sort(opt_cols)
    stop_rule <- stop_table(des$rule, pvs[[1]])

    set.seed(11)
    idx_ga <- sample.int(n_full, n_ga)
    t0 <- elapsed()
    bank <- lapply(kfull, function(km) km[idx_ga, opt_cols, drop = FALSE])
    fixed_kmat <- kfull[[1]][idx_ga, , drop = FALSE]
    bank_secs <- elapsed() - t0
    fixed_same <- all(vapply(
        kfull,
        function(km) identical(km[, -opt_cols], kfull[[1]][, -opt_cols]),
        logical(1)
    ))
    stamp(
        "[bank %s] bank built %.2f s; fixed columns identical across nodes: %s",
        name,
        bank_secs,
        fixed_same
    )

    gc <- des$gc
    x0 <- multigrain:::.build_start_matrix(
        gc,
        list(list(hyp_weight = NULL, trans_matrix = NULL))
    )
    p_graph <- ncol(x0)
    q_start <- match(-2, nodes)
    u_start <- (q_start - 0.5) / M
    base_opts <- multigrain:::control_prepare_dims(
        multigrain_control(),
        nsim = n_full,
        m = m,
        verbose = "silent"
    )$global_opt
    pop_rule <- min(max(40L * (m + 1L), 200L), 500L)

    run_gene <- function(seed, pop) {
        fit <- make_bank_obj(
            gc,
            des$gain,
            fixed_kmat,
            bank,
            opt_cols,
            p_graph,
            K,
            stop_rule
        )
        immutable <- list(
            type = "real-valued",
            fitness = fit,
            lower = rep(0, p_graph + 1L),
            upper = rep(1, p_graph + 1L),
            population = multigrain:::.cauchy_population,
            mutation = multigrain:::.make_cauchy_mutation_multi(
                p_param_mutate = 0.1,
                scale = 1.0
            ),
            optim = TRUE,
            suggestions = cbind(x0, u_start)
        )
        opts <- base_opts
        opts$popSize <- pop
        ga_args <- utils::modifyList(immutable, opts)
        set.seed(100 + seed)
        t0 <- elapsed()
        ga <- do.call(GA::ga, ga_args)
        secs <- elapsed() - t0
        sol <- ga@solution[1, ]
        u <- min(max(sol[[p_graph + 1L]], 0), 1)
        q <- min(M, 1L + floor(u * M))
        gr <- multigrain:::param_to_solution(
            sol[seq_len(p_graph)],
            gc,
            process = TRUE
        )
        full <- eval_graph(
            kfull[[q]],
            K,
            des$gain,
            gr$hyp_weight,
            gr$trans_matrix,
            stop_rule
        )
        arc <- attr(fit, "archive")
        list(
            config = sprintf("gene, popSize %d", pop),
            seed = seed,
            iter = ga@iter,
            evals = arc$n_eval,
            kernel_evals = arc$n_kernel,
            secs = secs,
            node = nodes[[q]],
            w1 = gr$hyp_weight[[1]],
            best_fitness = ga@fitnessValue,
            full_gain = full,
            nodes_visited = sum(arc$visits > 0),
            visits = arc$visits,
            archive_best = arc$best
        )
    }
    run_plain <- function(seed) {
        gain_c <- counting_gain(des$gain)
        opts <- base_opts
        set.seed(100 + seed)
        t0 <- elapsed()
        ga <- multigrain:::.graph_optimise_ga_gsd(
            pvals = pvs[[q_start]],
            graph_constraint = gc,
            trial_success = gain_c,
            alpha = alpha,
            nsim = n_ga,
            start_graph = list(list(hyp_weight = NULL, trans_matrix = NULL)),
            global_opts = opts,
            num_threads = n_threads,
            stop_rule = stop_rule,
            verbose = "silent"
        )
        secs <- elapsed() - t0
        list(
            config = sprintf("no gene (gamma -2), popSize %d", opts$popSize),
            seed = seed,
            iter = ga$ga_output@iter,
            evals = count_of(gain_c),
            kernel_evals = NA,
            secs = secs,
            node = -2,
            w1 = ga$ga_hyp_weight[[1]],
            best_fitness = ga$ga_output@fitnessValue,
            full_gain = ga$ga_trial_success,
            nodes_visited = 1,
            visits = NULL,
            archive_best = NULL
        )
    }

    runs <- list()
    for (s in bank_seeds) {
        for (kind in c("gene", "plain")) {
            r <- if (kind == "gene") {
                run_gene(s, base_opts$popSize)
            } else {
                run_plain(s)
            }
            runs[[length(runs) + 1L]] <- r
            stamp(
                "[bank %s] %s seed %d: iter %d, evals %d, %.1f s, node %.2f, best %.6f, full %.6f",
                name,
                r$config,
                s,
                as.integer(r$iter),
                as.integer(r$evals),
                r$secs,
                r$node,
                r$best_fitness,
                r$full_gain
            )
            saveRDS(
                list(runs = runs),
                part_file(paste0("bank_", name, "_partial"))
            )
        }
    }
    for (s in bank_seeds_pop300) {
        r <- run_gene(s, 300L)
        runs[[length(runs) + 1L]] <- r
        stamp(
            "[bank %s] %s seed %d: iter %d, evals %d, %.1f s, node %.2f, best %.6f, full %.6f",
            name,
            r$config,
            s,
            as.integer(r$iter),
            as.integer(r$evals),
            r$secs,
            r$node,
            r$best_fitness,
            r$full_gain
        )
        saveRDS(list(runs = runs), part_file(paste0("bank_", name, "_partial")))
    }
    saveRDS(
        list(
            design = name,
            runs = runs,
            transform_s = tr_secs,
            bank_s = bank_secs,
            fixed_same = fixed_same,
            p_graph = p_graph,
            pop_rule = pop_rule,
            popSize_default = base_opts$popSize,
            u_start = u_start,
            seeds = list(idx_ga = 11, ga = 100 + bank_seeds)
        ),
        part_file(paste0("bank_", name))
    )
    unlink(part_file(paste0("bank_", name, "_partial")))
}

# ---- report --------------------------------------------------------------------

report <- function() {
    options(width = 200)
    cat("# d04_eval_counts_and_budget.R output\n")
    cat(sprintf(
        "R %s; multigrain %s; GA %s; nloptr %s; %s; num_threads = %d; nsim = 1e5\n",
        getRversion(),
        packageVersion("multigrain"),
        packageVersion("GA"),
        packageVersion("nloptr"),
        format(Sys.time()),
        n_threads
    ))
    cat(
        "tracemem()/copy numbers: not repeated, d03 is in (see d03_memory_and_copy_output.txt)\n"
    )
    out <- list()
    if (file.exists(part_file("counts"))) {
        cn <- readRDS(part_file("counts"))
        out$counts <- cn
        cat(
            "\n== Table D6: section-2 design (sfLDOF x 3), drivers in graph_optimise_gsd() order ==\n"
        )
        cat(sprintf(
            "transform (sfLDOF, G = 1024): %.2f s\n",
            cn$transform_ldof_s
        ))
        for (lab in names(cn$d6)) {
            x <- cn$d6[[lab]]
            tab <- x$st$table
            tab$secs <- sprintf("%.2f", tab$secs)
            cat(sprintf("\n-- %s --\n", lab))
            print(tab, row.names = FALSE)
            cat(sprintf(
                "staged total %.1f s; one-call graph_optimise_gsd() %s s; identical(hyp_weight) %s; identical(trans_matrix too) %s\n",
                x$staged_secs,
                if (is.na(x$one_secs)) {
                    "not run"
                } else {
                    sprintf("%.1f", x$one_secs)
                },
                x$identical_hyp_weight,
                x$identical_graph
            ))
            cat(sprintf(
                "gain staged %.6f, one-call %s; source %s; GA full-sample %s, local full-sample %.6f\n",
                x$st$gain,
                if (is.na(x$one_gain)) "-" else sprintf("%.6f", x$one_gain),
                x$st$source,
                if (is.na(x$st$local_ga[["ga"]])) {
                    "-"
                } else {
                    sprintf("%.6f", x$st$local_ga[["ga"]])
                },
                x$st$local_ga[["local"]]
            ))
            ga_row <- tab[tab$stage == "GA", ]
            if (nrow(ga_row)) {
                cat(sprintf(
                    "GA: %.2f ms per evaluation (incl. GA overhead); GA share of the staged run %.0f %%\n",
                    1000 * as.numeric(ga_row$secs) / ga_row$evals,
                    100 * as.numeric(ga_row$secs) / x$staged_secs
                ))
            }
        }
        cat(
            "\n== Table D7: warm-started COBYLA counts E_C, section-2 design, from the gamma = -2 optimum ==\n"
        )
        cat("start graph:", cn$g2_source, "\n")
        d7 <- cn$d7
        d7p <- d7
        for (cc in c("gain_start", "gain_refined", "gain")) {
            d7p[[cc]] <- sprintf("%.6f", d7[[cc]])
        }
        for (cc in c("cobyla_s", "node_s", "transform_s", "ms_per_eval")) {
            d7p[[cc]] <- sprintf("%.2f", d7[[cc]])
        }
        print(d7p, row.names = FALSE)
        ec <- d7$evals[d7$node != -2]
        cat(sprintf(
            "E_C over the six off-start nodes: median %d, range %d-%d; [derived] design 4.2: 200-800\n",
            as.integer(stats::median(ec)),
            as.integer(min(ec)),
            as.integer(max(ec))
        ))
    }
    d02_file <- file.path(here, "d02_warm_start_profile_results.rds")
    d02 <- if (file.exists(d02_file)) readRDS(d02_file) else NULL
    for (nm in c("M0", "M1")) {
        f <- part_file(paste0("bank_", nm))
        if (!file.exists(f)) {
            next
        }
        b <- readRDS(f)
        out[[paste0("bank_", nm)]] <- b
        cat(sprintf(
            "\n== Table D8 (%s): GA with and without the parameter gene ==\n",
            nm
        ))
        cat(sprintf(
            "21 transforms %.1f s; bank (21 x 5e4 x 2 slices) %.2f s; fixed H1 columns identical across nodes: %s\n",
            b$transform_s,
            b$bank_s,
            b$fixed_same
        ))
        cat(sprintf(
            "graph genes p_graph = %d; popSize rule min(max(40 (m + d), 200), 500) = %d at m + d = 3 (inert: default is %d)\n",
            b$p_graph,
            b$pop_rule,
            b$popSize_default
        ))
        argmax <- NA
        if (!is.null(d02) && !is.null(d02$designs[[nm]])) {
            p <- d02$designs[[nm]]$profile_df
            argmax <- p$node[which.max(p$P)]
            pl <- plateau_of(p$node, p$P)
            cat(sprintf(
                "d02 profile argmax %.2f (P~ %.6f), plateau [%.2f, %.2f]\n",
                argmax,
                max(p$P),
                pl[1],
                pl[2]
            ))
        }
        df <- do.call(
            rbind,
            lapply(b$runs, function(r) {
                data.frame(
                    config = r$config,
                    seed = r$seed,
                    generations = r$iter,
                    evals = r$evals,
                    wall_s = round(r$secs, 1),
                    node = r$node,
                    w1 = round(r$w1, 4),
                    best_fitness_sub = sprintf("%.6f", r$best_fitness),
                    full_gain = sprintf("%.6f", r$full_gain),
                    nodes_visited = r$nodes_visited
                )
            })
        )
        print(df, row.names = FALSE)
        agg <- do.call(
            rbind,
            lapply(
                split(b$runs, vapply(b$runs, `[[`, character(1), "config")),
                function(rs) {
                    data.frame(
                        config = rs[[1]]$config,
                        n = length(rs),
                        median_generations = stats::median(vapply(
                            rs,
                            `[[`,
                            numeric(1),
                            "iter"
                        )),
                        median_evals = stats::median(vapply(
                            rs,
                            `[[`,
                            numeric(1),
                            "evals"
                        )),
                        median_wall_s = round(
                            stats::median(vapply(rs, `[[`, numeric(1), "secs")),
                            1
                        ),
                        nodes_found = paste(
                            sprintf(
                                "%.2f",
                                vapply(rs, `[[`, numeric(1), "node")
                            ),
                            collapse = " "
                        ),
                        max_full_gain = sprintf(
                            "%.6f",
                            max(vapply(rs, `[[`, numeric(1), "full_gain"))
                        )
                    )
                }
            )
        )
        cat("\nsummary:\n")
        print(agg, row.names = FALSE)
        gen_gene <- agg$median_generations[
            agg$config == sprintf("gene, popSize %d", b$popSize_default)
        ]
        gen_plain <- agg$median_generations[grepl("^no gene", agg$config)]
        cat(sprintf(
            "generation ratio gene / no gene at popSize %d (medians): %.2f (rule threshold 1.5)\n",
            b$popSize_default,
            gen_gene / gen_plain
        ))
        gene_runs <- Filter(function(r) !is.null(r$visits), b$runs)
        vis <- Reduce(`+`, lapply(gene_runs, `[[`, "visits"))
        cat("visits per node, summed over gene runs (nodes -4 .. 1):\n")
        print(stats::setNames(vis, sprintf("%.2f", lattice)))
    }
    saveRDS(out, file.path(here, "d04_eval_counts_and_budget_results.rds"))
}

if (sys.nframe() == 0L) {
    for (a in args) {
        if (a == "report") {
            report()
        } else if (a == "counts") {
            run_counts()
        } else if (a %in% c("bank_M0", "bank_M1")) {
            run_bank(sub("bank_", "", a))
        } else {
            stop("unknown part ", a)
        }
    }
}
