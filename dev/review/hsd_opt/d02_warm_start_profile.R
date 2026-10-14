# d02_warm_start_profile.R
# Phase D, item 2 (design section 4.3 item 2; brief D-b): the measurement the
# local stage rests on. For M0, M1 (README 7.4 at nsim = 1e5) and the README
# section-2 design of the brief, compare
#   - an outer loop: transform at gamma, then a full default
#     graph_optimise_gsd(), at gamma in {-4, -3, -2, -1, 0, 1} x seeds {1, 2, 3}
#   - a warm-started profile over the 21-node lattice seq(-4, 1, 0.25):
#     COBYLA at each node from the better of the previous node's refined
#     optimum and the overall best graph, walking outward from the node of the
#     best seed-1 outer result.
# Gate D4: max over seeds of U_full - P~ <= 5e-4 at every outer node inside
# the plateau of P~ (5e-4 around its argmax).
#
# Run from the repository root:
#   Rscript dev/review/hsd_opt/d02_warm_start_profile.R            # all + report
#   Rscript dev/review/hsd_opt/d02_warm_start_profile.R M0 M1      # parts
#   Rscript dev/review/hsd_opt/d02_warm_start_profile.R report     # tables only
# Parts write d02_warm_start_profile_results_<design>.rds; the report merges
# them into d02_warm_start_profile_results.rds (output to *_output.txt).

source(file.path(
    if (dir.exists("dev/review/hsd_opt")) "dev/review/hsd_opt" else ".",
    "d_helpers.R"
))

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) {
    args <- c("M0", "M1", "S2", "report")
}
part_file <- function(nm) {
    file.path(here, sprintf("d02_warm_start_profile_results_%s.rds", nm))
}
stamp <- function(...) {
    cat(format(Sys.time(), "%H:%M:%S"), " ", sprintf(...), "\n", sep = "")
    flush.console()
}

run_design <- function(name) {
    t_design <- proc.time()[["elapsed"]]
    des <- make_design(name)
    stamp("[%s] design built", name)

    # ---- 1. one transform per lattice node (cached) -----------------------
    nodes <- lattice
    pv <- vector("list", length(nodes))
    tr_secs <- numeric(length(nodes))
    for (j in seq_along(nodes)) {
        tn <- transform_node(des, nodes[[j]])
        pv[[j]] <- tn$pv
        tr_secs[[j]] <- tn$secs
    }
    stamp(
        "[%s] %d transforms: %.1f s (median %.2f s per node)",
        name,
        length(nodes),
        sum(tr_secs),
        stats::median(tr_secs)
    )
    stop_rule <- stop_table(des$rule, pv[[1]])
    local_opts <- local_opts_default(pv[[1]]$nsim, des$m)

    # ---- 2. outer-loop reference ------------------------------------------
    # Seeds per design: D02_SEEDS_<design>="1,2" cuts repetitions (S2 was cut
    # to seeds 1 and 2 for wall time; see the discovery report). A partial
    # file from an interrupted run is resumed: finished (node, seed) runs are
    # kept, not repeated.
    seeds <- as.integer(strsplit(Sys.getenv(paste0("D02_SEEDS_", name), "1,2,3"), ",")[[1]])
    outer <- list()
    pf <- part_file(paste0(name, "_partial"))
    if (file.exists(pf)) {
        outer <- readRDS(pf)$outer
        stamp("[%s] resuming: %d outer runs already done", name, length(outer))
    }
    done <- vapply(outer, function(o) paste(o$node, o$seed), character(1))
    for (g in outer_nodes) {
        j <- match(g, nodes)
        for (s in seeds) {
            if (paste(g, s) %in% done) next
            set.seed(s)
            t0 <- proc.time()[["elapsed"]]
            res <- graph_optimise_gsd(
                pv[[j]],
                des$gc,
                des$gain,
                stopping_rule = des$rule,
                num_threads = n_threads,
                verbose = "silent"
            )
            secs <- proc.time()[["elapsed"]] - t0
            outer[[length(outer) + 1L]] <- list(
                node = g,
                seed = s,
                gain = res$power$trial_success,
                w = unname(res$hyp_weight),
                G = unname(res$trans_matrix),
                secs = secs,
                opt_source = res$solution$opt_source,
                ga_iter = res$global_output@iter
            )
            stamp(
                "[%s] outer node %5.2f seed %d: gain %.6f w1 %.4f (%s) %.1f s",
                name,
                g,
                s,
                res$power$trial_success,
                res$hyp_weight[[1]],
                res$solution$opt_source,
                secs
            )
            saveRDS(list(outer = outer), part_file(paste0(name, "_partial")))
        }
    }
    outer_df <- data.frame(
        node = vapply(outer, `[[`, numeric(1), "node"),
        seed = vapply(outer, `[[`, numeric(1), "seed"),
        gain = vapply(outer, `[[`, numeric(1), "gain"),
        w1 = vapply(outer, function(o) o$w[[1]], numeric(1)),
        secs = vapply(outer, `[[`, numeric(1), "secs"),
        opt_source = vapply(outer, `[[`, character(1), "opt_source"),
        ga_iter = vapply(outer, `[[`, numeric(1), "ga_iter")
    )

    # ---- 3. warm-started profile ------------------------------------------
    s1 <- outer_df$seed == 1
    k_best <- which(s1)[which.max(outer_df$gain[s1])]
    start_node <- outer_df$node[[k_best]]
    j0 <- match(start_node, nodes)
    start_graph <- list(w = outer[[k_best]]$w, G = outer[[k_best]]$G)
    stamp("[%s] profile starts at node %.2f", name, start_node)

    prof <- vector("list", length(nodes))
    node_secs <- numeric(length(nodes))

    t0 <- proc.time()[["elapsed"]]
    prof[[j0]] <- refine_node(des, pv[[j0]], list(start_graph), stop_rule, local_opts)
    node_secs[[j0]] <- proc.time()[["elapsed"]] - t0
    overall <- prof[[j0]]$graph
    S <- rep(NA_real_, length(nodes))
    S[[j0]] <- eval_graph(
        multigrain:::.gsd_kernel_matrix(pv[[j0]]),
        pv[[j0]]$K,
        des$gain,
        overall$w,
        overall$G,
        stop_rule
    )

    walk <- function(idx_seq) {
        prev <- overall
        for (j in idx_seq) {
            t0 <- proc.time()[["elapsed"]]
            r <- refine_node(
                des,
                pv[[j]],
                list(prev = prev, overall = overall),
                stop_rule,
                local_opts
            )
            node_secs[[j]] <<- proc.time()[["elapsed"]] - t0
            prof[[j]] <<- r
            S[[j]] <<- r$cand_vals[[2]]
            prev <- r$graph
            stamp(
                "[%s] profile node %5.2f: P~ %.6f S %.6f n_eval %d %.1f s",
                name,
                nodes[[j]],
                r$gain,
                S[[j]],
                as.integer(r$n_eval),
                node_secs[[j]]
            )
        }
    }
    if (j0 < length(nodes)) walk(seq(j0 + 1L, length(nodes)))
    if (j0 > 1L) walk(seq(j0 - 1L, 1L))

    profile_df <- data.frame(
        node = nodes,
        P = vapply(prof, `[[`, numeric(1), "gain"),
        S = S,
        refined = vapply(prof, `[[`, numeric(1), "refined"),
        start_from = vapply(prof, `[[`, numeric(1), "start_from"),
        w1 = vapply(prof, function(r) r$graph$w[[1]], numeric(1)),
        n_eval = vapply(prof, `[[`, numeric(1), "n_eval"),
        transform_s = tr_secs,
        refine_s = node_secs,
        cobyla_s = vapply(prof, `[[`, numeric(1), "cobyla_secs")
    )
    graphs <- lapply(prof, `[[`, "graph")

    out <- list(
        design = name,
        nodes = nodes,
        outer = outer,
        outer_df = outer_df,
        start_node = start_node,
        start_graph = start_graph,
        profile_df = profile_df,
        graphs = graphs,
        overall = overall,
        local_opts = local_opts,
        seeds = list(sim = if (name == "S2") 1 else 20261007, outer = sort(unique(outer_df$seed))),
        secs = proc.time()[["elapsed"]] - t_design
    )
    saveRDS(out, part_file(name))
    unlink(part_file(paste0(name, "_partial")))
    stamp("[%s] done in %.1f s", name, out$secs)
    invisible(out)
}

report <- function() {
    options(width = 200)
    designs <- c("M0", "M1", "S2")
    have <- designs[file.exists(vapply(designs, part_file, character(1)))]
    res <- lapply(stats::setNames(have, have), function(nm) readRDS(part_file(nm)))
    cat("# d02_warm_start_profile.R output\n")
    cat(sprintf(
        "R %s; multigrain %s; %s; num_threads = %d; nsim = 1e5\n",
        getRversion(),
        packageVersion("multigrain"),
        format(Sys.time()),
        n_threads
    ))
    cat("Designs reported:", paste(have, collapse = ", "), "\n")
    gate <- list()
    for (nm in have) {
        r <- res[[nm]]
        p <- r$profile_df
        o <- r$outer_df
        cat(sprintf("\n== %s ==\n", nm))
        cat(sprintf(
            "simulation seed %s; outer runs per node and seed listed below; profile start node %.2f (best seed-1 outer result)\n",
            r$seeds$sim,
            r$start_node
        ))
        cat("\n-- outer-loop reference (full default graph_optimise_gsd per node and seed) --\n")
        o_print <- o
        o_print$gain <- sprintf("%.6f", o$gain)
        o_print$w1 <- sprintf("%.4f", o$w1)
        o_print$secs <- sprintf("%.1f", o$secs)
        print(o_print, row.names = FALSE)

        cat("\n-- warm-started profile --\n")
        p_print <- data.frame(
            node = sprintf("%.2f", p$node),
            P_tilde = sprintf("%.6f", p$P),
            S_fixed_graph = sprintf("%.6f", p$S),
            P_minus_S = sprintf("%.2e", p$P - p$S),
            start = c("seed1", "prev", "overall")[ifelse(p$node == r$start_node, 1, p$start_from + 1)],
            w1 = sprintf("%.4f", p$w1),
            n_eval = p$n_eval,
            transform_s = sprintf("%.2f", p$transform_s),
            refine_s = sprintf("%.2f", p$refine_s)
        )
        print(p_print, row.names = FALSE)

        pl <- plateau_of(p$node, p$P)
        pl_s <- plateau_of(p$node, p$S)
        am <- p$node[which.max(p$P)]
        cat(sprintf(
            "argmax P~ = %.2f (P~ = %.6f); plateau(P~, 5e-4) = [%.2f, %.2f]; plateau(S, 5e-4) = [%.2f, %.2f]\n",
            am,
            max(p$P),
            pl[1],
            pl[2],
            pl_s[1],
            pl_s[2]
        ))

        cat("\n-- Table D4: gate --\n")
        d4 <- do.call(rbind, lapply(outer_nodes, function(g) {
            u <- o$gain[o$node == g]
            us <- o$seed[o$node == g]
            data.frame(
                node = g,
                P_tilde = p$P[p$node == g],
                S = p$S[p$node == g],
                U_seed1 = if (1 %in% us) u[us == 1] else NA_real_,
                U_seed2 = if (2 %in% us) u[us == 2] else NA_real_,
                U_seed3 = if (3 %in% us) u[us == 3] else NA_real_,
                max_U = max(u),
                gap = max(u) - p$P[p$node == g],
                in_plateau = g >= pl[1] & g <= pl[2]
            )
        }))
        d4_print <- d4
        for (cc in c("P_tilde", "S", "U_seed1", "U_seed2", "U_seed3", "max_U")) {
            d4_print[[cc]] <- sprintf("%.6f", d4[[cc]])
        }
        d4_print$gap <- sprintf("%+.2e", d4$gap)
        print(d4_print, row.names = FALSE)
        inside <- d4$in_plateau
        max_gap_in <- if (any(inside)) max(d4$gap[inside]) else NA_real_
        max_gap_out <- if (any(!inside)) max(d4$gap[!inside]) else NA_real_
        pass <- isTRUE(all(d4$gap[inside] <= 5e-4)) && any(inside)
        cat(sprintf(
            "max gap inside plateau = %+.2e (%d outer nodes inside); max gap outside = %s\n",
            max_gap_in,
            sum(inside),
            if (is.na(max_gap_out)) "none" else sprintf("%+.2e", max_gap_out)
        ))
        cat(sprintf("GATE %s: %s\n", nm, if (pass) "PASS" else "FAIL"))
        gate[[nm]] <- list(d4 = d4, pass = pass, plateau = pl, plateau_S = pl_s,
                           argmax = am, max_gap_in = max_gap_in,
                           max_gap_out = max_gap_out)

        cat("\n-- cost --\n")
        cat(sprintf(
            "outer loop: %d runs, %.1f s in all (median %.1f s per run; transforms extra)\n",
            nrow(o),
            sum(o$secs),
            stats::median(o$secs)
        ))
        cat(sprintf(
            "profile: 21 transforms %.1f s (median %.2f s) + refinement %.1f s (median %.2f s, %d evaluations in all, median %d per node); total %.1f s\n",
            sum(p$transform_s),
            stats::median(p$transform_s),
            sum(p$refine_s),
            stats::median(p$refine_s),
            as.integer(sum(p$n_eval)),
            as.integer(stats::median(p$n_eval)),
            sum(p$transform_s) + sum(p$refine_s)
        ))
        cat(sprintf(
            "COBYLA per evaluation (full nsim): %.2f ms\n",
            1000 * sum(p$cobyla_s) / sum(p$n_eval)
        ))
    }

    if ("M1" %in% have) {
        p <- res$M1$profile_df
        am <- p$node[which.max(p$P)]
        cat("\n== README 7.5 reproduction (M1) ==\n")
        cat(sprintf(
            "M1: P~ at gamma -4 (LDOF-like H2) = %.6f with w1 = %.4f; argmax gamma = %.2f with P~ = %.6f and w1 = %.4f; gain moved %+.4f\n",
            p$P[p$node == -4],
            p$w1[p$node == -4],
            am,
            max(p$P),
            p$w1[which.max(p$P)],
            max(p$P) - p$P[p$node == -4]
        ))
        if ("M0" %in% have) {
            p0 <- res$M0$profile_df
            cat(sprintf(
                "M0 for contrast: P~ at gamma -4 = %.6f (w1 %.4f); argmax gamma = %.2f with P~ = %.6f (w1 %.4f)\n",
                p0$P[p0$node == -4],
                p0$w1[p0$node == -4],
                p0$node[which.max(p0$P)],
                max(p0$P),
                p0$w1[which.max(p0$P)]
            ))
        }
        cat("README 7.5 (nsim 2e5, LDOF vs sfLDPocock on H2, grid in w1): 0.6707 at w1 0.32 -> 0.6903 at w1 0.46\n")
    }

    cat("\n== Gate summary ==\n")
    for (nm in names(gate)) {
        cat(sprintf(
            "%s: %s (max gap inside plateau %+.2e; plateau [%.2f, %.2f])\n",
            nm,
            if (gate[[nm]]$pass) "PASS" else "FAIL",
            gate[[nm]]$max_gap_in,
            gate[[nm]]$plateau[1],
            gate[[nm]]$plateau[2]
        ))
    }
    all_pass <- length(gate) == 3 && all(vapply(gate, `[[`, logical(1), "pass"))
    cat(sprintf("OVERALL WARM-START GATE: %s\n", if (all_pass) "PASS" else "FAIL or incomplete"))
    saveRDS(list(designs = res, gate = gate), file.path(here, "d02_warm_start_profile_results.rds"))
}

if (sys.nframe() == 0L) {
    for (a in args) {
        if (a == "report") {
            report()
        } else {
            run_design(a)
        }
    }
}
