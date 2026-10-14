# d05_grid_step_loss.R
# Phase D, item 5 (design section 4.3 item 5; brief D-b): the loss from a
# coarse lattice step. At the d02 profile argmax g* on M1 and the section-2
# design: fine grid seq(g* - 0.5, g* + 0.5, by = 0.025) (41 nodes), one
# transform and a warm-started COBYLA from the g* optimum at each; fine max F.
# For each coarse step D in {0.5, 0.25, 0.125} the coarse nodes are the fine
# nodes at multiples of D; loss = F - best coarse node gain. Table D9.
# Decision (plan section 3): coarsest D with loss <= 4e-5 on both designs,
# else 0.25.
#
# Needs d02_warm_start_profile_results_<design>.rds (M1, S2).
# Run from the repository root:
#   Rscript dev/review/hsd_opt/d05_grid_step_loss.R > dev/review/hsd_opt/d05_grid_step_loss_output.txt 2>&1

source(file.path(
    if (dir.exists("dev/review/hsd_opt")) "dev/review/hsd_opt" else ".",
    "d_helpers.R"
))
options(width = 200)

designs <- c("M1", "S2")
steps <- c(0.5, 0.25, 0.125)
loss_tol <- 4e-5
fine_step <- 0.025
half_width <- 0.5

t_script <- proc.time()[["elapsed"]]
cat("# d05_grid_step_loss.R output\n")
cat(sprintf(
    "R %s; multigrain %s; %s; num_threads = %d; nsim = 1e5\n",
    getRversion(),
    packageVersion("multigrain"),
    format(Sys.time()),
    n_threads
))

results <- list()
for (nm in designs) {
    d02 <- readRDS(file.path(
        here,
        sprintf("d02_warm_start_profile_results_%s.rds", nm)
    ))
    p <- d02$profile_df
    j_star <- which.max(p$P)
    g_star <- p$node[[j_star]]
    graph_star <- d02$graphs[[j_star]]
    des <- make_design(nm)
    fine <- round(
        seq(g_star - half_width, g_star + half_width, by = fine_step),
        10
    )
    local_opts <- local_opts_default(1e5, des$m)
    rows <- list()
    stop_rule <- NULL
    for (g in fine) {
        tn <- transform_node(des, g)
        if (is.null(stop_rule) && !is.null(des$rule)) {
            stop_rule <- stop_table(des$rule, tn$pv)
        }
        t0 <- proc.time()[["elapsed"]]
        r <- refine_node(des, tn$pv, list(graph_star), stop_rule, local_opts)
        rows[[length(rows) + 1L]] <- data.frame(
            node = g,
            gain = r$gain,
            fixed_graph = r$cand_vals[[1]],
            w1 = r$graph$w[[1]],
            n_eval = r$n_eval,
            transform_s = tn$secs,
            refine_s = proc.time()[["elapsed"]] - t0
        )
        cat(sprintf(
            "%s [%s] node %.3f: gain %.6f (fixed graph %.6f), n_eval %d\n",
            format(Sys.time(), "%H:%M:%S"),
            nm,
            g,
            r$gain,
            r$cand_vals[[1]],
            as.integer(r$n_eval)
        ))
        flush.console()
    }
    df <- do.call(rbind, rows)
    F_max <- max(df$gain)
    # inside the default bounds [-4, 1]: nodes outside them are not reachable
    # by any lattice step, so the decision uses the bounded window; the
    # unbounded window is reported beside it
    inb <- df$node >= -4 - 1e-9 & df$node <= 1 + 1e-9
    F_in <- max(df$gain[inb])
    d9 <- do.call(
        rbind,
        lapply(steps, function(D) {
            on <- abs(df$node / D - round(df$node / D)) < 1e-9
            best <- max(df$gain[on])
            on_in <- on & inb
            best_in <- max(df$gain[on_in])
            data.frame(
                design = nm,
                step = D,
                n_coarse = sum(on_in),
                best_coarse_node = df$node[on_in][which.max(df$gain[on_in])],
                best_coarse = best_in,
                loss = F_in - best_in,
                loss_unbounded_window = F_max - best
            )
        })
    )
    results[[nm]] <- list(
        g_star = g_star,
        P_star_d02 = p$P[[j_star]],
        fine = df,
        F_max = F_max,
        F_node = df$node[which.max(df$gain)],
        F_in = F_in,
        F_in_node = df$node[inb][which.max(df$gain[inb])],
        d9 = d9
    )
}

cat("\n== Fine grids ==\n")
for (nm in designs) {
    r <- results[[nm]]
    cat(sprintf(
        "\n-- %s: g* = %.2f (d02 P~ %.6f); fine grid [%.3f, %.3f] step %.3f; F = %.6f at %.3f; inside [-4, 1]: F = %.6f at %.3f --\n",
        nm,
        r$g_star,
        r$P_star_d02,
        min(r$fine$node),
        max(r$fine$node),
        fine_step,
        r$F_max,
        r$F_node,
        r$F_in,
        r$F_in_node
    ))
    fp <- r$fine
    fp$gain <- sprintf("%.6f", fp$gain)
    fp$fixed_graph <- sprintf("%.6f", fp$fixed_graph)
    fp$w1 <- sprintf("%.4f", fp$w1)
    fp$transform_s <- sprintf("%.2f", fp$transform_s)
    fp$refine_s <- sprintf("%.2f", fp$refine_s)
    print(fp, row.names = FALSE)
    cat(sprintf(
        "cost: transforms %.1f s, refinement %.1f s\n",
        sum(r$fine$transform_s),
        sum(r$fine$refine_s)
    ))
}

cat(
    "\n== Table D9: loss of a coarse lattice step (F - best coarse node within g* +/- 0.5, intersected with [-4, 1]) ==\n"
)
d9 <- do.call(rbind, lapply(results, `[[`, "d9"))
d9p <- d9
d9p$best_coarse <- sprintf("%.6f", d9$best_coarse)
d9p$loss <- sprintf("%.2e", d9$loss)
d9p$loss_unbounded_window <- sprintf("%.2e", d9$loss_unbounded_window)
print(d9p, row.names = FALSE)

ok <- vapply(
    steps,
    function(D) all(d9$loss[d9$step == D] <= loss_tol),
    logical(1)
)
chosen <- if (any(ok)) max(steps[ok]) else 0.25
cat(sprintf(
    "\nDECISION: steps with loss <= %.0e on both designs: %s; default step = %s%s\n",
    loss_tol,
    if (any(ok)) paste(steps[ok], collapse = ", ") else "none",
    chosen,
    if (any(ok)) "" else " (fallback)"
))
cat(sprintf("Script runtime: %.1f s\n", proc.time()[["elapsed"]] - t_script))
saveRDS(
    list(results = results, d9 = d9, chosen_step = chosen, loss_tol = loss_tol),
    file.path(here, "d05_grid_step_loss_results.rds")
)
