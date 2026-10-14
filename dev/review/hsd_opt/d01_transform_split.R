# d01_transform_split.R
# Phase D discovery, item 1 of dev/hsd_opt_design_revised.md section 4.3.
# No package code is changed; internals are reached with `multigrain:::`.
#
# (a) Cost split of transform_pvalues_gsd(): the boundary-table build
#     (.gsd_boundary_table(), linear in grid_size) timed separately from the
#     inversion (.gsd_invert(), independent of grid_size), so that the
#     derived split of design section 4.1 (c ~ 1.06 ms per gsBound1() call,
#     I ~ 0) is replaced by a measurement. Plus the interpolation error of
#     coarse tables against G = 1024 (and of every table against an exact
#     root-find, the method of the record's Appendix B), and one full
#     transform at m = 3, K = 2, N = 1e5, G = 1024 to anchor the brief's 3.2 s.
# (b) Well-ordering and invertibility of HSD tables at every node of
#     seq(-10, 2, 0.25) on four schedules plus two corners (design 2.1), the
#     identity h(t) = sfHSD(1, t, g)$spend, exact proportionality in the
#     level, and acceptance of param = 0.
#
# Run from the repository root:
#   Rscript dev/review/hsd_opt/d01_transform_split.R

suppressPackageStartupMessages({
    library(multigrain)
    library(gsDesign)
})

here <- "dev/review/hsd_opt"
if (!dir.exists(here)) {
    here <- "."
}
out_txt <- file.path(here, "d01_transform_split_output.txt")
out_rds <- file.path(here, "d01_transform_split_results.rds")

t_script <- proc.time()
alpha <- 0.025
n_rep <- 5L
n_cores <- 2L

# Output goes to the console and to a buffer; the buffer is written to the
# text file at the end, under a header that carries the sweep verdict.
buffer <- character()
say <- function(...) {
    line <- paste0(...)
    cat(line, "\n", sep = "")
    buffer <<- c(buffer, line)
}
say_table <- function(df, digits = 4) {
    lines <- utils::capture.output(print(df, row.names = FALSE, digits = digits))
    for (ln in lines) say(ln)
}

hsd_fn <- function(g) {
    force(g)
    function(a, t) gsDesign::sfHSD(a, t, g)$spend
}
elapsed <- function(expr) {
    system.time(expr)[["elapsed"]]
}

say("R ", getRversion(), "; gsDesign ", packageVersion("gsDesign"),
    "; multigrain ", packageVersion("multigrain"), "; ", Sys.time())

# ---- (a1) boundary-table time by grid size and number of looks -------------
# Median of n_rep builds at HSD(-2). Table time / G is the per-gsBound1() cost.
say("")
say("== (a1) .gsd_boundary_table() time, HSD(-2), median of ", n_rep, " ==")
schedules_a <- list(`2` = c(0.5, 1), `3` = c(1 / 3, 2 / 3, 1))
grid_sizes <- c(128L, 256L, 512L, 1024L)
sp_m2 <- hsd_fn(-2)

tab_time <- expand.grid(G = grid_sizes, K = c(2L, 3L))
tab_time$table_s <- NA_real_
for (r in seq_len(nrow(tab_time))) {
    t_look <- schedules_a[[as.character(tab_time$K[r])]]
    times <- vapply(seq_len(n_rep), function(i) {
        elapsed(multigrain:::.gsd_boundary_table(
            t_look,
            spending = sp_m2,
            alpha = alpha,
            grid_size = tab_time$G[r]
        ))
    }, numeric(1))
    tab_time$table_s[r] <- stats::median(times)
}
tab_time$ms_per_gsBound1 <- 1000 * tab_time$table_s / tab_time$G
say_table(tab_time)

# Least-squares fit table_s = c * G + b0 per K (b0 is the per-table overhead).
fit_c <- lapply(c(2L, 3L), function(k) {
    d <- tab_time[tab_time$K == k, ]
    stats::coef(stats::lm(table_s ~ G, data = d))
})
names(fit_c) <- c("2", "3")
for (k in names(fit_c)) {
    say(sprintf(
        "K=%s: fit table_s = %.4f ms * G + %.4f s (slope = ms per gsBound1 call)",
        k, 1000 * fit_c[[k]][["G"]], fit_c[[k]][["(Intercept)"]]
    ))
}

# ---- (a2) inversion time by N and K ---------------------------------------
# .gsd_invert() per look on a prebuilt G = 1024 table; per hypothesis = all
# K looks. Median of 11 with bench::mark (sub-millisecond resolution).
# Raw p-values uniform.
say("")
say("== (a2) .gsd_invert() time on G = 1024 tables, median of 11 ==")
tabs_1024 <- lapply(schedules_a, function(tl) {
    multigrain:::.gsd_boundary_table(tl, spending = sp_m2, alpha = alpha,
                                     grid_size = 1024L)
})
inv_time <- expand.grid(N = c(1e4, 1e5), K = c(2L, 3L))
inv_time$per_hyp_ms <- NA_real_
for (r in seq_len(nrow(inv_time))) {
    k <- inv_time$K[r]
    tab <- tabs_1024[[as.character(k)]]
    set.seed(101 + r)
    raw <- matrix(stats::runif(inv_time$N[r] * k), ncol = k)
    bm <- bench::mark(
        for (l in seq_len(k)) {
            multigrain:::.gsd_invert(tab$bounds[, l], grid = tab$grid,
                                     p = raw[, l])
        },
        iterations = 11L, check = FALSE, filter_gc = FALSE
    )
    inv_time$per_hyp_ms[r] <- 1000 * as.numeric(bm$median)
}
inv_time$per_look_ms <- inv_time$per_hyp_ms / inv_time$K
inv_time$m2_ms <- 2 * inv_time$per_hyp_ms
inv_time$m3_ms <- 3 * inv_time$per_hyp_ms
inv_time$per_look_ms_per_1e5_rows <- inv_time$per_look_ms * 1e5 / inv_time$N
say_table(inv_time)

# ---- (a3) one full transform, m = 3, K = 2, N = 1e5, G = 1024 -------------
# One HSD(-2) function recycled to all three hypotheses (three tables built).
say("")
say("== (a3) full transform_pvalues_gsd(), m = 3, K = 2, N = 1e5, G = 1024 ==")
set.seed(20261006)
raw_full <- array(stats::runif(1e5 * 3 * 2), dim = c(1e5, 3L, 2L))
full_s <- elapsed(
    pv_full <- transform_pvalues_gsd(
        raw_full,
        info_frac = c(0.5, 1),
        spending = sp_m2,
        alpha = alpha,
        grid_size = 1024L
    )
)
pred_full <- 3 * tab_time$table_s[tab_time$G == 1024L & tab_time$K == 2L] +
    3 * inv_time$per_hyp_ms[inv_time$N == 1e5 & inv_time$K == 2L] / 1000
say(sprintf("measured %.2f s; predicted 3 tables + 3 inversions = %.2f s; brief 3.2 s",
            full_s, pred_full))
rm(raw_full, pv_full)

# ---- (a4) interpolation error ---------------------------------------------
# 200 log-uniform p-values on [1e-6, 0.025] per look, transformed through the
# K = 2 table at each G. Error of each coarse table against G = 1024, and of
# every table against the exact level (look 1: a = p / h(t1), since the first
# boundary is the first spend; look 2: uniroot on log a). Only p-values inside
# the table (result neither 1 nor the floor) enter.
say("")
say("== (a4) interpolation error, K = 2, t = (0.5, 1), 200 log-uniform p in [1e-6, 0.025] ==")
t2 <- c(0.5, 1)
thetas <- c(-4, -2, 1)
bound_at <- function(a, sp) {
    cs <- sp(a, t2)
    b <- gsDesign::gsBound1(theta = 0, I = t2, a = rep(-20, 2),
                            probhi = diff(c(0, cs)))$b
    stats::pnorm(b, lower.tail = FALSE)
}
interp <- list()
for (g in thetas) {
    sp <- hsd_fn(g)
    tabs <- lapply(grid_sizes, function(G) {
        multigrain:::.gsd_boundary_table(t2, spending = sp, alpha = alpha,
                                         grid_size = G)
    })
    names(tabs) <- grid_sizes
    set.seed(round(1000 + 10 * g))
    p <- exp(stats::runif(200, log(1e-6), log(0.025)))
    h1 <- gsDesign::sfHSD(1, t2, g)$spend[1]
    for (l in 1:2) {
        tr <- lapply(tabs, function(tab) {
            multigrain:::.gsd_invert(tab$bounds[, l], grid = tab$grid, p = p)
        })
        inside <- tr[["1024"]] < 1 & tr[["1024"]] > multigrain:::gsd_grid_min
        exact <- rep(NA_real_, length(p))
        if (l == 1) {
            exact[inside] <- p[inside] / h1
        } else {
            for (j in which(inside)) {
                f <- function(la) log(bound_at(exp(la), sp)[2]) - log(p[j])
                exact[j] <- exp(stats::uniroot(
                    f,
                    interval = log(c(multigrain:::gsd_grid_min, alpha)),
                    tol = 1e-13
                )$root)
            }
        }
        for (G in grid_sizes) {
            gk <- as.character(G)
            vs_1024 <- abs(tr[[gk]][inside] / tr[["1024"]][inside] - 1)
            vs_exact <- abs(tr[[gk]][inside] / exact[inside] - 1)
            interp[[length(interp) + 1L]] <- data.frame(
                theta = g, look = l, G = G, n = sum(inside),
                max_vs_1024 = max(vs_1024), med_vs_1024 = stats::median(vs_1024),
                max_vs_exact = max(vs_exact), med_vs_exact = stats::median(vs_exact)
            )
        }
    }
}
interp <- do.call(rbind, interp)
say_table(interp, digits = 3)
say("Worst case over theta and look, per G:")
interp_g <- stats::aggregate(
    cbind(max_vs_1024, max_vs_exact) ~ G, data = interp, FUN = max
)
say_table(interp_g, digits = 3)

# ---- (b) well-ordering and invertibility sweep ----------------------------
# Every node of seq(-10, 2, 0.25) on four schedules, plus the two corners.
# The builder aborts on a well-ordering failure: the message is recorded.
# Invertibility is checked as .gsd_invert() checks it (distinct positive
# boundaries per column >= 2) and by actually calling .gsd_invert() on each
# column. Run on 2 cores with parallel::mclapply (pure R, deterministic).
say("")
say("== (b) well-ordering and invertibility sweep, G = 1024 ==")
nodes <- seq(-10, 2, by = 0.25)
schedules_b <- list(
    `(0.5,1)` = c(0.5, 1),
    `(1/3,2/3,1)` = c(1 / 3, 2 / 3, 1),
    `(0.25,0.5,0.75,1)` = c(0.25, 0.5, 0.75, 1),
    `(0.6,1)` = c(0.6, 1)
)
cases <- expand.grid(node = nodes, schedule = names(schedules_b),
                     stringsAsFactors = FALSE)
corners <- data.frame(node = c(-10, 2), schedule = c("(0.1,1)", "(0.9,1)"),
                      stringsAsFactors = FALSE)
cases <- rbind(cases, corners)
all_t <- c(schedules_b, list(`(0.1,1)` = c(0.1, 1), `(0.9,1)` = c(0.9, 1)))

sweep_one <- function(r) {
    g <- cases$node[r]
    t_look <- all_t[[cases$schedule[r]]]
    sp <- hsd_fn(g)
    spend_ok <- tryCatch(
        multigrain:::.gsd_check_spending(sp, alpha = alpha, hyp = 1),
        error = function(e) conditionMessage(e)
    )
    tab <- tryCatch(
        multigrain:::.gsd_boundary_table(t_look, spending = sp, alpha = alpha,
                                         grid_size = 1024L, hyp = 1),
        error = function(e) e
    )
    if (inherits(tab, "error")) {
        return(data.frame(
            node = g, schedule = cases$schedule[r],
            spend_check = isTRUE(spend_ok), well_ordered = FALSE,
            message = conditionMessage(tab), min_distinct = NA_integer_,
            invert_ok = FALSE, min_bound = NA_real_, min_pos_bound = NA_real_,
            n_zero = NA_integer_
        ))
    }
    distinct <- apply(tab$bounds, 2, function(col) sum(!duplicated(col) & col > 0))
    inv <- vapply(seq_len(ncol(tab$bounds)), function(l) {
        res <- tryCatch(
            multigrain:::.gsd_invert(tab$bounds[, l], grid = tab$grid,
                                     p = c(1e-12, 1e-6, 1e-3, 0.02)),
            error = function(e) NULL
        )
        !is.null(res)
    }, logical(1))
    data.frame(
        node = g, schedule = cases$schedule[r],
        spend_check = isTRUE(spend_ok), well_ordered = TRUE, message = "",
        min_distinct = min(distinct), invert_ok = all(inv),
        min_bound = min(tab$bounds),
        min_pos_bound = min(tab$bounds[tab$bounds > 0]),
        n_zero = sum(tab$bounds <= 0)
    )
}
t_sweep <- elapsed(
    sweep <- do.call(rbind, parallel::mclapply(
        seq_len(nrow(cases)), sweep_one, mc.cores = n_cores
    ))
)
sweep$pass <- sweep$spend_check & sweep$well_ordered & sweep$invert_ok &
    !is.na(sweep$min_distinct) & sweep$min_distinct >= 2L
say(sprintf("%d tables built in %.1f s on %d cores", nrow(sweep), t_sweep, n_cores))

main <- sweep[sweep$schedule %in% names(schedules_b), ]
pass_mat <- matrix(
    ifelse(main$pass, "pass", "FAIL"),
    nrow = length(nodes),
    dimnames = list(format(nodes), names(schedules_b))
)
mindist_mat <- matrix(main$min_distinct, nrow = length(nodes),
                      dimnames = list(format(nodes), names(schedules_b)))
say("Pass/fail (node x schedule) and min over columns of distinct positive boundaries:")
both <- matrix(paste0(pass_mat, " ", mindist_mat), nrow = length(nodes),
               dimnames = dimnames(pass_mat))
for (ln in utils::capture.output(print(noquote(both)))) say(ln)
say("Corners:")
say_table(sweep[!sweep$schedule %in% names(schedules_b),
                c("node", "schedule", "pass", "min_distinct", "min_bound",
                  "min_pos_bound")], digits = 3)
say("Per schedule: smallest distinct count, smallest boundary, smallest positive boundary:")
say_table(stats::aggregate(
    cbind(min_distinct, min_bound, min_pos_bound) ~ schedule,
    data = sweep, FUN = min
), digits = 3)
say(sprintf("Boundaries <= 0 anywhere: %d", sum(sweep$n_zero, na.rm = TRUE)))
fails <- sweep[!sweep$pass, ]
if (nrow(fails) > 0) {
    say("FAILING cases:")
    say_table(fails[, c("node", "schedule", "spend_check", "well_ordered",
                        "invert_ok", "min_distinct", "message")])
}

# ---- (b2) sfHSD identities ------------------------------------------------
# param = 0 accepted (linear); h(t) = sfHSD(1, t, g)$spend against the closed
# form; exact proportionality in the level: spend(1e-14)/1e-14 vs
# spend(0.025)/0.025, worst relative difference over all nodes and schedules.
say("")
say("== (b2) sfHSD identities ==")
p0 <- tryCatch(gsDesign::sfHSD(0.025, c(0.5, 1), 0)$spend,
               error = function(e) conditionMessage(e))
say("sfHSD(0.025, c(0.5, 1), 0)$spend = ", toString(format(p0, digits = 15)))
ident <- do.call(rbind, lapply(nodes, function(g) {
    do.call(rbind, lapply(names(all_t), function(s) {
        t_look <- all_t[[s]]
        h <- gsDesign::sfHSD(1, t_look, g)$spend
        closed <- if (g == 0) t_look else (1 - exp(-g * t_look)) / (1 - exp(-g))
        lo <- gsDesign::sfHSD(1e-14, t_look, g)$spend / 1e-14
        hi <- gsDesign::sfHSD(0.025, t_look, g)$spend / 0.025
        data.frame(
            node = g, schedule = s,
            h_vs_closed = max(abs(h / closed - 1)),
            prop_rel_diff = max(abs(lo / hi - 1)),
            min_increment_1e14 = min(diff(c(0, lo * 1e-14)))
        )
    }))
}))
say(sprintf("max |sfHSD(1,t,g) / closed form - 1| = %.3g", max(ident$h_vs_closed)))
say(sprintf("max |[spend(1e-14)/1e-14] / [spend(0.025)/0.025] - 1| = %.3g",
            max(ident$prop_rel_diff)))
say(sprintf("smallest spend increment at level 1e-14 = %.3g (node %g, schedule %s)",
            min(ident$min_increment_1e14),
            ident$node[which.min(ident$min_increment_1e14)],
            ident$schedule[which.min(ident$min_increment_1e14)]))

# ---- verdict, header, save ------------------------------------------------
n_main_fail <- sum(!main$pass)
n_corner_fail <- sum(!sweep$pass) - n_main_fail
verdict <- if (sum(!sweep$pass) == 0) {
    "VERDICT: no node in [-10, 2] fails (well-ordering and invertibility, all 4 schedules and both corners)."
} else {
    sprintf(paste0(
        "!!! PAUSE CONDITION: %d node x schedule case(s) FAIL (%d on the four ",
        "schedules, %d at the corners); see 'FAILING cases' below. !!!"
    ), sum(!sweep$pass), n_main_fail, n_corner_fail)
}
runtime <- (proc.time() - t_script)[["elapsed"]]
say("")
say(verdict)
say(sprintf("Script runtime: %.1f s", runtime))

header <- c(
    "# d01_transform_split.R output",
    verdict,
    ""
)
writeLines(c(header, buffer), out_txt)
saveRDS(
    list(
        table_time = tab_time, fit_c = fit_c, invert_time = inv_time,
        full_transform_s = full_s, full_transform_pred_s = pred_full,
        interp = interp, interp_by_G = interp_g, sweep = sweep,
        pass_matrix = pass_mat, min_distinct_matrix = mindist_mat,
        param0_spend = p0, identities = ident, verdict = verdict,
        runtime_s = runtime, session = utils::sessionInfo()
    ),
    out_rds
)
