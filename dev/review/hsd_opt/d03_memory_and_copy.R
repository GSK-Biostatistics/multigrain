# d03_memory_and_copy.R
# Phase D discovery, item 3 of dev/hsd_opt_design_revised.md section 4.3.
# No package code is changed; internals are reached with `multigrain:::`.
#
# The bank closure of P2 keeps one N_GA x mK working matrix and, per GA
# evaluation, overwrites the K columns of each tagged hypothesis with the
# bank slice of the chosen node before calling the kernel. This script
# measures, at N_GA = 5e4, m = 3, K = 2 (README section-2 design, sfLDOF):
#   1. whether the kernel call or a column assignment copies the matrix
#      (tracemem), for the assignment patterns a closure can use;
#   2. the time of one K-column assignment against one kernel call, and their
#      ratio against the 15 % threshold of the plan (a C++ column copier is
#      added only above it);
#   3. object.size() of the matrix and the projected bank bytes
#      8 * K * N_GA * sum(M_g) for M = 21, 41, 97 at m_act = 1, 2, 4.
#
# Run from the repository root:
#   Rscript dev/review/hsd_opt/d03_memory_and_copy.R

suppressPackageStartupMessages({
    library(multigrain)
    library(gsDesign)
})

here <- "dev/review/hsd_opt"
if (!dir.exists(here)) {
    here <- "."
}
out_txt <- file.path(here, "d03_memory_and_copy_output.txt")
out_rds <- file.path(here, "d03_memory_and_copy_results.rds")

t_script <- proc.time()
buffer <- character()
say <- function(...) {
    line <- paste0(...)
    cat(line, "\n", sep = "")
    buffer <<- c(buffer, line)
}
say_table <- function(df, digits = 4) {
    lines <- utils::capture.output(print(
        df,
        row.names = FALSE,
        digits = digits
    ))
    for (ln in lines) {
        say(ln)
    }
}

say(
    "R ",
    getRversion(),
    "; multigrain ",
    packageVersion("multigrain"),
    "; profmem capability ",
    capabilities("profmem"),
    "; ",
    Sys.time()
)

# ---- design: README section 2 (m = 3, K = 2, sfLDOF), subset to N_GA -------
alpha <- 0.025
power_nom <- c(H1 = 0.90, H2 = 0.80, H3 = 0.70)
corr <- matrix(0.5, 3, 3)
diag(corr) <- 1
info_frac <- c(0.5, 1)
n_ga <- 5e4
n_threads <- 2L
K <- 2L
m <- 3L

set.seed(20261006)
raw <- simulate_pvalues_gsd(
    power_nominal = power_nom,
    alpha = alpha,
    corr_matrix = corr,
    info_frac = info_frac,
    nsim = 1e5
)
pv <- transform_pvalues_gsd(raw, spending = sfLDOF, alpha = alpha)
# the GA subsample, drawn as .sample_pvals_gsd() draws it; the same rows are
# used for the bank slice below
sample_rows <- function(obj) {
    set.seed(1)
    multigrain:::.sample_pvals_gsd(obj, n_ga)
}
pv_ga <- sample_rows(pv)

x <- multigrain:::.gsd_kernel_matrix(pv_ga)
w <- rep(1 / 3, 3)
G <- matrix(0.5, 3, 3)
diag(G) <- 0

# Replacement for hypothesis 2 (columns 2 and 5 = H2 at looks 1 and 2): a
# bank slice from another spending function (HSD(-2)) on the same rows.
pv_alt <- transform_pvalues_gsd(
    raw,
    spending = function(a, t) gsDesign::sfHSD(a, t, -2)$spend,
    alpha = alpha
)
pv_alt <- sample_rows(pv_alt)
bank_slice <- multigrain:::.gsd_kernel_matrix(pv_alt)[, c(2L, 5L)]
cols <- c(2L, 5L)
rm(raw, pv, pv_alt)
invisible(gc())

say(sprintf("Working matrix: %d x %d, storage %s", nrow(x), ncol(x), typeof(x)))

kernel <- function(mat) {
    multigrain:::graph_shortcut_gsd_parallel(
        mat,
        alpha,
        w,
        G,
        K,
        n_threads,
        -1L,
        NULL
    )
}

# ---- 1. tracemem ----------------------------------------------------------
# Each probe tracemem()s a fresh matrix and captures what R prints when the
# object is duplicated. "copies" is the number of tracemem[...] lines.
say("")
say("== 1. tracemem: does the kernel call or a column assignment copy? ==")
count_copies <- function(lines) sum(grepl("^tracemem\\[", lines))
probes <- list()

# (a) top level: kernel call, then a column assignment
x_a <- x + 0
invisible(tracemem(x_a))
log_kernel <- utils::capture.output(res_a <- kernel(x_a))
log_assign <- utils::capture.output(x_a[, cols] <- bank_slice)
untracemem(x_a)
probes$kernel_call <- log_kernel
probes$assign_top_level_after_kernel <- log_assign

# (b) closure that owns the matrix and assigns with <<- (the P2 pattern)
make_closure_super <- function(mat) {
    force(mat)
    list(
        trace = function() tracemem(mat),
        untrace = function() untracemem(mat),
        eval = function(slice) {
            mat[, cols] <<- slice
            kernel(mat)
        }
    )
}
# Evaluations logged one by one. When `mat` arrives as an argument its
# binding is a promise; the first <<- rebinds `mat` to the value while the
# promise still references it, so the second evaluation copies once. That is
# a one-off, not a per-call cost; (b2) shows it vanishes when the matrix is
# built inside the factory.
cl_b <- make_closure_super(x + 0)
invisible(cl_b$trace())
for (i in 1:3) {
    probes[[paste0("closure_superassign_eval_", i)]] <-
        utils::capture.output(invisible(cl_b$eval(bank_slice)))
}
cl_b$untrace()

# (b2) as (b), but the matrix is built inside the factory, as the P2 closure
# builds it from the spec (no argument promise holds it)
make_closure_built <- function() {
    mat <- multigrain:::.gsd_kernel_matrix(pv_ga)
    list(
        trace = function() tracemem(mat),
        untrace = function() untracemem(mat),
        eval = function(slice) {
            mat[, cols] <<- slice
            kernel(mat)
        }
    )
}
cl_b2 <- make_closure_built()
invisible(cl_b2$trace())
for (i in 1:3) {
    probes[[paste0("closure_built_inside_eval_", i)]] <-
        utils::capture.output(invisible(cl_b2$eval(bank_slice)))
}
cl_b2$untrace()

# (c) closure that assigns with <- (local copy of the enclosing binding)
make_closure_local <- function(mat) {
    force(mat)
    list(
        trace = function() tracemem(mat),
        untrace = function() untracemem(mat),
        eval = function(slice) {
            mat[, cols] <- slice
            kernel(mat)
        }
    )
}
cl_c <- make_closure_local(x + 0)
invisible(cl_c$trace())
log_c <- utils::capture.output(
    for (i in 1:3) {
        cl_c$eval(bank_slice)
    }
)
cl_c$untrace()
probes$closure_local_assign_3_evals <- log_c

# (d) matrix held in an environment, assigned through env$mat[, cols] <-
env_d <- new.env()
env_d$mat <- x + 0
invisible(tracemem(env_d$mat))
log_d <- utils::capture.output(
    for (i in 1:3) {
        env_d$mat[, cols] <- bank_slice
        kernel(env_d$mat)
    }
)
untracemem(env_d$mat)
probes$env_assign_3_evals <- log_d

# (e) a second binding to the same matrix (e.g. the caller keeps `x`):
# the first assignment must copy
x_e <- x + 0
x_e_alias <- x_e
invisible(tracemem(x_e))
log_e <- utils::capture.output({
    x_e[, cols] <- bank_slice
    x_e[, cols] <- bank_slice
})
untracemem(x_e)
probes$aliased_2_assigns <- log_e
rm(x_e_alias)

copy_tab <- data.frame(
    probe = names(probes),
    copies = vapply(probes, count_copies, integer(1))
)
say_table(copy_tab)
for (nm in names(probes)) {
    if (length(probes[[nm]]) > 0) {
        say(nm, ": ", paste(probes[[nm]], collapse = " | "))
    }
}
tracemem_verdict <- if (
    copy_tab$copies[copy_tab$probe == "kernel_call"] == 0 &&
        copy_tab$copies[copy_tab$probe == "assign_top_level_after_kernel"] ==
            0 &&
        sum(copy_tab$copies[grepl(
            "^closure_superassign_eval_",
            copy_tab$probe
        )]) <=
            1 &&
        all(
            copy_tab$copies[grepl("^closure_built_inside", copy_tab$probe)] == 0
        )
) {
    paste0(
        "TRACEMEM: no per-evaluation copy - the kernel call never copies; a K-column ",
        "assignment with <<- copies at most once in total when the matrix arrived as a ",
        "factory argument (",
        sum(copy_tab$copies[grepl(
            "^closure_superassign_eval_",
            copy_tab$probe
        )]),
        " copy over 3 evaluations), never when it is built inside the closure ",
        "factory, and on every evaluation with a plain <- inside the closure."
    )
} else {
    "TRACEMEM: a PER-EVALUATION COPY was triggered by the kernel call or the column assignment; see the probes above."
}
say(tracemem_verdict)

# ---- 2. assignment vs kernel time ----------------------------------------
# 200 K-column assignments vs 200 kernel calls, each timed individually with
# bench::hires_time(); medians reported. Interleaved so that drift hits both.
say("")
say(
    "== 2. K-column assignment vs kernel call, N_GA = 5e4, m = 3, K = 2, 2 threads =="
)
n_iter <- 200L
x_t <- x + 0
slices <- list(bank_slice, x[, cols])
t_assign <- numeric(n_iter)
t_kernel <- numeric(n_iter)
t_cl <- numeric(n_iter)
cl_t <- make_closure_super(x + 0)
invisible(kernel(x_t))
for (i in seq_len(n_iter)) {
    s <- slices[[1L + (i %% 2L)]]
    t0 <- bench::hires_time()
    x_t[, cols] <- s
    t1 <- bench::hires_time()
    invisible(kernel(x_t))
    t2 <- bench::hires_time()
    invisible(cl_t$eval(s))
    t3 <- bench::hires_time()
    t_assign[i] <- t1 - t0
    t_kernel[i] <- t2 - t1
    t_cl[i] <- t3 - t2
}
timing <- data.frame(
    operation = c(
        "assign K = 2 columns (top level)",
        "kernel call",
        "closure eval (<<- assign + kernel)"
    ),
    median_ms = 1000 *
        c(
            stats::median(t_assign),
            stats::median(t_kernel),
            stats::median(t_cl)
        ),
    p90_ms = 1000 *
        c(
            stats::quantile(t_assign, 0.9),
            stats::quantile(t_kernel, 0.9),
            stats::quantile(t_cl, 0.9)
        )
)
say_table(timing)
ratio <- stats::median(t_assign) / stats::median(t_kernel)
closure_overhead <- (stats::median(t_cl) - stats::median(t_kernel)) /
    stats::median(t_kernel)
say(sprintf(
    "copy / kernel = %.1f %% (threshold 15 %%): %s",
    100 * ratio,
    if (ratio > 0.15) {
        "ABOVE - C++ column copier warranted"
    } else {
        "below - no C++ column copier needed"
    }
))
say(sprintf(
    "closure eval overhead over a bare kernel call = %.1f %%",
    100 * closure_overhead
))
# Per tagged hypothesis count: the bank closure assigns K columns per active
# tagged hypothesis, so the cost scales with m_act.
for (m_act in c(1, 2, 4)) {
    say(sprintf(
        "  projected at m_act = %d: %.1f %% of a kernel call",
        m_act,
        100 * m_act * ratio
    ))
}

# ---- 3. memory ------------------------------------------------------------
say("")
say("== 3. memory ==")
say(sprintf(
    "object.size(working matrix) = %s (%.0f bytes; 8 * N_GA * mK = %.0f)",
    format(object.size(x), units = "MB"),
    as.numeric(object.size(x)),
    8 * n_ga * m * K
))
bank <- expand.grid(M = c(21L, 41L, 97L), m_act = c(1L, 2L, 4L), K = c(2L, 3L))
bank$MB <- 8 * bank$K * n_ga * bank$M * bank$m_act / 1e6
bank_wide <- stats::reshape(
    bank,
    idvar = c("K", "M"),
    timevar = "m_act",
    direction = "wide"
)
names(bank_wide) <- sub("MB.", "MB_m_act_", names(bank_wide), fixed = TRUE)
say(
    "Projected bank bytes 8 * K * N_GA * M * m_act, N_GA = 5e4, in MB (1e6 bytes):"
)
say_table(bank_wide[order(bank_wide$K, bank_wide$M), ], digits = 5)

runtime <- (proc.time() - t_script)[["elapsed"]]
say("")
say(sprintf("Script runtime: %.1f s", runtime))

writeLines(
    c("# d03_memory_and_copy.R output", tracemem_verdict, "", buffer),
    out_txt
)
saveRDS(
    list(
        copies = copy_tab,
        tracemem_logs = probes,
        tracemem_verdict = tracemem_verdict,
        timing = timing,
        t_assign = t_assign,
        t_kernel = t_kernel,
        t_closure = t_cl,
        ratio = ratio,
        closure_overhead = closure_overhead,
        object_size = as.numeric(object.size(x)),
        bank = bank,
        runtime_s = runtime,
        session = utils::sessionInfo()
    ),
    out_rds
)
