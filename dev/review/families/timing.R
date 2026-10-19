# F1 gate for families of hypotheses: cost of the within-family adjustment
# relative to one kernel call. Decides whether "raw p-values everywhere"
# (re-adjusting before each kernel call outside the objective) stays.
#
# Usage: Rscript dev/review/families/timing.R [<library with multigrain>]
#
# Design of the families record, section 1.4: m = 20, M = 6, Hommel in every
# family, the drawn graph. Serial timings, median of 7 repeats.
# Output: dev/review/families/results/timing.txt
#
# Local timings: keep the laptop on AC power and awake (see the local
# benchmark notes); the ratio is the useful number, not the absolute times.

args <- commandArgs(TRUE)
if (length(args)) {
    .libPaths(c(args[1], .libPaths()))
}
library(multigrain)

repo <- normalizePath(".")
out_dir <- file.path(repo, "dev", "review", "families", "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

env <- new.env(parent = asNamespace("multigrain"))
sys.source(file.path(repo, "tests", "testthat", "helper-families.R"), env)

sink(file.path(out_dir, "timing.txt"), split = TRUE)
on.exit(sink(), add = TRUE)

cat("multigrain", format(utils::packageVersion("multigrain")), "\n")
cat(format(Sys.time()), "\n")
cat(utils::sessionInfo()$running, "\n\n")

with(env, {
    members <- list(1, 2:5, 6:8, 9:12, 13:16, 17:20)
    fam_index <- family_index0(members)
    w <- c(1, 0, 0, 0, 0, 0)
    G <- matrix(0, 6, 6)
    G[1, 2:3] <- c(0.25, 0.75)
    G[2, 3] <- 1
    G[3, 4:6] <- 1 / 3
    G[4, c(2, 5, 6)] <- c(0.001, 0.4995, 0.4995)
    G[5, c(4, 6)] <- 0.5
    G[6, c(4, 5)] <- 0.5

    # Median over `reps` of the time per call, each rep timing enough calls
    # (about 2e6 rows in total) to be well above the 10 ms Windows timer.
    time_median <- function(expr_fun, n, reps = 7) {
        calls <- max(1, round(2e6 / n))
        per_call <- replicate(reps, {
            system.time(
                for (i in seq_len(calls)) {
                    expr_fun()
                }
            )[["elapsed"]] /
                calls
        })
        stats::median(per_call)
    }

    set.seed(1)
    cat(sprintf(
        "%-8s %12s %12s %8s\n",
        "N",
        "adjust (s)",
        "kernel (s)",
        "ratio"
    ))
    for (n in c(5e4, 2e5, 1e6)) {
        pvals <- matrix(
            stats::pnorm(stats::rnorm(n * 20, mean = 2.5), lower.tail = FALSE),
            ncol = 20
        )
        padj <- adjust_families_r(pvals, members, "hommel")
        t_adjust <- time_median(
            \() adjust_families_r(pvals, members, "hommel"),
            n = n
        )
        t_kernel <- time_median(
            \() graph_shortcut_families(padj, fam_index, 0.025, w, G),
            n = n
        )
        cat(sprintf(
            "%-8s %12.4f %12.4f %8.2f\n",
            format(n, scientific = TRUE),
            t_adjust,
            t_kernel,
            t_adjust / t_kernel
        ))
    }
})
