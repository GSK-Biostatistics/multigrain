# gsd-joint-spending-vignette.R
# Precomputes every slow result shown in the article
# vignettes/articles/gsd-joint-spending.Rmd and stores them in
# vignettes/articles/data/gsd-joint-spending.rds. The article shows the code
# of these steps with `eval = FALSE` and prints the stored objects, so it
# builds in seconds. Run from the package root after installing the package:
#
#   Rscript data-raw/gsd-joint-spending-vignette.R
#
# Main design: the running example of trial_success.qmd, "Optimising the
# spending parameter" (seed 20261007, nsim 1e5), so its numbers can be checked
# against dev/review/hsd_opt/article_example_output.txt. Second design: the
# manuscript's Example 5 (PFS and OS) with the spending of OS chosen jointly.

library(multigrain)
options(multigrain_verbosity = "silent", width = 80)

out_file <- "vignettes/articles/data/gsd-joint-spending.rds"
n_threads <- 2L

timed <- function(label, expr) {
    t0 <- proc.time()[["elapsed"]]
    value <- force(expr)
    cat(sprintf("[time] %s: %.1f s\n", label, proc.time()[["elapsed"]] - t0))
    value
}

t_start <- proc.time()[["elapsed"]]

# ---- main design: H1 primary, H2 key secondary, stop once H1 is rejected ----

set.seed(20261007)
raw <- simulate_pvalues_gsd(
    power_nominal = c(H1 = 0.90, H2 = 0.80),
    corr_matrix = matrix(c(1, 0.5, 0.5, 1), nrow = 2),
    info_frac = c(0.5, 1),
    nsim = 1e5
)

gain <- trial_success_gsd(0.6 * d(t1) + 0.4 * d(t2), d = c(1, 0.8))
rule <- stopping_rule_gsd(r1)

pv_ldof <- timed(
    "main: transform, sfLDOF on both",
    transform_pvalues_gsd(raw, spending = gsDesign::sfLDOF)
)
main_base <- timed(
    "main: graph_optimise_gsd(), fixed spending",
    graph_optimise_gsd(
        pv_ldof,
        graph_constraint_free(2),
        gain,
        stopping_rule = rule,
        num_threads = n_threads
    )
)

spec <- timed(
    "main: transform, H2 tagged (spec)",
    transform_pvalues_gsd(
        raw,
        spending = list(gsDesign::sfLDOF, sf_hsd_opt())
    )
)
main_spec_print <- utils::capture.output(print(spec))

main_joint <- timed(
    "main: graph_optimise_gsd(), joint",
    graph_optimise_gsd(
        spec,
        graph_constraint_free(2),
        gain,
        stopping_rule = rule,
        num_threads = n_threads
    )
)

# the chosen design on fresh draws
set.seed(20261008)
raw_new <- simulate_pvalues_gsd(
    power_nominal = c(H1 = 0.90, H2 = 0.80),
    corr_matrix = matrix(c(1, 0.5, 0.5, 1), nrow = 2),
    info_frac = c(0.5, 1),
    nsim = 1e5
)
spec_new <- timed(
    "main: fresh spec",
    transform_pvalues_gsd(
        raw_new,
        spending = list(gsDesign::sfLDOF, sf_hsd_opt())
    )
)
pv_new <- resolve_pvalues_gsd(spec_new, main_joint$spending_param)
check <- calc_power_pvals_gsd(
    pv_new,
    hyp_weight = main_joint$hyp_weight,
    trans_matrix = main_joint$trans_matrix,
    stopping_rule = rule,
    custom_power = list(gain = gain)
)
main_check <- list(
    spending = pv_new$spending,
    gain = check$gain,
    local_power = check$local_power,
    stop_distribution = check$stop_distribution
)

rm(raw, raw_new, spec, spec_new, pv_ldof, pv_new)
invisible(gc())

# ---- Example 5: PFS complete at the interim, OS at 70% and 100% --------------

info_frac_ex5 <- rbind(PFS = c(1, 1), OS = c(0.7, 1))

set.seed(20261014)
raw5 <- simulate_pvalues_gsd(
    power_nominal = c(PFS = 0.98, OS = 0.93),
    corr_matrix = matrix(c(1, 0.5, 0.5, 1), nrow = 2),
    info_frac = info_frac_ex5,
    nsim = 1e5
)

gain5 <- trial_success_gsd(
    0.4 * d_pfs(t1) + 0.6 * d_os(t2),
    d_pfs = c(1, 0.50),
    d_os = c(1, 0.75)
)

pv5_ldof <- timed(
    "ex5: transform, sfLDOF",
    transform_pvalues_gsd(raw5, spending = gsDesign::sfLDOF)
)
ex5_base <- timed(
    "ex5: graph_optimise_gsd(), fixed spending",
    graph_optimise_gsd(
        pv5_ldof,
        graph_constraint_free(2),
        gain5,
        num_threads = n_threads
    )
)

ex5_joint_run <- function(bounds) {
    spec5 <- timed(
        sprintf("ex5: spec, bounds [%g, %g]", bounds[1L], bounds[2L]),
        transform_pvalues_gsd(
            raw5,
            spending = list(gsDesign::sfLDOF, sf_hsd_opt(bounds = bounds))
        )
    )
    res <- timed(
        "ex5: graph_optimise_gsd(), joint",
        graph_optimise_gsd(
            spec5,
            graph_constraint_free(2),
            gain5,
            num_threads = n_threads
        )
    )
    list(spec_print = utils::capture.output(print(spec5)), result = res)
}

ex5_bounds <- c(-4, 1)
ex5 <- ex5_joint_run(ex5_bounds)
if (ex5$result$spending_param[[1L]] >= ex5_bounds[2L]) {
    cat("[ex5] optimum on the upper bound 1; widening to [-4, 2]\n")
    ex5_bounds <- c(-4, 2)
    ex5 <- ex5_joint_run(ex5_bounds)
}

# ---- store -------------------------------------------------------------------

out <- list(
    main_base = main_base,
    main_spec_print = main_spec_print,
    main_joint = main_joint,
    main_check = main_check,
    ex5_base = ex5_base,
    ex5_spec_print = ex5$spec_print,
    ex5_joint = ex5$result,
    ex5_bounds = ex5_bounds,
    versions = c(
        R = as.character(getRversion()),
        multigrain = as.character(utils::packageVersion("multigrain")),
        gsDesign = as.character(utils::packageVersion("gsDesign")),
        GA = as.character(utils::packageVersion("GA"))
    )
)
saveRDS(out, out_file)

# keep the stored file small: the genetic algorithm's object is not shown
if (file.size(out_file) > 2e6) {
    for (nm in c("main_base", "main_joint", "ex5_base", "ex5_joint")) {
        out[[nm]]$global_output <- NULL
    }
    saveRDS(out, out_file)
}
cat(sprintf("[size] %s: %.0f kB\n", out_file, file.size(out_file) / 1024))

# ---- report ------------------------------------------------------------------

report <- function(label, res) {
    cat(sprintf(
        "[%s] gain %.4f, w1 %.3f, spending %s\n",
        label,
        res$power$trial_success,
        res$hyp_weight[[1L]],
        paste(res$spending %||% "fixed", collapse = " | ")
    ))
}
report("main base", main_base)
report("main joint", main_joint)
cat(sprintf("[main fresh] gain %.4f\n", main_check$gain))
report("ex5 base", ex5_base)
report("ex5 joint", ex5$result)
print(main_joint)
print(ex5$result)

cat(sprintf("\n[time] total: %.1f s\n", proc.time()[["elapsed"]] - t_start))
