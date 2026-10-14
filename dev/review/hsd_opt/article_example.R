# article_example.R
# The running example of the "Optimising the spending parameter" subsection of
# vignettes/articles/trial_success.qmd (build plan, phase P5; brief
# dev/review/hsd_opt/briefs/P5_documentation.md, deliverable 3). The article
# chunks are `#| eval: false` because quarto is not installed here; every
# `#>` line pasted into them comes from one run of this script, whose output
# is kept beside it:
#
#   Rscript dev/review/hsd_opt/article_example.R
#
# Design: README section 7.5 of dev/review/gsd_user_testing/ under model M1
# (two hypotheses, powers 0.90 and 0.80, correlation 0.5, an interim at half
# information, gain 0.6 * d(t1) + 0.4 * d(t2) with d = c(1, 0.8), stop once
# H1 is rejected), sfLDOF on H1 and sf_hsd_opt() on H2, nsim 1e5, seed
# 20261007 (the draws of Phase D and of gate_p3_benchmarks.R), 2 threads.
# Then, on the same draws, the graph optimised with H2's spending fixed at
# sfLDOF and at sfLDPocock (the two choices of README 7.4 and 7.5), and the
# chosen design re-scored on fresh draws (seed 20261008).

Sys.setenv(NOT_CRAN = "true")
suppressMessages(devtools::load_all(quiet = TRUE))
options(multigrain_verbosity = "silent", width = 80)

if (!exists("%||%", mode = "function")) {
    `%||%` <- function(x, y) if (is.null(x)) y else x
}

out_file <- "dev/review/hsd_opt/article_example_output.txt"
con <- file(out_file, open = "wt")
sink(con, split = TRUE)

timed <- function(label, expr) {
    t0 <- proc.time()[["elapsed"]]
    value <- force(expr)
    cat(sprintf("[time] %s: %.1f s\n", label, proc.time()[["elapsed"]] - t0))
    value
}
show <- function(header, expr) {
    cat("\n## ", header, "\n", sep = "")
    print(expr)
}

t_start <- proc.time()[["elapsed"]]
cat(sprintf(
    "article_example: R %s; multigrain %s; gsDesign %s; GA %s; num_threads 2\n",
    getRversion(),
    utils::packageVersion("multigrain"),
    utils::packageVersion("gsDesign"),
    utils::packageVersion("GA")
))

# ---- chunk 1: simulate and build the spec ------------------------------------

set.seed(20261007)
raw <- simulate_pvalues_gsd(
    power_nominal = c(H1 = 0.90, H2 = 0.80),
    corr_matrix = matrix(c(1, 0.5, 0.5, 1), nrow = 2),
    info_frac = c(0.5, 1),
    nsim = 1e5
)

spec <- timed("transform_pvalues_gsd() spec", transform_pvalues_gsd(
    raw,
    spending = list(gsDesign::sfLDOF, sf_hsd_opt())
))
show("spec", spec)

# ---- chunk 2: optimise graph and parameter jointly ---------------------------

gain <- trial_success_gsd(
    0.6 * d(t1) + 0.4 * d(t2),
    d = c(1, 0.8),
    verbose = "silent"
)
rule <- stopping_rule_gsd(r1)

result <- timed("graph_optimise_gsd() on the spec", graph_optimise_gsd(
    spec,
    graph_constraint_free(2),
    gain,
    stopping_rule = rule,
    num_threads = 2
))
show("print(result)", result)
cat("\n## summary(result)\n")
summary(result)

cat("\n## fields\n")
show("result$spending_param", result$spending_param)
show("result$hyp_spending_param", result$hyp_spending_param)
show("result$spending", result$spending)
show("round(result$hyp_weight, 3)", round(result$hyp_weight, 3))
show("round(result$power$trial_success, 4)",
     round(result$power$trial_success, 4))
show("result$solution$opt_source", result$solution$opt_source)
show("result$profile_converged", result$profile_converged)
cat("\n## result$spending_profile (digits 6)\n")
print(result$spending_profile, digits = 6, row.names = FALSE)
cat("\n## profile graphs: w1 by node\n")
w1_by_node <- vapply(
    attr(result$spending_profile, "graphs"),
    function(g) g$hyp_weight[[1L]],
    double(1L)
)
print(data.frame(
    value = result$spending_profile$value,
    w1 = round(w1_by_node, 3)
), row.names = FALSE)
show("result$spending_search", result$spending_search)

# ---- chunk 3: the same draws with fixed spending on H2 -----------------------

fixed_run <- function(label, h2_spending) {
    pv <- transform_pvalues_gsd(
        raw,
        spending = list(gsDesign::sfLDOF, h2_spending)
    )
    res <- timed(paste("graph_optimise_gsd(), H2", label), graph_optimise_gsd(
        pv,
        graph_constraint_free(2),
        gain,
        stopping_rule = rule,
        num_threads = 2
    ))
    cat(sprintf(
        "[fixed] H2 %s: gain %.6f, w1 %.3f\n",
        label,
        res$power$trial_success,
        res$hyp_weight[[1L]]
    ))
    res
}
fixed_ldof <- fixed_run("sfLDOF", gsDesign::sfLDOF)
fixed_pocock <- fixed_run("sfLDPocock", gsDesign::sfLDPocock)

cat("\n## fixed spending on H2 vs the joint optimum (same draws)\n")
comparison <- data.frame(
    H2_spending = c("sfLDOF", "sfLDPocock", result$spending[[2L]]),
    w1 = round(c(
        fixed_ldof$hyp_weight[[1L]],
        fixed_pocock$hyp_weight[[1L]],
        result$hyp_weight[[1L]]
    ), 3),
    gain = round(c(
        fixed_ldof$power$trial_success,
        fixed_pocock$power$trial_success,
        result$power$trial_success
    ), 4)
)
print(comparison, row.names = FALSE)

# ---- chunk 4: the chosen design on fresh draws -------------------------------

set.seed(20261008)
raw_new <- simulate_pvalues_gsd(
    power_nominal = c(H1 = 0.90, H2 = 0.80),
    corr_matrix = matrix(c(1, 0.5, 0.5, 1), nrow = 2),
    info_frac = c(0.5, 1),
    nsim = 1e5
)
spec_new <- timed("transform_pvalues_gsd() fresh spec", transform_pvalues_gsd(
    raw_new,
    spending = list(gsDesign::sfLDOF, sf_hsd_opt())
))
pv_new <- resolve_pvalues_gsd(spec_new, result$spending_param)
show("pv_new$spending", pv_new$spending)
check <- calc_power_pvals_gsd(
    pv_new,
    hyp_weight = result$hyp_weight,
    trans_matrix = result$trans_matrix,
    stopping_rule = rule,
    custom_power = list(gain = gain)
)
show("round(check$gain, 4)", round(check$gain, 4))
show("round(check$local_power, 3)", round(check$local_power, 3))
show("round(check$stop_distribution, 3)", round(check$stop_distribution, 3))

cat(sprintf(
    "\n[time] total: %.1f s\n",
    proc.time()[["elapsed"]] - t_start
))

sink()
close(con)
