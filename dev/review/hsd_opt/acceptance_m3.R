# dev/review/hsd_opt/acceptance_m3.R -- the end-to-end acceptance run of the
# spending-parameter build (build plan section 6 item 5): the README section 2
# design at nsim 1e4 with one shared parameter on the two secondary
# hypotheses, a stopping rule, two threads and the default control.
#
# Run from the repository root:
#     Rscript dev/review/hsd_opt/acceptance_m3.R

suppressPackageStartupMessages(library(multigrain))

set.seed(1)
raw <- simulate_pvalues_gsd(
    power_nominal = c(0.9, 0.85, 0.8),
    corr_matrix = diag(3),
    info_frac = c(0.5, 1),
    nsim = 1e4
)

t0 <- Sys.time()
spec <- transform_pvalues_gsd(
    raw,
    spending = list(
        gsDesign::sfLDOF,
        sf_hsd_opt(id = "sec"),
        sf_hsd_opt(id = "sec")
    )
)
cat(
    "spec built in",
    round(as.numeric(Sys.time() - t0, units = "secs"), 1),
    "s\n\n"
)
print(spec)

gain <- trial_success_gsd(r1 + r2 + r3, verbose = "silent")
rule <- stopping_rule_gsd(r1)

t1 <- Sys.time()
set.seed(7)
res <- graph_optimise_gsd(
    spec,
    graph_constraint_free(3),
    gain,
    stopping_rule = rule,
    num_threads = 2L,
    verbose = "silent"
)
wall <- as.numeric(Sys.time() - t1, units = "secs")
cat("\noptimiser:", round(wall, 1), "s\n\n")

print(res)
cat("\n")
summary(res)

# the final evaluation is the resolved object's
pv_star <- resolve_pvalues_gsd(spec, res$spending_param)
check <- calc_power_pvals_gsd(
    pv_star,
    res$hyp_weight,
    res$trans_matrix,
    custom_power = list(trial_success = gain),
    stopping_rule = rule
)
on_lattice <- all(vapply(
    seq_along(res$spending_param),
    function(r) {
        any(abs(spec$param$lattice[[r]] - res$spending_param[[r]]) <= 1e-12)
    },
    logical(1L)
))
cat(
    "\nidentical(res$power, calc_power_pvals_gsd(resolve(...))):",
    identical(res$power, check),
    "\n"
)
cat("spending_param on the lattice:", on_lattice, "\n")
cat(
    "profile rows:",
    nrow(res$spending_profile),
    " converged:",
    res$profile_converged,
    "\n"
)
cat(
    "VERDICT:",
    if (identical(res$power, check) && on_lattice) "PASS" else "FAIL",
    "\n"
)
