# Optimised graph for the "Families of hypotheses" article
# (vignettes/articles/families.qmd).
#
# The p-values, families, graph constraint and trial success measure are the
# same as in the article; keep the two in sync. The result is saved to
# vignettes/articles/data/families-graph.rds and read by the article.
#
# Run from the package root with the development version installed:
#   Rscript data-raw/families-vignette.R

library(multigrain)

# Hypotheses: each is a comparison of a triple-therapy dose with a dual
# therapy (A or B) on one endpoint.
comparison <- c(
    "high_vs_A", "high_vs_B", "high_vs_B", "high_vs_A", "high_vs_A", "high_vs_A",
    "low_vs_A", "low_vs_B", "low_vs_B", "low_vs_A", "low_vs_A", "low_vs_A"
)
endpoint <- rep(
    c("exac", "exac", "fev1_auc", "fev1_trough", "time_exac", "sgrq"),
    times = 2
)

# Comparisons that share an arm (equal allocation) have correlation 0.5
corr_comparison <- matrix(
    c(
        1, 0.5, 0.5, 0,
        0.5, 1, 0, 0.5,
        0.5, 0, 1, 0.5,
        0, 0.5, 0.5, 1
    ),
    nrow = 4,
    dimnames = rep(list(c("high_vs_A", "high_vs_B", "low_vs_A", "low_vs_B")), 2)
)

# Endpoints: 0.6 between the two exacerbation endpoints, 0.7 between the two
# lung function endpoints, 0.2 otherwise
endpoints <- c("exac", "time_exac", "fev1_auc", "fev1_trough", "sgrq")
corr_endpoint <- matrix(0.2, 5, 5, dimnames = list(endpoints, endpoints))
diag(corr_endpoint) <- 1
corr_endpoint["exac", "time_exac"] <- 0.6
corr_endpoint["time_exac", "exac"] <- 0.6
corr_endpoint["fev1_auc", "fev1_trough"] <- 0.7
corr_endpoint["fev1_trough", "fev1_auc"] <- 0.7

corr_matrix <- corr_comparison[comparison, comparison] *
    corr_endpoint[endpoint, endpoint]
dimnames(corr_matrix) <- NULL

set.seed(2020)
pvals <- simulate_pvalues(
    power_nominal = c(
        0.85, 0.75, 0.95, 0.80, 0.70, 0.60, # high dose
        0.80, 0.70, 0.95, 0.80, 0.65, 0.55 # low dose
    ),
    corr_matrix = corr_matrix,
    nsim = 1e5
)

fams <- family_spec(
    list(
        exac_high = 1:2,
        lung_high = 3:4,
        sec_high = 5:6,
        exac_low = 7:8,
        lung_low = 9:10,
        sec_low = 11:12
    ),
    test = c("hommel", "hommel", "hochberg", "hommel", "hommel", "hochberg")
)

# All initial weight on the high-dose exacerbation family; NA = free edge
fam_constraint <- graph_constraint(
    hyp_constraint = c(1, 0, 0, 0, 0, 0),
    trans_constraint = matrix(
        c(
            0, NA, 0, NA, 0, 0, # exac_high -> lung_high, exac_low
            0, 0, NA, NA, 0, 0, # lung_high -> sec_high, exac_low
            0, 0, 0, NA, NA, NA, # sec_high  -> low-dose families
            0, NA, 0, 0, NA, 0, # exac_low  -> lung_high, lung_low
            0, NA, 0, 0, 0, NA, # lung_low  -> lung_high, sec_low
            0, NA, NA, 0, NA, 0 # sec_low   -> lung_high, sec_high, lung_low
        ),
        nrow = 6,
        byrow = TRUE
    ),
    names = names(fams$members)
)

gain <- trial_success(
    4 * (r1 && r2) + 2 * (r7 && r8) +
        r3 + r4 + r9 + r10 +
        0.5 * (r5 + r6 + r11 + r12)
)

set.seed(1)
graph_fam <- graph_optimise(
    pvals = pvals,
    graph_constraint = fam_constraint,
    trial_success = gain,
    families = fams,
    num_threads = 4
)

saveRDS(graph_fam, "vignettes/articles/data/families-graph.rds")
