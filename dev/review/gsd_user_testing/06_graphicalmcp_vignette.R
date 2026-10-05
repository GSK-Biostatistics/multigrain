# 06_graphicalmcp_vignette.R
# The two case studies of the graphicalMCP vignette "Group Sequential Design
# with Graphical Approaches" (graphicalMCP 0.3.0), run through multigrain's
# transform and kernel and compared with graphicalMCP's own answer.
#
#   A. Maurer and Bretz (2013) diabetes trial: 4 hypotheses, analyses at 1/3
#      and 2/3 of the information (the trial stops after the second), no
#      look-back.
#   B. Oncology trial: 6 hypotheses with 3, 3, 2, 2, 1 and 1 analyses, NA
#      padding, look-back on. Three sets of p-values: the vignette's first set,
#      its "look-back makes a difference" set, and that set without look-back.
#
# Conventions that differ between the two packages (design record 4.9):
#   * graphicalMCP pads a finished endpoint with NA and, without look-back,
#     never tests it again. multigrain carries its repeated p-value forward, so
#     it can still be rejected when alpha reaches it later. The two agree when
#     look-back is on, and can differ when it is off.
#   * multigrain's decision time is the analysis at which the rejection was
#     declared (graphicalMCP's `decision_at`), not the earliest analysis whose
#     boundary was crossed (`first_rejected_at`).
#
# Run from the repository root:
#   Rscript dev/review/gsd_user_testing/06_graphicalmcp_vignette.R

suppressPackageStartupMessages({
    library(multigrain)
    library(graphicalMCP)
})

here <- "dev/review/gsd_user_testing"
if (!dir.exists(here)) here <- "."

alpha <- 0.025
ldof  <- graphicalMCP::spending_of     # plain numeric cumulative spend

# One observed trial through multigrain: p is m x K (NA allowed).
run_multigrain <- function(p, info_frac, w, G, look_back) {
    m <- nrow(p); K <- ncol(p)
    raw <- array(p, dim = c(1L, m, K))
    info <- if (is.matrix(info_frac)) info_frac else matrix(info_frac, m, K, byrow = TRUE)
    # multigrain wants the information fraction of a finished endpoint carried
    # forward (it is still "at full information"), where graphicalMCP has NA
    for (i in seq_len(m)) {
        done <- which(info[i, ] >= 1)
        if (length(done)) {
            k0 <- done[1]
            if (k0 < K) {
                info[i, (k0 + 1):K] <- 1
                raw[1, i, (k0 + 1):K] <- raw[1, i, k0]
            }
        }
    }
    pv <- transform_pvalues_gsd(raw, info_frac = info, spending = ldof,
                                alpha = alpha, look_back = look_back)
    x <- pv$pvals
    dim(x) <- c(1L, m * K)
    out <- multigrain:::graph_shortcut_gsd(x, alpha, w, G, K)
    list(rejected = as.vector(out$rejected), time = as.vector(out$time),
         repeated = matrix(pv$pvals, m, K))
}

run_oracle <- function(p, info_frac, w, G, look_back) {
    g <- graphicalMCP::graph_create(w, G)
    pm <- p
    dimnames(pm) <- list(names(g$hypotheses), NULL)
    inf <- info_frac
    if (is.matrix(inf)) dimnames(inf) <- list(names(g$hypotheses), NULL)
    o <- graphicalMCP::graph_test_shortcut_gsd(
        graph = g, p = pm, alpha = alpha, info_frac = inf,
        spending_fn = ldof, look_back = look_back
    )
    list(rejected = unname(o$outputs$rejected),
         decision_at = as.integer(unname(o$outputs$decision_at)),
         first_rejected_at = as.integer(unname(o$outputs$first_rejected_at)))
}

compare <- function(label, p, info_frac, w, G, look_back, hyp) {
    mg <- run_multigrain(p, info_frac, w, G, look_back)
    or <- run_oracle(p, info_frac, w, G, look_back)
    tab <- data.frame(
        hypothesis = hyp,
        multigrain_rejected = mg$rejected,
        graphicalMCP_rejected = or$rejected,
        multigrain_time = mg$time,
        graphicalMCP_decision_at = ifelse(or$rejected, or$decision_at, 0L),
        graphicalMCP_first_rejected_at = ifelse(is.na(or$first_rejected_at), 0L,
                                                or$first_rejected_at)
    )
    cat("\n== ", label, "\n", sep = "")
    print(tab, row.names = FALSE)
    same_rej  <- identical(mg$rejected, or$rejected)
    same_time <- identical(mg$time[or$rejected & mg$rejected],
                           or$decision_at[or$rejected & mg$rejected])
    cat(sprintf("rejections agree: %s; decision times agree on jointly rejected hypotheses: %s\n",
                same_rej, same_time))
    invisible(list(table = tab, same_rejections = same_rej, same_times = same_time,
                   repeated = mg$repeated))
}

# ---- A. Maurer and Bretz diabetes case study --------------------------------
w_a <- c(0.5, 0.5, 0, 0)
G_a <- rbind(c(0, 0.5, 0.5, 0), c(0.5, 0, 0, 0.5), c(0, 1, 0, 0), c(1, 0, 0, 0))
p_a <- rbind(c(0.0062, 0.0002), c(0.017, 0.0035), c(0.009, 0.002), c(0.13, 0.06))
res_a <- compare("A. Diabetes trial (Maurer and Bretz 2013), analyses at 1/3 and 2/3, no look-back",
                 p_a, c(1 / 3, 2 / 3), w_a, G_a, FALSE, paste0("H", 1:4))
cat("vignette: H1, H2 and H3 rejected at analysis 2, H4 retained\n")

# Repeated p-values against the paper's Table 2. multigrain reports 1 above
# alpha, so the table is rebuilt up to 0.5 for this comparison only.
raw_a <- array(p_a, dim = c(1L, 4L, 2L))
wide  <- transform_pvalues_gsd(raw_a, info_frac = c(1 / 3, 2 / 3), spending = ldof,
                               alpha = 0.5)$pvals
paper <- cbind(c(0.1141, 0.1683, 0.1316, 0.382), c(0.0024, 0.0172, 0.0117, 0.1285))
cat("\nRepeated p-values, multigrain (table built to 0.5) against Maurer and Bretz Table 2:\n")
print(data.frame(hypothesis = paste0("H", 1:4),
                 look1 = signif(wide[1, , 1], 4), paper1 = paper[, 1],
                 look2 = signif(wide[1, , 2], 4), paper2 = paper[, 2]),
      row.names = FALSE)
cat(sprintf("largest absolute difference from the paper: %.1e\n",
            max(abs(matrix(wide, 4, 2) - paper))))

# ---- B. Oncology case study ---------------------------------------------------
hyp_b <- c("H1_OS_S", "H2_OS_A", "H3_PFS_S", "H4_PFS_A", "H5_ORR_S", "H6_ORR_A")
w_b <- c(0.01, 0.01, 0.004, 0, 0.0005, 0.0005) / alpha
G_b <- rbind(
    c(0, 1, 0, 0, 0, 0),
    c(0, 0, 0.5, 0.5, 0, 0),
    c(0, 0, 0, 1, 0, 0),
    c(0, 0, 0, 0, 0.5, 0.5),
    c(0, 0, 0, 0, 0, 1),
    c(0.5, 0.5, 0, 0, 0, 0)
)
p_b <- rbind(
    c(0.03, 0.0001, 0.000001),
    c(0.2, 0.15, 0.1),
    c(0.2, 0.001, NA),
    c(0.3, 0.2, NA),
    c(0.00001, NA, NA),
    c(0.1, NA, NA)
)
info_b <- rbind(
    c(185 / 295, 245 / 295, 1),
    c(529 / 800, 700 / 800, 1),
    c(265 / 310, 1, NA),
    c(675 / 750, 1, NA),
    c(1, NA, NA),
    c(1, NA, NA)
)
res_b1 <- compare("B1. Oncology trial, vignette p-values, look-back on",
                  p_b, info_b, w_b, G_b, TRUE, hyp_b)
cat("vignette: H5 rejected at analysis 1; H1 and H3 at analysis 2; H2, H4 and H6 retained\n")

p_b2 <- p_b
p_b2[2, ] <- c(0.003, 0.005, 0.1)
p_b2[4, ] <- c(0.0001, 0.02, NA)
p_b2[5, ] <- c(0.0008, NA, NA)
res_b2 <- compare("B2. Oncology trial, 'look-back makes a difference' p-values, look-back on",
                  p_b2, info_b, w_b, G_b, TRUE, hyp_b)
cat("vignette: with look-back H2, H4 and H5 are rejected as well\n")

res_b3 <- compare("B3. The same p-values, look-back off (conventions differ for finished endpoints)",
                  p_b2, info_b, w_b, G_b, FALSE, hyp_b)
cat("vignette: without look-back H4 and H5 are not rejected, H2 is rejected at analysis 2\n")

saveRDS(list(A = res_a, B1 = res_b1, B2 = res_b2, B3 = res_b3),
        file.path(here, "06_graphicalmcp_vignette_results.rds"))
