# F1 gate for families of hypotheses: checks too large for tests/.
#
# Usage: Rscript dev/review/families/verify.R [<library with multigrain>]
#
#   1. Kernel vs literal Algorithm 5.7 (tests/testthat/helper-families.R):
#      60 random cases per within-family test, lowest-index and random
#      rejection order, m 2-12, M 1-6, non-contiguous members, rows of G
#      summing to 1 or below 1, 400 trials per case.
#   2. All-singleton families identical() to graph_shortcut(): 40 random
#      graphs, m 2-10, 2000 trials each.
#   3. adjust_pvals_families() identical() to stats::p.adjust(): k 1-7,
#      10,000 rows per method and k.
#   4. Serial identical() to parallel: 20 random cases, 1e5 trials, 1-4
#      threads.
#
# Output: dev/review/families/results/verify.txt

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

log_file <- file.path(out_dir, "verify.txt")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

cat("multigrain", format(utils::packageVersion("multigrain")), "\n")
cat(format(Sys.time()), "\n\n")

with(env, {
    set.seed(20261019)
    tests <- c("hommel", "hochberg", "holm")

    # 1. Literal Algorithm 5.7 ----------------------------------------------
    cat("1. Kernel vs literal Algorithm 5.7 (400 trials per case)\n")
    for (test in tests) {
        for (order in c("lowest", "random")) {
            n_match <- 0L
            for (case in 1:60) {
                m <- sample(2:12, 1)
                num_fam <- sample.int(min(m, 6), 1)
                members <- random_partition(m, num_fam)
                gr <- random_family_graph(num_fam, full_rows = case %% 2 == 0)
                pvals <- matrix(stats::runif(400 * m)^4, ncol = m)

                kernel <- shortcut_families_r(pvals, members, test, gr$w, gr$G)
                oracle <- alg57_oracle(
                    pvals,
                    members,
                    test,
                    gr$w,
                    gr$G,
                    order = order
                )
                n_match <- n_match + identical(unname(kernel), oracle)
            }
            cat(sprintf(
                "   %-8s %-6s order: %d / 60 identical\n",
                test,
                order,
                n_match
            ))
        }
    }

    # 2. All-singleton families vs graph_shortcut() --------------------------
    n_match <- 0L
    for (case in 1:40) {
        m <- sample(2:10, 1)
        gr <- random_family_graph(m, full_rows = case %% 2 == 0)
        pvals <- matrix(stats::runif(2000 * m)^4, ncol = m)
        n_match <- n_match +
            identical(
                graph_shortcut_families(pvals, 0:(m - 1), 0.025, gr$w, gr$G),
                graph_shortcut(pvals, 0.025, gr$w, gr$G)
            )
    }
    cat(sprintf(
        "\n2. All-singleton vs graph_shortcut(): %d / 40 identical\n",
        n_match
    ))

    # 3. Within-family adjustment vs stats::p.adjust() -----------------------
    cat("\n3. adjust_pvals_families() vs stats::p.adjust() (10,000 rows)\n")
    for (test in tests) {
        res <- vapply(
            1:7,
            function(k) {
                pvals <- matrix(stats::runif(1e4 * k)^2, ncol = k)
                expected <- t(apply(pvals, 1, stats::p.adjust, method = test))
                if (k == 1) {
                    expected <- t(expected)
                }
                identical(
                    adjust_families_r(pvals, list(seq_len(k)), test),
                    unname(expected)
                )
            },
            logical(1)
        )
        cat(sprintf(
            "   %-8s k = 1..7 identical: %s\n",
            test,
            paste(res, collapse = " ")
        ))
    }

    # 4. Serial vs parallel --------------------------------------------------
    n_match <- 0L
    for (case in 1:20) {
        m <- sample(2:12, 1)
        num_fam <- sample.int(min(m, 6), 1)
        members <- random_partition(m, num_fam)
        gr <- random_family_graph(num_fam)
        padj <- adjust_families_r(
            matrix(stats::runif(1e5 * m)^4, ncol = m),
            members,
            sample(tests, 1)
        )
        fam_index <- family_index0(members)
        serial <- graph_shortcut_families(padj, fam_index, 0.025, gr$w, gr$G)
        n_match <- n_match +
            all(vapply(
                1:4,
                function(nt) {
                    identical(
                        graph_shortcut_families_parallel(
                            padj,
                            fam_index,
                            0.025,
                            gr$w,
                            gr$G,
                            num_threads = nt
                        ),
                        serial
                    )
                },
                logical(1)
            ))
    }
    cat(sprintf(
        "\n4. Serial vs parallel (1-4 threads): %d / 20 identical\n",
        n_match
    ))
})
