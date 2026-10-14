# dev/review/hsd_opt/gsd_fixed_path_identity.R -- the group sequential
# "nothing changed" gate of the spending-parameter build (build plan section
# 6 item 4). Installs the branch base 3e55fb5 (gsd-build before any of this
# work) and the current working tree into two private libraries, runs the
# same untagged `transform_pvalues_gsd()` and `graph_optimise_gsd()` calls
# under each in its own Rscript process, and compares the results with
# `identical()`. The fixed-sample gate is `dev/gsd_identity_check.R`; this
# script is its GSD counterpart.
#
# Run from the repository root:
#     Rscript dev/review/hsd_opt/gsd_fixed_path_identity.R
#
# Nothing is written inside the repository and nothing is installed into the
# user's own library.

repo <- normalizePath(getwd(), winslash = "/")
stopifnot(file.exists(file.path(repo, "DESCRIPTION")))

rscript <- file.path(R.home("bin"), "Rscript")
base_ref <- "3e55fb5"

base <- tempfile("gsd_fixed_path_")
worktree <- file.path(base, "base_worktree")
lib_base <- file.path(base, "lib_base")
lib_branch <- file.path(base, "lib_branch")
dir.create(lib_base, recursive = TRUE)
dir.create(lib_branch, recursive = TRUE)

cat("Scratch directory:", base, "\n")

install_to <- function(pkg_dir, lib) {
    script <- file.path(base, paste0("install_", basename(lib), ".R"))
    writeLines(
        c(
            sprintf(".libPaths(c(%s, .libPaths()))", deparse(lib)),
            sprintf(
                paste0(
                    "devtools::install(pkg = %s, quick = TRUE, ",
                    "dependencies = FALSE, upgrade = FALSE, ",
                    "reload = FALSE, quiet = TRUE)"
                ),
                deparse(pkg_dir)
            )
        ),
        script
    )
    r_libs <- paste(
        c(lib, .libPaths()), # nolint: undesirable_function_linter.
        collapse = .Platform$path.sep
    )
    status <- withr::with_envvar(
        c(R_LIBS = r_libs),
        system2(rscript, shQuote(script))
    )
    if (status != 0L) {
        stop("install of ", pkg_dir, " failed", call. = FALSE)
    }
    invisible(TRUE)
}

# The driver: the test-optimisation_gsd.R fixture, untagged spending, with a
# stopping rule so that the whole fixed GSD path is exercised.
driver <- file.path(base, "driver.R")
writeLines(
    r"(args <- commandArgs(trailingOnly = TRUE)
lib <- args[[1L]]
out <- args[[2L]]
.libPaths(c(lib, .libPaths()))
library(multigrain, lib.loc = lib)

set.seed(101)
raw <- simulate_pvalues_gsd(
    power_nominal = c(0.9, 0.85, 0.8),
    corr_matrix = diag(3),
    info_frac = c(0.5, 1),
    nsim = 2000
)
pvals <- transform_pvalues_gsd(
    raw,
    spending = gsDesign::sfLDOF,
    grid_size = 128L
)
gain <- trial_success_gsd(d(t1) + d(t2) + d(t3), d = c(1, 0.75),
                          verbose = "silent")
ctrl <- multigrain_control()
ctrl <- control_global(ctrl, maxiter = 20, run = 5, popSize = 30)
ctrl <- control_local(ctrl, maxeval = 200)
ctrl <- control_nsim_global(ctrl, 1000)

set.seed(202)
result <- graph_optimise_gsd(
    pvals = pvals,
    graph_constraint = graph_constraint_free(3),
    trial_success = gain,
    stopping_rule = stopping_rule_gsd(r1),
    control = ctrl,
    verbose = "silent"
)

saveRDS(
    list(
        pvals = pvals,
        result = result,
        seed = .Random.seed,
        version = as.character(packageVersion("multigrain")),
        lib = lib
    ),
    out
))",
    driver
)

run_driver <- function(lib, tag) {
    out <- file.path(base, paste0(tag, ".rds"))
    status <- system2(rscript, shQuote(c(driver, lib, out)))
    if (status != 0L) {
        stop("driver failed under ", lib, call. = FALSE)
    }
    readRDS(out)
}

build_and_run <- function() {
    cat("Creating a worktree of", base_ref, "...\n")
    system2("git", shQuote(c(
        "-C", repo, "worktree", "add", "--detach", worktree, base_ref
    )))
    on.exit(
        {
            cat("Removing the worktree ...\n")
            system2("git", shQuote(c(
                "-C", repo, "worktree", "remove", "--force", worktree
            )))
        },
        add = TRUE
    )
    cat("Installing", base_ref, "into", lib_base, "...\n")
    install_to(worktree, lib_base)
    cat("Installing the working tree into", lib_branch, "...\n")
    install_to(repo, lib_branch)
    cat("Running the driver under each build ...\n")
    list(
        base = run_driver(lib_base, "base"),
        branch = run_driver(lib_branch, "branch")
    )
}

runs <- build_and_run()
a <- runs$base
b <- runs$branch

cat("\n== group sequential fixed-path identity check ==\n")
cat("base    :", a$version, "at", a$lib, "\n")
cat("branch  :", b$version, "at", b$lib, "\n")
cat("identical(transform_pvalues_gsd() object):", identical(a$pvals, b$pvals), "\n")
cat("identical(.Random.seed):", identical(a$seed, b$seed), "\n")

differing <- names(a$result)[
    !vapply(
        names(a$result),
        function(nm) identical(a$result[[nm]], b$result[[nm]]),
        logical(1L)
    )
]
cat("Result components that are not identical():", toString(differing), "\n")
cat("Result components present on the branch only:",
    toString(setdiff(names(b$result), names(a$result))), "\n")

for (nm in c("hyp_weight", "trans_matrix", "power", "solution",
             "stopping_rule")) {
    cat(sprintf("  %-16s identical: %s\n", nm,
                identical(a$result[[nm]], b$result[[nm]])))
}
cat(sprintf("  %-16s identical: %s\n", "global@solution",
            identical(a$result$global_output@solution,
                      b$result$global_output@solution)))

ts_a <- a$result$trial_success
ts_b <- b$result$trial_success
cat("\n`trial_success`, field by field:\n")
for (nm in c("m", "objective", "cpp_code")) {
    cat(sprintf("  %-16s identical: %s\n", nm,
                identical(ts_a[[nm]], ts_b[[nm]])))
}
cat(sprintf("  %-16s identical: %s  (external pointer; may differ)\n",
            "func", identical(ts_a$func, ts_b$func)))

ok <- identical(a$pvals, b$pvals) &&
    identical(a$seed, b$seed) &&
    setequal(differing, "trial_success") &&
    length(setdiff(names(b$result), names(a$result))) == 0L
cat("\nVERDICT:", if (ok) "PASS" else "FAIL", "\n")
