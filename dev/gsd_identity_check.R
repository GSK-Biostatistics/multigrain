# dev/gsd_identity_check.R -- the "nothing changed" gate of design record
# section 4.10. Installs `main` (324b7ca) and the current working tree into two
# private libraries, runs the same fixed-sample `graph_optimise()` call under
# each in its own Rscript process, and compares the results with `identical()`.
#
# Run from the repository root:
#     Rscript dev/gsd_identity_check.R
#
# Nothing is written inside the repository and nothing is installed into the
# user's own library: both builds go into temporary directories, and the
# worktree is removed at the end.

repo <- normalizePath(getwd(), winslash = "/")
stopifnot(file.exists(file.path(repo, "DESCRIPTION")))

rscript <- file.path(R.home("bin"), "Rscript")
main_ref <- "324b7ca"

base <- tempfile("gsd_identity_")
worktree <- file.path(base, "main_worktree")
lib_main <- file.path(base, "lib_main")
lib_branch <- file.path(base, "lib_branch")
dir.create(lib_main, recursive = TRUE)
dir.create(lib_branch, recursive = TRUE)

cat("Scratch directory:", base, "\n")

# --- install a package directory into a private library --------------------
# `R CMD INSTALL` installs into the first library of the process that runs it,
# so the target library is put first on R_LIBS for the child process.
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
    if (!file.exists(file.path(lib, "multigrain", "DESCRIPTION"))) {
        stop("multigrain did not land in ", lib, call. = FALSE)
    }

    invisible(TRUE)
}

# --- the driver run under each library -------------------------------------
driver <- file.path(base, "driver.R")
writeLines(
    r"(args <- commandArgs(trailingOnly = TRUE)
lib <- args[[1L]]
out <- args[[2L]]
.libPaths(c(lib, .libPaths()))
library(multigrain, lib.loc = lib)

set.seed(11)
pvals <- simulate_pvalues(
    power_nominal = c(0.9, 0.85, 0.8),
    corr_matrix = diag(3),
    nsim = 2000
)
ts <- trial_success(r1 + r2 + r3, verbose = "silent")
ctrl <- multigrain_control()
ctrl <- control_global(ctrl, maxiter = 20, run = 5, popSize = 30)
ctrl <- control_local(ctrl, maxeval = 200)
ctrl <- control_nsim_global(ctrl, 1000)

result <- graph_optimise(
    pvals = pvals,
    graph_constraint = graph_constraint_free(3),
    trial_success = ts,
    control = ctrl,
    verbose = "silent"
)

saveRDS(
    list(
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

# --- do it -----------------------------------------------------------------
# Wrapped in a function so that `on.exit()` removes the worktree even when a
# step fails: `on.exit()` at the top level of a script never fires.
build_and_run <- function() {
    cat("Creating a worktree of", main_ref, "...\n")
    system2("git", shQuote(c(
        "-C", repo, "worktree", "add", "--detach", worktree, main_ref
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

    cat("Installing", main_ref, "into", lib_main, "...\n")
    install_to(worktree, lib_main)

    cat("Installing the working tree into", lib_branch, "...\n")
    install_to(repo, lib_branch)

    cat("Running the driver under each build ...\n")
    list(
        main = run_driver(lib_main, "main"),
        branch = run_driver(lib_branch, "branch")
    )
}

runs <- build_and_run()
a <- runs$main
b <- runs$branch

cat("\n== graph_optimise() identity check ==\n")
cat("main    :", a$version, "at", a$lib, "\n")
cat("branch  :", b$version, "at", b$lib, "\n")

same_object <- identical(a$result, b$result)
same_seed <- identical(a$seed, b$seed)
cat("identical(result):", same_object, "\n")
cat("identical(.Random.seed):", same_seed, "\n")

if (!same_object) {
    cat("\nFalling back to component comparison.\n")
    differing <- names(a$result)[
        !vapply(
            names(a$result),
            function(nm) identical(a$result[[nm]], b$result[[nm]]),
            logical(1L)
        )
    ]
    cat("Top-level components that are not identical():",
        toString(differing), "\n")

    report <- function(label, x, y) {
        cat(sprintf(
            "  %-22s identical: %-5s  all.equal: %s\n",
            label,
            identical(x, y),
            toString(isTRUE(all.equal(x, y)))
        ))
    }
    report("hyp_weight", a$result$hyp_weight, b$result$hyp_weight)
    report("trans_matrix", a$result$trans_matrix, b$result$trans_matrix)
    report("power", a$result$power, b$result$power)
    report("solution", a$result$solution, b$result$solution)
    if (!is.null(a$result$global_output)) {
        report(
            "global_output@solution",
            a$result$global_output@solution,
            b$result$global_output@solution
        )
    }

    # Only `func` can legitimately differ: it is compiled by
    # `Rcpp::sourceCpp()` and its environment holds an external pointer that
    # is process-specific. Everything else about the gain, including the
    # generated C++ and the body of the compiled function, must match, or the
    # gain itself has changed.
    cat("\n`trial_success`, field by field:\n")
    ts_a <- a$result$trial_success
    ts_b <- b$result$trial_success
    report("trial_success$m", ts_a$m, ts_b$m)
    report("trial_success$objective", ts_a$objective, ts_b$objective)
    report("trial_success$cpp_code", ts_a$cpp_code, ts_b$cpp_code)
    report("body(func)", body(ts_a$func), body(ts_b$func))
    report("trial_success$func", ts_a$func, ts_b$func)

    cat(
        "\nNote: only `trial_success$func` may differ. It is compiled by",
        "\n`Rcpp::sourceCpp()` and carries a process-specific external",
        "pointer, so\n`identical()` on the whole object cannot hold across",
        "processes. Every other\nfield of the gain, including the generated",
        "C++, is compared above.\n"
    )
}

cat("\nDone.\n")
