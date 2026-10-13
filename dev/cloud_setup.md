# Setting up a Claude Code cloud container for `multigrain`

The cloud container is Ubuntu 24.04 with R 4.3.3 (`r-base-core 4.3.3-2build2`), gcc, make
and pandoc, but no R packages. Installing every dependency from CRAN source takes 30 to
40 minutes on the 4-core container because stringi, ragg, igraph, RcppArmadillo, V8 and
friends compile from scratch. The script below gets the same result in a few minutes by
taking prebuilt binaries from Ubuntu's own `r-cran-*` packages (built for the container's
R 4.3.3) and compiling only the four CRAN packages Ubuntu does not ship, plus the three
development tools where the version matters.

Paste it into the environment's settings as the **Setup script** (cloud environment menu in
the session's title bar, then Edit). It runs when each new session starts, so a session is
ready as soon as it opens. The working directory of the script is the repository clone.

## What the proxy allows, measured 2026-10-13

| Source | Reachable | Used for |
|---|---|---|
| `archive.ubuntu.com` (apt) | yes | system libraries, `r-cran-*` binaries |
| `cloud.r-project.org` (CRAN source tarballs) | yes | packages not in Ubuntu |
| Posit Package Manager binaries (`rspm-sync.rstudio.com`) | no (403) | nothing |
| r2u (`r2u.stat.illinois.edu`) | no | nothing |
| GitHub release downloads (for `air`) | no (403) | nothing: format by hand, lint with lintr |

## The setup script

```bash
#!/usr/bin/env bash
# multigrain cloud container setup: Ubuntu binaries first, CRAN source for the rest.
set -euo pipefail
export DEBIAN_FRONTEND=noninteractive

# 1. System libraries that R packages link against (most are already in the image;
#    libuv1-dev and libtiff-dev were missing on 2026-10-13).
apt-get update -qq
apt-get install -y -qq \
    libuv1-dev libtiff-dev libcurl4-openssl-dev libssl-dev libxml2-dev \
    libfontconfig1-dev libfreetype-dev libharfbuzz-dev libfribidi-dev \
    libpng-dev libjpeg-dev libgit2-dev libglpk-dev libgmp-dev libicu-dev

# 2. Prebuilt R packages from Ubuntu universe (Imports, Suggests, dev tools, and their
#    compiled dependencies: stringi, ragg, systemfonts, RcppArmadillo, igraph, ...).
apt-get install -y -qq \
    r-cran-cli r-cran-dplyr r-cran-ggplot2 r-cran-ggraph r-cran-glue r-cran-igraph \
    r-cran-mvtnorm r-cran-nloptr r-cran-purrr r-cran-rcpp r-cran-rcppparallel \
    r-cran-rlang r-cran-tibble r-cran-tidyr r-cran-vctrs r-cran-withr r-cran-spelling \
    r-cran-svglite r-cran-testthat r-cran-vdiffr r-cran-devtools r-cran-roxygen2 \
    r-cran-lintr r-cran-cyclocomp r-cran-remotes r-cran-pkgload r-cran-pkgbuild \
    r-cran-knitr r-cran-rmarkdown r-cran-pkgdown

# 3. From CRAN source, upgrading what Ubuntu ships too old. Ubuntu's rlang, cli, vctrs,
#    withr, fs, R6, ps, processx, usethis and about thirty others are below the minimum
#    versions that current testthat, roxygen2, gsDesign and friends require, so a plain
#    install.packages() of the targets fails at byte-compilation. remotes with
#    upgrade = "always" resolves the dependency closure of DESCRIPTION and of the dev
#    tools and rebuilds only what is outdated or missing (86 packages, mostly pure R).
#    devtools itself stays at Ubuntu's 2.4.5: it delegates to the upgraded pkgload,
#    roxygen2 and testthat, and the CRAN version would also need usethis rebuilt.
Rscript -e '
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 4)
remotes::install_deps(".", dependencies = TRUE, upgrade = "always")
remotes::install_cran(c("roxygen2", "lintr", "cyclocomp", "testthat"), upgrade = "always")
for (p in c("Rcpp", "RcppParallel", "gsDesign", "GA", "nloptr", "testthat", "devtools",
            "roxygen2", "lintr", "cyclocomp", "graphicalMCP", "gMCPLite", "mvtnorm",
            "vdiffr", "svglite"))
    if (!requireNamespace(p, quietly = TRUE)) stop("missing package: ", p)
cat("roxygen2", as.character(packageVersion("roxygen2")),
    "testthat", as.character(packageVersion("testthat")),
    "RcppParallel", as.character(packageVersion("RcppParallel")), "\n")
'

# 4. Build and install the package once so the compiled kernels are ready.
if [ -f DESCRIPTION ]; then
    LC_ALL=C.UTF-8 Rscript -e '
    devtools::clean_dll()
    devtools::install(quick = TRUE, dependencies = FALSE, quiet = TRUE)
    cat("multigrain installed:", as.character(packageVersion("multigrain")), "\n")
    '
fi
```

## Notes for the session that uses it

- Set `Sys.setenv(NOT_CRAN = "true")` before `testthat::test_file()`, or every
  `skip_on_cran()` block silently skips. Use `devtools::load_all()` and `test_file()` on
  the files you touch; never `testthat::test_package()`.
- Run `devtools::document()` under `LC_ALL=C.UTF-8`: a POSIX locale rewrites
  `man/multigrain-package.Rd` with `<U+0219>` escapes.
- Ubuntu's `r-cran-*` versions are older than CRAN (for example `r-cran-roxygen2` is
  7.3.1 and `r-cran-lintr` 3.1.1). Step 3 installs the CRAN versions on top, into
  `/usr/local/lib/R/site-library`, which comes first on `.libPaths()`. If `devtools::document()`
  wants to change `Config/roxygen2/version` in `DESCRIPTION`, the roxygen2 in use is not
  the pinned one.
- `graphicalMCP (>= 0.3.0)` is a Suggests used only as a test oracle; if it ever fails to
  install, the oracle tests skip and the gate counts in `dev/gsd_progress.md` will not
  match.
- After any change under `src/`, run `devtools::clean_dll()` before `devtools::load_all()`.
- The CRAN source route for everything (what this branch's first session did) is the
  fallback if apt is ever blocked:
  `Rscript -e 'options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 4); install.packages(c("remotes", "devtools", "testthat", "roxygen2", "lintr", "cyclocomp")); remotes::install_deps(".", dependencies = TRUE)'`
  after `apt-get install -y libuv1-dev libtiff-dev`. Expect 30 to 40 minutes.

## Measured on 2026-10-13 (4 cores)

| Step | Time | Result |
|---|---|---|
| CRAN source for everything (the fallback) | about 35 min | 150 packages |
| 2. apt `r-cran-*` binaries | 2 min 27 s | 663 Debian packages (the closure of the list above) |
| 3. CRAN source upgrades on top of apt | 7 min 19 s | 86 packages rebuilt; roxygen2 8.1.1, testthat 3.3.2, lintr 3.4.0, RcppParallel 6.2.1, gsDesign 3.11.0, graphicalMCP 0.3.0, GA 3.2.5 |
| 4. build multigrain | under 1 min | |

So a session is ready in about ten minutes instead of thirty-five. Three earlier variants of
step 3 failed and are recorded here so nobody repeats them: `install.packages()` of just the
seven missing targets (rlang too old for gsDesign and testthat); the same with rlang, cli,
lifecycle, vctrs and withr added in the same call (parallel builds began before the new
rlang was in place); the same as two sequential calls (then fs and R6 were too old). Only
the upgrade-aware `remotes` install resolved the closure; its one failure, the CRAN
devtools needing usethis >= 3.2.1, is why devtools is left at Ubuntu's version.
