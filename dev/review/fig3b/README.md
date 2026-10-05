# Reproducible code for "Gain-function optimisation of graphical multiple testing procedures for confirmatory clinical trials"

Spiers, Grayling, Wheeler, Mander — submitted to *Statistics in Medicine*.

---

## 1. What this reproduces

| Example | Section | Optimise script | Figure script | Figures | Table |
|--------:|---------|-----------------|---------------|---------|-------|
| 1 | 3.1 Enriched subgroup | `Ex1-optimise.R` | `Ex1-figures.R` | `subgroup_fig.pdf` | Table 2 |
| 2 | 3.2 Sample size (FIBRONEER) | `Ex2-optimise.R` | `Ex2-figures.R` | `sample_size_fig.pdf` | inline (n\*, N_total) |
| 3 | 3.3 Elicitation (tralokinumab) | `Ex3-optimise.R` | `Ex3-figures.R` | `elicit_fig.pdf` | Table 4 |
| 4 | 3.4 Bayes gain | `Ex4-optimise.R` | `Ex4-figures.R` | `bayes_fig.pdf`, `tau_sensitivity_fig.pdf` | Table 5 |
| 5 | 3.5 GSD (PFS/OS) | `Ex5-optimise.R` | `Ex5-figures.R` | `two_node_graph_bw.pdf`, `gsd_w1_wide.pdf` | Table 6 |

`00-shared-settings.R` holds shared constants and is sourced by every script.

---

## 2. Requirements and installation

- **R** ≥ 4.2; a C++ toolchain (Rtools on Windows).

```r
remotes::install_github("GSK-Biostatistics/multigrain")
# or, from the supplied Supporting Information source package:
install.packages("multigrain-sample-size-pub.zip", repos = NULL, type = "source")
```

> **Note.** Example 2 requires `graph_optimise_n()`, available in the
> `sample-size-pub` build above.

---

## 3. How to run

Run from inside `submission/`.

### Reproduce figures and tables (fast path)

Precomputed results ship as `results/Ex{1..5}-results.rds`. Run only the figure scripts:

```bash
cd submission
Rscript Ex1-figures.R
Rscript Ex2-figures.R
Rscript Ex3-figures.R
Rscript Ex4-figures.R
Rscript Ex5-figures.R
```

### Re-run optimisations (optional, hours, requires HPC)

Run `Rscript Ex{N}-optimise.R` before the corresponding figure script. Examples 1, 3, and 4 can take hours; Example 5 takes minutes. Set `NUM_THREADS` in `00-shared-settings.R` to match your machine (default `60`). Runtime can be reduced by lowering `nsim_local` and `nsim_global`, at the cost of precision.
---