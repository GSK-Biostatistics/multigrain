# Brief P2: objective closure and slice bank

Branch `gsd-build-spending`, repo `/home/user/multigrain-spending`. You implement phase P2 of
`dev/hsd_opt_build_plan.md`. P1 is committed: read its code before anything else. **Report, do
not commit.** Complete files, not diffs.

## Read first (~40 minutes)

1. `dev/hsd_opt_build_plan.md` §1 (rows 2, 3, 6 bind you), §4 "Shared test fixture" and "P2".
2. `dev/hsd_opt_design_revised.md` §1.3 "Stage G", §2.3, §3 items 1, 2, 5, 10, §5.1
   (`R/objective_function_gsd_opt.R`), §5.2 "Clamping", "Archive", "Mutation on the parameter
   gene", "Memory guard", Appendix A.2.
3. `R/transform_pvalues_gsd_opt.R` (P1: the spec fields, `resolve_pvalues_gsd()`,
   `.gsd_invert_hyp()`, lattices, `node_tables` keys) and `tests/testthat/helper-hsd_opt.R`.
4. `R/objective_function_gsd.R` in full — you will copy its decode-and-penalise block verbatim.
   `R/objective_function.R` (`split_theta()`, `recover_full_weights()`,
   `recover_full_trans_matrix()`); `R/check_gsd.R` `.gsd_kernel_matrix()` (column layout
   `(k-1)*m + i`); `R/optimisation_gsd.R` `.graph_optimise_ga_gsd()` (:407-512) for how the
   closure is wired into `GA::ga()`, and `.sample_pvals_gsd()` (:388).
5. `tests/testthat/test-objective_function_gsd.R` (style of the oracle tests).
6. `dev/review/hsd_opt/d03_memory_and_copy_output.txt` (copy cost vs kernel; whether a column
   assignment copies the matrix).

## What to build — new file `R/objective_function_gsd_opt.R`

- `.hsd_node_index(u, M)`: `u <- pmin(pmax(u, 0), 1)` (vectorised), then
  `pmin(M, 1L + floor(u * M))` as integer. Documented inline: equal prior mass per node;
  `u = 1` maps to node `M`.
- `.hsd_bank(spec, idx)`: for each **active** hypothesis `i` (`!is.na(spec$group[i])`) and
  each node `q` of its lattice, the `length(idx) × K` matrix of repeated p-values obtained by
  `.gsd_invert_hyp()` on `spec$raw[idx, i, , drop = FALSE]` with the table
  `spec$node_tables[[key(i, node)]]` (reuse P1's key helper — do not reimplement the key).
  Returns `list(slices = list over active i of list over q of matrix, idx = idx, bytes = …)`.
  Rows are in `idx` order. Call `.check_bank_memory(bytes)` before building.
- `.check_bank_memory(bytes, max_bank_gb = 4, call)`: `cli::cli_warn` above 2 GB,
  `cli::cli_abort` above `max_bank_gb` GB, both naming the three remedies (coarser `step`,
  narrower `bounds`, smaller `nsim_global`). Bytes formula `8 * K * N_idx * sum over active i of
  M_{g(i)}`.
- `create_obj_func_gsd_opt(spec, bank, power_criterion, hyp_constraint, trans_constraint,
  alpha, num_threads = 1L, stop_rule = NULL, p_graph, start_param)`:
  - `force()` everything as `create_obj_func_gsd()` does.
  - Build the working matrix once, **inside the factory body** (not passed in as an argument):
    `.gsd_kernel_matrix()` of `resolve_pvalues_gsd(spec, start_param)` subset to `bank$idx` rows
    (`drop = FALSE` on the array before `dim<-`). Fixed hypotheses' columns never change
    afterwards. **Assign columns with `mat[, cols] <<- slice`** from the closure. d03 measured
    (`dev/review/hsd_opt/d03_memory_and_copy_output.txt`): a matrix built inside the factory
    and updated with `<<-` is never copied; a plain `<-` inside the closure copies the whole
    matrix on every evaluation; a matrix passed in as a factory argument is copied once. Do not
    bind the matrix to a second name anywhere.
  - Precompute for each active hypothesis the column indices `(k - 1) * m + i` for `k = 1..K`
    and its group `r`.
  - Returned closure `function(x)`: `xg <- x[seq_len(p_graph)]`; `u <- x[p_graph +
    seq_len(d)]`; `q <- .hsd_node_index(u, M_r)` per group; **then the decode-and-penalise block
    of `create_obj_func_gsd()` copied verbatim on `xg`** (returns the same penalties before any
    kernel work); then for each active `i`: `mat[, cols_i] <- bank$slices[[i]][[q[r_i]]]` (skip
    when `q` is unchanged from the previous call for that group — keep `last_q` in the closure);
    kernel call exactly as the fixed closure (`graph_shortcut_gsd()` or `_parallel()` with
    `grain_size = -1L`); `fit <- power_criterion(res$time)`; archive update; return `fit`.
  - Archive: an environment `arch` with, per axis `r`, `best` (numeric M_r, init `-Inf`),
    `visits` (integer M_r), `genes` (matrix M_r × p_graph, init NA): on each evaluation that
    reached the kernel, `visits[q_r] += 1`; if `fit > best[q_r]` then `best[q_r] <- fit;
    genes[q_r, ] <- xg`. Expose it as `attr(fn, "archive")` and `attr(fn, "d")`,
    `attr(fn, "p_graph")`.
  - The `d` and the per-group `M_r` come from `spec$param`.
- No RNG. No change to any existing file.

## Tests — `tests/testthat/test-objective_function_gsd_opt.R`

Use the P1 fixtures (`hsd_opt_spec()` d=1 with 5 nodes, the d=2 fixture). `idx <-
withr::with_seed(3, sample.int(spec$nsim, 1000))`. Tiny `graph_constraint_free(m)`, gain
`trial_success_gsd(r1 + r2 + r3)` (and the m=4 gain for the d=2 fixture), `p_graph` from
`ncol(multigrain:::.build_start_matrix(gc, <default start_graph>))`.
- [ ] **node oracle** (the gate): for 20 random graph gene vectors `x` (`runif(p_graph)`) and
  every node `q` of every axis (and, for d=2, every pair of nodes), `identical(fn(c(x, u_q)),
  ref_q(x))` where `ref_q <- create_obj_func_gsd(m, gain$func, hc, tc, pvals =
  .gsd_kernel_matrix(sub_q), K, alpha, 1L, stop_rule)` and `sub_q` is `resolve_pvalues_gsd(spec,
  node_q)` subset to `idx` (build a `multigrain_pvals_gsd` with `$pvals[idx, , , drop =
  FALSE]` and `$nsim <- length(idx)`). `u_q <- (q - 0.5) / M`. Also with a stopping rule table
  (`.gsd_stop_table()` of `stopping_rule_gsd(r1)`) and with `num_threads = 2L`.
- [ ] genes outside [0, 1]: `u = -0.3` → node 1, `u = 1.7` → node M, `u = 1` → node M, `u = 0`
  → node 1 (check through the oracle values).
- [ ] penalties: a gene vector with a negative weight gene returns the same penalty as the fixed
  closure (no kernel call; `visits` unchanged).
- [ ] archive: replay 30 random evaluations by hand and compare `best`, `visits`, `genes`.
- [ ] `last_q` short-cut does not change results: 50 evaluations alternating nodes equal the
  oracle.
- [ ] the working matrix is not copied per evaluation: `tracemem()` the matrix inside the
  closure's environment (reach it with `environment(fn)$mat`) and assert no copy output over
  10 evaluations (skip on platforms where `tracemem` is unavailable: `skip_if_not(
  capabilities("profmem"))`).
- [ ] bank: dims of every slice, `bytes` formula, `.check_bank_memory()` warns at 2.5e9 and
  aborts at 4.5e9 (snapshot of both messages, `expect_snapshot(error = TRUE)` for the abort).
- [ ] unchanged: `test-objective_function_gsd.R` 21.

## Gate (run before reporting; the orchestrator re-runs from a clean build)

```r
Sys.setenv(NOT_CRAN = "true"); devtools::load_all()
testthat::test_file("tests/testthat/test-objective_function_gsd_opt.R", package = "multigrain")
testthat::test_file("tests/testthat/test-objective_function_gsd.R", package = "multigrain")  # 21
testthat::test_file("tests/testthat/test-transform_pvalues_gsd_opt.R", package = "multigrain") # P1 count
lintr::lint("R/objective_function_gsd_opt.R"); lintr::lint("tests/testthat/test-objective_function_gsd_opt.R")
```
`load_all()` is fine; **no** `install()`/`clean_dll()` (shared library; a measurement agent may
be running). m ≤ 4, nsim ≤ 2000 in tests.

## Report

Files; gate counts; divergences from the brief/plan and why; what you could not verify; what
the design did not anticipate. Do not commit.
