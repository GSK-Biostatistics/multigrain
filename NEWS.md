# multigrain (development version)

# multigrain 0.3.0

## New functionality

* `graph_optimal_get_control()` can be used to extract the optimisation settings (i.e. the `multigrain_control` object) from the optimised graph.
* `graph_optimise()` (named changed - see "Changes" below) now accepts a `num_threads` argument for parallel execution of the shortcut algorithm. This replaces the previous `control_parallel()` workflow. Default is `1L` (serial).
* Users can supply names to `graph_random()`.
* Users can supply titles when plotting `multigrain_graph_constraint` or `multigrain_graph_optimal` objects.

## Bug fixes

* `trial_success()` now accepts expressions/strings containing only `r<digit>` symbols, numeric literals, and operators; this fixes failures when an object named `r<digit>` exists in the global environment.

## Changes

* `optimise_graph()` and `optimize_graph()` have been renamed to `graph_optimise()` and `graph_optimize()`
* the `verbose` argument to `graph_optimise()` and `trial_success()` is no longer a logical, but a character with the user being able to choose one of three verbosity levels:
  * `"info"`: will only show milestones / informational messages highlighting the progress of the optimisation at coarse-grained level. This corresponds to what was previously `TRUE`.
  * `"detail"`: will show milestones and information about fine-grained optimisation events. All local and global optimisation messages are displayed in the console. This is a newly introduced level.
  * `"silent"`: no information about the progress of the optimisation is printed to the console. This corresponds to what was previously `FALSE`.
  * verbosity can be control both at _function-_ and _package-_ level. At _package-level_ this is done with the `multigrain_verbosity` option which should be set to one of the three possible values ("detail" > "info" > "silent").
* in calls to the following functions all optional arguments must be named (they can no longer be passed by position): `calc_power_pvals()`, `graph_constraint()`, `graph_optimise()`, `simulate_pvalues()`, and `normalise_sum()`. These functions also gain an `...` (also know as the "ellipsis") argument to allow for future extensions.
* `random_graph()` has been renamed to `graph_random()` (for consistency with `graph_optimise()` and to avoid a conflict with `graphicalMCP::random_graph()`).
* Updates to the `multigrain_control` object's print method:
    * `optimArgs` settings for the global optimisation are now printed, not just the top level class.
* Changes to the `multigrain_graph_optimal` object:
    * the `GA_output` element has been renamed to `global_output`.
    * the `nloptr_output` element has been renamed to `local_output`.
    * the `opt_settings$global_search` element is replaced by `global_search`.
* Changes to parallelisation:
    * `optimise_graph_parallel()` has been removed.
    * use `graph_optimise()` with `num_threads` instead.
    * `control_parallel()` and the `parallel_opt` slot on `multigrain_control` are removed.
    * `grain_size` is no longer user-facing; it is always auto-tuned internally.
* `global_output` and `local_output` now store the raw global (`GA`) / local (`nloptr`) objects directly.
* `control_global_search()` has been removed. Reverted to `global_search` being a direct argument to `graph_optimise()`.
* the names of the multigrain S3 classes are now prefixed with `"multigrain_":
    * `graph_constraint` -> `multigrain_graph_constraint`.
    * `graph_optimal` -> `multigrain_graph_optimal`.
    * `trial_success` -> `multigrain_trial_success`.
    * no change for _control_, which has always been `multigrain_control`.
    * the `summary()` and `print()` methods for these objects have been updated.
* `graph_optimise()`'s renormalises hypothesis-weight vectors and transition-matrix rows
via `normalise_sum()`. Now each call will use a `tolerance` passed from `graph_constraint`.
If not supplied, it will fall back to `tolerance = sqrt(.Machine$double.eps)`


# multigrain 0.2.0

## New functionality

* Advanced optimisation configuration is now handled via a new `multigrain_control` object and a family of `control_*()` functions. The main optimisation functions (`optimise_graph()`, `optimise_graph_parallel()`) take a single control argument, providing a unified and extensible API for setting optimisation parameters.
* Dedicated function to find closest graph conditional to constraints (`closest_graph_to_constraints()`).
* Added `optimise_graph_parallel()`, which has identical behaviour to `optimise_graph()` but exploits multithreaded parallel execution of graphical test using RcppParallel package.
  * Users can set parallelisation options (`num_threads`, `grain_size`), with `multigrain_control()` + `control_parallel()`.
  * Default behaviour is to use a single thread.
* Both `optimise_graph()` and `optimise_graph_parallel()` now use a memory efficient implementation of the Bretz shortcut algorithm in the internal functions `graph_shortcut()` and `graph_shortcut_parallel()`. The full CTP is still available in `apply_ctp()` - this will be used for Simes and parametric graph optimisation in future releases.
* Users can supply names for unconstrained graphs via the `names` argument to `graph_constraint_free().

## Bug fixes

* `$power$trial_success` from `graph_optimal` object did not match the graph the user sees (eg when plotting), due to power being calculated before pruning. `multigrain` now stores pre-pruned powers as `$global_opt_power` / `$local_opt_power`, and set `$power` to the post-pruned evaluation.
* Improved argument checking and error messages for optimisation parameters.

## Changes

* The sample size optimisation functionality (`optimise_N()`) has been removed.
* Improvements to print and summary methods for key S3 classes (`graph_optimal`, `graph_constraint`, and `trial_success`).


# multigrain 0.1.1 - 2025-10-02

## New functionality

## Bug fixes

* Fixed bug that causes 2-hypothesis optimisation of graphs using `optimise_graphs()` to fail.
* Fixed the behaviour of `graph_constraint_free(2)`.

## Changes

* Fixed print/summary S3 methods for `graph_optimal` to refer to "trial success" instead of "power objective".

# multigrain 0.1.0 - 2025-09-20

## New functionality

* Sample size optimisation through `optimise_N()`
  * Minimises sample size, finding the optimal graph to do this with
  trial success / per-hypothesis marginal power targets.
* Added plotting methods for `graph_optimal` and `graph_constraint` objects.
  * Both can be plotted with `autoplot()` and `plot()` functions.
* Overhauled the `graph_constraint` functionality:
  * a `graph_constraint` is created with `graph_constraint()`
  * added built-in validation
  * diagnosis is optional and can be accessed by running
  `graph_constraint(..., diagnose = TRUE)`
  * create an unconstrained graph constraint with `graph_constraint_free()`.
  * users can get and set `graph_constraint` elements with `[` and `$`. A
  modified `graph_constraint` is then automatically (re-)validated.
  * `graph_constraint` validation happens with tolerance
* pkgwdown site and vignette
* `simulate_pvalues()` function simulates raw p-values from under the
alternative hypotheses and the assumption that the distribution of test
statistics is a multivariate normal distribution.
* Add graph validation functions

## Bug fixes
- Bug where regex translation of trial_success function does not create
correct code.
- Internal function `recover_G()` did not recover transition matrices from
genetic chromosome every time,
this has been fixed with its replacement `recover_full_trans_matrix()`

## Changes
- `calc_power()` function created based on `eval_power()` (the latter is deleted)
- Many functions, classes and arguments have been renamed. The most significant are:
  - The `graphOpt` S3 class is now named the `graph_optimal`.
  - Power objectives are referred to as "trial success measures"; `customPowerObj()` is renamed
  to `trial_success()` to reflect this.
  - Functions and function arguments have been changed to have snake-case convention, e.g.,
    - `graphOpt()` function is now `optimise_graph()`
    - `graphConstraint()` function is now `graph_constraint()`, etc.
  - Checks of `pvals`/`graph_constraint`/`start_graph` dimensions and arguments have been
  introduced for `optimise_graph()` to prevent downstream crashing.


# multigrain 0.0.2 - 2024-10-21

## New functionality
- GA hybrid-global optimisation route
- Documentation for core functions
- print and summary S3 methods for `graphOpt` class
- Unit tests for all core functions
- Proprietary license
- Changed function name `designConstraint` to `graphConstraint`
- Added m=2 option for optimisation
- `is_graph_valid()` function to check graph validity
- `normalise_sum()` (replaced `normalise_exactly()`) to fix floating point
errors causing w and G-rows to not sum to 1. Now uses tolerance-based
convergence and respects fixed-element constraints.


## Changes
- Confirmed R package dependencies (GA, nloptr, gMCPLite, etc.)
- Updated the version in the DESCRIPTION file from 0.0.1 to 0.0.2
- Removed inequality constraints from NLOPT routine

## Bug fixes
- Bug where local/global power comparison at end of `optimise_graph()` would not return
correct graph (should be highest-power valid graph)
