# sf_hsd_opt() validation messages name the limits

    Code
      sf_hsd_opt(bounds = c(-11, 1))
    Condition
      Error in `sf_hsd_opt()`:
      ! The lower bound in `bounds` must be at least -10, not -11.

---

    Code
      sf_hsd_opt(bounds = c(-4, 2.5))
    Condition
      Error in `sf_hsd_opt()`:
      ! The upper bound in `bounds` must be at most 2, not 2.5.

---

    Code
      sf_hsd_opt(bounds = c(-1, -2))
    Condition
      Error in `sf_hsd_opt()`:
      ! The lower bound in `bounds` must be below the upper bound, not -1 and -2.

---

    Code
      sf_hsd_opt(bounds = c(-1, -1))
    Condition
      Error in `sf_hsd_opt()`:
      ! The lower bound in `bounds` must be below the upper bound, not -1 and -1.

---

    Code
      sf_hsd_opt(bounds = c(-4, Inf))
    Condition
      Error in `sf_hsd_opt()`:
      ! `bounds` must be a finite numeric vector of length 2.

---

    Code
      sf_hsd_opt(bounds = -4)
    Condition
      Error in `sf_hsd_opt()`:
      ! `bounds` must be a finite numeric vector of length 2.

---

    Code
      sf_hsd_opt(step = 0.1)
    Condition
      Error in `sf_hsd_opt()`:
      ! `step` must be one of 0.5, 0.25, or 0.125.
      i Binary fractions keep every lattice node exact.

---

    Code
      sf_hsd_opt(step = c(0.5, 0.25))
    Condition
      Error in `sf_hsd_opt()`:
      ! `step` must be one of 0.5, 0.25, or 0.125.
      i Binary fractions keep every lattice node exact.

---

    Code
      sf_hsd_opt(start = 1.5)
    Condition
      Error in `sf_hsd_opt()`:
      ! `start` must lie inside `bounds` [-4, 1], not 1.5.

---

    Code
      sf_hsd_opt(start = NA_real_)
    Condition
      Error in `sf_hsd_opt()`:
      ! `start` must be `NULL` or a single finite number.

---

    Code
      sf_hsd_opt(id = c("a", "b"))
    Condition
      Error in `sf_hsd_opt()`:
      ! `id` must be `NULL` or a single non-missing string or number.

---

    Code
      sf_hsd_opt(id = NA)
    Condition
      Error in `sf_hsd_opt()`:
      ! `id` must be `NULL` or a single non-missing string or number.

# group members must agree on bounds, step and start

    Code
      transform_pvalues_gsd(raw, spending = list(gsDesign::sfLDOF, sf_hsd_opt(id = "sec",
        bounds = c(-4, -2), step = 0.5), sf_hsd_opt(id = "sec", bounds = c(-4, -1),
      step = 0.5)), grid_size = 16L)
    Condition
      Error in `transform_pvalues_gsd()`:
      ! Hypotheses 2 and 3 share a spending parameter but differ in `bounds` of `sf_hsd_opt()`.
      i Hypotheses with the same `id` must use identical `bounds`, `step` and `start`.

---

    Code
      transform_pvalues_gsd(raw, spending = list(gsDesign::sfLDOF, sf_hsd_opt(id = "sec",
        step = 0.5), sf_hsd_opt(id = "sec", step = 0.25)), grid_size = 16L)
    Condition
      Error in `transform_pvalues_gsd()`:
      ! Hypotheses 2 and 3 share a spending parameter but differ in `step` of `sf_hsd_opt()`.
      i Hypotheses with the same `id` must use identical `bounds`, `step` and `start`.

---

    Code
      transform_pvalues_gsd(raw, spending = list(gsDesign::sfLDOF, sf_hsd_opt(id = "gamma2",
        bounds = c(-4, -2), step = 0.5), sf_hsd_opt(bounds = c(-4, -2), step = 0.5)),
      grid_size = 16L)
    Condition
      Error in `transform_pvalues_gsd()`:
      ! Spending parameter names must be unique.
      x "gamma2" names more than one parameter.
      i Give each `sf_hsd_opt()` group a distinct `id`.

# an all-short-circuited spec is the fixed transform

    Code
      out <- transform_pvalues_gsd(raw, info_frac = info, spending = list(gsDesign::sfLDOF,
      sf_hsd_opt(bounds = c(-4, -2), step = 0.5)), grid_size = 64L)
    Message
      i Every hypothesis tagged with `sf_hsd_opt()` has a single analysis at full information, where spending has no effect.
        Returning a <multigrain_pvals_gsd> object; there is no spending parameter to optimise.

# off-lattice and malformed parameter values are refused

    Code
      resolve_pvalues_gsd(spec, -3.2)
    Condition
      Error in `resolve_pvalues_gsd()`:
      ! `spending_param` value -3.2 for sec is not a node of its lattice.
      i The nearest nodes are -3.5 and -3; the lattice runs from -4 to -2 in steps of 0.5.

---

    Code
      resolve_pvalues_gsd(spec, -5)
    Condition
      Error in `resolve_pvalues_gsd()`:
      ! `spending_param` value -5 for sec is not a node of its lattice.
      i The nearest node is -4; the lattice runs from -4 to -2 in steps of 0.5.

---

    Code
      resolve_pvalues_gsd(spec, c(-3, -2))
    Condition
      Error in `resolve_pvalues_gsd()`:
      ! `spending_param` must be 1 finite number, one per spending parameter ("sec").

---

    Code
      resolve_pvalues_gsd(spec, c(prim = -3))
    Condition
      Error in `resolve_pvalues_gsd()`:
      ! The names of `spending_param` must be the spending parameter names.
      x Expected "sec", got "prim".

---

    Code
      resolve_pvalues_gsd(spec, NA_real_)
    Condition
      Error in `resolve_pvalues_gsd()`:
      ! `spending_param` must be 1 finite number, one per spending parameter ("sec").

---

    Code
      resolve_pvalues_gsd(hsd_opt_raw(), -3)
    Condition
      Error in `resolve_pvalues_gsd()`:
      ! `spec` must be a multigrain_pvals_gsd_spec object (see sf_hsd_opt()), not a double 3D array.

# a node that is not well ordered aborts naming it

    Code
      transform_pvalues_gsd(hsd_opt_raw(), spending = list(gsDesign::sfLDOF,
      gsDesign::sfLDOF, not_well_ordered), grid_size = 64L)
    Condition
      Error in `transform_pvalues_gsd()`:
      ! Hypothesis 3 is not well ordered at sec = -2.
      i Narrow `bounds` of `sf_hsd_opt()`.
      Caused by error in `transform_pvalues_gsd()`:
      ! Spending function for hypothesis 3 is not well ordered.
      x Its nominal boundary at analysis 1 decreases as the allocated level increases.
      i The group sequential graphical procedure is not valid for such a spending function.

# a table without two distinct positive boundaries is refused

    Code
      .hsd_check_invertible(tab, hyp = 2L, node = -3, name = "sec")
    Condition
      Error:
      ! The boundary table of hypothesis 2 cannot be inverted at sec = -3.
      x Fewer than two distinct positive boundaries were found at analysis 2.
      i Narrow `bounds` of `sf_hsd_opt()`.

# check_pvals_gsd() refuses a spec with a hint

    Code
      check_pvals_gsd(spec)
    Condition
      Error:
      ! `spec` is a spending spec, not a transformed p-value object.
      i Pass it to `graph_optimise_gsd()` to optimise the spending parameters, or call `resolve_pvalues_gsd()` with chosen values.

# the spec prints

    Code
      print(hsd_opt_spec())
    Output
      <multigrain_pvals_gsd_spec>
      2000 simulation(s), 3 hypothes(es), 2 analys(es); alpha = 0.025
      Information fractions (NA = no data at that analysis):
         analysis 1 analysis 2
      H1        0.5          1
      H2        0.5          1
      H3        0.5          1
      Look-back: H1 = FALSE, H2 = FALSE, H3 = FALSE
      Spending:
      H1: gsDesign::sfLDOF
      H2: sfHSD(param = sec in [-4, -2])
      H3: sfHSD(param = sec in [-4, -2])
      Spending parameters:
      sec: bounds [-4, -2], step 0.5, 5 nodes, start -3
      Boundary tables built: 5

---

    Code
      print(hsd_opt_spec_d2())
    Output
      <multigrain_pvals_gsd_spec>
      1000 simulation(s), 4 hypothes(es), 2 analys(es); alpha = 0.025
      Information fractions (NA = no data at that analysis):
         analysis 1 analysis 2
      H1        0.5          1
      H2        0.5          1
      H3        0.5          1
      H4         NA          1
      Look-back: H1 = FALSE, H2 = FALSE, H3 = FALSE, H4 = FALSE
      Spending:
      H1: gsDesign::sfLDOF
      H2: sfHSD(param = gamma1 in [-3, -2])
      H3: sfHSD(param = b in [-3, -2])
      H4: sfHSD(param = -2.5) [short-circuited: one look at full information]
      Spending parameters:
      gamma1: bounds [-3, -2], step 0.5, 3 nodes, start -2.5
      b: bounds [-3, -2], step 0.5, 3 nodes, start -2.5
      Boundary tables built: 3

