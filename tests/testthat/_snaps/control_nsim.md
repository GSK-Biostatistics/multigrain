# control_nsim_local complains with incorrect nsim_local

    Code
      control_nsim_local(ctrl, "foo")
    Condition
      Error in `control_nsim_local()`:
      ! `nsim_local` must be a whole number, not the string "foo".

---

    Code
      control_nsim_local(ctrl, c(100, 200))
    Condition
      Error in `control_nsim_local()`:
      ! `nsim_local` must be a whole number, not a double vector.

# control_nsim_global complains with incorrect nsim_global

    Code
      control_nsim_global(ctrl, "foo")
    Condition
      Error in `control_nsim_global()`:
      ! `nsim_global` must be a whole number, not the string "foo".

---

    Code
      control_nsim_global(ctrl, c(100, 200))
    Condition
      Error in `control_nsim_global()`:
      ! `nsim_global` must be a whole number, not a double vector.

# adjust_nsim_local

    Number of simulations for trial success calculation (`nsim_local`) is greater than the number of rows in `pvals`.
    * Setting `nsim_local` = `nrows(pvals)`.

# adjust_nsim_global

    Number of simulations for trial success calculation within global optimisation algorithm (`nsim_global`) is greater than the number of rows in `pvals`.
    * Setting `nsim_global` = `nrows(pvals)`.

