# bank slices, byte count and the memory guard

    Code
      .check_bank_memory(2.5e+09)
    Condition
      Warning:
      The spending-parameter slice bank needs 2.5 GB.
      i Use a coarser `step` or narrower `bounds` in `sf_hsd_opt()`, or a smaller `nsim_global` (see `control_nsim_global()`).

---

    Code
      .check_bank_memory(4.5e+09)
    Condition
      Error:
      ! The spending-parameter slice bank would need 4.5 GB, above the limit of 4 GB.
      i Use a coarser `step` or narrower `bounds` in `sf_hsd_opt()`, or a smaller `nsim_global` (see `control_nsim_global()`).

---

    Code
      .check_bank_memory(1.5e+09, max_bank_gb = 1)
    Condition
      Error:
      ! The spending-parameter slice bank would need 1.5 GB, above the limit of 1 GB.
      i Use a coarser `step` or narrower `bounds` in `sf_hsd_opt()`, or a smaller `nsim_global` (see `control_nsim_global()`).

