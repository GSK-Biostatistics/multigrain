# powerFunc matches numeric and logical matrix inputs

    Code
      ts <- new_trial_success(e)
    Message
      v Trial success function compiled and sourced successfully.

---

    Code
      ts <- new_trial_success(e)
    Message
      v Trial success function compiled and sourced successfully.

---

    Code
      ts <- new_trial_success(e)
    Message
      v Trial success function compiled and sourced successfully.

---

    Code
      ts <- new_trial_success(e)
    Message
      v Trial success function compiled and sourced successfully.

---

    Code
      ts <- new_trial_success(e)
    Message
      v Trial success function compiled and sourced successfully.

---

    Code
      ts <- new_trial_success(e)
    Message
      v Trial success function compiled and sourced successfully.

# trial_success print and summary methods

    Code
      print(ts)
    Output
      <multigrain_trial_success>
      r1 + r2 + r3 + r4

---

    Code
      summary(ts)
    Output
      
      Trial success function:
      r1 + r2 + r3 + r4

# trial_success chatty

    Code
      ts <- trial_success(r1 + r2 + r3 + r4, verbose = TRUE)
    Message
      v Trial success function compiled and sourced successfully.

---

    Code
      ts <- trial_success(r1 + r2 + r3 + r4, verbose = FALSE)

---

    Code
      trial_success(r1 + r2 + r3 + r4, verbose = 2)
    Condition
      Error in `trial_success()`:
      ! `verbose` must be a logical vector, not the number 2.

