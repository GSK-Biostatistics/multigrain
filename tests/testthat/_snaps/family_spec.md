# family_spec() validates `members`

    Code
      family_spec(1:3)
    Condition
      Error in `family_spec()`:
      ! `members` must be a non-empty list of integer vectors, not an integer vector.

---

    Code
      family_spec(list())
    Condition
      Error in `family_spec()`:
      ! `members` must be a non-empty list of integer vectors, not an empty list.

---

    Code
      family_spec(list(1, 2.5))
    Condition
      Error in `family_spec()`:
      ! `members[[2]]` must be integer-like, not the number 2.5.

---

    Code
      family_spec(list(1, integer(0), 2))
    Condition
      Error in `family_spec()`:
      ! Every family must contain at least one hypothesis. Family 2 is empty.

---

    Code
      family_spec(list(0:1, 2))
    Condition
      Error in `family_spec()`:
      ! The hypotheses in `members` must be positive whole numbers (columns of the p-value matrix).

---

    Code
      family_spec(list(1:2, 2:3))
    Condition
      Error in `family_spec()`:
      ! Each hypothesis must belong to exactly one family. Hypothesis 2 appears in more than one family.

---

    Code
      family_spec(list(1, 3))
    Condition
      Error in `family_spec()`:
      ! `members` must contain every hypothesis from 1 to 3. Hypothesis 2 is not in any family.

---

    Code
      family_spec(list(1, 4))
    Condition
      Error in `family_spec()`:
      ! `members` must contain every hypothesis from 1 to 4. Hypotheses 2 and 3 are not in any family.

---

    Code
      family_spec(list(a = 1, 2))
    Condition
      Error in `family_spec()`:
      ! `members` must be fully named or not named at all.

---

    Code
      family_spec(list(a = 1, a = 2))
    Condition
      Error in `family_spec()`:
      ! Family names in `members` must be unique.

# family_spec() validates `test`

    Code
      family_spec(list(1, 2), test = "bonferroni")
    Condition
      Error in `family_spec()`:
      ! `test` must be one of "hommel", "hochberg", or "holm", not "bonferroni".

---

    Code
      family_spec(list(1, 2, 3), test = c("holm", "hommel"))
    Condition
      Error in `family_spec()`:
      ! `test` must have length 1 or 3 (one per family), not 2.

---

    Code
      family_spec(list(1, 2), test = 1)
    Condition
      Error in `family_spec()`:
      ! `test` must be a character vector, not the number 1.

# check_family_spec()

    Code
      check_family_spec(list(1, 2:3))
    Condition
      Error:
      ! `list(1, 2:3)` must be a multigrain family specification object, not a list.

# check_families_dims()

    Code
      check_families_dims(fams, pvals[, 1:2], num_nodes = 2)
    Condition
      Error:
      ! `pvals` has 2 columns but `families` covers 3 hypotheses.

---

    Code
      check_families_dims(fams, pvals, num_nodes = 3)
    Condition
      Error:
      ! The graph has 3 nodes but `families` has 2 families.
      i With `families`, the graph has one node per family.

# family_spec print and summary methods

    Code
      print(fams)
    Output
      <multigrain_family_spec>
      Hypotheses (columns of the p-value matrix) in each family:
      primary   (hommel): 1
      secondary (hochberg): 2 3 4

---

    Code
      summary(fams)
    Output
      
      Families of hypotheses:
      primary   (hommel): 1
      secondary (hochberg): 2 3 4

