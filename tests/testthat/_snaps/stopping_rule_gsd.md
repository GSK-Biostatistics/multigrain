# print and summary methods

    Code
      print(recycled)
    Output
      <multigrain_stopping_rule_gsd>
      Every interim: r1 && (r2 || r3)
      Hypotheses referenced: up to r3

---

    Code
      summary(recycled)
    Output
      
      Stopping rule (group sequential):
      Every interim: r1 && (r2 || r3)
      Hypotheses referenced: up to r3

---

    Code
      print(two)
    Output
      <multigrain_stopping_rule_gsd>
      Interim 1: r1
      Interim 2: r1 && (r2 || r3)
      Hypotheses referenced: up to r3

---

    Code
      summary(two)
    Output
      
      Stopping rule (group sequential):
      Interim 1: r1
      Interim 2: r1 && (r2 || r3)
      Hypotheses referenced: up to r3

