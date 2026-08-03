# create_obj_func: 3-m disjunctive power: improvedFallbackI

    Code
      cat("3-m disjunctive power: improvedFallbackI")
    Output
      3-m disjunctive power: improvedFallbackI
    Code
      cat(paste("multigrain power:", round(test_power, 4)))
    Output
      multigrain power: 0.9347
    Code
      cat(paste("gMCP power:", round(power_disj$disj_power, 4)))
    Output
      gMCP power: 0.9347

# create_obj_func: 4-m conjunctive pow: improvedParallelGatekeeping

    Code
      cat("4-m conjunctive power: improvedParallelGatekeeping")
    Output
      4-m conjunctive power: improvedParallelGatekeeping
    Code
      cat(paste("multigrain power:", round(test_power, 4)))
    Output
      multigrain power: 0.3979
    Code
      cat(paste("gMCP power:", round(power_conj$conj_power, 4)))
    Output
      gMCP power: 0.3979

# create_obj_func: 6-m average power: BretzEtAl2011

    Code
      cat("6-m average power: BretzEtAl2011")
    Output
      6-m average power: BretzEtAl2011
    Code
      cat(paste("multigrain power:", round(test_power, 4)))
    Output
      multigrain power: 3.5351
    Code
      cat(paste("gMCP power:", round(power_avg$exp_rejections, 4)))
    Output
      gMCP power: 3.5351

# create_obj_func handles near-zero transition weights

    Code
      cat(paste("multigrain power:", round(package_power, 4)))
    Output
      multigrain power: 0.3603
    Code
      cat(paste("gMCP power:", round(power_nearzero$conj_power, 4)))
    Output
      gMCP power: 0.3603

---

    Code
      cat(paste("multigrain power:", round(package_power, 4)))
    Output
      multigrain power: 0.3978
    Code
      cat(paste("gMCP power:", round(power_noadjust$conj_power, 4)))
    Output
      gMCP power: 0.3978

