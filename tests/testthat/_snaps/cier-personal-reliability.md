# print renders the locked cli summary (deterministic PR)

    Code
      print(cier_personal_reliability(x, it, resample = FALSE))
    Output
      -- cier_personal_reliability ---------------------------------------------------
      Direction: upper -- higher values flag carelessness.
      Cutoff: 1 -- respondents with value >= 1 are flagged.
      Flagged: 3 of 11 scored respondents (27.3%).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

# print reports abstaining respondents on their own line (PR)

    Code
      print(cier_personal_reliability(x, it, resample = FALSE))
    Output
      -- cier_personal_reliability ---------------------------------------------------
      Direction: upper -- higher values flag carelessness.
      Cutoff: 1 -- respondents with value >= 1 are flagged.
      Flagged: 3 of 11 scored respondents (27.3%).
      Abstained: 1 (no score).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

