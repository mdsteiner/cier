# print renders the locked cli summary (deterministic PR)

    Code
      print(out)
    Output
      -- cier_personal_reliability ---------------------------------------------------
      Direction: upper -- higher values flag carelessness.
      Cutoff: 1 -- respondents with value >= 1 are flagged.
      Flagged: 12 of 23 scored respondents (52.2%).
      Abstained: 1 (no score).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

# print reports abstaining respondents on their own line (PR)

    Code
      print(out)
    Output
      -- cier_personal_reliability ---------------------------------------------------
      Direction: upper -- higher values flag carelessness.
      Cutoff: 1 -- respondents with value >= 1 are flagged.
      Flagged: 12 of 23 scored respondents (52.2%).
      Abstained: 2 (no score).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

