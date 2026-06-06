# print renders the locked cli summary (no abstention)

    Code
      print(cier_longstring(x))
    Output
      -- cier_longstring -------------------------------------------------------------
      Direction: upper -- higher values flag carelessness.
      Cutoff: 3 -- respondents with value >= 3 are flagged.
      Flagged: 3 of 6 scored respondents (50.0%).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

# print reports abstaining respondents on their own line

    Code
      print(cier_longstring(x))
    Output
      -- cier_longstring -------------------------------------------------------------
      Direction: upper -- higher values flag carelessness.
      Cutoff: 3 -- respondents with value >= 3 are flagged.
      Flagged: 3 of 5 scored respondents (60.0%).
      Abstained: 1 (no responses).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

