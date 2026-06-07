# print renders the locked cli summary (lower direction)

    Code
      print(cier_irv(x))
    Output
      -- cier_irv --------------------------------------------------------------------
      Direction: lower -- lower values flag carelessness.
      Cutoff: ~0.781 -- respondents with value <= ~0.781 are flagged.
      Flagged: 2 of 30 scored respondents (6.7%).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

# print reports abstaining respondents on their own line

    Code
      print(cier_irv(x))
    Output
      -- cier_irv --------------------------------------------------------------------
      Direction: lower -- lower values flag carelessness.
      Cutoff: ~0.871 -- respondents with value <= ~0.871 are flagged.
      Flagged: 2 of 29 scored respondents (6.9%).
      Abstained: 1 (no responses).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

