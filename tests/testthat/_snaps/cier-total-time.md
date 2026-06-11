# print renders the locked cli summary (lower direction)

    Code
      print(cier_total_time(spread_seconds()))
    Output
      -- cier_total_time -------------------------------------------------------------
      Direction: lower -- lower values flag carelessness.
      Cutoff: ~34.1 -- respondents with value <= ~34.1 are flagged.
      Flagged: 2 of 30 scored respondents (6.7%).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

# print reports abstaining respondents on their own line

    Code
      print(cier_total_time(x))
    Output
      -- cier_total_time -------------------------------------------------------------
      Direction: lower -- lower values flag carelessness.
      Cutoff: ~34 -- respondents with value <= ~34 are flagged.
      Flagged: 2 of 29 scored respondents (6.9%).
      Abstained: 1 (no score).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

