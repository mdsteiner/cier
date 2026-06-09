# print renders the locked cli summary (lower direction)

    Code
      print(cier_ht(m, poly_items()))
    Output
      -- cier_ht ---------------------------------------------------------------------
      Direction: lower -- lower values flag carelessness.
      Cutoff: ~-0.0822 -- respondents with value <= ~-0.0822 are flagged.
      Flagged: 2 of 30 scored respondents (6.7%).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

# print reports abstaining respondents as '(no score)'

    Code
      print(cier_ht(m, poly_items(12L)))
    Output
      -- cier_ht ---------------------------------------------------------------------
      Direction: lower -- lower values flag carelessness.
      Cutoff: ~-0.0736 -- respondents with value <= ~-0.0736 are flagged.
      Flagged: 2 of 29 scored respondents (6.9%).
      Abstained: 1 (no score).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

