# print renders the locked cli summary (lower direction)

    Code
      print(cier_person_total(rand_matrix(n = 30L, p = 12L, seed = 11L)))
    Output
      -- cier_person_total -----------------------------------------------------------
      Direction: lower -- lower values flag carelessness.
      Cutoff: ~-0.183 -- respondents with value <= ~-0.183 are flagged.
      Cutoff method: 5th sample percentile (fpr = 0.05).
      Flagged: 2 of 30 scored respondents (6.7%).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

# print reports abstaining respondents on their own line

    Code
      print(cier_person_total(x))
    Output
      -- cier_person_total -----------------------------------------------------------
      Direction: lower -- lower values flag carelessness.
      Cutoff: ~-0.284 -- respondents with value <= ~-0.284 are flagged.
      Cutoff method: 5th sample percentile (fpr = 0.05).
      Flagged: 2 of 29 scored respondents (6.9%).
      Abstained: 1 (no score).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

