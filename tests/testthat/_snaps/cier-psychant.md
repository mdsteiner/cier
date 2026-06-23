# print renders the locked cli summary (upper direction)

    Code
      print(cier_psychant(ant_matrix(n = 30L, seed = 11L)))
    Output
      -- cier_psychant ---------------------------------------------------------------
      Direction: upper -- higher values flag carelessness.
      Cutoff: ~-0.583 -- respondents with value >= ~-0.583 are flagged.
      Cutoff method: 95th sample percentile (fpr = 0.05).
      Flagged: 2 of 30 scored respondents (6.7%).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

# print reports abstaining respondents as '(no score)'

    Code
      print(cier_psychant(x))
    Output
      -- cier_psychant ---------------------------------------------------------------
      Direction: upper -- higher values flag carelessness.
      Cutoff: ~-0.548 -- respondents with value >= ~-0.548 are flagged.
      Cutoff method: 95th sample percentile (fpr = 0.05).
      Flagged: 2 of 29 scored respondents (6.9%).
      Abstained: 1 (no score).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

