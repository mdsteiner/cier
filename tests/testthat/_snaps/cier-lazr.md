# print renders the locked cli summary (upper direction)

    Code
      print(cier_lazr(lazr_fixture(n = 30L, p = 12L)))
    Output
      -- cier_lazr -------------------------------------------------------------------
      Direction: upper -- higher values flag carelessness.
      Cutoff: ~0.652 -- respondents with value >= ~0.652 are flagged.
      Cutoff method: 95th sample percentile (fpr = 0.05).
      Flagged: 3 of 30 scored respondents (10.0%).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

# print reports abstaining respondents on their own line

    Code
      print(cier_lazr(x))
    Output
      -- cier_lazr -------------------------------------------------------------------
      Direction: upper -- higher values flag carelessness.
      Cutoff: ~0.606 -- respondents with value >= ~0.606 are flagged.
      Cutoff method: 95th sample percentile (fpr = 0.05).
      Flagged: 5 of 29 scored respondents (17.2%).
      Abstained: 1 (no score).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

