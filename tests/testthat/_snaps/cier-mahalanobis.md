# print renders the locked cli summary (upper direction)

    Code
      print(cier_mahalanobis(x, alpha = 0.2))
    Output
      -- cier_mahalanobis ------------------------------------------------------------
      Direction: upper -- higher values flag carelessness.
      Cutoff: ~8.56 -- respondents with value >= ~8.56 are flagged.
      Flagged: 3 of 30 scored respondents (10.0%).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

# print reports abstaining respondents on their own line

    Code
      print(cier_mahalanobis(x, alpha = 0.2))
    Output
      -- cier_mahalanobis ------------------------------------------------------------
      Direction: upper -- higher values flag carelessness.
      Cutoff: ~8.56 -- respondents with value >= ~8.56 are flagged.
      Flagged: 5 of 29 scored respondents (17.2%).
      Abstained: 1 (no score).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

