# print renders the locked cli summary (upper direction)

    Code
      print(cier_autocorrelation(ac_fixture(n = 30L, p = 12L), max_lag = 5L))
    Output
      -- cier_autocorrelation --------------------------------------------------------
      Direction: upper -- higher values flag carelessness.
      Cutoff: ~0.725 -- respondents with value >= ~0.725 are flagged.
      Flagged: 2 of 30 scored respondents (6.7%).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

# print reports abstaining respondents on their own line

    Code
      print(cier_autocorrelation(x, max_lag = 5L))
    Output
      -- cier_autocorrelation --------------------------------------------------------
      Direction: upper -- higher values flag carelessness.
      Cutoff: ~0.798 -- respondents with value >= ~0.798 are flagged.
      Flagged: 2 of 29 scored respondents (6.9%).
      Abstained: 1 (no score).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

