# print renders the locked cli summary (upper direction)

    Code
      print(cier_attention(small_checks(), small_pass()))
    Output
      -- cier_attention --------------------------------------------------------------
      Direction: upper -- higher values flag carelessness.
      Cutoff: 1 -- respondents with value >= 1 are flagged.
      Cutoff method: fixed count.
      Flagged: 3 of 4 scored respondents (75.0%).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

# print reports abstaining respondents on their own line

    Code
      print(cier_attention(ch, list(c(1, 2), 0)))
    Output
      -- cier_attention --------------------------------------------------------------
      Direction: upper -- higher values flag carelessness.
      Cutoff: 1 -- respondents with value >= 1 are flagged.
      Cutoff method: fixed count.
      Flagged: 1 of 2 scored respondents (50.0%).
      Abstained: 1 (no score).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

