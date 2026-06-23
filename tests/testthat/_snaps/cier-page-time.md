# print renders the locked cli summary (upper direction)

    Code
      print(cier_page_time(small_pages(), small_ipp()))
    Output
      -- cier_page_time --------------------------------------------------------------
      Direction: upper -- higher values flag carelessness.
      Cutoff: 1 -- respondents with value >= 1 are flagged.
      Cutoff method: fixed count.
      Flagged: 2 of 4 scored respondents (50.0%).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

# print reports abstaining respondents on their own line

    Code
      print(cier_page_time(ps, c(2L, 2L, 2L)))
    Output
      -- cier_page_time --------------------------------------------------------------
      Direction: upper -- higher values flag carelessness.
      Cutoff: 1 -- respondents with value >= 1 are flagged.
      Cutoff method: fixed count.
      Flagged: 2 of 3 scored respondents (66.7%).
      Abstained: 1 (no score).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

