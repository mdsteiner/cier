# print renders the locked cli summary (upper direction)

    Code
      print(cier_gnormed(m, poly_items(), cutoff = 0.5))
    Output
      -- cier_gnormed ----------------------------------------------------------------
      Direction: upper -- higher values flag carelessness.
      Cutoff: 0.5 -- respondents with value >= 0.5 are flagged.
      Flagged: 2 of 30 scored respondents (6.7%).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

# print reports abstaining respondents as '(no score)'

    Code
      print(cier_gnormed(m, poly_items(12L), cutoff = 0.5))
    Output
      -- cier_gnormed ----------------------------------------------------------------
      Direction: upper -- higher values flag carelessness.
      Cutoff: 0.5 -- respondents with value >= 0.5 are flagged.
      Flagged: 2 of 29 scored respondents (6.9%).
      Abstained: 1 (no score).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

