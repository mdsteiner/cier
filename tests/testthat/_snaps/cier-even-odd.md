# print renders the locked cli summary (upper direction)

    Code
      print(cier_even_odd(rand_matrix(11L, 12L, 11L), it))
    Output
      -- cier_even_odd ---------------------------------------------------------------
      Direction: upper -- higher values flag carelessness.
      Cutoff: 1 -- respondents with value >= 1 are flagged.
      Flagged: 3 of 11 scored respondents (27.3%).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

# print reports abstaining respondents on their own line

    Code
      print(cier_even_odd(x, it))
    Output
      -- cier_even_odd ---------------------------------------------------------------
      Direction: upper -- higher values flag carelessness.
      Cutoff: 1 -- respondents with value >= 1 are flagged.
      Flagged: 3 of 11 scored respondents (27.3%).
      Abstained: 1 (no score).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

