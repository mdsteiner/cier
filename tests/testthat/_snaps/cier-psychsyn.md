# print renders the locked cli summary (lower direction)

    Code
      print(cier_psychsyn(syn_matrix(n = 30L, seed = 11L)))
    Output
      -- cier_psychsyn ---------------------------------------------------------------
      Direction: lower -- lower values flag carelessness.
      Cutoff: ~-0.239 -- respondents with value <= ~-0.239 are flagged.
      Flagged: 2 of 30 scored respondents (6.7%).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

# print reports abstaining respondents as '(no score)'

    Code
      print(cier_psychsyn(x))
    Output
      -- cier_psychsyn ---------------------------------------------------------------
      Direction: lower -- lower values flag carelessness.
      Cutoff: ~-0.132 -- respondents with value <= ~-0.132 are flagged.
      Flagged: 2 of 29 scored respondents (6.9%).
      Abstained: 1 (no score).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

