# print renders the locked cli summary with the Cutoff method line

    Code
      print(cier_irv(prov_matrix()))
    Output
      -- cier_irv --------------------------------------------------------------------
      Direction: lower -- lower values flag carelessness.
      Cutoff: ~0.579 -- respondents with value <= ~0.579 are flagged.
      Cutoff method: 5th sample percentile (fpr = 0.05).
      Flagged: 2 of 30 scored respondents (6.7%).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

# summary.cier_index renders the locked snapshot (scored + abstain)

    Code
      summary(cier_irv(prov_matrix()))
    Output
      -- cier_irv (summary) ----------------------------------------------------------
      Direction: lower -- lower values flag carelessness.
      Cutoff: ~0.579 (5th sample percentile; fpr = 0.05).
      Scored 30 of 30 respondents (0 abstained).
      Score quartiles: min 0.426 | Q1 0.899 | median 1.05 | Q3 1.19 | max 1.46.
      At or below the cutoff: 2 of 30 scored (6.7%).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

---

    Code
      summary(suppressWarnings(cier_irv(x)))
    Output
      -- cier_irv (summary) ----------------------------------------------------------
      Direction: lower -- lower values flag carelessness.
      Cutoff: ~0.871 (5th sample percentile; fpr = 0.05).
      Scored 29 of 30 respondents (1 abstained).
      Score quartiles: min 0.632 | Q1 1.17 | median 1.47 | Q3 1.72 | max 2.19.
      At or below the cutoff: 2 of 29 scored (6.9%).
      i Per-respondent scores in `$value`, flags in `$flagged`.
      

