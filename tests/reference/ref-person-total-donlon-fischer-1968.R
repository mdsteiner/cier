# Paper-faithful reference implementation of the person-total correlation
# (r_pbis) index per Donlon & Fischer (1968) / Curran (2016, pp. 12-13).
#
# Definition (item-total form, as in PerFit::r.pbis):
#   1. Compute the per-item mean across all respondents,
#      m[j] = mean over respondents of x[, j] (na.rm).
#   2. For each respondent, take the Pearson correlation between their
#      answered responses and m at those items.
#   3. Return that per-respondent correlation; low values flag careless
#      responding (their pattern does not track the group's item-mean
#      profile).
#
# This is an independent per-row stats::cor() loop; the production kernel
# `kernel_person_total()` is the vectorised masked-sum equivalent. The two
# are cross-checked to 1e-12 in test-cier-person-total.R, so the kernel's
# vectorisation is validated against a transparent reference rather than
# against itself.

ref_person_total_row <- function(row, item_means) {
  ok <- !is.na(row)
  if (sum(ok) < 3L) {
    return(NA_real_)
  }
  suppressWarnings(stats::cor(as.numeric(row[ok]), item_means[ok]))
}

ref_person_total <- function(x) {
  if (!is.matrix(x)) x <- as.matrix(x)
  item_means <- colMeans(x, na.rm = TRUE)
  vapply(seq_len(nrow(x)),
         function(i) ref_person_total_row(x[i, ], item_means),
         numeric(1L))
}
