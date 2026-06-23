# Shared pure kernels for the direct family: attention-check indices scoring answers
# against author-declared pass values.

# Per-respondent count of failed attention checks among answered ones (Meade & Craig
# 2012; Goldammer et al. 2024). A check fails when its answered (non-NA) response is not
# in that column's pass-set (pass[[j]] to column j). Returns numeric length-n in [0, k];
# a respondent with no answered check abstains (NA). Unanswered cells give no evidence:
# the mask makes `answered & not_in_pass` FALSE, not NA.
kernel_attention <- function(checks, pass) {
  k <- ncol(checks)
  answered <- !is.na(checks)
  failed <- matrix(FALSE, nrow = nrow(checks), ncol = k)
  for (j in seq_len(k)) {
    failed[, j] <- answered[, j] & !(checks[, j] %in% pass[[j]])
  }
  value <- as.numeric(rowSums(failed))
  value[rowSums(answered) == 0L] <- NA_real_  # every check unanswered -> abstain
  unname(value)
}
