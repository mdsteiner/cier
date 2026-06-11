# Purpose: Shared pure kernels for the direct family -- attention-check indices
#          that score the respondent's answers against author-declared correct /
#          non-endorsing values, rather than the response pattern. cier_attention
#          is the first (and, in v0.2, only) direct index, so this file opens with
#          it; further direct indices would join here under the single-kernel rule.
# Args:    See per-kernel signatures.
# Returns: Documented per kernel.
# Invariants:
#   - Kernels are pure (no side effects, no I/O, no registry read).
#   - No statistical computation in user-facing files; this file holds the
#     direct-family math consumed by R/cier-attention.R.

# Purpose: Per-respondent count of FAILED attention checks among ANSWERED ones --
#          the Meade & Craig (2012) / Goldammer et al. (2024) direct attention
#          index. A check is failed when the respondent's answered (non-NA)
#          response is NOT a member of that check's pass-set.
# Args:    checks - numeric matrix (n x k); raw check-item responses, one column
#                   per check. NA marks an unanswered check; all cells are finite
#                   or NA (validated upstream by check_responses()).
#          pass   - list length k; pass[[j]] is the passing value(s) of check j,
#                   each a non-empty finite numeric vector (validated upstream by
#                   check_pass()). Positional: pass[[j]] applies to column j.
# Returns: numeric length-n; the count of failed checks per respondent, NA where
#          every check is unanswered (no evidence). value[i] is in [0, k] when the
#          respondent answered at least one check.
# Invariants:
#   - A failure is NON-membership (`!(x %in% pass[[j]])`) among answered cells, so
#     a passing response is never counted. An NA (unanswered) check contributes no
#     evidence: `!is.na(col)` is FALSE there, so `answered & not_in_pass` is FALSE
#     (not NA) and the row sums are plain integer arithmetic with no NA propagation
#     (`NA %in% pass[[j]]` is FALSE because pass[[j]] carries no NA).
#   - The pass-set is applied PER COLUMN (pass[[j]] to column j), never pooled.
#   - A respondent with no answered check (every cell NA) scores NA (abstains).
kernel_attention <- function(checks, pass) {
  k <- ncol(checks)
  answered <- !is.na(checks)                  # materialise the mask once and reuse
  failed <- matrix(FALSE, nrow = nrow(checks), ncol = k)
  for (j in seq_len(k)) {
    failed[, j] <- answered[, j] & !(checks[, j] %in% pass[[j]])
  }
  value <- as.numeric(rowSums(failed))
  value[rowSums(answered) == 0L] <- NA_real_  # every check unanswered -> abstain
  unname(value)
}
