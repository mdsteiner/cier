# Independent reference for cier_attention() (direct family; Meade & Craig 2012,
# <doi:10.1037/a0028085>; Goldammer et al. 2024).
#
# attention is a counting rule, not an estimator: per respondent, count the
# attention checks that were FAILED among those ANSWERED. A check is failed when
# its answered (non-NA) response is NOT a member of that check's pass-set. An NA
# (unanswered) check contributes no evidence -- it is neither passed nor failed --
# and a respondent who answered NO check at all abstains (NA).
#
# These helpers re-derive each quantity by an independent path -- a per-respondent
# ROW loop with %in%, never the production kernel's column-wise vectorised
# membership matrix -- so any divergence is attributable to the wrapper. attention
# has no CRAN parity partner (no package implements this attention-check counting
# rule as a C/IER index), so this oracle is the sole parity
# check -- oracle-only trust, like page_time / total_time / PR / RPR. Same base-R
# membership arithmetic as the wrapper, so the target tolerance is 0 (exact integer
# counts).

# Per-respondent failed-check count. `pass` is a positional list of length
# ncol(checks); pass[[j]] gives the passing value(s) of check j. For respondent i
# the count is the number of columns j where checks[i, j] is observed (non-NA) and
# checks[i, j] is NOT in pass[[j]]. A respondent with every check NA scores NA
# (abstains). A mutant that counts an NA check as a failure, counts a missing row
# as 0 instead of abstaining, pools the pass-sets across columns, applies pass to
# the wrong column, or inverts the membership (counts passes) diverges from this.
ref_attention <- function(checks, pass) {
  checks <- as.matrix(checks)
  n <- nrow(checks)
  vapply(seq_len(n), function(i) {
    row <- checks[i, ]
    answered <- !is.na(row)
    if (!any(answered)) {
      return(NA_real_)
    }
    failed <- 0
    for (j in seq_along(row)) {
      if (answered[[j]] && !(row[[j]] %in% pass[[j]])) {
        failed <- failed + 1
      }
    }
    as.numeric(failed)
  }, numeric(1L))
}

# Upper-tail flag rule: a respondent is flagged iff their failed-check count is at
# or above the cutoff (ties flag). An NA cutoff flags nobody; an abstaining (NA
# value) respondent is never flagged -- and, mirroring the production
# new_cier_index() one-abstention rule, their flag is NA, not FALSE, so this oracle
# stays valid if reused on a fixture that contains abstaining rows.
ref_attention_flags <- function(value, cutoff) {
  if (is.na(cutoff)) {
    return(rep(FALSE, length(value)))
  }
  flagged <- value >= cutoff
  flagged[is.na(value)] <- NA
  flagged
}
