# Purpose: Low-level numerical kernels shared by the indirect-index public
#          wrappers (R/cier-longstring.R, ...). Each kernel takes a numeric
#          matrix of responses and returns the per-respondent scores it is
#          responsible for. Single-kernel rule: one production implementation
#          per statistic.
# Args:    See per-kernel documentation below.
# Returns: Numeric vectors; never raises typed errors (the wrappers validate).
# Invariants:
#   - Kernels are pure (no I/O, no global state) and never mutate inputs.

# ---- Longstring -------------------------------------------------------------

# Maximum run length of consecutive identical responses per respondent, over
# the *raw* row (no scale blocking). Bytewise compatible with
# careless::longstring(): base::rle() treats each NA as its own run (NA == NA
# is NA, not TRUE), so identical values separated by NA do not merge and an
# all-NA row yields a max run length of 1. The wrapper applies NA-abstention for
# rows with no present responses; the kernel itself stays pure.
kernel_longstring <- function(responses) {
  vapply(
    seq_len(nrow(responses)),
    function(i) max(rle(responses[i, ])$lengths),
    numeric(1L)
  )
}

# ---- IRV (intra-individual response variability) ----------------------------

# Per-respondent SAMPLE standard deviation (denominator n - 1) across the items
# a respondent answered. `matrixStats::rowSds()` is the C-vectorised equivalent
# of the original `apply(x, 1, stats::sd, na.rm = TRUE)`: it shares stats::sd's
# NA semantics under `na.rm = TRUE`, so a row with fewer than two present values
# yields NA (never NaN) and abstention falls out of the summary itself with no
# extra guard (unlike longstring, whose all-NA row scores 1). Its summation
# order differs from a hand-rolled two-pass formula at ulp level, so the
# reference oracle and careless::irv() are held to 1e-10, not bytewise (see
# tests/reference/TOLERANCES.md).
kernel_irv <- function(responses) {
  matrixStats::rowSds(responses, na.rm = TRUE)
}

# ---- Mahalanobis ------------------------------------------------------------

# Per-respondent SQUARED Mahalanobis distance D^2 of each response vector from
# the column-mean centroid, using pairwise covariance. Numerically equivalent to
# careless::mahad() and psych::outlier(plot = FALSE, bad = 0, na.rm = TRUE) --
# both of which return D^2 -- at 1e-10 (see tests/reference/TOLERANCES.md):
#   * sx <- stats::cov(x, use = "pairwise"); sx_inv <- solve(sx).
#   * Centre by column means via scale(scale = FALSE), then zero-fill the
#     remaining NAs in the CENTRED matrix. Zeros contribute nothing to the
#     bilinear sum, so the missing terms drop exactly as psych::outlier's
#     na.rm = TRUE skips them. This is NOT na.rm on the product: one NA there
#     turns a whole sx_inv %*% xc product NA and collapses the row to 0. The
#     zero-fill lets one BLAS matrix multiply replace the per-row apply() loop:
#     D^2 = diag(X sx_inv X^T) = rowSums((X sx_inv) * X).
# Returns list(value, status), status a scalar:
#   "ok"                     - D^2 computed (a row is NA only when it is all-NA).
#   "insufficient_responses" - fewer than two respondents carry any data.
#   "singular_covariance"    - the pairwise covariance is not invertible: solve()
#                              failed, OR it returned a non-finite inverse. The
#                              latter is the pairwise-specific case where a cov
#                              cell is NA (two items never co-answered, or an item
#                              answered by nobody): solve() does NOT error on an
#                              NA matrix, it returns an all-NA inverse, so an
#                              `anyNA()` check is required on top of tryCatch.
# On the two non-"ok" statuses every value is NA; the wrapper raises the typed
# warning. The kernel stays pure (it raises no conditions).
kernel_mahalanobis <- function(responses) {
  value <- rep(NA_real_, nrow(responses))
  keep <- rowSums(!is.na(responses)) > 0L
  if (sum(keep) < 2L) {
    return(list(value = value, status = "insufficient_responses"))
  }
  x_f <- responses[keep, , drop = FALSE]
  sx_inv <- tryCatch(solve(stats::cov(x_f, use = "pairwise")),
                     error = function(e) NULL)
  if (is.null(sx_inv) || anyNA(sx_inv)) {
    return(list(value = value, status = "singular_covariance"))
  }
  x_centered <- scale(x_f, scale = FALSE)
  x_centered[is.na(x_centered)] <- 0
  value[keep] <- matrixStats::rowSums2((x_centered %*% sx_inv) * x_centered)
  list(value = value, status = "ok")
}
