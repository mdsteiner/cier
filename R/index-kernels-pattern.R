# Purpose: Low-level numerical kernels for the pattern-family indirect indices
#          (autocorrelation now; Laz.R lands in a later slice). Pure functions on
#          a numeric response matrix; the wrappers (R/cier-autocorrelation.R)
#          validate input. Single-kernel rule: one production implementation per
#          statistic. Pattern indices read the as-clicked response sequence, so
#          they are computed on the RAW responses with no reverse-keying.
# Args:    See per-kernel documentation below.
# Returns: Numeric vectors; never raises typed errors (the wrappers validate).
# Invariants:
#   - Kernels are pure (no I/O, no global state) and never mutate inputs.

# ---- Autocorrelation --------------------------------------------------------

# One lag of the per-respondent autocorrelation, vectorised across respondents.
# `row1_mat` / `row2_mat` are the n x (p - lag) lag slices (columns 1..(p-lag)
# and (1+lag)..p of the response matrix). Returns the per-respondent
# autocorrelation for this lag following the responsePatterns::rp.acors()
# decision tree:
#   * a slice with fewer than two non-NA elements -> NA;
#   * a zero-variance slice (constant over its non-NA elements) -> 1 (NOT the
#     acf() NaN; a straightliner is deliberately sent to the top of the tail);
#   * fewer than two complete (both-present) pairs -> NA;
#   * else the pairwise-complete Pearson correlation.
# The zero-variance test uses the FULL non-NA slice (matching
# var(slice, na.rm = TRUE)), distinct from the paired-complete subset the
# correlation is summed over. The correlation is a masked-sum Pearson (the
# kernel_person_total technique): one set of rowSums replaces a per-row cor()
# call. It sums in a different order than stats::cor(), so it matches the
# independent oracle and responsePatterns to 1e-10, not bytewise (see
# tests/reference/TOLERANCES.md).
autocorrelation_lag <- function(row1_mat, row2_mat) {
  both_ok <- !is.na(row1_mat) & !is.na(row2_mat)
  n_pairs <- rowSums(both_ok)
  v1 <- matrixStats::rowVars(row1_mat, na.rm = TRUE)
  v2 <- matrixStats::rowVars(row2_mat, na.rm = TRUE)
  v1[rowSums(!is.na(row1_mat)) < 2L] <- NA_real_   # var() is NA on < 2 non-NA
  v2[rowSums(!is.na(row2_mat)) < 2L] <- NA_real_
  r1 <- row1_mat
  r2 <- row2_mat
  r1[!both_ok] <- 0                                # drop non-paired cells
  r2[!both_ok] <- 0
  s1 <- rowSums(r1)
  s2 <- rowSums(r2)
  denom <- n_pairs * (n_pairs - 1)
  cov12 <- (n_pairs * rowSums(r1 * r2) - s1 * s2) / denom
  vp1 <- (n_pairs * rowSums(r1 * r1) - s1 * s1) / denom
  vp2 <- (n_pairs * rowSums(r2 * r2) - s2 * s2) / denom
  ac <- cov12 / sqrt(vp1 * vp2)
  na_slice <- is.na(v1) | is.na(v2)
  ac[na_slice] <- NA_real_
  ac[!na_slice & (v1 == 0 | v2 == 0)] <- 1         # zero-variance -> 1
  ac[!is.finite(ac)] <- NA_real_                   # too-few-pairs / degenerate
  ac
}

# Reduce the per-lag autocorrelation matrix (n x n_lags) to the per-respondent
# maximum absolute autocorrelation. A respondent with no finite lag (every lag
# NA) abstains (NA); rowMaxs is therefore taken only over the rows that carry at
# least one finite value, so it never returns -Inf.
reduce_max_abs_ac <- function(per_lag) {
  abs_lag <- abs(per_lag)
  any_finite <- rowSums(!is.na(abs_lag)) > 0L
  out <- rep(NA_real_, nrow(per_lag))
  if (any(any_finite)) {
    out[any_finite] <- matrixStats::rowMaxs(
      abs_lag[any_finite, , drop = FALSE], na.rm = TRUE
    )
  }
  out
}

# Per-respondent maximum absolute lag autocorrelation (Gottfried, Jezek, Kralova
# & Rihacek 2022). For each lag in min_lag:max_lag, correlate the row with its
# lag-shifted self; the per-respondent value is the maximum absolute correlation
# over lags (NA when every lag abstains). With `na_rm = FALSE` (default)
# missingness is handled PAIRWISE within each lag; with `na_rm = TRUE` each row's
# NAs are stripped before lagging (which collapses administration-order spacing,
# so it is rarely appropriate for an order-sensitive index). The wrapper resolves
# the lag range and the percentile cutoff; the kernel returns the bare value
# vector (the lean-schema shape).
kernel_autocorrelation <- function(responses, min_lag, max_lag, na_rm) {
  if (isTRUE(na_rm)) {
    return(kernel_autocorrelation_na_rm(responses, min_lag, max_lag))
  }
  p <- ncol(responses)
  lags <- seq.int(min_lag, max_lag)
  per_lag <- matrix(NA_real_, nrow = nrow(responses), ncol = length(lags))
  for (k in seq_along(lags)) {
    lag <- lags[[k]]
    if (p - lag < 3L) {                            # slice too short for a lag
      next
    }
    row1 <- responses[, seq_len(p - lag), drop = FALSE]
    row2 <- responses[, seq.int(1L + lag, p), drop = FALSE]
    per_lag[, k] <- autocorrelation_lag(row1, row2)
  }
  reduce_max_abs_ac(per_lag)
}

# `na_rm = TRUE` variant: each respondent's NAs are stripped first, so the
# compacted length differs across rows and the vectorised lag slicing does not
# apply. The rare path (the default is `na_rm = FALSE`); a per-row loop on the
# stripped (NA-free) row, where var() / cor() need no na.rm. Matches the oracle's
# na_rm = TRUE branch.
kernel_autocorrelation_na_rm <- function(responses, min_lag, max_lag) {
  lags <- seq.int(min_lag, max_lag)
  vapply(seq_len(nrow(responses)), function(i) {
    row <- responses[i, ]
    row <- row[!is.na(row)]
    p_eff <- length(row)
    best <- NA_real_
    for (lag in lags) {
      if (p_eff - lag < 3L) {
        next
      }
      x <- row[seq_len(p_eff - lag)]
      y <- row[seq.int(1L + lag, p_eff)]
      ac <- if (stats::var(x) == 0 || stats::var(y) == 0) 1 else stats::cor(x, y)
      if (is.na(best) || abs(ac) > best) {
        best <- abs(ac)
      }
    }
    best
  }, numeric(1L))
}
