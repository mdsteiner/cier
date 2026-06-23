# Independent reference for the per-respondent maximum-absolute-lag
# autocorrelation index of Gottfried, Jezek, Kralova & Rihacek (2022,
# *Practical Assessment, Research & Evaluation*, 27(2),
# <doi:10.7275/vyxb-gt24>), as packaged in `responsePatterns::rp.acors()`.
#
# This oracle re-derives the statistic from base R (`stats::var` +
# `stats::cor`) and NEVER calls the production kernel, so any divergence is
# attributable to the kernel, not a shared helper. It is the genuine
# independent check because the only CRAN partner (`responsePatterns`) is the
# authors' own implementation.
#
# Definition (per respondent, per lag in min_lag:max_lag):
#   row  <- response row
#   if (na_rm) row <- row[!is.na(row)]            # strip-then-lag variant
#   if (length(row) - lag < 3)        ac <- NA    # slice too short for a lag
#   else:
#     row1 <- row[1:(p_eff - lag)]; row2 <- row[(1 + lag):p_eff]
#     v1 <- var(row1, na.rm = TRUE);  v2 <- var(row2, na.rm = TRUE)
#     if (is.na(v1) || is.na(v2))     ac <- NA    # < 2 non-NA in a slice
#     else if (v1 == 0 || v2 == 0)    ac <- 1     # zero-variance -> 1 (NOT acf)
#     else if (#complete pairs < 3)   ac <- NA    # a 2-pair cor is +/-1
#     else ac <- cor(row1, row2, use = "pairwise.complete.obs")
# Per respondent: value = max(|ac|) over lags; NA when every lag is NA.
#
# Default lag range: `max_lag = NULL` resolves to
# `min(p - 3, 10)` -- the low lags Gottfried et al. recommend -- mirroring the
# wrapper's default so default-vs-oracle parity tests compare at the same range
# (the exact formula is pinned independently by the wrapper's default test).
#
# NA convention: `na_rm = FALSE` handles missingness
# PAIRWISE (a missing cell drops only the lagged pairs it touches), matching
# `responsePatterns` 0.3.x -- NOT the harsher / buggy 0.1.1 behaviour. The
# zero-variance check uses the FULL non-NA slice (`na.rm = TRUE`), distinct from
# the paired-complete subset the correlation is computed over. A lag whose
# pairwise-complete subset has fewer than 3 pairs
# abstains (NA) -- a 2-pair correlation is +/-1 by construction -- but a
# zero-variance slice still scores 1 first (the straightliner convention wins).
# Only the `na_rm = FALSE` pairwise path can reach < 3 complete pairs: under
# `na_rm = TRUE` the row is compacted to non-NA cells first, so complete pairs
# equal the slice length, already >= 3 by the slice-too-short guard.

ref_autocorrelation_row <- function(row, min_lag = 1L, max_lag = NULL,
                                    na_rm = FALSE) {
  p <- length(row)
  if (is.null(max_lag)) {
    max_lag <- min(p - 3L, 10L)
  }
  min_lag <- as.integer(min_lag)
  max_lag <- as.integer(max_lag)
  if (na_rm) {
    row <- row[!is.na(row)]
  }
  p_eff <- length(row)
  lags <- seq.int(min_lag, max_lag)
  per_lag <- rep(NA_real_, length(lags))
  for (k in seq_along(lags)) {
    lag <- lags[[k]]
    if (p_eff - lag < 3L) {
      next
    }
    row1 <- row[seq_len(p_eff - lag)]
    row2 <- row[seq.int(1L + lag, p_eff)]
    v1 <- stats::var(row1, na.rm = TRUE)
    v2 <- stats::var(row2, na.rm = TRUE)
    if (is.na(v1) || is.na(v2)) {
      per_lag[[k]] <- NA_real_
    } else if (v1 == 0 || v2 == 0) {
      per_lag[[k]] <- 1
    } else if (sum(!is.na(row1) & !is.na(row2)) < 3L) {
      per_lag[[k]] <- NA_real_       # < 3 complete pairs (cor would be +/-1)
    } else {
      per_lag[[k]] <- suppressWarnings(
        stats::cor(row1, row2, method = "pearson",
                   use = "pairwise.complete.obs")
      )
    }
  }
  names(per_lag) <- paste0("lag", lags)
  abs_lag <- abs(per_lag)
  max_abs_ac <- if (any(!is.na(abs_lag))) max(abs_lag, na.rm = TRUE) else NA_real_
  list(per_lag = per_lag, max_abs_ac = max_abs_ac)
}

ref_autocorrelation <- function(x, min_lag = 1L, max_lag = NULL,
                                na_rm = FALSE) {
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (is.null(max_lag)) {
    max_lag <- min(ncol(x) - 3L, 10L)
  }
  lapply(seq_len(nrow(x)),
         function(i) ref_autocorrelation_row(x[i, ], min_lag, max_lag, na_rm))
}

# Convenience: the per-respondent max-absolute-lag autocorrelation vector (the
# production `value`), so a test can compare in one line.
ref_autocorrelation_value <- function(x, min_lag = 1L, max_lag = NULL,
                                      na_rm = FALSE) {
  vapply(ref_autocorrelation(x, min_lag, max_lag, na_rm),
         function(r) r$max_abs_ac, numeric(1L))
}
