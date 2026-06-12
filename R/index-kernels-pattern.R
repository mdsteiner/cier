# Purpose: Low-level numerical kernels for the pattern-family indirect indices
#          (autocorrelation and Laz.R). Pure functions on a numeric response
#          matrix; the wrappers (R/cier-autocorrelation.R, R/cier-lazr.R)
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
# decision tree, with the D5 minimum-pairs guard:
#   * a slice with fewer than two non-NA elements -> NA;
#   * a zero-variance slice (constant over its non-NA elements) -> 1 (NOT the
#     acf() NaN; a straightliner is deliberately sent to the top of the tail);
#   * else fewer than THREE complete (both-present) pairs -> NA (D5: a 2-pair
#     correlation is +/-1 by construction, which was deterministically flagging
#     early dropouts; raised from rp.acors's 2-pair floor, recorded in
#     TOLERANCES.md). The zero-variance branch above takes PRECEDENCE, so a
#     constant slice still scores 1 even with fewer than three complete pairs;
#   * else the pairwise-complete Pearson correlation.
# Only the pairwise (na_rm = FALSE) path can reach < 3 complete pairs; on
# complete data every lag has p - lag >= 3 pairs (the slice-too-short guard), so
# D5 is a no-op there and the rp.acors / oracle parity on complete data is
# unchanged. The zero-variance test uses the FULL non-NA slice (matching
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
  zero_var <- !na_slice & (v1 == 0 | v2 == 0)      # clean logical (no NA)
  ac[na_slice] <- NA_real_
  ac[zero_var] <- 1                                # zero-variance -> 1 (precedence)
  ac[!is.finite(ac)] <- NA_real_                   # 0/1-pair degenerate
  ac[n_pairs < 3L & !zero_var] <- NA_real_         # D5: < 3 complete pairs -> NA
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

# ---- Laz.R ------------------------------------------------------------------

# Per-respondent Laz.R (Biemann, Koch-Bayram, Meier-Barthold & Aguinis 2025,
# Eq. 3): the average probability with which the previous answer predicts the
# next. Over each respondent's consecutive (cur, nxt) pairs with both endpoints
# present (an NA endpoint drops that transition -- the drop-NA convention),
# tabulate the transition counts T, row-normalise to P, and return
#   value = sum(P * T) / n_trans = sum_i (sum_j T_ij^2) / rowsum_i / n_trans.
# Higher = more predictable = more careless; the value lies in (~1/s, 1].
#
# The statistic is invariant to the assumed anchor count, so the kernel needs no
# `s`: it rank-maps the observed values to a dense 1..s index over the DISTINCT
# anchors actually used (so any base -- 0-based, bipolar -- and any unobserved
# anchor score identically). Ranking over distinct values, not the raw value
# span, keeps the bin space bounded by the number of distinct responses
# (<= ncol), so a stray large value (an un-recoded numeric missing code, a
# sentinel, a leaked id) cannot inflate `n * s * s` into a 32-bit integer
# overflow or a huge allocation. A respondent with fewer than two valid
# transitions abstains (NA): a single transition gives P_ij * T_ij / 1 = 1
# regardless of the pattern. Bare value vector (lean schema); the wrapper
# resolves the percentile cutoff. The flat-index tabulate + reshape is one
# vectorised pass, no per-row R call.
kernel_lazr <- function(responses) {
  n <- nrow(responses)
  p <- ncol(responses)
  if (p < 2L) {
    return(rep(NA_real_, n))
  }
  cur <- responses[, -p, drop = FALSE]
  nxt <- responses[, -1L, drop = FALSE]
  valid <- !is.na(cur) & !is.na(nxt)
  n_trans <- rowSums(valid)
  cells <- which(valid)
  if (length(cells) == 0L) {
    return(rep(NA_real_, n))
  }
  cur_v <- cur[cells]
  nxt_v <- nxt[cells]
  anchors <- sort(unique(c(cur_v, nxt_v)))         # distinct used anchors, ranked
  s <- length(anchors)
  cur_i <- match(cur_v, anchors)
  nxt_i <- match(nxt_v, anchors)
  row_idx <- ((cells - 1L) %% n) + 1L              # which() is column-major
  flat <- (row_idx - 1L) * (s * s) + (cur_i - 1L) * s + nxt_i
  counts <- tabulate(flat, nbins = n * s * s)
  # block[k, r] = T_r entry k (k = (cur - 1) * s + nxt, cur outer / nxt inner).
  # The s x (s * n) view groups the s within-(r, i) cells into one column, so
  # column sums give the per-(r, i) row sums of T.
  block <- matrix(as.numeric(counts), nrow = s, ncol = s * n)
  row_sum <- colSums(block)
  sq_sum <- colSums(block^2)
  contrib <- numeric(length(row_sum))
  nz <- row_sum > 0
  contrib[nz] <- sq_sum[nz] / row_sum[nz]
  per_row_num <- colSums(matrix(contrib, nrow = s, ncol = n))
  value <- per_row_num / n_trans
  value[n_trans < 2L] <- NA_real_                  # also rewrites 0/0 = NaN
  value
}
