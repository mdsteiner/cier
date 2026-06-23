# Numerical kernels for the pattern-family indirect indices (autocorrelation and Laz.R).
# Pure functions on the raw (as-clicked, no reverse-keying) response matrix; the wrappers
# validate input. Exception to the pure-value contract: kernel_lazr() raises a typed
# cier_error_input when the pooled distinct-anchor count exceeds the survey-plausible
# ceiling; the wrapper threads `call`.

# ---- Autocorrelation --------------------------------------------------------

# One lag of the per-respondent autocorrelation, vectorised across respondents.
# `row1_mat` / `row2_mat` are the n x (p - lag) lag slices. Per-respondent decision tree:
#   * < 2 non-NA elements -> NA;
#   * zero-variance slice (constant over non-NA elements) -> 1 (straightliner sent to the
#     top of the tail, not acf()'s NaN). Takes precedence: a constant slice scores 1 even
#     with < 3 complete pairs;
#   * < 3 complete (both-present) pairs -> NA;
#   * else pairwise-complete Pearson correlation, computed as a masked rowSums (no per-row
#     cor() call). The zero-variance test uses the full non-NA slice, not the paired subset.
# Only this (na_rm = FALSE) path can reach < 3 complete pairs; on complete data every lag
# has p - lag >= 3 pairs, so the minimum is a no-op there.
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
  ac[n_pairs < 3L & !zero_var] <- NA_real_         # < 3 complete pairs -> NA
  ac
}

# Reduce the per-lag autocorrelation matrix (n x n_lags) to the per-respondent max absolute
# autocorrelation. No finite lag -> NA; rowMaxs runs only over rows with a finite value, so
# it never returns -Inf.
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

# Per-respondent maximum absolute lag autocorrelation (Gottfried, Jezek, Kralova & Rihacek
# 2022). For each lag in min_lag:max_lag, correlate the row with its lag-shifted self; value
# is the max absolute correlation over lags (NA when every lag abstains). `na_rm = FALSE`
# (default) handles missingness pairwise within each lag; `na_rm = TRUE` strips each row's
# NAs before lagging (collapses administration-order spacing, rarely appropriate).
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

# `na_rm = TRUE` variant: NAs stripped per row first, so compacted lengths differ across rows
# and vectorised lag slicing does not apply -- a per-row loop on the stripped row, where
# var() / cor() need no na.rm.
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

# Row ranges over which kernel_lazr() chunks its transition tabulate. Each chunk tabulates
# chunk_n * s^2 bins, bounding peak allocation by cell_budget regardless of n; per-row counts
# are independent, so the chunked result is byte-identical to a single-pass tabulate. Returns
# a list of c(start, stop) ranges tiling 1..n exactly.
lazr_row_chunks <- function(n, ss, cell_budget) {
  chunk_n <- max(1L, as.integer(cell_budget %/% ss))
  starts <- seq.int(1L, n, by = chunk_n)
  lapply(starts, function(start) c(start, min(start + chunk_n - 1L, n)))
}

# Per-respondent Laz.R (Biemann, Koch-Bayram, Meier-Barthold & Aguinis 2025, Eq. 3): the
# average probability with which the previous answer predicts the next. Over each
# respondent's consecutive (cur, nxt) pairs with both endpoints present (an NA endpoint
# drops that transition), tabulate the transition counts T, row-normalise to P, and return
#   value = sum(P * T) / n_trans.
# Higher = more predictable = more careless; value lies in (~1/s, 1]. Fewer than two valid
# transitions abstains (NA), since a single transition gives 1 regardless of pattern.
#
# The statistic is invariant to the assumed anchor count, so the kernel rank-maps observed
# values to a dense 1..s index over the distinct used anchors. Bin space is bounded by the
# count of distinct pooled integer values, not ncol: a 0-100 slider or stray unique-integer
# column makes s large. Two guards keep `n * s * s` from exploding: (1) the tabulate is
# chunked over rows via lazr_row_chunks(), so peak allocation is chunk_n * s^2; (2) when s^2
# exceeds cell_budget (default 2^22 => s <= 2048) the kernel aborts typed (stray non-item
# column the usual cause). cell_budget is an internal seam so the chunk boundary and ceiling
# are testable.
kernel_lazr <- function(responses, cell_budget = 4194304, call = rlang::caller_env()) {
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
  ss <- as.numeric(s) * s                          # s^2 in double (overflow-proof)
  if (ss > cell_budget) {
    max_anchors <- as.integer(sqrt(cell_budget))   # floor(sqrt(.)); 2048 at default
    cier_abort(
      "cier_error_input",
      c("{.fn cier_lazr} found too many distinct response values to tabulate.",
        "x" = "Pooled distinct anchors: {s}, above the limit of {max_anchors} \\
               (the transition table needs {.code s^2} cells per respondent).",
        "i" = "A response scale rarely has this many anchors. A stray ID, \\
               timestamp, free-numeric, or un-recoded missing-code column is the \\
               usual cause -- pass only the item columns.",
        "i" = "For a genuinely fine-grained scale, bin the responses to fewer \\
               anchors first."),
      data = list(arg = "responses", observed = s, expected = max_anchors),
      call = call
    )
  }
  cur_i <- match(cur_v, anchors)
  nxt_i <- match(nxt_v, anchors)
  row_idx <- ((cells - 1L) %% n) + 1L              # which() is column-major
  per_row_num <- numeric(n)
  for (rng in lazr_row_chunks(n, ss, cell_budget)) {
    in_chunk <- row_idx >= rng[[1L]] & row_idx <= rng[[2L]]
    if (!any(in_chunk)) next                        # a chunk of wholly abstaining rows
    local_n <- rng[[2L]] - rng[[1L]] + 1L
    local_row <- row_idx[in_chunk] - rng[[1L]] + 1L
    ci <- cur_i[in_chunk]
    ni <- nxt_i[in_chunk]
    flat <- (local_row - 1) * ss + (ci - 1) * s + ni   # double: never overflows
    counts <- tabulate(flat, nbins = local_n * ss)
    # block[k, r] = T_r entry k (k = (cur - 1) * s + nxt, cur outer / nxt inner).
    # The s x (s * local_n) view groups the s within-(r, i) cells into one column,
    # so column sums give the per-(r, i) row sums of T.
    block <- matrix(as.numeric(counts), nrow = s, ncol = s * local_n)
    row_sum <- colSums(block)
    sq_sum <- colSums(block^2)
    contrib <- numeric(length(row_sum))
    nz <- row_sum > 0
    contrib[nz] <- sq_sum[nz] / row_sum[nz]
    per_row_num[rng[[1L]]:rng[[2L]]] <-
      colSums(matrix(contrib, nrow = s, ncol = local_n))
  }
  value <- per_row_num / n_trans
  value[n_trans < 2L] <- NA_real_                  # also rewrites 0/0 = NaN
  # unname(): n_trans = rowSums(valid) carries the input row names into value; the
  # cier_index value / flagged must be bare positional vectors.
  unname(value)
}
