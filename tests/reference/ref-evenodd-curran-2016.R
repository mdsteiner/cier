# Paper-faithful reference implementation of the even-odd consistency index per
# Curran (2016, *JESP*) and Jackson (1976). Used by
# tests/testthat/test-cier-even-odd.R as an INDEPENDENT cross-check of the
# production kernel: it re-derives the statistic step by step and never calls
# cier's kernel_even_odd() / kernel_split_half_row().
#
# Definition:
#   For each respondent r and scale block k, compute:
#     even_rk = mean over items at within-scale EVEN positions
#     odd_rk  = mean over items at within-scale ODD  positions
#   Then for each respondent:
#     r_r     = cor(even_r, odd_r) across scales (pairwise-complete)
#     sb_r    = 2*r_r / (1 + r_r), clamped at >= -1
#     value_r = -sb_r
#   (so HIGHER values flag carelessness, matching careless v1.2.0+).
#
# The oracle assumes trait-aligned input: callers that exercise reverse-keying
# must reverse-score the matrix independently ((categories + 1) - x) before
# passing it here, so the keying and the Curran (2016) statistic are checked from
# separate sources.

ref_even_odd_row <- function(row, blocks) {
  evens <- rep(NA_real_, length(blocks))
  odds  <- rep(NA_real_, length(blocks))
  for (k in seq_along(blocks)) {
    cols <- blocks[[k]]
    if (length(cols) < 2L) next
    vals <- row[cols]
    pos <- seq_along(vals)
    even_vals <- vals[pos %% 2L == 0L]
    odd_vals  <- vals[pos %% 2L == 1L]
    e_mean <- mean(even_vals, na.rm = TRUE)
    o_mean <- mean(odd_vals,  na.rm = TRUE)
    if (is.nan(e_mean)) e_mean <- NA_real_
    if (is.nan(o_mean)) o_mean <- NA_real_
    evens[[k]] <- e_mean
    odds[[k]]  <- o_mean
  }
  finite <- !is.na(evens) & !is.na(odds)
  if (sum(finite) < 2L) {
    return(NA_real_)
  }
  r <- suppressWarnings(stats::cor(evens, odds,
                                   use = "pairwise.complete.obs"))
  if (is.na(r)) {
    return(NA_real_)
  }
  sb <- (2 * r) / (1 + r)
  if (is.na(sb) || sb < -1) sb <- -1
  -sb
}

ref_even_odd <- function(x, blocks) {
  if (!is.matrix(x)) x <- as.matrix(x)
  vapply(seq_len(nrow(x)), function(i) {
    ref_even_odd_row(x[i, ], blocks)
  }, numeric(1L))
}
