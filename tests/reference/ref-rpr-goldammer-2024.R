# Independent reference oracle for Resampled Personal Reliability (RPR) per
# Goldammer et al. (2024, *Behavior Research Methods*, 56(8), 8831-8851).
# Re-derives the statistic independently and NEVER calls the production kernel.
# Used by tests/testthat/test-cier-personal-reliability.R.
#
# Definition:
#   For each respondent r:
#     Repeat B times:
#       - For each scale block k, draw a uniform random split of the items into a
#         first half of ceil(k / 2) items and the rest; the per-half mean is the
#         mean of that half's present values.
#       - Across scales, correlate first-half and second-half means.
#       - Apply Spearman-Brown correction `2r / (1 + r)`, clamp at `-1`, negate.
#     The respondent's RPR is the mean of the B per-iteration values
#     (na.rm = TRUE); a respondent whose every iteration is NA stays NA.
#
# Reproducibility (a deliberate WHITE-BOX constraint, not a tautology): the
# statistic above is re-derived from scratch, but the random-draw ORDER is
# coordinated with production so a fixed seed matches (to 1e-10). Production draws
# the uniform split by RANK-OF-UNIFORMS -- one runif(n * len) matrix per >= 2-item
# scale per iteration, the ceil(len / 2) smallest-uniform positions forming each
# respondent's first half (kernel_rpr()); the outer loop is over iterations, the
# scales drawn in block order, one runif() per scale with two or more items and
# NONE for a one-item scale. This oracle mirrors that draw order EXACTLY (so the
# seeded uniforms agree), but scores INDEPENDENTLY: it recovers each respondent's
# half from the same uniforms with base `order()` (not matrixStats::rowRanks),
# forms half-means with base `mean()` (not the masked-sum kernel), and correlates
# with `stats::cor()` (not the vectorised masked-sum Pearson). The coordinated
# uniforms pin the split; the independent scoring is the parity check.
#
# Note: the rank-of-uniforms split is statistically identical to drawing a random
# permutation and taking its first half (every size-ceil subset equally likely);
# production switched to it purely to vectorise the draw. The per-seed values
# therefore differ from a sample.int-based mechanism, but the STATISTIC -- the
# expected split-half consistency over uniform splits -- is unchanged.

# One RPR iteration's per-respondent values, drawing the uniform splits in
# production's order (one runif(n * len) per >= 2-item scale, in block order) and
# scoring each respondent independently.
ref_rpr_iter <- function(x, blocks) {
  n <- nrow(x)
  n_blocks <- length(blocks)
  # Draw the per-scale uniform matrices up front, in block order, exactly as
  # production does -- so the consumed RNG stream matches.
  draws <- vector("list", n_blocks)
  for (k in seq_len(n_blocks)) {
    len <- length(blocks[[k]])
    if (len < 2L) next
    draws[[k]] <- matrix(stats::runif(n * len), nrow = n, ncol = len)
  }
  vapply(seq_len(n), function(i) {
    firsts <- rep(NA_real_, n_blocks)
    seconds <- rep(NA_real_, n_blocks)
    for (k in seq_len(n_blocks)) {
      cols <- blocks[[k]]
      len <- length(cols)
      if (len < 2L) next
      # The ceil(len / 2) smallest-uniform positions are this respondent's first
      # half (order() recovers the same subset rowRanks() selects in production).
      first_idx <- order(draws[[k]][i, ])[seq_len(ceiling(len / 2L))]
      vals <- x[i, cols]
      f_mean <- mean(vals[first_idx], na.rm = TRUE)
      s_mean <- mean(vals[-first_idx], na.rm = TRUE)
      if (is.nan(f_mean)) f_mean <- NA_real_
      if (is.nan(s_mean)) s_mean <- NA_real_
      firsts[[k]] <- f_mean
      seconds[[k]] <- s_mean
    }
    finite <- !is.na(firsts) & !is.na(seconds)
    if (sum(finite) < 2L) {
      return(NA_real_)
    }
    r <- suppressWarnings(stats::cor(firsts, seconds,
                                     use = "pairwise.complete.obs"))
    if (is.na(r)) return(NA_real_)
    sb <- (2 * r) / (1 + r)
    if (is.na(sb) || sb < -1) sb <- -1
    -sb
  }, numeric(1L))
}

ref_rpr <- function(x, blocks, n_resamples, seed) {
  if (!is.matrix(x)) x <- as.matrix(x)
  set.seed(seed)
  per_iter <- matrix(NA_real_, nrow = n_resamples, ncol = nrow(x))
  for (b in seq_len(n_resamples)) {
    per_iter[b, ] <- ref_rpr_iter(x, blocks)
  }
  vapply(seq_len(nrow(x)), function(i) {
    vals <- per_iter[, i]
    if (all(is.na(vals))) NA_real_ else mean(vals, na.rm = TRUE)
  }, numeric(1L))
}
