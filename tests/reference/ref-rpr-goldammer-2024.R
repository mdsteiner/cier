# Independent reference oracle for Resampled Personal Reliability (RPR) per
# Goldammer et al. (2024, *Behavior Research Methods*, 56(8), 8831-8851).
# Re-derives the statistic independently and NEVER calls the production kernel.
# Used by tests/testthat/test-cier-personal-reliability.R.
#
# Definition:
#   For each respondent r:
#     Repeat B times:
#       - For each scale block k, draw a uniform random permutation of the
#         items; the first ceil(k / 2) form the first half, the rest the second
#         half.
#       - Compute the per-half mean for each scale.
#       - Across scales, correlate first-half and second-half means.
#       - Apply Spearman-Brown correction `2r / (1 + r)`, clamp at `-1`, negate.
#     The respondent's RPR is the mean of the B per-iteration values
#     (na.rm = TRUE); a respondent whose every iteration is NA stays NA.
#
# Reproducibility (a deliberate WHITE-BOX constraint, not a tautology): the
# statistic above is re-derived from scratch, but the random-draw ORDER is
# coordinated with production so a fixed seed matches bytewise. Production
# resolves RPR with a plain `set.seed(seed)` + `for` over iterations
# (kernel_rpr()), each iteration looping respondents then -- via
# compute_block_splits() -- scales, drawing one permutation per scale with two
# or more items and NONE for a one-item scale. The reference mirrors that order:
# outer `for` over iterations, inner `vapply` over respondents, inner `for` over
# scales, `next` (no draw) on a one-item scale.

ref_rpr_iter <- function(x, blocks) {
  vapply(seq_len(nrow(x)), function(i) {
    firsts <- rep(NA_real_, length(blocks))
    seconds <- rep(NA_real_, length(blocks))
    row <- x[i, ]
    for (k in seq_along(blocks)) {
      cols <- blocks[[k]]
      k_n <- length(cols)
      if (k_n < 2L) next
      vals <- row[cols]
      perm <- sample.int(k_n, k_n)
      first_idx <- perm[seq_len(ceiling(k_n / 2L))]
      second_idx <- setdiff(perm, first_idx)
      f_mean <- mean(vals[first_idx], na.rm = TRUE)
      s_mean <- mean(vals[second_idx], na.rm = TRUE)
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
