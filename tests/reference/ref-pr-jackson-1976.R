# Independent reference oracle for Personal Reliability (PR) per Jackson (1976).
# Re-derives each respondent's -SB(r) step by step and NEVER calls the
# production kernel. Used by tests/testthat/test-cier-personal-reliability.R.
#
# Definition:
#   For each respondent and each scale block k, split items into a first half
#   (first ceil(n_k / 2) within-scale positions) and a second half (the
#   remaining floor(n_k / 2)). Compute the per-half mean. Across scales,
#   correlate the first-half-means and second-half-means. Apply the
#   Spearman-Brown correction `2r / (1 + r)`, clamp at >= -1, and return `-r`
#   so that higher values flag carelessness.
#
# `blocks` is a list of within-matrix column-index vectors, one per scale,
# built independently of production (blocks_from_scale() in the test file).

ref_pr_row <- function(row, blocks) {
  firsts  <- rep(NA_real_, length(blocks))
  seconds <- rep(NA_real_, length(blocks))
  for (k in seq_along(blocks)) {
    cols <- blocks[[k]]
    if (length(cols) < 2L) next
    vals <- row[cols]
    n_k <- length(vals)
    first_idx  <- seq_len(ceiling(n_k / 2L))
    second_idx <- setdiff(seq_len(n_k), first_idx)
    f_mean <- mean(vals[first_idx],  na.rm = TRUE)
    s_mean <- mean(vals[second_idx], na.rm = TRUE)
    if (is.nan(f_mean)) f_mean <- NA_real_
    if (is.nan(s_mean)) s_mean <- NA_real_
    firsts[[k]]  <- f_mean
    seconds[[k]] <- s_mean
  }
  finite <- !is.na(firsts) & !is.na(seconds)
  if (sum(finite) < 2L) {
    return(NA_real_)
  }
  r <- suppressWarnings(stats::cor(firsts, seconds,
                                   use = "pairwise.complete.obs"))
  if (is.na(r)) {
    return(NA_real_)
  }
  sb <- (2 * r) / (1 + r)
  if (is.na(sb) || sb < -1) sb <- -1
  -sb
}

ref_pr <- function(x, blocks) {
  if (!is.matrix(x)) x <- as.matrix(x)
  vapply(seq_len(nrow(x)),
         function(i) ref_pr_row(x[i, ], blocks),
         numeric(1L))
}
