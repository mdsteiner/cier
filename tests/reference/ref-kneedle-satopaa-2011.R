# Independent reference for the Kneedle elbow detector of Satopaa, Albrecht,
# Irwin & Raghavan (2011, *2011 31st International Conference on Distributed
# Computing Systems Workshops*, <doi:10.1109/ICDCSW.2011.20>). Biemann et al.
# (2025) recommend it for a sample-specific Laz.R cutoff; cier exposes it as the
# `kneedle = TRUE` cutoff on cier_lazr().
#
# This oracle re-derives the parameter-free elbow from scratch and NEVER calls
# the production kernel (kneedle_knee / resolve_kneedle_cutoff in R/cutoff.R), so
# any divergence is attributable to the kernel, not a shared helper. There is no
# CRAN parity partner: the only R implementation (etam4260/kneedle) is
# GitHub-only and cannot be a CRAN Suggests, so this oracle is the sole parity
# check -- oracle-only trust like PR / RPR (see tests/reference/TOLERANCES.md).
#
# Definition (Satopaa et al. 2011, pp. 168-169, Section 3.1), the parameter-free
# form cier ships (no sensitivity S, no smoothing -- the raw sorted scores):
#   1. Sort the n finite input values v_(1) <= ... <= v_(n) ascending.
#   2. Build the normalised curve (x_i, y_i) where x_i = (i - 1) / (n - 1) and
#      y_i = (v_(i) - min(v)) / (max(v) - min(v)).
#   3. The knee is the point of greatest distance from the diagonal y = x (the
#      line joining the normalised curve's endpoints). For a convex + increasing
#      curve -- the right-skewed C/IER index distribution, where a moderate bulk
#      gives way to a high-scoring careless tail -- that is the point of maximum
#      *negative* deviation, argmin(y_i - x_i); for a concave + increasing curve
#      it is argmax(y_i - x_i). A decreasing curve is mapped onto its increasing
#      counterpart by reversing the normalised y vector, with the index mapped
#      back to the original sort order.
#   4. Return the value at the chosen sorted-vector index.
#
# cier_lazr only ever uses the convex + increasing case (high Laz.R = careless =
# upper tail), but the oracle carries the curve / direction arguments so it is a
# faithful, fully independent re-derivation of the cited algorithm.

ref_kneedle <- function(values,
                        curve     = c("convex", "concave"),
                        direction = c("increasing", "decreasing")) {
  curve     <- match.arg(curve)
  direction <- match.arg(direction)
  values <- values[is.finite(values)]
  n <- length(values)
  if (n < 3L) {
    return(list(value = NA_real_, knee_index = NA_integer_,
                status = "insufficient_items"))
  }
  if (length(unique(values)) == 1L) {
    return(list(value = as.numeric(values[[1L]]), knee_index = 1L,
                status = "all_equal"))
  }
  sorted_values <- sort(values)
  x_norm <- (seq_len(n) - 1L) / (n - 1L)
  y_min  <- min(sorted_values)
  y_max  <- max(sorted_values)
  y_norm <- (sorted_values - y_min) / (y_max - y_min)
  if (identical(direction, "decreasing")) {
    y_norm <- rev(y_norm)
  }
  diff_yx <- y_norm - x_norm
  knee_idx_increasing <- if (identical(curve, "convex")) {
    which.min(diff_yx)
  } else {
    which.max(diff_yx)
  }
  knee_idx <- if (identical(direction, "decreasing")) {
    n - knee_idx_increasing + 1L
  } else {
    knee_idx_increasing
  }
  list(value = as.numeric(sorted_values[[knee_idx]]),
       knee_index = as.integer(knee_idx), status = "ok")
}
