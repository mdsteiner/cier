# Independent reference for the careless pattern mutators + extent/onset that
# cier_simulate() uses to plant C/IER.
#
# These oracles re-derive each deterministic mutator from scratch and NEVER call
# the production kernels in R/sim-patterns.R / R/sim-plan.R, so any divergence is
# attributable to the kernel, not a shared helper. The simulator is a generator,
# not an index, and has no CRAN parity partner (verified 2026-06-11: no package
# simulates C/IER with a planted truth), so these closed-form oracles plus
# paper-anchored pins through the shipped indices (cier_lazr / cier_longstring /
# cier_autocorrelation) are the sole trust, oracle-only like personal reliability.

# Round half UP: floor(x + 0.5). The convention for
# the position map and the midpoint mutator, so a 4-option midpoint is the upper
# middle (3), matching the archive ceiling((K+1)/2). Re-derived here independently
# of sim_round_half_up() -- note this is NOT base R's round() (round half to even),
# which a mutant might use (round(2.5) == 2 vs ref_round_half_up(2.5) == 3).
ref_round_half_up <- function(x) {
  floor(x + 0.5)
}

# Per-item position value: round_half_up(min + q * (max - min)). `q` is a relative
# position in [0, 1] (scalar or length-p); `mins` / `maxs` are per-item. The closed
# form the straightline-position and the midpoint (q = 0.5) mutators emit.
ref_position_value <- function(q, mins, maxs) {
  ref_round_half_up(mins + q * (maxs - mins))
}

# Cyclic diagonal POSITION sequence over K categories: pos_j = ((start - 1) +
# step * (j - 1)) mod K + 1, j = 1..p. Returns a length-p vector of positions in
# 1..K; the item value is offset to min_j + pos - 1. A diagonal-liner that wraps.
ref_diagonal_cyclic <- function(start, step, k, p) {
  ((start - 1L) + step * (seq_len(p) - 1L)) %% k + 1L
}

# Bounce (triangle-wave) diagonal POSITION sequence: m = ((start - 1) +
# step * (j - 1)) mod 2(K - 1); pos = if (m < K) m + 1 else 2K - 1 - m. The
# up-then-back-down zigzag of Biemann et al.'s (2025) footnote 2.
ref_diagonal_bounce <- function(start, step, k, p) {
  period <- 2L * (k - 1L)
  m <- ((start - 1L) + step * (seq_len(p) - 1L)) %% period
  ifelse(m < k, m + 1L, 2L * k - 1L - m)
}

# Alternating POSITION sequence: values[((offset + j - 1) mod period) + 1], a
# cyclic walk through `values` (positions in 1..K). Period-2 with values c(1, K)
# is the classic high/low seesaw.
ref_alternating <- function(values, offset, period, p) {
  values[((offset + seq_len(p) - 1L) %% period) + 1L]
}

# The Biemann et al. (2025) footnote-2 reference sequence: K = 5, p = 10, start at
# category position 1, step 1, bounce. A literal so the production bounce is pinned
# byte-for-byte AND through the shipped cier_lazr (Laz.R = 2/3, while the cyclic
# diagonal 1,2,3,4,5,1,2,3,4,5 gives Laz.R = 1).
ref_biemann_bounce <- c(1L, 2L, 3L, 4L, 5L, 4L, 3L, 2L, 1L, 2L)

# Rank (Mann-Whitney) AUC of a careless-high `score` for the careless rows over
# the attentive rows, re-derived independently of the package: P(score_careless >
# score_attentive) with ties at 0.5. The recovery smoke tests assert this clears a
# modest per-pattern floor (the matched index ranks planted rows toward the top).
ref_rank_auc <- function(score, careless) {
  ok <- is.finite(score) & !is.na(careless)
  s <- score[ok]
  l <- as.logical(careless[ok])
  pos <- s[l]
  neg <- s[!l]
  if (!length(pos) || !length(neg)) {
    return(NA_real_)
  }
  r <- rank(c(pos, neg))
  (sum(r[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) /
    (length(pos) * length(neg))
}
