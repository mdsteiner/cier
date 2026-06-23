# Independent reference for the per-respondent Laz.R index of Biemann,
# Koch-Bayram, Meier-Barthold & Aguinis (2025, *Organizational Research
# Methods*, <doi:10.1177/10944281251334778>).
#
# This oracle re-derives Eq. 3 from scratch with an explicit transition
# double-loop and NEVER calls the production kernel, so any divergence is
# attributable to the kernel, not a shared helper. Laz.R has no CRAN partner
# (verified 2026-06-10: neither `careless` nor any other package implements it),
# so this oracle is the sole parity check, oracle-only trust like PR / RPR.
#
# Definition (Eq. 3): for a response row r = (r_1, ..., r_N) of integer anchors,
#   T[i, j] = #{t in 1:(N-1) : r_t = i and r_{t+1} = j}   over consecutive pairs
#             with both endpoints present,
#   P[i, j] = T[i, j] / sum_j' T[i, j']   (rows whose total is zero stay zero),
#   Laz.R   = sum(P * T) / n_trans
# where `*` is the Hadamard product and n_trans is the count of valid
# transitions. Higher values indicate more predictable (more careless)
# responding; the value lies in (~1/s, 1].
#
# Conventions encoded here, matching the production wrapper:
#   - drop-NA: a transition is dropped whenever either endpoint is missing, so
#     the denominator is the count of VALID transitions, not N - 1.
#   - a row with fewer than TWO valid transitions abstains (NA): a single
#     transition gives P_ij * T_ij / 1 = 1 regardless of the underlying pattern.
#   - matrix-only: the anchor set is inferred from the observed integer values
#     (the Laz.R value is invariant to the assumed anchor count s), so the
#     transition matrix is built over sort(unique(observed)) and the value is
#     unchanged by any constant integer shift of the responses.

ref_lazr_row <- function(row) {
  row <- as.numeric(row)
  n_total <- length(row)
  if (n_total < 2L) {
    return(list(value = NA_real_, transition_matrix = matrix(0L, 0L, 0L)))
  }
  cur <- row[seq_len(n_total - 1L)]
  nxt <- row[seq.int(2L, n_total)]
  valid <- !is.na(cur) & !is.na(nxt)
  n_trans <- sum(valid)
  if (n_trans < 2L) {
    return(list(value = NA_real_, transition_matrix = matrix(0L, 0L, 0L)))
  }
  cur_v <- as.integer(round(cur[valid]))
  nxt_v <- as.integer(round(nxt[valid]))
  # Anchors are the observed integer values, ordered; the matrix is square over
  # them. This is the s-invariant form: an unobserved or higher anchor would add
  # an all-zero row/column and leave the value unchanged.
  anchors <- sort(unique(c(cur_v, nxt_v)))
  s <- length(anchors)
  tmat <- matrix(0L, nrow = s, ncol = s)
  for (k in seq_len(n_trans)) {
    i <- match(cur_v[[k]], anchors)
    j <- match(nxt_v[[k]], anchors)
    tmat[i, j] <- tmat[i, j] + 1L
  }
  pmat <- matrix(0, nrow = s, ncol = s)
  for (i in seq_len(s)) {
    rs <- sum(tmat[i, ])
    if (rs > 0L) {
      pmat[i, ] <- tmat[i, ] / rs
    }
  }
  value <- sum(pmat * tmat) / n_trans
  list(value = value, transition_matrix = tmat)
}

ref_lazr <- function(x) {
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  vapply(seq_len(nrow(x)),
         function(i) ref_lazr_row(x[i, ])$value, numeric(1L))
}

# John's worked example (Eq. 4 of the paper). 50 items over the four anchors
# 1..4. Paper text: "John's sequence of answers is '1-2-3-4-3-2-1-2...1-2'".
# The triangle-wave repeat unit is c(1, 2, 3, 4, 3, 2) of length 6; rep_len(., 50)
# is 8 full cycles (positions 1..48) followed by 1, 2 (positions 49..50). That
# yields the transition counts T_12 = 9 and T_21 = T_23 = T_32 = T_34 = T_43 = 8,
# so Laz.R = 33 / 49 = 0.673.
ref_lazr_john_sequence <- function() {
  pattern <- c(1L, 2L, 3L, 4L, 3L, 2L)
  rep_len(pattern, 50L)
}
