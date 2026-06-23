# Independent reference for the attentive graded-response model (GRM) that
# cier_simulate() uses to generate careful-respondent data.
#
# These oracles re-derive the attentive layer from scratch and NEVER call the
# production kernels in R/sim-attentive.R / R/sim-marginals.R, so any divergence
# is attributable to the kernel, not a shared helper. The simulator is a
# generator, not an index, and has no CRAN parity partner (verified 2026-06-11:
# no package simulates C/IER with a planted truth), so these closed-form oracles
# plus paper-anchored pins through the shipped indices are the sole
# trust, oracle-only like personal reliability.
#
# Three pieces:
#   1. ref_grm_pmf(preset, K, strength)  -- the target marginal pmf each named
#      preset stands for, re-derived from the agreed closed forms. This is the
#      EXPECTED category-probability vector for an item with K categories.
#   2. ref_grm_categorise(eta, thresholds) -- categorise continuous latent
#      values with an explicit threshold-count double loop (no findInterval()).
#   3. ref_marginal_implied(thresholds) -- the pmf a threshold vector implies on
#      a unit-variance normal latent, P(cat = c) = Phi(tau_c) - Phi(tau_{c-1}).

# The five named marginal-shape presets. Every
# preset returns a strictly-positive pmf of length K (no zero-probability
# category, so tau_k = qnorm(cumsum(p)) stays finite and strictly increasing).
#   - uniform: flat, 1 / K.
#   - peaked: symmetric triangular w_k = min(k, K + 1 - k); K = 5 -> (1,2,3,2,1)/9.
#   - skewed_left: geometric climb w_k = strength^(k - 1), so mass piles on the
#     HIGH (agreement) categories; default strength 1.5 (>= 1, per-spec overridable).
#   - skewed_right: the mirror of skewed_left (mass on the low categories).
#   - bimodal: deep reflected-triangular U, w_k = (max(tri) + 1) - tri_k; K = 5 ->
#     (3,2,1,2,3)/11 (centre one third of the ends).
ref_grm_pmf <- function(preset, k, strength = 1.5) {
  k <- as.integer(k)
  tri <- pmin(seq_len(k), rev(seq_len(k)))          # min(j, K + 1 - j)
  w <- switch(
    preset,
    uniform      = rep(1, k),
    peaked       = tri,
    skewed_left  = strength^(seq_len(k) - 1L),
    skewed_right = rev(strength^(seq_len(k) - 1L)),
    bimodal      = (max(tri) + 1) - tri,
    stop("unknown preset: ", preset)
  )
  w / sum(w)
}

# Categorise an n x p matrix of latent values against per-item thresholds, with
# an explicit count loop instead of findInterval(): category = 1 + #{thresholds
# at or below eta}. The "at or below" (`<=`) boundary matches findInterval()'s
# default left-closed convention, so a latent value exactly on a threshold falls
# in the HIGHER category -- pinned here so a mutant flipping the boundary shows.
ref_grm_categorise <- function(eta, thresholds) {
  n <- nrow(eta)
  p <- ncol(eta)
  out <- matrix(0L, nrow = n, ncol = p)
  for (j in seq_len(p)) {
    tau <- thresholds[[j]]
    for (i in seq_len(n)) {
      out[i, j] <- 1L + sum(tau <= eta[i, j])
    }
  }
  out
}

# The marginal pmf a threshold vector implies on a standard-normal latent: the
# probability mass between consecutive cut-points, P(cat = c) = Phi(tau_c) -
# Phi(tau_{c-1}) with Phi(tau_0) = 0 and Phi(tau_K) = 1. This is the closed-form
# inverse of tau_k = qnorm(cumsum(p)): feeding the resolver's thresholds back
# through it must return the target pmf (to ~1e-12), the exact marginals oracle.
ref_marginal_implied <- function(thresholds) {
  diff(c(0, stats::pnorm(thresholds), 1))
}
