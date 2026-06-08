# Paper-faithful reference oracle for the normed polytomous Guttman-error
# person-fit statistic (Gnormed) wrapped by cier_gnormed() (registry citation
# key: niessen2016). Consumed by tests/testthat/test-cier-gnormed.R.
#
# This oracle is INDEPENDENT of the upstream routine (PerFit::Gnormed.poly): it
# recomputes the statistic from its originating definition rather than from the
# package's person-fit function. For the nonparametric Guttman-error index the
# oracle is a pure function of the response matrix (no IRT model is fitted).
# Ht's oracle (mokken / PerFit) joins this file when its slice lands.
#
# References:
#   Molenaar (1991); Emons (2008, APM) -- polytomous normed Guttman (Gnormed).
#   Niessen, Meijer & Tendeiro (2016, JRP) -- survey C/IER orientation.

# ---- PerFit: Gnormed (normed polytomous Guttman errors) --------------------

# Gnormed (Emons 2008; Molenaar 1991) is the polytomous normed Guttman-error
# index Gp / max(Gp | score). Independent closed form (a different formulation
# from PerFit's ISD-matrix + G() + anti-diagonal denominator DP): expand each
# response to item-steps, give every item-step a popularity position
# (1 = most popular). For a respondent endorsing NC total steps,
#   Gp  = sum(position of endorsed steps) - NC(NC + 1)/2,
# because NC(NC + 1)/2 is the error-free Guttman position-sum and each excess
# unit is one inversion. The maximum for that NC is
#   maxG = (max prefix-constrained position-sum over allocations summing to NC)
#          - NC(NC + 1)/2,
# a max-plus knapsack over items (each item contributes a prefix of its steps).
# Gnormed = Gp / maxG, with 0/0 -> 0 (all-min / all-max respondents). Reverse-
# keyed items are reverse-scored first (cat + 1 - x) and respondents with any
# missing cell are NA, exactly as the production bridge; a block with fewer
# than three items abstains (PerFit's recursive denominator needs >= 3 items).
# Matches PerFit::Gnormed.poly to its 4-decimal output rounding (< 1e-4; see
# TOLERANCES.md) and reduces to the dichotomous PerFit::Gnormed when ncat = 2.
ref_personfit_gnormed_poly <- function(data, ncat = NULL) {
  resp <- data$responses
  rk <- data$items$reverse_keyed
  if (any(rk)) {
    reflect <- data$items$categories[rk] + 1L
    resp[, rk] <- rep(reflect, each = nrow(resp)) - resp[, rk]
  }
  if (is.null(ncat)) ncat <- as.integer(data$items$categories[[1L]])
  n <- nrow(resp)
  i <- ncol(resp)
  mm <- ncat - 1L
  out <- rep(NA_real_, n)
  complete <- stats::complete.cases(resp)
  if (sum(complete) < 2L || i < 3L) {
    return(out)
  }
  m <- resp[complete, , drop = FALSE] - 1L          # zero-base to 0..M
  storage.mode(m) <- "double"
  # Item-step popularities (item-major: item 1 steps 1..M, item 2 ...).
  probs <- vapply(seq_len(mm), function(s) colMeans(m >= s), numeric(i))
  if (mm == 1L) probs <- matrix(probs, ncol = 1L)
  pos_vect <- rank(-as.vector(t(probs)), ties.method = "first")
  posmat <- matrix(pos_vect, nrow = i, ncol = mm, byrow = TRUE)
  cs <- apply(posmat, 1L, cumsum)
  if (mm == 1L) cs <- matrix(cs, nrow = 1L)
  cumrank <- cbind(0, t(cs))                         # i x (M + 1); steps 0..M
  nc <- rowSums(m)
  sumpos <- vapply(seq_len(nrow(m)), function(r)
    sum(cumrank[cbind(seq_len(i), m[r, ] + 1L)]), numeric(1L))
  num <- sumpos - nc * (nc + 1) / 2
  # T[k + 1] = max position-sum over allocations of k steps (max-plus knapsack).
  total <- i * mm
  dp <- rep(-Inf, total + 1L)
  dp[[1L]] <- 0
  for (it in seq_len(i)) {
    ndp <- rep(-Inf, total + 1L)
    cr <- cumrank[it, ]
    for (k in seq_along(dp)) {
      if (is.finite(dp[[k]])) {
        for (s in 0:mm) {
          idx <- (k - 1L) + s + 1L
          if (idx <= total + 1L) {
            ndp[[idx]] <- max(ndp[[idx]], dp[[k]] + cr[[s + 1L]])
          }
        }
      }
    }
    dp <- ndp
  }
  max_g <- dp - (0:total) * (1:(total + 1L)) / 2
  res <- num / max_g[nc + 1L]
  res[is.nan(res)] <- 0
  out[complete] <- res
  out
}
