# Paper-faithful reference oracles for the nonparametric person-fit C/IER
# statistics wrapped by cier_gnormed() and cier_ht() (registry citation key:
# niessen2016). Consumed by tests/testthat/test-cier-gnormed.R and
# tests/testthat/test-cier-ht.R.
#
# These oracles are INDEPENDENT of the upstream routines (PerFit::Gnormed.poly,
# mokken::coefH): they recompute each statistic from its originating definition
# rather than from the package's scorer. Both are nonparametric -- a pure
# function of the response matrix (no IRT model is fitted).
#
# References:
#   Molenaar (1991); Emons (2008, APM) -- polytomous normed Guttman (Gnormed),
#     and Molenaar's weighted Loevinger H underlying polytomous Ht.
#   Sijtsma (1986); Sijtsma & Meijer (1992) -- Ht (person H / person scalability).
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
# Matches PerFit::Gnormed.poly to its 4-decimal output rounding
# and reduces to the dichotomous PerFit::Gnormed when ncat = 2.
ref_personfit_gnormed_poly <- function(data, ncat = NULL) {
  resp <- data$responses
  rk <- data$items$reverse_keyed
  if (any(rk)) {
    # The oracle's fixtures are base-1 (1..max), so (min + max) - x is
    # (1 + max) - x and the category count Ncat equals max.
    reflect <- data$items$max[rk] + 1L
    resp[, rk] <- rep(reflect, each = nrow(resp)) - resp[, rk]
  }
  if (is.null(ncat)) ncat <- as.integer(data$items$max[[1L]])
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

# ---- mokken: Ht polytomous (transposed Mokken scalability) -----------------

# Polytomous person Ht is the item scalability of the transposed response
# matrix (mokken::coefH(t(X))$Hi). Independent closed form: for the
# complete-case, reverse-scored block,
#   Ht_v = cov(x_v, Tot - x_v) / cov(sort(x_v), SortedTot - sort(x_v))
# over items, with Tot = colSums and SortedTot = colSums of the
# per-respondent-sorted rows. This is the bilinearity collapse of
# sum_{w!=v} cov(X_v, X_w) / sum_{w!=v} cov_max(X_v, X_w), where cov_max is
# the Frechet / comonotonic maximum (Hardy-Littlewood-Polya rearrangement:
# sort both, pair). It matches mokken::coefH to ~1e-14 and reduces to the
# dichotomous Ht (PerFit::Ht, to that package's 4-dp output rounding). Zero-
# variance (straightline) rows and respondents with any missing cell are NA, as
# in the production bridge; reverse-keyed items are reverse-scored first
# (cat + 1 - x), exactly as apply_split_half_keying does on a base-1 scale.
ref_personfit_ht_poly <- function(data, ncat = NULL) {
  resp <- data$responses
  rk <- data$items$reverse_keyed
  if (any(rk)) {
    # Base-1 fixtures: (min + max) - x is (1 + max) - x.
    reflect <- data$items$max[rk] + 1L
    resp[, rk] <- rep(reflect, each = nrow(resp)) - resp[, rk]
  }
  n <- nrow(resp)
  out <- rep(NA_real_, n)
  complete <- stats::complete.cases(resp)
  if (sum(complete) < 2L) {
    return(out)
  }
  z <- resp[complete, , drop = FALSE]
  # Mirror the bridge's all-constant guard: < 2 non-constant complete
  # respondents -> nothing scorable (coefH would error or return all NaN).
  if (sum(apply(z, 1L, stats::var) > 0) < 2L) {
    return(out)
  }
  z <- z - min(z)
  s <- t(apply(z, 1L, sort))
  tot <- colSums(z)
  sorted_tot <- colSums(s)
  pcov <- function(a, b) mean(a * b) - mean(a) * mean(b)
  v <- rep(NA_real_, nrow(z))
  for (r in seq_len(nrow(z))) {
    if (stats::var(z[r, ]) == 0) next
    den <- pcov(s[r, ], sorted_tot - s[r, ])
    if (den != 0) v[r] <- pcov(z[r, ], tot - z[r, ]) / den
  }
  out[complete] <- v
  out
}
