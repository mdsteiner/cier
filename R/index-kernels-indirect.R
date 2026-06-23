# Numerical kernels for the matrix-only indirect indices: longstring, IRV, Mahalanobis,
# person-total, and the synonym / antonym pairing kernels. Each maps a numeric response
# matrix to a per-respondent score vector -- except kernel_mahalanobis(), which returns
# list(value, status). Kernels are pure and never raise typed conditions.

# ---- Longstring -------------------------------------------------------------

# Maximum run length of consecutive identical responses per respondent, over the raw row.
# One columnwise pass (p - 1 vectorised transitions updating an n-length current-run and
# running-max) replaces n R-level rle() calls. The `eq[is.na(eq)] <- FALSE` reset
# reproduces rle() NA semantics: each NA is its own run, so values separated by NA do not
# merge and an all-NA row scores a max run of 1.
kernel_longstring <- function(responses) {
  n <- nrow(responses)
  current <- rep(1, n)
  best <- rep(1, n)
  for (j in seq_len(ncol(responses))[-1L]) {
    eq <- responses[, j] == responses[, j - 1L]
    eq[is.na(eq)] <- FALSE                 # an NA cell starts a new run
    current <- current * eq + 1            # extend the run, or reset to 1
    best <- pmax(best, current)
  }
  best
}

# ---- IRV (intra-individual response variability) ----------------------------

# Per-respondent sample standard deviation (denominator n - 1) across answered items.
# matrixStats::rowSds() shares stats::sd na.rm semantics: a row with fewer than two present
# values yields NA (never NaN), so abstention falls out with no extra guard.
kernel_irv <- function(responses) {
  # unname(): rowSds carries input row names into the score; the cier_index value /
  # flagged must be bare positional vectors.
  unname(matrixStats::rowSds(responses, na.rm = TRUE))
}

# ---- Mahalanobis ------------------------------------------------------------

# Per-respondent squared Mahalanobis distance D^2 from the column-mean centroid, using
# pairwise covariance. Centre by column means via scale(scale = FALSE), then zero-fill the
# centred NAs: zeros contribute nothing to the bilinear sum, so missing terms drop (NOT
# na.rm on the product, where one NA would collapse the row). The zero-fill lets one BLAS
# multiply replace the per-row apply() loop: D^2 = rowSums((X sx_inv) * X). Returns
# list(value, status), status a scalar:
#   "ok"                     - D^2 computed (a row is NA only when it is all-NA).
#   "insufficient_responses" - fewer than two respondents carry any data.
#   "singular_covariance"    - cov not invertible: solve() failed, or returned a non-finite
#                              inverse (a cov cell is NA because two items never
#                              co-answered; solve() does not error on an NA matrix, so
#                              anyNA() is needed atop tryCatch).
#   "indefinite_covariance"  - invertible but not positive definite: pairwise estimation
#                              assembles each cov cell from a different subsample, so under
#                              heavy/structured missingness Sigma can gain a negative
#                              eigenvalue while solve() still succeeds. The bilinear form is
#                              then signed (negative "squared distance"), invalid for every
#                              row, so the kernel abstains wholesale. chol() is the test;
#                              the inverse still comes from solve(). Repairing Sigma to the
#                              nearest positive-definite matrix was rejected: the smoothed
#                              statistic would silently mask the data problem.
# On non-"ok" statuses every value is NA; the wrapper raises the typed warning.
kernel_mahalanobis <- function(responses) {
  value <- rep(NA_real_, nrow(responses))
  keep <- rowSums(!is.na(responses)) > 0L
  if (sum(keep) < 2L) {
    return(list(value = value, status = "insufficient_responses"))
  }
  x_f <- responses[keep, , drop = FALSE]
  sx <- stats::cov(x_f, use = "pairwise")
  sx_inv <- tryCatch(solve(sx), error = function(e) NULL)
  if (is.null(sx_inv) || anyNA(sx_inv)) {
    return(list(value = value, status = "singular_covariance"))
  }
  if (is.null(tryCatch(chol(sx), error = function(e) NULL))) {
    return(list(value = value, status = "indefinite_covariance"))
  }
  x_centered <- scale(x_f, scale = FALSE)
  x_centered[is.na(x_centered)] <- 0
  value[keep] <- matrixStats::rowSums2((x_centered %*% sx_inv) * x_centered)
  list(value = value, status = "ok")
}

# ---- Person-total correlation (r_pbis) --------------------------------------

# Per-respondent Pearson correlation of each respondent's raw answered responses with the
# whole-sample per-item mean vector (item-total form: the respondent's own value is
# included in the item mean). Low values flag carelessness. Keying-insensitive and the
# means are whole-sample, so the kernel needs no item metadata.
#
# Vectorised masked-sum Pearson: present items mask out missing cells, the item mean m
# carries only over answered cells, one pass of masked sums covers the whole battery. Items
# answered by nobody give a 0/0 = NaN mean forced to 0, contributing nothing. Returns a
# bare numeric vector, NA where a respondent answered fewer than three items or either side
# has zero variance (a straightliner, or a flat item-mean profile).
#
# Zero-variance robustness: a constant row's deviation sum-of-squares is 0 in exact
# arithmetic, but floating cancellation lands it a few ulp either side -- tiny-negative
# sends sqrt() to NaN, tiny-positive leaks a spurious finite score. So the variance terms
# are clamped at 0 under the sqrt, and a constant row is detected exactly (masked min ==
# max) and forced to NA. suppressWarnings() guards the all-NA-row reduction.
kernel_person_total <- function(responses) {
  present <- !is.na(responses)
  k <- rowSums(present)
  xf <- responses
  xf[!present] <- 0
  col_n <- colSums(present)
  m <- colSums(xf) / col_n
  m[!is.finite(m)] <- 0
  sx  <- rowSums(xf)
  sxx <- rowSums(xf * xf)
  sxm <- as.numeric(xf %*% m)
  sm  <- as.numeric(present %*% m)
  smm <- as.numeric(present %*% (m * m))
  num <- sxm - sx * sm / k
  den <- sqrt(pmax(sxx - sx * sx / k, 0) * pmax(smm - sm * sm / k, 0))
  value <- num / den
  constant <- suppressWarnings(
    matrixStats::rowMins(responses, na.rm = TRUE) ==
      matrixStats::rowMaxs(responses, na.rm = TRUE)
  )
  value[k < 3L | constant | !is.finite(value)] <- NA_real_
  as.numeric(value)
}

# ---- Psychometric synonyms / antonyms ---------------------------------------

# Whole-sample inter-item Pearson correlation used to discover synonym / antonym item
# pairs. Pairwise-complete, so a missing cell drops only the pairs it touches.
# suppressWarnings: a constant item column yields NA correlations (forms no pair) and
# stats::cor() would emit a "standard deviation is zero" warning that must not leak.
pairing_cor <- function(responses) {
  suppressWarnings(stats::cor(responses, use = "pairwise.complete.obs"))
}

# Discover the item pairs whose whole-sample inter-item correlation clears the `critical_r`
# magnitude. Returns a 2-column integer matrix with pairs[, 1] > pairs[, 2] (lower-triangle
# indices, larger item first), or a 0-row matrix when none qualify. `pairing` selects the
# tail: "syn" keeps r > critical_r, "ant" keeps r < -critical_r. `cor_mat` optionally
# injects a precomputed pairing_cor() so a caller sweeping thresholds builds the p x p
# matrix once; NULL (default) computes it here.
find_item_pairs <- function(responses, critical_r, pairing, cor_mat = NULL) {
  critical_r <- abs(critical_r)
  if (is.null(cor_mat)) {
    cor_mat <- pairing_cor(responses)
  }
  p <- ncol(cor_mat)
  diag(cor_mat) <- NA_real_
  cor_mat[upper.tri(cor_mat, diag = FALSE)] <- NA_real_
  flat_idx <- if (identical(pairing, "syn")) {
    which(cor_mat > critical_r)
  } else {
    which(cor_mat < -critical_r)
  }
  if (length(flat_idx) == 0L) {
    return(matrix(integer(0L), ncol = 2L))
  }
  row_idx <- ((flat_idx - 1L) %% p) + 1L
  col_idx <- ((flat_idx - 1L) %/% p) + 1L
  matrix(c(row_idx, col_idx), ncol = 2L)
}

# Psychometric-synonyms / antonyms kernel: discover the qualifying pairs on the whole
# sample once, then score every respondent in one vectorised pass. The score is the
# within-respondent Pearson correlation across the qualifying pairs between stacked
# first-item values (A) and second-item values (B), dropping pairs where either side is
# missing. Computed as a masked-sum Pearson accumulated over pair chunks, without keeping
# all qualifying pairs resident at once. Returns a bare numeric vector, NA where a
# respondent has fewer than three complete pairs or a pair side has zero variance. The two
# variance terms are clamped at 0 so floating noise cannot send sqrt() to NaN; a constant
# pair side is additionally detected exactly (masked min == max) and forced to NA.
# suppressWarnings() guards the no-complete-pair all-NA reduction. `pairing` is the
# discovery tail ("syn" / "ant"); pairing and scoring both use raw responses.
kernel_psychsyn <- function(responses, critical_r, pairing, cor_mat = NULL,
                            pairs = NULL, chunk_cells = 5e6) {
  # `pairs` optionally injects the precomputed find_item_pairs() set so a wrapper that
  # already discovered it, or a threshold sweep, does not re-run the lower-triangle scan.
  if (is.null(pairs)) {
    pairs <- find_item_pairs(responses, critical_r, pairing, cor_mat = cor_mat)
  }
  n <- nrow(responses)
  if (nrow(pairs) == 0L) {
    return(rep(NA_real_, n))
  }

  k <- integer(n)
  sa <- sb <- saa <- sbb <- sab <- numeric(n)
  a_min <- b_min <- rep(Inf, n)
  a_max <- b_max <- rep(-Inf, n)
  chunk_pairs <- max(1L, min(nrow(pairs), floor(chunk_cells / max(n, 1L))))
  starts <- seq.int(1L, nrow(pairs), by = chunk_pairs)

  for (start in starts) {
    end <- min(nrow(pairs), start + chunk_pairs - 1L)
    idx <- start:end
    a_cols <- responses[, pairs[idx, 1L], drop = FALSE]
    b_cols <- responses[, pairs[idx, 2L], drop = FALSE]
    mask <- !is.na(a_cols) & !is.na(b_cols)

    a_complete <- a_cols
    b_complete <- b_cols
    a_complete[!mask] <- NA_real_
    b_complete[!mask] <- NA_real_
    a_min <- pmin(
      a_min,
      suppressWarnings(matrixStats::rowMins(a_complete, na.rm = TRUE)),
      na.rm = TRUE
    )
    a_max <- pmax(
      a_max,
      suppressWarnings(matrixStats::rowMaxs(a_complete, na.rm = TRUE)),
      na.rm = TRUE
    )
    b_min <- pmin(
      b_min,
      suppressWarnings(matrixStats::rowMins(b_complete, na.rm = TRUE)),
      na.rm = TRUE
    )
    b_max <- pmax(
      b_max,
      suppressWarnings(matrixStats::rowMaxs(b_complete, na.rm = TRUE)),
      na.rm = TRUE
    )

    a_cols[!mask] <- 0
    b_cols[!mask] <- 0
    k <- k + rowSums(mask)
    sa <- sa + rowSums(a_cols)
    sb <- sb + rowSums(b_cols)
    saa <- saa + rowSums(a_cols * a_cols)
    sbb <- sbb + rowSums(b_cols * b_cols)
    sab <- sab + rowSums(a_cols * b_cols)
  }

  constant <- (a_min == a_max) | (b_min == b_max)
  num <- sab - sa * sb / k
  den <- sqrt(pmax(saa - sa * sa / k, 0) * pmax(sbb - sb * sb / k, 0))
  value <- num / den
  value[k <= 2L | constant | !is.finite(value)] <- NA_real_
  as.numeric(value)
}
