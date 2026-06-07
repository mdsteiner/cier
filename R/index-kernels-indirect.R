# Purpose: Low-level numerical kernels shared by the indirect-index public
#          wrappers (R/cier-longstring.R, ...). Each kernel takes a numeric
#          matrix of responses and returns the per-respondent scores it is
#          responsible for. Single-kernel rule: one production implementation
#          per statistic.
# Args:    See per-kernel documentation below.
# Returns: Numeric vectors; never raises typed errors (the wrappers validate).
# Invariants:
#   - Kernels are pure (no I/O, no global state) and never mutate inputs.

# ---- Longstring -------------------------------------------------------------

# Maximum run length of consecutive identical responses per respondent, over
# the *raw* row (no scale blocking). Bytewise compatible with
# careless::longstring(): base::rle() treats each NA as its own run (NA == NA
# is NA, not TRUE), so identical values separated by NA do not merge and an
# all-NA row yields a max run length of 1. The wrapper applies NA-abstention for
# rows with no present responses; the kernel itself stays pure.
kernel_longstring <- function(responses) {
  vapply(
    seq_len(nrow(responses)),
    function(i) max(rle(responses[i, ])$lengths),
    numeric(1L)
  )
}

# ---- IRV (intra-individual response variability) ----------------------------

# Per-respondent SAMPLE standard deviation (denominator n - 1) across the items
# a respondent answered. `matrixStats::rowSds()` is the C-vectorised equivalent
# of the original `apply(x, 1, stats::sd, na.rm = TRUE)`: it shares stats::sd's
# NA semantics under `na.rm = TRUE`, so a row with fewer than two present values
# yields NA (never NaN) and abstention falls out of the summary itself with no
# extra guard (unlike longstring, whose all-NA row scores 1). Its summation
# order differs from a hand-rolled two-pass formula at ulp level, so the
# reference oracle and careless::irv() are held to 1e-10, not bytewise (see
# tests/reference/TOLERANCES.md).
kernel_irv <- function(responses) {
  matrixStats::rowSds(responses, na.rm = TRUE)
}

# ---- Mahalanobis ------------------------------------------------------------

# Per-respondent SQUARED Mahalanobis distance D^2 of each response vector from
# the column-mean centroid, using pairwise covariance. Numerically equivalent to
# careless::mahad() and psych::outlier(plot = FALSE, bad = 0, na.rm = TRUE) --
# both of which return D^2 -- at 1e-10 (see tests/reference/TOLERANCES.md):
#   * sx <- stats::cov(x, use = "pairwise"); sx_inv <- solve(sx).
#   * Centre by column means via scale(scale = FALSE), then zero-fill the
#     remaining NAs in the CENTRED matrix. Zeros contribute nothing to the
#     bilinear sum, so the missing terms drop exactly as psych::outlier's
#     na.rm = TRUE skips them. This is NOT na.rm on the product: one NA there
#     turns a whole sx_inv %*% xc product NA and collapses the row to 0. The
#     zero-fill lets one BLAS matrix multiply replace the per-row apply() loop:
#     D^2 = diag(X sx_inv X^T) = rowSums((X sx_inv) * X).
# Returns list(value, status), status a scalar:
#   "ok"                     - D^2 computed (a row is NA only when it is all-NA).
#   "insufficient_responses" - fewer than two respondents carry any data.
#   "singular_covariance"    - the pairwise covariance is not invertible: solve()
#                              failed, OR it returned a non-finite inverse. The
#                              latter is the pairwise-specific case where a cov
#                              cell is NA (two items never co-answered, or an item
#                              answered by nobody): solve() does NOT error on an
#                              NA matrix, it returns an all-NA inverse, so an
#                              `anyNA()` check is required on top of tryCatch.
# On the two non-"ok" statuses every value is NA; the wrapper raises the typed
# warning. The kernel stays pure (it raises no conditions).
kernel_mahalanobis <- function(responses) {
  value <- rep(NA_real_, nrow(responses))
  keep <- rowSums(!is.na(responses)) > 0L
  if (sum(keep) < 2L) {
    return(list(value = value, status = "insufficient_responses"))
  }
  x_f <- responses[keep, , drop = FALSE]
  sx_inv <- tryCatch(solve(stats::cov(x_f, use = "pairwise")),
                     error = function(e) NULL)
  if (is.null(sx_inv) || anyNA(sx_inv)) {
    return(list(value = value, status = "singular_covariance"))
  }
  x_centered <- scale(x_f, scale = FALSE)
  x_centered[is.na(x_centered)] <- 0
  value[keep] <- matrixStats::rowSums2((x_centered %*% sx_inv) * x_centered)
  list(value = value, status = "ok")
}

# ---- Person-total correlation (r_pbis) --------------------------------------

# Per-respondent Pearson correlation of each respondent's RAW answered responses
# with the whole-sample per-item mean vector (item-TOTAL form, as in
# PerFit::r.pbis: the respondent's own value is included in the item mean). Low
# values flag carelessness -- a careless pattern decouples from the group's
# item-mean profile. Reverse-keying plays no part (keying-insensitive by design)
# and the means are whole-sample, never scale-blocked, so the kernel needs no
# item metadata.
#
# Vectorised masked-sum Pearson, identical (<= ~1.6e-14) to a per-row
# stats::cor(x_i, m) loop (the reference oracle ref_person_total()): for each
# respondent the present items mask out missing cells, the item mean m carries
# only over answered cells, and the correlation is assembled from masked sums so
# one pass covers the whole battery. Items answered by nobody give a 0/0 = NaN
# mean that is forced to 0 -- they then contribute nothing to the masked products
# rather than poisoning every row with NaN.
#
# Returns a bare numeric vector the length of nrow(responses): the per-respondent
# correlation, NA where a respondent answered fewer than three items or either
# side of the correlation has zero variance (a straightliner row, or a flat
# item-mean profile). The wrapper needs no per-row warning -- the only condition
# is the shared percentile-cutoff abstention when every respondent is NA -- so
# the kernel returns the score vector alone, not a list. It stays pure.
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
  den <- sqrt((sxx - sx * sx / k) * (smm - sm * sm / k))
  value <- num / den
  value[k < 3L | !is.finite(value)] <- NA_real_
  as.numeric(value)
}

# ---- Even-odd / split-half family -------------------------------------------

# Group column indices by their scale label, in first-appearance order, so a
# heterogeneously-ordered battery still yields a reproducible block sequence.
# `items$scale` is aligned to the columns of `responses`. Shared by the even-odd
# kernel here and the personal-reliability split-half kernels.
scale_block_indices <- function(items) {
  scales <- items$scale
  uniq <- unique(scales)
  out <- vector("list", length(uniq))
  names(out) <- uniq
  for (lab in uniq) {
    out[[lab]] <- which(scales == lab)
  }
  out
}

# Spearman-Brown correction 2r/(1+r) with the careless clamp at -1. NA in -> NA
# out; a perfect inverse (r = -1) sends the ratio to -Inf, which the clamp pulls
# to -1 (so the negated index value tops out at +1, the most-careless score).
spearman_brown_clamp <- function(r) {
  if (is.na(r)) {
    return(NA_real_)
  }
  val <- (2 * r) / (1 + r)
  if (is.na(val) || val < -1) {
    return(-1)
  }
  val
}

# Reverse-score reverse-keyed items so the split-half family forms its half-means
# on trait-aligned responses. A PER-METHOD keying step (not a global rescore:
# longstring / irv / person-total need the raw responses), so only the split-half
# wrappers call it. Strict no-op (returns `responses` unchanged) when no item is
# reverse-keyed -- preserving the bytewise careless parity of no-reverse data --
# and NA-preserving. Each reverse item is reflected with the self-inverse
# (min + max) - x, where max = min + categories - 1 and min is the scale base
# (items$min, default 1 -> the classic (categories + 1) - x for 1..categories
# coding); declaring a 0-based or bipolar base keeps the reflection on the same
# range. In the wrapper path `items` has been validated by check_items()
# (categories present and non-NA, min finite/whole on every reverse item, min
# defaulted to 1); the categories and min guards below are defensive backstops for
# direct callers, keeping the single-kernel reuse path safe.
apply_split_half_keying <- function(responses, items, call = rlang::caller_env()) {
  rk <- items$reverse_keyed
  if (!any(rk)) {
    return(responses)
  }
  cats <- items$categories
  if (is.null(cats) || anyNA(cats[rk])) {
    offending <- if (is.null(cats)) which(rk) else which(rk & is.na(cats))
    cier_abort(
      "cier_error_input",
      c("Reverse-keyed items need a known {.field categories} to reverse-score.",
        "x" = "Reverse-keyed item(s) with no category count: {.val {offending}}."),
      data = list(arg = "items$categories", observed = offending), call = call
    )
  }
  mins <- items$min                       # scale base; default 1 for direct callers
  if (is.null(mins)) {
    mins <- rep(1L, length(rk))
  } else if (anyNA(mins[rk])) {
    offending <- which(rk & is.na(mins))
    cier_abort(
      "cier_error_input",
      c("Reverse-keyed items need a known {.field min} (scale base) to reverse-score.",
        "x" = "Reverse-keyed item(s) with no scale base: {.val {offending}}."),
      data = list(arg = "items$min", observed = offending), call = call
    )
  }
  rev_cols <- responses[, rk, drop = FALSE]
  rev_max  <- mins[rk] + cats[rk] - 1L    # scale maximum per reverse item
  reflect  <- mins[rk] + rev_max          # (min + max) self-inverse reflection
  responses[, rk] <- rep(reflect, each = nrow(responses)) - rev_cols
  responses
}

# Even-odd split: EVEN within-scale positions form the "first" half, ODD the
# "second". cor() is symmetric, so this reproduces the even-vs-odd correlation
# exactly (and bytewise with careless::evenodd()). Only called with k >= 2
# (single-item scales are skipped in kernel_split_half_row).
even_odd_split_fn <- function(k) {
  pos <- seq_len(k)
  list(first_idx  = pos[pos %% 2L == 0L],
       second_idx = pos[pos %% 2L == 1L])
}

# Per-respondent split-half consistency over the scale blocks. `split_fn(k)`
# returns list(first_idx, second_idx) of within-scale positions. For each block
# of >= 2 items it forms the two half-means; across blocks it correlates the
# first-half vs second-half mean vectors (pairwise-complete), Spearman-Brown
# corrects, and returns the NEGATED value (high = careless). Returns a single
# numeric, NA where the row abstains: fewer than two blocks yield a finite
# half-mean pair, or the across-block correlation is undefined (a straightliner /
# flat profile -> zero variance). Shared kernel for even-odd here and personal
# reliability, per the single-kernel rule.
kernel_split_half_row <- function(row, blocks, split_fn) {
  n_blocks <- length(blocks)
  first_means  <- rep(NA_real_, n_blocks)
  second_means <- rep(NA_real_, n_blocks)
  for (k in seq_len(n_blocks)) {
    cols <- blocks[[k]]
    if (length(cols) < 2L) {
      next
    }
    vals <- row[cols]
    split <- split_fn(length(cols))
    first_means[[k]]  <- mean(vals[split$first_idx],  na.rm = TRUE)
    second_means[[k]] <- mean(vals[split$second_idx], na.rm = TRUE)
  }
  first_means[is.nan(first_means)]   <- NA_real_
  second_means[is.nan(second_means)] <- NA_real_
  if (sum(!is.na(first_means) & !is.na(second_means)) < 2L) {
    return(NA_real_)
  }
  r <- suppressWarnings(stats::cor(first_means, second_means,
                                   use = "pairwise.complete.obs"))
  sb <- spearman_brown_clamp(r)
  if (is.na(sb)) NA_real_ else -sb
}

# Even-odd consistency kernel: kernel_split_half_row with the even/odd split, one
# row at a time. Returns a per-respondent numeric vector (NA where the row
# abstains). The wrapper reverse-scores keyed items before calling and resolves
# the percentile abstention when every row is NA.
kernel_even_odd <- function(responses, blocks) {
  vapply(
    seq_len(nrow(responses)),
    function(i) kernel_split_half_row(responses[i, ], blocks, even_odd_split_fn),
    numeric(1L)
  )
}
