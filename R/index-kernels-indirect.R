# Purpose: Low-level numerical kernels shared by the indirect-index public
#          wrappers (R/cier-longstring.R, ...). Each kernel takes a numeric
#          matrix of responses and returns the per-respondent scores it is
#          responsible for. Single-kernel rule: one production implementation
#          per statistic.
# Args:    See per-kernel documentation below.
# Returns: The KERNELS return numeric vectors and never raise typed conditions
#          (the wrappers validate). The shared pre-score helpers in this file are
#          the documented exception: apply_split_half_keying() aborts on a bad
#          reverse-keying declaration and warn_two_scale_consistency() warns on a
#          degenerate two-scale design -- both are validation / diagnostic steps,
#          not numeric kernels.
# Invariants:
#   - Kernels are pure (no I/O, no global state) and never mutate inputs. The one
#     exception is kernel_rpr(): with a non-NULL seed it sets a local RNG seed and
#     restores the caller's `.Random.seed` on exit, so it is net state-preserving.

# ---- Longstring -------------------------------------------------------------

# Maximum run length of consecutive identical responses per respondent, over
# the *raw* row (no scale blocking). One columnwise pass over the whole matrix
# (p - 1 vectorised column transitions updating an n-length current-run and
# running-max) instead of n separate R-level rle() calls; integer run counting
# is exact, so the result is bytewise compatible with careless::longstring().
# Its rle() NA semantics are reproduced by the `eq[is.na(eq)] <- FALSE` reset:
# base::rle() treats each NA as its own run (NA == NA is NA, not TRUE), so
# identical values separated by NA do not merge and an all-NA row yields a max
# run length of 1. The wrapper applies NA-abstention for rows with no present
# responses; the kernel itself stays pure.
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
#   "indefinite_covariance"  - the pairwise covariance is invertible but NOT
#                              positive definite. Pairwise estimation assembles
#                              each cov cell from a different subsample, so under
#                              heavy or structured missingness the cells can be
#                              mutually inconsistent and Sigma gains a negative
#                              eigenvalue while solve() still succeeds. The
#                              bilinear form is then SIGNED -- a row can score a
#                              negative "squared distance" the upper-tail
#                              chi-square flag can never reach, and the ranking
#                              among the positive rows is distorted too -- so the
#                              distance is invalid for EVERY row and the kernel
#                              abstains wholesale. chol() is the test (it errors
#                              iff Sigma is not positive definite); the inverse
#                              still comes from solve(), keeping every
#                              positive-definite input byte-identical to
#                              careless::mahad / psych::outlier. Repairing Sigma
#                              instead (e.g. Matrix::nearPD) was deliberately
#                              rejected: the smoothed statistic traces to no
#                              cited paper or trusted package and would silently
#                              mask the data problem (see ADR.md, "Mahalanobis
#                              degenerate covariance: warn and abstain").
# On the non-"ok" statuses every value is NA; the wrapper raises the typed
# warning. The kernel stays pure (it raises no conditions).
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
#
# Zero-variance robustness: for a constant (straightliner) row the deviation
# sum-of-squares sxx - sx^2/k is 0 in exact arithmetic, but on a NON-integer
# constant floating cancellation lands it a few ulp on EITHER side of zero --
# tiny-negative would send sqrt() to NaN with a leaked base-R warning, and
# tiny-positive would leak a spurious finite ~1e-7 score into the percentile
# pool. So (a) the variance terms are clamped at 0 under the sqrt (the
# kernel_psychsyn technique, silencing the warning) and (b) a constant row is
# detected EXACTLY (masked min == max) and forced to NA, independent of which
# way the cancellation fell. The item-mean side needs no exact detection: a flat
# masked mean profile is exactly equal floats, cancelling to exactly 0 -> den
# 0 -> non-finite -> NA. suppressWarnings() guards the all-NA-row reduction
# (those rows abstain via k < 3 regardless).
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

# Warn when exactly two scale blocks are scorable (>= 2 items each; a one-item
# block yields a NULL split downstream and is skipped). With only two blocks the
# split-half consistency correlation is taken across two block means, so it is
# +/-1 by construction: a degenerate point mass for even-odd / PR, and a coarse
# average of +/-1 iterations for RPR. The trigger is a structural property of the
# item design (block sizes), so the warning fires once at the wrapper -- before
# scoring, independent of the responses and of the resample mode -- WITHOUT
# changing the statistic (careless parity preserved; the user may still want the
# coarse consistent / inconsistent split, or override with >= 3 scales). Shared by
# the even-odd and personal-reliability wrappers (single source). See ADR.md
# ("Split-half two-scale degeneracy", D6).
warn_two_scale_consistency <- function(blocks, call = rlang::caller_env()) {
  if (sum(lengths(blocks) >= 2L) == 2L) {
    cier_warn(
      "cier_warning_two_scale_consistency",
      c("Only two scorable scale blocks: the split-half consistency correlation \\
         is taken across just two block means, so it can only be {.val -1} or \\
         {.val 1}.",
        "i" = "The consistency score is therefore degenerate (clustered at the \\
               {.val -1} / {.val 1} extremes, with no ranking inside the careless \\
               tail). Provide >= 3 multi-item scales for a graded consistency \\
               score."),
      data = list(n_scorable = 2L), call = call
    )
  }
  invisible(NULL)
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
# (min + max) - x, where max is the largest response option (items$max) and min
# the scale base (items$min, default 1 -> the classic (max + 1) - x for 1..max
# coding); declaring a 0-based or bipolar base keeps the reflection on the same
# range. In the wrapper path `items` has been validated by check_items()
# (max present and non-NA, min finite/whole on every reverse item, min
# defaulted to 1); the max and min guards below are defensive backstops for
# direct callers, keeping the single-kernel reuse path safe. The observed-range
# cross-check below is NOT a backstop: type-valid metadata can still contradict
# the data, which only the data can reveal.
apply_split_half_keying <- function(responses, items, call = rlang::caller_env()) {
  rk <- items$reverse_keyed
  if (!any(rk)) {
    return(responses)
  }
  maxs <- items$max
  if (is.null(maxs) || anyNA(maxs[rk])) {
    offending <- if (is.null(maxs)) which(rk) else which(rk & is.na(maxs))
    cier_abort(
      "cier_error_input",
      c("Reverse-keyed items need a known {.field max} to reverse-score.",
        "x" = "Reverse-keyed item(s) with no scale maximum: {.val {offending}}."),
      data = list(arg = "items$max", observed = offending), call = call
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
  rev_max  <- maxs[rk]                    # scale maximum per reverse item
  # Cross-check the DECLARED range against the OBSERVED responses before
  # reflecting: a type-valid but wrong declaration (the classic case: 0-based
  # 0..4 data declared max = 5 with the default min = 1) would otherwise
  # reflect to off-scale values and silently corrupt the consistency score.
  # The person-fit bridges catch the same mistake in personfit_zero_base(); this
  # gives the split-half family (and Ht's reverse items) the equivalent guard.
  # An all-NA reverse column reduces to Inf/-Inf, which fails neither comparison
  # and so is (correctly) not an offender.
  obs_min <- suppressWarnings(matrixStats::colMins(rev_cols, na.rm = TRUE))
  obs_max <- suppressWarnings(matrixStats::colMaxs(rev_cols, na.rm = TRUE))
  bad <- obs_min < mins[rk] | obs_max > rev_max
  if (any(bad)) {
    offending <- which(rk)[bad]
    cier_abort(
      "cier_error_input",
      c("Reverse-keyed item responses must lie within the declared scale range \\
         {.val [min, max]}.",
        "x" = "Out-of-range reverse-keyed item(s): {.val {offending}}.",
        "i" = "Check {.field max} / {.field min} against the data (a \\
               0-based scale needs {.field min} = 0); reflecting an off-range \\
               value would silently corrupt the score."),
      data = list(arg = "items", observed = offending), call = call
    )
  }
  reflect  <- mins[rk] + rev_max          # (min + max) self-inverse reflection
  responses[, rk] <- rep(reflect, each = nrow(responses)) - rev_cols
  responses
}

# Even-odd split: EVEN within-scale positions form the "first" half, ODD the
# "second". cor() is symmetric, so this reproduces the even-vs-odd correlation
# exactly (and bytewise with careless::evenodd()). Only called with k >= 2
# (single-item scales yield a NULL split in compute_block_splits).
even_odd_split_fn <- function(k) {
  pos <- seq_len(k)
  list(first_idx  = pos[pos %% 2L == 0L],
       second_idx = pos[pos %% 2L == 1L])
}

# Precompute the within-scale split positions for every scale block: a per-block
# list aligned to `blocks`, each entry list(first_idx, second_idx) for a block of
# >= 2 items or NULL for a one-item block (skipped downstream). `split_fn(k)`
# maps a block length to its split. Computing the splits here -- once for the
# deterministic even-odd / PR callers (outside the per-respondent loop), or once
# per respondent for RPR's random splits -- keeps split construction out of the
# hot row loop in kernel_split_half_row(). RPR relies on this being one
# sample.int() draw per >= 2-item block and NONE for a one-item block, so its RNG
# stream stays in step with the reference oracle.
compute_block_splits <- function(blocks, split_fn) {
  lapply(blocks, function(cols) {
    if (length(cols) < 2L) NULL else split_fn(length(cols))
  })
}

# Per-respondent split-half consistency over the scale blocks. `splits` is the
# precomputed per-block list from compute_block_splits(): each entry is
# list(first_idx, second_idx) of within-scale positions, or NULL for a one-item
# block (skipped). For each scored block it forms the two half-means; across
# blocks it correlates the first-half vs second-half mean vectors
# (pairwise-complete), Spearman-Brown corrects, and returns the NEGATED value
# (high = careless). Returns a single numeric, NA where the row abstains: fewer
# than two blocks yield a finite half-mean pair, or the across-block correlation
# is undefined (a straightliner / flat profile -> zero variance). Shared kernel
# for even-odd and personal reliability, per the single-kernel rule.
kernel_split_half_row <- function(row, blocks, splits) {
  n_blocks <- length(blocks)
  first_means  <- rep(NA_real_, n_blocks)
  second_means <- rep(NA_real_, n_blocks)
  for (k in seq_len(n_blocks)) {
    split <- splits[[k]]
    if (is.null(split)) {
      next
    }
    vals <- row[blocks[[k]]]
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

# Build a deterministic split-half kernel: a `function(responses, blocks)` that
# precomputes the per-block split ONCE via `split_fn` (outside the per-respondent
# loop) and scores every respondent through kernel_split_half_row(). Even-odd and
# personal reliability (PR) are the SAME row loop differing only in their
# deterministic split, so they share this factory (single-kernel rule). RPR is
# deliberately NOT built here: its random splits must be drawn per respondent, so
# kernel_rpr() keeps its own loop.
make_split_half_kernel <- function(split_fn) {
  function(responses, blocks) {
    splits <- compute_block_splits(blocks, split_fn)
    vapply(
      seq_len(nrow(responses)),
      function(i) kernel_split_half_row(responses[i, ], blocks, splits),
      numeric(1L)
    )
  }
}

# Even-odd consistency kernel: the even/odd within-scale split. The wrapper
# reverse-scores keyed items before calling and resolves the percentile
# abstention when every row is NA.
kernel_even_odd <- make_split_half_kernel(even_odd_split_fn)

# ---- Personal reliability (PR / RPR) ----------------------------------------

# First/second-half split for Personal Reliability (Jackson 1976): the first
# ceil(k / 2) within-scale positions form the "first" half, the remaining
# floor(k / 2) the "second". Only called with k >= 2 (one-item scales yield a
# NULL split in compute_block_splits and are skipped).
split_half_indices <- function(k) {
  first <- seq_len(ceiling(k / 2L))
  # The second half is the suffix after the first ceil(k/2) positions; an index
  # slice avoids a set difference (first is always a prefix of seq_len(k)).
  list(first_idx = first, second_idx = seq_len(k)[-first])
}

# Personal Reliability (PR) kernel: the deterministic first/second-half split.
kernel_personal_reliability <- make_split_half_kernel(split_half_indices)

# Random within-scale split for RPR: a uniform permutation of the k positions,
# its first ceil(k / 2) entries forming the "first" half and the rest the
# "second". Consumes exactly one sample.int(k, k) draw per call. Only called with
# k >= 2 (one-item scales yield a NULL split and draw nothing -- see
# compute_block_splits -- so the RNG stream stays in step with the oracle).
random_split_indices <- function(k) {
  perm <- sample.int(k, k)
  # First ceil(k/2) of the shuffled permutation vs the rest; an index slice of
  # the already-distinct permutation avoids a set difference.
  half <- seq_len(ceiling(k / 2L))
  list(first_idx = perm[half], second_idx = perm[-half])
}

# Restore the global RNG state captured before a local set.seed(): re-assign the
# saved `.Random.seed`, or remove it when none existed, so the session returns to
# its pre-call state. Lets a seeded kernel be reproducible WITHOUT disturbing the
# caller's random stream.
restore_random_seed <- function(saved) {
  global <- globalenv()
  if (is.null(saved)) {
    if (exists(".Random.seed", envir = global, inherits = FALSE)) {
      rm(".Random.seed", envir = global)
    }
  } else {
    global[[".Random.seed"]] <- saved
  }
}

# Resampled Personal Reliability (RPR) kernel (Goldammer et al. 2024): the mean
# over `n_resamples` random within-scale split-half iterations of the per-row PR
# statistic. With a non-NULL seed the draw order reproduces the reference oracle
# (to 1e-10): iterations (outer), then respondents, then scales, with
# compute_block_splits() drawing one permutation per >= 2-item scale per
# respondent and none for a one-item scale. The seed is applied LOCALLY -- the
# caller's global RNG state is saved and restored on exit -- so a seeded call is
# reproducible WITHOUT disturbing the caller's random stream. A NULL seed draws
# from the ambient stream (results then vary per call). Per respondent the value
# is the mean of the finite per-iteration values (na.rm = TRUE), NA when every
# iteration abstains. The seed / for RNG logic lives here, not in the thin
# wrapper.
kernel_rpr <- function(responses, blocks, n_resamples, seed) {
  if (!is.null(seed)) {
    saved <- globalenv()[[".Random.seed"]]   # NULL when no RNG has been drawn yet
    on.exit(restore_random_seed(saved), add = TRUE)
    set.seed(seed)
  }
  n <- nrow(responses)
  per_iter <- matrix(NA_real_, nrow = n_resamples, ncol = n)
  for (b in seq_len(n_resamples)) {
    per_iter[b, ] <- vapply(seq_len(n), function(i) {
      splits <- compute_block_splits(blocks, random_split_indices)
      kernel_split_half_row(responses[i, ], blocks, splits)
    }, numeric(1L))
  }
  # Per-respondent mean over iterations; an all-NA column yields NaN -> NA.
  means <- matrixStats::colMeans2(per_iter, na.rm = TRUE)
  means[is.nan(means)] <- NA_real_
  means
}

# ---- Psychometric synonyms / antonyms ---------------------------------------

# Whole-sample inter-item Pearson correlation used to discover synonym / antonym
# item pairs. Pairwise-complete, so a missing cell drops only the pairs it
# touches. Factored out so the synonym and antonym wrappers build the pairing
# correlation with identical arguments. suppressWarnings: a constant
# (zero-variance) item column yields NA correlations -- it simply forms no pair --
# and stats::cor() would emit a base-R, locale-dependent "standard deviation is
# zero" warning that must not leak (the package signals only typed cli conditions).
pairing_cor <- function(responses) {
  suppressWarnings(stats::cor(responses, use = "pairwise.complete.obs"))
}

# Discover the item pairs whose whole-sample inter-item correlation clears the
# `critical_r` magnitude. Returns a 2-column integer matrix with pairs[, 1] >
# pairs[, 2] (lower-triangle column indices, the larger item index first), or a
# 0-row matrix when none qualify. `pairing` selects the tail: "syn" keeps
# r > critical_r (strong positive), "ant" keeps r < -critical_r (strong
# negative). Column-major traversal of the lower triangle makes the pair order
# and orientation byte-identical to careless:::get_item_pairs(), preserving the
# bytewise careless::psychsyn() parity. Shared by psychsyn (pairing = "syn") and
# the antonyms index (pairing = "ant"); the `pairing` tail is distinct from the
# registry flag direction. `cor_mat` optionally injects a precomputed
# pairing_cor(responses) so a caller evaluating several thresholds (the
# critical_r sweep) builds the p x p matrix once; NULL (the default) computes it
# here, byte-identically.
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

# Psychometric-synonyms / antonyms kernel: discover the qualifying pairs on the
# whole sample once, then score every respondent in a single vectorised pass. The
# per-respondent score is the Pearson correlation, within that respondent, across
# the K qualifying pairs between the stacked first-item values (A) and second-item
# values (B), dropping pairs where either side is missing. It is computed as a
# masked-sum Pearson over the n x K matrices A and B (the kernel_person_total
# technique): one set of rowSums replaces a per-row stats::cor() loop, a 4-6x
# speedup that scales to large samples without a compiled backend. Returns a bare
# numeric vector (the lean-schema shape), NA where a respondent has fewer than
# three complete pairs (careless's `> 2` guard) or a pair side has zero variance
# (correlation undefined); a respondent who answered nothing has no complete pairs
# and so abstains. The two variance terms are clamped at 0 so floating noise on a
# zero-variance side cannot send sqrt() to NaN with a warning -- such rows fall to
# NA through the finite check. The clamp alone covers only the tiny-NEGATIVE side
# of the cancellation: on a NON-integer constant the deviation sum-of-squares can
# also land tiny-POSITIVE, leaking a spurious finite score (~1.0) instead of the
# documented abstention -- so a constant pair side is additionally detected
# EXACTLY (masked min == max over the complete-pair values, NA-aligned before the
# zero-fill) and forced to NA. suppressWarnings() guards the no-complete-pair
# all-NA reduction (those rows abstain via k <= 2 regardless). `pairing` is the
# pair-discovery tail ("syn" / "ant"); pairing and scoring both use the raw
# responses, with no reverse-keying.
#
# Parity note: this masked-sum form is the same Pearson correlation as the per-row
# cor() loop it replaces, but sums in a different order, so it matches
# careless::psychsyn(resample_na = FALSE) to 1e-12 rather than bytewise (the loop
# was exact). The independent oracle parity (1e-12) is unaffected. See
# tests/reference/TOLERANCES.md and the ADR entry.
kernel_psychsyn <- function(responses, critical_r, pairing, cor_mat = NULL) {
  pairs <- find_item_pairs(responses, critical_r, pairing, cor_mat = cor_mat)
  n <- nrow(responses)
  if (nrow(pairs) == 0L) {
    return(rep(NA_real_, n))
  }
  a_cols <- responses[, pairs[, 1L], drop = FALSE]
  b_cols <- responses[, pairs[, 2L], drop = FALSE]
  mask <- !is.na(a_cols) & !is.na(b_cols)
  a_cols[!mask] <- NA_real_                # align both sides to complete pairs
  b_cols[!mask] <- NA_real_
  constant <- suppressWarnings(
    matrixStats::rowMins(a_cols, na.rm = TRUE) ==
      matrixStats::rowMaxs(a_cols, na.rm = TRUE) |
      matrixStats::rowMins(b_cols, na.rm = TRUE) ==
        matrixStats::rowMaxs(b_cols, na.rm = TRUE)
  )
  a_cols[!mask] <- 0
  b_cols[!mask] <- 0
  k   <- rowSums(mask)
  sa  <- rowSums(a_cols)
  sb  <- rowSums(b_cols)
  saa <- rowSums(a_cols * a_cols)
  sbb <- rowSums(b_cols * b_cols)
  sab <- rowSums(a_cols * b_cols)
  num <- sab - sa * sb / k
  den <- sqrt(pmax(saa - sa * sa / k, 0) * pmax(sbb - sb * sb / k, 0))
  value <- num / den
  value[k <= 2L | constant | !is.finite(value)] <- NA_real_
  as.numeric(value)
}
