# Numerical kernels for the split-half consistency family -- even-odd, personal
# reliability (PR), and resampled personal reliability (RPR) -- plus the reverse-keying
# and declared-range helpers shared with cier_gnormed and cier_ht.
#
# Kernels return numeric vectors and never raise typed conditions (the wrappers validate);
# the exceptions are apply_split_half_keying() (aborts on a bad reverse-keying declaration)
# and warn_two_scale_consistency() (warns on a degenerate two-scale design). Kernels are
# pure, except kernel_rpr(): with a non-NULL seed it draws under a local seed and restores
# the caller's `.Random.seed` on exit.

# ---- Even-odd / split-half family -------------------------------------------

# Group column indices by scale label, in first-appearance order, so a
# heterogeneously-ordered battery yields a reproducible block sequence.
# `items$scale` is aligned to the columns of `responses`.
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

# Warn when exactly two scale blocks are scorable (>= 2 items each; a one-item block is
# skipped downstream). Across two block means the consistency correlation is +/-1 by
# construction: a degenerate point mass for even-odd / PR, a coarse average of +/-1
# iterations for RPR. A structural property of the item design, so the warning fires once
# at the wrapper, before scoring, without changing the statistic.
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

# Reduce the per-respondent first-/second-half scale-mean matrices (each n x n_blocks, NA
# column where a scale yields no half-mean) to the negated, Spearman-Brown-corrected
# split-half consistency score. The shared scoring tail of even-odd, PR, and every RPR
# iteration: per respondent, the across-scale Pearson correlation of first- vs second-half
# means over pairwise-complete blocks, Spearman-Brown corrected (2r/(1+r), clamped at -1)
# and negated (high = careless); NA when fewer than two blocks pair or either side has zero
# variance.
#
# Vectorised masked-sum Pearson: one set of rowSums scores every respondent at once. Two
# details reproduce stats::cor() exactly:
#   * a constant (zero-variance) complete-profile side is detected by masked min == max and
#     forced to NA; the deviation sum-of-squares alone can land a few ulp either side of 0
#     on a non-integer constant and leak a spurious finite score;
#   * r is clamped to [-1, 1] before Spearman-Brown, since the masked-sum form can drift a
#     few ulp past -1 on a two-block (+/-1) profile and r = -1-eps would explode 2r/(1+r).
split_half_pearson_neg <- function(first_means, second_means) {
  both <- !is.na(first_means) & !is.na(second_means)
  k <- rowSums(both)
  a_masked <- first_means
  a_masked[!both] <- NA_real_
  b_masked <- second_means
  b_masked[!both] <- NA_real_
  constant <- suppressWarnings(
    matrixStats::rowMins(a_masked, na.rm = TRUE) ==
      matrixStats::rowMaxs(a_masked, na.rm = TRUE) |
      matrixStats::rowMins(b_masked, na.rm = TRUE) ==
        matrixStats::rowMaxs(b_masked, na.rm = TRUE)
  )
  a <- first_means
  a[!both] <- 0
  b <- second_means
  b[!both] <- 0
  sa  <- rowSums(a)
  sbv <- rowSums(b)
  saa <- rowSums(a * a)
  sbb <- rowSums(b * b)
  sab <- rowSums(a * b)
  num <- sab - sa * sbv / k
  den <- sqrt(pmax(saa - sa * sa / k, 0) * pmax(sbb - sbv * sbv / k, 0))
  r <- num / den
  r[k < 2L | constant | !is.finite(r)] <- NA_real_
  r <- pmax(pmin(r, 1), -1)
  corrected <- (2 * r) / (1 + r)
  clamp <- !is.na(r) & (is.na(corrected) | corrected < -1)
  corrected[clamp] <- -1
  -corrected
}

# Cross-check the declared range against observed responses for the columns in `in_scope` (a
# logical mask). A type-valid but wrong declaration silently corrupts the score: a 0-based
# 0..4 scale declared max = 5 reflects a reverse item off-scale; a forward item with a
# 99-style sentinel shifts its half-mean. An all-NA column reduces to Inf / -Inf, so it is
# not an offender. Caller decides scope (every item for split-half / Gnormed; reverse items
# only for Ht).
#
# Bounds are gated independently: the upper bound for every in-scope item declaring a `max`,
# the lower bound additionally needs a non-NA `min`. No `max` declared -> exempt entirely.
assert_items_in_range <- function(responses, mins, maxs, in_scope, call) {
  hi <- in_scope & !is.na(maxs)            # upper bound: in scope and a max declared
  lo <- hi & !is.na(mins)                  # lower bound: also a min declared
  if (!any(hi)) {
    return(invisible(NULL))
  }
  obs_min <- suppressWarnings(matrixStats::colMins(responses, na.rm = TRUE))
  obs_max <- suppressWarnings(matrixStats::colMaxs(responses, na.rm = TRUE))
  bad <- (hi & obs_max > maxs) | (lo & obs_min < mins)
  if (any(bad)) {
    offending <- which(bad)
    cier_abort(
      "cier_error_input",
      c("Item responses must lie within the declared scale range \\
         {.field [min, max]}.",
        "x" = "{cli::qty(length(offending))}Out-of-range item{?s}: \\
               {.val {offending}}.",
        "i" = "Check {.field max} / {.field min} against the data (a 0-based \\
               scale needs {.field min} = 0); a value outside the declared range \\
               corrupts the score."),
      data = list(arg = "items", observed = offending), call = call
    )
  }
  invisible(NULL)
}

# Reverse-score reverse-keyed items so the split-half family forms its half-means on
# trait-aligned responses. A per-method step, not a global rescore (longstring / irv /
# person-total need the raw responses), so only the split-half and person-fit wrappers call
# it. It range-checks the data (assert_items_in_range), then reflects each reverse item with
# the self-inverse (min + max) - x (min default 1). NA-preserving; returns `responses`
# unchanged on valid data with no reverse item. `forward_range = TRUE` (split-half / Gnormed)
# range-checks every item declaring a min/max; Ht passes FALSE because its forward `max` is
# unused metadata, so only reverse items' range is a contract there. The reverse max / min
# presence guards below are defensive backstops for direct callers (the wrapper path
# validates via check_items()).
apply_split_half_keying <- function(responses, items, call = rlang::caller_env(),
                                    forward_range = TRUE) {
  rk <- items$reverse_keyed
  mins <- items$min                       # scale base; default 1 for direct callers
  if (is.null(mins)) {
    mins <- rep(1L, length(rk))
  }
  maxs <- items$max
  if (is.null(maxs)) {
    maxs <- rep(NA_real_, length(rk))     # no declared max -> nothing to range-check
  }
  # Ht restricts the range check to reverse items (its forward `max` is unused);
  # the split-half family and Gnormed check every item.
  in_scope <- if (forward_range) rep(TRUE, length(rk)) else rk
  assert_items_in_range(responses, mins, maxs, in_scope, call)
  if (!any(rk)) {
    return(responses)
  }
  # Reverse items must declare a max (and min) to be reflected.
  if (anyNA(maxs[rk])) {
    offending <- which(rk & is.na(maxs))
    cier_abort(
      "cier_error_input",
      c("Reverse-keyed items need a known {.field max} to reverse-score.",
        "x" = "Reverse-keyed item(s) with no scale maximum: {.val {offending}}."),
      data = list(arg = "items$max", observed = offending), call = call
    )
  }
  if (anyNA(mins[rk])) {
    offending <- which(rk & is.na(mins))
    cier_abort(
      "cier_error_input",
      c("Reverse-keyed items need a known {.field min} (scale base) to reverse-score.",
        "x" = "Reverse-keyed item(s) with no scale base: {.val {offending}}."),
      data = list(arg = "items$min", observed = offending), call = call
    )
  }
  rev_cols <- responses[, rk, drop = FALSE]
  reflect  <- mins[rk] + maxs[rk]         # (min + max) self-inverse reflection
  responses[, rk] <- rep(reflect, each = nrow(responses)) - rev_cols
  responses
}

# Even-odd split: even within-scale positions form the "first" half, odd the "second".
# cor() is symmetric, so this reproduces the even-vs-odd correlation exactly. Only called
# with k >= 2 (one-item scales yield a NULL split in compute_block_splits).
even_odd_split_fn <- function(k) {
  pos <- seq_len(k)
  list(first_idx  = pos[pos %% 2L == 0L],
       second_idx = pos[pos %% 2L == 1L])
}

# Precompute the within-scale split positions for every block, once outside the scoring pass:
# a per-block list aligned to `blocks`, each entry list(first_idx, second_idx) for a >= 2-item
# block or NULL for a one-item block (skipped downstream). Used by the deterministic even-odd
# / PR kernels; RPR draws its per-respondent random splits inline.
compute_block_splits <- function(blocks, split_fn) {
  lapply(blocks, function(cols) {
    if (length(cols) < 2L) NULL else split_fn(length(cols))
  })
}

# Row-wise mean of present values across a column submatrix, as masked sum / count: one
# scale's half-mean for every respondent in one pass. NA (never NaN) where a row has no
# present value, so unlike base mean(na.rm = TRUE) no is.nan() post-pass is needed.
masked_col_mean <- function(sub) {
  cnt <- rowSums(!is.na(sub))
  m <- rowSums(sub, na.rm = TRUE) / cnt
  m[cnt == 0L] <- NA_real_
  m
}

# Per-respondent split-half consistency over the scale blocks, all respondents in one pass.
# `splits` is the precomputed per-block list from compute_block_splits(). Each scored block's
# half-means come from masked_col_mean(); split_half_pearson_neg() reduces the two mean
# matrices to the score. Shared by even-odd and PR, reused by kernel_rpr() per iteration.
kernel_split_half <- function(responses, blocks, splits) {
  n <- nrow(responses)
  n_blocks <- length(blocks)
  first_means  <- matrix(NA_real_, n, n_blocks)
  second_means <- matrix(NA_real_, n, n_blocks)
  for (k in seq_len(n_blocks)) {
    split <- splits[[k]]
    if (is.null(split)) {
      next
    }
    cols <- blocks[[k]]
    first_means[, k]  <- masked_col_mean(responses[, cols[split$first_idx],
                                                   drop = FALSE])
    second_means[, k] <- masked_col_mean(responses[, cols[split$second_idx],
                                                   drop = FALSE])
  }
  split_half_pearson_neg(first_means, second_means)
}

# Build a deterministic split-half kernel: a function(responses, blocks) that precomputes the
# per-block split via `split_fn` and scores via kernel_split_half(). Even-odd and PR differ
# only in the deterministic split, so they share this factory. RPR is not built here: its
# random splits are drawn per respondent, so kernel_rpr() keeps its own draw loop.
make_split_half_kernel <- function(split_fn) {
  function(responses, blocks) {
    splits <- compute_block_splits(blocks, split_fn)
    kernel_split_half(responses, blocks, splits)
  }
}

# Even-odd consistency kernel: the even/odd within-scale split. The wrapper reverse-scores
# keyed items before calling and resolves the percentile abstention when every row is NA.
kernel_even_odd <- make_split_half_kernel(even_odd_split_fn)

# ---- Personal reliability (PR / RPR) ----------------------------------------

# First/second-half split for Personal Reliability (Jackson 1976): the first ceil(k / 2)
# within-scale positions form the "first" half, the remaining floor(k / 2) the "second".
# Only called with k >= 2 (one-item scales yield a NULL split and are skipped).
split_half_indices <- function(k) {
  first <- seq_len(ceiling(k / 2L))
  # The second half is the suffix after the first ceil(k/2) positions; an index
  # slice avoids a set difference (first is always a prefix of seq_len(k)).
  list(first_idx = first, second_idx = seq_len(k)[-first])
}

# Personal Reliability (PR) kernel: the deterministic first/second-half split.
kernel_personal_reliability <- make_split_half_kernel(split_half_indices)

# Resampled Personal Reliability (RPR) kernel (Goldammer et al. 2024): per respondent, the
# mean (na.rm = TRUE; NA when every iteration abstains) over `n_resamples` iterations of the
# PR statistic, each drawing a fresh uniform within-scale split (not the shared even-odd / PR
# split). The whole iteration -- draws, half-means, across-scale correlation -- is vectorised
# across respondents.
#
# Uniform split via rank-of-uniforms: per scale of `len` items, one n x len matrix of U(0, 1)
# draws is ranked per row (matrixStats::rowRanks) and the ceil(len / 2) smallest-draw
# positions form that respondent's first half. Ranking i.i.d. uniforms makes every size-ceil
# subset equally likely, so the split is uniform -- identical to a random permutation's first
# half -- while replacing n * n_resamples * scales sample.int() calls with one
# runif()/rowRanks() per scale per iteration (the dominant speedup).
#
# Reproducibility: the draw order (iteration -> scale, one runif(n * len) per >= 2-item scale)
# is fixed, and the seed is applied locally via with_local_seed() -- the caller's global RNG
# state is restored on exit -- so a seeded call reproduces without disturbing the caller's
# stream; a NULL seed draws from the ambient stream.
kernel_rpr <- function(responses, blocks, n_resamples, seed) {
  with_local_seed(seed, function() {
    n <- nrow(responses)
    n_blocks <- length(blocks)
    # Per-block quantities independent of the random split: the response submatrix, its
    # present mask, the NA -> 0 fill, the per-row half size, and the whole-block sum / count
    # (the second half is the block total minus the first half). One-item block -> NULL.
    blk <- lapply(blocks, function(cols) {
      len <- length(cols)
      if (len < 2L) {
        return(NULL)
      }
      vals <- responses[, cols, drop = FALSE]
      present <- !is.na(vals)
      filled <- vals
      filled[!present] <- 0
      list(len = len, half = ceiling(len / 2L), present = present,
           filled = filled, tot_sum = rowSums(filled), tot_cnt = rowSums(present))
    })
    scorable <- which(!vapply(blk, is.null, logical(1L)))
    per_iter <- matrix(NA_real_, nrow = n_resamples, ncol = n)
    for (b in seq_len(n_resamples)) {
      first_means  <- matrix(NA_real_, n, n_blocks)
      second_means <- matrix(NA_real_, n, n_blocks)
      for (k in scorable) {
        bk <- blk[[k]]
        # Random first-half membership: the ceil(len/2) smallest-uniform positions
        # per respondent.
        ranks <- matrixStats::rowRanks(
          matrix(stats::runif(n * bk$len), n, bk$len), ties.method = "first"
        )
        in_first <- ranks <= bk$half
        f_sum <- rowSums(bk$filled * in_first)
        f_cnt <- rowSums(bk$present & in_first)
        f_mean <- f_sum / f_cnt
        f_mean[f_cnt == 0L] <- NA_real_
        s_cnt <- bk$tot_cnt - f_cnt
        s_mean <- (bk$tot_sum - f_sum) / s_cnt
        s_mean[s_cnt == 0L] <- NA_real_
        first_means[, k]  <- f_mean
        second_means[, k] <- s_mean
      }
      per_iter[b, ] <- split_half_pearson_neg(first_means, second_means)
    }
    # Per-respondent mean over iterations; an all-NA column yields NaN -> NA.
    means <- matrixStats::colMeans2(per_iter, na.rm = TRUE)
    means[is.nan(means)] <- NA_real_
    means
  })
}
