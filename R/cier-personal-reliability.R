# Purpose: cier_personal_reliability() -- the public personal-reliability C/IER
#          index, in both the classical (PR; Jackson 1976) and resampled (RPR;
#          Goldammer et al. 2024) variants.
# Args:    See documentation below.
# Returns: A light `cier_index` (see new_cier_index()).
# Invariants:
#   - Second metadata index: validates `items` via the shared check_items() and
#     reverse-scores keyed items via apply_split_half_keying() before the kernel.
#   - PR is deterministic; RPR is reproducible only when `seed` is supplied
#     (a non-NULL seed is applied locally and the caller's RNG state restored, so
#     a seeded call does not disturb the session RNG stream).
#   - The cutoff routes through the single resolve_cutoff() path (percentile
#     method, upper direction, fpr = 0.05 by default).

#' Personal reliability (PR / RPR) C/IER index
#'
#' Computes each respondent's **personal reliability**. Within each scale block
#' (given by `items$scale`) the items are split into two halves; the mean of each
#' half is formed, the two half-mean vectors are correlated across scale blocks,
#' Spearman-Brown corrected (`2r / (1 + r)`, clamped at `-1`), and **negated**. A
#' careful respondent answers the two halves of each scale consistently, so the
#' halves correlate highly and the negated score is low; a careless respondent's
#' halves decouple, so the score is high. **High** values therefore flag
#' carelessness.
#'
#' Two variants share this definition and differ only in how each scale is
#' halved:
#' - **PR** (Jackson, 1976; `resample = FALSE`) uses a single deterministic
#'   split: the first half of each scale's items against the second half.
#' - **RPR** (Goldammer et al., 2024; `resample = TRUE`, the default) averages
#'   the statistic over `n_resamples` **random** within-scale splits. Resampling
#'   removes the arbitrariness of a single fixed split and is the recommended
#'   reliability index.
#'
#' @details
#' **Reproducibility (RPR).** RPR draws random splits, so its result depends on
#' the random-number stream. With the default `seed = NULL` the ambient session
#' RNG is used and **repeated calls give different scores**; pass an integer
#' `seed` to make a run reproducible. A non-`NULL` `seed` is applied locally and
#' the previous RNG state is restored on return, so a seeded call reproduces
#' **without** disturbing the caller's random stream. PR (`resample = FALSE`) is
#' deterministic and ignores `seed` and `n_resamples`.
#'
#' **Reverse-keying.** Reverse-keyed items (`items$reverse_keyed`) are
#' reverse-scored with the self-inverse reflection `(min + max) - x` (with `max`
#' the largest response option, `items$max`) *before* the half-means are formed,
#' so the index always sees trait-aligned responses. The scale base `min`
#' defaults to `1` (the `1..max` coding), giving the classic `(max + 1) - x`;
#' declare `items$min` for a 0-based or bipolar scale. Supply the raw responses
#' and declare reverse items through `items`; pre-scoring them and leaving
#' `reverse_keyed = FALSE` yields the same result.
#'
#' **Cutoff.** The default flags the lowest-reliability respondents: the cutoff
#' is the empirical `1 - fpr` quantile of the observed scores (the 95th
#' percentile by default) and respondents at or above it are flagged. This is a
#' **ranking** convention -- a sample percentile flags **at least** `fpr` of
#' respondents (more when scores tie at the cutoff) -- not a calibrated
#' false-positive rate. Adjust the target tail
#' with `fpr`, or pass an absolute `cutoff` in `[-1, 1]` to flag on a literal
#' score threshold (e.g. one carried over from a calibration sample). `fpr` and
#' `cutoff` are mutually exclusive.
#'
#' **Abstention.** A respondent for whom fewer than two scale blocks yield a
#' finite half-mean pair, or whose half-mean vectors have no variance across
#' scales, has no defined reliability correlation: both `value` and `flagged` are
#' `NA` and the row is excluded from the flag count and rate. A scale with only
#' one item cannot be split and is skipped. For RPR a respondent is `NA` only
#' when **every** resampled iteration abstains; otherwise the score averages the
#' finite iterations.
#'
#' @section What this catches:
#' Random, alternating, and opposite-pattern responding whose within-scale halves
#' diverge. On a forward-keyed (or scale-midpoint) battery it **misses** pure
#' straightlining and midpoint locking -- the constant row's half-means are
#' identical, so it abstains -- so straightlining is best caught directly with
#' [cier_longstring()] and [cier_irv()].
#'
#' @param responses A numeric matrix (or a data.frame / tibble coerced
#'   internally) of responses, one row per respondent and one column per item.
#'   `NA` marks a missing response.
#' @param items A data.frame of item metadata, one row per item, aligned to the
#'   columns of `responses`. Must carry a `scale` column with at least two
#'   distinct scale labels. An optional logical `reverse_keyed` column marks
#'   reverse-keyed items (default: none); when any item is reverse-keyed an
#'   integer `max` column (the largest response option; at least `min + 1`) is
#'   required so those items can be reverse-scored. An optional integer `min`
#'   column gives the smallest response option (the scale base; default `1`).
#' @param resample Logical: `TRUE` (default) computes the resampled variant
#'   (RPR); `FALSE` computes the deterministic classical variant (PR).
#' @param n_resamples Positive integer: the number of random split-halves
#'   averaged for RPR (default `25`). Validated but otherwise unused when
#'   `resample = FALSE`.
#' @param seed Optional integer seed for RPR's random splits. `NULL` (default)
#'   draws from the ambient session RNG (results are not reproducible across
#'   calls); supplying a seed makes a run reproducible without disturbing the
#'   caller's RNG stream. Validated but otherwise unused when `resample = FALSE`.
#' @param fpr Optional target false-positive tail mass for the percentile
#'   cutoff. `NULL` (default) uses the registry default `0.05`. A finite number
#'   in the open interval `(0, 1)`; the cutoff is that upper-tail (`1 - fpr`)
#'   quantile of the observed scores. Mutually exclusive with `cutoff`.
#' @param cutoff Optional **literal** cutoff on the score, supplied instead of
#'   `fpr`. A single finite number in `[-1, 1]`; respondents whose score is at or
#'   above it are flagged (the upper-tail direction). Use it to apply an absolute
#'   threshold rather than a sample percentile.
#'
#' @return A `cier_index`: a list with per-respondent `value` (numeric, `NA` on
#'   abstention) and `flagged` (logical) vectors plus the `method`, `cutoff`, and
#'   `direction` metadata. Use [as.data.frame()][as.data.frame.cier_index] for a
#'   tidy data frame and `print()` for a summary.
#'
#' @references
#' Jackson, D. N. (1976). The appraisal of personal reliability. Paper presented
#' at the meeting of the Society of Multivariate Experimental Psychology,
#' University Park, PA.
#'
#' Goldammer, P., Stöckli, P. L., Annen, H., & Schmitz-Wilhelmy, A. (2024). A
#' comparison of conventional and resampled personal reliability in detecting
#' careless responding. *Behavior Research Methods*, 56, 8831–8851.
#' \doi{10.3758/s13428-024-02506-0}
#'
#' @seealso [cier_even_odd()], [cier_longstring()], [cier_irv()]
#' @family indirect indices
#' @export
#' @examples
#' # Build item metadata from the BFI-44 column names: the scale is the letter
#' # run after "v_BFI_", a trailing "_R" marks a reverse-keyed item, and the
#' # items are coded 1..5.
#' nm <- names(bfi_careless)[1:44]
#' items <- data.frame(
#'   scale = sub("^v_BFI_([A-Za-z]+)[0-9].*$", "\\1", nm),
#'   reverse_keyed = grepl("_R$", nm),
#'   max = 5L
#' )
#' # Resampled personal reliability (the default); set a seed for reproducibility.
#' out <- cier_personal_reliability(bfi_careless[, 1:44], items, seed = 1)
#' out
#' head(as.data.frame(out))
cier_personal_reliability <- function(responses, items, resample = TRUE,
                                      n_resamples = 25L, seed = NULL,
                                      fpr = NULL, cutoff = NULL) {
  call <- rlang::caller_env()
  responses <- check_responses(responses, call = call)
  items <- check_items(items, ncol(responses), min_scales = 2L, call = call)
  check_flag(resample, "resample", call = call)
  check_count(n_resamples, "n_resamples", call = call)
  if (!is.null(seed)) check_int(seed, "seed", call = call)
  # `fpr` is a tail mass in (0, 1); a literal `cutoff` is a score threshold in
  # [-1, 1]; the two are mutually exclusive.
  check_percentile_overrides(fpr, cutoff, lower = -1, upper = 1, call = call)
  row <- cier_method_row("cier_personal_reliability")
  responses <- apply_split_half_keying(responses, items, call = call)
  blocks <- scale_block_indices(items)
  value <- if (resample) {
    kernel_rpr(responses, blocks, n_resamples, seed)
  } else {
    kernel_personal_reliability(responses, blocks)
  }
  resolve_index_cutoff(value, row, fpr, cutoff, call = call)
}
