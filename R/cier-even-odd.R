# Purpose: cier_even_odd() -- the public even-odd consistency C/IER index.
# Args:    See documentation below.
# Returns: A light `cier_index` (see new_cier_index()).
# Invariants:
#   - Deterministic (no RNG). First metadata index: validates `items` via the
#     shared check_items() and reverse-scores keyed items via the shared
#     apply_split_half_keying() before the kernel.
#   - Bytewise compatible with careless::evenodd() on data with no reverse-keyed
#     items (careless does not reverse-key).
#   - The cutoff routes through the single resolve_cutoff() path (percentile
#     method, upper direction, fpr = 0.05 by default).

#' Even-odd consistency C/IER index
#'
#' Computes each respondent's **even-odd consistency**. Within each scale block
#' (given by `items$scale`) the mean of the even-positioned items and the mean of
#' the odd-positioned items are formed; the two are correlated across scale
#' blocks, Spearman-Brown corrected (`2r / (1 + r)`, clamped at `-1`), and
#' **negated**. A careful respondent answers the two halves of each scale
#' consistently, so the halves correlate highly and the negated score is low; a
#' careless respondent's halves decouple, so the score is high. **High** values
#' therefore flag carelessness.
#'
#' @details
#' **Reverse-keying.** Reverse-keyed items (`items$reverse_keyed`) are
#' reverse-scored with the self-inverse reflection `(min + max) - x` (with `max`
#' the largest response option, `items$max`) *before* the half-means are formed,
#' so the index always sees trait-aligned responses. The scale base `min`
#' defaults to `1` (the `1..max` coding), giving the classic `(max + 1) - x`;
#' declare `items$min` for a 0-based or bipolar scale so reverse items reflect
#' onto the same range. Supply the raw responses and declare reverse items
#' through `items`; pre-scoring them and leaving `reverse_keyed = FALSE` yields
#' the same result.
#'
#' **Cutoff.** The default flags the highest-inconsistency respondents: the cutoff
#' is the empirical `1 - fpr` quantile of the observed scores (the 95th percentile
#' by default) and respondents at or above it are flagged. This is a **ranking**
#' convention -- a sample percentile flags `fpr` of respondents by construction --
#' not a calibrated false-positive rate. Adjust the target tail with `fpr`, or
#' pass an absolute `cutoff` in `[-1, 1]` to flag on a literal score threshold
#' (e.g. one carried over from a calibration sample). `fpr` and `cutoff` are
#' mutually exclusive.
#'
#' **Abstention.** A respondent for whom fewer than two scale blocks yield a
#' finite even-and-odd mean pair, or whose (reverse-scored) half-mean vectors have
#' no variance across scales, has no defined consistency correlation: both `value`
#' and `flagged` are `NA` and the row is excluded from the flag count and rate. A
#' scale with only one item cannot be split and is skipped. Note that a constant
#' (raw straightliner) row abstains only when reverse-scoring leaves it constant --
#' a forward-keyed battery, or a value at the scale midpoint. With reverse-keyed
#' items an off-midpoint constant row is reflected into a non-constant one and so
#' receives a finite score rather than abstaining.
#'
#' @section What this catches:
#' Random, alternating, and opposite-pattern responding whose within-scale halves
#' diverge. On a forward-keyed (or scale-midpoint) battery it **misses** pure
#' straightlining and midpoint locking -- the constant row's even and odd means
#' are identical, so it abstains. With reverse-keyed items an off-midpoint
#' straightliner is instead reflected into a non-constant row and scored, but the
#' score is not a reliable careless signal, so straightlining is still best caught
#' directly with [cier_longstring()] and [cier_irv()].
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
#'   column gives the smallest response option (the scale base; default `1`, i.e.
#'   `1..max` coding) -- declare it for 0-based or bipolar scales.
#' @param fpr Optional target false-positive tail mass for the percentile
#'   cutoff. `NULL` (default) uses the registry default `0.05`. A finite number
#'   in the open interval `(0, 1)`; the cutoff is that upper-tail (`1 - fpr`)
#'   quantile of the observed scores. Mutually exclusive with `cutoff`.
#' @param cutoff Optional **literal** cutoff on the score, supplied instead of
#'   `fpr`. A single finite number in `[-1, 1]`; respondents whose even-odd score
#'   is at or above it are flagged (the upper-tail direction). Use it to apply an
#'   absolute threshold rather than a sample percentile.
#'
#' @return A `cier_index`: a list with per-respondent `value` (numeric, `NA` on
#'   abstention) and `flagged` (logical) vectors plus the `method`, `cutoff`, and
#'   `direction` metadata. Use [as.data.frame()][as.data.frame.cier_index] for a
#'   tidy data frame and `print()` for a summary.
#'
#' @references
#' Curran, P. G. (2016). Methods for the detection of carelessly invalid
#' responses in survey data. *Journal of Experimental Social Psychology*, 66,
#' 4–19. \doi{10.1016/j.jesp.2015.07.006}
#'
#' Jackson, D. N. (1976). The appraisal of personal reliability. Paper presented
#' at the meeting of the Society of Multivariate Experimental Psychology,
#' University Park, PA.
#'
#' @seealso [careless::evenodd()], [cier_longstring()], [cier_irv()]
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
#' out <- cier_even_odd(bfi_careless[, 1:44], items)
#' out
#' head(as.data.frame(out))
cier_even_odd <- function(responses, items, fpr = NULL, cutoff = NULL) {
  call <- rlang::caller_env()
  responses <- check_responses(responses, call = call)
  items <- check_items(items, ncol(responses), min_scales = 2L, call = call)
  # `fpr` is a tail mass in (0, 1); a literal `cutoff` is an even-odd score
  # threshold in [-1, 1]; the two are mutually exclusive.
  check_percentile_overrides(fpr, cutoff, lower = -1, upper = 1, call = call)
  row <- cier_method_row("cier_even_odd")
  # Reverse-score keyed items, then score the even-odd consistency. The kernel
  # returns NA where a row abstains, so abstention needs no separate guard.
  responses <- apply_split_half_keying(responses, items, call = call)
  value <- kernel_even_odd(responses, scale_block_indices(items))
  resolve_index_cutoff(value, row, fpr, cutoff, call = call)
}
