#' Personal reliability (PR / RPR) C/IER index
#'
#' Computes each respondent's personal reliability. Within each scale block (given by
#' `items$scale`) the items are split into two halves, each half's mean is formed, and
#' the two half-mean vectors are correlated across scale blocks, Spearman-Brown
#' corrected (`2r / (1 + r)`, clamped at `-1`), and negated. A careful respondent's
#' halves correlate highly (low negated score); a careless respondent's decouple, so
#' high values flag carelessness.
#'
#' The two variants differ only in how each scale is halved:
#' - **PR** (Jackson, 1976; `resample = FALSE`): a single deterministic split, first
#'   half of each scale against the second.
#' - **RPR** (Goldammer et al., 2024; `resample = TRUE`, default, recommended): averages
#'   the statistic over `n_resamples` random within-scale splits.
#'
#' @details
#' **Reproducibility.** RPR draws random splits, so the default `seed = NULL` gives
#' different scores per call; pass an integer `seed` for a reproducible run (applied
#' locally, restoring the caller's RNG state, and assuming the same [RNGkind()]). PR is
#' deterministic and ignores `seed` and `n_resamples`.
#'
#' **Reverse-keying.** Reverse-keyed items (`items$reverse_keyed`) are reverse-scored
#' with `(min + max) - x` before the half-means, so the index sees trait-aligned
#' responses. `min` defaults to `1`; declare `items$min` for a 0-based or bipolar scale.
#'
#' **Cutoff.** The default flags at the `1 - fpr` quantile (95th percentile), a ranking
#' convention that flags at least `fpr` of respondents rather than a calibrated
#' false-positive rate. Adjust with `fpr`, or pass an absolute `cutoff` in `[-1, 1]`;
#' the two are mutually exclusive.
#'
#' **Score range.** The Spearman-Brown clamp maps every across-block correlation at or
#' below `-1/3` to exactly `+1`, giving the negated per-iteration score a point mass at
#' `+1` (PR returns it directly, RPR averages it). A percentile cutoff landing on the
#' atom flags more than `fpr` of respondents.
#'
#' **Abstention.** A respondent with fewer than two scorable scale blocks, or whose
#' half-mean vectors have no across-scale variance, has no defined correlation (`value`,
#' `flagged` are `NA`); for RPR, only when every iteration abstains. With exactly two
#' scorable blocks the correlation is `+1`/`-1` by construction -- a degenerate score
#' returned with a typed warning; `>= 3` multi-item scales are recommended.
#'
#' @section What this catches:
#' Random, alternating, and opposite-pattern responding whose within-scale halves
#' diverge. It misses pure straightlining and midpoint locking (the constant row
#' abstains) -- catch those with [cier_longstring()] and [cier_irv()].
#'
#' @template responses
#' @template items-scale-keyed
#' @param resample Logical; `TRUE` (default) for RPR, `FALSE` for deterministic PR.
#' @param n_resamples Positive integer; random split-halves averaged for RPR (default
#'   `25`). Unused when `resample = FALSE`.
#' @param seed Optional integer seed for RPR's splits; `NULL` (default) uses the ambient
#'   session RNG. Unused when `resample = FALSE`.
#' @param fpr Optional false-positive tail mass in `(0, 1)`; `NULL` (default) uses
#'   `0.05`, giving the upper-tail (`1 - fpr`) quantile. Mutually exclusive with
#'   `cutoff`.
#' @param cutoff Optional literal score cutoff in `[-1, 1]`; respondents at or above it
#'   are flagged. Mutually exclusive with `fpr`.
#'
#' @template return-cier-index
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
#' # Item metadata from BFI-44 column names: scale = letters after "v_BFI_",
#' # "_R" suffix = reverse-keyed, items coded 1..5.
#' nm <- names(bfi_careless)[1:44]
#' items <- data.frame(
#'   scale = sub("^v_BFI_([A-Za-z]+)[0-9].*$", "\\1", nm),
#'   reverse_keyed = grepl("_R$", nm),
#'   max = 5L
#' )
#' # Resampled personal reliability (default); seed for reproducibility.
#' out <- cier_personal_reliability(bfi_careless[, 1:44], items, seed = 1)
#' out
#' head(as.data.frame(out))
cier_personal_reliability <- function(responses, items, resample = TRUE,
                                      n_resamples = 25L, seed = NULL,
                                      fpr = NULL, cutoff = NULL) {
  call <- rlang::caller_env()
  responses <- check_responses(responses, call = call)
  inform_if_unkeyed(items, ncol(responses))
  items <- check_items(items, ncol(responses), min_scales = 2L, call = call,
                       response_names = colnames(responses))
  check_flag(resample, "resample", call = call)
  check_count(n_resamples, "n_resamples", call = call)
  if (!is.null(seed)) check_int(seed, "seed", call = call)
  # cutoff bounds: scores live in [-1, 1].
  check_percentile_overrides(fpr, cutoff, lower = -1, upper = 1, call = call)
  row <- cier_method_spec("cier_personal_reliability")
  responses <- apply_split_half_keying(responses, items, call = call)
  blocks <- scale_block_indices(items)
  warn_two_scale_consistency(blocks, call = call)
  value <- if (resample) {
    kernel_rpr(responses, blocks, n_resamples, seed)
  } else {
    kernel_personal_reliability(responses, blocks)
  }
  resolve_index_cutoff(value, row, fpr, cutoff, call = call)
}
