#' Even-odd consistency C/IER index
#'
#' Computes each respondent's even-odd consistency. Within each scale block (given
#' by `items$scale`) the even-position and odd-position item means are formed,
#' correlated across scale blocks, Spearman-Brown corrected (`2r / (1 + r)`, clamped
#' at `-1`), and negated. A careful respondent's halves correlate highly (low
#' negated score); a careless respondent's decouple, so high values flag
#' carelessness.
#'
#' @details
#' Reverse-keyed items (`items$reverse_keyed`) are reverse-scored with `(min + max) - x`
#' before the half-means, so the index sees trait-aligned responses. `min` defaults to
#' `1`; declare `items$min` for a 0-based or bipolar scale.
#'
#' The default cutoff flags the upper `1 - fpr` quantile (95th percentile) -- a ranking
#' convention flagging at least `fpr` of respondents, not a calibrated false-positive
#' rate. The Spearman-Brown clamp maps every across-block correlation `<= -1/3` to `+1`,
#' so the negated score has a point mass at `+1`; a percentile cutoff landing on that
#' atom flags more than `fpr`.
#'
#' Abstention (`value` and `flagged` are `NA`): fewer than two scorable scale blocks, or
#' half-mean vectors with no across-scale variance (a straightliner row abstains only if
#' reverse-scoring leaves it constant). With exactly two scorable blocks the correlation
#' is `+1`/`-1` by construction -- a degenerate binary score returned with a typed
#' warning; `>= 3` multi-item scales are recommended.
#'
#' @section What this catches:
#' Random, alternating, and opposite-pattern responding whose within-scale halves
#' diverge. It misses pure straightlining and midpoint locking (the constant row
#' abstains) -- catch those with [cier_longstring()] and [cier_irv()].
#'
#' @template responses
#' @template items-scale-keyed
#' @param fpr Optional tail mass in `(0, 1)` for the upper-tail (`1 - fpr`) quantile
#'   cutoff; `NULL` (default) uses `0.05`. Mutually exclusive with `cutoff`.
#' @param cutoff Optional literal cutoff in `[-1, 1]`; respondents at or above it are
#'   flagged.
#'
#' @template return-cier-index
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
#' @seealso [cier_longstring()], [cier_irv()]
#' @family indirect indices
#' @export
#' @examples
#' # Item metadata from BFI-44 column names: scale = letter run after "v_BFI_",
#' # trailing "_R" = reverse-keyed, items coded 1..5.
#' nm <- names(bfi_careless)[1:44]
#' items <- data.frame(
#'   scale = sub("^v_BFI_([A-Za-z]+)[0-9].*$", "\\1", nm),
#'   reverse_keyed = grepl("_R$", nm),
#'   max = 5L
#' )
#' # Many respondents tie at the +1 clamp maximum, so the default 95th-percentile
#' # cutoff lands on it and flags >5% -- an expected `cier_warning_saturated_cutoff`.
#' # Set an explicit `cutoff` to avoid it.
#' out <- cier_even_odd(bfi_careless[, 1:44], items)
#' out
#' head(as.data.frame(out))
cier_even_odd <- function(responses, items, fpr = NULL, cutoff = NULL) {
  call <- rlang::caller_env()
  responses <- check_responses(responses, call = call)
  inform_if_unkeyed(items, ncol(responses))
  items <- check_items(items, ncol(responses), min_scales = 2L, call = call,
                       response_names = colnames(responses))
  # `fpr` is a tail mass in (0, 1); a literal `cutoff` is a score in [-1, 1].
  check_percentile_overrides(fpr, cutoff, lower = -1, upper = 1, call = call)
  row <- cier_method_spec("cier_even_odd")
  # Reverse-score keyed items before scoring; kernel returns NA on abstention.
  responses <- apply_split_half_keying(responses, items, call = call)
  blocks <- scale_block_indices(items)
  warn_two_scale_consistency(blocks, call = call)
  value <- kernel_even_odd(responses, blocks)
  resolve_index_cutoff(value, row, fpr, cutoff, call = call)
}
