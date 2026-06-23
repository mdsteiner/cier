#' Person-total correlation (r_pbis) C/IER index
#'
#' Computes each respondent's person-total correlation -- the Pearson correlation
#' between their item responses and the per-item means across the whole sample (the
#' "group-determined item difficulties" of Donlon & Fischer, 1968). A careless
#' pattern decouples from the group profile, driving the correlation toward zero, so
#' low values flag carelessness. Computed over the raw response row, so unlike the
#' split-half consistency indices it also works on a single scale.
#'
#' @details
#' Item-total form (per-item mean taken over all respondents, including the scored
#' one): numerically near-identical to the item-rest leave-one-out variant on real
#' data and far cheaper to vectorise. Keying-insensitive because the per-item means
#' are empirical, so it takes no item metadata.
#'
#' The default flags the lowest-correlation respondents at the `fpr` quantile (a
#' ranking convention that flags at least `fpr` of respondents, not a calibrated
#' false-positive rate). Adjust with `fpr`, or pass an absolute `cutoff` in `[-1, 1]`.
#'
#' A respondent answering fewer than three items, or whose answered responses have
#' zero variance (a straightliner), has no defined correlation (`value` and `flagged`
#' are `NA`). When all item means are equal (e.g. a perfectly balanced design) every
#' respondent abstains with the standard insufficient-scores warning.
#'
#' @section What this catches:
#' Random and partial / changepoint responding and idiosyncratic patterns (diagonal,
#' alternating) whose profile diverges from the sample's item-mean profile. It misses
#' pure straightlining and midpoint locking (a constant row abstains); pair it with
#' [cier_longstring()] and [cier_irv()].
#'
#' @template responses
#' @param fpr Optional lower-tail mass for the percentile cutoff, in `(0, 1)`; `NULL`
#'   (default) uses `0.05`. Mutually exclusive with `cutoff`.
#' @param cutoff Optional literal cutoff on the correlation, in `[-1, 1]`; respondents
#'   at or below it are flagged. Mutually exclusive with `fpr`.
#'
#' @template return-cier-index
#'
#' @references
#' Donlon, T. F., & Fischer, F. E. (1968). An index of an individual's agreement
#' with group-determined item difficulties. *Educational and Psychological
#' Measurement*, 28(1), 105–113. \doi{10.1177/001316446802800110}
#'
#' Curran, P. G. (2016). Methods for the detection of carelessly invalid
#' responses in survey data. *Journal of Experimental Social Psychology*, 66,
#' 4–19. \doi{10.1016/j.jesp.2015.07.006}
#'
#' @seealso [cier_mahalanobis()], [cier_irv()]
#' @family indirect indices
#' @export
#' @examples
#' # The 44 BFI items are the first 44 columns of the example data.
#' out <- cier_person_total(bfi_careless[, 1:44])
#' out
#' head(as.data.frame(out))
cier_person_total <- function(responses, fpr = NULL, cutoff = NULL) {
  call <- rlang::caller_env()
  responses <- check_responses(responses, call = call)
  # `fpr` is a tail mass in (0, 1); a literal `cutoff` a correlation in [-1, 1].
  check_percentile_overrides(fpr, cutoff, lower = -1, upper = 1, call = call)
  row <- cier_method_spec("cier_person_total")
  # Kernel returns NA for <3 answered items or zero variance, so no abstention guard.
  value <- kernel_person_total(responses)
  resolve_index_cutoff(value, row, fpr, cutoff, call = call)
}
