#' Intra-individual response variability (IRV) C/IER index
#'
#' Each respondent's IRV is the sample SD of their answers across items
#' (Marjanovic et al., 2015), computed over the raw response row (no scale
#' blocking). Low values flag carelessness.
#'
#' @details
#' The default flags the lowest-variability respondents at the `fpr` quantile (5th
#' percentile), a ranking convention that flags at least `fpr` of respondents, not a
#' calibrated false-positive rate. Blank items are dropped; a respondent with fewer
#' than two answered items has no defined SD (`value` and `flagged` are `NA`), while a
#' constant respondent with at least two items scores `0` and is flagged.
#'
#' @section What this catches:
#' Straightlining and midpoint responding, which compress the response range. Not
#' useful against random responding, which inflates variability.
#'
#' @template responses
#' @param fpr Optional lower-tail mass in `(0, 1)` for the percentile cutoff; `NULL`
#'   (default) uses `0.05`. Mutually exclusive with `cutoff`.
#' @param cutoff Optional literal IRV cutoff; respondents at or below it are flagged.
#'   Mutually exclusive with `fpr`.
#'
#' @template return-cier-index
#'
#' @references
#' Marjanovic, Z., Holden, R., Struthers, W., Cribbie, R., & Greenglass, E.
#' (2015). The inter-item standard deviation (ISD): An index that discriminates
#' between conscientious and random responders. *Personality and Individual
#' Differences*, 84, 79–83. \doi{10.1016/j.paid.2014.08.021}
#'
#' @family indirect indices
#' @export
#' @examples
#' # The 44 BFI items are the first 44 columns of the bundled example data.
#' out <- cier_irv(bfi_careless[, 1:44])
#' out
#' head(as.data.frame(out))
cier_irv <- function(responses, fpr = NULL, cutoff = NULL) {
  call <- rlang::caller_env()
  responses <- check_responses(responses, call = call)
  check_percentile_overrides(fpr, cutoff, lower = 0, call = call)
  row <- cier_method_spec("cier_irv")
  # rowSds(na.rm = TRUE) returns NA below two answered items, so abstention needs
  # no separate guard.
  value <- kernel_irv(responses)
  resolve_index_cutoff(value, row, fpr, cutoff, call = call)
}
