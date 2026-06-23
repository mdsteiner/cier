#' Psychometric synonyms (within-person synonym consistency) C/IER index
#'
#' Each respondent's psychometric-synonyms score is the Pearson correlation, within
#' that respondent, across the item pairs that are strongly positively correlated in
#' the whole sample (the "synonyms"). Careless answers decouple near-synonymous items,
#' driving the correlation toward zero, so low values flag carelessness (Meade &
#' Craig, 2012).
#'
#' @details
#' Synonym pairs are item pairs whose whole-sample inter-item correlation exceeds
#' `critical_r` (default `0.60`). Pairing uses raw responses with no reverse-scoring,
#' so a reverse-keyed item (negative correlations) is never a synonym (it may be an
#' antonym instead). The score is the single correlation across all qualifying pairs'
#' stacked values.
#'
#' `critical_r` is a property of the inventory: on broad personality inventories the
#' strongest inter-item correlations can fall below `0.60`, so no pairs qualify, every
#' value is `NA`, and the index abstains with a warning. Lower `critical_r`, or use a
#' more homogeneous scale.
#'
#' The default cutoff flags the lowest-consistency respondents at the `fpr` quantile
#' (5th percentile), flagging at least `fpr` of respondents rather than a calibrated
#' false-positive rate. A respondent with fewer than three complete synonym pairs, or
#' zero-variance pair values (a straightliner), has no defined correlation
#' (`value`/`flagged` are `NA`).
#'
#' @section What this catches:
#' Random and partial / changepoint responding. It misses pure straightlining and
#' midpoint locking (a constant row abstains); pair it with [cier_longstring()] and
#' [cier_irv()].
#'
#' @section Pairing reference on contaminated samples:
#' When much of the sample responds carelessly the inter-item correlations shrink and
#' no pair clears `critical_r`, so the index abstains exactly when carelessness is
#' highest. Supplying `reference` discovers pairs on a clean source while every
#' respondent is still scored on the full analysis sample. Because `critical_r = 0.60`
#' is strict, this workflow usually supplies a clean `reference` and relaxes
#' `critical_r` (e.g. to `0.40`). `reference = NULL` (default) is whole-sample
#' self-pairing.
#'
#' @template responses
#' @param critical_r Inter-item correlation a pair must exceed to count as a synonym;
#'   a number in `(0, 1)`, default `0.60` (Meade & Craig, 2012).
#' @param fpr Optional lower-tail mass for the percentile cutoff; `NULL` (default) uses
#'   `0.05`, else a number in `(0, 1)`. Mutually exclusive with `cutoff`.
#' @param cutoff Optional literal cutoff on the correlation in `[-1, 1]`; respondents
#'   at or below it are flagged. Mutually exclusive with `fpr`.
#' @param reference Optional clean source for pair discovery (scoring always uses the
#'   full `responses`). `NULL` (default) discovers pairs on `responses`. Otherwise a
#'   logical mask (length `nrow(responses)`) or integer row indices selecting a clean
#'   subset, or a separate clean sample (numeric matrix / data.frame, same items, >= 3
#'   rows, aligned by column name when both are named else positionally). See the
#'   *Pairing reference on contaminated samples* section.
#'
#' @template return-cier-index
#'
#' @references
#' Meade, A. W., & Craig, S. B. (2012). Identifying careless responses in survey
#' data. *Psychological Methods*, 17(3), 437–455. \doi{10.1037/a0028085}
#'
#' @seealso [cier_person_total()], [cier_irv()]
#' @family indirect indices
#' @export
#' @examples
#' # BFI items (first 44 columns); lower critical_r since BFI correlations fall
#' # below the 0.60 default.
#' out <- cier_psychsyn(bfi_careless[, 1:44], critical_r = 0.5)
#' out
#' head(as.data.frame(out))
#'
#' # Discover pairs on a clean reference subset while scoring on the full sample
#' # (contaminated-sample workflow; relax critical_r too).
#' out2 <- cier_psychsyn(bfi_careless[, 1:44], critical_r = 0.4,
#'                       reference = 1:200)
cier_psychsyn <- function(responses, critical_r = 0.60, fpr = NULL,
                          cutoff = NULL, reference = NULL) {
  call <- rlang::caller_env()
  responses <- check_responses(responses, call = call)
  check_open_unit(critical_r, "critical_r", call = call)
  # `fpr` is a tail mass in (0, 1); a literal `cutoff` a correlation in [-1, 1].
  check_percentile_overrides(fpr, cutoff, lower = -1, upper = 1, call = call)
  row <- cier_method_spec("cier_psychsyn")
  # `reference` (default NULL = self-pairing) moves only pair DISCOVERY to a clean
  # sample; scoring stays on full `responses`. When no pair qualifies, the shared tail
  # swaps the generic percentile abstention for the actionable no-pairs warning.
  cor_mat <- resolve_pairing_cor(responses, reference, call)
  pairs <- find_item_pairs(responses, critical_r, "syn", cor_mat = cor_mat)
  value <- kernel_psychsyn(responses, critical_r, "syn", pairs = pairs)
  resolve_pair_index_cutoff(value, row, fpr, cutoff, nrow(pairs) == 0L, "syn",
                            critical_r, cor_mat, call)
}
