# cier_psychant(): psychometric-antonyms index. Shares the psychsyn kernel with the
# "ant" tail. Pairs on raw responses (reverse-keying would collapse the antonym
# correlations). Percentile cutoff, upper direction.

#' Psychometric antonyms (within-person antonym opposition) C/IER index
#'
#' Each respondent's psychometric-antonyms score is the within-respondent Pearson
#' correlation across the item pairs that are strongly negatively correlated in the
#' whole sample (the "antonyms"). A careful respondent answers genuine antonyms (e.g.
#' a trait item and its reverse-keyed twin) in opposition; a careless respondent's
#' answers decouple toward zero, so high values flag carelessness -- the mirror of
#' [cier_psychsyn()] (Meade & Craig, 2012).
#'
#' @details
#' Antonym pairs are item pairs whose whole-sample inter-item correlation is below
#' `-critical_r`. Pairing uses raw responses (no reverse-scoring, which would collapse
#' those negative correlations toward zero); the score is the single correlation across
#' all qualifying pairs' stacked values. Supply `critical_r` as a positive magnitude in
#' `(0, 1)` applied to the negative tail, matching [cier_psychsyn()] and the pairing
#' diagnostics ([cier_synonym_pairs()], [cier_psychsyn_critval()] with `antonym = TRUE`)
#' so one threshold means the same pair strength for synonyms and antonyms; a negative
#' value is an input error.
#'
#' `critical_r` is a property of the inventory: if no two items oppose each other
#' strongly enough, no pairs qualify, every value is `NA`, and the index abstains with a
#' warning -- lower `critical_r`, or use a scale pairing forward and reverse items
#' ([cier_synonym_pairs()] with `antonym = TRUE` lists qualifying pairs). A respondent
#' with fewer than three complete pairs, or zero-variance pair values (a straightliner),
#' also has `NA` `value` and `flagged`.
#'
#' The default cutoff flags the lowest-opposition respondents at the `1 - fpr` quantile
#' (95th percentile). This ranking flags at least `fpr` of respondents, not a calibrated
#' false-positive rate; adjust with `fpr`, or pass an absolute `cutoff` in `[-1, 1]`
#' (mutually exclusive).
#'
#' @section What this catches:
#' Random and partial / changepoint responding. It misses pure straightlining and
#' midpoint locking (a constant row abstains); pair it with [cier_longstring()] and
#' [cier_irv()]. It also needs genuine antonym pairs to exist -- usually forward/reverse
#' item pairs.
#'
#' @section Pairing reference on contaminated samples:
#' When much of the sample responds carelessly the antonym correlations shrink toward
#' zero and no pair clears `-critical_r`, so the index abstains exactly when carelessness
#' is highest. `reference` estimates the pairs on a clean source -- a screened subset of
#' the analysis rows, or a separate clean sample -- while every respondent is still
#' scored on the full sample. Such workflows usually supply a clean `reference` and relax
#' the strict default `critical_r` (e.g. to `0.40`). `reference = NULL` (default) uses
#' whole-sample self-pairing.
#'
#' @template responses
#' @param critical_r Antonym-pair threshold: positive magnitude in `(0, 1)` applied to
#'   the negative tail (a pair must have inter-item r below `-critical_r`). Default
#'   `0.60` follows Meade & Craig (2012).
#' @param fpr Optional upper-tail mass for the percentile cutoff; a finite number in
#'   `(0, 1)`. `NULL` (default) uses `0.05`. Mutually exclusive with `cutoff`.
#' @param cutoff Optional literal correlation cutoff, a finite number in `[-1, 1]`;
#'   respondents at or above it are flagged. Mutually exclusive with `fpr`.
#' @param reference Optional clean source for pair discovery only (scoring always uses
#'   the full `responses`). `NULL` (default) self-pairs on `responses`. Otherwise a
#'   logical mask (length `nrow(responses)`) or integer row indices selecting a clean
#'   subset, or a separate clean sample (numeric matrix / data.frame, same items, at
#'   least 3 rows; aligned by column name when both carry names, else positionally).
#'   See the *Pairing reference on contaminated samples* section.
#'
#' @template return-cier-index
#'
#' @references
#' Meade, A. W., & Craig, S. B. (2012). Identifying careless responses in survey
#' data. *Psychological Methods*, 17(3), 437–455. \doi{10.1037/a0028085}
#'
#' @seealso [cier_psychsyn()], [cier_synonym_pairs()], [cier_irv()]
#' @family indirect indices
#' @export
#' @examples
#' # The BFI's strongest oppositions are milder than the default 0.60 (most negative
#' # inter-item r is about -0.48), so lower critical_r to surface antonym pairs.
#' out <- cier_psychant(bfi_careless[, 1:44], critical_r = 0.40)
#' out
#' head(as.data.frame(out))
#'
#' # Discover pairs on a clean reference subset while scoring every respondent:
#' out2 <- cier_psychant(bfi_careless[, 1:44], critical_r = 0.40,
#'                       reference = 1:200)
cier_psychant <- function(responses, critical_r = 0.60, fpr = NULL,
                          cutoff = NULL, reference = NULL) {
  call <- rlang::caller_env()
  responses <- check_responses(responses, call = call)
  check_open_unit(critical_r, "critical_r", call = call)
  # `cutoff` is a literal correlation in [-1, 1].
  check_percentile_overrides(fpr, cutoff, lower = -1, upper = 1, call = call)
  row <- cier_method_spec("cier_psychant")
  # One correlation serves discovery and scoring; the "ant" tail keeps pairs with
  # r < -critical_r, and `reference` (NULL = self-pairing) moves only DISCOVERY to a
  # clean sample. When no pair qualifies, the shared tail emits the no-pairs warning.
  cor_mat <- resolve_pairing_cor(responses, reference, call)
  pairs <- find_item_pairs(responses, critical_r, "ant", cor_mat = cor_mat)
  value <- kernel_psychsyn(responses, critical_r, "ant", pairs = pairs)
  resolve_pair_index_cutoff(value, row, fpr, cutoff, nrow(pairs) == 0L, "ant",
                            critical_r, cor_mat, call)
}
