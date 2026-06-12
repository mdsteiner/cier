# Purpose: cier_psychsyn() -- the public psychometric-synonyms C/IER index.
# Args:    See documentation below.
# Returns: A light `cier_index` (see new_cier_index()).
# Invariants:
#   - Deterministic (no RNG); whole-sample pairing on RAW responses, with no
#     reverse-keying (item keying plays no part in synonym pairing).
#   - Scores cor() on the full stacked pair vectors, matching
#     careless::psychsyn(resample_na = FALSE) bytewise.
#   - The cutoff routes through the single resolve_cutoff() path (percentile
#     method, lower direction, fpr = 0.05 by default).

#' Psychometric synonyms (within-person synonym consistency) C/IER index
#'
#' Computes each respondent's **psychometric-synonyms** score -- the Pearson
#' correlation, within that respondent, across the item pairs that are strongly
#' positively correlated in the whole sample (the "synonyms"). A careful
#' respondent answers near-synonymous items consistently, so their within-person
#' pattern tracks the sample's and the correlation is high; a careless
#' respondent's answers to synonyms decouple, driving the correlation toward
#' zero. **Low** values therefore flag carelessness (Meade & Craig, 2012).
#'
#' @details
#' **Pairing on raw responses.** The synonym pairs are the item pairs whose
#' whole-sample inter-item correlation exceeds `critical_r` (default `0.60`).
#' Pairing uses the **raw** responses, with no reverse-scoring: a reverse-keyed
#' item simply forms negative correlations and so is never paired as a synonym
#' (it may instead be an *antonym*; see the companion antonyms index). The score
#' is the single correlation across all qualifying pairs' stacked values, matching
#' `careless::psychsyn(resample_na = FALSE)` to the bit.
#'
#' **No pairs found.** `critical_r` is a property of the *inventory*, not the
#' respondent: on broad personality inventories the strongest inter-item
#' correlations can fall below `0.60`, so no synonym pairs qualify, every value is
#' `NA`, and the index abstains with a warning. Lower `critical_r` to surface
#' pairs, or compute the index on a more homogeneous scale.
#'
#' **Cutoff.** The default flags the lowest-consistency respondents: the cutoff is
#' the empirical `fpr` quantile of the observed scores (the 5th percentile by
#' default) and respondents at or below it are flagged. This is a **ranking**
#' convention -- a sample percentile flags **at least** `fpr` of respondents (more
#' when scores tie at the cutoff) -- not a calibrated false-positive rate. Adjust
#' the target tail with `fpr`, or
#' pass an absolute `cutoff` in `[-1, 1]` to flag on a literal correlation
#' threshold. `fpr` and `cutoff` are mutually exclusive.
#'
#' **Abstention.** A respondent with fewer than three complete synonym pairs (an
#' all-missing row included), or whose answered pair values have zero variance (a
#' straightliner), has no defined correlation: both `value` and `flagged` are `NA`
#' and the row is excluded from the flag count and rate.
#'
#' @section What this catches:
#' Random and partial / changepoint responding, whose answers to near-synonymous
#' items lose the sample's consistency. It **misses** pure straightlining and
#' midpoint locking (a constant row has zero variance and abstains); pair it with
#' [cier_longstring()] and [cier_irv()].
#'
#' @param responses A numeric matrix (or a data.frame / tibble coerced
#'   internally) of responses, one row per respondent and one column per item.
#'   `NA` marks a missing response.
#' @param critical_r The whole-sample inter-item correlation a pair must exceed to
#'   count as a synonym. A single number in the open interval `(0, 1)`; the
#'   default `0.60` follows Meade & Craig (2012).
#' @param fpr Optional target false-positive tail mass for the percentile cutoff.
#'   `NULL` (default) uses the registry default `0.05`. A finite number in the
#'   open interval `(0, 1)`; the cutoff is that lower-tail quantile of the
#'   observed scores. Mutually exclusive with `cutoff`.
#' @param cutoff Optional **literal** cutoff on the correlation, supplied instead
#'   of `fpr`. A single finite number in `[-1, 1]`; respondents whose synonym
#'   correlation is at or below it are flagged (the lower-tail direction).
#'
#' @return A `cier_index`: a list with per-respondent `value` (numeric, `NA` on
#'   abstention) and `flagged` (logical) vectors plus the `method`, `cutoff`, and
#'   `direction` metadata. Use [as.data.frame()][as.data.frame.cier_index] for a
#'   tidy data frame and `print()` for a summary.
#'
#' @references
#' Meade, A. W., & Craig, S. B. (2012). Identifying careless responses in survey
#' data. *Psychological Methods*, 17(3), 437–455. \doi{10.1037/a0028085}
#'
#' @seealso [careless::psychsyn()], [cier_person_total()], [cier_irv()]
#' @family indirect indices
#' @export
#' @examples
#' # The 44 BFI items are the first 44 columns of the bundled example data. The
#' # BFI's strongest inter-item correlations fall below the default critical_r of
#' # 0.60, so lower it to surface synonym pairs on this inventory.
#' out <- cier_psychsyn(bfi_careless[, 1:44], critical_r = 0.5)
#' out
#' head(as.data.frame(out))
cier_psychsyn <- function(responses, critical_r = 0.60, fpr = NULL,
                          cutoff = NULL) {
  call <- rlang::caller_env()
  responses <- check_responses(responses, call = call)
  check_open_unit(critical_r, "critical_r", call = call)
  # `fpr` is a tail mass in (0, 1); a literal `cutoff` is a correlation threshold
  # in [-1, 1] (outside it flags everyone or no one); mutually exclusive.
  check_percentile_overrides(fpr, cutoff, lower = -1, upper = 1, call = call)
  row <- cier_method_row("cier_psychsyn")
  # One pairing correlation serves discovery and scoring. The kernel returns NA
  # where a respondent has fewer than three complete synonym pairs; when NO pair
  # qualifies at all, the shared tail swaps the generic percentile abstention
  # for the actionable no-pairs warning (cause + remedy).
  cor_mat <- pairing_cor(responses)
  pairs <- find_item_pairs(responses, critical_r, "syn", cor_mat = cor_mat)
  value <- kernel_psychsyn(responses, critical_r, "syn", cor_mat = cor_mat)
  resolve_pair_index_cutoff(value, row, fpr, cutoff, nrow(pairs) == 0L, "syn",
                            critical_r, cor_mat, call)
}
