# Purpose: cier_psychant() -- the public psychometric-antonyms C/IER index.
# Args:    See documentation below.
# Returns: A light `cier_index` (see new_cier_index()).
# Invariants:
#   - Deterministic (no RNG); whole-sample pairing on RAW responses, with no
#     reverse-keying (reverse-scoring would collapse the genuine antonym
#     correlations toward zero and so must NOT precede pairing).
#   - Shares the single psychsyn/antonym kernel with the "ant" tail: scores cor()
#     on the full stacked pair vectors, matching careless's antonym scorer
#     (psychsyn(anto = TRUE, resample_na = FALSE)) to 1e-12.
#   - The cutoff routes through the single resolve_cutoff() path (percentile
#     method, upper direction, fpr = 0.05 by default).

#' Psychometric antonyms (within-person antonym opposition) C/IER index
#'
#' Computes each respondent's **psychometric-antonyms** score -- the Pearson
#' correlation, within that respondent, across the item pairs that are strongly
#' *negatively* correlated in the whole sample (the "antonyms"). A careful
#' respondent answers genuine antonyms (for instance a trait item and its
#' reverse-keyed twin) in opposition, so their within-person pattern mirrors the
#' sample's and the correlation stays strongly negative; a careless respondent's
#' answers to antonyms decouple, driving the correlation up toward zero. **High**
#' values therefore flag carelessness -- the mirror of [cier_psychsyn()], whose
#' low values flag (Meade & Craig, 2012).
#'
#' @details
#' **Pairing on raw responses.** The antonym pairs are the item pairs whose
#' whole-sample inter-item correlation is below `-critical_r` (more negative than
#' `-0.60` by default). Pairing uses the **raw** responses, with no reverse-
#' scoring: reverse-keying would collapse the very negative correlations that
#' define an antonym pair toward zero. The score is the single correlation across
#' all qualifying pairs' stacked values, matching careless's antonym scorer
#' (`careless::psychsyn(anto = TRUE, resample_na = FALSE)`) to 1e-12.
#'
#' **`critical_r` is a positive magnitude.** Supply `critical_r` as a **positive**
#' number in `(0, 1)` (default `0.60`); the index applies it to the negative tail
#' internally, keeping pairs whose inter-item correlation is below `-critical_r`.
#' This matches
#' [cier_psychsyn()] and the pairing diagnostics
#' ([cier_synonym_pairs()], [cier_psychsyn_critval()] with `antonym = TRUE`), so a
#' single positive threshold means the same pair strength for synonyms and
#' antonyms alike. Note that `careless::psychant()` instead takes a **negative**
#' `critval` (e.g. `-0.60`); here pass the positive magnitude `0.60`. A negative
#' `critical_r` is rejected as an input error.
#'
#' **No pairs found.** `critical_r` is a property of the *inventory*, not the
#' respondent: a survey may contain no two items that oppose each other strongly
#' enough, so no antonym pairs qualify, every value is `NA`, and the index
#' abstains with a warning. Lower `critical_r` to surface pairs, or compute the
#' index on a scale that pairs forward and reverse items. Use
#' [cier_synonym_pairs()] (`antonym = TRUE`) to see which pairs qualify.
#'
#' **Cutoff.** The default flags the lowest-opposition respondents: the cutoff is
#' the empirical `1 - fpr` quantile of the observed scores (the 95th percentile by
#' default) and respondents at or above it are flagged. This is a **ranking**
#' convention -- a sample percentile flags `fpr` of respondents by construction --
#' not a calibrated false-positive rate. Adjust the target tail with `fpr`, or
#' pass an absolute `cutoff` in `[-1, 1]` to flag on a literal correlation
#' threshold. `fpr` and `cutoff` are mutually exclusive.
#'
#' **Abstention.** A respondent with fewer than three complete antonym pairs (an
#' all-missing row included), or whose answered pair values have zero variance (a
#' straightliner), has no defined correlation: both `value` and `flagged` are `NA`
#' and the row is excluded from the flag count and rate.
#'
#' @section What this catches:
#' Random and partial / changepoint responding, whose answers to opposed items
#' lose the sample's consistency. It **misses** pure straightlining and midpoint
#' locking (a constant row has zero variance and abstains); pair it with
#' [cier_longstring()] and [cier_irv()]. It also needs genuine antonym pairs to
#' exist in the inventory -- forward/reverse item pairs are the usual source.
#'
#' @param responses A numeric matrix (or a data.frame / tibble coerced
#'   internally) of responses, one row per respondent and one column per item.
#'   `NA` marks a missing response.
#' @param critical_r The whole-sample inter-item correlation magnitude a pair must
#'   exceed (on the negative tail) to count as an antonym. A single **positive**
#'   number in the open interval `(0, 1)`; the default `0.60` follows Meade &
#'   Craig (2012) and keeps pairs whose inter-item correlation is below `-0.60`.
#'   (Unlike `careless::psychant()`, which takes a negative `critval`.)
#' @param fpr Optional target false-positive tail mass for the percentile cutoff.
#'   `NULL` (default) uses the registry default `0.05`. A finite number in the
#'   open interval `(0, 1)`; the cutoff is that upper-tail quantile of the
#'   observed scores. Mutually exclusive with `cutoff`.
#' @param cutoff Optional **literal** cutoff on the correlation, supplied instead
#'   of `fpr`. A single finite number in `[-1, 1]`; respondents whose antonym
#'   correlation is at or above it are flagged (the upper-tail direction).
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
#' @seealso [careless::psychant()], [cier_psychsyn()], [cier_synonym_pairs()],
#'   [cier_irv()]
#' @family indirect indices
#' @export
#' @examples
#' # The 44 BFI items are the first 44 columns of the bundled example data. A
#' # broad personality inventory's strongest item oppositions are milder than the
#' # default critical_r of 0.60 (the BFI's most negative inter-item r is about
#' # -0.48), so lower it to surface antonym pairs -- typically forward/reverse
#' # item pairs within a scale.
#' out <- cier_psychant(bfi_careless[, 1:44], critical_r = 0.40)
#' out
#' head(as.data.frame(out))
cier_psychant <- function(responses, critical_r = 0.60, fpr = NULL,
                          cutoff = NULL) {
  call <- rlang::caller_env()
  responses <- check_responses(responses, call = call)
  check_open_unit(critical_r, "critical_r", call = call)
  # `fpr` is a tail mass in (0, 1); a literal `cutoff` is a correlation threshold
  # in [-1, 1] (outside it flags everyone or no one); mutually exclusive.
  check_percentile_overrides(fpr, cutoff, lower = -1, upper = 1, call = call)
  row <- cier_method_row("cier_psychant")
  # The kernel returns NA where a respondent has fewer than three complete antonym
  # pairs (or no pair qualifies at all), so abstention needs no separate guard --
  # the shared percentile cutoff abstains when every respondent is NA. The "ant"
  # tail keeps pairs with r < -critical_r (psychsyn passes "syn" for r > critical_r).
  value <- kernel_psychsyn(responses, critical_r, "ant")
  resolve_index_cutoff(value, row, fpr, cutoff, call = call)
}
