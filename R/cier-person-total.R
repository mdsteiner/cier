# Purpose: cier_person_total() -- the public person-total correlation (r_pbis)
#          index.
# Args:    See documentation below.
# Returns: A light `cier_index` (see new_cier_index()).
# Invariants:
#   - Deterministic (no RNG); whole-sample, raw responses, keying-insensitive.
#   - Item-total form (per-item mean over all respondents, as in PerFit::r.pbis);
#     equal to PerFit::r.pbis() within 1e-4 (its 4-dp output rounding).
#   - The cutoff routes through the single resolve_cutoff() path (percentile
#     method, lower direction, fpr = 0.05 by default).

#' Person-total correlation (r_pbis) C/IER index
#'
#' Computes each respondent's **person-total correlation** -- the Pearson
#' correlation between their item responses and the per-item means computed
#' across the whole sample (the "group-determined item difficulties" of Donlon &
#' Fischer, 1968). A careful respondent endorses items roughly in line with the
#' group, so their pattern tracks the item-mean profile and the correlation is
#' high; a careless respondent's pattern decouples from the group, driving the
#' correlation toward zero. **Low** values therefore flag carelessness. The score
#' is computed over the raw response row (no scale blocking), so -- unlike the
#' split-half consistency indices -- it also works on a single scale.
#'
#' @details
#' **Item-total, raw, keying-insensitive.** This is the item-**total** form (the
#' per-item mean is taken over *all* respondents, including the respondent being
#' scored, as in `PerFit::r.pbis`); it is numerically near-identical to the
#' item-rest (leave-one-out) variant on real data and vectorises far more
#' cheaply. Because the per-item means are empirical, the index is
#' **keying-insensitive**: it is computed on the raw responses with no
#' reverse-scoring, so it takes no item metadata and a reverse-keyed item does
#' not need to be declared.
#'
#' **Cutoff.** The default flags the lowest-correlation respondents: the cutoff is
#' the empirical `fpr` quantile of the observed scores (the 5th percentile by
#' default) and respondents at or below it are flagged. This is a **ranking**
#' convention -- a sample percentile flags **at least** `fpr` of respondents (more
#' when scores tie at the cutoff) -- not a calibrated false-positive rate. Adjust
#' the target tail with `fpr`, or
#' pass an absolute `cutoff` in `[-1, 1]` to flag on a literal correlation
#' threshold (e.g. one carried over from a calibration sample). `fpr` and `cutoff`
#' are mutually exclusive.
#'
#' **Abstention.** A respondent who answered fewer than three items, or whose
#' answered responses have zero variance (a straightliner), has no defined
#' correlation: both `value` and `flagged` are `NA` and the row is excluded from
#' the flag count and rate. An item answered by nobody drops out rather than
#' poisoning the computation.
#'
#' @section What this catches:
#' Random and partial / changepoint responding and idiosyncratic patterns
#' (diagonal, alternating) whose profile diverges from the sample's item-mean
#' profile. It **misses** pure straightlining and midpoint locking (a constant
#' row has zero variance and abstains); pair it with [cier_longstring()] and
#' [cier_irv()].
#'
#' @param responses A numeric matrix (or a data.frame / tibble coerced
#'   internally) of responses, one row per respondent and one column per item.
#'   `NA` marks a missing response.
#' @param fpr Optional target false-positive tail mass for the percentile
#'   cutoff. `NULL` (default) uses the registry default `0.05`. A finite number
#'   in the open interval `(0, 1)`; the cutoff is that lower-tail quantile of the
#'   observed scores. Mutually exclusive with `cutoff`.
#' @param cutoff Optional **literal** cutoff on the correlation, supplied instead
#'   of `fpr`. A single finite number in `[-1, 1]`; respondents whose
#'   person-total correlation is at or below it are flagged (the lower-tail
#'   direction). Use it to apply an absolute threshold rather than a sample
#'   percentile.
#'
#' @return A `cier_index`: a list with per-respondent `value` (numeric, `NA` on
#'   abstention) and `flagged` (logical) vectors plus the `method`, `cutoff`, and
#'   `direction` metadata. Use [as.data.frame()][as.data.frame.cier_index] for a
#'   tidy data frame and `print()` for a summary.
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
#' @seealso [PerFit::r.pbis()], [cier_mahalanobis()], [cier_irv()]
#' @family indirect indices
#' @export
#' @examples
#' # The 44 BFI items are the first 44 columns of the bundled example data.
#' out <- cier_person_total(bfi_careless[, 1:44])
#' out
#' head(as.data.frame(out))
cier_person_total <- function(responses, fpr = NULL, cutoff = NULL) {
  call <- rlang::caller_env()
  responses <- check_responses(responses, call = call)
  # `fpr` is a tail mass in (0, 1); a literal `cutoff` is a correlation threshold
  # in [-1, 1] (outside it flags everyone or no one); mutually exclusive.
  check_percentile_overrides(fpr, cutoff, lower = -1, upper = 1, call = call)
  row <- cier_method_row("cier_person_total")
  # The kernel returns NA where a respondent answered fewer than three items or
  # either side of the correlation has zero variance, so abstention needs no
  # separate guard (cf. irv).
  value <- kernel_person_total(responses)
  resolve_index_cutoff(value, row, fpr, cutoff, call = call)
}
