# Purpose: cier_ht() -- the public polytomous person-scalability (Ht) C/IER index,
#          the second external-backend bridge (mokken).
# Args:    See documentation below.
# Returns: A light `cier_index` (see new_cier_index()).
# Invariants:
#   - The production scorer is mokken::coefH on the transposed scale
#     (single-kernel rule); the wrapper reverse-keys, then delegates to kernel_ht.
#   - mokken is required (a typed cier_error_input when absent): the standalone
#     index cannot compute its statistic without it.
#   - Deterministic: the cutoff routes through the single resolve_cutoff() path
#     (percentile method, lower direction, fpr = 0.05 by default). Unlike Gnormed
#     there is no Monte-Carlo null -- none exists for the mokken-backed polytomous
#     Ht -- so there is no `seed` argument.

#' Polytomous person-scalability (Ht) C/IER index
#'
#' Computes each respondent's **person-scalability** coefficient Ht (Sijtsma &
#' Meijer, 1992; Molenaar, 1991) -- Loevinger's scalability coefficient H computed
#' in the transposed (person) space via [mokken::coefH()] on the transposed
#' response matrix. A careful respondent answers in a way that scales consistently
#' with the sample's item ordering, giving a high Ht; an aberrant (careless)
#' response pattern scales poorly, so **low** values flag carelessness. The
#' statistic is **nonparametric** (no item response model is fitted).
#'
#' @details
#' **mokken is required.** Ht is scored by [mokken::coefH()], so the `mokken`
#' package must be installed; the index aborts with an informative error
#' otherwise. Reverse-keyed items (`items$reverse_keyed`) are reverse-scored first
#' with the self-inverse reflection `(min + max) - x` (where
#' `max = min + categories - 1`, and `min` defaults to `1`, i.e. the classic
#' `(categories + 1) - x`); scalability is then computed on the transposed scale.
#' Supply the raw responses and declare reverse items through `items`. Only
#' reverse-keyed items need a `categories` value (used to reverse-score them);
#' `mokken` itself accepts a mix of category counts, so no common number of
#' response categories is required.
#'
#' **Cutoff -- an empirical percentile (a ranking convention).** Ht's default
#' cutoff is the `fpr` quantile of the observed scores (the 5th percentile by
#' default); respondents at or below it are flagged. This is a **ranking**
#' convention -- a sample percentile flags `fpr` of respondents by construction --
#' not a calibrated false-positive rate. Unlike [cier_gnormed()], Ht has **no
#' Monte-Carlo null**: a model-conforming simulated null is available only for
#' `PerFit`'s own dichotomous Ht, not for the polytomous `mokken` statistic
#' computed here. Adjust the target tail with `fpr`, or pass an absolute `cutoff`
#' in `[-1, 1]` to flag on a literal scalability threshold (e.g. one carried over
#' from a calibration sample). `fpr` and `cutoff` are mutually exclusive.
#'
#' **Abstention.** Because [mokken::coefH()] cannot take missing values, a
#' respondent with **any** missing cell is excluded: both `value` and `flagged`
#' are `NA` and the row is dropped from the flag count and rate. A **straightliner**
#' (a fully-answered but zero-variance row) is structurally unscalable, so it too
#' abstains (`NA`) -- the complement of [cier_gnormed()], which *scores*
#' straightliners. If fewer than two respondents are complete, or fewer than two of
#' the complete respondents vary, every `value` is `NA`, a warning is emitted, and
#' no one is flagged.
#'
#' @section What this catches:
#' Aberrant, inconsistent response patterns -- random and patterned responding --
#' that scale poorly against the sample's item ordering. Unlike most indirect
#' indices it **abstains on** rather than misses pure straightlining (a
#' zero-variance row is unscalable); pair it with [cier_longstring()] and
#' [cier_irv()], which target the straightlining it cannot score.
#'
#' @param responses A numeric matrix (or a data.frame / tibble coerced
#'   internally) of responses, one row per respondent and one column per item.
#'   `NA` marks a missing response.
#' @param items A data.frame of item metadata, one row per item, aligned to the
#'   columns of `responses`. An optional logical `reverse_keyed` column marks
#'   reverse-keyed items (default: none); each reverse-keyed item additionally
#'   needs an integer `categories` value (the number of response options, `>= 2`)
#'   and an optional integer `min` (the scale base; default `1`, declare it for a
#'   0-based or bipolar scale) so it can be reverse-scored. Forward-keyed items
#'   need no metadata, and `categories` may be heterogeneous across items
#'   (`mokken` does not require a common number of response categories).
#' @param fpr Optional target false-positive tail mass for the percentile cutoff.
#'   `NULL` (default) uses the registry default `0.05`. A finite number in the open
#'   interval `(0, 1)`; the cutoff is that lower-tail quantile of the observed
#'   scores. Mutually exclusive with `cutoff`.
#' @param cutoff Optional **literal** cutoff on the scalability coefficient,
#'   supplied instead of `fpr`. A single finite number in `[-1, 1]`; respondents
#'   whose Ht is at or below it are flagged (the lower-tail direction). Use it to
#'   apply an absolute threshold rather than a sample percentile.
#'
#' @return A `cier_index`: a list with per-respondent `value` (numeric, `NA` on
#'   abstention) and `flagged` (logical) vectors plus the `method`, `cutoff`, and
#'   `direction` metadata. Use [as.data.frame()][as.data.frame.cier_index] for a
#'   tidy data frame and `print()` for a summary.
#'
#' @references
#' Niessen, A. S. M., Meijer, R. R., & Tendeiro, J. N. (2016). Detecting careless
#' respondents in web-based questionnaires: Which method to use? *Journal of
#' Research in Personality*, 63, 1–11. \doi{10.1016/j.jrp.2016.04.010}
#'
#' Sijtsma, K., & Meijer, R. R. (1992). A method for investigating the
#' intersection of item response functions in Mokken's nonparametric IRT model.
#' *Applied Psychological Measurement*, 16(2), 149–157.
#'
#' Molenaar, I. W. (1991). A weighted Loevinger H-coefficient extending Mokken
#' scaling to multicategory items. *Kwantitatieve Methoden*, 12(37), 97–117.
#'
#' @seealso [mokken::coefH()], [cier_gnormed()], [cier_longstring()], [cier_irv()]
#' @family person-fit indices
#' @export
#' @examplesIf requireNamespace("mokken", quietly = TRUE)
#' # The 44 BFI items are the first 44 columns of the bundled example data; they
#' # are 5-point items. Build item metadata from the column names: a trailing "_R"
#' # marks a reverse-keyed item.
#' nm <- names(bfi_careless)[1:44]
#' items <- data.frame(reverse_keyed = grepl("_R$", nm), categories = 5L)
#' out <- cier_ht(bfi_careless[, 1:44], items)
#' out
#' head(as.data.frame(out))
cier_ht <- function(responses, items, fpr = NULL, cutoff = NULL) {
  call <- rlang::caller_env()
  responses <- check_responses(responses, call = call)
  items <- check_items_ht(items, ncol(responses), call = call)
  # `fpr` is a tail mass in (0, 1); a literal `cutoff` is a scalability value in
  # [-1, 1] (outside it flags everyone or no one); the two are mutually exclusive.
  check_percentile_overrides(fpr, cutoff, lower = -1, upper = 1, call = call)
  require_suggested("mokken", "cier_ht", call = call)
  row <- cier_method_row("cier_ht")
  # Reverse-score keyed items, then score on the transposed complete-case block.
  responses <- apply_split_half_keying(responses, items, call = call)
  value <- kernel_ht(responses, call = call)
  resolve_index_cutoff(value, row, fpr, cutoff, call = call)
}
