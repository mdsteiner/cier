# Purpose: cier_gnormed() -- the public normed polytomous Guttman-error (Gnormed)
#          C/IER index, the first external-backend bridge (PerFit).
# Args:    See documentation below.
# Returns: A light `cier_index` (see new_cier_index()).
# Invariants:
#   - The production scorer is PerFit::Gnormed.poly (single-kernel rule); the
#     wrapper reverse-keys, complete-cases, and zero-bases, then delegates.
#   - PerFit is required (a typed cier_error_input when absent): the standalone
#     index cannot compute its statistic without it.
#   - The default cutoff is the PerFit Monte-Carlo null (perfit_null), resolved at
#     the bridge from the fitted object; it is reproducible via `seed`.

#' Normed polytomous Guttman-error (Gnormed) C/IER index
#'
#' Computes each respondent's **normed polytomous Guttman-error** score (Emons,
#' 2008; Molenaar, 1991) via [PerFit::Gnormed.poly()]. Ordering the item-steps by
#' descending popularity, a Guttman error is an endorsement of a less popular step
#' while a more popular one is missed; the count is normed by the maximum possible
#' given the respondent's total step score. A careful respondent produces few
#' such reversals, so a careless (aberrant) response pattern scores **high** --
#' high values flag carelessness. The statistic is **nonparametric** (no item
#' response model is fitted).
#'
#' @details
#' **PerFit is required.** Gnormed is scored by [PerFit::Gnormed.poly()], so the
#' `PerFit` package must be installed; the index aborts with an informative error
#' otherwise. Reverse-keyed items (`items$reverse_keyed`) are reverse-scored first
#' with the self-inverse reflection `(min + max) - x` (`min` defaults to `1`,
#' i.e. the classic `(max + 1) - x`); the responses are then recoded to PerFit's
#' documented `0..(Ncat - 1)` contract. Supply the raw responses and declare
#' reverse items through `items`. PerFit scores a **single number of response
#' categories** (`Ncat = max - min + 1`), so every item must have the same
#' `max - min` span -- items may still differ in base (a `1..5` and a `0..4`
#' item both have five options and score together). Genuinely mixed-format
#' data (say four- and seven-option items) is a **backend limit**, not a
#' metadata mistake: the typed error carries a dedicated subclass, and
#' [cier_screen()] records Gnormed as skipped with that reason instead of
#' aborting the battery -- score homogeneous item subsets separately if you
#' need Gnormed on such data.
#'
#' **Cutoff -- the PerFit Monte-Carlo null.** Unlike the sample-percentile
#' indices, Gnormed's default cutoff is referenced to a **simulated null**: PerFit
#' resamples model-conforming response vectors and takes the nominal-rate
#' (`Blvl = fpr`) quantile of the statistic ([PerFit::cutoff()]). Respondents
#' beyond the upper tail are flagged. Because this is a real null rather than a
#' ranking convention, the flag rate is **not** pinned at `fpr`: it is whatever
#' share of respondents exceed the simulated quantile, and can be informatively
#' higher under contamination. The simulation is randomised; pass an integer
#' `seed` for a reproducible cutoff (it is applied locally and does not disturb
#' your session's random stream). Override the cutoff with **one** of two mutually
#' exclusive arguments: `fpr` (the nominal level) or `cutoff` (a literal threshold
#' on the score).
#'
#' **The sample must span the declared scale.** PerFit's item-step popularities
#' are undefined when a scale end never occurs, so the complete-case block must
#' **attain both extremes**: the declared `min` and `max` must each be chosen
#' at least once somewhere in the sample. A perfectly valid dataset in which,
#' say, nobody ever picked the top category raises a typed error rather than
#' scoring -- check `max` / `min` against the data, or expect this on small
#' samples with rarely-endorsed extreme categories. Like the mixed-format case
#' above, this is a **backend limit** on otherwise-valid data, not a metadata
#' mistake: the typed error carries the same dedicated subclass, so
#' [cier_screen()] records Gnormed as skipped with that reason instead of
#' aborting the battery.
#'
#' **Abstention.** Because [PerFit::Gnormed.poly()] needs complete data, a
#' respondent with **any** missing cell is excluded: both `value` and `flagged`
#' are `NA` and the row is dropped from the flag count and rate. If fewer than two
#' respondents are complete, or fewer than three items are present (PerFit's
#' recursive denominator needs `>= 3`), every `value` is `NA`, a warning is
#' emitted, and no one is flagged.
#'
#' @section What this catches:
#' Aberrant, inconsistent response patterns -- random and some patterned
#' responding -- that produce many Guttman reversals against the sample's
#' item-step popularities. It has a **documented blind spot**: a straightliner is
#' scored (it receives a *low*, non-aberrant value) rather than abstaining, so
#' straightlining evades this index -- pair it with [cier_longstring()] and
#' [cier_irv()].
#'
#' @param responses A numeric matrix (or a data.frame / tibble coerced
#'   internally) of responses, one row per respondent and one column per item.
#'   `NA` marks a missing response.
#' @param items A data.frame of item metadata, one row per item, aligned to the
#'   columns of `responses`. Must carry an integer `max` column (the largest
#'   response option) on every item, with the **same** span `max - min` on every
#'   item (PerFit's single `Ncat`; at least two options, `max >= min + 1`). An
#'   optional logical `reverse_keyed` column marks reverse-keyed items (default:
#'   none). An optional integer `min` column gives the smallest response option
#'   (the scale base; default `1`) -- declare it for a 0-based or bipolar scale.
#' @param fpr Optional nominal level for the Monte-Carlo null cutoff (PerFit's
#'   `Blvl`). `NULL` (default) uses the registry default `0.05`. A finite number
#'   in the open interval `(0, 1)`. Mutually exclusive with `cutoff`.
#' @param cutoff Optional **literal** cutoff on the score, supplied instead of
#'   `fpr`. A single finite number in `[0, 1]`; respondents whose Gnormed value is
#'   at or above it are flagged (the upper-tail direction).
#' @param seed Optional integer making the Monte-Carlo null cutoff reproducible.
#'   `NULL` (default) draws from the ambient random stream (the cutoff then varies
#'   per call); a non-`NULL` seed is applied locally and restores your session's
#'   random state on exit. Ignored when a literal `cutoff` is supplied.
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
#' Emons, W. H. M. (2008). Nonparametric person-fit analysis of polytomous item
#' scores. *Applied Psychological Measurement*, 32(3), 224–247.
#'
#' Molenaar, I. W. (1991). A weighted Loevinger H-coefficient extending Mokken
#' scaling to multicategory items. *Kwantitatieve Methoden*, 12(37), 97–117.
#'
#' @seealso [PerFit::Gnormed.poly()], [cier_longstring()], [cier_irv()]
#' @family person-fit indices
#' @export
#' @examplesIf requireNamespace("PerFit", quietly = TRUE)
#' # The 44 BFI items are the first 44 columns of the bundled example data; they
#' # are 5-point items coded 1..5. Build item metadata from the column names: a
#' # trailing "_R" marks a reverse-keyed item.
#' nm <- names(bfi_careless)[1:44]
#' items <- data.frame(reverse_keyed = grepl("_R$", nm), max = 5L)
#' out <- cier_gnormed(bfi_careless[, 1:44], items, seed = 1)
#' out
#' head(as.data.frame(out))
cier_gnormed <- function(responses, items, fpr = NULL, cutoff = NULL,
                         seed = NULL) {
  call <- rlang::caller_env()
  responses <- check_responses(responses, call = call)
  items <- check_items_personfit(items, ncol(responses), call = call)
  # `fpr` is the Monte-Carlo nominal level in (0, 1); a literal `cutoff` is a
  # normed-Guttman value in [0, 1]; the two are mutually exclusive. `seed` makes
  # the simulated null reproducible.
  check_percentile_overrides(fpr, cutoff, lower = 0, upper = 1, call = call)
  if (!is.null(seed)) check_int(seed, "seed", call = call)
  require_suggested("PerFit", "cier_gnormed", call = call)
  row <- cier_method_row("cier_gnormed")
  # Reverse-score keyed items, then score on the complete-case block. The kernel
  # returns the per-respondent values (NA where abstaining) and the fitted PerFit
  # object that the Monte-Carlo null cutoff reuses.
  responses <- apply_split_half_keying(responses, items, call = call)
  res <- kernel_gnormed(responses, items$ncat, items$min, call = call)
  cutoff_value <- resolve_gnormed_cutoff(res$fit, fpr, cutoff,
                                         row$default_cutoff_value, seed,
                                         call = call)
  flagged <- apply_flag(res$value, cutoff_value, row$flag_direction, call = call)
  new_cier_index(res$value, flagged, row$method, cutoff_value, row$flag_direction)
}
