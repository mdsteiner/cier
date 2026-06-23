#' Normed polytomous Guttman-error (Gnormed) C/IER index
#'
#' Computes each respondent's normed polytomous Guttman-error score (Emons, 2008;
#' Molenaar, 1991). Ordering the item-steps by descending popularity, a Guttman error
#' endorses a less popular step while missing a more popular one; the count is normed
#' by the maximum possible given the respondent's total step score. Aberrant patterns
#' score high -- high values flag carelessness. Nonparametric (no item response model
#' is fitted).
#'
#' @details
#' **Computation.** Item-step popularity positions normed by an n-independent max-plus
#' knapsack maximum. Reverse-keyed items (`items$reverse_keyed`) are reverse-scored
#' with `(min + max) - x` (`min` defaults to `1`), then recoded to the `0..(Ncat - 1)`
#' contract; supply raw responses and declare reverse items through `items`. Scores a
#' single category count (`Ncat = max - min + 1`), so every item must share the same
#' `max - min` span (the base may differ: a `1..5` and a `0..4` item score together).
#' Mixed-format data raises a typed error with a dedicated subclass, so [cier_screen()]
#' skips Gnormed rather than aborting the battery -- score homogeneous subsets
#' separately.
#'
#' **Cutoff -- a Monte-Carlo null.** The default references a sum-score-conditional
#' nonparametric resample (resample observed sum scores, then draw each item from its
#' category frequencies among respondents at that score), scored by the same closed
#' form and summarised by the upper-tail bootstrap quantile at level `fpr`. Being a
#' real null, the flag rate is not pinned at `fpr`. With `seed = NULL` the cutoff and
#' flags can vary run to run; pass an integer `seed` to fix them (applied locally).
#' Override with one of `fpr` or `cutoff`.
#'
#' **The sample must span the declared scale.** Item-step popularities are undefined
#' when a scale end never occurs, so the complete-case block must attain both the
#' declared `min` and `max`. Otherwise a typed error is raised (same subclass as the
#' mixed-format case, so [cier_screen()] skips rather than aborts).
#'
#' **Abstention.** Computed on complete cases, so any missing cell excludes a
#' respondent (`value`, `flagged` are `NA`). With fewer than two complete respondents
#' or fewer than three items, every `value` is `NA`, a warning is emitted, and no one
#' is flagged.
#'
#' @section What this catches:
#' Aberrant, inconsistent patterns -- random and some patterned responding -- with
#' many Guttman reversals. Blind spot: straightlining scores a low, non-aberrant value
#' rather than abstaining, so it evades this index -- pair with [cier_longstring()] and
#' [cier_irv()].
#'
#' @template responses
#' @param items Item-metadata data.frame, one row per item aligned to the columns of
#'   `responses`. Needs an integer `max` column (largest response option) on every
#'   item, all sharing the span `max - min` (`Ncat = max - min + 1 >= 2`). Optional
#'   logical `reverse_keyed` (default none) and integer `min` (scale base, default
#'   `1`). When `responses` has column names, an optional `item` column (or matching
#'   row names) is cross-checked, so a reordered frame is a typed error.
#' @param fpr Optional nominal level for the Monte-Carlo null cutoff, a number in
#'   `(0, 1)`; `NULL` (default) uses `0.05`. Mutually exclusive with `cutoff`.
#' @param cutoff Optional literal cutoff on the score, a number in `[0, 1]`;
#'   respondents at or above it are flagged.
#' @param seed Optional integer making the Monte-Carlo null reproducible; applied
#'   locally and restores the session's random state. `NULL` (default) draws from the
#'   ambient random stream, so cutoff and flags can vary. Ignored when `cutoff` is
#'   supplied.
#'
#' @template return-cier-index
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
#' @seealso [cier_longstring()], [cier_irv()]
#' @family person-fit indices
#' @export
#' @examples
#' # Columns 1:44 are the BFI items (5-point, 1..5); a trailing "_R" marks reverse.
#' nm <- names(bfi_careless)[1:44]
#' items <- data.frame(reverse_keyed = grepl("_R$", nm), max = 5L)
#' out <- cier_gnormed(bfi_careless[, 1:44], items, seed = 1)
#' out
#' head(as.data.frame(out))
cier_gnormed <- function(responses, items, fpr = NULL, cutoff = NULL,
                         seed = NULL) {
  call <- rlang::caller_env()
  responses <- check_responses(responses, call = call)
  inform_if_unkeyed(items, ncol(responses))
  items <- check_items_personfit(items, ncol(responses), call = call,
                                 response_names = colnames(responses))
  check_percentile_overrides(fpr, cutoff, lower = 0, upper = 1, call = call)
  if (!is.null(seed)) check_int(seed, "seed", call = call)
  row <- cier_method_spec("cier_gnormed")
  # Reverse-score keyed items, then score the complete-case block. The kernel returns
  # per-respondent values (NA on abstention) plus the zero-based block the
  # Monte-Carlo null regenerates from.
  responses <- apply_split_half_keying(responses, items, call = call)
  res <- kernel_gnormed(responses, items$ncat, items$min, call = call)
  cutoff_value <- resolve_gnormed_cutoff(res, fpr, cutoff,
                                         row$default_cutoff_value, seed,
                                         call = call)
  flag_and_assemble(res$value, cutoff_value, row, cutoff, "mc_null",
                    if (is.null(fpr)) row$default_cutoff_value else fpr,
                    call = call)
}
