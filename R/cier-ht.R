# cier_ht(): polytomous person-scalability (Ht) index. Reverse-keys, then delegates
# to kernel_ht.

#' Polytomous person-scalability (Ht) C/IER index
#'
#' Computes each respondent's person-scalability coefficient Ht (Sijtsma & Meijer,
#' 1992; Molenaar, 1991) -- Loevinger's scalability coefficient H in the transposed
#' (person) space. Aberrant patterns scale poorly against the sample's item ordering,
#' so low Ht flags carelessness. Nonparametric (no item response model is fitted).
#'
#' @details
#' **Computation.** Ht is a closed form -- the Frechet / rearrangement collapse that
#' reduces the transposed Mokken scalability to a per-respondent covariance ratio in
#' `O(n * p)` time and memory (one pass plus a per-respondent sort), without
#' materialising the n-by-n person matrix. Reverse-keyed items (`items$reverse_keyed`)
#' are reverse-scored first with `(min + max) - x` (`min` defaults to `1`). Only
#' reverse-keyed items need a `max`; heterogeneous category counts are accepted with
#' no ceiling, so wide scales (11-point, 0-100, mixed-width) score directly.
#'
#' **Cutoff.** By default the cutoff is the `fpr` quantile of the observed scores (5th
#' percentile); respondents at or below it are flagged. This ranking convention flags
#' at least `fpr` of respondents, not a calibrated false-positive rate. Setting
#' `method = "mc_null"` instead resolves the cutoff against a model-conforming
#' simulated null -- the same mechanism [cier_gnormed()] uses (a sum-score-conditional
#' resample of the scored block, excluding constant vectors), scored by the Ht kernel
#' and summarised by the lower-tail bootstrap quantile at level `fpr`; its flag rate
#' is not pinned at `fpr`. With `seed = NULL` the null cutoff and flags vary run to
#' run; pass an integer `seed` to fix them (applied locally). `fpr` and `cutoff` (an
#' absolute value in `[-1, 1]`) are mutually exclusive, and a literal `cutoff`
#' overrides `method` (then `seed` is unused).
#'
#' **Abstention.** Computed on complete cases: any missing cell excludes a respondent
#' (`value` and `flagged` are `NA`). A straightliner (fully-answered zero-variance row)
#' is structurally unscalable, so it too abstains -- the complement of [cier_gnormed()],
#' which scores straightliners. If fewer than two respondents are complete, or fewer
#' than two vary, every `value` is `NA`, a warning is emitted, and no one is flagged.
#'
#' @section What this catches:
#' Aberrant, inconsistent response patterns (random and patterned responding) that
#' scale poorly against the sample's item ordering. It abstains on rather than misses
#' pure straightlining; pair it with [cier_longstring()] and [cier_irv()].
#'
#' @template responses
#' @param items A data.frame of item metadata, one row per item, aligned to the
#'   columns of `responses`. Optional logical `reverse_keyed` marks reverse-keyed items
#'   (default none); each such item needs an integer `max` (largest option) and
#'   optional integer `min` (scale base, default `1`) for reverse-scoring. Forward-keyed
#'   items need no metadata, and `max` may be heterogeneous. When `responses` has column
#'   names, an optional `item` column (or matching row names) is cross-checked so a
#'   reordered frame is a typed error.
#' @param fpr Optional target tail mass -- the lower-tail quantile of the observed
#'   scores (percentile) or the nominal level of the `mc_null` null. `NULL` (default)
#'   uses `0.05`; otherwise a finite number in `(0, 1)`. Mutually exclusive with
#'   `cutoff`.
#' @param cutoff Optional literal cutoff on the scalability coefficient, a finite
#'   number in `[-1, 1]`; respondents at or below it are flagged. Overrides `method`.
#' @param method Cutoff strategy, `"percentile"` (default) or `"mc_null"` (see
#'   Details). Ignored when a literal `cutoff` is supplied.
#' @param seed Optional integer making the `mc_null` null reproducible. `NULL`
#'   (default) draws from the ambient random stream (null cutoff and flags vary call to
#'   call); a non-`NULL` seed is applied locally and restores the session's random
#'   state. No effect under the percentile method or a literal `cutoff`.
#'
#' @template return-cier-index
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
#' @seealso [cier_gnormed()], [cier_longstring()], [cier_irv()]
#' @family person-fit indices
#' @export
#' @examples
#' # First 44 columns are the BFI items; a trailing "_R" marks a reverse item.
#' nm <- names(bfi_careless)[1:44]
#' items <- data.frame(reverse_keyed = grepl("_R$", nm), max = 5L)
#' out <- cier_ht(bfi_careless[, 1:44], items)
#' out
#' head(as.data.frame(out))
#'
#' # Opt in to the Monte-Carlo null (reproducible via seed):
#' cier_ht(bfi_careless[, 1:44], items, method = "mc_null", seed = 1)
cier_ht <- function(responses, items, fpr = NULL, cutoff = NULL,
                    method = "percentile", seed = NULL) {
  call <- rlang::caller_env()
  responses <- check_responses(responses, call = call)
  inform_if_unkeyed(items, ncol(responses))
  items <- check_items_ht(items, ncol(responses), call = call,
                          response_names = colnames(responses))
  check_percentile_overrides(fpr, cutoff, lower = -1, upper = 1, call = call)
  check_choice(method, "method", c("percentile", "mc_null"), call = call)
  if (!is.null(seed)) check_int(seed, "seed", call = call)
  row <- cier_method_spec("cier_ht")
  # forward_range = FALSE: kernel zero-bases globally and tolerates heterogeneous
  # category counts, so only reverse items' declared range is range-checked.
  responses <- apply_split_half_keying(responses, items, call = call,
                                       forward_range = FALSE)
  res <- kernel_ht(responses, call = call)
  # Percentile resolves via resolve_index_cutoff (override = NULL); mc_null and a
  # literal `cutoff` are passed as explicit overrides, the latter taking precedence.
  override <- NULL
  provenance <- NULL
  if (!is.null(cutoff)) {
    override <- cutoff
  } else if (identical(method, "mc_null")) {
    fpr_level <- if (is.null(fpr)) row$default_cutoff_value else fpr
    override <- resolve_ht_cutoff(res, fpr_level, seed, call = call)
    provenance <- list(method = "mc_null", rate = fpr_level)
  }
  resolve_index_cutoff(res$value, row, fpr, override, provenance = provenance,
                       call = call)
}
