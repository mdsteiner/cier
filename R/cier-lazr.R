# Purpose: cier_lazr() -- the public first-order Markov "Lazy Respondents"
#          (Laz.R) predictability index (Biemann, Koch-Bayram, Meier-Barthold &
#          Aguinis 2025, Organizational Research Methods).
# Args:    See documentation below.
# Returns: A light `cier_index` (see new_cier_index()).
# Invariants:
#   - Reproduces Eq. 3 of Biemann et al. (2025): John worked example = 33/49,
#     the one-liner Laz.R(c(1,2,3,4,5,4,3,2,1,2)) = 2/3, straightliner = 1.
#   - Oracle-only trust (no CRAN partner); see tests/reference/TOLERANCES.md.
#   - The cutoff routes through the single resolve_index_cutoff() path
#     (percentile method, upper direction, fpr = 0.05 by default).

#' Laz.R first-order Markov predictability C/IER index
#'
#' Computes each respondent's **Laz.R** score (Biemann, Koch-Bayram,
#' Meier-Barthold & Aguinis, 2025): the average probability with which their
#' previous answer predicts the next. Over each respondent's response sequence
#' the function builds the transition-count matrix \eqn{T} and its row-normalised
#' probability matrix \eqn{P}, and returns
#' \deqn{\mathrm{Laz.R} = \frac{\sum_{i,j} P_{ij}\,T_{ij}}{n_{\mathrm{trans}}}}
#' the count-weighted average transition probability, bounded by roughly
#' \eqn{1/s} and \eqn{1} for an \eqn{s}-anchor scale. A respondent who repeats a
#' predictable pattern (straightlining, midpoint or extreme locking,
#' diagonal-lining `1-2-3-4-5-1-2-...`, or a seesaw `1-5-1-5-...`) has a high
#' Laz.R, so **high** values flag carelessness. The score is computed on the
#' **raw** response row in administration order -- there is no reverse-keying and
#' no scale blocking.
#'
#' @details
#' **Cutoff.** The default flags the most predictable respondents: the cutoff is
#' the empirical upper `fpr` quantile of the observed scores (the 95th percentile
#' by default). This is a deliberate divergence from Biemann et al.'s
#' sample-specific Kneedle elbow -- it follows cier's single-`fpr`-knob ranking
#' convention, and the Laz.R score's documented dependence on sequence length
#' makes an absolute cutoff indefensible. Adjust the target tail with `fpr`, or
#' pass an absolute `cutoff` in `[0, 1]` to flag on a literal Laz.R magnitude;
#' `fpr` and `cutoff` are mutually exclusive.
#'
#' **Integer anchors.** Laz.R scores a first-order Markov chain over discrete
#' response anchors, so responses must be integer-coded; non-integer (averaged or
#' POMP-rescaled) input is a typed error rather than being silently dropped. The
#' value is invariant to the anchor count and to the scale base: a 0-based or
#' bipolar coding, and any unused higher anchor, score identically.
#'
#' **Missing data.** A transition is counted only when both of its endpoints are
#' present; a missing response drops the two transitions it touches, and the
#' denominator is the count of valid transitions (not the item count minus one).
#' A respondent with fewer than **two** valid transitions abstains -- both
#' `value` and `flagged` are `NA` and the row is excluded from the flag count and
#' rate -- because a single transition gives the uninformative `1` regardless of
#' the pattern. This drops missing transitions rather than treating a gap as a
#' predictable state, so a careful but incomplete respondent is not flagged.
#' Numeric missing codes (`99`, `-9`, and the like) must be recoded to `NA`
#' first; otherwise they are scored as genuine response anchors and add spurious
#' transitions.
#'
#' **Assumptions.** The index reads the response columns as an ordered sequence
#' on a common response scale, so it is meaningful only when the columns are in
#' administration order and share the same response options. Biemann et al.
#' recommend at least 20 items and a uniform anchor count across items; on a
#' mixed-format survey, score each same-format block in a separate call.
#'
#' @section What this catches:
#' Predictable / repetitive responding -- straightlining, midpoint and extreme
#' locking, diagonal-lining, and seesaw / alternating patterns -- whose next
#' answer is well predicted by the previous. It is **not** useful against uniform
#' random responding, whose Laz.R is low by design; pair it with [cier_irv()],
#' [cier_mahalanobis()], or [cier_psychsyn()], which catch random responding by
#' other means.
#'
#' @param responses A numeric matrix (or a data.frame / tibble coerced
#'   internally) of integer-coded responses, one row per respondent and one
#'   column per item, in administration order. `NA` marks a missing response.
#' @param fpr Optional target false-positive tail mass for the percentile cutoff.
#'   `NULL` (default) uses the registry default `0.05`. A finite number in the
#'   open interval `(0, 1)`; the cutoff is that upper-tail quantile of the
#'   observed scores. Mutually exclusive with `cutoff`.
#' @param cutoff Optional **literal** cutoff on the Laz.R score, supplied instead
#'   of `fpr`. A single finite number in `[0, 1]`; respondents whose score is at
#'   or above it are flagged.
#'
#' @return A `cier_index`: a list with per-respondent `value` (numeric, `NA` on
#'   abstention) and `flagged` (logical) vectors plus the `method`, `cutoff`, and
#'   `direction` metadata. Use [as.data.frame()][as.data.frame.cier_index] for a
#'   tidy data frame and `print()` for a summary.
#'
#' @references
#' Biemann, T., Koch-Bayram, I. F., Meier-Barthold, M., & Aguinis, H. (2025).
#' Using Markov chains to detect careless responding in survey research.
#' *Organizational Research Methods*, 28(4), 543-568.
#' \doi{10.1177/10944281251334778}
#'
#' @seealso [cier_autocorrelation()], [cier_longstring()]
#' @family indirect indices
#' @export
#' @examples
#' # The 44 BFI items are the first 44 columns of the bundled example data.
#' out <- cier_lazr(bfi_careless[, 1:44])
#' out
#' head(as.data.frame(out))
cier_lazr <- function(responses, fpr = NULL, cutoff = NULL) {
  call <- rlang::caller_env()
  responses <- check_responses(responses, call = call)
  check_integer_responses(responses, call = call)
  # `value` is a Laz.R magnitude in (~1/s, 1]; `fpr` is a tail mass in (0, 1), a
  # literal `cutoff` a magnitude in [0, 1]; they are mutually exclusive.
  check_percentile_overrides(fpr, cutoff, lower = 0, upper = 1, call = call)
  row <- cier_method_row("cier_lazr")
  value <- kernel_lazr(responses)
  resolve_index_cutoff(value, row, fpr, cutoff, call = call)
}
