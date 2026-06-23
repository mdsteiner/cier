#' Laz.R first-order Markov predictability C/IER index
#'
#' Computes each respondent's Laz.R score (Biemann, Koch-Bayram, Meier-Barthold &
#' Aguinis, 2025): the count-weighted average probability with which the previous
#' answer predicts the next. Over each response sequence it builds the
#' transition-count matrix \eqn{T} and its row-normalised probability matrix \eqn{P}:
#' \deqn{\mathrm{Laz.R} = \frac{\sum_{i,j} P_{ij}\,T_{ij}}{n_{\mathrm{trans}}}}
#' A respondent who repeats a predictable pattern (straightlining, locking,
#' diagonal-lining, seesaw) scores high, so high values flag carelessness. Computed on
#' the raw response row in administration order -- no reverse-keying or scale blocking.
#'
#' @details
#' **Cutoff.** Respondents at or above the resolved cutoff are flagged. The default is
#' the upper `fpr` quantile of the scores (95th percentile, the Biemann et al.
#' convention); this flags at least `fpr` of respondents, not a calibrated
#' false-positive rate. Override mutually exclusively with `fpr`, a literal `cutoff`,
#' or `kneedle = TRUE` (the data-driven Kneedle elbow). On a strongly bimodal sample
#' the knee can flag a majority -- legitimate for a low-quality panel, but cier warns
#' when it flags more than half; for a fixed top-tail share prefer `fpr`.
#'
#' **Integer anchors.** Laz.R is a first-order Markov chain over discrete anchors, so
#' responses must be integer-coded (non-integer input errors). Invariant to anchor
#' count and scale base.
#'
#' **Missing data.** A transition counts only when both endpoints are present; a
#' missing response drops the two transitions it touches, and the denominator is the
#' valid-transition count. A respondent with fewer than two valid transitions abstains
#' (`value`, `flagged` are `NA`), since a lone transition gives the uninformative `1`.
#' Recode numeric missing codes (`99`, `-9`) to `NA` first, else they score as anchors.
#'
#' **Assumptions.** Columns are read as an ordered sequence on a common response scale,
#' so the index is meaningful only in administration order with the same options
#' throughout. Biemann et al. recommend at least 20 items; score each same-format block
#' of a mixed-format survey in a separate call.
#'
#' Standalone: [cier_screen()] does not run it, so call it directly.
#'
#' @section What this catches:
#' Predictable / repetitive responding -- straightlining, locking, diagonal-lining,
#' seesaw / alternating patterns. Not useful against uniform random responding (low
#' Laz.R by design); pair it with [cier_irv()], [cier_mahalanobis()], or
#' [cier_psychsyn()].
#'
#' @param responses Numeric matrix (or data.frame / tibble, coerced internally) of
#'   integer-coded responses, one row per respondent, one column per item, in
#'   administration order. `NA` marks a missing response.
#' @param fpr Optional upper-tail mass for the percentile cutoff; `NULL` (default) uses
#'   `0.05`, else a finite number in `(0, 1)`. Mutually exclusive with `cutoff`,
#'   `kneedle`.
#' @param cutoff Optional literal cutoff on the Laz.R score, finite in `[0, 1]`.
#'   Mutually exclusive with `fpr`, `kneedle`.
#' @param kneedle Logical (default `FALSE`); when `TRUE` the cutoff is the
#'   sample-specific Kneedle elbow (Satopaa et al., 2011) of the sorted scores instead
#'   of the `fpr` percentile. Abstains (`NA` cutoff, with a warning) when fewer than
#'   three finite scores remain or all scores are equal. Mutually exclusive with
#'   `fpr`, `cutoff`.
#'
#' @template return-cier-index
#'
#' @references
#' Biemann, T., Koch-Bayram, I. F., Meier-Barthold, M., & Aguinis, H. (2025).
#' Using Markov chains to detect careless responding in survey research.
#' *Organizational Research Methods*, 28(4), 543-568.
#' \doi{10.1177/10944281251334778}
#'
#' Satopaa, V., Albrecht, J., Irwin, D., & Raghavan, B. (2011). Finding a
#' "kneedle" in a haystack: Detecting knee points in system behavior. *2011 31st
#' International Conference on Distributed Computing Systems Workshops*, 166-171.
#'
#' @seealso [cier_autocorrelation()], [cier_longstring()]
#' @family indirect indices
#' @export
#' @examples
#' # The 44 BFI items are the first 44 columns of the bundled example data.
#' out <- cier_lazr(bfi_careless[, 1:44])
#' out
#' head(as.data.frame(out))
cier_lazr <- function(responses, fpr = NULL, cutoff = NULL, kneedle = FALSE) {
  call <- rlang::caller_env()
  responses <- check_responses(responses, call = call)
  check_integer_responses(responses, call = call)
  check_flag(kneedle, "kneedle", call = call)
  use_kneedle <- isTRUE(kneedle)
  # Validate each cutoff knob, then reject any pair (the three are mutually exclusive).
  if (!is.null(fpr)) check_open_unit(fpr, "fpr", call = call)
  if (!is.null(cutoff)) check_number(cutoff, "cutoff", lower = 0, upper = 1, call = call)
  knobs <- list(fpr = fpr, cutoff = cutoff,
                kneedle = if (use_kneedle) TRUE else NULL)
  assert_single_cutoff(knobs, call = call)
  row <- cier_method_spec("cier_lazr")
  value <- kernel_lazr(responses, call = call)
  # Literal `cutoff` or Kneedle elbow is an explicit override; `NULL` lets the shared
  # percentile tail resolve the `fpr` default.
  override <- if (!is.null(cutoff)) {
    cutoff
  } else if (use_kneedle) {
    resolve_kneedle_cutoff(value, call = call)
  } else {
    NULL
  }
  # The Kneedle elbow arrives as `cutoff`, so tag its provenance explicitly.
  provenance <- if (use_kneedle) list(method = "kneedle", rate = NA_real_) else NULL
  resolve_index_cutoff(value, row, fpr, override, provenance = provenance,
                       call = call)
}
