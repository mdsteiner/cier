# Purpose: cier_autocorrelation() -- the public maximum-absolute-lag
#          autocorrelation index (Gottfried, Jezek, Kralova & Rihacek 2022),
#          packaged in responsePatterns::rp.acors().
# Args:    See documentation below.
# Returns: A light `cier_index` (see new_cier_index()).
# Invariants:
#   - Equal to responsePatterns::rp.acors()$indices$max.abs.ac within 1e-10 on
#     complete data (the masked-sum kernel vs their per-row cor()).
#   - The cutoff routes through the single resolve_index_cutoff() path
#     (percentile method, upper direction, fpr = 0.05 by default).

# Resolve and validate the lag range against the item count. `max_lag = NULL`
# defaults to `n_items - 3` (the rp.acors() default: the largest lag whose slice
# still holds three observations). An explicit lag range must satisfy
# `1 <= min_lag <= max_lag <= n_items - 3`, else a typed input error -- so a
# battery with fewer than four items (no valid lag) also fails here, at the
# boundary, before the kernel runs.
resolve_autocorrelation_lags <- function(min_lag, max_lag, n_items, call) {
  check_count(min_lag, "min_lag", call = call)
  if (is.null(max_lag)) {
    max_lag <- n_items - 3L
  }
  check_count(max_lag, "max_lag", call = call)
  min_lag <- as.integer(min_lag)
  max_lag <- as.integer(max_lag)
  if (min_lag > max_lag || max_lag > n_items - 3L) {
    cier_abort(
      "cier_error_input",
      c("The lag range cannot be reconciled with the item count.",
        "x" = "Need 1 <= min_lag <= max_lag <= ncol(responses) - 3 = \\
               {n_items - 3L}.",
        "i" = "Got min_lag = {min_lag}, max_lag = {max_lag}."),
      data = list(arg = "max_lag", observed = c(min_lag, max_lag),
                  expected = n_items - 3L),
      call = call
    )
  }
  list(min_lag = min_lag, max_lag = max_lag)
}

#' Autocorrelation C/IER index
#'
#' Computes each respondent's **maximum absolute lag autocorrelation** -- the
#' largest absolute Pearson correlation between their response vector and its
#' lag-shifted self over a range of lags (Gottfried, Jezek, Kralova & Rihacek,
#' 2022), as packaged in `responsePatterns::rp.acors()`. A respondent whose
#' answers follow a predictable repeating pattern (seesaw / alternating,
#' diagonal, or straightlining) has a high autocorrelation, so **high** values
#' flag carelessness. The score is computed over the **raw** response row in
#' administration order -- there is no reverse-keying and no scale blocking.
#'
#' @details
#' **Cutoff.** The default flags the most autocorrelated respondents: the cutoff
#' is the empirical upper `fpr` quantile of the observed scores (the 95th
#' percentile by default) and respondents at or above it are flagged. This
#' empirical-ranking convention is the paper's own recommendation ("rely on
#' relative, not absolute criteria"). Adjust the target tail with `fpr`, or pass
#' an absolute `cutoff` in `[0, 1]` to flag on a literal autocorrelation
#' magnitude; `fpr` and `cutoff` are mutually exclusive.
#'
#' **Lag range.** `max_lag = NULL` resolves to `ncol(responses) - 3`, the
#' `rp.acors()` default (the largest lag whose slice still holds three
#' observations). On a long instrument this includes many short, high-lag slices
#' whose correlation saturates near 1; Gottfried et al. recommend restricting
#' attention to a handful of low lags (roughly 5-12) on runs of 10-40 items, so
#' set `max_lag` (and `min_lag`) deliberately for a long battery.
#'
#' **Zero-variance convention.** A lag whose slice is constant (a straightliner)
#' has an undefined correlation; following `rp.acors()` it scores `1` -- the top
#' of the tail -- rather than abstaining. This differs from `stats::acf()` (which
#' returns `NaN`). A straightliner is therefore **flagged**, not dropped.
#'
#' **Missing data.** With `na_rm = FALSE` (default) missingness is handled
#' pairwise within each lag (a missing cell drops only the lagged pairs it
#' touches). With `na_rm = TRUE` each row's missing responses are removed before
#' lagging; this collapses the administration-order spacing the index depends on
#' and is rarely appropriate. A respondent with no scorable lag (every lag has
#' too few observations) abstains: both `value` and `flagged` are `NA` and the
#' row is excluded from the flag count and rate.
#'
#' **Assumptions.** The index reads the response columns as an ordered sequence
#' on a common response scale, so it is meaningful only when the columns are in
#' administration order and share the same number of response options. This is
#' the source paper's own scale-format guidance: "only questions with the same
#' answer scales should be analyzed at one time", since mixing scales "with
#' vastly different number ranges ... can bias the results to a great extent",
#' and questions with unique scales or answer options (say gender or education)
#' "should be excluded prior to screening" (Gottfried et al., 2022). On a
#' mixed-format survey, score each same-format block in a separate call.
#'
#' @section What this catches:
#' Repetitive / periodic careless responding -- seesaw, alternating, diagonal,
#' and (via the zero-variance convention) straightlining -- whose lagged response
#' is predictable from the previous. It is **not** useful against random
#' responding, whose autocorrelation is near zero.
#'
#' @param responses A numeric matrix (or a data.frame / tibble coerced
#'   internally) of responses, one row per respondent and one column per item, in
#'   administration order. `NA` marks a missing response.
#' @param min_lag Positive whole number. The smallest lag to evaluate. Default
#'   `1`.
#' @param max_lag Positive whole number or `NULL`. The largest lag to evaluate.
#'   `NULL` (default) resolves to `ncol(responses) - 3`. An explicit value must
#'   satisfy `min_lag <= max_lag <= ncol(responses) - 3`. The default is sized
#'   from the full column count even when `na_rm = TRUE`, so a heavily missing
#'   row may have its higher lags abstain once its responses are stripped.
#' @param na_rm Single `TRUE` / `FALSE`. When `TRUE`, each row's missing
#'   responses are stripped before lagging; when `FALSE` (default) missingness is
#'   handled pairwise within each lag.
#' @param fpr Optional target false-positive tail mass for the percentile cutoff.
#'   `NULL` (default) uses the registry default `0.05`. A finite number in the
#'   open interval `(0, 1)`; the cutoff is that upper-tail quantile of the
#'   observed scores. Mutually exclusive with `cutoff`.
#' @param cutoff Optional **literal** cutoff on the autocorrelation magnitude,
#'   supplied instead of `fpr`. A single finite number in `[0, 1]`; respondents
#'   whose score is at or above it are flagged.
#'
#' @return A `cier_index`: a list with per-respondent `value` (numeric, `NA` on
#'   abstention) and `flagged` (logical) vectors plus the `method`, `cutoff`, and
#'   `direction` metadata. Use [as.data.frame()][as.data.frame.cier_index] for a
#'   tidy data frame and `print()` for a summary.
#'
#' @references
#' Gottfried, J., Jezek, S., Kralova, M., & Rihacek, T. (2022). Autocorrelation
#' screening: A potentially efficient method for detecting repetitive response
#' patterns in questionnaire data. *Practical Assessment, Research & Evaluation*,
#' 27, Article 2. \doi{10.7275/vyxb-gt24}
#'
#' @seealso [responsePatterns::rp.acors()]
#' @family indirect indices
#' @export
#' @examples
#' # The 44 BFI items are the first 44 columns of the bundled example data.
#' # Restrict to low lags (the paper's recommendation): the default
#' # max_lag = ncol - 3 includes short high-lag slices whose correlation
#' # saturates near 1 on a long battery, flagging far more than the target rate.
#' out <- cier_autocorrelation(bfi_careless[, 1:44], max_lag = 8L)
#' out
#' head(as.data.frame(out))
cier_autocorrelation <- function(responses, min_lag = 1L, max_lag = NULL,
                                 na_rm = FALSE, fpr = NULL, cutoff = NULL) {
  call <- rlang::caller_env()
  responses <- check_responses(responses, call = call)
  check_flag(na_rm, "na_rm", call = call)
  # `value` is a max |Pearson autocorrelation| in [0, 1]; `fpr` is a tail mass in
  # (0, 1), a literal `cutoff` a magnitude in [0, 1]; they are mutually exclusive.
  check_percentile_overrides(fpr, cutoff, lower = 0, upper = 1, call = call)
  lags <- resolve_autocorrelation_lags(min_lag, max_lag, ncol(responses), call)
  row <- cier_method_row("cier_autocorrelation")
  value <- kernel_autocorrelation(responses, lags$min_lag, lags$max_lag, na_rm)
  resolve_index_cutoff(value, row, fpr, cutoff, call = call)
}
