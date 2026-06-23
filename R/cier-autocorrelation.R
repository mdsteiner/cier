# Resolve and validate the lag range against the item count. `max_lag = NULL`
# defaults to `min(n_items - 3, 10)`: the low lags Gottfried et al. recommend,
# capped (default-only; an explicit `max_lag` may reach `n_items - 3`) so a long
# battery's short high-lag slices do not saturate the score. Fewer than four items
# has no valid lag, so it fails first with the item-count error.
resolve_autocorrelation_lags <- function(min_lag, max_lag, n_items, call) {
  check_count(min_lag, "min_lag", call = call)
  if (n_items < 4L) {
    cier_abort(
      "cier_error_input",
      c("{.fn cier_autocorrelation} needs at least 4 items.",
        "x" = "Got {n_items} item{?s}; the shortest lag (lag 1) then leaves a \\
               slice of {n_items - 1L}, but an autocorrelation needs at least 3.",
        "i" = "Score each same-format block separately, or use a different index \\
               on a very short scale."),
      data = list(arg = "responses", observed = n_items, expected = 4L),
      call = call
    )
  }
  if (is.null(max_lag)) {
    max_lag <- min(n_items - 3L, 10L)
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
#' Computes each respondent's maximum absolute lag autocorrelation -- the largest
#' absolute Pearson correlation between their raw response row (administration order,
#' no reverse-keying or scale blocking) and its lag-shifted self over a range of lags
#' (Gottfried, Jezek, Kralova & Rihacek, 2022). A predictable repeating pattern
#' (seesaw, diagonal, straightlining) yields high autocorrelation, so high values
#' flag carelessness.
#'
#' @details
#' **Cutoff.** The default flags the upper `fpr` quantile of observed scores (95th
#' percentile), the paper's relative-ranking recommendation; or pass an absolute
#' `cutoff` in `[0, 1]` instead.
#'
#' **Zero-variance convention.** A constant lag slice (a straightliner) has an
#' undefined correlation and scores `1` rather than abstaining, so it is flagged.
#'
#' **Missing data.** With `na_rm = FALSE` (default) missingness is pairwise within
#' each lag; a lag with fewer than three complete pairs abstains (a two-pair
#' correlation is `+/-1` by construction), and a row with too few scorable lags
#' abstains entirely (`value` and `flagged` `NA`). `na_rm = TRUE` strips each row's
#' missing responses before lagging, collapsing administration-order spacing, and is
#' rarely appropriate.
#'
#' **Assumptions.** Columns are read as an ordered sequence on a common response
#' scale, so the index is meaningful only when they are in administration order with
#' the same number of options; score each same-format block separately on a
#' mixed-format survey.
#'
#' Standalone: [cier_screen()] does not run it, so call it directly.
#'
#' @section What this catches:
#' Repetitive / periodic careless responding -- seesaw, alternating, diagonal, and
#' (via the zero-variance convention) straightlining. Not useful against random
#' responding, whose autocorrelation is near zero.
#'
#' @param responses Numeric matrix (or data.frame / tibble coerced internally), one
#'   row per respondent, one column per item, in administration order; `NA` marks a
#'   missing response.
#' @param min_lag Positive whole number; smallest lag to evaluate. Default `1`.
#' @param max_lag Positive whole number or `NULL` (default, resolves to
#'   `min(ncol(responses) - 3, 10)`); largest lag to evaluate. An explicit value must
#'   satisfy `min_lag <= max_lag <= ncol(responses) - 3`.
#' @param na_rm Single `TRUE` / `FALSE`. `TRUE` strips each row's missing responses
#'   before lagging; `FALSE` (default) handles missingness pairwise within each lag.
#' @param fpr Optional false-positive tail mass for the percentile cutoff; finite
#'   number in `(0, 1)`, `NULL` (default) uses `0.05`. Mutually exclusive with
#'   `cutoff`.
#' @param cutoff Optional **literal** cutoff on the autocorrelation magnitude, in
#'   `[0, 1]`, supplied instead of `fpr`; scores at or above it are flagged.
#'
#' @template return-cier-index
#'
#' @references
#' Gottfried, J., Jezek, S., Kralova, M., & Rihacek, T. (2022). Autocorrelation
#' screening: A potentially efficient method for detecting repetitive response
#' patterns in questionnaire data. *Practical Assessment, Research & Evaluation*,
#' 27, Article 2. \doi{10.7275/vyxb-gt24}
#'
#' @family indirect indices
#' @export
#' @examples
#' # First 44 columns are the BFI items; max_lag defaults to min(ncol - 3, 10) = 10.
#' out <- cier_autocorrelation(bfi_careless[, 1:44])
#' out
#' head(as.data.frame(out))
#' # Restrict to a handful of low lags:
#' cier_autocorrelation(bfi_careless[, 1:44], max_lag = 6L)
cier_autocorrelation <- function(responses, min_lag = 1L, max_lag = NULL,
                                 na_rm = FALSE, fpr = NULL, cutoff = NULL) {
  call <- rlang::caller_env()
  responses <- check_responses(responses, call = call)
  check_flag(na_rm, "na_rm", call = call)
  check_percentile_overrides(fpr, cutoff, lower = 0, upper = 1, call = call)
  lags <- resolve_autocorrelation_lags(min_lag, max_lag, ncol(responses), call)
  row <- cier_method_spec("cier_autocorrelation")
  value <- kernel_autocorrelation(responses, lags$min_lag, lags$max_lag, na_rm)
  resolve_index_cutoff(value, row, fpr, cutoff, call = call)
}
