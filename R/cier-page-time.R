# Purpose: cier_page_time() -- the per-page rapid-responding C/IER index, the
#          within-survey member of the timing family.
# Args:    See documentation below.
# Returns: A light `cier_index` (see new_cier_index()).
# Invariants:
#   - Value = per-respondent count of pages whose mean per-item time (page total
#     / items on the page) is strictly below `min_seconds`; direction upper.
#   - Oracle-only trust (no CRAN partner); see tests/reference/TOLERANCES.md.
#   - Cutoff routes through the single resolve_cutoff() path (fixed method): the
#     default `1` (any rapid page), a `frac` fraction of the page count, or a
#     literal `cutoff` count -- the last two mutually exclusive.

#' Per-page rapid-responding C/IER index (Bowling et al. 2023)
#'
#' Counts, per respondent, the number of survey pages answered too fast to be
#' plausible -- pages whose **mean per-item time** falls strictly below
#' `min_seconds` (default `2` seconds per item, following Bowling, Huang, Brower,
#' & Bragg, 2023). The mean per-item time on a page is the respondent's **total**
#' time on that page divided by the number of items on it, so **more** rapid
#' pages flag carelessness (`direction = "upper"`).
#'
#' @details
#' **Input.** `page_seconds` is an n x pages matrix (a data.frame / tibble is
#' coerced) holding each respondent's **total** time on each page, one column per
#' page -- the shape survey platforms export (e.g. a Qualtrics page-submit timer).
#' `items_per_page` gives the number of items on each page, in the same column
#' order; the page total is divided by it to get the per-item rate. (If your
#' input is already a per-item rate, pass `items_per_page = rep(1, ncol)`.) Each
#' observed page time must be strictly positive; a zero or negative duration is a
#' typed error, and a missing time (`NA`) marks an untimed page that contributes
#' no evidence.
#'
#' **Cutoff.** The default flags a respondent with **any** rapid page (the cited
#' Bowling et al. rule: cutoff `1`). On a long survey a single rapid page can be
#' an over-sensitive trigger, so override the count with **one** of two mutually
#' exclusive arguments: `frac`, a fraction of the **total page count** in
#' `(0, 1]` resolving to `ceiling(frac * pages)` (e.g. `frac = 0.1` on 40 pages
#' requires four rapid pages), or `cutoff`, a literal rapid-page count in
#' `[1, pages]`. Respondents whose rapid-page count is at or above the resolved
#' cutoff are flagged.
#'
#' **Abstention.** A page with no recorded time (`NA`) contributes no evidence:
#' it is neither counted rapid nor counted toward the page total, and the per-item
#' denominator is the declared `items_per_page`, not a count of answered cells. A
#' respondent for whom **every** page is untimed abstains: both `value` and
#' `flagged` are `NA` and the row is excluded from the flag count and rate.
#'
#' **This is an absolute-threshold heuristic by design.** The `2` s/item floor is
#' a cited fixed rule, not an empirical percentile -- a respondent below it is
#' called rapid regardless of the rest of the sample. Page time is a *floor* on
#' effort, not a content check.
#'
#' @section What this catches:
#' Within-survey speeding bursts that the total-time index misses -- a respondent
#' who attends to the first part of the survey but races one or more later pages
#' registers those pages as rapid even when the whole-survey total looks plausible.
#'
#' @section What this misses:
#' Slow careless responding (a respondent who clicks deliberately but
#' inattentively never trips the rapid-page threshold). Pair it with a
#' response-content index from the indirect or pattern families, and with
#' [cier_total_time()] for the whole-survey floor.
#'
#' @param page_seconds A numeric matrix (or a data.frame / tibble coerced
#'   internally) of per-page total times in seconds, one row per respondent and
#'   one column per page. Must be strictly positive; `NA` marks an untimed page.
#' @param items_per_page A vector of positive whole numbers, length
#'   `ncol(page_seconds)`, giving the number of items on each page (in column
#'   order).
#' @param min_seconds Single finite number `>= 0`. The per-item time threshold in
#'   seconds; a page whose mean per-item time is strictly below it is counted
#'   rapid. Default `2` per Bowling et al. (2023).
#' @param frac Optional cutoff as a **fraction** of the page count: a single
#'   finite number in `(0, 1]`, resolving to `ceiling(frac * pages)`. Mutually
#'   exclusive with `cutoff`. `NULL` (default) uses the registry default count
#'   of `1`.
#' @param cutoff Optional **literal** rapid-page-count cutoff: a single finite
#'   number in `[1, pages]`, used verbatim. Supplied instead of `frac`.
#'
#' @return A `cier_index`: a list with per-respondent `value` (numeric, `NA` on
#'   abstention) and `flagged` (logical) vectors plus the `method`, `cutoff`, and
#'   `direction` metadata. Use [as.data.frame()][as.data.frame.cier_index] for a
#'   tidy data frame and `print()` for a summary.
#'
#' @references
#' Bowling, N. A., Huang, J. L., Brower, C. K., & Bragg, C. B. (2023). The quick
#' and the careless: The construct validity of page time as a measure of
#' insufficient effort responding to surveys. *Organizational Research Methods*,
#' 26(2), 323-352. \doi{10.1177/10944281211056520}
#'
#' Huang, J. L., Curran, P. G., Keeney, J., Poposki, E. M., & DeShon, R. P.
#' (2012). Detecting and deterring insufficient effort responding to surveys.
#' *Journal of Business and Psychology*, 27(1), 99-114.
#'
#' @seealso [cier_total_time()] for the whole-survey timing floor;
#'   [cier_methods()] for the registry defaults.
#' @family timing indices
#' @export
#' @examples
#' # Five respondents, three two-item pages, per-page total times in seconds.
#' page_seconds <- matrix(c(20, 18, 22,
#'                          19, 21, 20,
#'                           3,  2, 24,
#'                          22, 20, 19,
#'                           2,  2,  2),
#'                        nrow = 5L, byrow = TRUE)
#' cier_page_time(page_seconds, items_per_page = c(2, 2, 2))
#'
#' # Flag only respondents with at least half their pages rapid (here 2 of 3),
#' # rather than any single rapid page:
#' cier_page_time(page_seconds, items_per_page = c(2, 2, 2), frac = 0.5)
cier_page_time <- function(page_seconds, items_per_page, min_seconds = 2,
                           frac = NULL, cutoff = NULL) {
  call <- rlang::caller_env()
  page_seconds <- check_page_seconds(page_seconds, call = call)
  n_pages <- ncol(page_seconds)
  items_per_page <- check_items_per_page(items_per_page, n_pages, call = call)
  check_number(min_seconds, "min_seconds", lower = 0, call = call)
  # The two cutoff overrides mirror longstring: a fraction of the page count
  # (`frac`) and a literal count (`cutoff`), mutually exclusive and validated
  # before the kernel runs.
  if (!is.null(frac)) check_fraction(frac, "frac", call = call)
  if (!is.null(cutoff)) {
    check_number(cutoff, "cutoff", lower = 1, upper = n_pages, call = call)
  }
  assert_single_override(frac, "frac", cutoff, call = call)
  row <- cier_method_row("cier_page_time")
  value <- kernel_page_time(page_seconds, items_per_page, min_seconds)
  # A literal `cutoff` passes through verbatim; `frac` resolves a fraction of the
  # page count; otherwise the registry default count of 1 (any rapid page flags).
  cutoff_value <- if (!is.null(cutoff)) {
    cutoff
  } else if (!is.null(frac)) {
    resolve_cutoff(method = "fixed", value = frac, n_items = n_pages, call = call)
  } else {
    resolve_cutoff(method = "fixed", value = row$default_cutoff_value, call = call)
  }
  flagged <- apply_flag(value, cutoff_value, row$flag_direction, call = call)
  new_cier_index(value, flagged, row$method, cutoff_value, row$flag_direction)
}
