# Warn when every scored respondent is flagged rapid: usually a unit mismatch (times
# in minutes/another unit), or a too-high threshold when `min_seconds` is off its
# default. `min_seconds_default` selects the hint; computed over scored rows only.
warn_page_time_all_flagged <- function(value, flagged, min_seconds,
                                       min_seconds_default,
                                       call = rlang::caller_env()) {
  scored <- !is.na(value)
  if (!any(scored) || !all(flagged[scored])) {
    return(invisible(NULL))
  }
  hint <- if (min_seconds_default) {
    "This usually means the times are not in seconds -- minutes or another \\
     unit make every per-item mean fall below {.arg min_seconds}. Check the \\
     export's unit; {.fun cier_page_time} expects per-page totals in seconds."
  } else {
    "With {.arg min_seconds} set to {min_seconds}, either the times are not in \\
     seconds (a unit mismatch pushes every per-item mean below the threshold) \\
     or {min_seconds} s per item is high relative to your page times. Check the \\
     export's unit -- {.fun cier_page_time} expects per-page totals in seconds \\
     -- and that {.arg min_seconds} is the threshold you intend."
  }
  cier_warn(
    "cier_warning_all_flagged",
    c("Every scored respondent was flagged as rapid.", "i" = hint),
    data = list(n_scored = sum(scored)), call = call
  )
  invisible(NULL)
}

#' Per-page rapid-responding C/IER index (Bowling et al. 2023)
#'
#' Counts, per respondent, the survey pages answered too fast to be plausible: pages
#' whose mean per-item time (page total divided by item count) falls below
#' `min_seconds` (default `2` s/item, following Bowling, Huang, Brower, & Bragg, 2023).
#' More rapid pages flag carelessness (`direction = "upper"`).
#'
#' @details
#' `page_seconds` is an n x pages matrix (data.frame / tibble coerced) of per-page
#' total times in seconds; `items_per_page` gives the declared item count per page
#' (same column order, not the number of answered cells) and divides the page total to a
#' per-item rate (pass `rep(1, ncol)` for an already per-item rate). A page time must be
#' non-negative (`0` is maximal speeding evidence; a negative errors); `NA` marks an
#' untimed page contributing no evidence, and a respondent with every page untimed
#' abstains (`value`, `flagged` `NA`).
#'
#' Times must be in **seconds**: a millisecond export flags no one, a minutes export
#' flags everyone (with a warning). The `2` s/item floor is an absolute rule, not an
#' empirical percentile -- a respondent below it is rapid regardless of the sample.
#' This is a floor on effort, not a content check.
#'
#' This index is standalone: [cier_screen()] does not run it, so call it directly.
#'
#' The default cutoff flags any rapid page (`cutoff = 1`), which can be over-sensitive
#' on a long survey; override with `frac` (a fraction of the page count) or `cutoff`
#' (a literal count).
#'
#' @section What this catches:
#' Within-survey speeding bursts that [cier_total_time()]'s whole-survey floor misses
#' -- racing later pages while the whole-survey total still looks plausible.
#'
#' @section What this misses:
#' Slow careless responding, which never trips the rapid-page threshold. Pair it with a
#' response-content index and [cier_total_time()] for the whole-survey floor.
#'
#' @param page_seconds Numeric n x pages matrix (data.frame / tibble coerced) of
#'   per-page total times in seconds. Non-negative (`0` is valid speeding evidence,
#'   negative errors); `NA` marks an untimed page.
#' @param items_per_page Positive whole numbers, length `ncol(page_seconds)`, giving
#'   the item count per page in column order.
#' @param min_seconds Single finite number `>= 0`; per-item time threshold below which
#'   a page is counted rapid. Default `2` per Bowling et al. (2023).
#' @param frac Optional cutoff as a fraction of the page count, in `(0, 1]`, resolving
#'   to `ceiling(frac * pages)`. `NULL` (default) uses count `1`. Mutually exclusive
#'   with `cutoff`.
#' @param cutoff Optional literal rapid-page count, a whole number in `[1, pages]`.
#'   Mutually exclusive with `frac`.
#'
#' @template return-cier-index
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
#' @seealso [cier_total_time()] for the whole-survey timing floor.
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
  # `frac` and `cutoff` are mutually exclusive; validate before resolving.
  if (!is.null(frac)) check_fraction(frac, "frac", call = call)
  if (!is.null(cutoff)) {
    check_number(cutoff, "cutoff", lower = 1, upper = n_pages, whole = TRUE,
                 call = call)
  }
  assert_single_cutoff(list(frac = frac, cutoff = cutoff), call = call)
  row <- cier_method_spec("cier_page_time")
  value <- kernel_page_time(page_seconds, items_per_page, min_seconds)
  # Literal `cutoff` verbatim; `frac` as a page-count fraction; else default count 1.
  cutoff_value <- if (!is.null(cutoff)) {
    cutoff
  } else if (!is.null(frac)) {
    resolve_cutoff(method = "fixed", value = frac, n_items = n_pages, call = call)
  } else {
    resolve_cutoff(method = "fixed", value = row$default_cutoff_value, call = call)
  }
  # Flag explicitly (not flag_and_assemble) so the all-flagged warning can read
  # `flagged`; `min_seconds == 2` signals the default hint.
  flagged <- apply_flag(value, cutoff_value, row$flag_direction, call = call)
  warn_page_time_all_flagged(value, flagged, min_seconds, min_seconds == 2, call = call)
  new_cier_index_provenance(value, flagged, row$method, cutoff_value,
                            row$flag_direction, cutoff,
                            if (is.null(frac)) "fixed_count" else "fixed_fraction",
                            if (is.null(frac)) NA_real_ else frac)
}
