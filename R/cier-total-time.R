# cier_total_time(): the total survey completion-time index. The per-respondent value
# is the validated seconds vector itself (no kernel math); direction lower -- low
# totals flag speeding.

#' Total survey completion-time C/IER index
#'
#' Flags respondents with the fastest total completion times as likely speeders
#' (Ward & Meade, 2023; Huang et al., 2012). The value is the completion time in
#' seconds itself, so low totals flag carelessness (`direction = "lower"`); there is
#' no response matrix.
#'
#' @details
#' `seconds` is a bare numeric vector of whole-survey durations, one per respondent
#' (a 2-D input is rejected). Each must be strictly positive (zero/negative is a typed
#' error; `NA` abstains, yielding `NA` `value`/`flagged`). Percentile and
#' median-relative cutoffs are scale-invariant, but a literal `cutoff` is in seconds,
#' so confirm the export's unit. This index is standalone: [cier_screen()] does not
#' run it.
#'
#' The default flags the fastest respondents at the lower `fpr` quantile (5th
#' percentile) -- a ranking convention flagging at least `fpr` of respondents, not a
#' calibrated false-positive rate; it needs at least 20 finite times to resolve a
#' tail. `frac_median` instead flags respondents at or below that fraction of the
#' sample median (Leiner's 2019 RSI uses `0.5`; Greszki et al. 2015 use 0.5 / 0.4 /
#' 0.3), robust to up to half the sample responding carelessly. If no respondent has a
#' finite time the cutoff is `NA`, nobody is flagged, and a typed warning is raised.
#'
#' @section What this catches:
#' Respondents who rushed the whole survey -- the canonical speeder profile and the
#' strongest single signal on fully-careless protocols (Ward & Meade, 2023).
#'
#' @section What this misses:
#' Within-survey speeding bursts (a respondent who races only the second half keeps a
#' plausible total) and careless respondents who took their time -- completion time is
#' a floor on effort, not a content check. Biemann et al. (2025) report 46-80% of
#' whole-survey straightliners pass a 2 s/item screen, so pair completion time with a
#' response-content index ([cier_irv()] or [cier_mahalanobis()]) and a page-level
#' timing index when available.
#'
#' @param seconds A numeric vector of total completion times in seconds, one per
#'   respondent. Must be strictly positive; `NA` abstains.
#' @param fpr Optional percentile tail mass; `NULL` (default) uses `0.05`, else a
#'   finite number in `(0, 1)`. Mutually exclusive with `cutoff` and `frac_median`.
#' @param cutoff Optional literal cutoff in seconds, a finite number `>= 0`; flags
#'   respondents at or below it. Mutually exclusive with `fpr` and `frac_median`.
#' @param frac_median Optional median-relative cutoff in `(0, 1]`, resolving to
#'   `frac_median * median` of the finite times; flags respondents at or below it.
#'   Mutually exclusive with `fpr` and `cutoff`.
#'
#' @template return-cier-index
#'
#' @references
#' Ward, M. K., & Meade, A. W. (2023). Dealing with careless responding in survey
#' data. *Annual Review of Psychology*, 74, 577-596.
#' \doi{10.1146/annurev-psych-040422-045007}
#'
#' Huang, J. L., Curran, P. G., Keeney, J., Poposki, E. M., & DeShon, R. P.
#' (2012). Detecting and deterring insufficient effort responding to surveys.
#' *Journal of Business and Psychology*, 27(1), 99-114.
#'
#' Leiner, D. J. (2019). Too fast, too straight, too weird: Non-reactive
#' indicators for meaningless data in internet surveys. *Survey Research Methods*,
#' 13(3), 229-248.
#'
#' Greszki, R., Meyer, M., & Schoen, H. (2015). Exploring the effects of removing
#' "too fast" responses and respondents from web surveys. *Public Opinion
#' Quarterly*, 79(2), 471-503.
#'
#' @seealso [cier_irv()] and [cier_mahalanobis()] for response-content indices to
#'   pair with completion time.
#' @family timing indices
#' @export
#' @examples
#' # A careful bulk of several-minute completions with a fast tail of speeders.
#' secs <- c(420, 412, 505, 388, 470, 401, 440, 415, 380, 460,
#'           405, 390, 430, 365, 480, 445, 398, 455, 425, 372,
#'           95, 47, 88, 110)
#' cier_total_time(secs)
#'
#' # Flag respondents faster than half the sample median (Leiner 2019 RSI):
#' cier_total_time(secs, frac_median = 0.5)
cier_total_time <- function(seconds, fpr = NULL, cutoff = NULL,
                            frac_median = NULL) {
  call <- rlang::caller_env()
  value <- check_seconds(seconds, call = call)
  # Validate each cutoff knob, then reject any pair (mutually exclusive).
  if (!is.null(fpr)) check_open_unit(fpr, "fpr", call = call)
  if (!is.null(frac_median)) check_fraction(frac_median, "frac_median", call = call)
  if (!is.null(cutoff)) check_number(cutoff, "cutoff", lower = 0, call = call)
  knobs <- list(fpr = fpr, frac_median = frac_median, cutoff = cutoff)
  assert_single_cutoff(knobs, call = call)
  row <- cier_method_spec("cier_total_time")
  # Literal or median-relative cutoff is an explicit override; else the percentile
  # tail resolves the `fpr` default.
  override <- if (!is.null(cutoff)) {
    cutoff
  } else if (!is.null(frac_median)) {
    resolve_median_cutoff(value, frac_median, call = call)
  } else {
    NULL
  }
  # A median-relative override arrives as `cutoff`; pass provenance explicitly so the
  # resolver does not label it "literal".
  provenance <- if (!is.null(frac_median)) {
    list(method = "median_relative", rate = frac_median)
  } else {
    NULL
  }
  resolve_index_cutoff(value, row, fpr, override, provenance = provenance,
                       call = call)
}
