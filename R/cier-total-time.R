# Purpose: cier_total_time() -- the public total survey completion-time C/IER
#          index, the coarsest member of the timing family.
# Args:    See documentation below.
# Returns: A light `cier_index` (see new_cier_index()).
# Invariants:
#   - The per-respondent value is the validated seconds vector itself (an
#     identity; there is no kernel math), direction lower -- low totals flag
#     speeding.
#   - Oracle-only trust (no CRAN partner); see tests/reference/TOLERANCES.md.
#   - Three mutually-exclusive cutoff knobs: the percentile `fpr` default, the
#     median-relative `frac_median`, and a literal `cutoff`.

#' Total survey completion-time C/IER index
#'
#' Takes each respondent's **total completion time** in seconds -- one duration
#' per respondent, as survey platforms export it -- and flags the fastest as
#' likely speeders (Ward & Meade, 2023; Huang et al., 2012). The value is the
#' completion time itself, so **low** totals flag carelessness
#' (`direction = "lower"`); there is no response matrix and no transformation of
#' the input.
#'
#' @details
#' **Input.** `seconds` is a bare numeric vector, one total per respondent. A
#' two-dimensional input (a matrix or data frame) is rejected -- sum your
#' per-cell response times to one total per respondent first. Each time must be
#' strictly positive; a zero or negative duration is a typed error, and a missing
#' duration (`NA`) abstains.
#'
#' **Cutoff.** The default flags the fastest respondents: the cutoff is the
#' empirical lower `fpr` quantile of the observed times (the 5th percentile by
#' default), and respondents at or below it are flagged. This is a **ranking**
#' convention -- a sample percentile flags `fpr` of respondents by construction --
#' not a calibrated false-positive rate. Override it with **one** of three
#' mutually exclusive arguments:
#' - `fpr`, a different target tail mass for the percentile;
#' - `frac_median`, a **median-relative** rule flagging respondents faster than
#'   that fraction of the sample median (the published instrument-free speeding
#'   rules: Leiner's 2019 Relative Speed Index flags those slower than half the
#'   median, i.e. `frac_median = 0.5`; Greszki et al. 2015 use 0.5 / 0.4 / 0.3).
#'   Because it is anchored to the median it is robust to up to half the sample
#'   responding carelessly, unlike the empirical percentile;
#' - `cutoff`, a literal threshold in seconds.
#' Respondents whose time is at or below the resolved cutoff are flagged.
#'
#' **Abstention.** A respondent with a missing time (`NA`) abstains: both `value`
#' and `flagged` are `NA` and the row is excluded from the flag count and rate.
#' If no respondent has a finite time the cutoff cannot be resolved: it is `NA`,
#' nobody is flagged, and a typed warning is raised.
#'
#' **Assumptions.** The index reads `seconds` as whole-survey durations on a
#' common clock; it cannot see *within*-survey speeding (a respondent who answers
#' the first half attentively and races the second half keeps a plausible total).
#' Biemann et al. (2025) report that 46-80% of whole-survey straightliners pass a
#' 2 s/item screen, so completion time must be paired with a response-content
#' index. Pair it with [cier_irv()] or [cier_mahalanobis()], and with a
#' page-level timing index when page times are available.
#'
#' @section What this catches:
#' Respondents who rushed the whole survey -- the canonical speeder profile, the
#' strongest single signal on fully-careless protocols (Ward & Meade, 2023).
#'
#' @section What this misses:
#' Within-survey speeding bursts, and careless respondents who took their time.
#' Completion time is a *floor* on effort, not a content check.
#'
#' @param seconds A numeric vector of total completion times in seconds, one per
#'   respondent. Must be strictly positive; `NA` marks a missing duration (which
#'   abstains).
#' @param fpr Optional target false-positive tail mass for the percentile cutoff.
#'   `NULL` (default) uses the registry default `0.05`. A finite number in the
#'   open interval `(0, 1)`; the cutoff is that lower-tail quantile of the
#'   observed times. Mutually exclusive with `cutoff` and `frac_median`.
#' @param cutoff Optional **literal** cutoff in seconds, supplied instead of a
#'   rate. A single finite number `>= 0`; respondents whose time is at or below it
#'   are flagged. Mutually exclusive with `fpr` and `frac_median`.
#' @param frac_median Optional **median-relative** cutoff: a single finite number
#'   in `(0, 1]`, resolving to `frac_median * median` of the finite times.
#'   Respondents at or below it are flagged. Mutually exclusive with `fpr` and
#'   `cutoff`.
#'
#' @return A `cier_index`: a list with per-respondent `value` (numeric, `NA` on
#'   abstention) and `flagged` (logical) vectors plus the `method`, `cutoff`, and
#'   `direction` metadata. Use [as.data.frame()][as.data.frame.cier_index] for a
#'   tidy data frame and `print()` for a summary.
#'
#' @references
#' Ward, M. K., & Meade, A. W. (2023). Dealing with careless responding in survey
#' data. *Annual Review of Psychology*, 74, 577-606.
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
#'   pair with completion time; [cier_methods()] for the registry defaults.
#' @family timing indices
#' @export
#' @examples
#' # One total completion time in seconds per respondent.
#' secs <- c(420, 380, 95, 510, 47, 460, 405, 88, 390, 430)
#' cier_total_time(secs)
#'
#' # Flag respondents faster than half the sample median (Leiner 2019 RSI):
#' cier_total_time(secs, frac_median = 0.5)
cier_total_time <- function(seconds, fpr = NULL, cutoff = NULL,
                            frac_median = NULL) {
  call <- rlang::caller_env()
  value <- check_seconds(seconds, call = call)
  # Three mutually-exclusive cutoff knobs: a target tail mass (`fpr`), a
  # median-relative fraction (`frac_median`), and a literal seconds threshold
  # (`cutoff`). Validate each up front, then reject any combination.
  if (!is.null(fpr)) check_open_unit(fpr, "fpr", call = call)
  if (!is.null(frac_median)) check_fraction(frac_median, "frac_median", call = call)
  if (!is.null(cutoff)) check_number(cutoff, "cutoff", lower = 0, call = call)
  knobs <- list(fpr = fpr, frac_median = frac_median, cutoff = cutoff)
  assert_single_cutoff(knobs, call = call)
  row <- cier_method_row("cier_total_time")
  # A literal `cutoff`, or a resolved median-relative cutoff, is an explicit
  # override; otherwise (`NULL`) the shared percentile tail resolves the `fpr`
  # default. Either way the cutoff -> flag -> assemble step is the one shared
  # resolve_index_cutoff() the other percentile indices use.
  override <- if (!is.null(cutoff)) {
    cutoff
  } else if (!is.null(frac_median)) {
    resolve_median_cutoff(value, frac_median, call = call)
  } else {
    NULL
  }
  resolve_index_cutoff(value, row, fpr, override, call = call)
}
