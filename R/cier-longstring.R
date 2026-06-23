#' Longest-run-length C/IER index
#'
#' Computes the maximum run length of consecutive identical responses for every
#' respondent (Johnson, 2005; Meade & Craig, 2012). Long runs flag straightlining.
#' The score is computed over the raw response row (no scale blocking).
#'
#' @details
#' Default cutoff is `ceiling(0.5 * p)` (`p` = item count) -- the "greater than
#' half the scale length" heuristic of Meade & Craig (2012), flagging a longest
#' run `>=` it. It is not a calibrated false-positive rate and is lenient on a
#' multi-scale instrument. Override with `frac` or `cutoff`.
#'
#' The half-length default is over-sensitive below ~10 items, where a moderate run
#' is a large fraction of the scale; when it resolves to 2 or below (`p <= 4`) the
#' index warns. An all-`NA` row abstains (`value` and `flagged` are `NA`).
#'
#' @section What this catches:
#' Straightlining (locking on one option). Not random or alternating responding,
#' which give short runs.
#'
#' @template responses
#' @param frac Optional cutoff as a fraction of the item count: a finite number in
#'   `(0, 1]`, resolving to `ceiling(frac * p)`. `NULL` (default) uses `0.5`.
#'   Mutually exclusive with `cutoff`.
#' @param cutoff Optional literal run-length cutoff: a finite number in `[1, p]`.
#'   Mutually exclusive with `frac`.
#'
#' @template return-cier-index
#'
#' @references
#' Johnson, J. A. (2005). Ascertaining the validity of individual protocols from
#' web-based personality inventories. *Journal of Research in Personality*,
#' 39(1), 103–129. \doi{10.1016/j.jrp.2004.09.009}
#'
#' Meade, A. W., & Craig, S. B. (2012). Identifying careless responses in survey
#' data. *Psychological Methods*, 17(3), 437–455. \doi{10.1037/a0028085}
#'
#' @family indirect indices
#' @export
#' @examples
#' # 44 BFI items; default frac = 0.5 is a run of 22 -- lenient on a multi-scale
#' # instrument, so compare a tighter `frac`.
#' out <- cier_longstring(bfi_careless[, 1:44])
#' out
#' head(as.data.frame(out))
#'
#' cier_flagged_cases(out)
#'
#' cier_longstring(bfi_careless[, 1:44], frac = 0.25) # cutoff = 11
cier_longstring <- function(responses, frac = NULL, cutoff = NULL) {
  call <- rlang::caller_env()
  responses <- check_responses(responses, call = call)
  p <- ncol(responses)
  # `frac` in (0, 1]; `cutoff` in [1, p]; mutually exclusive.
  if (!is.null(frac)) check_fraction(frac, "frac", call = call)
  if (!is.null(cutoff)) check_number(cutoff, "cutoff", lower = 1, upper = p,
                                     call = call)
  assert_single_cutoff(list(frac = frac, cutoff = cutoff), call = call)
  row <- cier_method_spec("cier_longstring")
  value <- kernel_longstring(responses)
  value[rowSums(!is.na(responses)) == 0L] <- NA_real_   # abstain on all-NA rows
  # Literal cutoff passes through; else resolve the fraction of items.
  cutoff_value <- if (!is.null(cutoff)) {
    cutoff
  } else {
    f <- if (is.null(frac)) row$default_cutoff_value else frac
    resolve_cutoff(method = "fixed", value = f, n_items = p, call = call)
  }
  if (is.null(frac) && is.null(cutoff) && cutoff_value <= 2) {
    warn_short_battery(p, cutoff_value, call = call)
  }
  flag_and_assemble(value, cutoff_value, row, cutoff, "fixed_fraction",
                    if (is.null(frac)) row$default_cutoff_value else frac,
                    call = call)
}


# Warn when the half-length default cutoff lands at 2 or below (<= 4 items), where
# the "greater than half" heuristic flags many attentive respondents. Default path
# only; score unchanged.
warn_short_battery <- function(p, cutoff_value, call = rlang::caller_env()) {
  cier_warn(
    "cier_warning_short_battery",
    c("The default longstring cutoff resolved to {cutoff_value} run-length on \\
       {p} item{?s} (`ceiling(0.5 * p)`).",
      "i" = "On short batteries the half-length default flags a large share of \\
             attentive respondents (about 20-56% below ~10 items). Set \\
             {.arg cutoff} or {.arg frac} deliberately, or inspect the score \\
             distribution."),
    data = list(n_items = p, cutoff = cutoff_value), call = call
  )
}
