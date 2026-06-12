# Purpose: cier_longstring() -- the public maximum-run-length C/IER index.
# Args:    See documentation below.
# Returns: A light `cier_index` (see new_cier_index()).
# Invariants:
#   - Bytewise compatible with careless::longstring() on present data.
#   - The cutoff routes through the single resolve_cutoff() path (fixed method):
#     `frac` resolves a fraction of the item count, `cutoff` a literal count.

#' Longest-run-length C/IER index
#'
#' Computes the **maximum** run length of consecutive identical responses for
#' every respondent — the canonical longstring statistic of Johnson (2005) and
#' Meade & Craig (2012). Long runs flag straightlining. The score is computed
#' over the raw response row (no scale blocking).
#'
#' @details
#' **Cutoff.** The default cutoff is `ceiling(0.5 * p)` where `p` is the number
#' of items — the "greater than half the scale length" starting point suggested
#' by Meade & Craig (2012), applied here to the full item set. It is a heuristic,
#' **not** a calibrated false-positive rate: longstring has no universal cutoff
#' ([careless::longstring()] leaves the choice to inspection), and half of *all*
#' items is lenient for a multi-scale instrument — consider a smaller `frac` or
#' `cutoff` and check the score distribution. Override the default with **one** of
#' two mutually exclusive arguments: `frac`, a fraction of the item count in
#' `(0, 1]` (resolving to `ceiling(frac * p)`), or `cutoff`, a literal run-length
#' count in `[1, p]`. Respondents whose longest run is greater than or equal to
#' the resolved cutoff are flagged.
#'
#' **Abstention.** A respondent who answered nothing (an all-`NA` row) abstains:
#' both `value` and `flagged` are `NA` and the row is excluded from the flag
#' count and rate. (`careless::longstring()` returns `1` for such a row; `cier`
#' reports it as missing rather than as the least-careless score.) The flag rate
#' shown by `print()` is taken over respondents who produced a usable score.
#'
#' @section What this catches:
#' Straightlining, midpoint responding, and extreme responding (when a careless
#' respondent locks on one option). It is **not** useful against random
#' responding, which produces short runs.
#'
#' @param responses A numeric matrix (or a data.frame / tibble coerced
#'   internally) of responses, one row per respondent and one column per item.
#'   `NA` marks a missing response.
#' @param frac Optional cutoff as a **fraction** of the item count: a single
#'   finite number in `(0, 1]`, resolving to `ceiling(frac * p)`. `NULL`
#'   (default) uses the registry default `0.5`. Mutually exclusive with `cutoff`.
#' @param cutoff Optional **literal** run-length cutoff: a single finite number
#'   in `[1, p]`, used verbatim. Supplied instead of `frac`.
#'
#' @return A `cier_index`: a list with per-respondent `value` (numeric, `NA` on
#'   abstention) and `flagged` (logical) vectors plus the `method`, `cutoff`, and
#'   `direction` metadata. Use [as.data.frame()][as.data.frame.cier_index] for a
#'   tidy data frame and `print()` for a summary.
#'
#' @references
#' Johnson, J. A. (2005). Ascertaining the validity of individual protocols from
#' web-based personality inventories. *Journal of Research in Personality*,
#' 39(1), 103–129. \doi{10.1016/j.jrp.2004.09.009}
#'
#' Meade, A. W., & Craig, S. B. (2012). Identifying careless responses in survey
#' data. *Psychological Methods*, 17(3), 437–455. \doi{10.1037/a0028085}
#'
#' @seealso [careless::longstring()]
#' @family indirect indices
#' @export
#' @examples
#' # The 44 BFI items are the first 44 columns of the bundled example data.
#' out <- cier_longstring(bfi_careless[, 1:44])
#' out
#' head(as.data.frame(out))
cier_longstring <- function(responses, frac = NULL, cutoff = NULL) {
  call <- rlang::caller_env()
  responses <- check_responses(responses, call = call)
  p <- ncol(responses)
  # Validate every input up front so a bad argument fails before the kernel runs.
  # `frac` is a fraction of the item count in (0, 1]; a literal `cutoff` is a
  # run-length count in [1, p]; the two overrides are mutually exclusive.
  if (!is.null(frac)) check_fraction(frac, "frac", call = call)
  if (!is.null(cutoff)) check_number(cutoff, "cutoff", lower = 1, upper = p,
                                     call = call)
  assert_single_override(frac, "frac", cutoff, call = call)
  row <- cier_method_row("cier_longstring")
  value <- kernel_longstring(responses)
  value[rowSums(!is.na(responses)) == 0L] <- NA_real_   # abstain on all-NA rows
  # A literal cutoff passes through verbatim (already validated); otherwise the
  # default is a fraction of the item count, ceiling(frac * p).
  cutoff_value <- if (!is.null(cutoff)) {
    cutoff
  } else {
    f <- if (is.null(frac)) row$default_cutoff_value else frac
    resolve_cutoff(method = "fixed", value = f, n_items = p, call = call)
  }
  flagged <- apply_flag(value, cutoff_value, row$flag_direction, call = call)
  new_cier_index(value, flagged, row$method, cutoff_value, row$flag_direction)
}
