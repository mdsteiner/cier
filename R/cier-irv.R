# Purpose: cier_irv() -- the public intra-individual response variability index.
# Args:    See documentation below.
# Returns: A light `cier_index` (see new_cier_index()).
# Invariants:
#   - Equal to careless::irv() within 1e-10 (matrixStats::rowSds vs stats::sd).
#   - The cutoff routes through the single resolve_cutoff() path (percentile
#     method, lower direction, fpr = 0.05 by default).

#' Intra-individual response variability (IRV) C/IER index
#'
#' Computes each respondent's **intra-individual response variability** -- the
#' sample standard deviation of their answers across items (Marjanovic et al.,
#' 2015). A respondent who barely varies their answers (straightlining, midpoint
#' responding) has a low IRV, so **low** values flag carelessness. The score is
#' computed over the raw response row (no scale blocking); items left blank are
#' dropped, and a respondent who answered fewer than two items abstains.
#'
#' @details
#' **Cutoff.** The default flags the lowest-variability respondents: the cutoff
#' is the empirical `fpr` quantile of the observed IRV scores (the 5th
#' percentile by default) and respondents at or below it are flagged. This is a
#' **ranking** convention -- a sample percentile flags **at least** `fpr` of
#' respondents (more when scores tie at the cutoff) -- not a calibrated
#' false-positive rate. Adjust the target tail
#' with `fpr`, or pass an absolute `cutoff` instead to flag on a literal IRV
#' threshold (e.g. one carried over from a calibration sample). `fpr` and
#' `cutoff` are mutually exclusive.
#'
#' **Abstention.** A respondent with fewer than two answered items has no
#' defined sample standard deviation: both `value` and `flagged` are `NA` and
#' the row is excluded from the flag count and rate. A constant respondent who
#' answered at least two items scores `0` -- the strongest low-variability
#' signal -- and is flagged; that is not abstention.
#'
#' @section What this catches:
#' Straightlining and midpoint responding, which compress the response range. It
#' is **not** useful against random responding, which inflates variability.
#'
#' @param responses A numeric matrix (or a data.frame / tibble coerced
#'   internally) of responses, one row per respondent and one column per item.
#'   `NA` marks a missing response.
#' @param fpr Optional target false-positive tail mass for the percentile
#'   cutoff. `NULL` (default) uses the registry default `0.05`. A finite number
#'   in the open interval `(0, 1)`; the cutoff is that lower-tail quantile of the
#'   observed scores. Mutually exclusive with `cutoff`.
#' @param cutoff Optional **literal** cutoff on the IRV score, supplied instead
#'   of `fpr`. A single finite number; respondents whose IRV is at or below it
#'   are flagged (the lower-tail direction). Use it to apply an absolute
#'   threshold rather than a sample percentile.
#'
#' @return A `cier_index`: a list with per-respondent `value` (numeric, `NA` on
#'   abstention) and `flagged` (logical) vectors plus the `method`, `cutoff`, and
#'   `direction` metadata. Use [as.data.frame()][as.data.frame.cier_index] for a
#'   tidy data frame and `print()` for a summary.
#'
#' @references
#' Marjanovic, Z., Holden, R., Struthers, W., Cribbie, R., & Greenglass, E.
#' (2015). The inter-item standard deviation (ISD): An index that discriminates
#' between conscientious and random responders. *Personality and Individual
#' Differences*, 84, 79–83. \doi{10.1016/j.paid.2014.08.021}
#'
#' @seealso [careless::irv()]
#' @family indirect indices
#' @export
#' @examples
#' # The 44 BFI items are the first 44 columns of the bundled example data.
#' out <- cier_irv(bfi_careless[, 1:44])
#' out
#' head(as.data.frame(out))
cier_irv <- function(responses, fpr = NULL, cutoff = NULL) {
  call <- rlang::caller_env()
  responses <- check_responses(responses, call = call)
  # `fpr` is a tail mass in (0, 1); a literal IRV `cutoff` is a non-negative SD
  # threshold (0 flags only exact straightliners); they are mutually exclusive.
  check_percentile_overrides(fpr, cutoff, lower = 0, call = call)
  row <- cier_method_row("cier_irv")
  # rowSds(na.rm = TRUE) returns NA where fewer than two items were answered, so
  # abstention needs no separate guard (cf. longstring's all-NA row).
  value <- kernel_irv(responses)
  resolve_index_cutoff(value, row, fpr, cutoff, call = call)
}
