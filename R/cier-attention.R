# Purpose: cier_attention() -- the attention-check C/IER index, the only member
#          of the direct family in v0.2.
# Args:    See documentation below.
# Returns: A light `cier_index` (see new_cier_index()).
# Invariants:
#   - Value = per-respondent count of attention checks FAILED among ANSWERED ones
#     (an answered response not in that check's pass-set); direction upper.
#   - NA on a check is no evidence; a respondent with no answered check abstains.
#   - Oracle-only trust (no CRAN partner); see tests/reference/TOLERANCES.md.
#   - Cutoff is the fixed default 1 (any failed check) or a literal `cutoff` count.

#' Attention-check C/IER index (Meade & Craig 2012)
#'
#' Counts, per respondent, how many **attention checks** were failed among those
#' answered. Attention checks are items with an author-known correct or
#' non-endorsing response -- instructed-response items ("to show you are reading,
#' select Strongly Disagree"), bogus items (a nonsensical statement no attentive
#' respondent endorses), and infrequency items (a statement almost everyone
#' answers one way). A respondent **fails** a check when their answered response
#' is **not** one of that check's declared passing values, so **more** failed
#' checks flag carelessness (`direction = "upper"`).
#'
#' @details
#' **Input.** `checks` is an n x k matrix (a data.frame / tibble is coerced)
#' holding the raw responses to the k attention-check items, one column per check.
#' `pass` is a list of length k, in the same column order; element `j` gives the
#' value(s) that **pass** check `j` -- the directed option for an instructed item
#' (e.g. `0`), the non-endorsing options for a bogus item (e.g. `c(1, 2)` on a
#' 1-5 scale), or every value but the infrequent one for an infrequency item. Any
#' answered response outside `pass[[j]]` is a failure.
#'
#' **No reverse-keying.** The pass values are read in the **literal observed-response
#' coding**: declare them as the survey platform delivers them. `cier_attention()`
#' applies no reverse-keying and takes no `items` metadata. It also makes no
#' same-scale assumption across checks -- each check has its own pass-set, so
#' checks with different response scales mix freely.
#'
#' **Cutoff.** The default flags a respondent with **any** failed check (the cited
#' Meade & Craig rule: cutoff `1`). To require more failures, pass `cutoff`, a
#' literal failed-check count in `[1, k]` (e.g. `cutoff = 2` flags only respondents
#' who failed at least two checks). Respondents whose failed-check count is at or
#' above the cutoff are flagged.
#'
#' **Abstention.** An unanswered check (`NA`) contributes no evidence: it is
#' counted as neither a pass nor a failure. A respondent who answered **at least
#' one** check scores a finite count; a respondent who answered **none** abstains
#' (`value` and `flagged` are `NA`, and the row is excluded from the flag count and
#' rate). This is a deliberate choice -- a missing check is treated as absent
#' information, not as a failure.
#'
#' @section What this catches:
#' Random responders (likely to miss an instructed key or endorse a bogus item);
#' agreement / extreme straightliners whose constant answer crosses a check's
#' non-endorsing set; and any respondent inattentive enough to overlook a planted
#' instruction. Direct checks are the bundle anchor when a survey includes them.
#'
#' @section What this misses:
#' A careless respondent who happens to answer every check correctly (a
#' straightliner whose single value matches each pass-set; someone who spots and
#' complies with the checks while answering the rest carelessly). Direct checks are
#' easy to pass with low effort once noticed, and must be designed into the survey
#' in advance. Pair them with at least one response-pattern (indirect / pattern)
#' or timing index.
#'
#' @param checks A numeric matrix (or a data.frame / tibble coerced internally) of
#'   raw attention-check responses, one row per respondent and one column per
#'   check. `NA` marks an unanswered check.
#' @param pass A list of length `ncol(checks)`, in column order, giving the
#'   passing value(s) of each check: `pass[[j]]` is a non-empty numeric vector of
#'   the response value(s) that pass check `j`. An answered response outside it is
#'   a failure.
#' @param cutoff Optional **literal** failed-check-count cutoff: a single finite
#'   number in `[1, ncol(checks)]`, used verbatim. `NULL` (default) uses the
#'   registry default count of `1` (any failed check flags).
#'
#' @return A `cier_index`: a list with per-respondent `value` (numeric, `NA` on
#'   abstention) and `flagged` (logical) vectors plus the `method`, `cutoff`, and
#'   `direction` metadata. Use [as.data.frame()][as.data.frame.cier_index] for a
#'   tidy data frame and `print()` for a summary.
#'
#' @references
#' Meade, A. W., & Craig, S. B. (2012). Identifying careless responses in survey
#' data. *Psychological Methods*, 17(3), 437-455. \doi{10.1037/a0028085}
#'
#' Goldammer P., Stöckli, P. L., Escher, Y. A., Annen, H., Jonas, K., &
#' Antonakis, J. (2024) Careless responding detection revisited: Accuracy of
#' direct and indirect measures. *Behavior Research Methods*, 56, 8422-8449.
#'
#' @seealso [cier_methods()] for the registry defaults; [bfi_careless] for the
#'   bundled bogus and instructed-response check columns.
#' @family direct indices
#' @export
#' @examples
#' # Five respondents, two checks: a bogus item (pass = the non-endorsing 1 or 2)
#' # and an instructed item (pass = the directed option 0).
#' checks <- matrix(c(1, 0,
#'                    3, 0,
#'                    2, 5,
#'                    4, 7,
#'                    1, 0),
#'                  nrow = 5L, byrow = TRUE)
#' cier_attention(checks, pass = list(c(1, 2), 0))
#'
#' # The bundled Bruhlmann data ships two attention checks: a bogus item (passing
#' # responses 1 or 2) and an instructed item v_IRI (passing response 0).
#' att <- cier_attention(bfi_careless[, c("v_Bogus_Item", "v_IRI")],
#'                       pass = list(c(1, 2), 0))
#' sum(att$flagged)   # respondents failing at least one of the two checks
cier_attention <- function(checks, pass, cutoff = NULL) {
  call <- rlang::caller_env()
  checks <- check_responses(checks, arg = "checks", call = call)
  n_checks <- ncol(checks)
  pass <- check_pass(pass, n_checks, call = call)
  if (!is.null(cutoff)) {
    check_number(cutoff, "cutoff", lower = 1, upper = n_checks, call = call)
  }
  row <- cier_method_row("cier_attention")
  value <- kernel_attention(checks, pass)
  # A literal `cutoff` passes through verbatim; otherwise the registry default
  # count of 1 (any failed check flags).
  cutoff_value <- if (!is.null(cutoff)) {
    cutoff
  } else {
    resolve_cutoff(method = "fixed", value = row$default_cutoff_value, call = call)
  }
  flagged <- apply_flag(value, cutoff_value, row$flag_direction, call = call)
  new_cier_index(value, flagged, row$method, cutoff_value, row$flag_direction)
}
