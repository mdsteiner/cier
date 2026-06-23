# cier_attention(): the attention-check index. Per-respondent count of checks
# failed among those answered; abstains when no check was answered.

# Warn (not error: a never-endorsed value is legal) when a check's pass-set shares
# no value with the column's observed responses, so every answered respondent fails
# it -- almost always a mis-specified pass-set.
warn_pass_disjoint <- function(checks, pass, call = rlang::caller_env()) {
  offending <- vapply(seq_len(ncol(checks)), function(j) {
    col <- checks[, j]
    obs <- col[!is.na(col)]
    length(obs) > 0L && !any(obs %in% pass[[j]])
  }, logical(1L))
  if (!any(offending)) {
    return(invisible(NULL))
  }
  idx <- which(offending)
  cn <- colnames(checks)
  labels <- if (is.null(cn)) {
    paste("column", idx)
  } else {
    ifelse(is.na(cn[idx]) | !nzchar(cn[idx]), paste("column", idx), cn[idx])
  }
  cier_warn(
    "cier_warning_disjoint_pass",
    c("{cli::qty(length(labels))}The pass-set{?s} for check{?s} {.val {labels}} \\
       {?matches/match} none of {?that column's/those columns'} observed \\
       responses, so every answered respondent fails {?it/them}.",
      "i" = "A pass-set lists the PASSING values in the column's own coding -- \\
             e.g. {.code pass = list(0)} on 1-5 data endorses nothing and flags \\
             everyone. Check each pass-set against its column's coding \\
             (all-NA columns are unaffected)."),
    data = list(arg = "pass", observed = labels), call = call
  )
}

# Warn when every scored respondent is flagged: at 100% a mis-specified pass-set
# is the usual cause. Computed over scored rows only.
warn_attention_all_flagged <- function(value, flagged,
                                       call = rlang::caller_env()) {
  scored <- !is.na(value)
  if (any(scored) && all(flagged[scored])) {
    cier_warn(
      "cier_warning_all_flagged",
      c("Every scored respondent was flagged on the attention checks.",
        "i" = "This usually means a mis-specified {.arg pass} -- the wrong \\
               passing values, the wrong column order, or values outside a \\
               column's coding make every respondent fail. Check {.arg pass} \\
               matches each check's observed responses."),
      data = list(n_scored = sum(scored)), call = call
    )
  }
  invisible(NULL)
}

#' Attention-check C/IER index (Meade & Craig 2012)
#'
#' Counts, per respondent, how many attention checks were failed among those
#' answered. A check is failed when the answered response is not one of its
#' declared passing values, so more failures flag carelessness
#' (`direction = "upper"`). Attention checks are items with an author-known
#' passing response: instructed-response, bogus, and infrequency items.
#'
#' @details
#' Any answered response outside `pass[[j]]` is a failure. Pass values use the
#' literal observed-response coding -- no reverse-keying, no `items` metadata -- so
#' checks with different response scales mix freely. An unanswered check (`NA`) is
#' neither pass nor failure; a respondent who answered none abstains (`value` and
#' `flagged` are `NA`).
#'
#' This index is standalone: [cier_screen()] does not run it, so call it directly.
#'
#' @section What this catches:
#' Random responders, and agreement / extreme straightliners whose constant answer
#' crosses a check's non-endorsing set.
#'
#' @section What this misses:
#' A careless respondent who answers every check correctly (a straightliner whose
#' value matches each pass-set, or someone who spots the checks). Pair direct checks
#' with a response-pattern or timing index.
#'
#' @param checks A numeric matrix (data.frame / tibble coerced internally) of raw
#'   attention-check responses, one row per respondent, one column per check. `NA`
#'   marks an unanswered check.
#' @param pass A length-`ncol(checks)` list in column order; `pass[[j]]` is a
#'   non-empty numeric vector of the value(s) that pass check `j`.
#' @param cutoff Optional whole-number failed-check count in `[1, ncol(checks)]`;
#'   `NULL` (default) flags any failed check.
#'
#' @template return-cier-index
#'
#' @references
#' Meade, A. W., & Craig, S. B. (2012). Identifying careless responses in survey
#' data. *Psychological Methods*, 17(3), 437-455. \doi{10.1037/a0028085}
#'
#' Goldammer P., Stöckli, P. L., Escher, Y. A., Annen, H., Jonas, K., &
#' Antonakis, J. (2024) Careless responding detection revisited: Accuracy of
#' direct and indirect measures. *Behavior Research Methods*, 56, 8422-8449.
#'
#' @seealso [bfi_careless] for the bundled bogus and instructed-response check
#'   columns.
#' @family direct indices
#' @export
#' @examples
#' # Two checks: a bogus item (pass = 1 or 2) and an instructed item (pass = 0).
#' checks <- matrix(c(1, 0,
#'                    3, 0,
#'                    2, 5,
#'                    4, 7,
#'                    1, 0),
#'                  nrow = 5L, byrow = TRUE)
#' cier_attention(checks, pass = list(c(1, 2), 0))
#'
#' # The bundled Bruhlmann data ships two attention checks: bogus (pass 1 or 2)
#' # and instructed v_IRI (pass 0).
#' att <- cier_attention(bfi_careless[, c("v_Bogus_Item", "v_IRI")],
#'                       pass = list(c(1, 2), 0))
#' sum(att$flagged)   # respondents failing at least one of the two checks
cier_attention <- function(checks, pass, cutoff = NULL) {
  call <- rlang::caller_env()
  checks <- check_responses(checks, arg = "checks", call = call)
  n_checks <- ncol(checks)
  pass <- check_pass(pass, n_checks, check_names = colnames(checks), call = call)
  if (!is.null(cutoff)) {
    check_number(cutoff, "cutoff", lower = 1, upper = n_checks, whole = TRUE,
                 call = call)
  }
  warn_pass_disjoint(checks, pass, call = call)
  row <- cier_method_spec("cier_attention")
  value <- kernel_attention(checks, pass)
  # Literal `cutoff` verbatim; otherwise the package default count of 1.
  cutoff_value <- if (!is.null(cutoff)) {
    cutoff
  } else {
    resolve_cutoff(method = "fixed", value = row$default_cutoff_value, call = call)
  }
  # Flag explicitly (not flag_and_assemble): all-flagged warning reads `flagged`
  # before assembly.
  flagged <- apply_flag(value, cutoff_value, row$flag_direction, call = call)
  warn_attention_all_flagged(value, flagged, call = call)
  new_cier_index_provenance(value, flagged, row$method, cutoff_value,
                            row$flag_direction, cutoff, "fixed_count", NA_real_)
}
