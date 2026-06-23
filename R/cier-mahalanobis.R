# Compose the cli message for the degenerate-covariance warning. The cause is pasted
# in (not cli-interpolated) because it is local to this helper, outside the
# environment cli sees when the wrapper raises the warning.
mahalanobis_abstain_message <- function(status) {
  cause <- switch(status,
    singular_covariance   = "the covariance matrix is singular",
    indefinite_covariance = "the pairwise covariance is not positive definite",
    "fewer than two respondents have data"
  )
  msg <- c(paste0("Mahalanobis distance is undefined: ", cause, "."),
           "i" = "All respondents abstain (value = {.val NA}); no one is flagged.")
  if (identical(status, "indefinite_covariance")) {
    msg <- c(msg,
             "i" = paste0("Heavy or structured missingness can make the pairwise ",
                          "covariance inconsistent; consider scoring complete ",
                          "cases or reducing missingness."))
  }
  msg
}

# Warn once when partial (mean-imputed) rows are scored, so a cier_screen() user who
# never reads the help page learns of it. Inert on abstaining statuses.
warn_partial_rows <- function(value, responses, call = rlang::caller_env()) {
  scored <- !is.na(value)
  n_partial <- sum(scored & rowSums(is.na(responses)) > 0L)
  if (n_partial == 0L) {
    return(invisible(NULL))
  }
  n_scored <- sum(scored)
  cier_warn(
    "cier_warning_partial_rows",
    c("{n_partial} of {n_scored} scored respondent{?s} had one or more missing \\
       items, mean-imputed before scoring.",
      "i" = "The observed terms are still weighted by the full inverse covariance \\
             and referred to the same chi-square reference, so partial rows are \\
             over-flagged under inter-item correlation.",
      "i" = "Re-run on complete cases as a sensitivity check."),
    data = list(n_partial = n_partial, n_scored = n_scored), call = call
  )
  invisible(NULL)
}

#' Mahalanobis-distance C/IER index
#'
#' Computes each respondent's squared Mahalanobis distance D² of their response
#' vector from the sample centroid (Curran, 2016), with covariance estimated
#' pairwise. A profile far from the multivariate centre has a large D², so high
#' values flag carelessness.
#'
#' @details
#' **Cutoff.** D² is referred to a chi-square reference (a real null, not a sample
#' percentile) with `p` degrees of freedom (`p` = number of items): the cutoff is
#' `qchisq(1 - alpha, df = p)` with `alpha = 0.001` (Tabachnick & Fidell's
#' conservative `p < .001` criterion). The flag rate is therefore not pinned at
#' `alpha` and can run higher under contamination. Override with one of `alpha`
#' (chi-square tail probability) or `cutoff` (literal threshold on D²); `df = p`
#' regardless.
#'
#' **Partial rows.** A respondent's missing items are mean-imputed so the row is
#' scored, but the observed terms are still weighted by the full inverse covariance
#' and referred to the same `χ²(p)`, which over-flags partial rows under inter-item
#' correlation. When any scored row carries a mean-imputed item the index warns
#' (`cier_warning_partial_rows`); re-run on complete cases as a sensitivity check.
#'
#' **Abstention.** An all-`NA` respondent abstains (`value` and `flagged` are `NA`).
#' If no distance is computable for anyone -- fewer than two respondents carry data,
#' singular covariance, or pairwise covariance indefinite under heavy/structured
#' missingness -- every `value` is `NA`, a warning is emitted, and no one is flagged.
#'
#' @section What this catches:
#' Multivariate outliers: response profiles inconsistent with the sample's covariance
#' structure, including some random and patterned responding. Being relative to this
#' sample's centroid, a largely careless sample shifts the centroid and weakens it.
#'
#' @template responses
#' @param alpha Optional chi-square upper-tail probability; `NULL` (default) uses
#'   `0.001`, else a number in `(0, 1)` giving cutoff `qchisq(1 - alpha, df = p)`.
#'   Mutually exclusive with `cutoff`.
#' @param cutoff Optional literal cutoff on D², a finite number `>= 0`; respondents
#'   at or above it are flagged.
#'
#' @template return-cier-index
#'
#' @references
#' Curran, P. G. (2016). Methods for the detection of carelessly invalid
#' responses in survey data. *Journal of Experimental Social Psychology*, 66,
#' 4–19. \doi{10.1016/j.jesp.2015.07.006}
#'
#' Tabachnick, B. G., & Fidell, L. S. (2014). *Using Multivariate Statistics*
#' (6th ed.). Pearson New International Edition. (Conservative `p < .001`
#' criterion for multivariate outliers via Mahalanobis distance.)
#'
#' @family indirect indices
#' @export
#' @examples
#' # The 44 BFI items are the first 44 columns of the bundled example data.
#' out <- cier_mahalanobis(bfi_careless[, 1:44])
#' out
#' head(as.data.frame(out))
cier_mahalanobis <- function(responses, alpha = NULL, cutoff = NULL) {
  call <- rlang::caller_env()
  responses <- check_responses(responses, call = call)
  if (!is.null(alpha)) check_open_unit(alpha, "alpha", call = call)
  if (!is.null(cutoff)) check_number(cutoff, "cutoff", lower = 0, call = call)
  assert_single_cutoff(list(alpha = alpha, cutoff = cutoff), call = call)
  row <- cier_method_spec("cier_mahalanobis")
  res <- kernel_mahalanobis(responses)
  if (!identical(res$status, "ok")) {
    cier_warn("cier_warning_singular_covariance",
              mahalanobis_abstain_message(res$status),
              data = list(reason = res$status), call = call)
  }
  warn_partial_rows(res$value, responses, call = call)
  cutoff_value <- if (!is.null(cutoff)) {
    cutoff
  } else {
    resolve_cutoff(method = row$default_cutoff_method, df = ncol(responses),
                   alpha = if (is.null(alpha)) row$default_cutoff_value else alpha,
                   call = call)
  }
  flag_and_assemble(res$value, cutoff_value, row, cutoff, "chisq",
                    if (is.null(alpha)) row$default_cutoff_value else alpha,
                    call = call)
}
