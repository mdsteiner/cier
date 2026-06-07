# Purpose: cier_mahalanobis() -- the public squared-Mahalanobis-distance index,
#          plus the small message helper for its degenerate-covariance warning.
# Args:    See documentation below.
# Returns: A light `cier_index` (see new_cier_index()).
# Invariants:
#   - Equal to careless::mahad() and psych::outlier() within 1e-10 (both D^2).
#   - The cutoff routes through the single resolve_cutoff() path (chisq method,
#     upper direction, df = item count, alpha = 0.001 by default).

# Compose the cli message for the degenerate-covariance warning. Both causes
# share the cier_warning_singular_covariance class; the leading line names which
# one occurred so the message (not only the class) is informative. The cause is
# pasted in (not cli-interpolated) because it is local to this helper, outside
# the environment cli sees when the wrapper raises the warning.
mahalanobis_abstain_message <- function(status) {
  cause <- if (identical(status, "singular_covariance")) {
    "the covariance matrix is singular"
  } else {
    "fewer than two respondents have data"
  }
  c(paste0("Mahalanobis distance is undefined: ", cause, "."),
    "i" = "All respondents abstain (value = {.val NA}); no one is flagged.")
}

#' Mahalanobis-distance C/IER index
#'
#' Computes each respondent's **squared Mahalanobis distance** D² of their
#' response vector from the sample centroid (Curran, 2016). A response profile
#' far from the multivariate centre -- an inconsistent or atypical pattern across
#' items -- has a large D², so **high** values flag carelessness. The covariance
#' is estimated **pairwise**, and a respondent's missing items drop out of the
#' distance rather than discarding the whole row.
#'
#' @details
#' **Cutoff.** The default uses a **proper reference distribution** rather than a
#' sample percentile: D² is referred to a chi-square distribution with `p`
#' degrees of freedom (`p` = the number of items), and respondents beyond the
#' upper tail are flagged -- the cutoff is `qchisq(1 - alpha, df = p)` with
#' `alpha = 0.001`. Because this is a real null (not a ranking convention), the
#' flag rate is **not** pinned at `alpha`: it is whatever share of respondents
#' exceed the χ²(`p`) quantile, and can be informatively higher under
#' contamination. The `0.001` default follows Tabachnick and Fidell's (2019)
#' recommendation of a conservative `p < .001` criterion for identifying
#' multivariate outliers by Mahalanobis distance. It is **stricter** than
#' `careless::mahad()`, whose default flags the upper 1% tail (its
#' `confidence = 0.99`, i.e. `alpha = 0.01`); pass `alpha = 0.01` to reproduce
#' the `careless` default. Override the cutoff with **one** of two mutually
#' exclusive arguments: `alpha` (the chi-square tail probability) or `cutoff`
#' (a literal threshold on D²).
#'
#' **Abstention.** A respondent who answered nothing (an all-`NA` row) abstains:
#' both `value` and `flagged` are `NA`. If **no** distance can be computed for
#' anyone -- fewer than two respondents carry any data, or the covariance matrix
#' is singular (for example more items than respondents, or perfectly collinear
#' items) -- every `value` is `NA`, a warning is emitted, and no one is flagged.
#'
#' @section What this catches:
#' Multivariate outliers: response profiles inconsistent with the sample's
#' covariance structure, including some random and patterned responding. It is a
#' relative measure (distance from *this* sample's centroid), so a sample that is
#' largely careless shifts the centroid and weakens it.
#'
#' @param responses A numeric matrix (or a data.frame / tibble coerced
#'   internally) of responses, one row per respondent and one column per item.
#'   `NA` marks a missing response.
#' @param alpha Optional chi-square upper-tail probability for the cutoff. `NULL`
#'   (default) uses the registry default `0.001`. A single finite number in
#'   `(0, 1)`; the cutoff is `qchisq(1 - alpha, df = p)`. Mutually exclusive with
#'   `cutoff`.
#' @param cutoff Optional **literal** cutoff on the squared distance, supplied
#'   instead of `alpha`. A single finite number `>= 0`; respondents whose D² is
#'   at or above it are flagged.
#'
#' @return A `cier_index`: a list with per-respondent `value` (numeric, `NA` on
#'   abstention) and `flagged` (logical) vectors plus the `method`, `cutoff`, and
#'   `direction` metadata. Use [as.data.frame()][as.data.frame.cier_index] for a
#'   tidy data frame and `print()` for a summary.
#'
#' @references
#' Curran, P. G. (2016). Methods for the detection of carelessly invalid
#' responses in survey data. *Journal of Experimental Social Psychology*, 66,
#' 4–19. \doi{10.1016/j.jesp.2015.07.006}
#'
#' Tabachnick, B. G., & Fidell, L. S. (2019). *Using Multivariate Statistics*
#' (7th ed.). Pearson. (Conservative `p < .001` criterion for multivariate
#' outliers via Mahalanobis distance.)
#'
#' @seealso [careless::mahad()], [psych::outlier()]
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
  # Validate every input up front so a bad argument fails before the kernel runs.
  # `alpha` is a chi-square tail probability in the open interval (0, 1) -- 0
  # flags nobody, 1 flags everybody; `cutoff` is a literal threshold on D^2
  # (>= 0); the two overrides are mutually exclusive.
  if (!is.null(alpha)) check_open_unit(alpha, "alpha", call = call)
  if (!is.null(cutoff)) check_number(cutoff, "cutoff", lower = 0, call = call)
  assert_single_override(alpha, "alpha", cutoff, call = call)
  row <- cier_method_row("cier_mahalanobis")
  res <- kernel_mahalanobis(responses)
  if (!identical(res$status, "ok")) {
    cier_warn("cier_warning_singular_covariance",
              mahalanobis_abstain_message(res$status),
              data = list(reason = res$status), call = call)
  }
  # A literal cutoff passes through verbatim (already validated); otherwise the
  # chi-square default is qchisq(1 - alpha, df = item count).
  cutoff_value <- if (!is.null(cutoff)) {
    cutoff
  } else {
    resolve_cutoff(method = row$default_cutoff_method, df = ncol(responses),
                   alpha = if (is.null(alpha)) row$default_cutoff_value else alpha,
                   call = call)
  }
  flagged <- apply_flag(res$value, cutoff_value, row$flag_direction, call = call)
  new_cier_index(res$value, flagged, row$method, cutoff_value, row$flag_direction)
}
