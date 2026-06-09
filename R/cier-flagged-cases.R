# Purpose: cier_flagged_cases() -- extract the row positions of the flagged
#          respondents from a cier_index or a cier_screen. A pure accessor over
#          existing output: no statistics, no cutoffs.
# Args:    See documentation below.
# Returns: An integer vector of 1-based respondent row positions (ascending).
# Invariants:
#   - cier_index: which(flagged) -- abstaining (NA) respondents are excluded.
#   - cier_screen: which(rowSums(votes) >= min_votes) on the COLLAPSED votes, so
#     it agrees with the agreement diagnostic and never double-counts the
#     consistency construct.

#' Extract the flagged respondents from a cier index or screen
#'
#' Returns the **row positions** (1-based, into the original `responses`) of the
#' respondents a `cier_index` flagged, or -- for a [cier_screen()] -- of the
#' respondents flagged by at least `min_votes` *constructs*.
#'
#' @details
#' For a single index, a respondent is flagged when their value falls past the
#' index's cutoff; abstaining respondents (value `NA`) are never flagged.
#'
#' For a screen, the threshold counts **collapsed votes**, not raw per-index
#' flags: correlated indices that share a construct (even-odd and personal
#' reliability fuse into one `consistency` vote) count once. `min_votes` therefore
#' matches the cross-index agreement diagnostic exactly --
#' `length(cier_flagged_cases(screen, min_votes = k))` equals the screen's
#' "flagged by >= k votes" count. The function returns positions only; it imposes
#' no careless/not label and excludes no one -- the researcher chooses `min_votes`
#' and subsets `responses` themselves (e.g.
#' `responses[cier_flagged_cases(screen, min_votes = 3), ]`).
#'
#' @param x A `cier_index` (from an index function such as [cier_irv()]) or a
#'   [cier_screen()].
#' @param min_votes For a `cier_screen`, the minimum number of construct-votes a
#'   respondent must receive to be returned. A positive whole number; the default
#'   `1` returns everyone flagged by at least one construct. Ignored for a single
#'   `cier_index`.
#' @param ... Ignored.
#'
#' @return An integer vector of respondent row positions (ascending), empty
#'   (`integer(0)`) when no respondent qualifies.
#'
#' @seealso [cier_screen()]; the index functions [cier_irv()], [cier_longstring()]
#' @export
#' @examples
#' # From a single index: the respondents IRV flagged.
#' irv <- cier_irv(bfi_careless[, 1:44])
#' head(cier_flagged_cases(irv))
#'
#' # From a screen: respondents flagged by at least 3 constructs.
#' nm <- names(bfi_careless)[1:44]
#' items <- data.frame(scale = gsub("^v_BFI_|[0-9_R]+$", "", nm),
#'                     reverse_keyed = grepl("_R$", nm), categories = 5L)
#' screen <- cier_screen(bfi_careless[, 1:44], items,
#'                       control = list(cier_personal_reliability = list(seed = 1)))
#' cier_flagged_cases(screen, min_votes = 3)
cier_flagged_cases <- function(x, ...) {
  UseMethod("cier_flagged_cases")
}

#' @rdname cier_flagged_cases
#' @export
cier_flagged_cases.cier_index <- function(x, ...) {
  # NA (abstaining) rows are dropped by which(); unname() guarantees bare
  # positional indices even when the input carried respondent row names.
  unname(which(x$flagged))
}

#' @rdname cier_flagged_cases
#' @export
cier_flagged_cases.cier_screen <- function(x, min_votes = 1L, ...) {
  check_count(min_votes, "min_votes", call = rlang::caller_env())
  # rowSums of the collapsed votes is the per-respondent construct count. An empty
  # (0-vote) screen has a 0-column votes frame, whose rowSums is all-zero, so no
  # respondent reaches min_votes (>= 1) -- no special case needed. unname() keeps
  # the result bare positional indices regardless of input row names.
  counts <- rowSums(as.matrix(x$votes))
  unname(which(counts >= min_votes))
}

#' @rdname cier_flagged_cases
#' @export
cier_flagged_cases.default <- function(x, ...) {
  cier_abort(
    "cier_error_input",
    c("{.arg x} must be a {.cls cier_index} or {.cls cier_screen} object.",
      "x" = "Got {.cls {class(x)}}."),
    data = list(arg = "x"), call = rlang::caller_env()
  )
}
