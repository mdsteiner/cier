#' Extract the flagged respondents from a cier index or screen
#'
#' Returns the 1-based row positions (into the original `responses`) of the
#' respondents a `cier_index` flagged, or -- for a [cier_screen()] -- of those
#' flagged by at least `min_votes` *constructs*.
#'
#' @details
#' For a single index, a respondent is flagged when their value falls past the
#' cutoff; abstaining respondents (value `NA`) are never flagged.
#'
#' For a screen the threshold counts collapsed votes, not raw per-index flags:
#' correlated indices sharing a construct (even-odd and personal reliability fuse
#' into one `consistency` vote) count once.
#'
#' @param x A `cier_index` (e.g. from [cier_irv()]) or a [cier_screen()].
#' @param min_votes For a `cier_screen`, the minimum construct-votes a respondent
#'   must receive to be returned. Positive whole number, default `1`. Ignored for
#'   a single `cier_index`.
#' @param ... Ignored.
#'
#' @return Integer vector of respondent row positions (ascending), `integer(0)`
#'   when none qualify.
#'
#' @seealso [cier_screen()], [cier_irv()], [cier_longstring()]
#' @export
#' @examples
#' irv <- cier_irv(bfi_careless[, 1:44])
#' head(cier_flagged_cases(irv))
#'
#' # Seed both randomised pieces (RPR, Gnormed's Monte-Carlo null) for reproducibility.
#' nm <- names(bfi_careless)[1:44]
#' items <- data.frame(scale = sub("^v_BFI_([A-Za-z]+)[0-9].*$", "\\1", nm),
#'                     reverse_keyed = grepl("_R$", nm), max = 5L)
#' screen <- cier_screen(bfi_careless[, 1:44], items,
#'                       control = list(cier_personal_reliability = list(seed = 1),
#'                                      cier_gnormed = list(seed = 1)))
#' cier_flagged_cases(screen, min_votes = 3)
cier_flagged_cases <- function(x, ...) {
  UseMethod("cier_flagged_cases")
}

#' @rdname cier_flagged_cases
#' @export
cier_flagged_cases.cier_index <- function(x, ...) {
  # min_votes is cier_screen-only; warn rather than silently absorb a likely mistake.
  if ("min_votes" %in% names(list(...))) {
    cier_warn(
      "cier_warning_ignored_min_votes",
      c("{.arg min_votes} is ignored for a single {.cls cier_index}.",
        "i" = "It applies only to a {.fun cier_screen} (the cross-construct vote \\
               threshold); this index returns every flagged respondent."),
      call = rlang::caller_env()
    )
  }
  # which() drops NA (abstaining) rows.
  unname(which(x$flagged))
}

#' @rdname cier_flagged_cases
#' @export
cier_flagged_cases.cier_screen <- function(x, min_votes = 1L, ...) {
  check_count(min_votes, "min_votes", call = rlang::caller_env())
  # rowSums of collapsed votes = per-respondent construct count. An empty (0-column)
  # screen yields all-zero counts, so no one reaches min_votes -- no special case.
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
