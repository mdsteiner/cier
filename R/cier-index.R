# Purpose: The shared light `cier_index` object returned by every index wrapper,
#          plus its print / as.data.frame methods. `cier_index` is a list-based
#          S3 record (value / flagged vectors + method / cutoff / direction
#          metadata) -- the canonical S3 shape, so the metadata cannot be
#          silently dropped or desynced by subsetting (unlike attributes on a
#          data.frame). The flag count and rate are DERIVED on print from
#          `flagged`, never stored, so they cannot go stale.
# Args:    See per-function documentation.
# Returns: new_cier_index() -> a cier_index; print -> x invisibly;
#          as.data.frame -> a data.frame(value, flagged).
# Invariants:
#   - A respondent who abstains has `value = NA`; `flagged` is forced to NA
#     wherever `value` is NA (one rule for every index).
#   - print() output goes to stdout via cli::cli_format_method() + cat().

# Purpose: Assemble a cier_index from a per-respondent value + flagged vector and
#          the resolved cutoff metadata. The single constructor for all indices.
# Args:    value     - numeric per-respondent index value (NA where abstaining).
#          flagged   - logical per-respondent flag from apply_flag().
#          method    - the index id (e.g. "cier_longstring").
#          cutoff    - the resolved numeric cutoff (NA when unresolved).
#          direction - "upper" or "lower".
# Returns: a `cier_index` (list of value, flagged, method, cutoff, direction).
new_cier_index <- function(value, flagged, method, cutoff, direction) {
  flagged[is.na(value)] <- NA               # one abstention rule for every index
  structure(
    list(value = value, flagged = flagged, method = method,
         cutoff = cutoff, direction = direction),
    class = "cier_index"
  )
}

#' Coerce a cier index to a data frame
#'
#' @param x A `cier_index` object.
#' @param ... Ignored.
#' @return A data frame with one `value` and one `flagged` column per respondent.
#' @export
as.data.frame.cier_index <- function(x, ...) {
  data.frame(value = x$value, flagged = x$flagged)
}

#' Print a cier index
#'
#' The flag count and rate are computed over respondents who produced a usable
#' score; abstaining respondents (all-missing rows, `value = NA`) are excluded
#' from both the count and the rate and reported on their own line.
#'
#' @param x A `cier_index` object, as returned by an index function such as
#'   [cier_longstring()].
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.cier_index <- function(x, ...) {
  n         <- length(x$value)
  n_abstain <- sum(is.na(x$value))
  n_scored  <- n - n_abstain
  n_flagged <- sum(x$flagged, na.rm = TRUE)
  dir_word   <- if (identical(x$direction, "upper")) "higher" else "lower"
  comparator <- if (identical(x$direction, "upper")) ">=" else "<="
  pct <- if (n_scored > 0L) sprintf("%.1f", 100 * n_flagged / n_scored) else "--"
  # Round for display only (percentile cutoffs are arbitrary-precision quantiles;
  # integer fixed cutoffs are unchanged: signif(5, 3) == 5). The "~" marks a
  # rounded display so the printed threshold never claims false precision against
  # the exact cutoff that actually drives the flags.
  cutoff <- signif(x$cutoff, 3)
  mark <- if (!is.na(x$cutoff) && !isTRUE(all.equal(cutoff, x$cutoff))) "~" else ""
  lines <- cli::cli_format_method({
    cli::cli_rule(left = x$method)
    cli::cli_text("Direction: {x$direction} -- {dir_word} values flag carelessness.")
    if (is.na(cutoff)) {
      cli::cli_text("Cutoff: none (unresolved).")
    } else {
      cli::cli_text(
        "Cutoff: {mark}{cutoff} -- respondents with value {comparator} {mark}{cutoff} are flagged."
      )
    }
    cli::cli_text("Flagged: {n_flagged} of {n_scored} scored respondent{?s} ({pct}%).")
    if (n_abstain > 0L) {
      cli::cli_text("Abstained: {n_abstain} (no responses).")
    }
    cli::cli_alert_info("Per-respondent scores in {.code $value}, flags in {.code $flagged}.")
  })
  cat(lines, sep = "\n")
  cat("\n")
  invisible(x)
}
