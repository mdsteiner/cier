# The shared `cier_index` S3 record (value/flagged vectors + cutoff metadata) and its
# print / as.data.frame methods. Flag count and rate are derived on print, never stored.

# Assemble a cier_index. `cutoff_method` is the resolution strategy (a
# cier_index_cutoff_provenance() label, or NA); `cutoff_rate` is its rate (fpr / alpha /
# nominal / fraction), NA for rate-free strategies.
new_cier_index <- function(value, flagged, method, cutoff, direction,
                           cutoff_method = NA_character_,
                           cutoff_rate = NA_real_) {
  validate_new_cier_index(value, flagged, method, cutoff, direction,
                          cutoff_method, cutoff_rate)
  flagged[is.na(value)] <- NA               # abstention: NA value => NA flag
  structure(
    list(value = value, flagged = flagged, method = method,
         cutoff = cutoff, direction = direction,
         cutoff_method = cutoff_method, cutoff_rate = cutoff_rate),
    class = "cier_index"
  )
}

# Contract guard for new_cier_index(): `flagged` must be logical and length(value),
# checked BEFORE the abstention assignment, which would otherwise grow a short `flagged`.
validate_new_cier_index <- function(value, flagged, method, cutoff, direction,
                                    cutoff_method = NA_character_,
                                    cutoff_rate = NA_real_) {
  ok <- all(c(
    is.numeric(value),
    is.logical(flagged),
    length(flagged) == length(value),
    rlang::is_string(method),
    is.numeric(cutoff) && length(cutoff) == 1L,
    identical(direction, "upper") || identical(direction, "lower"),
    # cutoff_method: NA or a vocabulary label; cutoff_rate: length-1 numeric.
    (length(cutoff_method) == 1L && is.na(cutoff_method)) ||
      (rlang::is_string(cutoff_method) &&
         cutoff_method %in% cier_index_cutoff_provenance()),
    is.numeric(cutoff_rate) && length(cutoff_rate) == 1L
  ))
  if (!ok) {
    cier_abort(
      "cier_error_state",
      c("Invalid {.cls cier_index} construction (internal contract).",
        "i" = "Expected numeric {.field value}, a logical {.field flagged} of the \\
               same length, a scalar {.field method} / numeric {.field cutoff}, \\
               and {.field direction} in {.val {c(\"upper\", \"lower\")}}.",
        "i" = "The {.field cutoff_method} must be NA or one of \\
               {.val {cier_index_cutoff_provenance()}}, and {.field cutoff_rate} \\
               a scalar numeric."),
      data = list(cutoff_method = cutoff_method, cutoff_rate = cutoff_rate)
    )
  }
  invisible(NULL)
}

#' Coerce a cier index to a data frame
#'
#' @param x A `cier_index` object.
#' @param ... Ignored.
#' @return A data frame with one `value` and one `flagged` column per respondent.
#' @export
as.data.frame.cier_index <- function(x, ...) {
  out <- data.frame(value = x$value, flagged = x$flagged)
  # A named value vector would become row labels; force the plain row index.
  rownames(out) <- NULL
  out
}

# Cutoff rounded for display (signif 3) plus a "~" mark when rounding dropped precision.
# Shared by print and summary. NA cutoff -> NA, no mark.
cutoff_display <- function(cutoff) {
  value <- signif(cutoff, 3)
  mark <- if (!is.na(cutoff) && !isTRUE(all.equal(value, cutoff))) "~" else ""
  list(value = value, mark = mark)
}

# Whether the object carries cutoff provenance (a non-NA vocabulary label); print and
# summary omit it when NA.
has_provenance <- function(x) {
  length(x$cutoff_method) == 1L && !is.na(x$cutoff_method)
}

# English ordinal for a small positive integer (1st, 2nd, ...), with the 11/12/13
# exception. Phrases a percentile cutoff's position.
ordinal <- function(n) {
  n <- as.integer(round(n))
  suffix <- if (n %% 100L %in% 11:13) {
    "th"
  } else {
    switch(as.character(n %% 10L),
           "1" = "st", "2" = "nd", "3" = "rd", "th")
  }
  paste0(n, suffix)
}

# Render the cutoff provenance into a human phrase, split into a `position` clause and an
# optional `rate_phrase` so print (separate line) and summary (inline) punctuate it
# differently. The round() before integer rounding neutralises IEEE noise on a clean
# fraction.
cutoff_descriptor <- function(method, rate, direction) {
  none <- function(position) list(position = position, rate_phrase = "")
  switch(
    method,
    percentile = {
      tail <- if (identical(direction, "upper")) 1 - rate else rate
      q <- round(round(100 * tail, 9L))
      list(position = paste0(ordinal(q), " sample percentile"),
           rate_phrase = paste0("fpr = ", format(rate)))
    },
    chisq = list(position = "chi-square tail",
                 rate_phrase = paste0("alpha = ", format(rate))),
    mc_null = list(position = "Monte-Carlo null",
                   rate_phrase = paste0("nominal = ", format(rate))),
    fixed_fraction = none(paste0("fixed fraction ", format(rate),
                                 " of the item count")),
    median_relative = none(paste0(format(rate), " x the sample median")),
    fixed_count = none("fixed count"),
    kneedle = none("Kneedle elbow (parameter-free)"),
    literal = none("user-supplied threshold"),
    # switch must stay exhaustive over cier_index_cutoff_provenance(); fail loudly.
    cier_abort(
      "cier_error_state",
      "No display descriptor for cutoff_method {.val {method}}.",
      data = list(arg = "cutoff_method", observed = method)
    )
  )
}

# Provenance as one plain phrase for print's "Cutoff method:" line; NULL when no
# provenance (cutoff_method = NA).
cutoff_method_phrase <- function(x) {
  if (!has_provenance(x)) {
    return(NULL)
  }
  d <- cutoff_descriptor(x$cutoff_method, x$cutoff_rate, x$direction)
  if (nzchar(d$rate_phrase)) paste0(d$position, " (", d$rate_phrase, ")") else d$position
}

# Cutoff rendered with provenance inline in parentheses, for summary's single Cutoff
# line. `disp` is the already-rounded display value with its "~" mark.
cutoff_inline_phrase <- function(x, disp) {
  base <- paste0("Cutoff: ", disp)
  if (!has_provenance(x)) {
    return(base)
  }
  d <- cutoff_descriptor(x$cutoff_method, x$cutoff_rate, x$direction)
  inner <- if (nzchar(d$rate_phrase)) {
    paste0(d$position, "; ", d$rate_phrase)
  } else {
    d$position
  }
  paste0(base, " (", inner, ")")
}

#' Print a cier index
#'
#' Flag count and rate are computed over scored respondents; abstaining respondents
#' (`value = NA`) are excluded and reported separately.
#'
#' @param x A `cier_index` object, as returned by e.g. [cier_longstring()].
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
  # round() before sprintf keeps %.1f platform-stable (mirrors cier_screen()).
  pct <- if (n_scored > 0L) {
    sprintf("%.1f", round(100 * n_flagged / n_scored, 1))
  } else {
    "--"
  }
  disp <- cutoff_display(x$cutoff)        # signif(3) value + "~" mark
  cutoff <- disp$value
  mark <- disp$mark
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
    method_phrase <- cutoff_method_phrase(x)         # NULL when no provenance
    if (!is.null(method_phrase)) {
      cli::cli_text("Cutoff method: {method_phrase}.")
    }
    cli::cli_text("Flagged: {n_flagged} of {n_scored} scored respondent{?s} ({pct}%).")
    if (n_abstain > 0L) {
      cli::cli_text("Abstained: {n_abstain} (no score).")
    }
    cli::cli_alert_info("Per-respondent scores in {.code $value}, flags in {.code $flagged}.")
  })
  cat(lines, sep = "\n")
  cat("\n")
  invisible(x)
}

#' Summarise a cier index
#'
#' A fuller view than [print()][print.cier_index]: the header, scored / abstained counts,
#' a five-number summary of the scored values, and the share on the flagged side (equal
#' to the flag rate). Abstaining respondents (`value = NA`) are excluded.
#'
#' @param object A `cier_index` object, as returned by e.g. [cier_irv()].
#' @param ... Ignored.
#' @return `object`, invisibly.
#' @export
summary.cier_index <- function(object, ...) {
  x <- object
  scored    <- x$value[!is.na(x$value)]
  n         <- length(x$value)
  n_scored  <- length(scored)
  n_abstain <- n - n_scored
  n_flagged <- sum(x$flagged, na.rm = TRUE)
  dir_word  <- if (identical(x$direction, "upper")) "higher" else "lower"
  side      <- if (identical(x$direction, "upper")) "above" else "below"
  disp      <- cutoff_display(x$cutoff)   # signif(3) value + "~" mark
  cutoff    <- disp$value
  mark      <- disp$mark
  lines <- cli::cli_format_method({
    cli::cli_rule(left = paste0(x$method, " (summary)"))
    cli::cli_text("Direction: {x$direction} -- {dir_word} values flag carelessness.")
    if (is.na(cutoff)) {
      cli::cli_text("Cutoff: none (unresolved).")
    } else {
      cli::cli_text(paste0(cutoff_inline_phrase(x, paste0(mark, cutoff)), "."))
    }
    cli::cli_text("Scored {n_scored} of {n} respondent{?s} ({n_abstain} abstained).")
    if (n_scored > 0L) {
      qs <- signif(stats::quantile(scored, c(0, .25, .5, .75, 1),
                                   names = FALSE, type = 7L), 3)
      cli::cli_text("Score quartiles: min {qs[[1L]]} | Q1 {qs[[2L]]} | \\
                     median {qs[[3L]]} | Q3 {qs[[4L]]} | max {qs[[5L]]}.")
      if (!is.na(cutoff)) {
        # round() before sprintf keeps %.1f platform-stable (mirrors cier_screen()).
        pct <- sprintf("%.1f", round(100 * n_flagged / n_scored, 1))
        cli::cli_text("At or {side} the cutoff: {n_flagged} of {n_scored} \\
                       scored ({pct}%).")
      }
    }
    cli::cli_alert_info("Per-respondent scores in {.code $value}, flags in {.code $flagged}.")
  })
  cat(lines, sep = "\n")
  cat("\n")
  invisible(x)
}
