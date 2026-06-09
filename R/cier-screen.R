# Purpose: cier_screen() -- the public transparent flag-table combiner. It runs
#          the registry's screenable indices over one dataset, skips (with a
#          reason) the ones whose inputs are absent, and returns a `cier_screen`:
#          the per-index flag table, the per-construct collapsed votes, and the
#          cross-index agreement diagnostic. The orchestration / collapse helpers
#          live in R/screen-combiner.R; this file holds the wrapper, the
#          constructor, and the print / as.data.frame methods.
# Args:    See documentation below.
# Returns: a `cier_screen` (see new_cier_screen()).
# Invariants:
#   - Each `cier_index` in `$indices` is byte-identical to a direct index call.
#   - even-odd + personal_reliability collapse to ONE vote (one construct); the
#     count of flagged constructs never double-counts them.

# Assemble the cier_screen list. `methods` is derived from the indices that ran;
# the flag count and rate are derived on print, never stored.
new_cier_screen <- function(indices, flags, vote_group, votes, agreement,
                            skipped, n_respondents) {
  ran <- names(indices)
  if (is.null(ran)) {
    ran <- character(0L)
  }
  structure(
    list(indices = indices, flags = flags, vote_group = vote_group,
         votes = votes, agreement = agreement, skipped = skipped,
         methods = ran, n_respondents = n_respondents),
    class = "cier_screen"
  )
}

#' Screen a dataset with the indirect C/IER index battery
#'
#' Runs the package's default detection indices over one response dataset and
#' returns a transparent flag table plus a cross-index **agreement** diagnostic.
#' `cier_screen()` is an orchestrator, not a new statistic: every index it runs
#' is identical to calling that index directly. It deliberately produces **no
#' single careless/not label** -- the count of flagged constructs and the
#' agreement diagnostic are reported, and the researcher chooses how to threshold.
#'
#' @details
#' **What runs.** The registry's `screenable` indices (see `cier_methods()`),
#' optionally restricted by `methods`, are run in registry order. The four
#' indices that need item metadata ([cier_even_odd()], [cier_personal_reliability()],
#' [cier_gnormed()], [cier_ht()]) are **skipped with a recorded reason** when
#' `items` is `NULL`; the two backed by a `Suggests` package ([cier_gnormed()] via
#' `PerFit`, [cier_ht()] via `mokken`) are skipped when that package is not
#' installed. The skipped indices and their reasons are in `$skipped`. A genuinely
#' malformed `items` frame is **not** skipped -- the index's own typed error is
#' left to surface so you can fix the metadata.
#'
#' **Selecting indices.** Goldammer et al. (2024) report resampled personal
#' reliability as the single strongest indirect indicator, so weaker indices can
#' be **off-selected** (`methods = "cier_personal_reliability"`) to keep them from
#' diluting it. Tune any index through `control`: a named list keyed by method id
#' whose entries are argument lists spliced into that index's call -- e.g.
#' `control = list(cier_irv = list(fpr = 0.10), cier_gnormed = list(seed = 1))`.
#' Pass a `seed` for the random indices ([cier_personal_reliability()]'s RPR,
#' [cier_gnormed()]'s Monte-Carlo null) for a reproducible screen.
#'
#' **The combiner (transparent, no model).** Each index contributes one
#' per-respondent flag (`$flags`, a respondent x index table, `NA` where the index
#' abstains). Correlated indices then **collapse to one vote**: indices that share
#' a registry `vote_group` fuse into a single vote that fires when **any** member
#' flagged the respondent (an abstaining member counts as not flagged). Only
#' even-odd and personal reliability share a group (`consistency`) -- they measure
#' one construct, so weighting them as two votes would double-count. The collapsed
#' votes are in `$votes`; the count of flagged constructs per respondent is
#' `rowSums(votes)`.
#'
#' **Agreement, not a per-index rate.** Because an empirical-percentile cutoff
#' flags its target rate by construction, a per-index flag rate is tautological.
#' `$agreement` (from the internal cross-index agreement diagnostic) instead
#' reports, for each level k, the
#' observed share of respondents flagged by at least k votes against the share
#' expected if the votes fired independently; observed far above expected makes a
#' clustered careless subgroup visible. The Mahalanobis chi-square and Gnormed
#' Monte-Carlo votes carry their calibrated null nominal (the rest are
#' percentile, hence `NA`).
#'
#' @param responses A numeric matrix (or a data.frame / tibble coerced
#'   internally) of responses, one row per respondent and one column per item.
#'   `NA` marks a missing response.
#' @param items Optional item metadata, one row per item aligned to the columns
#'   of `responses` (`scale`, `reverse_keyed`, `categories`, optional `min`).
#'   `NULL` (default) runs only the six matrix-only indices and skips the four
#'   metadata indices with a recorded reason.
#' @param methods Optional character vector of method ids to run (a subset of the
#'   screenable registry methods; see `cier_methods()`). `NULL` (default) runs
#'   every screenable index. An unknown id is a typed error.
#' @param control Optional named list of per-index argument overrides, keyed by
#'   method id; each entry is a list of arguments forwarded to that index (for
#'   example `fpr`, `alpha`, `cutoff`, `seed`, `critical_r`, `n_resamples`). Names
#'   must be among the selected methods.
#'
#' @return A `cier_screen`: a list with `indices` (the per-index `cier_index`
#'   objects that ran), `flags` (a respondent x index flag table), `vote_group`
#'   (the index-to-construct map), `votes` (the respondent x construct collapsed
#'   votes), `agreement` (the cross-index agreement diagnostic, or `NULL` when
#'   nothing ran), `skipped` (a `data.frame(method, reason)`), `methods` (the
#'   indices that ran), and `n_respondents`. Use
#'   [as.data.frame()][as.data.frame.cier_screen] for a tidy per-respondent table
#'   and `print()` for a summary.
#'
#' @references
#' Curran, P. G. (2016). Methods for the detection of carelessly invalid
#' responses in survey data. *Journal of Experimental Social Psychology*, 66,
#' 4–19. \doi{10.1016/j.jesp.2015.07.006}
#'
#' Meade, A. W., & Craig, S. B. (2012). Identifying careless responses in survey
#' data. *Psychological Methods*, 17(3), 437–455. \doi{10.1037/a0028085}
#'
#' @seealso The index functions [cier_longstring()], [cier_irv()],
#'   [cier_personal_reliability()], [cier_gnormed()], [cier_ht()]
#' @family orchestration
#' @export
#' @examples
#' # The 44 BFI items are the first 44 columns of the bundled example data; a
#' # trailing "_R" in the column name marks a reverse-keyed item, and the scale
#' # is the BFI domain (the letters between the "v_BFI_" prefix and the item
#' # number, e.g. EX, AG, CON, NEU, OP).
#' nm <- names(bfi_careless)[1:44]
#' items <- data.frame(
#'   scale = gsub("^v_BFI_|[0-9_R]+$", "", nm),
#'   reverse_keyed = grepl("_R$", nm),
#'   categories = 5L
#' )
#' screen <- cier_screen(bfi_careless[, 1:44], items,
#'                       control = list(cier_personal_reliability = list(seed = 1)))
#' screen
#' head(as.data.frame(screen))
cier_screen <- function(responses, items = NULL, methods = NULL,
                        control = list()) {
  call <- rlang::caller_env()
  responses <- check_responses(responses, call = call)
  selected <- screen_resolve_methods(methods, call)
  screen_check_control(control, selected, call)
  reg <- load_method_registry()
  indices <- list()
  skipped_rows <- list()
  for (m in selected) {
    reason <- screen_skip_reason(m, items, reg[reg$method == m, , drop = FALSE])
    if (!is.na(reason)) {
      skipped_rows[[length(skipped_rows) + 1L]] <-
        data.frame(method = m, reason = reason, stringsAsFactors = FALSE)
      next
    }
    indices[[m]] <- screen_call_index(m, responses, items, control[[m]])
  }
  build_cier_screen(indices, skipped_rows, reg, control, nrow(responses))
}

# ---- print ------------------------------------------------------------------

# The header line: how many indices ran, on how many respondents, and how many
# votes remain after the construct collapse.
screen_header_line <- function(x) {
  ni <- length(x$indices)
  nv <- ncol(x$votes)
  sprintf(
    "%d %s on %d respondent%s; %d vote%s after collapsing shared constructs.",
    ni, if (ni == 1L) "index" else "indices", x$n_respondents,
    if (x$n_respondents == 1L) "" else "s", nv, if (nv == 1L) "" else "s"
  )
}

# One "method: n_flagged / n_scored (pct%)" line per index that ran, with the
# consistency members marked so the collapse is visible.
screen_index_lines <- function(x) {
  if (length(x$indices) == 0L) {
    return(character(0L))
  }
  body <- vapply(names(x$indices), function(m) {
    n_flagged <- sum(x$flags[[m]], na.rm = TRUE)
    n_scored <- sum(!is.na(x$indices[[m]]$value))
    pct <- if (n_scored > 0L) sprintf("%.1f", 100 * n_flagged / n_scored) else "--"
    grp <- if (identical(x$vote_group[[m]], "consistency")) "  [consistency]" else ""
    sprintf("  %-26s %d / %d (%s%%)%s", m, n_flagged, n_scored, pct, grp)
  }, character(1L), USE.NAMES = FALSE)
  c("", "Per-index flags (transparent; not independent votes):", body)
}

# The agreement table: observed share flagged by >= k votes vs the independence
# baseline, with the counts derived from the collapsed votes.
screen_agreement_lines <- function(agreement, n) {
  if (is.null(agreement)) {
    return(character(0L))
  }
  ag <- agreement$agreement
  body <- vapply(seq_len(nrow(ag)), function(i) {
    obs_n <- round(ag$observed[i] * n)
    mark <- if (ag$observed[i] > ag$expected[i] + 1e-9) "  <- excess" else ""
    sprintf("  flagged by >= %d vote%s: %d / %d (%.1f%%); expected %.1f%%%s",
            ag$k[i], if (ag$k[i] == 1L) "" else "s", obs_n, n,
            100 * ag$observed[i], 100 * ag$expected[i], mark)
  }, character(1L))
  c("", "Cross-index agreement (observed vs independence baseline):", body)
}

# The skip summary: a header count and one "method: reason" line per skip.
screen_skipped_lines <- function(skipped) {
  head <- sprintf("Skipped: %d", nrow(skipped))
  if (nrow(skipped) == 0L) {
    return(c("", head))
  }
  c("", head, sprintf("  %s: %s", skipped$method, skipped$reason))
}

#' Print a cier screen
#'
#' @param x A `cier_screen` object, as returned by [cier_screen()].
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.cier_screen <- function(x, ...) {
  body <- c(
    screen_header_line(x),
    screen_index_lines(x),
    screen_agreement_lines(x$agreement, x$n_respondents),
    screen_skipped_lines(x$skipped)
  )
  # Render the rule through cli (honours cli.width / cli.unicode), then cat the
  # body directly: cli_verbatim() drops empty-string elements, so the blank-line
  # section separators in `body` are emitted via cat(), which preserves them.
  rule <- cli::cli_format_method(cli::cli_rule(left = "cier_screen"))
  cat(rule, body, sep = "\n")
  cat("\n")
  invisible(x)
}

#' Coerce a cier screen to a data frame
#'
#' @param x A `cier_screen` object.
#' @param ... Ignored.
#' @return A tidy long data frame with one row per respondent x index: columns
#'   `respondent`, `method`, `value`, `flagged`, and `vote_group`.
#' @export
as.data.frame.cier_screen <- function(x, ...) {
  if (length(x$indices) == 0L) {
    return(data.frame(respondent = integer(0L), method = character(0L),
                      value = numeric(0L), flagged = logical(0L),
                      vote_group = character(0L), stringsAsFactors = FALSE))
  }
  n <- x$n_respondents
  pieces <- lapply(names(x$indices), function(m) {
    ix <- x$indices[[m]]
    data.frame(respondent = seq_len(n), method = m, value = ix$value,
               flagged = ix$flagged, vote_group = unname(x$vote_group[[m]]),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, pieces)
}
