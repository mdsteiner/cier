# cier_screen() combiner: wrapper, constructor, and print / as.data.frame methods.
# Orchestration and collapse helpers live in R/screen-combiner.R.

# `methods` is derived from the indices that ran; flag count and rate are derived on
# print, never stored.
new_cier_screen <- function(indices, flags, vote_group, votes, agreement,
                            skipped, notes, n_respondents) {
  ran <- names(indices)
  if (is.null(ran)) {
    ran <- character(0L)
  }
  structure(
    list(indices = indices, flags = flags, vote_group = vote_group,
         votes = votes, agreement = agreement, skipped = skipped,
         notes = notes, methods = ran, n_respondents = n_respondents),
    class = "cier_screen"
  )
}

#' Screen a dataset with the indirect C/IER index battery
#'
#' Runs the package's default detection indices over one response dataset and returns
#' a transparent flag table plus a cross-index agreement diagnostic. It is an
#' orchestrator, not a new statistic: every index it runs is identical to calling that
#' index directly. It produces no single careless/not label -- it reports the count of
#' flagged constructs and the agreement diagnostic, and the researcher thresholds.
#'
#' @details
#' **What runs.** Ten screenable indices run in package order. Six need only
#' `responses` ([cier_longstring()], [cier_irv()], [cier_psychsyn()],
#' [cier_psychant()], [cier_mahalanobis()], [cier_person_total()]); four also need
#' `items` ([cier_even_odd()], [cier_personal_reliability()], [cier_gnormed()],
#' [cier_ht()]). Standalone indices (e.g. [cier_autocorrelation()]) are not run.
#' `$skipped` records, with a reason, each metadata index when `items` is `NULL` and
#' any index hitting a typed backend limit on otherwise-valid data (the single-`Ncat`
#' contract for [cier_gnormed()] on mixed formats), so one index's limit never aborts
#' the battery. A genuinely malformed `items` frame is not skipped -- its typed error
#' surfaces.
#'
#' **Mixed response formats.** Only [cier_gnormed()] requires a single number of
#' response categories across items (Niessen et al., 2016), so the screen skips it on
#' mixed-format data. Other indices run on differing ranges, but published validation
#' (Curran, 2016) is almost entirely on uniform formats, so interpret mixed-format
#' screens with care. Caveat: a respondent who always picks the same option position
#' produces varying raw values across blocks of differing range, so the consistency
#' indices' zero-variance abstention no longer catches such straightliners -- pair them
#' with [cier_longstring()] and [cier_irv()].
#'
#' **Selecting indices.** Goldammer et al. (2024) report resampled personal reliability
#' as the single strongest indirect indicator, so weaker indices can be off-selected
#' (`methods = "cier_personal_reliability"`). Tune any index through `control`: a named
#' list keyed by method id whose entries are argument lists spliced into that index's
#' call -- e.g. `control = list(cier_irv = list(fpr = 0.10), cier_gnormed =
#' list(seed = 1))`. Pass a `seed` for the random indices (RPR, Gnormed's Monte-Carlo
#' null) for reproducibility.
#'
#' **The combiner (transparent, no model).** Each index contributes one per-respondent
#' flag (`$flags`, a respondent x index table, `NA` where the index abstains).
#' Correlated indices then collapse to one vote: those sharing a vote group fuse into a
#' single vote firing when any member flagged the respondent. Only even-odd and personal
#' reliability share a group (`consistency`), avoiding double-counting that construct.
#' Collapsed votes are in `$votes`; flagged constructs per respondent is
#' `rowSums(votes)`.
#'
#' **The vote count under missing data.** An index that abstains for a respondent counts
#' as not flagged, never as a missing vote. The denominator is the number of vote groups
#' that ran, so a respondent who abstained on several indices has fewer chances to reach
#' `k` votes; read a low vote count alongside `$flags` (which shows the `NA` abstentions)
#' and `$skipped`.
#'
#' **Agreement, not a per-index rate.** An empirical-percentile cutoff flags its target
#' rate by construction, so a per-index flag rate is tautological. `$agreement` instead
#' reports, for each level k, the observed share flagged by at least k votes against the
#' share expected if votes fired independently; observed far above expected makes a
#' clustered careless subgroup visible. The `<- excess` marker appears only where the
#' observed count is unlikely under vote independence (one-sided binomial tail below
#' 0.05) -- a descriptive guard, not a formal test. The per-vote table
#' (`$agreement$per_vote`) measures `excess` against a null defined only over scored
#' respondents, so when a null-referenced index (Mahalanobis, Gnormed) abstains it is
#' biased low.
#'
#' **Interpreting and reporting.** The percentile cutoff is a ranking convention, not a
#' calibrated false-positive rate: it cuts a tail of the size you choose (`fpr`) and
#' flags at least that share of this sample. Treat exclusion as a researcher decision:
#' report results before and after removing flagged respondents, and show the flag rate
#' across `fpr` (e.g. 0.01, 0.05, 0.10). `vignette("cier")` walks through this.
#'
#' @param responses Numeric matrix (or data.frame / tibble coerced internally), one row
#'   per respondent and one column per item; `NA` marks a missing response. Pass only
#'   item columns -- an ID, label, or free-text column makes the matrix non-numeric (or
#'   silently corrupts every index).
#' @param items Optional item metadata, one row per item aligned to the columns of
#'   `responses` (`scale`, `reverse_keyed`, `max`, optional `min`; see e.g.
#'   [cier_even_odd()]). When `responses` has column names you may give the metadata an
#'   `item` column (or row names) of matching item ids; a mismatch is a typed error.
#'   `NULL` (default) runs only the six matrix-only indices and skips the four metadata
#'   indices with a recorded reason.
#' @param methods Optional character vector of method ids to run. `NULL` (default) runs
#'   every screenable index. An unknown id is a typed error listing the screenable ids.
#' @param control Optional named list of per-index argument overrides, keyed by method
#'   id; each entry is a list of arguments forwarded to that index (e.g. `fpr`, `alpha`,
#'   `cutoff`, `seed`, `critical_r`, `n_resamples`). Names must be among the selected
#'   methods.
#' @param fpr Optional battery-wide target false-positive tail mass in `(0, 1)`,
#'   forwarded to the percentile-cutoff indices only so you can sweep the flag rate in
#'   one call. `NULL` (default) leaves each index at its own default. The Mahalanobis
#'   chi-square `alpha`, longstring `frac`, and Gnormed Monte-Carlo nominal live on their
#'   own scales and are not moved (tune via `control`). A per-index `control` entry
#'   setting `fpr` or a literal `cutoff` wins.
#'
#' @return A `cier_screen`: a list with `indices` (the per-index `cier_index` objects
#'   that ran), `flags` (respondent x index flag table), `vote_group` (index-to-construct
#'   map), `votes` (respondent x construct collapsed votes), `agreement` (cross-index
#'   agreement diagnostic, or `NULL` when nothing ran), `skipped`
#'   (`data.frame(method, reason)`), `notes` (`data.frame(method, note)` carrying the
#'   reason behind a `0 / 0` line), `methods` (the indices that ran), and
#'   `n_respondents`. Use [as.data.frame()][as.data.frame.cier_screen] for a tidy table
#'   and `print()` for a summary.
#'
#' @references
#' Curran, P. G. (2016). Methods for the detection of carelessly invalid
#' responses in survey data. *Journal of Experimental Social Psychology*, 66,
#' 4-19. \doi{10.1016/j.jesp.2015.07.006}
#'
#' Goldammer, P., Stöckli, P. L., Escher, Y. A., Annen, H., Jonas, K., &
#' Antonakis, J. (2024). Careless responding detection revisited: Accuracy of
#' direct and indirect measures. *Behavior Research Methods*, 56, 8422-8449.
#'
#' Meade, A. W., & Craig, S. B. (2012). Identifying careless responses in survey
#' data. *Psychological Methods*, 17(3), 437-455. \doi{10.1037/a0028085}
#'
#' Niessen, A. S. M., Meijer, R. R., & Tendeiro, J. N. (2016). Detecting
#' careless respondents in web-based questionnaires: Which method to use?
#' *Journal of Research in Personality*, 63, 1-11.
#' \doi{10.1016/j.jrp.2016.04.010}
#'
#' @seealso The index functions [cier_longstring()], [cier_irv()],
#'   [cier_personal_reliability()], [cier_gnormed()], [cier_ht()]
#' @family orchestration
#' @export
#' @examples
#' # First 44 columns are the BFI items, coded 1..5; "_R" marks a reverse-keyed item;
#' # scale is the BFI domain (letters between the "v_BFI_" prefix and the item number).
#' # Seed both randomised pieces (RPR resampling, Gnormed's Monte-Carlo null).
#' nm <- names(bfi_careless)[1:44]
#' items <- data.frame(
#'   scale = sub("^v_BFI_([A-Za-z]+)[0-9].*$", "\\1", nm),
#'   reverse_keyed = grepl("_R$", nm),
#'   max = 5L
#' )
#' screen <- cier_screen(bfi_careless[, 1:44], items,
#'                       control = list(cier_personal_reliability = list(seed = 1),
#'                                      cier_gnormed = list(seed = 1)))
#' screen
#' head(as.data.frame(screen))
cier_screen <- function(responses, items = NULL, methods = NULL,
                        control = list(), fpr = NULL) {
  call <- rlang::caller_env()
  responses <- check_responses(responses, call = call)
  selected <- screen_resolve_methods(methods, call)
  screen_check_control(control, selected, call)
  # Splice a battery-wide fpr into the percentile family; per-index control or a
  # literal cutoff wins (see screen_apply_fpr).
  control <- screen_apply_fpr(control, selected, fpr, call)
  # Emit the unkeyed warning once here (per-index repeats muffled) when any
  # keying-sensitive index runs.
  if (any(selected %in% cier_methods_requiring_items())) {
    inform_if_unkeyed(items, ncol(responses))
  }
  indices <- list()
  skipped_rows <- list()
  notes_rows <- list()
  for (m in selected) {
    reason <- screen_skip_reason(m, items)
    if (is.na(reason)) {
      # screen_call_index returns list(result, note): `result` is the index, or a
      # character skip reason on a typed backend limit (any other error propagates);
      # `note` is the captured no-pairs annotation for the "0 / 0" line.
      out <- screen_call_index(m, responses, items, control[[m]])
      if (!is.character(out$result)) {
        indices[[m]] <- out$result
        if (!is.null(out$note)) {
          notes_rows[[length(notes_rows) + 1L]] <-
            data.frame(method = m, note = out$note, stringsAsFactors = FALSE)
        }
        next
      }
      reason <- out$result
    }
    skipped_rows[[length(skipped_rows) + 1L]] <-
      data.frame(method = m, reason = reason, stringsAsFactors = FALSE)
  }
  build_cier_screen(indices, skipped_rows, notes_rows, control, nrow(responses))
}

# ---- print ------------------------------------------------------------------

# Header: indices that ran, respondent count, and votes remaining after the collapse.
screen_header_line <- function(x) {
  ni <- length(x$indices)
  nv <- ncol(x$votes)
  sprintf(
    "%d %s on %d respondent%s; %d vote%s after collapsing shared constructs.",
    ni, if (ni == 1L) "index" else "indices", x$n_respondents,
    if (x$n_respondents == 1L) "" else "s", nv, if (nv == 1L) "" else "s"
  )
}

# One "method: n_flagged / n_scored (pct%)" line per index, with consistency members
# marked so the collapse is visible.
screen_index_lines <- function(x) {
  if (length(x$indices) == 0L) {
    return(character(0L))
  }
  noted <- x$notes$method
  body <- vapply(names(x$indices), function(m) {
    n_flagged <- sum(x$flags[[m]], na.rm = TRUE)
    n_scored <- sum(!is.na(x$indices[[m]]$value))
    # round() before sprintf keeps the %.1f digit platform-stable (see format_share_pct).
    pct <- if (n_scored > 0L) {
      sprintf("%.1f", round(100 * n_flagged / n_scored, 1))
    } else {
      "--"
    }
    grp <- if (identical(x$vote_group[[m]], "consistency")) "  [consistency]" else ""
    # A "*" marker ties a "0 / 0" pair-index line to its Notes entry below.
    mark <- if (m %in% noted) "  *" else ""
    sprintf("  %-26s %d / %d (%s%%)%s%s", m, n_flagged, n_scored, pct, grp, mark)
  }, character(1L), USE.NAMES = FALSE)
  c("", "Per-index flags (transparent; not independent votes):", body)
}

# One "* method: reason" line per captured annotation (today the psychsyn / psychant
# no-pairs reason), parallel to the skip summary.
screen_notes_lines <- function(notes) {
  if (nrow(notes) == 0L) {
    return(character(0L))
  }
  c("", "Notes:", sprintf("  * %s: %s", notes$method, notes$note))
}

# Pre-round with round(., 1) BEFORE sprintf: a Poisson-binomial `expected` can land one
# ulp from a %.1f boundary where C's sprintf rounding diverges across platforms and
# breaks the snapshot bytewise. A nonzero share too small to survive rounding prints
# "<0.1%", not a misleading "0.0%".
format_share_pct <- function(p) {
  rounded <- round(100 * p, 1)
  if (p > 0 && rounded == 0) {
    "<0.1%"
  } else {
    sprintf("%.1f%%", rounded)
  }
}

# Observed share flagged by >= k votes vs the independence baseline. The "<- excess"
# marker fires only when the one-sided binomial tail P(count >= observed) < 0.05 -- a
# descriptive guard, not a calibrated test.
screen_agreement_lines <- function(agreement, n) {
  if (is.null(agreement)) {
    return(character(0L))
  }
  ag <- agreement$agreement
  body <- vapply(seq_len(nrow(ag)), function(i) {
    obs_n <- round(ag$observed[i] * n)
    tail_p <- stats::pbinom(obs_n - 1L, n, ag$expected[i], lower.tail = FALSE)
    mark <- if (tail_p < 0.05) "  <- excess" else ""
    sprintf("  flagged by >= %d vote%s: %d / %d (%s); expected %s%s",
            ag$k[i], if (ag$k[i] == 1L) "" else "s", obs_n, n,
            format_share_pct(ag$observed[i]), format_share_pct(ag$expected[i]),
            mark)
  }, character(1L))
  c("", "Cross-index agreement (observed vs independence baseline):", body)
}

# Header count and one "method: reason" line per skip.
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
    screen_notes_lines(x$notes),
    screen_agreement_lines(x$agreement, x$n_respondents),
    screen_skipped_lines(x$skipped)
  )
  # Render the rule through cli, then cat the body: cli_verbatim() drops empty strings,
  # so the blank-line section separators go through cat(), which keeps them.
  rule <- cli::cli_format_method(cli::cli_rule(left = "cier_screen"))
  cat(rule, body, sep = "\n")
  cat("\n")
  invisible(x)
}

#' Summarise a cier screen
#'
#' The screen's [print()][print.cier_screen] output is its summary, so `summary()`
#' delegates to it.
#'
#' @param object A `cier_screen` object, as returned by [cier_screen()].
#' @param ... Ignored.
#' @return `object`, invisibly.
#' @export
summary.cier_screen <- function(object, ...) {
  print(object)
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
  out <- do.call(rbind, pieces)
  # Index value vectors can carry the input's row names; reset to plain integers.
  rownames(out) <- NULL
  out
}
