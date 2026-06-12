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
#' **What runs.** Ten indices are screenable, run in registry order. Six need
#' only `responses`: [cier_longstring()], [cier_irv()], [cier_psychsyn()],
#' [cier_psychant()], [cier_mahalanobis()], [cier_person_total()]. Four also
#' need `items` metadata: [cier_even_odd()], [cier_personal_reliability()],
#' [cier_gnormed()] (requires `PerFit`), [cier_ht()] (requires `mokken`).
#' Standalone indices (e.g. [cier_autocorrelation()]) are not run by the screen
#' -- call them directly; [cier_methods()] lists every index with its registry
#' properties. The four metadata indices are **skipped with a recorded reason**
#' when `items` is `NULL`, and the two backend indices when their `Suggests`
#' package is not installed. An index that hits a typed **backend limit** on
#' otherwise-valid data (`mokken`'s 10-category ceiling for [cier_ht()];
#' `PerFit`'s single-`Ncat` contract for [cier_gnormed()] on mixed response
#' formats) is likewise recorded as skipped with the limit as the reason, so
#' one index's ceiling never aborts the battery. The skipped indices and their
#' reasons are in `$skipped`. A genuinely
#' malformed `items` frame is **not** skipped -- the index's own typed error is
#' left to surface so you can fix the metadata.
#'
#' **Mixed response formats.** Only [cier_gnormed()] requires a single number
#' of response categories across all items (`PerFit`'s polytomous statistics
#' score one `Ncat`; see also Niessen et al., 2016) -- on mixed-format data the
#' screen records it as skipped, as described above. Every other index runs on
#' items with differing response ranges, but the published validation of these
#' indices is almost entirely confined to uniform formats: Curran (2016)
#' describes even-odd consistency "across a series of similar scales using the
#' same response format", and `mokken` warns (since its 3.0.3 release) when
#' items have different numbers of response categories. Interpret mixed-format
#' screens with corresponding care. One mechanical caveat: a respondent who
#' always picks the same option *position* produces varying raw values across
#' blocks whose scale range differs, so the zero-variance abstention of the
#' consistency indices no longer applies to such straightliners -- pair them
#' with [cier_longstring()] and [cier_irv()], which still catch that pattern.
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
#' clustered careless subgroup visible. In the printed table the `<- excess`
#' marker appears only where the observed count would be unlikely under vote
#' independence (a one-sided binomial tail below 0.05) -- a descriptive guard so
#' ordinary sampling noise is not advertised as contamination, not a formal
#' test. The Mahalanobis chi-square and Gnormed
#' Monte-Carlo votes carry their calibrated null nominal (the rest are
#' percentile, hence `NA`).
#'
#' **Interpreting and reporting.** The percentile cutoff is a **ranking
#' convention**, not a calibrated false-positive rate: it cuts a tail of the size
#' you choose (`fpr`) and flags about that share of *this* sample by construction.
#' It is **not** Goldammer et al.'s simulated-null Sen95 operating point and makes
#' no claim about the true false-positive rate. Treat the exclusion as a
#' researcher decision: **report results before and after** removing flagged
#' respondents, and show the flag rate across `fpr` (for example 0.01, 0.05, 0.10)
#' so readers can see its sensitivity to a threshold you set. `vignette("cier")`
#' walks through this.
#'
#' @param responses A numeric matrix (or a data.frame / tibble coerced
#'   internally) of responses, one row per respondent and one column per item.
#'   `NA` marks a missing response.
#' @param items Optional item metadata, one row per item aligned to the columns
#'   of `responses` (`scale`, `reverse_keyed`, `max`, optional `min`; see e.g.
#'   [cier_even_odd()] for the column definitions).
#'   `NULL` (default) runs only the six matrix-only indices and skips the four
#'   metadata indices with a recorded reason.
#' @param methods Optional character vector of method ids to run.
#'   `NULL` (default) runs every screenable index. An unknown id is a typed
#'   error; see [cier_methods()] for the available set.
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
#' # The 44 BFI items are the first 44 columns of the bundled example data; a
#' # trailing "_R" in the column name marks a reverse-keyed item, the items are
#' # coded 1..5, and the scale is the BFI domain (the letters between the
#' # "v_BFI_" prefix and the item number, e.g. EX, AG, CON, NEU, OP). Seed BOTH
#' # randomised pieces (RPR's resampling, Gnormed's Monte-Carlo null) for a
#' # reproducible screen.
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
    if (is.na(reason)) {
      # screen_call_index returns a character skip reason instead of an index
      # when the call hit a typed backend limit (e.g. mokken's 10-category
      # ceiling); any other error propagates.
      res <- screen_call_index(m, responses, items, control[[m]])
      if (!is.character(res)) {
        indices[[m]] <- res
        next
      }
      reason <- res
    }
    skipped_rows[[length(skipped_rows) + 1L]] <-
      data.frame(method = m, reason = reason, stringsAsFactors = FALSE)
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
    # round() in R before sprintf: keep the %.1f digit platform-stable on a
    # value that could sit one ulp from a rounding boundary (see
    # screen_agreement_lines()).
    pct <- if (n_scored > 0L) {
      sprintf("%.1f", round(100 * n_flagged / n_scored, 1))
    } else {
      "--"
    }
    grp <- if (identical(x$vote_group[[m]], "consistency")) "  [consistency]" else ""
    sprintf("  %-26s %d / %d (%s%%)%s", m, n_flagged, n_scored, pct, grp)
  }, character(1L), USE.NAMES = FALSE)
  c("", "Per-index flags (transparent; not independent votes):", body)
}

# The agreement table: observed share flagged by >= k votes vs the independence
# baseline, with the counts derived from the collapsed votes. The "<- excess"
# marker is gated on chance, not on a strict point comparison: under
# independence each respondent is flagged by >= k votes with probability
# `expected[k]` (the exact Poisson-binomial tail), so the observed COUNT is
# Binomial(n, expected[k]) and would exceed the expectation about half the time
# on clean data. The marker therefore fires only when the one-sided binomial
# tail P(count >= observed) falls below 0.05 -- a descriptive guard against
# advertising pure sampling noise as contamination, not a formal calibrated
# test (the per-k rows are also not independent of each other).
#
# The displayed percentages are pre-rounded with `round(., 1)` (R's
# platform-independent rounding) BEFORE sprintf: a Poisson-binomial `expected`
# can land one ulp from a `%.1f` boundary (e.g. 0.2500000000000001%), where the
# C library's own sprintf rounding diverges across platforms (Windows -> "0.2",
# glibc -> "0.3") and breaks the print snapshot bytewise. Rounding in R first
# makes the formatted digit identical everywhere.
screen_agreement_lines <- function(agreement, n) {
  if (is.null(agreement)) {
    return(character(0L))
  }
  ag <- agreement$agreement
  body <- vapply(seq_len(nrow(ag)), function(i) {
    obs_n <- round(ag$observed[i] * n)
    tail_p <- stats::pbinom(obs_n - 1L, n, ag$expected[i], lower.tail = FALSE)
    mark <- if (tail_p < 0.05) "  <- excess" else ""
    sprintf("  flagged by >= %d vote%s: %d / %d (%.1f%%); expected %.1f%%%s",
            ag$k[i], if (ag$k[i] == 1L) "" else "s", obs_n, n,
            round(100 * ag$observed[i], 1), round(100 * ag$expected[i], 1), mark)
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
  out <- do.call(rbind, pieces)
  # An index value vector can carry the input's row names (kernels do not strip
  # them), which data.frame() would adopt and rbind() then de-duplicate into
  # mangled labels; a tidy long table always shows the plain integer row index.
  rownames(out) <- NULL
  out
}
