# Independent reference for cier_page_time() (timing family; Bowling, Huang,
# Brower & Bragg 2023, <doi:10.1177/10944281211056520>).
#
# page_time is a counting rule, not an estimator: per respondent, count the
# survey pages whose mean per-item time falls strictly below `min_seconds`. The
# mean per-item time on a page is the respondent's total time on that page
# divided by the number of items on it. The per-respondent value is that count;
# a page with no recorded time (NA) contributes no evidence, and a respondent
# with no timed page at all abstains (NA).
#
# These helpers re-derive each quantity by an independent path -- a per-
# respondent row loop, never the production kernel's vectorised masked rowSums --
# so any divergence is attributable to the wrapper. page_time has no CRAN parity
# partner (no package implements page time as a C/IER index; verified
# 2026-06-10), so this oracle is the sole parity check -- oracle-only trust, like
# total_time / PR / RPR. Same base-R arithmetic as the wrapper, so the target
# tolerance is 0 (exact integer counts).

# Per-respondent rapid-page count. For each respondent the per-item time on page
# j is page_seconds[i, j] / items_per_page[j]; a page is rapid when that mean is
# strictly below `min_seconds`. NA pages are dropped (no evidence); a respondent
# with every page NA scores NA (abstains). A mutant that compares the page TOTAL
# (forgets to divide), uses <= instead of <, counts non-rapid pages, or recycles
# items_per_page down rows instead of across columns diverges from this.
ref_page_time <- function(page_seconds, items_per_page, min_seconds = 2) {
  page_seconds <- as.matrix(page_seconds)
  n <- nrow(page_seconds)
  vapply(seq_len(n), function(i) {
    per_item <- page_seconds[i, ] / items_per_page
    observed <- !is.na(per_item)
    if (!any(observed)) {
      return(NA_real_)
    }
    sum(per_item[observed] < min_seconds)
  }, numeric(1L))
}

# Fixed-cutoff resolution. A literal count passes through verbatim (n_pages
# NULL); a fraction of the page count resolves to ceiling(frac * n_pages), with
# the same round-to-9-dp guard the production resolver uses so e.g. 0.28 * 25
# resolves to 7, not 8 (0.28 * 25 == 7.0000000000000009 in IEEE-754).
ref_page_time_fixed_cutoff <- function(value, n_pages = NULL) {
  if (is.null(n_pages)) {
    return(as.numeric(value))
  }
  ceiling(round(value * n_pages, 9L))
}

# Upper-tail flag rule: a respondent is flagged iff their (non-NA) count is at or
# above the cutoff (ties flag). An NA cutoff flags nobody; an NA value (abstain)
# is never flagged.
ref_page_time_flags <- function(value, cutoff) {
  if (is.na(cutoff)) {
    return(rep(FALSE, length(value)))
  }
  !is.na(value) & value >= cutoff
}
