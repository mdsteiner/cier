# Purpose: Shared pure kernels for the timing family. cier_total_time has no
#          kernel (its value is the validated seconds vector itself, an
#          identity), so this file opens with cier_page_time -- the first timing
#          index that actually computes something (the single-kernel rule governs
#          statistics, not identities; see ADR.md "Timing family").
# Args:    See per-kernel signatures.
# Returns: Documented per kernel.
# Invariants:
#   - Kernels are pure (no side effects, no I/O, no registry read).
#   - No statistical computation in user-facing files; this file holds the
#     timing math consumed by R/cier-page-time.R.

# Purpose: Per-respondent count of pages whose MEAN PER-ITEM time falls strictly
#          below `min_seconds` -- the Bowling, Huang, Brower & Bragg (2023)
#          page-time rapid-responding index. The mean per-item time on page j is
#          the respondent's total time on that page divided by the page's item
#          count, so a page total is normalised by `items_per_page[j]` before the
#          threshold comparison.
# Args:    page_seconds   - numeric matrix (n x pages); per-page TOTAL time in
#                           seconds, one column per page. NA marks an untimed
#                           page; all cells are finite and strictly positive
#                           (validated upstream by check_page_seconds()).
#          items_per_page - integer length-(pages); items on each page, all >= 1
#                           (validated upstream by check_items_per_page()).
#          min_seconds    - single non-negative number; the per-item threshold
#                           (default 2 per Bowling et al. 2023).
# Returns: numeric length-n; the count of rapid pages per respondent, NA where
#          every page is untimed (no evidence). value[i] is in [0, pages] when
#          the respondent has at least one timed page.
# Invariants:
#   - The comparison is strictly below (`<`), so a per-item mean exactly at
#     `min_seconds` is NOT rapid.
#   - An NA (untimed) page contributes no evidence: it is neither counted rapid
#     nor counted toward the timed-page total. Because `observed` is FALSE there,
#     `observed & per_item < min_seconds` is FALSE (not NA), so the row sums are
#     plain integer arithmetic with no NA propagation.
#   - A respondent with no timed page (every cell NA) scores NA (abstains).
kernel_page_time <- function(page_seconds, items_per_page, min_seconds) {
  # Divide each column by its item count to get the per-item mean time per page.
  per_item <- sweep(page_seconds, 2L, items_per_page, "/")
  observed <- !is.na(per_item)
  rapid <- observed & per_item < min_seconds
  value <- as.numeric(rowSums(rapid))
  value[rowSums(observed) == 0L] <- NA_real_
  unname(value)
}
