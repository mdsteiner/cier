# Purpose: Read-only pairing diagnostics for the psychometric synonyms / antonyms
#          indices -- cier_synonym_pairs() (which item pairs qualify, and their
#          inter-item correlations) and cier_psychsyn_critval() (a critical_r sweep
#          showing how many pairs qualify and how many respondents become scorable
#          at each candidate threshold). They help a user choose critical_r, which
#          is a property of the inventory, not a calibrated cutoff.
# Args:    See per-function documentation below.
# Returns: Plain data.frames (no new S3 class).
# Invariants:
#   - Both reuse the production pairing kernels verbatim (pairing_cor /
#     find_item_pairs / kernel_psychsyn) and derive no new statistic, so the
#     reported pairs and counts are exactly what the index uses (single-kernel
#     rule). One `antonym` switch routes the synonym (positive) or antonym
#     (negative) tail, serving both cier_psychsyn and cier_psychant.

# Validate the critical_r grid for the sweep: a non-empty numeric vector, every
# element a finite number in the open interval (0, 1).
check_critical_r_grid <- function(grid, call = rlang::caller_env()) {
  ok <- is.numeric(grid) && length(grid) >= 1L && all(is.finite(grid)) &&
    all(grid > 0 & grid < 1)
  if (!ok) {
    cier_abort(
      "cier_error_input",
      c("{.arg critical_r} must be a non-empty numeric vector of values in (0, 1).",
        "i" = "Each value is a candidate pairing threshold to evaluate."),
      data = list(arg = "critical_r"), call = call
    )
  }
  invisible(grid)
}

# Build the sorted lower-triangle pair table from the whole-sample pairing
# correlation. `antonym` selects the tail and the sort: synonyms keep r >
# critical_r and sort descending (strongest positive first); antonyms keep
# r < -critical_r and sort ascending (most negative first). `critical_r = NULL`
# lists every non-NA lower-triangle pair. Pairs touching a zero-variance (constant)
# item have an NA correlation and drop out. Column order matches find_item_pairs
# (item_i > item_j), so the filtered set is exactly the index's pairing.
build_pair_table <- function(responses, critical_r, antonym) {
  cor_mat <- pairing_cor(responses)
  low <- lower.tri(cor_mat)
  idx <- which(low, arr.ind = TRUE)
  r <- cor_mat[low]
  keep <- !is.na(r)
  if (!is.null(critical_r)) {
    keep <- keep & (if (antonym) r < -critical_r else r > critical_r)
  }
  item_i <- as.integer(idx[keep, 1L])
  item_j <- as.integer(idx[keep, 2L])
  r <- r[keep]
  ord <- order(r, decreasing = !antonym)
  item_i <- item_i[ord]
  item_j <- item_j[ord]
  nm <- colnames(responses)
  data.frame(
    item_i = item_i,
    item_j = item_j,
    name_i = if (is.null(nm)) rep(NA_character_, length(item_i)) else nm[item_i],
    name_j = if (is.null(nm)) rep(NA_character_, length(item_j)) else nm[item_j],
    r = r[ord],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

# The strongest in-tail inter-item correlation over the lower triangle: the
# maximum for synonyms, the minimum (most negative) for antonyms, ignoring the NA
# correlations of constant items; NA when no finite correlation exists.
strongest_pairing_cor <- function(cor_mat, antonym) {
  low <- cor_mat[lower.tri(cor_mat)]
  low <- low[!is.na(low)]
  if (length(low) == 0L) {
    return(NA_real_)
  }
  if (antonym) min(low) else max(low)
}

# Guidance message when the strictest threshold tested still finds no pairs.
inform_no_pairs <- function(strictest, strongest, antonym) {
  noun <- if (antonym) "antonyms" else "synonyms"
  strongest_disp <- round(strongest, 3)
  cier_inform(
    "cier_message_no_pairs",
    c("No {noun} clear {.arg critical_r} = {strictest} \\
       (strongest inter-item r = {strongest_disp}).",
      "i" = "Lower {.arg critical_r} to surface pairs and score more respondents."),
    data = list(critical_r = strictest, strongest_r = strongest)
  )
}

#' List the synonym (or antonym) item pairs and their inter-item correlations
#'
#' Surfaces the whole-sample inter-item correlations that [cier_psychsyn()] (and
#' the antonyms index) use to discover item pairs, so you can see which pairs
#' qualify at a candidate `critical_r` -- or the full distribution -- and choose a
#' threshold. On broad inventories the strongest inter-item correlations can fall
#' below the default `0.60`, leaving no synonym pairs; this diagnostic makes that
#' visible.
#'
#' @details
#' The pairs are the lower triangle of the whole-sample correlation matrix
#' (`stats::cor(responses, use = "pairwise.complete.obs")`), so each pair appears
#' once with the larger item index first (`item_i > item_j`). With `critical_r`
#' set, the returned set is exactly the pairs [cier_psychsyn()] scores at that
#' threshold (synonyms: `r > critical_r`; antonyms, `antonym = TRUE`:
#' `r < -critical_r`). With `critical_r = NULL` every non-missing pair is listed.
#' A pair that involves a constant (zero-variance) item has an undefined
#' correlation and is omitted. A threshold that no pair clears returns a zero-row
#' frame (not an error).
#'
#' @param responses A numeric matrix (or a data.frame / tibble coerced
#'   internally) of responses, one row per respondent and one column per item.
#' @param critical_r Either `NULL` (default) to list every pair, or a single
#'   number in the open interval `(0, 1)` to keep only the pairs whose inter-item
#'   correlation clears that magnitude.
#' @param antonym Single logical. `FALSE` (default) lists positively correlated
#'   synonym pairs, sorted by descending correlation; `TRUE` lists negatively
#'   correlated antonym pairs, sorted by ascending (most negative) correlation.
#'
#' @return A data frame with one row per pair and columns `item_i`, `item_j`
#'   (integer column indices, `item_i > item_j`), `name_i`, `name_j` (the item
#'   names, or `NA` when `responses` has no column names), and `r` (the signed
#'   inter-item correlation).
#'
#' @references
#' Meade, A. W., & Craig, S. B. (2012). Identifying careless responses in survey
#' data. *Psychological Methods*, 17(3), 437–455. \doi{10.1037/a0028085}
#'
#' @seealso [cier_psychsyn_critval()] for the coverage sweep, [cier_psychsyn()]
#'   for the index, and [careless::psychsyn_critval()].
#' @family pairing diagnostics
#' @export
#' @examples
#' # The BFI's strongest inter-item correlations fall below 0.60, so list the
#' # pairs that qualify at a lower threshold:
#' head(cier_synonym_pairs(bfi_careless[, 1:44], critical_r = 0.5))
cier_synonym_pairs <- function(responses, critical_r = NULL, antonym = FALSE) {
  call <- rlang::caller_env()
  responses <- check_responses(responses, call = call)
  check_flag(antonym, "antonym", call = call)
  if (!is.null(critical_r)) {
    check_open_unit(critical_r, "critical_r", call = call)
  }
  build_pair_table(responses, critical_r, antonym)
}

#' Sweep `critical_r` to weigh synonym pairs against respondent coverage
#'
#' Evaluates a grid of candidate `critical_r` values and reports, for each, how
#' many synonym (or antonym) item pairs qualify and how many respondents become
#' scorable by [cier_psychsyn()]. Lowering `critical_r` admits more pairs, which
#' usually scores more respondents (more reliable person estimates) at the cost of
#' weaker pairs -- this table makes that trade legible in one call, and directly
#' answers "how many more datapoints if I lower the cutoff?".
#'
#' @details
#' `critical_r` is a property of the **inventory**, not a calibrated cutoff: it
#' sets which item pairs count as near-synonyms. The sweep reuses the production
#' pairing and scoring (`find_item_pairs()` for `n_pairs`, [cier_psychsyn()]'s
#' scorer for `n_scored`), so the coverage shown is exactly what the index would
#' produce. `strongest_r` -- the strongest in-tail inter-item correlation -- is
#' constant across the grid and tells you the lowest threshold that can find any
#' pair. When even the strictest grid value finds no pairs, a typed message
#' (`cier_message_no_pairs`) points you to lower `critical_r`.
#'
#' @param responses A numeric matrix (or a data.frame / tibble coerced
#'   internally) of responses, one row per respondent and one column per item.
#' @param critical_r A numeric vector of candidate thresholds, each in the open
#'   interval `(0, 1)`. Defaults to a grid from `0.70` down to `0.30`.
#' @param antonym Single logical. `FALSE` (default) sweeps the synonym (positive)
#'   tail; `TRUE` sweeps the antonym (negative) tail.
#'
#' @return A data frame with one row per `critical_r` value and columns
#'   `critical_r`, `n_pairs`, `n_scored`, `n_abstain`, `abstention_rate` (the
#'   share of respondents with no usable score, rounded to 3 dp), and
#'   `strongest_r` (the strongest in-tail inter-item correlation, rounded to 3 dp;
#'   constant across rows).
#'
#' @references
#' Meade, A. W., & Craig, S. B. (2012). Identifying careless responses in survey
#' data. *Psychological Methods*, 17(3), 437–455. \doi{10.1037/a0028085}
#'
#' @seealso [cier_synonym_pairs()] to list the pairs themselves, [cier_psychsyn()]
#'   for the index, and [careless::psychsyn_critval()].
#' @family pairing diagnostics
#' @export
#' @examples
#' # On the BFI the default 0.60 finds no pairs; the sweep shows where coverage
#' # appears as the threshold drops:
#' cier_psychsyn_critval(bfi_careless[, 1:44])
cier_psychsyn_critval <- function(responses,
                                  critical_r = c(0.70, 0.65, 0.60, 0.55, 0.50,
                                                 0.45, 0.40, 0.35, 0.30),
                                  antonym = FALSE) {
  call <- rlang::caller_env()
  responses <- check_responses(responses, call = call)
  check_flag(antonym, "antonym", call = call)
  check_critical_r_grid(critical_r, call = call)
  pairing <- if (antonym) "ant" else "syn"
  n <- nrow(responses)
  strongest <- strongest_pairing_cor(pairing_cor(responses), antonym)
  n_pairs <- vapply(critical_r,
                    function(cr) nrow(find_item_pairs(responses, cr, pairing)),
                    integer(1L))
  n_scored <- vapply(critical_r,
                     function(cr) sum(!is.na(kernel_psychsyn(responses, cr, pairing))),
                     integer(1L))
  out <- data.frame(
    critical_r = critical_r,
    n_pairs = n_pairs,
    n_scored = n_scored,
    n_abstain = n - n_scored,
    abstention_rate = round((n - n_scored) / n, 3),
    strongest_r = round(strongest, 3),
    row.names = NULL
  )
  if (n_pairs[[which.max(critical_r)]] == 0L) {
    inform_no_pairs(max(critical_r), strongest, antonym)
  }
  out
}
