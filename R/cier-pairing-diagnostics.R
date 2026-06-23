# Read-only pairing diagnostics for the synonyms / antonyms indices: which item pairs
# qualify (cier_synonym_pairs) and a critical_r sweep of pairs vs. scorable respondents
# (cier_psychsyn_critval). They reuse the production kernels on the whole sample (the
# index default reference = NULL); neither takes a `reference`, so to mirror a
# cier_psychsyn(reference = ...) run pass that subset/sample as `responses`. One `antonym`
# switch routes the synonym (positive) or antonym (negative) tail.

# Validate the sweep grid: non-empty numeric, every element finite and in (0, 1).
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

# Sorted lower-triangle pair table from the whole-sample pairing correlation. `antonym`
# selects tail and sort: synonyms keep r > critical_r (descending), antonyms keep
# r < -critical_r (ascending). critical_r = NULL lists every non-NA pair; constant-item
# (NA) pairs drop out. Column order matches find_item_pairs (item_i > item_j).
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

# Strongest in-tail inter-item correlation over the lower triangle: max for synonyms,
# min (most negative) for antonyms, ignoring NA; NA when none.
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

# Typed warning when an index call (cier_psychsyn / cier_psychant) discovers no qualifying
# pairs, naming the cause (critical_r too strict for this inventory) and remedy that the
# shared percentile abstention omits. Carries cier_warning_insufficient_items alongside
# its own subclass so cier_screen()'s targeted muffler still silences it; a direct call
# sees the actionable message.
warn_no_pairs <- function(pairing, critical_r, cor_mat,
                          call = rlang::caller_env()) {
  antonym <- identical(pairing, "ant")
  noun <- if (antonym) "antonym" else "synonym"
  strongest <- strongest_pairing_cor(cor_mat, antonym)
  strongest_disp <- round(strongest, 3)
  # Tail-aware sweep hint: the antonym index sweeps the negative tail. Branched because cli
  # does not re-process inline markup inside an interpolated value.
  hint <- if (antonym) {
    "Strongest in-tail inter-item r = {strongest_disp}. Lower {.arg critical_r}, \\
     or sweep candidate thresholds with \\
     {.code cier_psychsyn_critval(antonym = TRUE)}."
  } else {
    "Strongest in-tail inter-item r = {strongest_disp}. Lower {.arg critical_r}, \\
     or sweep candidate thresholds with {.fun cier_psychsyn_critval}."
  }
  cier_warn(
    c("cier_warning_no_pairs", "cier_warning_insufficient_items"),
    c("No {noun} pairs clear {.arg critical_r} = {critical_r}; every \\
       respondent abstains.",
      "i" = hint),
    # `antonym` in the payload lets cier_screen's no-pairs note phrase the noun + sweep
    # hint from the condition itself, without re-deriving the pairing.
    data = list(critical_r = critical_r, strongest_r = strongest,
                antonym = antonym, n_used = 0L),
    call = call
  )
}

# Shared cutoff tail for the two pair-based indices. With qualifying pairs it is the
# ordinary percentile tail; with none it swaps the generic percentile abstention for the
# actionable no-pairs warning above (muffling the redundant generic one).
resolve_pair_index_cutoff <- function(value, row, fpr, cutoff, no_pairs,
                                      pairing, critical_r, cor_mat, call) {
  if (!no_pairs) {
    return(resolve_index_cutoff(value, row, fpr, cutoff, call = call))
  }
  warn_no_pairs(pairing, critical_r, cor_mat, call = call)
  withCallingHandlers(
    resolve_index_cutoff(value, row, fpr, cutoff, call = call),
    cier_warning_insufficient_items = function(w) invokeRestart("muffleWarning")
  )
}

# ---- Clean-reference pairing ------------------------------------------------

# Resolve the opt-in `reference` of cier_psychsyn() / cier_psychant() into the p x p
# inter-item correlation matrix that DISCOVERS the pairs, decoupling discovery from
# per-respondent scoring (which always uses the full `responses`). A clean reference
# restores the pairs on a heavily contaminated sample where whole-sample correlations
# shrink. `reference`:
#   * NULL (default) -> pairing_cor(responses), the self-pairing path.
#   * atomic, no `dim` (logical mask length n, or integer row indices in 1:n) -> subset
#     of the analysis rows.
#   * has `dim` (matrix / data.frame) -> external clean sample (m x p), checked to match
#     the item count.
# Either path needs >= 3 rows. Pairs are estimated on the raw responses (no reverse-keying).
resolve_pairing_cor <- function(responses, reference,
                                call = rlang::caller_env()) {
  if (is.null(reference)) {
    return(pairing_cor(responses))
  }
  ref_resp <- if (is.atomic(reference) && is.null(dim(reference))) {
    mask <- resolve_reference_mask(reference, nrow(responses), call)
    responses[mask, , drop = FALSE]
  } else {
    # An external clean sample. `arg = "reference"` so a malformed sample blames the
    # reference, not `responses`.
    ext <- check_responses(reference, arg = "reference", call = call)
    if (ncol(ext) != ncol(responses)) {
      cier_abort(
        "cier_error_input",
        c("{.arg reference} sample must measure the same items as \\
           {.arg responses}.",
          "x" = "It has {ncol(ext)} item{?s}; {.arg responses} has \\
                 {ncol(responses)}."),
        data = list(arg = "reference", observed = ncol(ext)), call = call
      )
    }
    align_reference_columns(ext, responses, call)
  }
  if (nrow(ref_resp) < 3L) {
    cier_abort(
      "cier_error_input",
      c("{.arg reference} must provide at least 3 rows to estimate item pairs.",
        "x" = "It provides {nrow(ref_resp)} row{?s}."),
      data = list(arg = "reference", observed = nrow(ref_resp)), call = call
    )
  }
  pairing_cor(ref_resp)
}

# Align an external reference sample to the analysis column order. Pairs are discovered
# positionally, so a different column order would silently pair the wrong columns. With
# unique names on both, items must match as a set and are reordered (name-set mismatch is
# a typed error); otherwise alignment is positional.
align_reference_columns <- function(ext, responses, call) {
  rn <- colnames(responses)
  en <- colnames(ext)
  by_name <- !is.null(rn) && !is.null(en) &&
    !anyDuplicated(rn) && !anyDuplicated(en)
  if (!by_name) {
    return(ext)
  }
  if (!setequal(rn, en)) {
    cier_abort(
      "cier_error_input",
      c("{.arg reference} sample must measure the same items as {.arg responses}.",
        "x" = "Its column names name different items than {.arg responses}.",
        "i" = "A same-items sample in a different column order is aligned by \\
               name; here the item sets differ."),
      data = list(arg = "reference"), call = call
    )
  }
  ext[, match(rn, en), drop = FALSE]
}

# Turn a subset-selector `reference` into a logical mask over the analysis rows: accepts a
# logical mask (length n, no NA, >= 1 TRUE) or integer row indices in 1:n; any other
# atomic vector is a typed input error.
resolve_reference_mask <- function(reference, n, call) {
  if (is.logical(reference)) {
    if (length(reference) != n || anyNA(reference) || !any(reference)) {
      cier_abort(
        "cier_error_input",
        c("{.arg reference} as a logical mask must have length {n}, no \\
           {.val NA}, and select at least one row.",
          "x" = "Got a logical of length {length(reference)} selecting \\
                 {sum(reference, na.rm = TRUE)} row{?s}."),
        data = list(arg = "reference"), call = call
      )
    }
    return(reference)
  }
  ok_idx <- is.numeric(reference) &&
    checkmate::test_integerish(reference, lower = 1L, upper = n,
                               any.missing = FALSE, min.len = 1L)
  if (ok_idx) {
    mask <- logical(n)
    mask[as.integer(reference)] <- TRUE
    return(mask)
  }
  cier_abort(
    "cier_error_input",
    c("{.arg reference} must be {.code NULL}, a logical mask of length {n}, \\
       integer row indices in 1:{n}, or a clean reference sample (a matrix or \\
       data.frame measuring the same items).",
      "x" = "Got {.obj_type_friendly {reference}}."),
    data = list(arg = "reference"), call = call
  )
}

#' List the synonym (or antonym) item pairs and their inter-item correlations
#'
#' Surfaces the whole-sample inter-item correlations that [cier_psychsyn()] (and the
#' antonyms index) use to discover item pairs, so you can see which pairs qualify at a
#' candidate `critical_r` -- or the full distribution -- and choose a threshold. On broad
#' inventories the strongest correlations can fall below the default `0.60`, leaving no
#' synonym pairs; this makes that visible.
#'
#' @details
#' Pairs are the lower triangle of the whole-sample correlation matrix, each appearing once
#' with the larger item index first (`item_i > item_j`). With `critical_r` set, the result
#' is exactly the pairs [cier_psychsyn()] discovers from the whole sample at that threshold;
#' `critical_r = NULL` lists every non-missing pair. No `reference` argument, so to inspect a
#' `cier_psychsyn(reference = ...)` run pass that subset/sample as `responses`. A constant
#' (zero-variance) item is omitted; a threshold no pair clears returns a zero-row frame.
#'
#' @param responses Numeric matrix (or data.frame / tibble coerced internally),
#'   one row per respondent, one column per item.
#' @param critical_r `NULL` (default) to list every pair, or a single number in
#'   `(0, 1)` to keep only pairs whose inter-item correlation clears that magnitude.
#' @param antonym Single logical. `FALSE` (default) lists positive synonym pairs
#'   descending by correlation; `TRUE` lists negative antonym pairs ascending (most
#'   negative first).
#'
#' @return A data frame, one row per pair: `item_i`, `item_j` (integer column
#'   indices, `item_i > item_j`), `name_i`, `name_j` (item names, or `NA` when
#'   `responses` has no column names), and `r` (signed inter-item correlation).
#'
#' @references
#' Meade, A. W., & Craig, S. B. (2012). Identifying careless responses in survey
#' data. *Psychological Methods*, 17(3), 437–455. \doi{10.1037/a0028085}
#'
#' @seealso [cier_psychsyn_critval()] for the coverage sweep, [cier_psychsyn()] for the
#'   index.
#' @family pairing diagnostics
#' @export
#' @examples
#' # BFI correlations fall below 0.60, so list pairs qualifying at a lower threshold:
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
#' Evaluates a grid of candidate `critical_r` values, reporting for each how many synonym
#' (or antonym) item pairs qualify and how many respondents become scorable by
#' [cier_psychsyn()]. Lowering `critical_r` admits more pairs, usually scoring more
#' respondents at the cost of weaker pairs -- this table makes that trade legible.
#'
#' @details
#' `critical_r` is a property of the inventory, not a calibrated cutoff: it sets which item
#' pairs count as near-synonyms. The sweep reuses the production pairing and scoring on the
#' whole sample, so coverage is exactly what the index produces under its default
#' `reference = NULL`. No `reference` argument; to weigh coverage for a
#' `cier_psychsyn(reference = ...)` run, pass that subset/sample as `responses`.
#' `strongest_r` is the strongest in-tail inter-item correlation, constant across the grid
#' (most positive for synonyms, most negative for antonyms), so it gives the lowest
#' threshold that can find any pair. When even the strictest grid value finds no pairs, a
#' typed message (`cier_message_no_pairs`) points you to lower `critical_r`.
#'
#' @param responses Numeric matrix (or data.frame / tibble coerced internally),
#'   one row per respondent, one column per item.
#' @param critical_r Numeric vector of candidate thresholds, each in `(0, 1)`.
#'   Defaults to a grid from `0.70` down to `0.30`.
#' @param antonym Single logical. `FALSE` (default) sweeps the synonym (positive)
#'   tail; `TRUE` sweeps the antonym (negative) tail.
#'
#' @return A data frame, one row per `critical_r`: `critical_r`, `n_pairs`,
#'   `n_scored`, `n_abstain`, `abstention_rate` (share of respondents with no usable
#'   score, rounded to 3 dp), and `strongest_r` (strongest in-tail inter-item
#'   correlation, rounded to 3 dp; constant across rows, signed -- most positive for
#'   synonyms, most negative for antonyms).
#'
#' @references
#' Meade, A. W., & Craig, S. B. (2012). Identifying careless responses in survey
#' data. *Psychological Methods*, 17(3), 437–455. \doi{10.1037/a0028085}
#'
#' @seealso [cier_synonym_pairs()] to list the pairs themselves, [cier_psychsyn()] for
#'   the index.
#' @family pairing diagnostics
#' @export
#' @examples
#' # On the BFI the default 0.60 finds no pairs; the sweep shows where coverage appears:
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
  # One p x p pairing correlation serves the whole sweep; the grid only moves the
  # threshold filter over it.
  cor_mat <- pairing_cor(responses)
  strongest <- strongest_pairing_cor(cor_mat, antonym)
  # Discover each threshold's pairs once, reused for the count and the scoring pass.
  pairs_by_cr <- lapply(critical_r, function(cr) {
    find_item_pairs(responses, cr, pairing, cor_mat = cor_mat)
  })
  n_pairs <- vapply(pairs_by_cr, nrow, integer(1L))
  n_scored <- vapply(seq_along(critical_r),
                     function(i) {
                       sum(!is.na(kernel_psychsyn(responses, critical_r[[i]],
                                                  pairing,
                                                  pairs = pairs_by_cr[[i]])))
                     },
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
