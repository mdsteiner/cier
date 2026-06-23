# The cier_sim object from cier_simulate(): a classed list whose slots feed the wrappers
# (cier_screen(x$responses, x$items), cier_total_time(x$seconds), ...) scored against
# x$truth, plus its validator and print method. Generator metadata (GRM structure, timing,
# failure probabilities, seed) rides along as the "generator" attribute.

# The single cier_sim constructor.
new_cier_sim <- function(responses, items, seconds, page_seconds,
                         items_per_page, checks, pass, truth, generator) {
  structure(
    list(responses = responses, items = items, seconds = seconds,
         page_seconds = page_seconds, items_per_page = items_per_page,
         checks = checks, pass = pass, truth = truth),
    generator = generator,
    class = "cier_sim"
  )
}

# Shared abort for the cier_sim validator helpers below.
sim_invalid <- function(what, call) {
  cier_abort(
    "cier_error_input",
    c("Invalid {.cls cier_sim} object: {what}.",
      "i" = "cier_simulate() builds this object; do not edit slots by hand."),
    data = list(arg = "x"), call = call
  )
}

# The slot set and the response / items / seconds shapes. Returns c(n, p).
validate_cier_sim_base <- function(x, call) {
  slots <- c("responses", "items", "seconds", "page_seconds",
             "items_per_page", "checks", "pass", "truth")
  if (!inherits(x, "cier_sim") || !is.list(x) || !identical(names(x), slots)) {
    sim_invalid("the slot set changed", call)
  }
  if (!is.matrix(x$responses) || !is.integer(x$responses)) {
    sim_invalid("$responses must be an integer matrix", call)
  }
  n <- nrow(x$responses)
  p <- ncol(x$responses)
  if (NROW(x$items) != p) {
    sim_invalid("$items must describe one row per item", call)
  }
  # Finite, not merely numeric: cier_total_time() rejects Inf / NaN, so catch it here.
  if (!is.numeric(x$seconds) || length(x$seconds) != n ||
        any(!is.finite(x$seconds))) {
    sim_invalid("$seconds must hold one finite total per respondent", call)
  }
  c(n, p)
}

# Page slots: items_per_page sums to p; page_seconds is a finite n x n_pages matrix
# (finite, not merely numeric -- cier_page_time() rejects Inf / NaN).
validate_cier_sim_pages <- function(x, n, p, call) {
  if (!is.integer(x$items_per_page) || sum(x$items_per_page) != p) {
    sim_invalid("$items_per_page must sum to the item count", call)
  }
  if (!is.matrix(x$page_seconds) ||
        !identical(dim(x$page_seconds), c(n, length(x$items_per_page))) ||
        any(!is.finite(x$page_seconds))) {
    sim_invalid("$page_seconds must be a finite respondents x pages matrix", call)
  }
  invisible(NULL)
}

# Attention-check slots: $checks and $pass travel together; when present their shapes
# agree (n rows, one pass-set per check column).
validate_cier_sim_checks <- function(x, n, call) {
  if (is.null(x$checks) != is.null(x$pass)) {
    sim_invalid("$checks and $pass must travel together", call)
  }
  if (!is.null(x$checks) &&
        (!is.matrix(x$checks) || nrow(x$checks) != n ||
           !is.list(x$pass) || length(x$pass) != ncol(x$checks))) {
    sim_invalid("$checks / $pass shapes disagree", call)
  }
  invisible(NULL)
}

# The truth schema and column types.
validate_cier_sim_truth <- function(truth, n, call) {
  truth_cols <- c("careless", "pattern", "extent", "onset_item",
                  "offset_item", "speeded", "params")
  if (!is.data.frame(truth) || nrow(truth) != n ||
        !identical(names(truth), truth_cols)) {
    sim_invalid("the $truth schema changed", call)
  }
  if (!is.logical(truth$careless) || !is.logical(truth$speeded) ||
        !is.list(truth$params)) {
    sim_invalid("the $truth column types changed", call)
  }
  invisible(NULL)
}

# Cross-check a cier_sim's slots (shapes, truth schema, checks/pass pairing). Raises
# cier_error_input on a violation; returns the object invisibly for validate-and-return.
validate_cier_sim <- function(x, call = rlang::caller_env()) {
  np <- validate_cier_sim_base(x, call)
  validate_cier_sim_pages(x, np[[1L]], np[[2L]], call)
  validate_cier_sim_checks(x, np[[1L]], call)
  validate_cier_sim_truth(x$truth, np[[1L]], call)
  invisible(x)
}

# Per-pattern counts of careless rows, label then count, in radix order so the printed
# order is locale-independent.
sim_print_patterns <- function(truth) {
  counts <- table(truth$pattern[truth$careless])
  counts <- counts[order(names(counts), method = "radix")]
  paste(names(counts), as.integer(counts), collapse = ", ")
}

# Extent counts in the fixed full / partial / temporary order, present only.
sim_print_extents <- function(truth) {
  counts <- table(truth$extent[truth$careless])
  levels_order <- c("full", "partial", "temporary")
  counts <- counts[levels_order[levels_order %in% names(counts)]]
  paste(as.integer(counts), names(counts), collapse = ", ")
}

#' Print a simulated C/IER dataset
#'
#' Summarizes a [cier_simulate()] object: planted careless share, per-pattern and
#' per-extent counts (from `$truth`), the timing and attention-check slots.
#'
#' @param x A `cier_sim` object from [cier_simulate()].
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.cier_sim <- function(x, ...) {
  n <- nrow(x$responses)
  p <- ncol(x$responses)
  n_scales <- length(unique(x$items$scale))
  n_careless <- sum(x$truth$careless)
  # round() before sprintf keeps the %.1f digit platform-stable (mirrors cier_screen()).
  pct <- sprintf("%.1f", round(if (n > 0L) 100 * n_careless / n else 0, 1))
  pattern_line <- if (n_careless > 0L) sim_print_patterns(x$truth) else "none"
  n_pages <- length(x$items_per_page)
  pages_word <- paste(x$items_per_page, collapse = "+")
  n_checks <- if (is.null(x$checks)) 0L else ncol(x$checks)
  lines <- cli::cli_format_method({
    cli::cli_rule(left = "cier_simulate")
    cli::cli_text("Respondents: {n} x {p} items ({n_scales} scale{?s}) -- \\
                   {n_careless} careless ({pct}%).")
    cli::cli_text("Patterns: {pattern_line}.")
    if (n_careless > 0L) {
      extent_line <- sim_print_extents(x$truth)
      cli::cli_text("Extent: {extent_line}.")
    }
    cli::cli_text("Timing: $seconds (totals) + $page_seconds \\
                   ({n_pages} page{?s}: {pages_word} items).")
    if (n_checks > 0L) {
      cli::cli_text("Checks: {n_checks} attention check{?s} in $checks \\
                     (pass sets in $pass).")
    } else {
      cli::cli_text("Checks: none.")
    }
    cli::cli_text("Truth: $truth -- careless, pattern, extent, onset_item, \\
                   offset_item, speeded, params.")
    cli::cli_text("Simulated data (power analysis / method comparison / \\
                   recovery tests), not evidence of real-world validity.")
  })
  cat(lines, sep = "\n")
  cat("\n")
  invisible(x)
}

#' Summarise a simulated C/IER dataset
#'
#' Delegates to [print()][print.cier_sim], whose output already is the object's summary.
#'
#' @param object A `cier_sim` object from [cier_simulate()].
#' @param ... Ignored.
#' @return `object`, invisibly.
#' @export
summary.cier_sim <- function(object, ...) {
  print(object)
}
