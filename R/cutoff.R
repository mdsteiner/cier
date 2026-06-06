# Purpose: The single cutoff path for the whole package and the flag
#          applicator that consumes it. `resolve_cutoff()` turns a vector of
#          per-respondent index values into one numeric cutoff; `apply_flag()`
#          turns that cutoff into per-respondent logical flags. Pure math, no
#          I/O and no registry read -- each index wrapper reads its own registry
#          row and passes the method / direction / fpr it needs.
# Args:    See per-function documentation below.
# Returns: resolve_cutoff() -> numeric scalar; apply_flag() -> logical vector.
# Invariants:
#   - The percentile branch flips on direction exactly ONCE (upper -> the
#     1 - fpr quantile, lower -> the fpr quantile). The registry stores each
#     method's literal directional quantile, so it must never be fed back
#     through the flip; the flag comparator is applied separately by
#     apply_flag(). quantile type is 7 (R's default), pinned against base R in
#     tests/reference/TOLERANCES.md.
#   - An abstaining percentile cutoff is NA_real_ with a typed warning, not an
#     error; an NA cutoff flags nobody.

# Empirical-percentile cutoff at the target false-positive tail mass. Single
# direction flip: upper flags the high tail (probs = 1 - fpr), lower flags the
# low tail (probs = fpr). Non-finite values are dropped first; with none left
# the cutoff abstains (NA + typed warning) rather than erroring.
resolve_percentile_cutoff <- function(values, direction, fpr, call) {
  if (!is.numeric(fpr) || length(fpr) != 1L || is.na(fpr) ||
        fpr <= 0 || fpr >= 1) {
    cier_abort(
      "cier_error_input",
      "{.arg fpr} must be a single number in the open interval (0, 1).",
      data = list(arg = "fpr", observed = fpr), call = call
    )
  }
  finite <- values[is.finite(values)]
  if (length(finite) == 0L) {
    cier_warn(
      "cier_warning_insufficient_items",
      c("Cannot resolve a percentile cutoff: no finite values remain.",
        "i" = "Returning {.val NA} as the cutoff."),
      data = list(n_used = 0L), call = call
    )
    return(NA_real_)
  }
  probs <- if (identical(direction, "upper")) 1 - fpr else fpr
  as.numeric(stats::quantile(finite, probs = probs, names = FALSE, type = 7L))
}

# Resolve a fixed cutoff that may be expressed as a fraction of the item count.
# With `n_items` supplied, a finite scalar in (0, 1] is a fraction of the item
# count (`ceiling(value * n_items)`, so `1` resolves to all items) and a scalar
# > 1 is an absolute count (passed through literally, including non-integers);
# `<= 0` or `> n_items` abort. Reuses check_number for the single-finite-scalar
# guard so the diagnostics match the rest of the package. The longstring default
# 0.5 and any user `cutoff` both flow through here.
resolve_fixed_cutoff <- function(value, n_items, call) {
  if (is.null(n_items)) {
    return(as.numeric(value))
  }
  check_number(value, "cutoff", call = call)
  if (value <= 0) {
    cier_abort("cier_error_input", "{.arg cutoff} must be positive.",
               data = list(arg = "cutoff", observed = value), call = call)
  }
  if (value <= 1) {
    return(ceiling(value * n_items))
  }
  if (value > n_items) {
    cier_abort(
      "cier_error_input",
      c("{.arg cutoff} cannot exceed the number of items.",
        "x" = "Got {value}; the data has {n_items} item{?s}."),
      data = list(arg = "cutoff", observed = value, expected = n_items),
      call = call
    )
  }
  as.numeric(value)
}

# Purpose: Resolve a flagging cutoff from per-respondent index values.
# Args:
#   values    - numeric vector of index values (NA / NaN / Inf dropped before
#               the percentile quantile); ignored for "fixed" / "chisq".
#   direction - "upper" (flag the high tail) or "lower" (the low tail); used
#               only by the percentile method's single flip.
#   method    - "percentile" (default), "fixed", or "chisq".
#   fpr       - false-positive tail mass for "percentile" (default 0.05); must
#               lie in the open interval (0, 1).
#   df        - chi-square degrees of freedom (required when method = "chisq").
#   alpha     - chi-square upper-tail probability (default 0.001).
#   value     - the cutoff for method = "fixed" (a literal count, or a
#               fraction-or-count when `n_items` is supplied).
#   n_items   - number of items; when supplied, `value` is interpreted as a
#               fraction-of-items-or-absolute-count for method = "fixed".
#   call      - calling environment for typed conditions.
# Returns: a numeric scalar cutoff; NA_real_ when the percentile method abstains.
resolve_cutoff <- function(values = NULL, direction = "upper",
                           method = "percentile", fpr = 0.05,
                           df = NULL, alpha = 0.001, value = NULL,
                           n_items = NULL, call = rlang::caller_env()) {
  check_choice(method, "method", cier_cutoff_methods(), call = call)
  if (identical(method, "fixed")) {
    if (is.null(value)) {
      cier_abort("cier_error_input",
                 "{.arg value} is required for the fixed cutoff method.",
                 data = list(arg = "value"), call = call)
    }
    return(resolve_fixed_cutoff(value, n_items, call))
  }
  if (identical(method, "chisq")) {
    if (is.null(df)) {
      cier_abort("cier_error_input",
                 "{.arg df} is required for the chi-square cutoff method.",
                 data = list(arg = "df"), call = call)
    }
    check_number(df, "df", lower = 1, call = call)
    check_number(alpha, "alpha", lower = 0, upper = 1, call = call)
    return(as.numeric(stats::qchisq(1 - alpha, df = df)))
  }
  check_choice(direction, "direction", cier_flag_directions(), call = call)
  resolve_percentile_cutoff(values, direction, fpr, call)
}

# Purpose: Apply a resolved cutoff to per-respondent index values.
# Args:
#   values    - numeric vector of index values.
#   cutoff    - numeric scalar cutoff (NA_real_ means the index abstained).
#   direction - "upper" (flag values >= cutoff) or "lower" (flag <= cutoff).
#   call      - calling environment for typed conditions.
# Returns: a logical vector the length of `values`. NA cutoff -> all FALSE;
#   NA values are never flagged; ties at the cutoff are flagged.
apply_flag <- function(values, cutoff, direction, call = rlang::caller_env()) {
  check_choice(direction, "direction", cier_flag_directions(), call = call)
  if (is.na(cutoff)) {
    return(rep(FALSE, length(values)))
  }
  comparator <- if (identical(direction, "upper")) `>=` else `<=`
  !is.na(values) & comparator(values, cutoff)
}
