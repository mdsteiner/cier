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

# Resolve a fixed cutoff. With `n_items` supplied the value is a fraction of the
# item count (`ceiling(value * n_items)`); without it the value is a literal
# threshold on the score and passes through verbatim. The two modes have
# distinct wrapper arguments (`frac` vs `cutoff`), so the wrapper -- which knows
# the argument name and its valid range -- validates before calling; this kernel
# only does the arithmetic. The product is rounded to remove floating-point
# noise before `ceiling()` so e.g. a fraction of 0.28 over 25 items resolves to
# 7 rather than 8 (0.28 * 25 == 7.0000000000000009 in IEEE-754).
resolve_fixed_cutoff <- function(value, n_items) {
  if (is.null(n_items)) {
    return(as.numeric(value))
  }
  ceiling(round(as.numeric(value) * n_items, 9L))
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
#   value     - the cutoff for method = "fixed": a literal threshold when
#               `n_items` is NULL, or a fraction of the item count when
#               `n_items` is supplied. Must be pre-validated by the caller (the
#               wrapper or resolve_index_cutoff); the fixed arithmetic does no
#               checking.
#   n_items   - number of items; when supplied, `value` is a fraction of the
#               item count (`ceiling(value * n_items)`) for method = "fixed".
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
    return(resolve_fixed_cutoff(value, n_items))
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

# Purpose: Resolve an index wrapper's flagging cutoff, honouring the
#   mutually-exclusive rate-vs-literal override every index exposes. Keeps the
#   rate-or-literal dispatch in one place so each index wrapper (and the eight
#   still to come) does not re-derive it. The wrapper validates `rate` and
#   `cutoff` (check_fraction / check_number) before calling; this helper only
#   enforces the mutual exclusion and routes a literal cutoff through the single
#   resolver.
# Args:
#   rate      - the index's rate knob (`fpr` / `alpha` / `frac`), or NULL.
#   rate_name - its argument name, for the typed mutual-exclusion message.
#   cutoff    - the literal-threshold knob (already validated), or NULL.
#   rate_fn   - a no-argument closure that resolves the rate-based default when
#               no literal `cutoff` is supplied (the index's own percentile /
#               chi-square / fraction resolution).
#   call      - calling environment for typed conditions.
# Returns: a numeric scalar cutoff.
resolve_index_cutoff <- function(rate, rate_name, cutoff, rate_fn,
                                 call = rlang::caller_env()) {
  assert_single_override(rate, rate_name, cutoff, call)
  if (is.null(cutoff)) {
    return(rate_fn())
  }
  resolve_cutoff(method = "fixed", value = cutoff, call = call)
}
