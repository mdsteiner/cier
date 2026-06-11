# Purpose: The single cutoff path for the whole package and the flag
#          applicator that consumes it. `resolve_cutoff()` turns a vector of
#          per-respondent index values into one numeric cutoff; `apply_flag()`
#          turns that cutoff into per-respondent logical flags;
#          `resolve_index_cutoff()` composes the two with `new_cier_index()` into
#          the shared cutoff -> flag -> assemble tail used by the percentile
#          indices. The cutoff math is pure, with no I/O and no registry read --
#          each index wrapper reads its own registry row and passes the method /
#          direction / rate it needs.
# Args:    See per-function documentation below.
# Returns: resolve_cutoff() -> numeric scalar; apply_flag() -> logical vector;
#          resolve_index_cutoff() -> a cier_index.
# Invariants:
#   - These are INTERNAL resolvers that TRUST their inputs. Every public index
#     wrapper validates the user-supplied rate (`fpr` / `alpha` / `frac`) and
#     literal `cutoff` before calling (early fail), and `method` / `direction`
#     come from the registry, so no input re-checking happens here. The only
#     condition raised below is the runtime percentile abstention.
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

# Median-relative cutoff for the timing family: flag respondents faster than a
# fraction `frac` of the sample median (Leiner 2019 Relative Speed Index; Greszki
# et al. 2015). The cutoff is `frac * median(finite values)`. Like the percentile
# resolver this drops non-finite values first and abstains -- NA plus the same
# typed `cier_warning_insufficient_items` -- when none remain, so the all-missing
# case mirrors the percentile path rather than returning a silent NaN. This is an
# OVERRIDE resolver (median-anchored, data-dependent), dispatched inline by the
# total-time wrapper, not a registry default routed through resolve_cutoff().
# Internal: the wrapper has validated `frac` in (0, 1], so it does no input check.
resolve_median_cutoff <- function(value, frac, call = rlang::caller_env()) {
  finite <- value[is.finite(value)]
  if (length(finite) == 0L) {
    cier_warn(
      "cier_warning_insufficient_items",
      c("Cannot resolve a median-relative cutoff: no finite values remain.",
        "i" = "Returning {.val NA} as the cutoff."),
      data = list(n_used = 0L), call = call
    )
    return(NA_real_)
  }
  frac * stats::median(finite)
}

# Pure parameter-free Kneedle elbow (Satopaa, Albrecht, Irwin & Raghavan 2011)
# for a convex, increasing curve -- the right-skewed C/IER index distribution,
# where a moderate bulk gives way to a high-scoring careless tail. The knee is
# the point of greatest distance BELOW the line joining the sorted curve's
# endpoints: sort the values ascending, normalise both axes to [0, 1], and take
# the value where `y_norm - x_norm` is most negative (the argmin). `values` must
# be finite, length >= 3, and non-constant -- resolve_kneedle_cutoff() guarantees
# this, so the kernel does no input checking. There is no sensitivity parameter
# and no smoothing (the raw sorted scores), so it is byte-identical to the
# independent oracle (tests/reference/ref-kneedle-satopaa-2011.R) at tolerance 0.
kneedle_knee <- function(values) {
  sorted <- sort(values)
  n <- length(sorted)
  x_norm <- (seq_len(n) - 1L) / (n - 1L)
  y_norm <- (sorted - sorted[[1L]]) / (sorted[[n]] - sorted[[1L]])
  sorted[[which.min(y_norm - x_norm)]]
}

# Kneedle elbow cutoff for the upper-tail Laz.R index (cier_lazr, kneedle = TRUE):
# the sample-specific cutoff Biemann et al. (2025) recommend in their companion
# app -- the parameter-free Satopaa et al. (2011) elbow of the sorted scores.
# Drops non-finite values first, then abstains -- NA plus the shared
# cier_warning_insufficient_items, so an abstaining cutoff flags nobody -- when
# fewer than three finite values remain OR they are all equal (a knee is
# undefined on a degenerate distribution; the signed-off deviation from the
# archived version, which returned the constant value). Like resolve_median_cutoff
# this is an OVERRIDE resolver (data-driven, dispatched inline by the wrapper and
# passed to resolve_index_cutoff() as the cutoff), not a registry default routed
# through resolve_cutoff(). It computes the upper-tail elbow only -- the sole
# wired use (cier_lazr is upper); the help page documents the convention.
resolve_kneedle_cutoff <- function(value, call = rlang::caller_env()) {
  finite <- value[is.finite(value)]
  # `||` short-circuits, so min/max are never taken on a < 3-length vector; for
  # length >= 3, `min == max` is the all-equal test and also states kneedle_knee's
  # precondition that its `max - min` denominator is non-zero.
  if (length(finite) < 3L || min(finite) == max(finite)) {
    cier_warn(
      "cier_warning_insufficient_items",
      c("Cannot resolve a Kneedle cutoff: fewer than three finite scores, or \\
         all scores equal -- the sorted curve has no knee.",
        "i" = "Returning {.val NA} as the cutoff."),
      data = list(n_used = length(finite)), call = call
    )
    return(NA_real_)
  }
  kneedle_knee(finite)
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

# Purpose: Resolve a flagging cutoff from per-respondent index values. Internal
#   resolver; it TRUSTS its inputs (the calling wrapper has validated the
#   user-supplied rate / `cutoff`, and `method` / `direction` come from the
#   registry), so it does no input checking -- it only does the math and signals
#   the runtime percentile abstention.
# Args:
#   values    - numeric vector of index values (NA / NaN / Inf dropped before
#               the percentile quantile); ignored for "fixed" / "chisq".
#   direction - "upper" (flag the high tail) or "lower" (the low tail); used
#               only by the percentile method's single flip.
#   method    - "percentile" (default), "fixed", or "chisq".
#   fpr       - false-positive tail mass for "percentile" (default 0.05).
#   df        - chi-square degrees of freedom (method = "chisq").
#   alpha     - chi-square upper-tail probability (default 0.001).
#   value     - the cutoff value for method = "fixed": a literal threshold when
#               `n_items` is NULL, or a fraction of the item count when
#               `n_items` is supplied.
#   n_items   - number of items; when supplied, `value` is a fraction of the
#               item count (`ceiling(value * n_items)`) for method = "fixed".
#   call      - calling environment for the abstention warning.
# Returns: a numeric scalar cutoff; NA_real_ when the percentile method abstains.
resolve_cutoff <- function(values = NULL, direction = "upper",
                           method = "percentile", fpr = 0.05,
                           df = NULL, alpha = 0.001, value = NULL,
                           n_items = NULL, call = rlang::caller_env()) {
  if (identical(method, "fixed")) {
    return(resolve_fixed_cutoff(value, n_items))
  }
  if (identical(method, "chisq")) {
    return(as.numeric(stats::qchisq(1 - alpha, df = df)))
  }
  if (identical(method, "percentile")) {
    return(resolve_percentile_cutoff(values, direction, fpr, call))
  }
  # Any other registry method (e.g. "perfit_null") is model-referenced and must
  # be resolved at its bridge from the fitted object -- it cannot be reached here.
  # Erroring (rather than silently falling through to the percentile branch)
  # guards a future caller that wires a non-value cutoff method into this path.
  cier_abort(
    "cier_error_state",
    c("Cutoff method {.val {method}} is not resolved by {.fn resolve_cutoff}.",
      "i" = "Model-referenced cutoffs (e.g. {.val perfit_null}) are resolved at \\
             the bridge from the fitted object, not the value-only path."),
    data = list(arg = "method", observed = method), call = call
  )
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

# Purpose: Resolve the flagging cutoff and assemble the cier_index for a
#   percentile index -- the shared cutoff -> flag -> assemble tail of cier_irv,
#   cier_even_odd, cier_person_total, and cier_personal_reliability.
# Args:
#   value  - numeric per-respondent index values (NA where the row abstains).
#   row    - the method's registry row (carries flag_direction,
#            default_cutoff_method, default_cutoff_value, method).
#   fpr    - the user's `fpr` override or NULL (NULL uses the registry default).
#   cutoff - the user's literal `cutoff` override or NULL. Already validated by
#            the wrapper; when supplied it passes through verbatim, otherwise the
#            default is the empirical `fpr` percentile in the flag direction.
#   call   - calling environment for the percentile-abstention warning.
# Returns: a `cier_index` (via new_cier_index()).
resolve_index_cutoff <- function(value, row, fpr, cutoff,
                                 call = rlang::caller_env()) {
  cutoff_value <- if (!is.null(cutoff)) {
    cutoff
  } else {
    resolve_cutoff(values = value, direction = row$flag_direction,
                   method = row$default_cutoff_method,
                   fpr = if (is.null(fpr)) row$default_cutoff_value else fpr,
                   call = call)
  }
  flagged <- apply_flag(value, cutoff_value, row$flag_direction, call = call)
  new_cier_index(value, flagged, row$method, cutoff_value, row$flag_direction)
}
