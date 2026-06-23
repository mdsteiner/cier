# The cutoff -> flag -> assemble tail of the percentile indices: resolve_cutoff() turns
# index values into one numeric cutoff, apply_flag() into logical flags,
# resolve_index_cutoff() composes both with new_cier_index().
#
# These resolvers TRUST their inputs (the wrapper validated the rate and literal
# `cutoff`; method/direction come from the method spec). The only conditions are the
# runtime abstention (cier_warning_insufficient_items; cutoff is NA_real_ and flags
# nobody) and the saturation diagnostic (cier_warning_saturated_cutoff; warns but
# resolves).

# Empirical-percentile cutoff at the target false-positive tail mass. Upper flags the
# high tail (probs = 1 - fpr), lower the low tail (probs = fpr). Drops non-finite
# values, then three degeneracy guards abstain.
resolve_percentile_cutoff <- function(values, direction, fpr, call) {
  finite <- values[is.finite(values)]
  n_used <- length(finite)
  # round() before ceiling() neutralises sub-1e-9 IEEE noise in 1 / fpr so a
  # reciprocal a hair above an integer is not pushed up.
  min_n <- ceiling(round(1 / fpr, 9L))
  if (n_used < min_n) {
    cier_warn(
      "cier_warning_insufficient_items",
      c("Cannot resolve a percentile cutoff at fpr = {fpr}: only {n_used} \\
         finite score{?s} (at least {min_n} are needed for a {fpr} tail).",
        "i" = "Returning {.val NA} as the cutoff; no respondent is flagged."),
      data = list(n_used = n_used, n_required = min_n, fpr = fpr), call = call
    )
    return(NA_real_)
  }
  upper <- identical(direction, "upper")
  probs <- if (upper) 1 - fpr else fpr
  cutoff <- as.numeric(stats::quantile(finite, probs = probs, names = FALSE,
                                       type = 7L))
  flags_all <- if (upper) cutoff <= min(finite) else cutoff >= max(finite)
  if (flags_all) {
    cier_warn(
      "cier_warning_insufficient_items",
      c("Cannot resolve a percentile cutoff: the scores are effectively \\
         constant, so every scored respondent would be flagged.",
        "i" = "Returning {.val NA} as the cutoff; no respondent is flagged."),
      data = list(n_used = n_used), call = call
    )
    return(NA_real_)
  }
  on_extreme <- if (upper) cutoff == max(finite) else cutoff == min(finite)
  if (on_extreme) {
    extreme <- if (upper) "maximum" else "minimum"
    cier_warn(
      "cier_warning_saturated_cutoff",
      c("The percentile cutoff ({.val {cutoff}}) equals the score {extreme}; \\
         a tie mass at the cutoff flags more than fpr = {fpr} of respondents.",
        "i" = "Ties at the cutoff can substantially exceed the target tail \\
               (the documented ranking convention). Inspect the realised flag \\
               rate, or set an explicit {.arg cutoff}."),
      data = list(cutoff = cutoff, fpr = fpr, direction = direction),
      call = call
    )
  }
  cutoff
}

# Median-relative cutoff for the timing family: flag respondents faster than `frac` of
# the sample median (Leiner 2019 Relative Speed Index; Greszki et al. 2015). Drops
# non-finite values, abstains when none remain. Override resolver dispatched inline by
# the total-time wrapper, not routed through resolve_cutoff().
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

# Parameter-free Kneedle elbow (Satopaa et al. 2011) for a convex increasing curve --
# the right-skewed C/IER distribution. The knee is the point of greatest distance below
# the line joining the sorted endpoints: normalise both axes to [0, 1], take the value
# where `y_norm - x_norm` is most negative. Precondition: `values` finite, length >= 3,
# non-constant (resolve_kneedle_cutoff guarantees this).
kneedle_knee <- function(values) {
  sorted <- sort(values)
  n <- length(sorted)
  x_norm <- (seq_len(n) - 1L) / (n - 1L)
  y_norm <- (sorted - sorted[[1L]]) / (sorted[[n]] - sorted[[1L]])
  sorted[[which.min(y_norm - x_norm)]]
}

# Integer flag-rate percent for the Kneedle saturation message; reaches 100 only when
# EVERY scored respondent is flagged (plain round() would show 100% for 999/1000).
saturation_pct <- function(n_flagged, n_scored) {
  pct <- round(100 * n_flagged / n_scored)
  if (pct >= 100 && n_flagged < n_scored) 99 else pct
}

# Kneedle elbow cutoff for the upper-tail Laz.R index (kneedle = TRUE). Drops non-finite
# values, abstains when fewer than three remain or all are equal. Otherwise resolves the
# elbow, warning (cier_warning_saturated_cutoff; warn but resolve) when it flags more
# than half the scored respondents -- the bimodal case where the knee sits at or below
# the careful bulk. Override resolver dispatched inline by the wrapper, not routed
# through resolve_cutoff().
resolve_kneedle_cutoff <- function(value, call = rlang::caller_env()) {
  finite <- value[is.finite(value)]
  # `||` short-circuits, so min/max are never taken on a < 3-length vector; for
  # length >= 3, `min == max` is the all-equal test.
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
  knee <- kneedle_knee(finite)
  # Realised flag rate with apply_flag's upper `>=`. Strict `> 0.5` leaves the exact-half
  # boundary unwarned; a bimodal sample lands the elbow on the lower cluster, flagging a
  # majority, which warns.
  n_scored <- length(finite)
  n_flagged <- sum(finite >= knee)
  rate <- n_flagged / n_scored
  if (rate > 0.5) {
    pct <- saturation_pct(n_flagged, n_scored)
    cier_warn(
      "cier_warning_saturated_cutoff",
      c("The Kneedle elbow ({.val {knee}}) flags {n_flagged} of {n_scored} \\
         scored respondents ({pct}%) -- more than half.",
        "i" = "A majority flag rate can reflect a genuinely low-quality sample, \\
               but also arises when the elbow sits at or below the careful bulk \\
               of a strongly bimodal distribution. Inspect the realised rate, or \\
               set a fixed top-tail share with {.arg fpr}."),
      data = list(cutoff = knee, n_flagged = n_flagged, n_scored = n_scored,
                  rate = rate),
      call = call
    )
  }
  knee
}

# Resolve a fixed cutoff. With `n_items` the value is a fraction of the item count
# (`ceiling(value * n_items)`); without it a literal threshold. Rounds before `ceiling()`
# to remove IEEE-754 noise (0.28 * 25 == 7.0000000000000009 -> 7).
resolve_fixed_cutoff <- function(value, n_items) {
  if (is.null(n_items)) {
    return(as.numeric(value))
  }
  ceiling(round(as.numeric(value) * n_items, 9L))
}

# Resolve a flagging cutoff from per-respondent index values. Dispatches on `method`:
# "percentile" (quantile at `fpr` tail in `direction`), "fixed" (literal threshold, or a
# fraction of the item count when `n_items` is given), or "chisq" (qchisq tail at
# `alpha`, df = `df`). Returns a numeric scalar, NA_real_ when percentile abstains.
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
  # Other methods (e.g. "mc_null") are null-referenced and resolved at their bridge from
  # the scored response block; erroring guards a future caller that wires one here.
  cier_abort(
    "cier_error_state",
    c("Cutoff method {.val {method}} is not resolved by {.fn resolve_cutoff}.",
      "i" = "Null-referenced cutoffs (e.g. {.val mc_null}) are resolved at \\
             their bridge from the scored response block, not the value-only path."),
    data = list(arg = "method", observed = method), call = call
  )
}

# Apply a resolved cutoff to index values: "upper" flags values >= cutoff, "lower"
# flags <= cutoff. Returns a logical vector the length of `values`; an NA cutoff
# (the index abstained) -> all FALSE, NA values are never flagged, ties are flagged.
apply_flag <- function(values, cutoff, direction, call = rlang::caller_env()) {
  check_choice(direction, "direction", cier_flag_directions(), call = call)
  if (is.na(cutoff)) {
    return(rep(FALSE, length(values)))
  }
  comparator <- if (identical(direction, "upper")) `>=` else `<=`
  !is.na(values) & comparator(values, cutoff)
}

# The cutoff-provenance pair (method label + rate): a literal `cutoff` override is
# "literal" (no rate), otherwise the index's default strategy + rate.
index_provenance <- function(cutoff, default_method, default_rate) {
  if (!is.null(cutoff)) {
    list(method = "literal", rate = NA_real_)
  } else {
    list(method = default_method, rate = default_rate)
  }
}

# Assemble a cier_index for a wrapper that resolved its own cutoff: derive the
# provenance (literal override vs default strategy) and construct.
new_cier_index_provenance <- function(value, flagged, method, cutoff_value,
                                      direction, cutoff, default_method,
                                      default_rate) {
  prov <- index_provenance(cutoff, default_method, default_rate)
  new_cier_index(value, flagged, method, cutoff_value, direction,
                 cutoff_method = prov$method, cutoff_rate = prov$rate)
}

# Flag via apply_flag() then assemble with provenance: the shared flag -> assemble tail
# of the wrappers that resolve their own cutoff. (cier_page_time and cier_attention keep
# an explicit apply_flag(): their all-flagged warning reads `flagged` before assembly.)
flag_and_assemble <- function(value, cutoff_value, row, cutoff,
                              default_method, default_rate,
                              call = rlang::caller_env()) {
  flagged <- apply_flag(value, cutoff_value, row$flag_direction, call = call)
  new_cier_index_provenance(value, flagged, row$method, cutoff_value,
                            row$flag_direction, cutoff, default_method,
                            default_rate)
}

# Resolve the cutoff and assemble the cier_index for a percentile index -- the shared
# cutoff -> flag -> assemble tail of cier_irv, cier_even_odd, cier_person_total,
# cier_personal_reliability. A literal `cutoff` passes through verbatim, otherwise the
# empirical `fpr` percentile. `provenance` optionally overrides the recorded label/rate
# (total_time / lazr pass it for their median-relative / Kneedle overrides, which arrive
# as `cutoff`).
resolve_index_cutoff <- function(value, row, fpr, cutoff, provenance = NULL,
                                 call = rlang::caller_env()) {
  rate <- if (is.null(fpr)) row$default_cutoff_value else fpr
  cutoff_value <- if (!is.null(cutoff)) {
    cutoff
  } else {
    resolve_cutoff(values = value, direction = row$flag_direction,
                   method = row$default_cutoff_method, fpr = rate, call = call)
  }
  prov <- if (!is.null(provenance)) {
    provenance
  } else {
    index_provenance(cutoff, row$default_cutoff_method, rate)
  }
  flagged <- apply_flag(value, cutoff_value, row$flag_direction, call = call)
  new_cier_index(value, flagged, row$method, cutoff_value, row$flag_direction,
                 cutoff_method = prov$method, cutoff_rate = prov$rate)
}
