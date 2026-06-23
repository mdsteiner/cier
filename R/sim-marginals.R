# Resolve a target marginal response distribution into per-item GRM thresholds for the
# attentive generator (sim-attentive.R), via tau_k = qnorm(cumsum(p)) so under normal traits
# the categorised marginal equals the target pmf. Pure; mutually exclusive with raw
# `thresholds`.
# Every pmf is strictly positive and re-normalised to sum exactly to 1, so a raw sum a few
# ulp above 1 cannot push qnorm to +Inf / NaN.

# The five named marginal-shape presets, as a strictly-positive pmf of length k. `strength`
# (>= 1, default 1.5) steepens the geometric skew and is ignored for the non-skewed presets.
#   - uniform: flat 1 / k.
#   - peaked: symmetric triangular w_j = min(j, k + 1 - j); k = 5 -> (1,2,3,2,1)/9.
#   - skewed_left: geometric w_j = strength^(j - 1), mass on the high (agreement) categories;
#     skewed_right is its mirror.
#   - bimodal: deep reflected-triangular U; k = 5 -> (3,2,1,2,3)/11.
sim_preset_pmf <- function(preset, k, strength = 1.5,
                           call = rlang::caller_env()) {
  check_choice(preset, "preset",
               c("uniform", "peaked", "skewed_left", "skewed_right", "bimodal"),
               call = call)
  k <- as.integer(k)
  tri <- pmin(seq_len(k), rev(seq_len(k)))
  w <- switch(preset,
    uniform      = rep(1, k),
    peaked       = tri,
    skewed_left  = strength^(seq_len(k) - 1L),
    skewed_right = rev(strength^(seq_len(k) - 1L)),
    bimodal      = (max(tri) + 1) - tri
  )
  w / sum(w)
}

# Validate an explicit pmf spec for an item with k categories: a strictly-positive
# numeric vector of length k summing to 1 (a zero entry would push a threshold to
# +/- Inf and break monotonicity; a negative entry is not a probability).
sim_validate_pmf <- function(p, k, call) {
  ok <- is.numeric(p) && length(p) == k && all(is.finite(p)) &&
    all(p > 0) && abs(sum(p) - 1) <= 1e-8
  if (!ok) {
    cier_abort(
      "cier_error_input",
      c("A {.arg marginals} pmf must be {k} strictly-positive probabilities \\
         summing to 1.",
        "x" = "Got length {length(p)}, sum {round(sum(p), 6)}."),
      data = list(arg = "marginals", expected = k), call = call
    )
  }
  invisible(p)
}

# Convert one spec (a preset string, an explicit pmf vector, or
# list(preset = , strength = )) to a pmf for an item with k categories.
sim_spec_to_pmf <- function(spec, k, call) {
  if (is.character(spec) && length(spec) == 1L) {
    return(sim_preset_pmf(spec, k, call = call))
  }
  if (is.numeric(spec)) {
    sim_validate_pmf(spec, k, call)
    return(as.numeric(spec))
  }
  if (is.list(spec) && !is.null(spec[["preset"]])) {
    extra <- setdiff(names(spec), c("preset", "strength"))
    if (length(extra) > 0L) {
      cier_abort(
        "cier_error_input",
        c("A {.code list(preset = , strength = )} spec has unknown key(s): \\
           {.val {extra}}.",
          "i" = "Only {.field preset} and {.field strength} are read."),
        data = list(arg = "marginals", observed = extra), call = call
      )
    }
    strength <- spec[["strength"]]
    if (is.null(strength)) {
      return(sim_preset_pmf(spec[["preset"]], k, call = call))
    }
    check_number(strength, "strength", lower = 1, call = call)
    return(sim_preset_pmf(spec[["preset"]], k, strength = strength, call = call))
  }
  cier_abort(
    "cier_error_input",
    "Each {.arg marginals} spec must be a preset name, a pmf vector, or \\
     {.code list(preset = , strength = )}.",
    data = list(arg = "marginals"), call = call
  )
}

# Validate the names of a scale-named `marginals` list: every entry named (no silent
# positional fallback), no repeated scale (no silent first-wins), and exactly the
# battery's scales covered.
sim_validate_scale_names <- function(nm, scales, call) {
  if (!all(nzchar(nm))) {
    cier_abort(
      "cier_error_input",
      c("A named {.arg marginals} list must name every entry (one per scale).",
        "x" = "Some entries are unnamed.",
        "i" = "Name all entries by scale, or pass an unnamed list of per-item \\
               specs."),
      data = list(arg = "marginals"), call = call
    )
  }
  if (anyDuplicated(nm) > 0L) {
    cier_abort(
      "cier_error_input",
      c("A scale-named {.arg marginals} list must not repeat a scale.",
        "x" = "Duplicated: {.val {unique(nm[duplicated(nm)])}}."),
      data = list(arg = "marginals"), call = call
    )
  }
  if (!setequal(nm, unique(scales))) {
    cier_abort(
      "cier_error_input",
      c("A scale-named {.arg marginals} list must name exactly the scales.",
        "x" = "Scales: {.val {unique(scales)}}; named: {.val {nm}}."),
      data = list(arg = "marginals"), call = call
    )
  }
  invisible(NULL)
}

# Expand `marginals` into a per-item list of specs (length p). Disambiguation: NULL ->
# peaked everywhere; a single preset / pmf / `preset`-keyed list -> one spec for every item;
# a fully scale-named list -> per scale; a fully unnamed length-p list -> per item. A
# partially-named list is rejected.
sim_marginal_specs <- function(marginals, scales, p, call) {
  if (is.null(marginals)) {
    return(rep(list("peaked"), p))
  }
  if (is.character(marginals) && length(marginals) == 1L) {
    return(rep(list(marginals), p))
  }
  if (is.numeric(marginals)) {
    return(rep(list(marginals), p))
  }
  if (is.list(marginals)) {
    if (!is.null(marginals[["preset"]])) {
      return(rep(list(marginals), p))
    }
    nm <- names(marginals)
    if (is.null(nm)) {
      if (length(marginals) == p) {
        return(marginals)
      }
      cier_abort(
        "cier_error_input",
        c("An unnamed {.arg marginals} list needs one spec per item.",
          "x" = "Got {length(marginals)} spec{?s} for {p} item{?s}."),
        data = list(arg = "marginals", expected = p), call = call
      )
    }
    sim_validate_scale_names(nm, scales, call)
    return(lapply(scales, function(s) marginals[[s]]))
  }
  cier_abort(
    "cier_error_input",
    c("{.arg marginals} must be a preset, a pmf vector, a scale-named list, \\
       or an unnamed list with one spec per item.",
      "x" = "Could not interpret the value for {p} item{?s}."),
    data = list(arg = "marginals", expected = p), call = call
  )
}

# Resolve `marginals` to the per-item target pmf list (length p). `items` is the
# validated check_items_simulate list (carries `scale` and `categories`).
sim_resolve_marginals <- function(marginals, items, call = rlang::caller_env()) {
  cats <- items$categories
  p <- length(cats)
  specs <- sim_marginal_specs(marginals, items$scale, p, call)
  out <- vector("list", p)
  for (j in seq_len(p)) {
    out[[j]] <- sim_spec_to_pmf(specs[[j]], cats[[j]], call)
  }
  out
}

# Convert a per-item pmf list to per-item GRM thresholds: tau_k = qnorm(cumsum(p)) for
# k = 1..K-1 (the final cut qnorm(1) = Inf is dropped). Normalised to sum exactly to 1 first,
# so a raw sum a few ulp above 1 cannot send qnorm to +Inf / NaN.
sim_marginals_to_thresholds <- function(pmf_list) {
  lapply(pmf_list, function(p) {
    p <- p / sum(p)
    k <- length(p)
    stats::qnorm(cumsum(p)[seq_len(k - 1L)])
  })
}

# Dispatch the attentive thresholds: raw `thresholds` and `marginals` are mutually
# exclusive. An explicit `thresholds` list passes through after validation,
# otherwise the (possibly NULL, meaning peaked) `marginals` is resolved to
# thresholds. `items` is the validated check_items_simulate list.
sim_resolve_thresholds <- function(items, marginals, thresholds,
                                   call = rlang::caller_env()) {
  if (!is.null(thresholds) && !is.null(marginals)) {
    cier_abort(
      "cier_error_input",
      c("Supply only one of {.arg marginals} and {.arg thresholds}.",
        "i" = "They are two ways to set the per-item category distribution."),
      data = list(args = c("marginals", "thresholds")), call = call
    )
  }
  if (!is.null(thresholds)) {
    sim_validate_thresholds(thresholds, items, call)
    return(thresholds)
  }
  sim_marginals_to_thresholds(sim_resolve_marginals(marginals, items, call))
}
