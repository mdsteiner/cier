# Validators for the per-item `items` metadata frame (scale, reverse-keying, min/max,
# span homogeneity, positional alignment to response columns). Raise `cier_error_input`
# on a malformed frame (heterogeneous-span case also carries `cier_error_backend_limit`);
# return a normalized list/column on success. Scalar/argument validators are in
# R/checks-args.R; payload coercers in R/checks-data.R.

# `scale`: required, character-coercible, every item labelled, >= min_scales distinct
# labels (the split-half indices correlate across scales).
check_items_scale <- function(items, min_scales, arg, call) {
  if (is.null(items$scale)) {
    cier_abort(
      "cier_error_input",
      "{.arg {arg}} must have a {.field scale} column (one label per item).",
      data = list(arg = arg), call = call
    )
  }
  scale <- as.character(items$scale)
  if (anyNA(scale) || any(!nzchar(scale))) {
    cier_abort(
      "cier_error_input",
      "{.field scale} must be a non-missing label for every item.",
      data = list(arg = arg), call = call
    )
  }
  n_scales <- length(unique(scale))
  if (n_scales < min_scales) {
    cier_abort(
      "cier_error_input",
      c("{.arg {arg}} must have at least {min_scales} distinct {.field scale}s.",
        "x" = "Found {n_scales}.",
        "i" = "The split-half consistency indices correlate across scales."),
      data = list(arg = arg, observed = n_scales, expected = min_scales),
      call = call
    )
  }
  scale
}

# `reverse_keyed`: optional logical, length n_items, no NA; defaults all-FALSE.
check_items_reverse <- function(items, n_items, arg, call) {
  rk <- items$reverse_keyed
  if (is.null(rk)) {
    return(rep(FALSE, n_items))
  }
  if (!is.logical(rk) || anyNA(rk)) {
    cier_abort(
      "cier_error_input",
      "{.field reverse_keyed} must be a logical value (no {.val NA}) per item.",
      data = list(arg = arg), call = call
    )
  }
  rk
}

# One-time typed message when an `items` frame omits reverse_keyed entirely (all items
# treated forward-keyed -- a silent mis-scoring risk on a keyed instrument). An explicit
# reverse_keyed = FALSE stays silent. Kept out of check_items_reverse(), which
# cier_simulate() also uses.
inform_if_unkeyed <- function(items, n_items) {
  if (is.data.frame(items) && is.null(items[["reverse_keyed"]])) {
    cier_inform(
      "cier_message_forward_keyed",
      c("{.arg items} has no {.field reverse_keyed} column; treating all \\
         {n_items} item{?s} as forward-keyed.",
        "i" = "Add a logical {.field reverse_keyed} column if the instrument has \\
               reverse-keyed items, or set {.code reverse_keyed = FALSE} to \\
               silence this."),
      data = list(arg = "items", n_items = n_items)
    )
  }
  invisible(NULL)
}

# `max`: largest response option per item. Required (finite whole, >= min + 1) only on
# reverse-keyed items, which reflect via (min + max) - x (`mins` is the resolved base,
# default 1). NA permitted on forward items; column may be absent when nothing is
# reverse-keyed. Returns the column unchanged (or NULL when absent).
check_items_max <- function(items, reverse_keyed, mins, arg, call) {
  maxs <- items$max
  if (!any(reverse_keyed)) {
    return(maxs)
  }
  rev_maxs <- if (is.null(maxs)) NA_real_ else maxs[reverse_keyed]
  # is_finite_whole() rejects NA/NaN/Inf (a non-finite max would poison the reverse
  # columns) and non-integers. Bound is min + 1, not an absolute 2: a 0/1 item is the
  # smallest valid scale.
  ok <- is_finite_whole(rev_maxs) && all(rev_maxs >= mins[reverse_keyed] + 1)
  if (!ok) {
    cier_abort(
      "cier_error_input",
      c("Reverse-keyed items need an integer {.field max} of at least \\
         {.field min} + 1.",
        "i" = "Set {.field max} (the largest response option) for every \\
               reverse-keyed item; {.field min} (the smallest, default 1) is \\
               read alongside it for the reflection (min + max) - x."),
      data = list(arg = arg), call = call
    )
  }
  maxs
}

# `min`: the response-scale base. Optional; defaults to 1. Generalises the reflection
# (min + max) - x so a 0-based or bipolar scale reflects onto itself. Finite whole on
# every reverse-keyed item (any integer base); NA permitted on forward items. Resolved
# before max, whose >= min + 1 bound reads it.
check_items_min <- function(items, reverse_keyed, n_items, arg, call) {
  mins <- items$min
  if (is.null(mins)) {
    return(rep(1L, n_items))
  }
  rev_mins <- mins[reverse_keyed]
  ok <- length(rev_mins) == 0L || is_finite_whole(rev_mins)
  if (!ok) {
    cier_abort(
      "cier_error_input",
      c("Reverse-keyed items need a finite integer {.field min} (the scale base).",
        "i" = "{.field min} is the smallest response option (default 1); set it \\
               for every reverse-keyed item on a 0-based or bipolar scale."),
      data = list(arg = arg), call = call
    )
  }
  mins
}

# Shared precondition every items validator runs first: `items` must be a data.frame with
# one row per item. When `response_names` is supplied and the frame carries explicit row
# identifiers (an `item` column, else character rownames), the two are cross-checked
# positionally, so a reordered metadata frame is a typed error rather than a silent
# reshuffle (the simulator passes response_names = NULL).
check_items_frame <- function(items, n_items, arg, call, response_names = NULL) {
  if (!is.data.frame(items)) {
    cier_abort(
      "cier_error_input",
      c("{.arg {arg}} must be a data.frame (one row per item).",
        "x" = "Got {.cls {class(items)}}."),
      data = list(arg = arg), call = call
    )
  }
  if (nrow(items) != n_items) {
    cier_abort(
      "cier_error_input",
      c("{.arg {arg}} must have one row per item (column of {.arg responses}).",
        "x" = "Got {nrow(items)} item row{?s} for {n_items} column{?s}."),
      data = list(arg = arg, observed = nrow(items), expected = n_items),
      call = call
    )
  }
  check_items_alignment(items, response_names, arg, call)
  invisible(NULL)
}

# Cross-check item row identifiers against response column names, by position. Identifier
# is an explicit `item` column if present, else character rownames -- automatic integer
# rownames are not an identifier, so the check is a no-op there.
check_items_alignment <- function(items, response_names, arg, call) {
  ids <- if (!is.null(items[["item"]])) {
    as.character(items[["item"]])
  } else if (is.character(.row_names_info(items, 0L))) {
    rownames(items)
  } else {
    NULL
  }
  if (is.null(ids) || is.null(response_names)) {
    return(invisible(NULL))
  }
  # NA on either side is a mismatch; `ids != response_names` returns NA there, which
  # which() would silently drop -- so test NA-ness explicitly.
  mism <- which(is.na(ids) | is.na(response_names) | ids != response_names)
  if (length(mism) > 0L) {
    cier_abort(
      "cier_error_input",
      c("{.arg {arg}} row identifiers do not match the columns of \\
         {.arg responses}.",
        "x" = "Mismatch at {cli::qty(length(mism))}position{?s}: {.val {mism}}.",
        "i" = "{.arg {arg}}: {.val {ids[mism]}}; {.arg responses}: \\
               {.val {response_names[mism]}}.",
        "i" = "Order {.arg {arg}} so each row matches its {.arg responses} \\
               column (same order)."),
      data = list(arg = arg, observed = mism), call = call
    )
  }
  invisible(NULL)
}

# Validate the `items` frame the split-half family uses (even-odd, personal reliability),
# aligned to `responses` columns. Returns a normalized list(scale, reverse_keyed, min,
# max); `max` is NULL when its column is absent, `min` defaults to all-1.
check_items <- function(items, n_items, min_scales = 2L,
                        arg = "items", call = rlang::caller_env(),
                        response_names = NULL) {
  check_items_frame(items, n_items, arg, call, response_names)
  scale <- check_items_scale(items, min_scales, arg, call)
  reverse_keyed <- check_items_reverse(items, n_items, arg, call)
  minimum <- check_items_min(items, reverse_keyed, n_items, arg, call)
  maximum <- check_items_max(items, reverse_keyed, minimum, arg, call)
  list(scale = scale, reverse_keyed = reverse_keyed,
       min = minimum, max = maximum)
}

# Validate the `items` frame cier_ht() uses. Ht needs metadata only to reverse-score
# keyed items: the kernel accepts mixed category counts and never reads `max` directly, so
# ranges need not be homogeneous, `max` is required only on reverse-keyed items, and
# `scale` is not required. Returns a normalized list(reverse_keyed, min, max).
check_items_ht <- function(items, n_items, arg = "items",
                           call = rlang::caller_env(),
                           response_names = NULL) {
  check_items_frame(items, n_items, arg, call, response_names)
  reverse_keyed <- check_items_reverse(items, n_items, arg, call)
  minimum <- check_items_min(items, reverse_keyed, n_items, arg, call)
  maximum <- check_items_max(items, reverse_keyed, minimum, arg, call)
  list(reverse_keyed = reverse_keyed, min = minimum, max = maximum)
}

# `max` for the person-fit (Gnormed) bridge: required on EVERY item (finite whole,
# >= min + 1), not only reverse-keyed ones, because reverse-keying, per-item zero-basing,
# and the category count Ncat = max - min + 1 all read it. Per-item validity is a plain
# input error; span homogeneity is classified separately in
# check_items_span_homogeneous(). Returns the per-item vector.
check_items_max_personfit <- function(items, mins, arg, call) {
  maxs <- items$max
  # is_finite_whole() returns FALSE for an absent (NULL) column, so the missing
  # `max` case is caught here without a separate is.null() guard.
  ok <- is_finite_whole(maxs) && all(maxs >= mins + 1)
  if (!ok) {
    cier_abort(
      "cier_error_input",
      c("{.arg {arg}} needs an integer {.field max} of at least \\
         {.field min} + 1 on every item.",
        "i" = "{.field max} is the largest response option; {.field min} (the \\
               smallest, default 1) is read alongside it."),
      data = list(arg = arg), call = call
    )
  }
  maxs
}

# The span homogeneity Gnormed's single-Ncat closed form requires: one category count
# Ncat = max - min + 1 across all items (bases may differ: 1..5 and 0..4 both have five
# options). A heterogeneous span on otherwise-valid metadata is accurate metadata for
# genuinely mixed-format data Gnormed cannot score, not a malformed frame -- so the abort
# carries cier_error_backend_limit and cier_screen() records the index as
# skipped-with-reason. Takes the validated `maxs`, so a malformed max cannot be
# misclassified as a backend limit. Returns the shared Ncat.
check_items_span_homogeneous <- function(maxs, mins, arg, call) {
  spans <- maxs - mins
  if (length(unique(spans)) != 1L) {
    cier_abort(
      c("cier_error_backend_limit", "cier_error_input"),
      c("Gnormed cannot score items with different numbers of response \\
         categories.",
        "x" = "Observed numbers of options ({.field max} - {.field min} + 1): \\
               {.val {sort(unique(spans)) + 1}}.",
        "i" = "Gnormed scores one Ncat across all items; screen homogeneous item \\
               subsets separately, or use another index."),
      data = list(arg = arg, observed = sort(unique(spans)) + 1,
                  reason = "mixed response-category counts (Gnormed scores a single Ncat)"),
      call = call
    )
  }
  as.integer(spans[[1L]] + 1L)
}

# `min` for the person-fit backends: the scale base, used to reverse-key and zero-base
# every item (unlike the split-half family, which keys only reverse items), so it must be
# a finite whole number on every item. Optional; defaults to all-1.
check_items_min_personfit <- function(items, n_items, arg, call) {
  mins <- items$min
  if (is.null(mins)) {
    return(rep(1L, n_items))
  }
  ok <- is_finite_whole(mins)
  if (!ok) {
    cier_abort(
      "cier_error_input",
      c("{.field min} (the scale base) must be a finite integer on every item.",
        "i" = "{.field min} is the smallest response option (default 1); set it \\
               for a 0-based or bipolar scale."),
      data = list(arg = arg), call = call
    )
  }
  mins
}

# Validate the `items` frame the Gnormed bridge uses. Unlike the split-half family this
# needs `max` on every item and a homogeneous span (single-Ncat), and does not use
# `scale`. Returns a normalized list(reverse_keyed, min, max, ncat); the
# heterogeneous-span case also carries cier_error_backend_limit (see
# check_items_span_homogeneous). `ncat` is the shared category count, computed once here.
check_items_personfit <- function(items, n_items, arg = "items",
                                  call = rlang::caller_env(),
                                  response_names = NULL) {
  check_items_frame(items, n_items, arg, call, response_names)
  reverse_keyed <- check_items_reverse(items, n_items, arg, call)
  minimum <- check_items_min_personfit(items, n_items, arg, call)
  maximum <- check_items_max_personfit(items, minimum, arg, call)
  ncat <- check_items_span_homogeneous(maximum, minimum, arg, call)
  list(reverse_keyed = reverse_keyed, min = minimum, max = maximum, ncat = ncat)
}
