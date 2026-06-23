# Scalar / argument validators and cutoff-override guards for the public wrappers.
# Each raises `cier_error_input` on a bad argument and returns its input invisibly on
# success. Payload coercers live in R/checks-data.R; `items`-frame validators in
# R/checks-items.R. The check_* helpers wrap checkmate's predicates so diagnostics stay
# typed cier_error_input with the package's cli phrasing.

check_string <- function(x, arg, call = rlang::caller_env()) {
  if (!isTRUE(checkmate::check_string(x, min.chars = 1L))) {
    cier_abort(
      "cier_error_input",
      "{.arg {arg}} must be a non-empty character string.",
      data = list(arg = arg), call = call
    )
  }
  invisible(x)
}

check_choice <- function(x, arg, choices, call = rlang::caller_env()) {
  if (!isTRUE(checkmate::check_choice(x, choices))) {
    cier_abort(
      "cier_error_input",
      "{.arg {arg}} must be one of {.val {choices}}.",
      data = list(arg = arg, observed = x), call = call
    )
  }
  invisible(x)
}

# Abort when more than one mutually-exclusive cutoff knob is supplied. `args` is a named
# list of knob values (NULL where unset); the message names the offenders.
assert_single_cutoff <- function(args, call = rlang::caller_env()) {
  supplied <- names(args)[!vapply(args, is.null, logical(1L))]
  if (length(supplied) > 1L) {
    cier_abort(
      "cier_error_input",
      c("Supply only one of {.arg {supplied}}.",
        "i" = "They are different ways to set the same cutoff."),
      data = list(args = supplied), call = call
    )
  }
  invisible(NULL)
}

check_number <- function(x, arg, lower = -Inf, upper = Inf, whole = FALSE,
                         call = rlang::caller_env()) {
  # whole = TRUE requires an integer (count cutoffs of cier_page_time /
  # cier_attention), rejecting a fractional 2.5.
  ok <- if (whole) {
    isTRUE(checkmate::check_int(x, lower = lower, upper = upper))
  } else {
    isTRUE(checkmate::check_number(x, lower = lower, upper = upper, finite = TRUE))
  }
  if (!ok) {
    unit <- if (whole) "whole number" else "number"
    cier_abort(
      "cier_error_input",
      "{.arg {arg}} must be a single finite {unit} in [{lower}, {upper}].",
      data = list(arg = arg, observed = x), call = call
    )
  }
  invisible(x)
}

# A tail probability / significance level (`alpha`, `fpr`): a finite number strictly
# inside (0, 1). The `&& x > 0 && x < 1` opens both ends checkmate's closed [0, 1]
# leaves, rejecting 0 (flags nobody) and 1 (everybody).
check_open_unit <- function(x, arg, call = rlang::caller_env()) {
  is_open_unit <- isTRUE(
    checkmate::check_number(x, lower = 0, upper = 1, finite = TRUE)
  ) && x > 0 && x < 1
  if (!is_open_unit) {
    cier_abort(
      "cier_error_input",
      "{.arg {arg}} must be a single number in the open interval (0, 1).",
      data = list(arg = arg, observed = x), call = call
    )
  }
  invisible(x)
}

# A fraction of the item count: a finite number in (0, 1] (0 would flag a zero-length
# run; >1 exceeds the item count).
check_fraction <- function(x, arg, call = rlang::caller_env()) {
  is_fraction <- isTRUE(
    checkmate::check_number(x, lower = 0, upper = 1, finite = TRUE)
  ) && x > 0
  if (!is_fraction) {
    cier_abort(
      "cier_error_input",
      "{.arg {arg}} must be a single number in the interval (0, 1].",
      data = list(arg = arg, observed = x), call = call
    )
  }
  invisible(x)
}

# A single non-missing TRUE/FALSE, for boolean wrapper switches (e.g. `resample`).
check_flag <- function(x, arg, call = rlang::caller_env()) {
  if (!isTRUE(checkmate::check_flag(x))) {
    cier_abort(
      "cier_error_input",
      "{.arg {arg}} must be a single {.val TRUE} or {.val FALSE}.",
      data = list(arg = arg), call = call
    )
  }
  invisible(x)
}

# A positive whole number (>= 1, no NA), for the resampling iteration count
# (`n_resamples`). Integerish tolerance rejects a non-whole 2.5.
check_count <- function(x, arg, call = rlang::caller_env()) {
  if (!isTRUE(checkmate::check_count(x, positive = TRUE))) {
    cier_abort(
      "cier_error_input",
      "{.arg {arg}} must be a single positive whole number.",
      data = list(arg = arg, observed = x), call = call
    )
  }
  invisible(x)
}

# A whole number >= 0, for a count that may legitimately be zero (`n_checks`, where a
# simulation may include no attention checks). Shared by cier_simulate() and
# sim_direct_checks().
check_count_nonneg <- function(x, arg, call = rlang::caller_env()) {
  if (!isTRUE(checkmate::check_count(x))) {
    cier_abort(
      "cier_error_input",
      "{.arg {arg}} must be a single non-negative whole number.",
      data = list(arg = arg, observed = x), call = call
    )
  }
  invisible(x)
}

# A single integer (integerish, finite, any sign, no NA), for an RNG `seed`: set.seed()
# would silently truncate a fractional value.
check_int <- function(x, arg, call = rlang::caller_env()) {
  if (!isTRUE(checkmate::check_int(x))) {
    cier_abort(
      "cier_error_input",
      "{.arg {arg}} must be a single integer.",
      data = list(arg = arg, observed = x), call = call
    )
  }
  invisible(x)
}

# `nms` non-NULL, NA-free, all non-empty, no duplicates. Shared front-end for the
# simulator's keyed override lists (`timing`, `pattern_params`).
check_unique_names <- function(nms, arg, call = rlang::caller_env()) {
  named_ok <- !is.null(nms) && !anyNA(nms) && all(nzchar(nms)) &&
    anyDuplicated(nms) == 0L
  if (!named_ok) {
    cier_abort(
      "cier_error_input",
      "{.arg {arg}} entries must be uniquely named.",
      data = list(arg = arg), call = call
    )
  }
  invisible(nms)
}

# Reject any of `nms` not in `allowed` so a typo cannot silently fall back to a default;
# the error lists offenders and permitted names. `label` is the message noun. Shared
# front-end for the keyed override lists. No-op when all allowed.
reject_unknown_keys <- function(nms, allowed, arg, label,
                                call = rlang::caller_env()) {
  extra <- setdiff(nms, allowed)
  if (length(extra) == 0L) {
    return(invisible(NULL))
  }
  cier_abort(
    "cier_error_input",
    c("Unknown {label}{?s}: {.val {extra}}.",
      "i" = if (length(allowed) == 0L) {
        "None is accepted here."
      } else {
        "Allowed: {.val {allowed}}."
      }),
    data = list(arg = arg, observed = extra), call = call
  )
}

# Validate the two mutually-exclusive cutoff-override knobs every percentile index
# exposes: `fpr` (tail mass in the open unit interval) and a literal `cutoff` on the
# score (a finite number in `[lower, upper]`). Run before the kernel to fail early.
check_percentile_overrides <- function(fpr, cutoff, lower = -Inf, upper = Inf,
                                       call = rlang::caller_env()) {
  if (!is.null(fpr)) check_open_unit(fpr, "fpr", call = call)
  if (!is.null(cutoff)) {
    check_number(cutoff, "cutoff", lower = lower, upper = upper, call = call)
  }
  assert_single_cutoff(list(fpr = fpr, cutoff = cutoff), call = call)
  invisible(NULL)
}

# Predicate: a numeric vector of finite whole numbers (no NA / NaN / Inf, every element
# integer-valued). Shared core of the items-column validators and the data-side count /
# integer checks. `is.numeric()` is first to short-circuit a non-numeric input.
is_finite_whole <- function(v) {
  is.numeric(v) && all(is.finite(v)) && all(v == round(v))
}
