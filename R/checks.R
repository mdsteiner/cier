# Purpose: Minimal user-input checks used by the public wrappers. Raise
#          `cier_error_input` (a bad user argument) rather than the
#          `cier_error_data` / `cier_error_state` raised by the schema
#          validators.
# Args:    The value under inspection plus its argument name.
# Returns: The input (invisibly) on success.
# Invariants: A failure always raises a typed cier condition.

# The check_* helpers wrap checkmate's predicate functions (which return TRUE or
# a message string) so the diagnostics stay typed cier_error_input conditions
# with the package's cli phrasing, rather than checkmate's plain assert errors.

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

# Guard the two cutoff-override knobs an index exposes: a rate (`fpr` / `alpha` /
# `frac`) and a literal `cutoff`. They are mutually exclusive -- two ways to set
# the same cutoff -- so accepting both would be ambiguous. Abort when both are
# supplied. (Message stays generic so it reads correctly for every rate name.)
assert_single_override <- function(rate, rate_name, cutoff,
                                   call = rlang::caller_env()) {
  if (!is.null(rate) && !is.null(cutoff)) {
    cier_abort(
      "cier_error_input",
      c("Supply only one of {.arg {rate_name}} and {.arg cutoff}.",
        "i" = "They are two ways to set the same cutoff."),
      data = list(args = c(rate_name, "cutoff")), call = call
    )
  }
  invisible(NULL)
}

check_number <- function(x, arg, lower = -Inf, upper = Inf,
                         call = rlang::caller_env()) {
  if (!isTRUE(checkmate::check_number(x, lower = lower, upper = upper,
                                      finite = TRUE))) {
    cier_abort(
      "cier_error_input",
      "{.arg {arg}} must be a single finite number in [{lower}, {upper}].",
      data = list(arg = arg, observed = x), call = call
    )
  }
  invisible(x)
}

# A fraction of the item count: a single finite number in the half-open
# interval (0, 1] (0 would flag a zero-length run; values above 1 exceed the
# item count). checkmate guards the closed [0, 1]; `&& x > 0` opens the lower end.
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

# Coerce a user response payload to a validated numeric matrix. A data.frame or
# tibble is accepted and coerced (so users never call as.matrix()); `NA` is the
# only allowed missing marker -- NaN and infinities are rejected. This is the
# function-first input contract shared by every index wrapper.
check_responses <- function(responses, arg = "responses",
                            call = rlang::caller_env()) {
  if (!is.matrix(responses) && !is.data.frame(responses)) {
    cier_abort(
      "cier_error_input",
      c("{.arg {arg}} must be a matrix or data.frame (one row per respondent).",
        "x" = "Got {.cls {class(responses)}}; a bare vector or array is rejected."),
      data = list(arg = arg), call = call
    )
  }
  m <- as.matrix(responses)
  if (!is.numeric(m)) {
    cier_abort(
      "cier_error_input",
      "{.arg {arg}} must be numeric (a matrix or data.frame of numbers).",
      data = list(arg = arg), call = call
    )
  }
  if (nrow(m) == 0L || ncol(m) == 0L) {
    cier_abort(
      "cier_error_input",
      "{.arg {arg}} must have at least one row and one column.",
      data = list(arg = arg, nrow = nrow(m), ncol = ncol(m)), call = call
    )
  }
  if (any(is.nan(m) | is.infinite(m))) {
    cier_abort(
      "cier_error_input",
      c("{.arg {arg}} may contain only finite numbers or {.val NA}.",
        "x" = "Found {.val NaN} or infinite values."),
      data = list(arg = arg), call = call
    )
  }
  m
}
