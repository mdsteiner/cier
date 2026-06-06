# Purpose: Minimal user-input checks used by the public wrappers. Raise
#          `cier_error_input` (a bad user argument) rather than the
#          `cier_error_data` / `cier_error_state` raised by the schema
#          validators.
# Args:    The value under inspection plus its argument name.
# Returns: The input (invisibly) on success.
# Invariants: A failure always raises a typed cier condition.

check_string <- function(x, arg, call = rlang::caller_env()) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    cier_abort(
      "cier_error_input",
      "{.arg {arg}} must be a non-empty character string.",
      data = list(arg = arg), call = call
    )
  }
  invisible(x)
}

check_choice <- function(x, arg, choices, call = rlang::caller_env()) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !(x %in% choices)) {
    cier_abort(
      "cier_error_input",
      "{.arg {arg}} must be one of {.val {choices}}.",
      data = list(arg = arg, observed = x), call = call
    )
  }
  invisible(x)
}

is_number_in_range <- function(x, lower, upper) {
  is.numeric(x) && length(x) == 1L && is.finite(x) && x >= lower && x <= upper
}

check_number <- function(x, arg, lower = -Inf, upper = Inf,
                         call = rlang::caller_env()) {
  if (!is_number_in_range(x, lower, upper)) {
    cier_abort(
      "cier_error_input",
      "{.arg {arg}} must be a single finite number in [{lower}, {upper}].",
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
