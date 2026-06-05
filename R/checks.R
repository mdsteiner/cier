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
