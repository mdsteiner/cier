# Purpose: Typed condition system for cier (errors, warnings, messages).
# Args:    See per-function documentation below.
# Returns: cier_abort() never returns; cier_warn()/cier_inform() return
#          invisibly.
# Invariants:
#   - Every condition raised through these helpers carries a class chain
#     with a documented parent (cier_error / cier_warning / cier_message)
#     so downstream callers can catch by parent or by sub-class.
#   - `data` is always a (possibly empty) named list attached to the
#     condition.
#   - All formatting goes through `cli`.

# ---- Internal helpers -------------------------------------------------------

# Purpose: Build the full S3 class chain for a condition.
# Args:    class — sub-class(es); kind — "error", "warning", or "message".
# Returns: A character vector: sub-class(es), parent, base classes.
cier_condition_classes <- function(class, kind) {
  parent <- switch(kind,
    error   = "cier_error",
    warning = "cier_warning",
    message = "cier_message"
  )
  base <- switch(kind,
    error   = c("error", "condition"),
    warning = c("warning", "condition"),
    message = c("message", "condition")
  )
  c(unique(c(class, parent)), base)
}

# Purpose: Guard the `class` / `data` arguments shared by the three
#          condition helpers. These are base-R stops by design: a typed
#          cier condition cannot validate the inputs used to build a typed
#          cier condition.
# Args:    class, data — the helper arguments.
# Returns: Invisibly NULL; stops on a violation.
check_condition_args <- function(class, data) {
  if (!is.character(class) || length(class) == 0L || anyNA(class)) {
    stop("`class` must be a non-empty character vector with no NAs.",
         call. = FALSE)
  }
  if (!is.list(data) || (length(data) > 0L && is.null(names(data)))) {
    stop("`data` must be a (possibly empty) named list.", call. = FALSE)
  }
  invisible(NULL)
}

# ---- Package-internal API ---------------------------------------------------

#' Abort with a typed cier condition
#'
#' Wraps [cli::cli_abort()] with a class-stable error. Every typed error in
#' cier is raised through this helper so callers can catch by family
#' (`cier_error`) or by sub-class (e.g. `cier_error_input`).
#'
#' @param class Character vector. Sub-class(es) to attach; the parent
#'   `cier_error` and the base condition classes are appended automatically.
#' @param message Character vector. Forwarded to [cli::cli_abort()].
#' @param ... Forwarded to [cli::cli_abort()].
#' @param data Named list. Optional payload attached to the condition.
#' @param call Calling environment used for the error trace.
#' @param .envir Environment for `cli` inline-markup interpolation.
#' @return Never returns; raises a condition.
#' @keywords internal
#' @noRd
cier_abort <- function(class, message, ..., data = list(),
                       call = rlang::caller_env(),
                       .envir = rlang::caller_env()) {
  check_condition_args(class, data)
  cli::cli_abort(
    message = message, ...,
    class = cier_condition_classes(class, "error"),
    data = data, call = call, .envir = .envir
  )
}

#' Warn with a typed cier condition
#'
#' Wraps [cli::cli_warn()] with a class-stable warning (parent
#' `cier_warning`).
#'
#' @inheritParams cier_abort
#' @return Invisibly `NULL`; signals a warning.
#' @keywords internal
#' @noRd
cier_warn <- function(class, message, ..., data = list(),
                      call = rlang::caller_env(),
                      .envir = rlang::caller_env()) {
  check_condition_args(class, data)
  cli::cli_warn(
    message = message, ...,
    class = cier_condition_classes(class, "warning"),
    data = data, call = call, .envir = .envir
  )
  invisible(NULL)
}

#' Inform with a typed cier condition
#'
#' Wraps [cli::cli_inform()] with a class-stable message (parent
#' `cier_message`).
#'
#' @inheritParams cier_abort
#' @return Invisibly `NULL`; signals a message.
#' @keywords internal
#' @noRd
cier_inform <- function(class, message, ..., data = list(),
                        .envir = rlang::caller_env()) {
  check_condition_args(class, data)
  cli::cli_inform(
    message = message, ...,
    class = cier_condition_classes(class, "message"),
    data = data, .envir = .envir
  )
  invisible(NULL)
}

#' Pull the structured payload out of a cier condition
#'
#' Conditions raised through the helpers above carry a named list in the
#' `data` slot. This centralises extraction (production reads it in the
#' screen's backend-limit handler; tests read it to assert payloads).
#'
#' @param cond A condition object.
#' @return The structured payload, or an empty list if absent.
#' @keywords internal
#' @noRd
cier_condition_data <- function(cond) {
  data <- cond$data
  if (is.null(data)) list() else data
}
