# Typed conditions for cier. Each carries a class chain with a parent
# (cier_error / cier_warning / cier_message) so callers catch by parent or
# sub-class; `data` is a named list; formatting goes through cli.

# ---- Internal helpers -------------------------------------------------------

# Class chain: sub-class(es), parent, base classes.
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

# Guard the shared `class`/`data` arguments. Base-R stop() by design: a typed
# cier condition cannot validate the inputs used to build one.
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
#' Wraps [cli::cli_abort()] with a class-stable error; callers catch by family
#' (`cier_error`) or sub-class (e.g. `cier_error_input`).
#'
#' @param class Character vector of sub-class(es); parent `cier_error` and the
#'   base condition classes are appended automatically.
#' @param message Character vector, forwarded to [cli::cli_abort()].
#' @param ... Forwarded to [cli::cli_abort()].
#' @param data Optional named list payload attached to the condition.
#' @param call Calling environment for the error trace.
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
#' Wraps [cli::cli_warn()] with a class-stable warning (parent `cier_warning`).
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
#' Wraps [cli::cli_inform()] with a class-stable message (parent `cier_message`).
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
#' Centralises extraction of the named list in a condition's `data` slot.
#'
#' @param cond A condition object.
#' @return The structured payload, or an empty list if absent.
#' @keywords internal
#' @noRd
cier_condition_data <- function(cond) {
  data <- cond$data
  if (is.null(data)) list() else data
}
