# Purpose: Method-properties registry: loader, accessors, and validator.
#          The CSV (inst/extdata/method-properties.csv) is the single source
#          of truth for each index's cutoff default, flag direction, optional
#          backend package, source citation, and screen membership.
# Args:    See documentation.
# Returns: A cier_method_info data.frame.
# Invariants: Cached in a package-private environment; validated on load.

.cier_registry_cache <- new.env(parent = emptyenv())

cier_registry_path <- function() {
  path <- system.file("extdata", "method-properties.csv",
                      package = "cier", mustWork = FALSE)
  if (nzchar(path)) {
    return(path)
  }
  # Fallback during devtools::load_all before the package is installed.
  candidate <- file.path("inst", "extdata", "method-properties.csv")
  if (file.exists(candidate)) candidate else ""
}

read_registry_csv <- function(path, call = rlang::caller_env()) {
  if (!nzchar(path) || !file.exists(path)) {
    cier_abort(
      "cier_error_state",
      c("Method-properties CSV not found.",
        "x" = "Looked for {.path inst/extdata/method-properties.csv}."),
      data = list(arg = "path", observed = path), call = call
    )
  }
  utils::read.csv(path, stringsAsFactors = FALSE, na.strings = c("", "NA"))
}

# Coerce the non-character columns to their schema types.
coerce_registry_columns <- function(df) {
  if ("paper_year" %in% names(df)) {
    df$paper_year <- as.integer(df$paper_year)
  }
  if ("default_cutoff_value" %in% names(df)) {
    df$default_cutoff_value <- as.numeric(df$default_cutoff_value)
  }
  if ("screenable" %in% names(df)) {
    df$screenable <- as.logical(df$screenable)
  }
  df
}

new_cier_method_info <- function(registry) {
  structure(registry, class = c("cier_method_info", "data.frame"))
}

# ---- registry validation ----------------------------------------------------

registry_columns <- function() {
  c("method", "family", "paper_year", "paper_citation_key", "doi",
    "default_cutoff_method", "default_cutoff_value", "flag_direction",
    "companion_methods", "backend", "screenable", "vote_group", "notes")
}

# Catch the realistic hand-edit mistake: a typo in a numeric/logical column
# that coercion turned into NA (e.g. a malformed default cutoff value).
validate_registry_values <- function(x) {
  for (col in c("default_cutoff_value", "paper_year", "screenable")) {
    if (anyNA(x[[col]])) {
      cier_abort("cier_error_data",
                 "{.field {col}} has missing values; check the registry CSV.",
                 data = list(arg = col))
    }
  }
  invisible(NULL)
}

validate_registry_vocab <- function(x) {
  if (anyNA(x$method) || any(!nzchar(x$method)) ||
        anyDuplicated(x$method) > 0L) {
    cier_abort("cier_error_data", "{.field method} must be unique, non-missing.")
  }
  if (anyNA(x$family) ||
        length(setdiff(x$family, cier_family_levels())) > 0L) {
    cier_abort("cier_error_data", "{.field family} must be in the vocabulary.")
  }
  if (anyNA(x$default_cutoff_method) ||
        length(setdiff(x$default_cutoff_method, cier_cutoff_methods())) > 0L) {
    cier_abort("cier_error_data", "{.field default_cutoff_method} unknown.")
  }
  if (anyNA(x$flag_direction) ||
        length(setdiff(x$flag_direction, cier_flag_directions())) > 0L) {
    cier_abort("cier_error_data", "{.field flag_direction} must be upper/lower.")
  }
  # `vote_group` drives the screen's vote collapse; a missing label would silently
  # drop the index from the agreement count, so every row must carry one.
  if (anyNA(x$vote_group) || any(!nzchar(x$vote_group))) {
    cier_abort("cier_error_data",
               "{.field vote_group} must be a non-missing label on every row.")
  }
  invisible(NULL)
}

validate_registry_companions <- function(x) {
  for (i in seq_len(nrow(x))) {
    comp <- x$companion_methods[[i]]
    if (is.na(comp) || !nzchar(comp)) {
      next
    }
    parts <- trimws(strsplit(comp, ",", fixed = TRUE)[[1L]])
    bad <- setdiff(parts, x$method)
    if (length(bad) > 0L) {
      cier_abort("cier_error_data",
                 "{.field companion_methods} must reference registry methods.",
                 data = list(arg = "companion_methods", observed = bad))
    }
  }
  invisible(NULL)
}

validate_cier_method_info <- function(x) {
  if (!inherits(x, "cier_method_info") || !is.data.frame(x)) {
    cier_abort("cier_error_state", "Expected a {.cls cier_method_info} frame.")
  }
  if (!identical(names(x), registry_columns())) {
    missing_cols <- setdiff(registry_columns(), names(x))
    extra_cols <- setdiff(names(x), registry_columns())
    cier_abort(
      "cier_error_data",
      c("{.cls cier_method_info} columns must match the frozen schema.",
        "x" = "Missing: {.field {missing_cols}}; extra: {.field {extra_cols}}.")
    )
  }
  validate_registry_vocab(x)
  validate_registry_values(x)
  validate_registry_companions(x)
  invisible(x)
}

# ---- loader + accessors -----------------------------------------------------

load_method_registry <- function() {
  call <- rlang::caller_env()
  if (!is.null(.cier_registry_cache$registry)) {
    return(.cier_registry_cache$registry)
  }
  raw <- read_registry_csv(cier_registry_path(), call)
  reg <- new_cier_method_info(coerce_registry_columns(raw))
  validate_cier_method_info(reg)
  .cier_registry_cache$registry <- reg
  reg
}

#' List the detection indices and their registry properties
#'
#' Returns the package's **method registry**: one row per detection index, with
#' the citable defaults every wrapper reads -- the default cutoff method and
#' value, the flag direction, the source paper, the optional backend package,
#' and the screen membership. Use it to see the valid `methods` ids for
#' [cier_screen()] and to trace where each index's defaults come from.
#'
#' @details
#' The registry is the single source of truth behind the index functions; this
#' accessor is read-only. Columns:
#' - `method`: the index function name (also the id for `cier_screen(methods=)`).
#' - `family`: `"indirect"` (response-pattern) or `"personfit"`.
#' - `paper_year`, `paper_citation_key`, `doi`: the source paper the default is
#'   traced to.
#' - `default_cutoff_method` / `default_cutoff_value`: how the default cutoff is
#'   resolved (`"percentile"` at a target tail mass, `"fixed"` fraction,
#'   `"chisq"` tail probability, or `"perfit_null"` Monte-Carlo level) and the
#'   cited default rate.
#' - `flag_direction`: which tail flags carelessness (`"upper"` / `"lower"`).
#' - `companion_methods`: indices commonly paired with this one.
#' - `backend`: the `Suggests` package that scores the statistic (`NA` for the
#'   pure indices).
#' - `screenable`: whether [cier_screen()] runs it.
#' - `vote_group`: the construct label the screen collapses correlated indices
#'   under (even-odd and personal reliability share `"consistency"`).
#' - `notes`: a one-line description.
#'
#' @return A data frame (class `cier_method_info`), one row per index, with the
#'   columns listed under Details.
#' @seealso [cier_screen()]; the index functions, e.g. [cier_longstring()],
#'   [cier_irv()], [cier_personal_reliability()]
#' @family orchestration
#' @export
#' @examples
#' cier_methods()[, c("method", "family", "default_cutoff_method",
#'                    "flag_direction", "screenable")]
cier_methods <- function() {
  load_method_registry()
}

# Purpose: Return the single registry row for one method.
# Args:    method - a registry method id (character scalar).
# Returns: A one-row cier_method_info data.frame.
cier_method_row <- function(method) {
  check_string(method, "method")
  reg <- load_method_registry()
  idx <- match(method, reg$method)
  if (is.na(idx)) {
    cier_abort(
      "cier_error_input",
      c("Method {.val {method}} not found in the registry.",
        "i" = "See {.fun cier_methods} for the available set."),
      data = list(arg = "method", observed = method)
    )
  }
  reg[idx, , drop = FALSE]
}

