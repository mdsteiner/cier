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
    "companion_methods", "backend", "screenable", "notes")
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

# Purpose: Friendly internal accessor for the method-properties registry
#          (the single, citable source of each index's default cutoff, flag
#          direction, backend, and source paper). Internal for now; a curated
#          public accessor can be added at the docs/release slice if needed.
# Args:    None.
# Returns: A cier_method_info data.frame, one row per index.
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

# Purpose: Return the family of one method.
# Args:    method - a registry method id (character scalar).
# Returns: A length-1 character family label.
cier_method_family <- function(method) {
  cier_method_row(method)$family
}
