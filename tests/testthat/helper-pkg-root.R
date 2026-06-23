# Shared helpers for tests that need the package source tree on disk
# (the developer-vocabulary guard, the roxygen-up-to-date check, the
# references-DOI guard). Loaded automatically by testthat because of the
# `helper-` prefix.

# Resolve the package root from the test location. `testthat::test_path()`
# returns the test directory; two `..` jumps land in the package root.
cier_pkg_root <- function() {
  normalizePath(testthat::test_path("..", ".."), mustWork = FALSE)
}

# Skip the calling test when the package source tree is not present (the
# typical case under `R CMD check`, where only the installed package is
# available rather than the development tree).
skip_if_no_description <- function(reason = NULL) {
  desc <- file.path(cier_pkg_root(), "DESCRIPTION")
  if (!file.exists(desc)) {
    testthat::skip(reason %||% "DESCRIPTION not on disk (installed package)")
  }
  invisible(TRUE)
}

# Fallback for environments where `rlang::%||%` is not in scope.
if (!exists("%||%", mode = "function", inherits = TRUE)) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}
