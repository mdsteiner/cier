# Slow-test tiering helper. Slow tests (e.g. external-package cross-checks that
# are heavy to install or run, or large-n paths) are opt-in: a normal local or
# CRAN run skips them, and CI requests them by setting CIER_SLOW_TESTS=true.
# Pair it with skip_on_cran() on any slow test so the tier is enforced in both
# directions. The opt-in value is "true" (case-insensitive); anything else, or an
# unset variable, skips.

skip_if_slow <- function() {
  if (!identical(tolower(Sys.getenv("CIER_SLOW_TESTS")), "true")) {
    testthat::skip("slow test; set CIER_SLOW_TESTS=true to run")
  }
  invisible(TRUE)
}
