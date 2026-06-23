# Unit test for the slow-test tiering helper skip_if_slow(). The helper gates
# slow tests on the CIER_SLOW_TESTS environment variable (default off): a normal
# local / CRAN run skips them, and CI opts in by setting CIER_SLOW_TESTS=true.
# Both branches are asserted here by catching the `skip` condition the helper
# raises, so the assertions run rather than skipping this test itself.

skipped_or_value <- function() {
  tryCatch(skip_if_slow(), skip = function(cnd) "skipped")
}

test_that("skip_if_slow() skips when CIER_SLOW_TESTS is unset", {
  withr::with_envvar(c(CIER_SLOW_TESTS = NA), {
    expect_identical(skipped_or_value(), "skipped")
  })
})

test_that("skip_if_slow() skips when CIER_SLOW_TESTS is a non-'true' value", {
  withr::with_envvar(c(CIER_SLOW_TESTS = "false"), {
    expect_identical(skipped_or_value(), "skipped")
  })
  withr::with_envvar(c(CIER_SLOW_TESTS = "1"), {
    expect_identical(skipped_or_value(), "skipped")
  })
})

test_that("skip_if_slow() runs (does not skip) when CIER_SLOW_TESTS=true", {
  withr::with_envvar(c(CIER_SLOW_TESTS = "true"), {
    expect_true(isTRUE(skipped_or_value()))
  })
})

test_that("skip_if_slow() opt-in is case-insensitive", {
  withr::with_envvar(c(CIER_SLOW_TESTS = "TRUE"), {
    expect_true(isTRUE(skipped_or_value()))
  })
})

test_that("skip_if_slow() returns its TRUE value invisibly when it runs", {
  withr::with_envvar(c(CIER_SLOW_TESTS = "true"), {
    expect_invisible(skip_if_slow())
  })
})
