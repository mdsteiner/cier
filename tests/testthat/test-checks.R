# Tests for the function-first input helper check_responses(), which lands with
# the first index. It coerces a data.frame / tibble to a numeric matrix and
# rejects anything that is not numeric-and-finite-or-NA. All failures are typed
# cier_error_input (asserted by class, never by message text).

test_that("a numeric matrix passes through unchanged in value", {
  m <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2L)
  out <- check_responses(m)
  expect_true(is.matrix(out) && is.numeric(out))
  expect_equal(unname(out), matrix(c(1, 2, 3, 4, 5, 6), nrow = 2L))
})

test_that("a data.frame of numeric columns is coerced, values preserved", {
  df <- data.frame(a = c(1, 2), b = c(3, 4), c = c(5, 6))
  out <- check_responses(df)
  expect_true(is.matrix(out) && is.numeric(out))
  expect_equal(unname(out), matrix(c(1, 2, 3, 4, 5, 6), nrow = 2L))
})

test_that("NA is an allowed missing marker", {
  m <- matrix(c(1, NA, 3, 4), nrow = 2L)
  expect_no_error(check_responses(m))
})

test_that("NaN is rejected", {
  m <- matrix(c(1, NaN, 3, 4), nrow = 2L)
  expect_error(check_responses(m), class = "cier_error_input")
})

test_that("Inf and -Inf are rejected", {
  expect_error(check_responses(matrix(c(1, Inf, 3, 4), nrow = 2L)),
               class = "cier_error_input")
  expect_error(check_responses(matrix(c(1, -Inf, 3, 4), nrow = 2L)),
               class = "cier_error_input")
})

test_that("a non-numeric (character) column is rejected", {
  df <- data.frame(a = c("x", "y"), b = c(1, 2), stringsAsFactors = FALSE)
  expect_error(check_responses(df), class = "cier_error_input")
})

test_that("a factor column is rejected", {
  df <- data.frame(a = factor(c("x", "y")), b = c(1, 2))
  expect_error(check_responses(df), class = "cier_error_input")
})

test_that("a logical matrix is rejected (not numeric)", {
  expect_error(check_responses(matrix(c(TRUE, FALSE, TRUE, FALSE), nrow = 2L)),
               class = "cier_error_input")
})

test_that("zero-row and zero-column inputs are rejected", {
  expect_error(check_responses(matrix(numeric(0), nrow = 0L, ncol = 3L)),
               class = "cier_error_input")
  expect_error(check_responses(matrix(numeric(0), nrow = 3L, ncol = 0L)),
               class = "cier_error_input")
})

test_that("a bare numeric vector is rejected, not silently shaped into a column", {
  # as.matrix(c(3,3,3,3)) would yield a 4x1 matrix (4 'respondents' of 1 item);
  # a vector is ambiguous, so require an explicit 2-D matrix / data.frame.
  expect_error(check_responses(c(3, 3, 3, 3)), class = "cier_error_input")
})

test_that("a higher-dimensional array is rejected", {
  expect_error(check_responses(array(1:24, dim = c(2L, 3L, 4L))),
               class = "cier_error_input")
})
