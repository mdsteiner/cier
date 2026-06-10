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

# Tests for check_items(), which lands with the first metadata index (even-odd).
# It validates the per-item `items` frame the split-half family uses: `scale`
# (>= min_scales distinct), an optional logical `reverse_keyed` (defaults
# all-FALSE), and `max` (the largest response option) -- required (integer
# >= min + 1, non-NA) only on items that are actually reverse-keyed. All
# failures are typed cier_error_input.

test_that("check_items returns normalized scale / reverse_keyed / max", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, TRUE),
                   max = 5L)
  out <- check_items(it, n_items = 4L)
  expect_identical(out$scale, c("A", "A", "B", "B"))
  expect_identical(out$reverse_keyed, c(FALSE, TRUE, FALSE, TRUE))
  expect_identical(out$max, rep(5L, 4L))
})

test_that("check_items defaults reverse_keyed to all-FALSE", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L))
  out <- check_items(it, n_items = 4L)
  expect_identical(out$reverse_keyed, rep(FALSE, 4L))
})

test_that("check_items coerces a factor scale column to character", {
  it <- data.frame(scale = factor(rep(c("A", "B"), each = 2L)),
                   reverse_keyed = FALSE)
  out <- check_items(it, n_items = 4L)
  expect_type(out$scale, "character")
  expect_identical(out$scale, c("A", "A", "B", "B"))
})

test_that("check_items requires the items row-count to equal n_items", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L))
  expect_error(check_items(it, n_items = 5L), class = "cier_error_input")
})

test_that("check_items requires at least min_scales distinct scales", {
  it <- data.frame(scale = rep("A", 4L))
  expect_error(check_items(it, n_items = 4L, min_scales = 2L),
               class = "cier_error_input")
})

test_that("check_items requires a scale column", {
  it <- data.frame(reverse_keyed = rep(FALSE, 4L))
  expect_error(check_items(it, n_items = 4L), class = "cier_error_input")
})

test_that("check_items rejects a non-data-frame items argument", {
  expect_error(check_items(list(scale = c("A", "B")), n_items = 2L),
               class = "cier_error_input")
})

test_that("check_items rejects a non-logical reverse_keyed column", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(0L, 1L, 0L, 1L))
  expect_error(check_items(it, n_items = 4L), class = "cier_error_input")
})

test_that("check_items rejects NA in reverse_keyed", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, NA, FALSE, TRUE))
  expect_error(check_items(it, n_items = 4L), class = "cier_error_input")
})

test_that("check_items requires max when an item is reverse-keyed", {
  # No max column but a reverse item -> cannot reverse-score -> error.
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, FALSE))
  expect_error(check_items(it, n_items = 4L), class = "cier_error_input")
})

test_that("check_items rejects NA max on a reverse item", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, FALSE),
                   max = c(5L, NA, 5L, 5L))
  expect_error(check_items(it, n_items = 4L), class = "cier_error_input")
})

test_that("check_items allows NA max on non-reverse items", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, FALSE),
                   max = c(NA, 5L, NA, 3L))
  expect_no_error(check_items(it, n_items = 4L))
})

test_that("check_items rejects max at or below min on a reverse item", {
  # max == min is a one-option "scale" (nothing to reflect); with the default
  # min = 1, max = 1 is the malformed declaration the old categories >= 2
  # bound caught.
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, FALSE),
                   max = c(5L, 1L, 5L, 5L))
  expect_error(check_items(it, n_items = 4L), class = "cier_error_input")
  # Same violation with an explicit min ABOVE the absolute-2 line: max == min
  # == 2 must be rejected too, so a validator testing `max >= 2` (ignoring min)
  # is wrong.
  it2 <- data.frame(scale = rep(c("A", "B"), each = 2L),
                    reverse_keyed = c(FALSE, TRUE, FALSE, FALSE),
                    max = c(5L, 2L, 5L, 5L), min = c(1L, 2L, 1L, 1L))
  expect_error(check_items(it2, n_items = 4L), class = "cier_error_input")
})

test_that("check_items accepts a two-option scale (max == min + 1)", {
  # The smallest valid scale: two response options. With min = 0 this is a
  # 0/1 item; a validator demanding max >= min + 2 (or max >= 2 regardless of
  # min) would wrongly reject it.
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, TRUE),
                   max = 1L, min = 0L)
  expect_no_error(check_items(it, n_items = 4L))
})

test_that("check_items rejects non-integer max on a reverse item", {
  # The reflection (min + max) - x assumes whole-number response options; a
  # fractional value is a malformed item definition.
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, FALSE),
                   max = c(5, 2.5, 5, 5))
  expect_error(check_items(it, n_items = 4L), class = "cier_error_input")
})

test_that("check_items rejects a non-finite (Inf) max on a reverse item", {
  # Inf passes is.numeric / >= min + 1 / == round(); only the is.finite guard
  # rejects it, so without that guard the item would reflect to (min + Inf) - x
  # = Inf and silently poison the reverse columns instead of erroring.
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, FALSE),
                   max = c(5, Inf, 5, 5))
  expect_error(check_items(it, n_items = 4L), class = "cier_error_input")
})

test_that("check_items does not require max on a forward-keyed battery", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L), reverse_keyed = FALSE)
  expect_no_error(check_items(it, n_items = 4L))
})

# `min` -- the response-scale base. Optional; defaults to 1 (1..max coding).
# When supplied it generalises the reverse-keying reflection to (min + max) - x
# so a 0-based (or bipolar) scale reflects onto itself.

test_that("check_items defaults min to all-1 when the column is absent", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, TRUE), max = 5L)
  out <- check_items(it, n_items = 4L)
  expect_identical(out$min, rep(1L, 4L))
})

test_that("check_items returns a supplied min column", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, TRUE),
                   max = 5L, min = 0L)
  out <- check_items(it, n_items = 4L)
  expect_identical(out$min, rep(0L, 4L))
})

test_that("check_items allows a zero, negative, or bipolar min base on reverse items", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, TRUE),
                   max = 5L, min = c(0L, -2L, 0L, -2L))
  expect_no_error(check_items(it, n_items = 4L))
})

test_that("check_items rejects a non-finite min on a reverse item", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, FALSE),
                   max = 5L, min = c(1, Inf, 1, 1))
  expect_error(check_items(it, n_items = 4L), class = "cier_error_input")
})

test_that("check_items rejects a non-integer min on a reverse item", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, FALSE),
                   max = 5L, min = c(1, 0.5, 1, 1))
  expect_error(check_items(it, n_items = 4L), class = "cier_error_input")
})

test_that("check_items rejects NA min on a reverse item", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, FALSE),
                   max = 5L, min = c(1L, NA, 1L, 1L))
  expect_error(check_items(it, n_items = 4L), class = "cier_error_input")
})

test_that("check_items allows NA min on non-reverse items", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, FALSE),
                   max = 5L, min = c(NA, 0L, NA, 3L))
  expect_no_error(check_items(it, n_items = 4L))
})
