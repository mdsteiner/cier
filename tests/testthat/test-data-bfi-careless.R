# Tests for the bundled bfi_careless example dataset.
#
# Trust model: bfi_careless is shipped data, not a computed statistic, so there
# is no oracle or cross-package parity. Instead this file pins the object's
# contract -- the exact shape, names, reverse-key structure, scales, and per-
# region value bounds that R/data.R documents and that every help-page example
# and the getting-started vignette rely on. A silent rebuild that transposed the
# matrix, relabelled or zero-based a column, dropped the reverse-keyed items,
# flattened a column to a constant, or let an out-of-range / NA / NaN value in
# must fail here.

# The 44 BFI item columns are the first 44; the two independent attention-check
# labels are the last two, in this order.
bfi_cols <- function() seq_len(44L)
bfi_names <- function() names(bfi_careless)[bfi_cols()]

# ---- Existence + shape ------------------------------------------------------

test_that("bfi_careless is a 394 x 46 data frame", {
  expect_s3_class(bfi_careless, "data.frame")
  expect_identical(dim(bfi_careless), c(394L, 46L))
})

# ---- Column names + structure ----------------------------------------------

test_that("the first 44 columns are the BFI-44 items and the last two are the labels", {
  nm <- names(bfi_careless)
  # Exactly the first 44 are BFI items; none of the others are.
  expect_true(all(grepl("^v_BFI_", nm[bfi_cols()])))
  expect_false(any(grepl("^v_BFI_", nm[45:46])))
  # The two attention-check labels are the final two columns, in this order.
  expect_identical(nm[45:46], c("v_Bogus_Item", "v_IRI"))
})

test_that("exactly 16 BFI items are reverse-keyed (the _R suffix)", {
  expect_identical(sum(grepl("_R$", bfi_names())), 16L)
})

test_that("the five Big-Five scales are present and are the only scales", {
  scales <- sub("^v_BFI_([A-Za-z]+)[0-9].*$", "\\1", bfi_names())
  expect_identical(sort(unique(scales)), c("AG", "CON", "EX", "NEU", "OP"))
  # Every BFI name resolves to one of the five scales (no name slips the regex).
  expect_true(all(scales %in% c("AG", "CON", "EX", "NEU", "OP")))
})

# ---- Value integrity --------------------------------------------------------

# Whole-number-valued (robust to integer-vs-double storage after a rebuild),
# numeric, no NA in any of the three answered regions, and no NaN / infinity
# anywhere.
expect_whole_in_range <- function(x, lo, hi) {
  expect_true(is.numeric(x))
  expect_false(anyNA(x))
  expect_true(all(x == as.integer(x)))
  expect_gte(min(x), lo)
  expect_lte(max(x), hi)
}

test_that("the 44 BFI items are whole numbers in 1:5 with no missing values", {
  bfi <- as.matrix(bfi_careless[, bfi_cols()])
  expect_whole_in_range(bfi, 1L, 5L)
})

test_that("v_Bogus_Item is a whole number in 1:5 with no missing values", {
  expect_whole_in_range(bfi_careless$v_Bogus_Item, 1L, 5L)
})

test_that("v_IRI is a whole number in 0:7 with no missing values", {
  expect_whole_in_range(bfi_careless$v_IRI, 0L, 7L)
})

test_that("no column contains NaN or infinite values", {
  bad <- vapply(bfi_careless, function(col) any(is.nan(col) | is.infinite(col)),
                logical(1L))
  expect_false(any(bad))
})

# ---- No degenerate columns --------------------------------------------------

# Real survey responses vary; a column flattened to a single value (a plausible
# rebuild slip) would pass the range and missingness checks but is not the
# bundled data. Pin that every column carries more than one distinct response.
test_that("every column has at least two distinct values", {
  n_distinct <- vapply(bfi_careless,
                       function(col) length(unique(col)), integer(1L))
  expect_true(all(n_distinct >= 2L))
})
