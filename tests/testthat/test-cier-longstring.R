# Tests for cier_longstring() (Slice 2 — the first index).
#
# Trust model: the independent oracle (ref_longstring) re-derives the max-run
# statistic with a hand-rolled loop and never calls the production kernel; the
# cross-package check pins bytewise parity with careless::longstring (the
# documented contract, tolerance 0); the property / mutant-killer block targets
# each mutant named in dev/restart/index-specs.md card 1.

source(test_path("..", "reference", "ref-longstring-johnson2005.R"))

# A self-contained fixture of present-data rows (no all-NA row, so the kernel
# value is well defined for every row and the oracle agrees row for row).
present_fixture <- function() {
  matrix(
    c(
      1, 2, 3, 4, 5, 1, 2, 3,   # cycling          -> longest 1
      3, 3, 3, 3, 3, 3, 3, 3,   # straightliner    -> longest 8
      1, 5, 1, 5, 1, 5, 1, 5,   # alternating      -> longest 1
      2, 2, 2, 2, NA, NA, NA, NA # 4 twos, 4 NAs   -> longest 4 (NAs break)
    ),
    nrow = 4L, byrow = TRUE
  )
}

# ---- Schema -----------------------------------------------------------------

test_that("cier_longstring returns a list-based cier_index with the pinned schema", {
  out <- cier_longstring(present_fixture())
  expect_s3_class(out, "cier_index")
  expect_type(out, "list")
  expect_identical(names(out),
                   c("value", "flagged", "method", "cutoff", "direction"))
  expect_type(out$value, "double")
  expect_type(out$flagged, "logical")
  expect_identical(length(out$value), 4L)
  expect_identical(out$method, "cier_longstring")
  expect_identical(out$direction, "upper")
  expect_type(out$cutoff, "double")
  # Concrete, non-tautological: present_fixture values are 1, 8, 1, 4; at the
  # default cutoff ceiling(0.5 * 8) = 4 exactly rows 2 and 4 flag.
  expect_identical(sum(out$flagged, na.rm = TRUE), 2L)
})

test_that("as.data.frame.cier_index returns the tidy per-respondent frame", {
  df <- as.data.frame(cier_longstring(present_fixture()))
  expect_s3_class(df, "data.frame")
  expect_identical(names(df), c("value", "flagged"))
  expect_identical(nrow(df), 4L)
  expect_equal(df$value, c(1, 8, 1, 4))
})

# ---- Independent oracle parity (exact) --------------------------------------

test_that("cier_longstring$value equals the hand-rolled oracle on a fixture", {
  x <- present_fixture()
  expect_equal(cier_longstring(x)$value, ref_longstring(x)$longest)
})

test_that("cier_longstring$value equals the oracle on a random complete matrix", {
  withr::with_seed(2026, {
    x <- matrix(sample.int(5L, 40L * 12L, replace = TRUE), nrow = 40L)
  })
  storage.mode(x) <- "double"
  expect_equal(cier_longstring(x)$value, ref_longstring(x)$longest)
})

# ---- Cross-package parity: bytewise vs careless (tolerance 0) ----------------

test_that("cier_longstring matches careless::longstring on a random matrix", {
  skip_if_not_installed("careless")
  withr::with_seed(2026, {
    x <- matrix(sample.int(5L, 20L * 30L, replace = TRUE), nrow = 20L)
  })
  storage.mode(x) <- "double"
  expect_equal(cier_longstring(x)$value,
               as.numeric(careless::longstring(x)), tolerance = 0)
})

test_that("cier_longstring matches careless on careless_dataset (bytewise)", {
  skip_if_not_installed("careless")
  x <- unname(as.matrix(careless::careless_dataset))
  storage.mode(x) <- "double"
  expect_equal(cier_longstring(x)$value,
               as.numeric(careless::longstring(careless::careless_dataset)),
               tolerance = 0)
})

# ---- Property / invariant + mutant-killers ----------------------------------

test_that("straightliner -> p and strict alternation -> 1 (max, not run count)", {
  p <- 10L
  expect_identical(cier_longstring(matrix(rep(3, p), nrow = 1L))$value, 10)
  alt <- matrix(rep(c(1, 4), length.out = p), nrow = 1L)
  expect_identical(cier_longstring(alt)$value, 1)
})

test_that("value is the max run length, not a count of long runs (sum>1 mutant)", {
  # runs are (3, 2, 1): max = 3, count-of-runs = 3, sum(lengths>1) = 2.
  expect_identical(cier_longstring(matrix(c(1, 1, 1, 2, 2, 3), nrow = 1L))$value, 3)
})

test_that("value is the max run length, not the average run length", {
  # runs (3, 1, 1, 1): max = 3, average = 1.5.
  expect_identical(cier_longstring(matrix(c(5, 5, 5, 1, 2, 3), nrow = 1L))$value, 3)
})

test_that("NA breaks runs and is never treated as equal to NA (no merge mutant)", {
  # rle gives runs (1),(NA),(NA),(1) all length 1 -> max 1. A mutant in which
  # NA == NA merges the two NAs would return 2.
  expect_identical(cier_longstring(matrix(c(1, NA, NA, 1), nrow = 1L))$value, 1)
})

test_that("NAs are not dropped before the run scan (no na.rm mutant)", {
  # c(5, NA, 5, NA, 5): all runs length 1 -> max 1. A mutant that drops NAs
  # first would see c(5, 5, 5) -> max 3.
  expect_identical(cier_longstring(matrix(c(5, NA, 5, NA, 5), nrow = 1L))$value, 1)
})

test_that("longstring uses raw column order (no scale-block / sort mutant)", {
  # Raw adjacency: runs (1),(2),(1,1),(2),(1) -> max 2. A mutant that sorted or
  # grouped columns would see c(1,1,1,1,2,2) -> max 4.
  expect_identical(cier_longstring(matrix(c(1, 2, 1, 1, 2, 1), nrow = 1L))$value, 2)
})

test_that("values lie in [1, p] for every present-data row", {
  withr::with_seed(7, {
    x <- matrix(sample.int(5L, 50L * 30L, replace = TRUE), nrow = 50L)
  })
  storage.mode(x) <- "double"
  v <- cier_longstring(x)$value
  expect_true(all(v >= 1 & v <= ncol(x)))
})

# ---- Edge cases -------------------------------------------------------------

test_that("an all-NA row abstains: value NA, flagged NA, excluded from n_flagged", {
  x <- matrix(c(rep(3, 6L), rep(NA_real_, 6L)), nrow = 2L, byrow = TRUE)
  out <- cier_longstring(x)
  expect_true(is.na(out$value[[2L]]))
  expect_true(is.na(out$flagged[[2L]]))
  expect_true(out$flagged[[1L]])
  expect_identical(sum(out$flagged, na.rm = TRUE), 1L)
})

test_that("an abstaining row in the middle keeps value/flagged row-aligned", {
  # Guards against a row-indexing / apply mutant that mis-aligns value and
  # flagged when an all-NA row is not the last row. p = 6, cutoff 3.
  x <- matrix(
    c(3, 3, 3, 3, 3, 3,        # row 1: value 6 -> flag
      NA, NA, NA, NA, NA, NA,  # row 2: abstains
      1, 2, 3, 4, 5, 1),       # row 3: value 1 -> no flag
    nrow = 3L, byrow = TRUE
  )
  out <- cier_longstring(x)
  expect_equal(out$value, c(6, NA, 1))
  expect_identical(out$flagged, c(TRUE, NA, FALSE))
  expect_identical(sum(out$flagged, na.rm = TRUE), 1L)
})

test_that("a single-column matrix yields longstring 1 for every row", {
  expect_identical(cier_longstring(matrix(c(3, 4, 5), ncol = 1L))$value, c(1, 1, 1))
})

test_that("single-row and constant matrices behave", {
  expect_identical(cier_longstring(matrix(c(1, 1, 2, 3), nrow = 1L))$value, 2)
  expect_identical(cier_longstring(matrix(rep(2, 20L), nrow = 2L))$value, c(10, 10))
})

# ---- Cutoff: default, direction, boundary -----------------------------------

test_that("default cutoff is ceiling(0.5 * p) and flagging includes ties (>=)", {
  # p = 10 -> cutoff 5. Row A longest run 5 (flag), row B longest run 4 (no).
  x <- matrix(
    c(3, 3, 3, 3, 3, 1, 2, 3, 4, 5,
      3, 3, 3, 3, 1, 2, 3, 4, 5, 1),
    nrow = 2L, byrow = TRUE
  )
  out <- cier_longstring(x)
  expect_identical(out$cutoff, 5)
  expect_identical(out$flagged, c(TRUE, FALSE))
})

test_that("direction is upper: high values flag, low values do not (flag-lower mutant)", {
  x <- matrix(
    c(rep(3, 10L),                       # straightliner -> value 10 -> flag
      rep(c(1, 4), length.out = 10L)),   # alternating   -> value 1  -> no flag
    nrow = 2L, byrow = TRUE
  )
  expect_identical(cier_longstring(x)$flagged, c(TRUE, FALSE))
})

# ---- Cutoff: fraction-or-count override -------------------------------------

test_that("a fractional cutoff in (0, 1] is a fraction of item count", {
  x <- matrix(rep(3, 10L), nrow = 1L)            # p = 10, value 10
  expect_identical(cier_longstring(x, cutoff = 0.4)$cutoff, 4) # ceiling(4)
})

test_that("an integer cutoff > 1 is an absolute count that drives the flags (>=)", {
  # GAP 1: drive the flag through the count path, not just the stored attribute.
  # p = 10. Row A longest run 5 (flag at count 5), row B longest run 4 (no).
  x <- matrix(
    c(3, 3, 3, 3, 3, 1, 2, 3, 4, 5,
      3, 3, 3, 3, 1, 2, 3, 4, 5, 1),
    nrow = 2L, byrow = TRUE
  )
  out <- cier_longstring(x, cutoff = 5)
  expect_identical(out$cutoff, 5)
  expect_identical(out$flagged, c(TRUE, FALSE))  # value 5 >= 5 flags (tie)
})

test_that("a count cutoff equal to p is accepted; p + 1 is rejected (boundary)", {
  # GAP 3: pin the > p reject boundary exactly, not only the far-out 999.
  x <- matrix(rep(3, 10L), nrow = 1L)            # p = 10, value 10
  out <- cier_longstring(x, cutoff = 10)
  expect_identical(out$cutoff, 10)
  expect_true(out$flagged)                       # 10 >= 10
  expect_error(cier_longstring(x, cutoff = 11), class = "cier_error_input")
})

test_that("a non-integer count > 1 is used as the literal cutoff (>= semantics)", {
  # GAP 4: contract pinned -- a count > 1 passes through literally (archive
  # parity, no rounding); >= against integer run lengths means 5.5 flags >= 6.
  x <- matrix(
    c(3, 3, 3, 3, 3, 3, 1, 2, 3, 4,   # longest run 6
      3, 3, 3, 3, 3, 1, 2, 3, 4, 5),  # longest run 5
    nrow = 2L, byrow = TRUE
  )
  out <- cier_longstring(x, cutoff = 5.5)
  expect_identical(out$cutoff, 5.5)
  expect_identical(out$flagged, c(TRUE, FALSE)) # 6 >= 5.5; 5 < 5.5
})

test_that("cutoff = 1 is the fraction path (ceiling(1 * p) = p), not the count 1", {
  x <- matrix(
    c(rep(3, 10L),                       # value 10 == p -> flagged at cutoff p
      rep(c(1, 4), length.out = 10L)),   # value 1       -> not flagged
    nrow = 2L, byrow = TRUE
  )
  out <- cier_longstring(x, cutoff = 1)
  expect_identical(out$cutoff, 10)
  expect_identical(out$flagged, c(TRUE, FALSE))
})

test_that("invalid cutoff overrides are typed input errors", {
  x <- matrix(rep(3, 40L), nrow = 5L, ncol = 8L)  # p = 8
  expect_error(cier_longstring(x, cutoff = -1), class = "cier_error_input")
  expect_error(cier_longstring(x, cutoff = 0), class = "cier_error_input")
  expect_error(cier_longstring(x, cutoff = NA_real_), class = "cier_error_input")
  expect_error(cier_longstring(x, cutoff = c(1, 2)), class = "cier_error_input")
  expect_error(cier_longstring(x, cutoff = 999), class = "cier_error_input") # > p
  expect_error(cier_longstring(x, cutoff = "x"), class = "cier_error_input")
})

# ---- print snapshot (locked, design-first) ----------------------------------

test_that("print renders the locked cli summary (no abstention)", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    x <- matrix(
      c(3, 3, 3, 3, 3, 3,   # value 6 -> flag
        1, 2, 3, 4, 5, 1,   # value 1
        2, 2, 2, 1, 1, 1,   # value 3 -> flag
        1, 2, 1, 2, 1, 2,   # value 1
        4, 4, 4, 4, 1, 2,   # value 4 -> flag
        5, 1, 5, 1, 5, 1),  # value 1
      nrow = 6L, byrow = TRUE
    )
    expect_snapshot(print(cier_longstring(x)))
  })
})

test_that("print reports abstaining respondents on their own line", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    x <- matrix(
      c(3, 3, 3, 3, 3, 3,        # value 6 -> flag
        1, 2, 3, 4, 5, 1,        # value 1
        2, 2, 2, 1, 1, 1,        # value 3 -> flag
        1, 2, 1, 2, 1, 2,        # value 1
        4, 4, 4, 4, 1, 2,        # value 4 -> flag
        NA, NA, NA, NA, NA, NA), # abstains
      nrow = 6L, byrow = TRUE
    )
    expect_snapshot(print(cier_longstring(x)))
  })
})
