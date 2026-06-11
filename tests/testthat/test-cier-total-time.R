# Tests for cier_total_time() (timing family; Ward & Meade 2023, Huang et al.
# 2012).
#
# Trust model: total_time's per-respondent value is the validated seconds vector
# itself (an identity), so the independent oracle (ref-total-time.R) re-derives
# the CONTRACT -- the identity value and the three cutoff resolutions -- by base
# R and never calls the production wrapper. There is no CRAN parity partner (no
# package implements completion time as a C/IER index; verified 2026-06-10), so
# the oracle plus the base-R quantile / median primitives are the parity checks,
# oracle-only trust like PR / RPR. Conventions under test (signed off this
# slice): bare-numeric-vector input; strictly-positive seconds (zero / negative
# are input errors, NA abstains); empirical lower percentile (fpr 0.05) default;
# a third median-relative override frac_median in (0, 1]; the three cutoff knobs
# mutually exclusive.

source(test_path("..", "reference", "ref-total-time.R"))

# ---- Fixtures ---------------------------------------------------------------

# Five respondents with a hand-checkable median of 30 (so frac_median = 0.5 ->
# cutoff 15, flagging only the 10 s respondent) and a clear fastest / slowest.
small_seconds <- function() c(10, 20, 30, 40, 50)

# A realistic spread for the percentile-default and print snapshots.
spread_seconds <- function(n = 30L, seed = 11L) {
  withr::with_seed(seed, stats::runif(n, min = 30, max = 300))
}

# ---- Schema -----------------------------------------------------------------

test_that("cier_total_time returns a list-based cier_index with the schema", {
  out <- cier_total_time(small_seconds())
  expect_s3_class(out, "cier_index")
  expect_type(out, "list")
  expect_identical(names(out),
                   c("value", "flagged", "method", "cutoff", "direction"))
  expect_type(out$value, "double")
  expect_type(out$flagged, "logical")
  expect_identical(length(out$value), 5L)
  expect_identical(length(out$flagged), 5L)
  expect_identical(out$method, "cier_total_time")
  expect_identical(out$direction, "lower")
  expect_type(out$cutoff, "double")
})

test_that("as.data.frame.cier_index returns the tidy per-respondent frame", {
  df <- as.data.frame(cier_total_time(small_seconds()))
  expect_s3_class(df, "data.frame")
  expect_identical(names(df), c("value", "flagged"))
  expect_identical(nrow(df), 5L)
  expect_identical(df$value, small_seconds())
})

# ---- Identity oracle (tolerance 0) ------------------------------------------

test_that("the value is the seconds vector verbatim (no transform)", {
  # The mutant-killer for any silent sum / scale / rank / log: the value must be
  # byte-identical to the input, NA preserved.
  x <- spread_seconds()
  expect_identical(cier_total_time(x)$value, ref_total_time(x))
})

test_that("an integer seconds vector is coerced to double, values unchanged", {
  x <- c(40L, 80L, 120L, 160L)
  out <- cier_total_time(x)
  expect_type(out$value, "double")
  expect_identical(out$value, as.numeric(x))
})

test_that("a single respondent scores (no minimum-n floor)", {
  # total_time needs no pair / second item the way IRV and Laz.R do, so one
  # finite value is a valid score; the percentile cutoff over one value is that
  # value. Pins that a spurious n < 2 abstention guard was not added.
  out <- cier_total_time(120)
  expect_identical(out$value, 120)
  expect_false(is.na(out$cutoff))
  expect_equal(out$cutoff, 120, tolerance = 1e-12)
})

# ---- Cutoff: percentile default, fpr override, NO-FLIP -----------------------

test_that("default cutoff is the lower-tail 5th percentile (NO double-flip)", {
  x <- spread_seconds(n = 60L)
  out <- cier_total_time(x)
  # Lower direction takes the fpr quantile DIRECTLY (not 1 - fpr).
  expect_equal(out$cutoff,
               ref_total_time_percentile_cutoff(x, 0.05),
               tolerance = 1e-12)
  expect_equal(out$cutoff,
               as.numeric(stats::quantile(x, 0.05, names = FALSE, type = 7L)),
               tolerance = 1e-12)
})

test_that("the fpr argument moves the percentile target", {
  x <- spread_seconds(n = 60L)
  out <- cier_total_time(x, fpr = 0.10)
  expect_equal(out$cutoff,
               ref_total_time_percentile_cutoff(x, 0.10),
               tolerance = 1e-12)
})

# ---- Cutoff: median-relative override (frac_median) -------------------------

test_that("frac_median resolves frac * median of the finite values", {
  # median(c(10,20,30,40,50)) = 30; frac_median 0.5 -> cutoff 15.
  out <- cier_total_time(small_seconds(), frac_median = 0.5)
  expect_equal(out$cutoff, 15, tolerance = 1e-12)
  expect_equal(out$cutoff,
               ref_total_time_median_cutoff(small_seconds(), 0.5),
               tolerance = 1e-12)
})

test_that("frac_median flags via the lower comparator (value <= cutoff)", {
  # cutoff 15 -> only the 10 s respondent is at or below it.
  out <- cier_total_time(small_seconds(), frac_median = 0.5)
  expect_identical(out$flagged, c(TRUE, FALSE, FALSE, FALSE, FALSE))
  expect_identical(out$flagged,
                   ref_total_time_flags(small_seconds(), out$cutoff))
})

test_that("frac_median is frac * median, not (1-frac)*median or median/frac", {
  # An asymmetric frac (0.25 != 0.75) discriminates the formula on a fixture
  # written for that purpose: 0.25 * median(=30) = 7.5, whereas (1-frac)*median
  # = 22.5 and median/frac = 120. (The frac_median = 1 test alone only rules out
  # (1-frac) by coincidence.)
  out <- cier_total_time(small_seconds(), frac_median = 0.25)
  expect_equal(out$cutoff, 7.5, tolerance = 1e-12)
})

test_that("frac_median uses the median, not the mean (mutant-killer)", {
  # Right-skewed: mean(c(10,20,30,40,200)) = 60, median = 30. frac 0.5 ->
  # cutoff 15 (median) not 30 (mean); only the 10 s row flags either way, but
  # the CUTOFF value distinguishes the two.
  x <- c(10, 20, 30, 40, 200)
  out <- cier_total_time(x, frac_median = 0.5)
  expect_equal(out$cutoff, 15, tolerance = 1e-12)   # 0.5 * median(=30)
})

test_that("frac_median ignores NA rows when taking the median", {
  # median of the finite values c(10,20,30,40,50) = 30 despite the NA.
  x <- c(10, 20, NA, 30, 40, 50)
  out <- cier_total_time(x, frac_median = 0.5)
  expect_equal(out$cutoff, 15, tolerance = 1e-12)
})

# ---- Cutoff: literal override -----------------------------------------------

test_that("a literal cutoff passes through and flags via the lower direction", {
  out <- cier_total_time(small_seconds(), cutoff = 25)
  expect_identical(out$cutoff, 25)
  expect_identical(out$flagged, c(TRUE, TRUE, FALSE, FALSE, FALSE))  # 10, 20
})

test_that("a literal cutoff of 0 is accepted and flags nobody", {
  # Seconds are strictly positive, so value <= 0 is never true; 0 is a valid
  # (degenerate) threshold rather than an input error.
  out <- cier_total_time(small_seconds(), cutoff = 0)
  expect_identical(out$cutoff, 0)
  expect_identical(out$flagged, rep(FALSE, 5L))
})

# ---- Three-way mutual exclusion ---------------------------------------------

test_that("supplying more than one cutoff knob is a typed input error", {
  x <- small_seconds()
  expect_error(cier_total_time(x, fpr = 0.1, cutoff = 20),
               class = "cier_error_input")
  expect_error(cier_total_time(x, fpr = 0.1, frac_median = 0.5),
               class = "cier_error_input")
  expect_error(cier_total_time(x, frac_median = 0.5, cutoff = 20),
               class = "cier_error_input")
  expect_error(cier_total_time(x, fpr = 0.1, frac_median = 0.5, cutoff = 20),
               class = "cier_error_input")
})

# ---- Direction --------------------------------------------------------------

test_that("direction is lower: fast rows flag, slow ones do not", {
  x <- spread_seconds(n = 40L)
  out <- cier_total_time(x)
  expect_true(out$flagged[[which.min(out$value)]])    # fastest respondent
  expect_false(out$flagged[[which.max(out$value)]])   # slowest respondent
  expect_identical(out$flagged,
                   !is.na(out$value) & out$value <= out$cutoff)
})

# ---- Abstention / NA edges --------------------------------------------------

test_that("an NA row abstains and keeps rows aligned", {
  # The abstainer sits in the middle, so value and flagged must stay aligned to
  # their respondents.
  x <- c(40, NA, 120)
  out <- cier_total_time(x)
  expect_true(is.na(out$value[[2L]]))
  expect_true(is.na(out$flagged[[2L]]))
  expect_false(is.na(out$value[[1L]]))
  expect_false(is.na(out$value[[3L]]))
  expect_identical(out$value[[1L]], 40)
})

test_that("an all-NA vector abstains: percentile warns and flags nobody", {
  expect_warning(
    out <- cier_total_time(rep(NA_real_, 5L)),
    class = "cier_warning_insufficient_items"
  )
  expect_true(all(is.na(out$value)))
  expect_true(is.na(out$cutoff))
  expect_true(all(is.na(out$flagged)))
})

test_that("frac_median on an all-NA vector warns and flags nobody", {
  # The median-relative resolver mirrors the percentile abstention: no finite
  # value -> NA cutoff + the same typed warning, flagging no one.
  expect_warning(
    out <- cier_total_time(rep(NA_real_, 5L), frac_median = 0.5),
    class = "cier_warning_insufficient_items"
  )
  expect_true(is.na(out$cutoff))
  expect_true(all(is.na(out$flagged)))
})

# ---- Input validation -------------------------------------------------------

test_that("a non-vector / non-numeric payload is a typed input error", {
  # Bare numeric vector only: a matrix, data.frame, list, or character input is
  # rejected (2-D is ambiguous about which axis is respondents).
  expect_error(cier_total_time(matrix(c(40, 80, 120, 160), ncol = 1L)),
               class = "cier_error_input")
  expect_error(cier_total_time(matrix(c(40, 80, 120, 160), nrow = 2L)),
               class = "cier_error_input")
  expect_error(cier_total_time(data.frame(t = c(40, 80))),
               class = "cier_error_input")
  expect_error(cier_total_time(list(40, 80)), class = "cier_error_input")
  expect_error(cier_total_time(c("40", "80")), class = "cier_error_input")
  expect_error(cier_total_time(NULL), class = "cier_error_input")
})

test_that("an empty vector is a typed input error", {
  expect_error(cier_total_time(numeric(0)), class = "cier_error_input")
})

test_that("NaN / infinite seconds are typed input errors", {
  expect_error(cier_total_time(c(40, NaN, 120)), class = "cier_error_input")
  expect_error(cier_total_time(c(40, Inf, 120)), class = "cier_error_input")
})

test_that("zero and negative seconds are typed input errors", {
  # A completion time cannot be <= 0; recode missing durations to NA instead.
  expect_error(cier_total_time(c(40, 0, 120)), class = "cier_error_input")
  expect_error(cier_total_time(c(40, -5, 120)), class = "cier_error_input")
})

test_that("a small-but-positive time is accepted and preserved", {
  # The boundary is strictly > 0, OPEN below: a sub-second time is valid, not an
  # "implausibly fast" input error. Kills a `value < 1` (minimum-plausible-time)
  # guard that would silently reject genuine fast respondents.
  out <- expect_no_error(cier_total_time(c(0.001, 40, 120)))
  expect_identical(out$value[[1L]], 0.001)
})

test_that("bad fpr values are typed input errors", {
  x <- small_seconds()
  expect_error(cier_total_time(x, fpr = 0), class = "cier_error_input")
  expect_error(cier_total_time(x, fpr = 1), class = "cier_error_input")
  expect_error(cier_total_time(x, fpr = -0.1), class = "cier_error_input")
  expect_error(cier_total_time(x, fpr = NA_real_), class = "cier_error_input")
  expect_error(cier_total_time(x, fpr = c(0.05, 0.1)),
               class = "cier_error_input")
  expect_error(cier_total_time(x, fpr = "x"), class = "cier_error_input")
})

test_that("bad frac_median values are typed input errors", {
  x <- small_seconds()
  expect_error(cier_total_time(x, frac_median = 0), class = "cier_error_input")
  expect_error(cier_total_time(x, frac_median = 1.5),
               class = "cier_error_input")   # > 1 is out of (0, 1]
  expect_error(cier_total_time(x, frac_median = -0.1),
               class = "cier_error_input")
  expect_error(cier_total_time(x, frac_median = NA_real_),
               class = "cier_error_input")
  expect_error(cier_total_time(x, frac_median = c(0.4, 0.5)),
               class = "cier_error_input")
  expect_error(cier_total_time(x, frac_median = "x"),
               class = "cier_error_input")
})

test_that("frac_median = 1 is accepted (the (0, 1] upper bound)", {
  # The closed upper end: 1.0 flags the at-or-below-median half. Pins that the
  # bound is (0, 1], not (0, 1).
  out <- cier_total_time(small_seconds(), frac_median = 1)
  expect_equal(out$cutoff, 30, tolerance = 1e-12)   # 1 * median(=30)
  expect_identical(out$flagged, c(TRUE, TRUE, TRUE, FALSE, FALSE))
})

test_that("invalid literal cutoff values are typed input errors", {
  x <- small_seconds()
  expect_error(cier_total_time(x, cutoff = -1), class = "cier_error_input")
  expect_error(cier_total_time(x, cutoff = NA_real_),
               class = "cier_error_input")
  expect_error(cier_total_time(x, cutoff = c(20, 30)),
               class = "cier_error_input")
  expect_error(cier_total_time(x, cutoff = "x"), class = "cier_error_input")
})

# ---- print snapshots (locked, design-first; shared cier_index print) --------

test_that("print renders the locked cli summary (lower direction)", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    expect_snapshot(print(cier_total_time(spread_seconds())))
  })
})

test_that("print reports abstaining respondents on their own line", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    x <- c(spread_seconds(n = 29L), NA_real_)   # one abstaining respondent
    expect_snapshot(print(cier_total_time(x)))
  })
})
