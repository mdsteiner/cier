# Tests for cier_page_time() (timing family; Bowling, Huang, Brower & Bragg
# 2023).
#
# Trust model: page_time is a counting rule -- per respondent, the number of
# pages whose mean per-item time (page total / items on the page) falls strictly
# below `min_seconds`. The independent oracle (ref-page-time.R) re-derives that
# count by a per-respondent row loop and never calls the production kernel; there
# is no CRAN parity partner (no package implements page time as a C/IER index;
# verified 2026-06-10), so the oracle plus base-R primitives are the parity
# checks -- oracle-only trust, like total_time / PR / RPR. Conventions under test
# (signed off this slice): lean input (an n x pages matrix of page TOTAL times +
# an explicit items-per-page vector, replacing the archive's cier_data + per-cell
# times + inferred page boundaries); page-level NA (a missing page contributes no
# evidence; the denominator is the declared items_per_page, not a count of
# answered cells); strictly-positive page times (zero / negative are input
# errors, NA abstains); the cited fixed cutoff of 1 (any rapid page flags) with a
# proportion override `frac` (a fraction of the page count) and a literal `cutoff`
# override, mutually exclusive.

source(test_path("..", "reference", "ref-page-time.R"))

# ---- Fixtures ---------------------------------------------------------------

# Four respondents, three two-item pages, page TOTAL times. Hand-checkable per
# the comments; mean per-item = total / 2.
small_pages <- function() {
  matrix(c(10, 10, 10,    # per-item c(5, 5, 5)   -> 0 rapid
           2, 10, 10,     # per-item c(1, 5, 5)   -> 1 rapid
           4, 10, 10,     # per-item c(2, 5, 5)   -> page 1 == 2.0, NOT < 2 -> 0
           2, 2, 2),      # per-item c(1, 1, 1)   -> 3 rapid
         nrow = 4L, byrow = TRUE)
}
small_ipp <- function() c(2L, 2L, 2L)

# A realistic spread for the oracle-parity and print snapshots: random positive
# page totals, varying items per page, n != n_pages (so a transposed / wrong-axis
# division is caught).
spread_pages <- function(n = 30L, n_pages = 8L, seed = 11L) {
  ipp <- rep_len(c(2L, 3L, 4L, 1L), n_pages)
  m <- withr::with_seed(seed, matrix(stats::runif(n * n_pages, 1, 30),
                                     nrow = n, ncol = n_pages))
  list(page_seconds = m, items_per_page = ipp)
}

# ---- Schema -----------------------------------------------------------------

test_that("cier_page_time returns a list-based cier_index with the schema", {
  out <- cier_page_time(small_pages(), small_ipp())
  expect_s3_class(out, "cier_index")
  expect_type(out, "list")
  expect_identical(names(out),
                   c("value", "flagged", "method", "cutoff", "direction"))
  expect_type(out$value, "double")
  expect_type(out$flagged, "logical")
  expect_identical(length(out$value), 4L)
  expect_identical(length(out$flagged), 4L)
  expect_identical(out$method, "cier_page_time")
  expect_identical(out$direction, "upper")
  expect_type(out$cutoff, "double")
})

test_that("as.data.frame.cier_index returns the tidy per-respondent frame", {
  df <- as.data.frame(cier_page_time(small_pages(), small_ipp()))
  expect_s3_class(df, "data.frame")
  expect_identical(names(df), c("value", "flagged"))
  expect_identical(nrow(df), 4L)
  expect_identical(df$value, c(0, 1, 0, 3))
})

# ---- Oracle parity (tolerance 0, integer counts) ----------------------------

test_that("the value equals the independent counting-rule oracle exactly", {
  fx <- spread_pages()
  out <- cier_page_time(fx$page_seconds, fx$items_per_page)
  expect_identical(out$value,
                   ref_page_time(fx$page_seconds, fx$items_per_page, 2))
})

test_that("page-time counts pages whose mean per-item time < min_seconds", {
  out <- cier_page_time(small_pages(), small_ipp(), cutoff = 1)
  expect_equal(out$value, c(0, 1, 0, 3))
  expect_identical(out$value,
                   ref_page_time(small_pages(), small_ipp(), 2))
})

# ---- Normalization pin: divide the page total by items_per_page -------------

test_that("rapidity is the per-ITEM mean, not the page total (divide pin)", {
  # One respondent, one 4-item page, total 6 s -> per-item 1.5 < 2 -> rapid.
  # A mutant that compares the page TOTAL (6) to min_seconds (2) finds it NOT
  # rapid, scoring 0 instead of 1.
  out <- cier_page_time(matrix(6, nrow = 1L, ncol = 1L), 4L)
  expect_equal(out$value, 1)
})

test_that("items_per_page divides COLUMN-wise (wrong-axis recycle pin)", {
  # n = 2 respondents, 3 pages, items_per_page varies across pages c(2, 4, 1) and
  # n != n_pages, so a naive `page_seconds / items_per_page` (column-major
  # recycle) or a transpose mis-divides and the oracle parity diverges.
  ps <- matrix(c(2, 12, 0.5,    # per-item c(1, 3, 0.5)  -> pages 1, 3 rapid -> 2
                 10, 4, 5),     # per-item c(5, 1, 5)    -> page 2 rapid     -> 1
               nrow = 2L, byrow = TRUE)
  ipp <- c(2L, 4L, 1L)
  out <- cier_page_time(ps, ipp)
  expect_equal(out$value, c(2, 1))
  expect_identical(out$value, ref_page_time(ps, ipp, 2))
})

# ---- Bowling 2 s/item boundary: strictly below ------------------------------

test_that("Bowling boundary is strict: per-item 1.9 flags, 2.0 does not", {
  # Two-item pages; mean per-item 1.9 (total 3.8) is rapid, 2.0 (total 4.0) not.
  ps <- matrix(c(3.8, 10,    # per-item c(1.9, 5) -> 1 rapid
                 4.0, 10),   # per-item c(2.0, 5) -> 0 rapid
               nrow = 2L, byrow = TRUE)
  out <- cier_page_time(ps, c(2L, 2L), cutoff = 1)
  expect_equal(out$value, c(1, 0))
})

# ---- min_seconds default and override ---------------------------------------

test_that("min_seconds defaults to 2 (Bowling)", {
  # Per-item 1.9 is rapid only because the default threshold is 2; a mutant with
  # a different default (e.g. 1) would score it 0.
  ps <- matrix(c(3.8, 10), nrow = 1L)
  expect_equal(cier_page_time(ps, c(2L, 2L))$value, 1)
})

test_that("min_seconds override moves the threshold", {
  # Per-item c(2.5, 5): rapid at min_seconds 3, not at the default 2.
  ps <- matrix(c(5, 10), nrow = 1L)
  expect_equal(cier_page_time(ps, c(2L, 2L), min_seconds = 3)$value, 1)
  expect_equal(cier_page_time(ps, c(2L, 2L))$value, 0)
})

test_that("min_seconds = 0 flags no page (page times are strictly positive)", {
  out <- cier_page_time(small_pages(), small_ipp(), min_seconds = 0)
  expect_equal(out$value, c(0, 0, 0, 0))
})

test_that("a fractional min_seconds above 1 is accepted (check_number lower 0)", {
  # min_seconds 0 (above) and 1.5 (here) together pin the validator as
  # check_number(lower = 0): a count / open-unit / fraction check would reject
  # one of them (0, or a value > 1).
  out <- expect_no_error(
    cier_page_time(matrix(c(1, 10), nrow = 1L), c(1L, 1L), min_seconds = 1.5)
  )
  expect_equal(out$value, 1)   # per-item c(1, 10); 1 < 1.5 rapid, 10 not
})

# ---- NA pages and abstention ------------------------------------------------

test_that("an NA page contributes no evidence; all-NA pages abstain", {
  # 3 respondents, 2 two-item pages.
  # Row 1: c(2, NA)   -> per-item c(1, .) -> page 1 rapid, page 2 no evidence -> 1
  # Row 2: c(NA, NA)  -> no timed page -> NA (abstains)
  # Row 3: c(10, NA)  -> per-item c(5, .) -> 0 rapid
  ps <- matrix(c(2, NA,
                 NA, NA,
                 10, NA),
               nrow = 3L, byrow = TRUE)
  # page_time uses a FIXED cutoff, so an abstaining respondent must NOT route the
  # cutoff through the percentile abstention (no cier_warning_insufficient_items).
  out <- expect_no_warning(cier_page_time(ps, c(2L, 2L)))
  expect_equal(out$value[c(1L, 3L)], c(1, 0))
  expect_true(is.na(out$value[[2L]]))
  expect_true(is.na(out$flagged[[2L]]))        # abstainer's flag is NA, aligned
  expect_identical(out$value, ref_page_time(ps, c(2L, 2L), 2))
  # The fixed cutoff stays 1 despite the abstaining row, and the scored rows
  # still flag against it (kills a percentile-misroute mutant that would abstain
  # the cutoff to NA and flag nobody).
  expect_identical(out$cutoff, 1)
  expect_true(out$flagged[[1L]])
  expect_false(out$flagged[[3L]])
})

test_that("one NA page does NOT make the respondent abstain (all-NA pin)", {
  # A single timed (rapid) page is enough to score; only all-NA abstains. Kills a
  # mutant that abstains whenever ANY page is NA.
  ps <- matrix(c(2, NA, NA), nrow = 1L)   # 1 timed page (rapid), 2 untimed
  out <- cier_page_time(ps, c(2L, 2L, 2L))
  expect_equal(out$value, 1)
  expect_false(is.na(out$value))
})

# ---- Fixed-cutoff contract (any rapid page) ---------------------------------

test_that("default cutoff is fixed 1: any rapid page flags (>= pin)", {
  out <- cier_page_time(small_pages(), small_ipp())
  expect_identical(out$cutoff, 1)
  # value c(0, 1, 0, 3) -> flag where count >= 1. A count of exactly 1 MUST flag
  # (kills a `>` mutant).
  expect_identical(out$flagged, c(FALSE, TRUE, FALSE, TRUE))
  expect_identical(out$flagged,
                   ref_page_time_flags(out$value, out$cutoff))
})

test_that("a literal cutoff override flags via value >= cutoff", {
  out <- cier_page_time(small_pages(), small_ipp(), cutoff = 2)
  expect_identical(out$cutoff, 2)
  # value c(0, 1, 0, 3) -> only the 3-rapid-page respondent clears 2.
  expect_identical(out$flagged, c(FALSE, FALSE, FALSE, TRUE))
})

# ---- frac override (proportion of pages) ------------------------------------

test_that("frac resolves ceiling(frac * n_pages) as the rapid-page cutoff", {
  # 4 single-item pages; frac 0.5 -> ceiling(2) = 2.
  ps <- matrix(c(1, 1, 10, 10,    # 2 rapid pages -> meets cutoff 2
                 1, 10, 10, 10),  # 1 rapid page  -> below cutoff 2
               nrow = 2L, byrow = TRUE)
  out <- cier_page_time(ps, rep(1L, 4L), frac = 0.5)
  expect_equal(out$cutoff, 2)
  expect_identical(out$flagged, c(TRUE, FALSE))
  expect_equal(out$cutoff, ref_page_time_fixed_cutoff(0.5, 4L))
})

test_that("frac rounds up: 0.25 -> 1, 0.6 -> 3 over four pages", {
  ps <- matrix(c(1, 10, 10, 10), nrow = 1L)   # 1 rapid page
  expect_equal(cier_page_time(ps, rep(1L, 4L), frac = 0.25)$cutoff, 1)
  expect_equal(cier_page_time(ps, rep(1L, 4L), frac = 0.6)$cutoff, 3)
})

test_that("frac eases the single-fast-page sensitivity on a long survey", {
  # The motivating case: 40 single-item pages, exactly ONE rapid. The cited
  # default (cutoff 1) flags the respondent; frac = 0.1 -> ceiling(4) = 4
  # requires a tenth of pages to be rapid, so one rapid page no longer flags.
  ps <- matrix(c(1, rep(10, 39L)), nrow = 1L)
  expect_true(cier_page_time(ps, rep(1L, 40L))$flagged[[1L]])
  out <- cier_page_time(ps, rep(1L, 40L), frac = 0.1)
  expect_equal(out$cutoff, 4)
  expect_false(out$flagged[[1L]])
})

test_that("frac carries the round-to-9dp guard (0.28 * 25 -> 7, not 8)", {
  ps <- matrix(stats::runif(25L, 1, 30), nrow = 1L)   # value irrelevant here
  out <- cier_page_time(ps, rep(1L, 25L), frac = 0.28)
  expect_equal(out$cutoff, 7)
})

test_that("frac = 1 is accepted: every page must be rapid", {
  ps <- matrix(c(1, 1, 1,    # 3 rapid -> meets cutoff 3
                 1, 1, 10),  # 2 rapid -> below cutoff 3
               nrow = 2L, byrow = TRUE)
  out <- cier_page_time(ps, rep(1L, 3L), frac = 1)
  expect_equal(out$cutoff, 3)
  expect_identical(out$flagged, c(TRUE, FALSE))
})

test_that("frac is a fraction of the TOTAL page count, not timed pages", {
  # One all-NA column (page 3 untimed for everyone): the survey has 5 pages, so
  # frac 0.5 -> ceiling(0.5 * 5) = 3. A mutant that takes the fraction over the
  # 4 timed pages resolves ceiling(0.5 * 4) = 2 and survives every fully-timed
  # frac fixture above -- this NA-column fixture is what discriminates them.
  ps <- matrix(c(1, 1, NA, 10, 10,
                 1, 10, NA, 10, 10),
               nrow = 2L, byrow = TRUE)
  out <- cier_page_time(ps, rep(1L, 5L), frac = 0.5)
  expect_equal(out$cutoff, 3)
  expect_equal(out$cutoff, ref_page_time_fixed_cutoff(0.5, 5L))
})

test_that("frac and cutoff are mutually exclusive", {
  expect_error(cier_page_time(small_pages(), small_ipp(),
                              frac = 0.5, cutoff = 2),
               class = "cier_error_input")
})

# ---- Direction --------------------------------------------------------------

test_that("direction is upper: more rapid pages flag", {
  out <- cier_page_time(small_pages(), small_ipp())
  expect_identical(out$direction, "upper")
  expect_true(out$flagged[[which.max(out$value)]])   # most rapid pages
  expect_false(out$flagged[[which.min(out$value)]])  # fewest rapid pages
})

test_that("registry flag_direction matches the wrapper direction", {
  declared <- cier_method_row("cier_page_time")$flag_direction
  out <- cier_page_time(small_pages(), small_ipp())
  expect_identical(out$direction, declared)
})

# ---- page_seconds input validation ------------------------------------------

test_that("a non-matrix / non-numeric page_seconds is a typed input error", {
  expect_error(cier_page_time(c(2, 4, 6), c(2L)), class = "cier_error_input")
  expect_error(cier_page_time(list(2, 4), c(2L)), class = "cier_error_input")
  expect_error(cier_page_time(matrix("a", 1L, 1L), 1L),
               class = "cier_error_input")
  expect_error(cier_page_time(NULL, 1L), class = "cier_error_input")
})

test_that("a data.frame page_seconds is accepted (coerced)", {
  df <- as.data.frame(small_pages())
  expect_no_error(cier_page_time(df, small_ipp()))
})

test_that("an empty page_seconds is a typed input error", {
  expect_error(cier_page_time(matrix(numeric(0), nrow = 0L, ncol = 0L), 1L),
               class = "cier_error_input")
})

test_that("NaN / infinite page_seconds are typed input errors", {
  expect_error(cier_page_time(matrix(c(2, NaN), nrow = 1L), c(2L, 2L)),
               class = "cier_error_input")
  expect_error(cier_page_time(matrix(c(2, Inf), nrow = 1L), c(2L, 2L)),
               class = "cier_error_input")
})

test_that("zero and negative page_seconds are typed input errors", {
  # A page time cannot be <= 0; recode untimed pages to NA instead.
  expect_error(cier_page_time(matrix(c(2, 0), nrow = 1L), c(2L, 2L)),
               class = "cier_error_input")
  expect_error(cier_page_time(matrix(c(2, -5), nrow = 1L), c(2L, 2L)),
               class = "cier_error_input")
})

test_that("a small-but-positive page time is accepted", {
  # The boundary is strictly > 0, OPEN below: a sub-second page total is valid.
  out <- expect_no_error(cier_page_time(matrix(c(0.001, 10), nrow = 1L),
                                        c(1L, 1L)))
  expect_equal(out$value, 1)   # 0.001 / 1 < 2 -> rapid
})

# ---- items_per_page input validation ----------------------------------------

test_that("items_per_page of the wrong length is a typed input error", {
  expect_error(cier_page_time(small_pages(), c(2L, 2L)),   # need length 3
               class = "cier_error_input")
})

test_that("items_per_page with NA / fractional / non-positive is an error", {
  expect_error(cier_page_time(small_pages(), c(2L, NA, 2L)),
               class = "cier_error_input")
  expect_error(cier_page_time(small_pages(), c(2, 2.5, 2)),
               class = "cier_error_input")
  expect_error(cier_page_time(small_pages(), c(2L, 0L, 2L)),
               class = "cier_error_input")
  expect_error(cier_page_time(small_pages(), c(2L, -1L, 2L)),
               class = "cier_error_input")
})

test_that("a non-numeric / 2-D items_per_page is a typed input error", {
  expect_error(cier_page_time(small_pages(), c("a", "b", "c")),
               class = "cier_error_input")
  expect_error(cier_page_time(small_pages(), matrix(2L, 3L, 1L)),
               class = "cier_error_input")
})

# ---- min_seconds / cutoff / frac argument validation ------------------------

test_that("bad min_seconds values are typed input errors", {
  expect_error(cier_page_time(small_pages(), small_ipp(), min_seconds = -1),
               class = "cier_error_input")
  expect_error(cier_page_time(small_pages(), small_ipp(), min_seconds = NA_real_),
               class = "cier_error_input")
  expect_error(cier_page_time(small_pages(), small_ipp(), min_seconds = Inf),
               class = "cier_error_input")
  expect_error(cier_page_time(small_pages(), small_ipp(), min_seconds = c(1, 2)),
               class = "cier_error_input")
  expect_error(cier_page_time(small_pages(), small_ipp(), min_seconds = "a"),
               class = "cier_error_input")
})

test_that("invalid literal cutoff values are typed input errors", {
  np <- ncol(small_pages())   # 3 pages -> cutoff in [1, 3]
  expect_error(cier_page_time(small_pages(), small_ipp(), cutoff = 0),
               class = "cier_error_input")
  expect_error(cier_page_time(small_pages(), small_ipp(), cutoff = np + 1L),
               class = "cier_error_input")
  expect_error(cier_page_time(small_pages(), small_ipp(), cutoff = NA_real_),
               class = "cier_error_input")
  expect_error(cier_page_time(small_pages(), small_ipp(), cutoff = c(1, 2)),
               class = "cier_error_input")
  expect_error(cier_page_time(small_pages(), small_ipp(), cutoff = "a"),
               class = "cier_error_input")
})

test_that("invalid frac values are typed input errors", {
  expect_error(cier_page_time(small_pages(), small_ipp(), frac = 0),
               class = "cier_error_input")
  expect_error(cier_page_time(small_pages(), small_ipp(), frac = 1.5),
               class = "cier_error_input")
  expect_error(cier_page_time(small_pages(), small_ipp(), frac = -0.1),
               class = "cier_error_input")
  expect_error(cier_page_time(small_pages(), small_ipp(), frac = NA_real_),
               class = "cier_error_input")
  expect_error(cier_page_time(small_pages(), small_ipp(), frac = c(0.4, 0.5)),
               class = "cier_error_input")
  expect_error(cier_page_time(small_pages(), small_ipp(), frac = "a"),
               class = "cier_error_input")
})

# ---- print snapshots (shared cier_index print; upper direction) -------------

test_that("print renders the locked cli summary (upper direction)", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    expect_snapshot(print(cier_page_time(small_pages(), small_ipp())))
  })
})

test_that("print reports abstaining respondents on their own line", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    # Row 2 is all-NA -> abstains; the others score.
    ps <- matrix(c(2, 10, 10,
                   NA, NA, NA,
                   2, 2, 10,
                   10, 10, 10),
                 nrow = 4L, byrow = TRUE)
    expect_snapshot(print(cier_page_time(ps, c(2L, 2L, 2L))))
  })
})
