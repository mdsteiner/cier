# Tests for cier_page_time() (timing family; Bowling, Huang, Brower & Bragg 2023).
#
# Trust model: page_time counts, per respondent, the pages whose mean per-item time
# (page total / items on the page) falls strictly below `min_seconds`. The oracle
# (ref-page-time.R) re-derives the count by a per-respondent row loop and never
# calls the kernel; no CRAN parity partner (no package implements page time as a
# C/IER index; verified 2026-06-10), so oracle + base-R primitives are the parity
# checks -- oracle-only trust, like total_time / PR / RPR. Conventions under test:
# lean input (n x pages matrix of page TOTAL times + an explicit items-per-page
# vector, replacing the archive's cier_data + per-cell times + inferred page
# boundaries); page-level NA (a missing page is no evidence; denominator is the
# declared items_per_page, not answered cells); strictly-positive page times (zero
# / negative are input errors, NA abstains); the cited fixed cutoff of 1 (any rapid
# page flags) with a proportion override `frac` (a fraction of the page count) and
# a literal `cutoff` override, mutually exclusive.

source(test_path("..", "reference", "ref-page-time.R"))

# ---- Fixtures ---------------------------------------------------------------

# Four respondents, three two-item pages, page TOTAL times. Mean per-item =
# total / 2; hand-checkable per the trailing comments.
small_pages <- function() {
  matrix(c(10, 10, 10,    # per-item c(5, 5, 5)   -> 0 rapid
           2, 10, 10,     # per-item c(1, 5, 5)   -> 1 rapid
           4, 10, 10,     # per-item c(2, 5, 5)   -> page 1 == 2.0, NOT < 2 -> 0
           2, 2, 2),      # per-item c(1, 1, 1)   -> 3 rapid
         nrow = 4L, byrow = TRUE)
}
small_ipp <- function() c(2L, 2L, 2L)

# A realistic spread for oracle-parity and print snapshots: random positive page
# totals, varying items per page, n != n_pages (catches a transposed / wrong-axis
# division).
spread_pages <- function(n = 30L, n_pages = 8L, seed = 11L) {
  ipp <- rep_len(c(2L, 3L, 4L, 1L), n_pages)
  m <- withr::with_seed(seed, matrix(stats::runif(n * n_pages, 1, 30),
                                     nrow = n, ncol = n_pages))
  list(page_seconds = m, items_per_page = ipp)
}

# ---- Schema -----------------------------------------------------------------

test_that("cier_page_time returns the pinned cier_index schema", {
  expect_cier_index_schema(cier_page_time(small_pages(), small_ipp()),
                           "cier_page_time", "upper", 4L)
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
  # One respondent, one 4-item page, total 6 s -> per-item 1.5 < 2 -> rapid. A
  # mutant comparing the page TOTAL (6) to min_seconds (2) finds it NOT rapid,
  # scoring 0. 100% flagged here (n = 1) trips the unit-mismatch warning, incidental
  # to this per-item-mean pin -> muffled.
  out <- suppressWarnings(cier_page_time(matrix(6, nrow = 1L, ncol = 1L), 4L))
  expect_equal(out$value, 1)
})

test_that("items_per_page divides COLUMN-wise (wrong-axis recycle pin)", {
  # n = 2 respondents, 3 pages, items_per_page c(2, 4, 1) varies, n != n_pages, so a
  # naive `page_seconds / items_per_page` (column-major recycle) or a transpose
  # mis-divides and the oracle parity diverges.
  ps <- matrix(c(2, 12, 0.5,    # per-item c(1, 3, 0.5)  -> pages 1, 3 rapid -> 2
                 10, 4, 5),     # per-item c(5, 1, 5)    -> page 2 rapid     -> 1
               nrow = 2L, byrow = TRUE)
  ipp <- c(2L, 4L, 1L)
  out <- suppressWarnings(cier_page_time(ps, ipp))   # both flagged -> warn muffled
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
  expect_equal(suppressWarnings(cier_page_time(ps, c(2L, 2L)))$value, 1)
})

test_that("min_seconds override moves the threshold", {
  # Per-item c(2.5, 5): rapid at min_seconds 3, not at the default 2.
  ps <- matrix(c(5, 10), nrow = 1L)
  # min_seconds = 3 flags the lone respondent (100%) -> all-flagged warn muffled.
  out <- suppressWarnings(cier_page_time(ps, c(2L, 2L), min_seconds = 3))
  expect_equal(out$value, 1)
  expect_equal(cier_page_time(ps, c(2L, 2L))$value, 0)
})

test_that("min_seconds = 0 flags no page (page times are strictly positive)", {
  out <- cier_page_time(small_pages(), small_ipp(), min_seconds = 0)
  expect_equal(out$value, c(0, 0, 0, 0))
})

test_that("a fractional min_seconds above 1 is accepted (check_number lower 0)", {
  # min_seconds 0 (above) and 1.5 (here) together pin the validator as
  # check_number(lower = 0): a count / open-unit / fraction check rejects one of
  # them (0, or a value > 1).
  out <- expect_no_error(suppressWarnings(   # 100% flagged -> unit warn muffled
    cier_page_time(matrix(c(1, 10), nrow = 1L), c(1L, 1L), min_seconds = 1.5)
  ))
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
  # The fixed cutoff stays 1 despite the abstaining row, and the scored rows still
  # flag against it (kills a percentile-misroute mutant that abstains the cutoff to
  # NA and flags nobody).
  expect_identical(out$cutoff, 1)
  expect_true(out$flagged[[1L]])
  expect_false(out$flagged[[3L]])
})

test_that("one NA page does NOT make the respondent abstain (all-NA pin)", {
  # A single timed (rapid) page is enough to score; only all-NA abstains. Kills a
  # mutant abstaining whenever ANY page is NA.
  ps <- matrix(c(2, NA, NA), nrow = 1L)   # 1 timed page (rapid), 2 untimed
  out <- suppressWarnings(cier_page_time(ps, c(2L, 2L, 2L)))  # 100% -> warn muffled
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
  # frac 0.25 -> cutoff 1 -> the single respondent flags (100%) -> warn muffled.
  expect_equal(suppressWarnings(cier_page_time(ps, rep(1L, 4L), frac = 0.25))$cutoff, 1)
  expect_equal(cier_page_time(ps, rep(1L, 4L), frac = 0.6)$cutoff, 3)
})

test_that("frac eases the single-fast-page sensitivity on a long survey", {
  # The motivating case: 40 single-item pages, exactly ONE rapid. The cited default
  # (cutoff 1) flags the respondent; frac = 0.1 -> ceiling(4) = 4 needs a tenth of
  # pages rapid, so one rapid page no longer flags.
  ps <- matrix(c(1, rep(10, 39L)), nrow = 1L)
  # The default cutoff flags the lone respondent (100%) -> all-flagged warn muffled.
  flagged_default <- suppressWarnings(cier_page_time(ps, rep(1L, 40L)))$flagged[[1L]]
  expect_true(flagged_default)
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
  # frac 0.5 -> ceiling(0.5 * 5) = 3. A mutant taking the fraction over the 4 timed
  # pages resolves ceiling(0.5 * 4) = 2 and survives every fully-timed frac fixture
  # above -- this NA-column fixture discriminates them.
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

test_that("method spec flag_direction matches the wrapper direction", {
  declared <- cier_method_spec("cier_page_time")$flag_direction
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

test_that("a negative page time is a typed input error", {
  # A page time cannot be negative; recode untimed pages to NA. Zero is now valid
  # (see below), so only the negative is rejected.
  expect_error(cier_page_time(matrix(c(2, -5), nrow = 1L), c(2L, 2L)),
               class = "cier_error_input")
})

test_that("the negative-page-time error says 'zero or greater' / 'negative'", {
  # The reworded hint drops the old "strictly positive / <= 0" wording: zero is now
  # valid, only a negative is rejected.
  msg <- tryCatch(cier_page_time(matrix(c(2, -5), nrow = 1L), c(2L, 2L)),
                  error = function(e) gsub("\\s+", " ", rlang::cnd_message(e)))
  expect_match(msg, "zero or greater", fixed = TRUE)
  expect_match(msg, "cannot be negative", fixed = TRUE)
})

test_that("zero and small-positive page times are accepted; 0 counts rapid", {
  # A 0 page total is maximal speeding evidence (per-item 0 < min_seconds), not an
  # input error; a sub-second total is likewise valid. Three rows so not every
  # scored respondent flags (keeps the all-flagged warning out).
  ps <- matrix(c(0, 10,        # per-item c(0, 5)     -> page 1 rapid -> 1
                 0.001, 10,    # per-item c(.0005, 5) -> page 1 rapid -> 1
                 10, 10),      # per-item c(5, 5)     -> 0 rapid      -> 0
               nrow = 3L, byrow = TRUE)
  out <- cier_page_time(ps, c(2L, 2L))
  expect_equal(out$value, c(1, 1, 0))
  expect_identical(out$flagged, c(TRUE, TRUE, FALSE))
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

test_that("a fractional literal cutoff is a typed input error", {
  # A page-count cutoff is a whole number of rapid pages; 2.5 is meaningless.
  expect_error(cier_page_time(small_pages(), small_ipp(), cutoff = 2.5),
               class = "cier_error_input")
  expect_s3_class(cier_page_time(small_pages(), small_ipp(), cutoff = 2),
                  "cier_index")
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

# ---- all-flagged unit-mismatch warning --------------------------------------

test_that("every scored respondent flagged emits the unit-mismatch warning", {
  # 12 respondents, every first page rapid -> 100% flagged -> the heuristic that the
  # times may not be in seconds (minutes / another unit make every per-item mean
  # fall below min_seconds).
  ps <- cbind(rep(1, 12L), rep(10, 12L))   # per-item c(1, 10): page 1 always rapid
  expect_warning(out <- cier_page_time(ps, c(1L, 1L)),
                 class = "cier_warning_all_flagged")
  expect_true(all(out$flagged))
  w <- testthat::capture_warnings(cier_page_time(ps, c(1L, 1L)))
  expect_length(w, 1L)
  expect_match(w, "every scored respondent", ignore.case = TRUE)
  expect_match(w, "seconds")
  # At the default threshold the hint stays the crisp unit-mismatch wording and does
  # NOT broaden to the min_seconds-set-high case (kills an always-broadened mutant).
  expect_false(any(grepl("high relative to your page times",
                         gsub("\\s+", " ", w), fixed = TRUE)))
})

test_that("a non-default min_seconds broadens the all-flagged hint", {
  # Raising min_seconds above the data's per-item times is a legitimate sensitivity
  # choice that also flags everyone. The hint must then name min_seconds (with its
  # value) as an alternative cause, not assert a unit mismatch as the sole reason.
  # (Page 1 per-item 1 s < 5; page 2 10 s not.)
  ps <- cbind(rep(1, 12L), rep(10, 12L))
  expect_warning(out <- cier_page_time(ps, c(1L, 1L), min_seconds = 5),
                 class = "cier_warning_all_flagged")
  expect_true(all(out$flagged))
  # Strip ANSI so the bare-digit match below sees only the interpolated value, not a
  # colour escape (cli styling can carry digits under a coloured theme).
  w <- testthat::capture_warnings(cier_page_time(ps, c(1L, 1L), min_seconds = 5))
  w <- gsub("\\s+", " ", cli::ansi_strip(w))
  expect_length(w, 1L)
  expect_match(w, "every scored respondent", ignore.case = TRUE)
  expect_match(w, "min_seconds", fixed = TRUE)
  expect_match(w, "high relative to your page times", fixed = TRUE)
  expect_match(w, "5", fixed = TRUE)   # shows the chosen threshold value
  expect_match(w, "seconds")           # still names the unit-mismatch possibility
})

test_that("an explicit min_seconds = 2 keeps the crisp hint (default value)", {
  # The branch is keyed on the default VALUE (2), not on whether the user typed it:
  # passing min_seconds = 2 explicitly must still give the crisp unit hint (kills a
  # missing()-keyed mutant that broadens any supplied value).
  ps <- cbind(rep(1, 12L), rep(10, 12L))
  w <- gsub("\\s+", " ",
            testthat::capture_warnings(cier_page_time(ps, c(1L, 1L),
                                                      min_seconds = 2)))
  expect_match(w, "every scored respondent", ignore.case = TRUE)
  expect_match(w, "seconds")
  expect_false(any(grepl("high relative to your page times", w, fixed = TRUE)))
})

test_that("the broadened hint interpolates the ACTUAL min_seconds value", {
  # Two distinct non-default thresholds: each message shows its own value, not the
  # other's. Kills a mutant hardcoding a literal threshold in the template (a single
  # asserted value cannot tell interpolation from a constant). ANSI stripped so the
  # bare-digit negative assertions see only interpolated numbers; `min_seconds` is
  # the only number in the text (`n_scored` lives in the condition's data slot).
  ps <- cbind(rep(1, 12L), rep(10, 12L))
  w4 <- testthat::capture_warnings(cier_page_time(ps, c(1L, 1L), min_seconds = 4))
  w4 <- gsub("\\s+", " ", cli::ansi_strip(w4))
  expect_match(w4, "set to 4", fixed = TRUE)
  expect_false(any(grepl("7", w4, fixed = TRUE)))
  w7 <- testthat::capture_warnings(cier_page_time(ps, c(1L, 1L), min_seconds = 7))
  w7 <- gsub("\\s+", " ", cli::ansi_strip(w7))
  expect_match(w7, "set to 7", fixed = TRUE)
  expect_false(any(grepl("4", w7, fixed = TRUE)))
})

test_that("all-flagged warns even at tiny n (no minimum-n threshold)", {
  # The spec deliberately has no min-n guard; n = 2 both flagged still warns. Pins
  # "no threshold" from the smallest side -- kills any `n_scored >= k` mutant.
  expect_warning(cier_page_time(rbind(c(1, 10), c(1, 10)), c(1L, 1L)),
                 class = "cier_warning_all_flagged")
})

test_that("a mixed flag set does NOT warn (negative case)", {
  ps <- cbind(c(rep(1, 6L), rep(10, 6L)), rep(10, 12L))   # half rapid, half not
  expect_no_warning(out <- cier_page_time(ps, c(1L, 1L)))
  expect_false(all(out$flagged))
})

test_that("an all-abstaining (all-NA) page set does NOT warn (no scored)", {
  # No scored respondent -> the all-flagged condition is vacuously false; the rows
  # abstain rather than triggering the unit warning.
  ps <- matrix(NA_real_, nrow = 5L, ncol = 2L)
  expect_no_warning(out <- cier_page_time(ps, c(1L, 1L)))
  expect_true(all(is.na(out$flagged)))
})

test_that("all-flagged with one abstainer still warns on the scored set", {
  # The condition is over SCORED respondents: an abstaining row does not block the
  # warning when every scored respondent flagged.
  ps <- rbind(matrix(rep(c(1, 10), 5L), ncol = 2L, byrow = TRUE),  # 5 rapid
              c(NA, NA))                                            # 1 abstains
  expect_warning(cier_page_time(ps, c(1L, 1L)),
                 class = "cier_warning_all_flagged")
})
