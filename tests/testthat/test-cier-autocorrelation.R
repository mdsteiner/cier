# Tests for cier_autocorrelation().
#
# Trust model: the independent oracle (ref_autocorrelation) re-derives the
# per-respondent max-absolute-lag autocorrelation with base-R var() + cor() and
# never calls the production kernel; the cross-package check pins parity with
# responsePatterns::rp.acors() (the authors' own implementation) on complete
# data at the 1e-10 tolerance in tests/reference/TOLERANCES.md. NA-bearing parity
# is oracle-only: responsePatterns 0.1.1 has two verified NA bugs (a var(row1)
# checked-twice crash and a na.rm=TRUE over-indexing slice) that cier deliberately
# avoids, so the partner is exercised only on NA-free fixtures.

source(test_path("..", "reference", "ref-autocorrelation-gottfried2022.R"))

# ---- Fixtures ---------------------------------------------------------------

ac_fixture <- function(n = 20L, p = 30L, seed = 2026L) {
  withr::with_seed(seed, {
    m <- matrix(sample.int(5L, n * p, replace = TRUE), nrow = n, ncol = p)
  })
  storage.mode(m) <- "double"
  m
}

# Hand fixture with analytically known maxima:
#   row 1: constant 3              -> every slice zero-variance -> value 1
#   row 2: seesaw 1,5,1,5,...      -> lag-1 cor -1, lag-2 cor 1 -> value 1
#   row 3: diagonal 1..5 repeated  -> lag-5 cor 1 -> value 1
#   row 4: random                  -> also saturates to 1 at the default
#                                     max_lag = p - 3 (short high-lag slices),
#                                     so it is NOT asserted here; the
#                                     de-saturation case is its own test below
pattern_fixture <- function(p = 30L) {
  withr::with_seed(4L, rnd <- sample.int(5L, p, replace = TRUE))
  m <- rbind(
    rep(3L, p),
    rep(c(1L, 5L), length.out = p),
    rep(seq_len(5L), length.out = p),
    rnd
  )
  storage.mode(m) <- "double"
  m
}

# ---- Schema -----------------------------------------------------------------

test_that("cier_autocorrelation returns a list-based cier_index schema", {
  out <- cier_autocorrelation(ac_fixture(n = 12L))
  expect_s3_class(out, "cier_index")
  expect_type(out, "list")
  expect_identical(names(out),
                   c("value", "flagged", "method", "cutoff", "direction"))
  expect_type(out$value, "double")
  expect_type(out$flagged, "logical")
  expect_identical(length(out$value), 12L)
  expect_identical(length(out$flagged), 12L)
  expect_identical(out$method, "cier_autocorrelation")
  expect_identical(out$direction, "upper")
  expect_type(out$cutoff, "double")
})

test_that("as.data.frame.cier_index returns the tidy per-respondent frame", {
  df <- as.data.frame(cier_autocorrelation(ac_fixture(n = 8L)))
  expect_s3_class(df, "data.frame")
  expect_identical(names(df), c("value", "flagged"))
  expect_identical(nrow(df), 8L)
})

# ---- Independent oracle parity (1e-10) --------------------------------------

test_that("cier_autocorrelation$value equals the oracle on complete data", {
  x <- ac_fixture(n = 25L, p = 30L)
  expect_equal(cier_autocorrelation(x)$value,
               ref_autocorrelation_value(x), tolerance = 1e-10)
})

test_that("oracle parity holds with explicit min_lag / max_lag", {
  x <- ac_fixture(n = 20L, p = 30L)
  expect_equal(cier_autocorrelation(x, min_lag = 2L, max_lag = 10L)$value,
               ref_autocorrelation_value(x, min_lag = 2L, max_lag = 10L),
               tolerance = 1e-10)
})

test_that("max_lag defaults to ncol - 3 (off-by-one default is caught)", {
  # Non-saturating continuous data so the boundary lag actually moves the max:
  # dropping lag ncol-3 (a ncol-4 default) changes 14 of 40 rows here.
  withr::with_seed(5L, x <- matrix(stats::runif(40L * 15L, 1, 5), nrow = 40L))
  v_default <- cier_autocorrelation(x)$value
  expect_equal(v_default, cier_autocorrelation(x, max_lag = ncol(x) - 3L)$value,
               tolerance = 1e-12)                       # default == ncol - 3
  dropped <- cier_autocorrelation(x, max_lag = ncol(x) - 4L)$value
  expect_false(isTRUE(all.equal(v_default, dropped)))
  expect_equal(v_default, ref_autocorrelation_value(x), tolerance = 1e-10)
})

test_that("min_lag defaults to 1 (lag 1 is included)", {
  # min_lag 1 vs 2 differ on 5 of 40 rows at max_lag = 6, so the default is
  # pinned independent of the production cutoff.
  withr::with_seed(5L, x <- matrix(stats::runif(40L * 15L, 1, 5), nrow = 40L))
  v_default <- cier_autocorrelation(x, max_lag = 6L)$value
  expect_equal(v_default,
               cier_autocorrelation(x, min_lag = 1L, max_lag = 6L)$value,
               tolerance = 1e-12)                        # default == explicit 1
  from_lag2 <- cier_autocorrelation(x, min_lag = 2L, max_lag = 6L)$value
  expect_false(isTRUE(all.equal(v_default, from_lag2)))
})

test_that("na_rm = FALSE handles NAs pairwise, matching the oracle", {
  mat <- rbind(
    c(1, 2, 3, NA, 5, 4, 3, 2),   # one NA in the middle (pairwise drop)
    c(5, 4, 3, 2, 1, 2, 3, 4),
    c(2, 4, 1, 5, 3, 1, 4, 2)
  )
  storage.mode(mat) <- "double"
  expect_equal(cier_autocorrelation(mat)$value,
               ref_autocorrelation_value(mat), tolerance = 1e-10)
})

test_that("na_rm = TRUE strips per row before lagging, matching the oracle", {
  # This row scores differently under the two NA paths (strip 0.756 vs pairwise
  # 0.945), so it directly distinguishes na_rm = TRUE from the default and kills
  # an "always pairwise" mutant on the value, not just on abstention.
  mat <- matrix(c(3, 5, 2, 3, 2, NA, NA, 3, NA, 4, 4, 5), nrow = 1L)
  strip <- ref_autocorrelation_value(mat, na_rm = TRUE)[[1L]]
  pairwise <- ref_autocorrelation_value(mat, na_rm = FALSE)[[1L]]
  expect_false(isTRUE(all.equal(strip, pairwise)))     # the fixture is diagnostic
  expect_equal(cier_autocorrelation(mat, na_rm = TRUE)$value[[1L]],
               strip, tolerance = 1e-10)
})

test_that("scattered-NA rows match the oracle on both NA paths", {
  # Four rows with different NA patterns at max_lag = 4, so several lags reduce
  # to 1 or 2 complete pairs (the n_pairs * (n_pairs - 1) denominator boundary)
  # and the stripped lengths differ across rows -- pinning the few-complete-pair
  # handling AND the na_rm = TRUE per-row alignment, neither of which the
  # single-row fixtures above reach.
  m <- rbind(
    c(2, 4, NA, 5, 1, 3, NA, 2, 4, 1),
    c(NA, 3, 1, 5, 2, NA, 4, 1, NA, 3),
    c(1, NA, NA, 2, 5, 4, NA, NA, 3, 2),
    c(5, 4, 3, 2, 1, 2, 3, 4, 5, 4)
  )
  storage.mode(m) <- "double"
  expect_equal(cier_autocorrelation(m, max_lag = 4L)$value,
               ref_autocorrelation_value(m, max_lag = 4L), tolerance = 1e-10)
  expect_equal(cier_autocorrelation(m, max_lag = 4L, na_rm = TRUE)$value,
               ref_autocorrelation_value(m, max_lag = 4L, na_rm = TRUE),
               tolerance = 1e-10)
})

# ---- Cross-package parity: responsePatterns::rp.acors (1e-10, complete) -----

test_that("cier_autocorrelation matches rp.acors on complete data (na_rm FALSE)", {
  skip_if_not_installed("responsePatterns")
  x <- ac_fixture(n = 20L, p = 30L, seed = 7L)
  ours <- cier_autocorrelation(x)$value
  rp <- suppressMessages(
    responsePatterns::rp.acors(as.data.frame(x), na.rm = FALSE)
  )
  theirs <- methods::slot(rp, "indices")$max.abs.ac
  expect_equal(ours, theirs, tolerance = 1e-10)
})

# ---- Convention pins: zero-variance slice scores 1 (NOT abstention) ---------

test_that("a straightliner scores exactly 1 and is flagged, never NA", {
  out <- cier_autocorrelation(matrix(rep(3, 30L), nrow = 1L))
  expect_identical(out$value, 1)
  expect_false(is.na(out$value))
  expect_true(out$flagged[[1L]])
})

test_that("seesaw and diagonal patterns score 1 (max |ac|)", {
  pf <- pattern_fixture(p = 30L)
  v <- cier_autocorrelation(pf)$value
  expect_identical(v[[1L]], 1)               # constant: zero-var branch, exact 1
  expect_equal(v[[2L]], 1, tolerance = 1e-10)  # seesaw 1,5,1,5,... (lag-2 cor 1)
  expect_equal(v[[3L]], 1, tolerance = 1e-10)  # diagonal 1..5 repeated (lag-5 cor 1)
})

test_that("values are bounded in [0, 1]; a moderate max_lag de-saturates", {
  # At the default max_lag = p - 3 the short high-lag slices saturate |ac| to 1
  # on discrete Likert data (the paper recommends lag <= 12); a moderate max_lag
  # spreads the scores strictly inside the unit interval.
  withr::with_seed(5L, x <- matrix(stats::runif(40L * 15L, 1, 5), nrow = 40L))
  v <- cier_autocorrelation(x, max_lag = 6L)$value
  expect_true(all(v >= 0 & v <= 1 + 1e-12))
  expect_true(any(v < 0.95))
})

# ---- NA / abstention edges --------------------------------------------------

test_that("na_rm = TRUE: a row stripping to < 3 elements abstains", {
  mat <- rbind(
    c(1, 2, 3, 4, 5, 4, 3, 2, 1, 2),       # scored
    c(5, 4, 3, 2, 1, 2, 3, 4, 5, 4),       # scored
    c(2, 4, rep(NA_real_, 8L))             # strips to 2 -> abstains
  )
  storage.mode(mat) <- "double"
  out <- cier_autocorrelation(mat, na_rm = TRUE)
  expect_true(is.na(out$value[[3L]]))
  expect_true(is.na(out$flagged[[3L]]))
  expect_false(is.na(out$value[[1L]]))
})

test_that("an all-NA row abstains and keeps rows aligned", {
  mat <- rbind(
    c(1, 2, 3, 4, 5, 4, 3, 2, 1, 2),
    rep(NA_real_, 10L),                    # abstains, in the middle
    c(5, 4, 3, 2, 1, 2, 3, 4, 5, 4)
  )
  storage.mode(mat) <- "double"
  out <- cier_autocorrelation(mat)
  expect_true(is.na(out$value[[2L]]))
  expect_true(is.na(out$flagged[[2L]]))
  expect_false(is.na(out$value[[1L]]))
  expect_false(is.na(out$value[[3L]]))
})

test_that("a wholly abstaining matrix warns and flags nobody", {
  expect_warning(
    out <- cier_autocorrelation(matrix(NA_real_, nrow = 3L, ncol = 30L)),
    class = "cier_warning_insufficient_items"
  )
  expect_true(all(is.na(out$value)))
  expect_true(is.na(out$cutoff))
  expect_true(all(is.na(out$flagged)))
})

# ---- Lag-window edges -------------------------------------------------------

test_that("p = 4 works with the default max_lag (= 1) and matches the oracle", {
  x <- ac_fixture(n = 6L, p = 4L, seed = 3L)
  out <- cier_autocorrelation(x)
  expect_s3_class(out, "cier_index")
  expect_identical(length(out$value), 6L)
  # The narrowest lag window (a single lag-1 slice over 4 columns) is value-
  # verified, not just shape-checked: a wrong single-lag slice or a default
  # off-by-one would change these numbers.
  expect_equal(out$value, ref_autocorrelation_value(x), tolerance = 1e-10)
})

test_that("fewer than four columns is a typed input error", {
  x <- ac_fixture(n = 5L, p = 3L)
  expect_error(cier_autocorrelation(x), class = "cier_error_input")
})

test_that("an unreconcilable lag range is a typed input error", {
  x <- ac_fixture(n = 6L, p = 10L)
  # min_lag > max_lag
  expect_error(cier_autocorrelation(x, min_lag = 5L, max_lag = 2L),
               class = "cier_error_input")
  # max_lag beyond ncol - 3 (= 7)
  expect_error(cier_autocorrelation(x, max_lag = 8L),
               class = "cier_error_input")
  # ncol - 3 (= 7) is the largest ACCEPTED max_lag (guards a too-tight bound)
  expect_s3_class(cier_autocorrelation(x, max_lag = 7L), "cier_index")
})

test_that("bad min_lag / max_lag are typed input errors", {
  x <- ac_fixture(n = 6L, p = 10L)
  expect_error(cier_autocorrelation(x, min_lag = 0L), class = "cier_error_input")
  expect_error(cier_autocorrelation(x, min_lag = NA_integer_),
               class = "cier_error_input")
  expect_error(cier_autocorrelation(x, min_lag = c(1L, 2L)),
               class = "cier_error_input")
  expect_error(cier_autocorrelation(x, min_lag = 1.5),
               class = "cier_error_input")
  expect_error(cier_autocorrelation(x, max_lag = 0L), class = "cier_error_input")
})

# ---- Cutoff: default, fpr override, NO double-flip, direction ---------------

test_that("default cutoff is the upper-tail 95th percentile (NO double-flip)", {
  withr::with_seed(5L, x <- matrix(stats::runif(60L * 12L, 1, 5), nrow = 60L))
  out <- cier_autocorrelation(x)
  # Upper direction takes the 1 - fpr quantile directly (NOT fpr): the registry
  # stores the literal fpr tail mass and the resolver must not re-flip.
  expect_equal(out$cutoff,
               as.numeric(stats::quantile(out$value, 0.95,
                                          names = FALSE, type = 7L)),
               tolerance = 1e-12)
})

test_that("the fpr argument moves the percentile target", {
  withr::with_seed(5L, x <- matrix(stats::runif(60L * 12L, 1, 5), nrow = 60L))
  out <- cier_autocorrelation(x, fpr = 0.10)
  expect_equal(out$cutoff,
               as.numeric(stats::quantile(out$value, 0.90,
                                          names = FALSE, type = 7L)),
               tolerance = 1e-12)
})

test_that("direction is upper: high-ac rows flag, low ones do not", {
  withr::with_seed(11L, rnd <- matrix(sample.int(5L, 20L * 16L, replace = TRUE),
                                      nrow = 20L))
  x <- rbind(rep(3, 16L), rnd)            # constant row -> value 1 (the maximum)
  storage.mode(x) <- "double"
  out <- cier_autocorrelation(x)
  expect_true(out$flagged[[1L]])                          # constant, value 1
  expect_false(out$flagged[[which.min(out$value)]])       # least autocorrelated
  expect_identical(out$flagged,
                   !is.na(out$value) & out$value >= out$cutoff)
})

test_that("a literal cutoff passes through and flags via the upper direction", {
  withr::with_seed(7L, rnd <- matrix(sample.int(5L, 10L * 16L, replace = TRUE),
                                     nrow = 10L))
  x <- rbind(rep(3, 16L), rnd)
  storage.mode(x) <- "double"
  out <- cier_autocorrelation(x, cutoff = 0.95)
  expect_identical(out$cutoff, 0.95)
  expect_true(out$flagged[[1L]])                          # value 1 >= 0.95
  expect_identical(out$flagged, !is.na(out$value) & out$value >= 0.95)
})

test_that("supplying both fpr and cutoff is a typed input error", {
  expect_error(cier_autocorrelation(ac_fixture(n = 6L), fpr = 0.1, cutoff = 0.5),
               class = "cier_error_input")
})

test_that("an invalid literal cutoff is a typed input error", {
  x <- ac_fixture(n = 6L)
  expect_error(cier_autocorrelation(x, cutoff = 1.5), class = "cier_error_input")
  expect_error(cier_autocorrelation(x, cutoff = -0.1), class = "cier_error_input")
  expect_error(cier_autocorrelation(x, cutoff = NA_real_),
               class = "cier_error_input")
  expect_error(cier_autocorrelation(x, cutoff = c(0.1, 0.2)),
               class = "cier_error_input")
  expect_error(cier_autocorrelation(x, cutoff = "x"), class = "cier_error_input")
})

test_that("invalid fpr values are typed input errors", {
  x <- ac_fixture(n = 6L)
  expect_error(cier_autocorrelation(x, fpr = 0), class = "cier_error_input")
  expect_error(cier_autocorrelation(x, fpr = 1), class = "cier_error_input")
  expect_error(cier_autocorrelation(x, fpr = -0.1), class = "cier_error_input")
  expect_error(cier_autocorrelation(x, fpr = NA_real_), class = "cier_error_input")
  expect_error(cier_autocorrelation(x, fpr = c(0.05, 0.1)),
               class = "cier_error_input")
  expect_error(cier_autocorrelation(x, fpr = "x"), class = "cier_error_input")
})

# ---- Input validation -------------------------------------------------------

test_that("a non-matrix / non-numeric payload is a typed input error", {
  expect_error(cier_autocorrelation(1:10), class = "cier_error_input")
  expect_error(cier_autocorrelation(matrix(letters[1:8], nrow = 2L)),
               class = "cier_error_input")
  expect_error(cier_autocorrelation(NULL), class = "cier_error_input")
})

test_that("a non-flag na_rm is a typed input error", {
  x <- ac_fixture(n = 6L)
  expect_error(cier_autocorrelation(x, na_rm = NA), class = "cier_error_input")
  expect_error(cier_autocorrelation(x, na_rm = c(TRUE, FALSE)),
               class = "cier_error_input")
  expect_error(cier_autocorrelation(x, na_rm = 1), class = "cier_error_input")
})

# ---- print snapshot (locked; reuses the shared cier_index print) ------------

test_that("print renders the locked cli summary (upper direction)", {
  # A moderate max_lag so the cutoff is a representative quantile rather than the
  # saturated ceiling of 1 the default max_lag produces on discrete data.
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    expect_snapshot(
      print(cier_autocorrelation(ac_fixture(n = 30L, p = 12L), max_lag = 5L))
    )
  })
})

test_that("print reports abstaining respondents on their own line", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    x <- ac_fixture(n = 29L, p = 12L)
    x <- rbind(x, rep(NA_real_, 12L))     # one abstaining respondent
    expect_snapshot(print(cier_autocorrelation(x, max_lag = 5L)))
  })
})
