# Tests for cier_autocorrelation().
#
# Trust model: the oracle (ref_autocorrelation) re-derives the per-respondent
# max-absolute-lag autocorrelation with base-R var() + cor() and never calls the
# kernel; the cross-package check pins parity with responsePatterns::rp.acors()
# (the authors' own implementation) on complete data at 1e-10. NA-bearing parity
# is oracle-only: responsePatterns 0.1.1 has two verified NA bugs (a var(row1)
# checked-twice crash and a na.rm=TRUE over-indexing slice) that cier avoids, so
# the partner is exercised only on NA-free fixtures.

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
#                                     max_lag = p - 3 (short high-lag slices), so
#                                     NOT asserted here; de-saturation is its own
#                                     test below
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

test_that("cier_autocorrelation returns the pinned cier_index schema", {
  out <- suppressWarnings(cier_autocorrelation(ac_fixture(n = 12L)))
  expect_cier_index_schema(out, "cier_autocorrelation", "upper", 12L)
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

test_that("max_lag defaults to min(ncol - 3, 10) (cap and off-by-one caught)", {
  # The default caps the lag window at 10 (Gottfried et al.'s low-lag range); the
  # full ncol - 3 saturates short high-lag slices to 1 on long batteries. Continuous
  # non-saturating data so the boundary lags actually move the max. p = 15 ->
  # ncol - 3 = 12 > 10, so the cap binds: the default is 10, not 12.
  withr::with_seed(5L, x <- matrix(stats::runif(40L * 15L, 1, 5), nrow = 40L))
  v_default <- cier_autocorrelation(x)$value
  expect_equal(v_default, cier_autocorrelation(x, max_lag = 10L)$value,
               tolerance = 1e-12)                       # default == 10 (the cap)
  # Off-by-one either side moves the result: dropping lag 10 (max_lag = 9) and
  # reaching ncol - 3 = 12 (cap not applied) both differ from the capped default.
  expect_false(isTRUE(all.equal(v_default,
                                cier_autocorrelation(x, max_lag = 9L)$value)))
  expect_false(isTRUE(all.equal(v_default,
                                cier_autocorrelation(x, max_lag = 12L)$value)))
  expect_equal(v_default, ref_autocorrelation_value(x), tolerance = 1e-10)

  # The cap is a min(), not a floor: when ncol - 3 <= 10 it stays the default, so
  # p - 3 is still reachable. p = 12 -> default == ncol - 3 == 9.
  withr::with_seed(6L, x2 <- matrix(stats::runif(30L * 12L, 1, 5), nrow = 30L))
  expect_equal(cier_autocorrelation(x2)$value,
               cier_autocorrelation(x2, max_lag = 9L)$value, tolerance = 1e-12)
  expect_equal(cier_autocorrelation(x2)$value,
               ref_autocorrelation_value(x2), tolerance = 1e-10)
  # The default never resolves ABOVE ncol - 3: an explicit max_lag = 10 is illegal
  # at p = 12 (10 > 9), pinning out a min(ncol - 2, 10) off-by-one that would hide
  # behind the slice-too-short guard.
  expect_error(cier_autocorrelation(x2, max_lag = 10L),
               class = "cier_error_input")
})

test_that("min_lag defaults to 1 (lag 1 is included)", {
  # min_lag 1 vs 2 differ on 5 of 40 rows at max_lag = 6, so the default is pinned
  # independent of the production cutoff.
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
  expect_equal(suppressWarnings(cier_autocorrelation(mat))$value,
               ref_autocorrelation_value(mat), tolerance = 1e-10)
})

test_that("na_rm = TRUE strips per row before lagging, matching the oracle", {
  # This row scores differently under the two NA paths (strip 0.756 vs pairwise
  # 0.945), distinguishing na_rm = TRUE from the default and killing an "always
  # pairwise" mutant on the value, not just on abstention.
  mat <- matrix(c(3, 5, 2, 3, 2, NA, NA, 3, NA, 4, 4, 5), nrow = 1L)
  strip <- ref_autocorrelation_value(mat, na_rm = TRUE)[[1L]]
  pairwise <- ref_autocorrelation_value(mat, na_rm = FALSE)[[1L]]
  expect_false(isTRUE(all.equal(strip, pairwise)))     # the fixture is diagnostic
  expect_equal(suppressWarnings(cier_autocorrelation(mat, na_rm = TRUE))$value[[1L]],
               strip, tolerance = 1e-10)
})

test_that("scattered-NA rows match the oracle on both NA paths", {
  # Four rows with different NA patterns at max_lag = 4, so several lags reduce to
  # 1 or 2 complete pairs (the n_pairs * (n_pairs - 1) denominator boundary) and the
  # stripped lengths differ across rows -- pinning the few-complete-pair handling AND
  # the na_rm = TRUE per-row alignment, neither of which the single-row fixtures
  # above reach.
  m <- rbind(
    c(2, 4, NA, 5, 1, 3, NA, 2, 4, 1),
    c(NA, 3, 1, 5, 2, NA, 4, 1, NA, 3),
    c(1, NA, NA, 2, 5, 4, NA, NA, 3, 2),
    c(5, 4, 3, 2, 1, 2, 3, 4, 5, 4)
  )
  storage.mode(m) <- "double"
  expect_equal(suppressWarnings(cier_autocorrelation(m, max_lag = 4L))$value,
               ref_autocorrelation_value(m, max_lag = 4L), tolerance = 1e-10)
  expect_equal(
    suppressWarnings(cier_autocorrelation(m, max_lag = 4L, na_rm = TRUE))$value,
    ref_autocorrelation_value(m, max_lag = 4L, na_rm = TRUE),
    tolerance = 1e-10
  )
})

# ---- Cross-package parity: responsePatterns::rp.acors (1e-10, complete) -----

test_that("cier_autocorrelation matches rp.acors on complete data (na_rm FALSE)", {
  skip_if_not_installed("responsePatterns")
  x <- ac_fixture(n = 20L, p = 30L, seed = 7L)
  # rp.acors()'s own default max_lag is ncol - 3 (cier caps at 10), so the parity
  # pins the shared statistic at rp.acors's range explicitly.
  ours <- suppressWarnings(cier_autocorrelation(x, max_lag = ncol(x) - 3L))$value
  rp <- suppressMessages(
    responsePatterns::rp.acors(as.data.frame(x), na.rm = FALSE)
  )
  theirs <- methods::slot(rp, "indices")$max.abs.ac
  expect_equal(ours, theirs, tolerance = 1e-10)
})

# ---- Convention pins: zero-variance slice scores 1 (NOT abstention) ---------

test_that("a straightliner scores exactly 1 (zero-variance slice -> 1, never NA)", {
  # The value convention: a zero-variance row has maximal autocorrelation 1, not NA.
  # A lone respondent cannot resolve a percentile cutoff, so it is not flagged here;
  # flagging once the cutoff is resolvable is pinned by the direction test below.
  expect_warning(out <- cier_autocorrelation(matrix(rep(3, 30L), nrow = 1L)),
                 class = "cier_warning_insufficient_items")
  expect_identical(out$value, 1)
  expect_false(is.na(out$value))
  expect_false(out$flagged[[1L]])
})

test_that("seesaw and diagonal patterns score 1 (max |ac|)", {
  pf <- pattern_fixture(p = 30L)
  v <- suppressWarnings(cier_autocorrelation(pf))$value
  expect_identical(v[[1L]], 1)               # constant: zero-var branch, exact 1
  expect_equal(v[[2L]], 1, tolerance = 1e-10)  # seesaw 1,5,1,5,... (lag-2 cor 1)
  expect_equal(v[[3L]], 1, tolerance = 1e-10)  # diagonal 1..5 repeated (lag-5 cor 1)
})

test_that("values are bounded in [0, 1]; a moderate max_lag de-saturates", {
  # At the default max_lag = p - 3 the short high-lag slices saturate |ac| to 1 on
  # discrete Likert data (the paper recommends lag <= 12); a moderate max_lag
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
  out <- suppressWarnings(cier_autocorrelation(mat, na_rm = TRUE))
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
  out <- suppressWarnings(cier_autocorrelation(mat))
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

# ---- per-lag minimum is 3 complete pairs (dropouts abstain, no spurious 1) ---

test_that("a sparse dropout row abstains (every lag has < 3 complete pairs)", {
  # Regression: under pairwise NA handling a respondent answering only a few
  # scattered items has no lag with 3 complete pairs, so every lag abstains and the
  # row scores NA -- it must NOT be sent to 1.0 by a 2-pair slice. Six of 44 items
  # answered at {1, 2, 20, 21, 40, 43}: the largest complete-pair count over lags
  # 1..10 is 2 (lag 1), so under the 3-pair minimum every lag abstains.
  dropout <- rep(NA_real_, 44L)
  dropout[c(1L, 2L, 20L, 21L, 40L, 43L)] <- c(2, 5, 1, 4, 3, 5)   # varied, not const
  # Embed among healthy rows so the percentile cutoff still resolves (>= 20 scored)
  # and the dropout's non-flag is meaningful.
  healthy <- ac_fixture(n = 25L, p = 44L, seed = 14L)
  m <- rbind(dropout, healthy)
  out <- suppressWarnings(cier_autocorrelation(m, max_lag = 10L))
  expect_true(is.na(out$value[[1L]]))            # abstains, not 1.0
  expect_true(is.na(out$flagged[[1L]]))
  expect_false(anyNA(out$value[-1L]))            # healthy rows still score
  expect_equal(out$value, ref_autocorrelation_value(m, max_lag = 10L),
               tolerance = 1e-10)
})

test_that("a 2-complete-pair lag no longer saturates the score to 1", {
  # Regression: an early-dropout row answering the first 6 of 44 items has exactly
  # 2 complete pairs at lag 4 ((1,5) and (2,6)), whose correlation is +/-1 by
  # construction -- the spurious saturation. Lags 1-3 keep 5/4/3 complete pairs, so
  # the row still scores, but below 1.
  dropout <- rep(NA_real_, 44L)
  dropout[1:6] <- c(1, 4, 2, 5, 3, 1)
  naive_lag4 <- stats::cor(c(dropout[[1L]], dropout[[2L]]),
                           c(dropout[[5L]], dropout[[6L]]))
  expect_equal(abs(naive_lag4), 1)               # the 2-pair slice IS +/-1
  v <- suppressWarnings(
    cier_autocorrelation(matrix(dropout, nrow = 1L), max_lag = 10L)
  )$value[[1L]]
  expect_false(is.na(v))                         # lags 1-3 still score it
  expect_true(v < 1)                             # but the spurious 1.0 is gone
  expect_equal(
    v,
    ref_autocorrelation_value(matrix(dropout, nrow = 1L), max_lag = 10L)[[1L]],
    tolerance = 1e-10
  )
})

test_that("the 3-pair minimum leaves the straightliner zero-variance convention intact", {
  # Zero-variance precedence: a lag slice constant over its non-NA cells still scores
  # 1 even with only 2 complete pairs -- the straightliner convention wins over the
  # < 3-pairs guard. The lag-1 slice of c(5, 5, NA, 5, 7) has row1 = (5, 5, NA, 5)
  # (var 0) and 2 complete pairs ((1) and (4)); a "pairs-guard-first" mutant returns
  # NA.
  zv <- matrix(c(5, 5, NA, 5, 7), nrow = 1L)
  v <- suppressWarnings(cier_autocorrelation(zv, max_lag = 1L))$value[[1L]]
  expect_identical(v, 1)                            # zero-variance wins, NOT NA
  expect_equal(v, ref_autocorrelation_value(zv, max_lag = 1L)[[1L]],
               tolerance = 1e-12)
})

# ---- the default lag window controls the bundled-data flag rate -------------

test_that("the default max_lag flags ~5% of bfi_careless, not ~32%", {
  # Regression: at the old default ncol - 3 the short high-lag slices saturate |ac|
  # to 1, flagging 126/394 (32%) and tripping the saturation diagnostic; the new
  # default min(ncol - 3, 10) flags 20/394 (~5%) with a representative cutoff and no
  # saturation warning.
  bfi <- as.matrix(bfi_careless[, 1:44])
  storage.mode(bfi) <- "double"

  expect_no_warning(default <- cier_autocorrelation(bfi),  # max_lag = min(41, 10)
                    class = "cier_warning_saturated_cutoff")
  expect_equal(default$cutoff, 0.9061, tolerance = 1e-3)
  expect_identical(sum(default$flagged, na.rm = TRUE), 20L)        # ~5% of 394

  # Assign INSIDE expect_warning: it returns the captured condition, not value.
  expect_warning(old <- cier_autocorrelation(bfi, max_lag = 44L - 3L),
                 class = "cier_warning_saturated_cutoff")
  expect_identical(sum(old$flagged, na.rm = TRUE), 126L)          # ~32% of 394
})

# ---- Lag-window edges -------------------------------------------------------

test_that("p = 4 works with the default max_lag (= 1) and matches the oracle", {
  x <- ac_fixture(n = 6L, p = 4L, seed = 3L)
  out <- suppressWarnings(cier_autocorrelation(x))
  expect_s3_class(out, "cier_index")
  expect_identical(length(out$value), 6L)
  # The narrowest lag window (a single lag-1 slice over 4 columns) is value-verified,
  # not just shape-checked: a wrong single-lag slice or a default off-by-one would
  # change these numbers.
  expect_equal(out$value, ref_autocorrelation_value(x), tolerance = 1e-10)
})

test_that("fewer than four items errors on the item count, not max_lag", {
  # A battery with no usable lag must blame the item count, not a max_lag the user
  # never supplied (a lag-1 slice of < 4 columns holds < 3 responses). The friendly
  # guard fires before the default max_lag resolves.
  x <- ac_fixture(n = 5L, p = 3L)
  err <- tryCatch(cier_autocorrelation(x),
                  cier_error_input = function(e) e)
  expect_s3_class(err, "cier_error_input")
  msg <- cli::ansi_strip(conditionMessage(err))
  expect_match(msg, "at least 4 items", fixed = TRUE)
  expect_no_match(msg, "max_lag", fixed = TRUE)         # never blames the knob
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
  expect_s3_class(suppressWarnings(cier_autocorrelation(x, max_lag = 7L)),
                  "cier_index")
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
  out <- suppressWarnings(cier_autocorrelation(x))
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
  # saturated ceiling of 1 the default produces on discrete data.
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
