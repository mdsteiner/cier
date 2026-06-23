# Tests for cier_gnormed() -- normed polytomous Guttman-error person-fit (Emons
# 2008; Molenaar 1991), computed in-package by the closed-form scorer
# (gnormed_scores) with an in-package nonparametric Monte-Carlo null cutoff.
#
# Trust model (INVERTED -- cier is now the implementation, not a PerFit bridge): the
# scorer is the in-package closed form, so the PRIMARY check is the INDEPENDENT
# closed-form oracle (ref_personfit_gnormed_poly: popularity-rank numerator +
# max-plus-knapsack denominator, never calling the kernel) at 1e-12. The PerFit
# cross-package parity confirms our EXACT score rounds to PerFit::Gnormed.poly's 4-dp
# output (round(got, 4) == PerFit, bytewise) and reduces to the dichotomous
# PerFit::Gnormed at Ncat = 2. The cutoff is the in-package Monte-Carlo null: its
# MECHANISM is tested directly (sum-score-conditional generation, bootstrap tail
# logic, seeded reproducibility), and its value is checked to AGREE with
# PerFit::cutoff to Monte-Carlo tolerance (not exact equality -- different RNG
# streams). Because Gnormed no longer needs a backend, the oracle-parity and edge
# tests run with or without PerFit; only the cross-package parity and the
# Monte-Carlo-agreement checks skip when PerFit is absent.

source(test_path("..", "reference", "ref-personfit-niessen-2016.R"))

# poly_matrix() / poly_items() fixtures are shared with test-cier-ht.R via
# helper-personfit.R.

# Recode 1..ncat -> 0..(ncat-1) for a hand-built PerFit call.
zero_base <- function(m) {
  z <- m - 1L
  storage.mode(z) <- "integer"
  z
}

# ---- Schema -----------------------------------------------------------------

test_that("cier_gnormed returns the pinned cier_index schema", {
  expect_cier_index_schema(cier_gnormed(poly_matrix(), poly_items(), seed = 1L),
                           "cier_gnormed", "upper", 60L)
})

# ---- Independent oracle parity (1e-12) --------------------------------------
# The PRIMARY check now that cier owns the implementation: the kernel value must
# equal the independent closed-form oracle. PerFit is not needed for these (the
# closed form is exact, where PerFit rounds to 4 dp).

test_that("cier_gnormed$value equals the oracle on a complete matrix", {
  m <- poly_matrix(n = 60L, p = 12L, seed = 21L)
  items <- poly_items(12L)
  got <- cier_gnormed(m, items, cutoff = 0.5)$value
  ref <- ref_personfit_gnormed_poly(list(responses = m, items = items))
  expect_equal(got, ref, tolerance = 1e-12)
})

test_that("cier_gnormed$value equals the oracle when rows carry NAs", {
  m <- poly_matrix(n = 60L, p = 12L, seed = 24L)
  m[2L, 3L] <- NA                   # one missing cell -> respondent excluded
  m[10L, ] <- NA                    # all-NA respondent -> excluded
  items <- poly_items(12L)
  got <- cier_gnormed(m, items, cutoff = 0.5)$value
  ref <- ref_personfit_gnormed_poly(list(responses = m, items = items))
  expect_equal(got, ref, tolerance = 1e-12)
  expect_true(is.na(got[[2L]]))
  expect_true(is.na(got[[10L]]))
})

test_that("cier_gnormed reverse-scores keyed items (oracle parity, not a no-op)", {
  m <- poly_matrix(n = 60L, p = 12L, seed = 22L)
  rk <- rep(c(FALSE, TRUE), 6L)
  items <- poly_items(12L, reverse = rk)
  got <- cier_gnormed(m, items, cutoff = 0.5)$value
  ref <- ref_personfit_gnormed_poly(list(responses = m, items = items))
  expect_equal(got, ref, tolerance = 1e-12)
  # Reverse-scoring is not a no-op: ignoring the keying gives a different score.
  raw <- ref_personfit_gnormed_poly(list(responses = m, items = poly_items(12L)))
  expect_false(isTRUE(all.equal(got, raw)))
})

# ---- Cross-package parity: PerFit preprocessing + dichotomous reduction ------
# Our exact value rounds to PerFit's 4-dp output bytewise (the inverted trust model:
# PerFit is the cross-package oracle, no longer the runtime scorer).

test_that("cier_gnormed rounds to a hand-built PerFit::Gnormed.poly call (bytewise)", {
  # Pins the bridge preprocessing (zero-basing to 0..(Ncat-1), persons-as-rows) and
  # that the exact closed form matches PerFit to its 4-dp output rounding. n != p so a
  # missing-transpose mutant returns the wrong length; a raw-1..k coding mutant shifts
  # every score. round(got, 4) == PerFit (already 4 dp) bytewise (tol 0).
  skip_if_not_installed("PerFit")
  m <- poly_matrix(n = 50L, p = 12L, seed = 7L)
  got <- cier_gnormed(m, poly_items(12L), cutoff = 0.5)$value
  invisible(utils::capture.output(
    fit <- PerFit::Gnormed.poly(matrix = zero_base(m), Ncat = 5L)
  ))
  ref <- as.numeric(fit$PFscores$PFscores)
  expect_equal(round(got, 4), ref, tolerance = 0)
})

test_that("cier_gnormed reduces to the dichotomous PerFit::Gnormed at Ncat = 2", {
  skip_if_not_installed("PerFit")
  b <- poly_matrix(n = 80L, p = 12L, ncat = 2L, seed = 11L)
  got <- cier_gnormed(b, poly_items(12L, ncat = 2L), cutoff = 0.5)$value
  invisible(utils::capture.output(
    dich <- as.numeric(PerFit::Gnormed(zero_base(b))$PFscores$PFscores)
  ))
  ok <- is.finite(got) & is.finite(dich)
  expect_gt(sum(ok), 10L)
  # Our exact poly value at Ncat = 2 equals the dichotomous Gnormed; both round to the
  # same 4 dp, so round(got, 4) == dich bytewise.
  expect_lt(max(abs(round(got[ok], 4) - dich[ok])), 1e-9)
})

# ---- Cutoff: in-package Monte-Carlo null -------------------------------------

test_that("the default cutoff agrees with PerFit::cutoff to Monte-Carlo tolerance", {
  # The engine reproduces PerFit::cutoff's NonParametric default mechanism. The two
  # draw from different RNG streams, so individual cutoffs are NOT equal; their MEANS
  # over several seeds agree within Monte-Carlo error (the per-call sd is ~0.01). A
  # wrong tail shifts the mean far beyond tolerance; the sum-score conditioning is
  # pinned separately by the mechanism test below (an unconditional resample drifts
  # the mean only ~0.01 here, within tolerance, so this test does not guard it).
  skip_if_not_installed("PerFit")
  m <- poly_matrix(n = 60L, p = 12L, seed = 5L)
  items <- poly_items(12L)
  invisible(utils::capture.output(
    fit <- PerFit::Gnormed.poly(matrix = zero_base(m), Ncat = 5L)
  ))
  seeds <- 1:8
  ours <- vapply(seeds, function(s) cier_gnormed(m, items, seed = s)$cutoff,
                 numeric(1L))
  perf <- vapply(seeds, function(s) {
    set.seed(s)
    invisible(utils::capture.output(co <- PerFit::cutoff(fit, Blvl = 0.05)))
    as.numeric(co$Cutoff)
  }, numeric(1L))
  expect_lt(abs(mean(ours) - mean(perf)), 0.02)
  # Not a degenerate constant, and a plausible normed-Guttman cutoff.
  expect_gt(stats::sd(ours), 0)
  expect_true(all(ours > 0 & ours < 1))
})

test_that("a larger fpr lowers the upper-tail Monte-Carlo cutoff (same seed)", {
  # The single direction flip: Gnormed is UPPER-tail, so the cutoff is the 1 - fpr
  # quantile. With the null generation and bootstrap resampling seed-fixed and
  # independent of fpr, a larger fpr can only LOWER the cutoff. A lower-tail mutant
  # (Blvl.use = fpr) would invert this (c at 0.10 > c at 0.05).
  m <- poly_matrix(n = 60L, p = 12L, seed = 5L)
  items <- poly_items(12L)
  c05 <- cier_gnormed(m, items, fpr = 0.05, seed = 7L)$cutoff
  c10 <- cier_gnormed(m, items, fpr = 0.10, seed = 7L)$cutoff
  expect_gt(c05, c10)
})

test_that("a seeded null cutoff is reproducible and leaves the ambient RNG intact", {
  m <- poly_matrix(n = 50L, p = 10L, seed = 3L)
  items <- poly_items(10L)
  set.seed(999L)
  before <- get(".Random.seed", envir = globalenv())
  c1 <- cier_gnormed(m, items, seed = 42L)$cutoff
  after <- get(".Random.seed", envir = globalenv())
  c2 <- cier_gnormed(m, items, seed = 42L)$cutoff
  expect_identical(c1, c2)             # same seed -> identical cutoff
  expect_identical(after, before)      # seeded call restored the caller's RNG
})

test_that("a NULL seed draws from (and advances) the ambient stream", {
  # The NULL-seed path must NOT save/restore: it consumes the caller's RNG (the
  # complement of the seeded path's restore contract). A mutant that always
  # save/restores, or returns a fixed cutoff, is caught by the advanced-stream
  # assertion.
  m <- poly_matrix(n = 50L, p = 10L, seed = 3L)
  items <- poly_items(10L)
  set.seed(1L)
  before <- get(".Random.seed", envir = globalenv())
  a <- cier_gnormed(m, items)$cutoff
  after <- get(".Random.seed", envir = globalenv())
  expect_true(is.finite(a))
  expect_false(identical(before, after))   # ambient stream advanced (no restore)
})

# ---- Monte-Carlo null engine (mechanism) ------------------------------------
# The cutoff value is stochastic, so the mechanism is pinned directly: the
# sum-score-conditional generator, the bootstrap tail logic, and the composer's
# determinism. These reuse the SHARED engine (personfit_null_matrix /
# bootstrap_tail_cutoff / personfit_null_cutoff) that a future Ht null also reuses.

test_that("the null generator conditions on the sum score (not independent columns)", {
  # All-min and all-max rows: at sum-score level 0 every item is 0, at the top level
  # every item is its max. A SUM-SCORE-CONDITIONAL draw can then only produce constant
  # rows (all 0 or all p). An unconditional column resample (the trap the design
  # warns against) would mix categories within a row and destroy the inter-item
  # covariance the statistic measures.
  p <- 6L
  z <- rbind(matrix(0L, 50L, p), matrix(1L, 50L, p))   # 0-based levels 0 and p
  null <- with_local_seed(1L,
                          function() personfit_null_matrix(z, ncat = 2L, nreps = 500L))
  expect_identical(dim(null), c(500L, p))
  rs <- rowSums(null)
  expect_true(all(rs %in% c(0L, p)))                   # every null row is constant
  expect_true(any(rs == 0L) && any(rs == p))           # both observed levels drawn
})

test_that("the null generator draws each item from its per-level category frequencies", {
  # One sum-score level (every row sums to 1) with a known per-item category mix: item
  # 1 is 0 for 80 rows and 1 for 20, balanced by item 2 so the row sum is constant.
  # The generated null's item-1 proportion must track the empirical 0.2 (a uniform
  # mutant would give ~0.5).
  n0 <- 80L
  n1 <- 20L
  z <- rbind(cbind(rep(0L, n0), rep(1L, n0)),          # item1=0, item2=1 (sum 1)
             cbind(rep(1L, n1), rep(0L, n1)))          # item1=1, item2=0 (sum 1)
  expect_identical(length(unique(rowSums(z))), 1L)     # a single sum-score level
  null <- with_local_seed(2L,
                          function() personfit_null_matrix(z, ncat = 2L, nreps = 4000L))
  expect_equal(mean(null[, 1L]), 0.2, tolerance = 0.05)
  expect_equal(mean(null[, 2L]), 0.8, tolerance = 0.05)
})

test_that("the bootstrap tail cutoff flips with the tail direction", {
  # On a fixed null vector the upper tail takes the 1 - fpr quantile, the lower tail
  # the fpr quantile. With fpr = 0.05 on 1..100 these straddle the distribution: the
  # upper cutoff sits near 95, the lower near 5.
  pfs <- as.numeric(1:100)
  up <- with_local_seed(3L, function() bootstrap_tail_cutoff(pfs, "upper", 0.05))
  lo <- with_local_seed(3L, function() bootstrap_tail_cutoff(pfs, "lower", 0.05))
  expect_gt(up, 80)
  expect_lt(lo, 20)
  expect_gt(up, lo)
})

test_that("bootstrap tail cutoff chunking preserves seeded output", {
  pfs <- as.numeric(1:100)
  unchunked <- with_local_seed(4L, function() {
    bootstrap_tail_cutoff(pfs, "upper", 0.05, breps = 101L,
                          chunk_cells = 1e6)
  })
  chunked <- with_local_seed(4L, function() {
    bootstrap_tail_cutoff(pfs, "upper", 0.05, breps = 101L,
                          chunk_cells = 300L)
  })
  expect_identical(chunked, unchunked)
})

test_that("personfit_null_cutoff is seed-reproducible and scorer/tail driven", {
  # End-to-end engine determinism on a fixed block: same seed -> identical cutoff, and
  # the Gnormed (upper) null cutoff is a finite normed-Guttman value in [0, 1].
  m <- poly_matrix(n = 60L, p = 12L, seed = 5L)
  z <- zero_base(m)
  score <- function(null) gnormed_scores(null, 5L)
  run <- function() {
    personfit_null_cutoff(z, 5L, score, tail = "upper", perfect = "allowed",
                          fpr = 0.05)
  }
  c1 <- with_local_seed(9L, run)
  c2 <- with_local_seed(9L, run)
  expect_identical(c1, c2)
  expect_true(is.finite(c1) && c1 >= 0 && c1 <= 1)
})

test_that("the engine's perfect='excluded' branch runs on polytomous data", {
  # Dormant today (Gnormed is perfect-allowed); exercised here so the constant-row
  # exclusion (matrixStats::rowVars > 0, NOT PerFit's dichotomous rowSums %% p rule,
  # which would wrongly drop non-constant polytomous rows) is locked in for the future
  # Ht-null reuse. On 5-point data it must run without error -- confirming rowVars
  # handles the integer null and the branch does not over-prune to an empty matrix --
  # and return a finite normed-Guttman cutoff.
  m <- poly_matrix(n = 60L, p = 12L, seed = 5L)
  z <- zero_base(m)
  score <- function(null) gnormed_scores(null, 5L)
  cut <- with_local_seed(11L, function() {
    personfit_null_cutoff(z, 5L, score, tail = "upper", perfect = "excluded",
                          fpr = 0.05)
  })
  expect_true(is.finite(cut) && cut >= 0 && cut <= 1)
})

# ---- Direction (upper) ------------------------------------------------------

test_that("direction is upper: a high-Guttman-error row flags, a low one does not", {
  # Deterministic comparator with a literal mid-cutoff (independent of the null):
  # the highest score flags, the lowest does not, flag rule is value >= cutoff.
  # A flag-lower mutant inverts both.
  m <- poly_matrix(n = 60L, p = 12L, seed = 11L)
  items <- poly_items(12L)
  v <- cier_gnormed(m, items, cutoff = 0.5)$value
  mid <- stats::median(v)
  out <- cier_gnormed(m, items, cutoff = mid)
  expect_true(out$flagged[[which.max(out$value)]])
  expect_false(out$flagged[[which.min(out$value)]])
  expect_identical(out$flagged, !is.na(out$value) & out$value >= out$cutoff)
})

# ---- Edge cases -------------------------------------------------------------

test_that("fewer than three items abstains (recursive denominator needs >= 3)", {
  m <- poly_matrix(n = 40L, p = 2L, seed = 26L)
  expect_warning(out <- cier_gnormed(m, poly_items(2L)),
                 class = "cier_warning_insufficient_items")
  expect_true(all(is.na(out$value)))
  expect_true(is.na(out$cutoff))
  expect_true(all(is.na(out$flagged)))
  # Length pins the full-length abstain vectors -- `all(is.na(numeric(0)))` is
  # vacuously TRUE, so a zero-length / scalar result would slip past; the index
  # must stay one row per respondent even when it abstains.
  expect_identical(length(out$value), 40L)
  expect_identical(length(out$flagged), 40L)
})

test_that("exactly three items scores (the items knife-edge, scoring side)", {
  # Brackets the >= 3 boundary from below: p = 2 abstains (above), p = 3 must
  # score, catching an off-by-one (abstain at p <= 3). Parity (not just
  # finiteness) pins value correctness at the minimum item count, where a small-p
  # preprocessing/orientation bug could hide.
  m <- poly_matrix(n = 40L, p = 3L, seed = 28L)
  items <- poly_items(3L)
  out <- cier_gnormed(m, items, cutoff = 0.5)
  expect_true(any(is.finite(out$value)))
  ref <- ref_personfit_gnormed_poly(list(responses = m, items = items))
  expect_equal(out$value, ref, tolerance = 1e-12)
})

test_that("fewer than two complete respondents abstains", {
  m <- poly_matrix(n = 5L, p = 6L, seed = 27L)
  m[1:4, 1L] <- NA                  # only respondent 5 is complete
  expect_warning(out <- cier_gnormed(m, poly_items(6L)),
                 class = "cier_warning_insufficient_items")
  expect_true(all(is.na(out$value)))
  expect_identical(length(out$value), 5L)
  expect_identical(length(out$flagged), 5L)
})

test_that("exactly two complete respondents scores (the complete knife-edge)", {
  # Brackets the >= 2-complete boundary from below: 1 complete abstains (above),
  # 2 complete must score, catching a `< 3 complete` off-by-one. The two complete
  # rows are forced to span the 1..5 scale (the closed form needs both extreme
  # categories); the count threshold, not the span, is what this pins.
  m <- poly_matrix(n = 6L, p = 6L, seed = 29L)
  m[1L, 1L] <- 1                    # row 1 alone spans the 1..5 scale (min..max)
  m[1L, 2L] <- 5
  m[3:6, 1L] <- NA                  # only respondents 1 and 2 are complete
  items <- poly_items(6L)
  out <- cier_gnormed(m, items, cutoff = 0.5)
  expect_true(is.finite(out$value[[1L]]))
  expect_true(is.finite(out$value[[2L]]))
  expect_true(all(is.na(out$value[3:6])))
  # Parity pins the two scored values and the NA placement against the oracle at
  # the 2-complete boundary.
  ref <- ref_personfit_gnormed_poly(list(responses = m, items = items))
  expect_equal(out$value, ref, tolerance = 1e-12)
})

test_that("a single-respondent matrix abstains", {
  m <- poly_matrix(n = 1L, p = 5L, seed = 1L)
  expect_warning(out <- cier_gnormed(m, poly_items(5L)),
                 class = "cier_warning_insufficient_items")
  expect_true(all(is.na(out$value)))
  expect_identical(length(out$value), 1L)
  expect_identical(length(out$flagged), 1L)
})

test_that("a straightliner gets a finite (low) score, not an abstention", {
  # Documented blind spot: unlike Ht, the Guttman-error index scores a
  # straightliner -- straightlining evades detection rather than abstaining.
  m <- poly_matrix(n = 40L, p = 10L, seed = 23L)
  m[1L, ] <- 3                      # midpoint straightline (complete row)
  out <- cier_gnormed(m, poly_items(10L), cutoff = 0.5)
  expect_true(is.finite(out$value[[1L]]))
  expect_true(all(is.finite(out$value)))
})

test_that("a respondent with one missing cell abstains; rows stay aligned", {
  m <- poly_matrix(n = 40L, p = 10L, seed = 24L)
  m[2L, 3L] <- NA
  out <- cier_gnormed(m, poly_items(10L), cutoff = 0.5)
  expect_true(is.na(out$value[[2L]]))
  expect_true(is.na(out$flagged[[2L]]))
  expect_gt(sum(is.finite(out$value)), 2L)
})

test_that("a fully-constant item column matches the oracle (popularity-rank edge)", {
  # A zero-variance item ties its popularity rank; the popularity-rank numerator
  # and the knapsack denominator must still agree with the oracle. A constant
  # COLUMN, distinct from the straightliner ROW above.
  m <- poly_matrix(n = 50L, p = 10L, seed = 30L)
  m[, 4L] <- 3                      # one constant (straightlined) item column
  items <- poly_items(10L)
  got <- cier_gnormed(m, items, cutoff = 0.5)$value
  ref <- ref_personfit_gnormed_poly(list(responses = m, items = items))
  expect_equal(got, ref, tolerance = 1e-12)
})

test_that("the scale base (items$min) is honoured for keying and zero-basing", {
  # Min-invariance: the same responses in 1..ncat coding (min = 1, default) and in
  # 0..(ncat-1) coding (min = 0) must score identically, with a reverse-keyed item
  # present so the (min + max) - x reflection is exercised, not just the zero-base.
  # A bridge hardcoding `- 1` zero-basing or `(max + 1) - x` reflection (ignoring
  # min), or deriving Ncat as `max` not `max - min + 1` (the 0-based declaration
  # has max = 4 but five categories), diverges or errors on the 0-based scale. The
  # oracle hardcodes base 1, so this is an internal invariance.
  rk <- rep(c(FALSE, TRUE), 5L)
  m1 <- poly_matrix(n = 50L, p = 10L, ncat = 5L, seed = 31L)   # 1..5 coding
  items1 <- poly_items(10L, reverse = rk)                       # min defaults 1
  m0 <- m1 - 1L                                                 # 0..4 coding
  items0 <- data.frame(reverse_keyed = rk, max = 4L, min = 0L)
  expect_equal(cier_gnormed(m1, items1, cutoff = 0.5)$value,
               cier_gnormed(m0, items0, cutoff = 0.5)$value,
               tolerance = 0)
})

test_that("items with equal span but different bases score together (per-item base)", {
  # The single-Ncat contract is about the NUMBER of response options (max - min + 1),
  # not identical min/max pairs: a battery mixing 1..5 and 0..4 items (both five
  # options) is valid and must equal the all-1..5 scoring of the same responses. A
  # homogeneity check on `max` alone (or a global rather than per-item zero-base)
  # would reject or mis-score it.
  m1 <- poly_matrix(n = 50L, p = 6L, ncat = 5L, seed = 41L)     # all 1..5
  m_mixed <- m1
  m_mixed[, 4:6] <- m_mixed[, 4:6] - 1L                         # cols 4-6: 0..4
  items_mixed <- data.frame(reverse_keyed = FALSE,
                            max = c(5L, 5L, 5L, 4L, 4L, 4L),
                            min = c(1L, 1L, 1L, 0L, 0L, 0L))
  expect_equal(cier_gnormed(m_mixed, items_mixed, cutoff = 0.5)$value,
               cier_gnormed(m1, poly_items(6L), cutoff = 0.5)$value,
               tolerance = 0)
})

test_that("a two-option 0/1 declaration (min = 0, max = 1) scores (dichotomous base)", {
  # The smallest valid scale through the personfit validator: max >= min + 1 is the
  # bound, so a validator demanding max >= 2 regardless of min would wrongly reject
  # this. Equality with the same data in 1..2 coding pins the min/max plumbing at
  # Ncat = 2.
  m2 <- poly_matrix(n = 40L, p = 6L, ncat = 2L, seed = 43L)   # 1..2 coding
  m2[1L, 1L] <- 1                  # force both extremes somewhere in the block
  m2[2L, 1L] <- 2
  m01 <- m2 - 1L                                              # 0/1 coding
  items01 <- data.frame(reverse_keyed = rep(FALSE, 6L), max = 1L, min = 0L)
  expect_equal(cier_gnormed(m01, items01, cutoff = 0.5)$value,
               cier_gnormed(m2, poly_items(6L, ncat = 2L), cutoff = 0.5)$value,
               tolerance = 0)
})

test_that("a response outside the declared scale is a plain input error, not a backend limit", {
  # The bridge zero-bases to 0..(Ncat - 1) and range-checks; a code above `max`
  # (or below `min`) is a contract violation, not silently coerced. check_responses
  # only rejects NaN/Inf, so this pins the bridge's own range guard. It is a PLAIN
  # cier_error_input, NOT a backend limit: an out-of-range value is a genuine data
  # defect that must keep propagating through cier_screen(), so the backend-limit
  # subclass (reserved for otherwise-valid data the kernel cannot score) must NOT
  # attach here -- a mutant tagging both zero-base branches alike would let the
  # screen swallow corrupt data.
  m <- poly_matrix(n = 10L, p = 6L, seed = 32L)
  m[1L, 1L] <- 6                    # exceeds max = 5
  err <- tryCatch(cier_gnormed(m, poly_items(6L)), error = function(e) e)
  expect_s3_class(err, "cier_error_input")
  expect_false(inherits(err, "cier_error_backend_limit"))
})

test_that("data that does not span the declared scale is a screen-survivable backend limit", {
  # The item-step popularities are undefined when a declared extreme category never
  # occurs (min == 0, max == Ncat - 1 after zero-basing). Responses confined to 2..4
  # of a declared 1..5 scale never reach the extremes, so the kernel surfaces a typed
  # error. This is OTHERWISE-VALID sample-dependent data the closed form refuses to
  # score (the contract inherited from PerFit; like the heterogeneous-span case), NOT
  # a metadata defect -- so the abort carries the cier_error_backend_limit subclass
  # with a compact data$reason, and cier_screen() records Gnormed as
  # skipped-with-reason. It stays cier_error_input for direct callers.
  withr::with_seed(33L, {
    m <- matrix(sample(2:4, 30L * 6L, replace = TRUE), nrow = 30L, ncol = 6L)
  })
  storage.mode(m) <- "double"
  err <- tryCatch(cier_gnormed(m, poly_items(6L)), error = function(e) e)
  expect_s3_class(err, "cier_error_input")
  expect_s3_class(err, "cier_error_backend_limit")
  # The exact compact reason is the screen's skip text (read via data$reason);
  # pin it bytewise so a mutant with an empty / drifted reason is caught.
  expect_identical(
    cier_condition_data(err)$reason,
    "sample does not attain both scale extremes (every end category must be observed)"
  )
})

# ---- abstention / abort messages name the real cause ------------------------

# Flatten a typed condition's full cli message (header + x/i bullets) to a
# whitespace-normalised, ANSI-stripped string for styling-agnostic substring
# matching (the idiom from test-cier-mahalanobis.R). The count is interpolated
# plainly (no {.val} quoting), so "only 2 items present" matches verbatim.
flatten_cnd_message <- function(cond) {
  gsub("\\s+", " ", cli::ansi_strip(rlang::cnd_message(cond)))
}

test_that("the too-few-items abstention names items, not 'no respondent'", {
  # The wholesale-abstain warning must distinguish its two causes. With p < 3 but
  # many complete respondents, the old "no respondent could be scored" / n_used = 0L
  # was false: nothing is wrong with the respondents, there are just too few items.
  # Pin the items wording (plural and singular) and the honest item-count payload,
  # and assert the complete-respondents phrasing is NOT used here (a mutant emitting
  # one message for both causes, or hardcoding the plural noun, dies here).
  m <- poly_matrix(n = 40L, p = 2L, seed = 26L)   # 40 complete respondents, 2 items
  w <- expect_warning(cier_gnormed(m, poly_items(2L)),
                      class = "cier_warning_insufficient_items")
  msg <- flatten_cnd_message(w)
  expect_match(msg, "only 2 items present", fixed = TRUE)
  expect_match(msg, "needs at least three items", fixed = TRUE)
  expect_no_match(msg, "complete respondent", fixed = TRUE)
  expect_no_match(msg, "no respondent could be scored", fixed = TRUE)
  # Honest item count; the misleading n_used = 0L must be gone.
  expect_identical(cier_condition_data(w)$n_items, 2L)
  expect_null(cier_condition_data(w)$n_used)
  # Singular form: a single item must read "1 item present", not "1 items".
  m1 <- poly_matrix(n = 40L, p = 1L, seed = 26L)
  w1 <- expect_warning(cier_gnormed(m1, poly_items(1L)),
                       class = "cier_warning_insufficient_items")
  expect_match(flatten_cnd_message(w1), "only 1 item present", fixed = TRUE)
  expect_identical(cier_condition_data(w1)$n_items, 1L)
})

test_that("the too-few-complete abstention names complete respondents, honest count", {
  # The complementary cause: enough items (p >= 3) but only one complete respondent.
  # The message must name complete respondents and carry the ACTUAL count (1), not
  # the hard-coded 0L the old code always reported.
  m <- poly_matrix(n = 5L, p = 6L, seed = 27L)
  m[1:4, 1L] <- NA                   # only respondent 5 is complete
  w <- expect_warning(cier_gnormed(m, poly_items(6L)),
                      class = "cier_warning_insufficient_items")
  msg <- flatten_cnd_message(w)
  expect_match(msg, "only 1 complete respondent", fixed = TRUE)
  expect_match(msg, "at least two complete respondents", fixed = TRUE)
  expect_no_match(msg, "items present", fixed = TRUE)
  expect_no_match(msg, "no respondent could be scored", fixed = TRUE)
  expect_identical(cier_condition_data(w)$n_used, 1L)   # honest count, not 0L
  expect_null(cier_condition_data(w)$n_items)
})

test_that("the too-few-complete count is the real complete-case count, even at zero", {
  # When no respondent is complete the count is genuinely 0 -- but it must be the
  # COMPUTED count (pluralised "respondents"), not the old hard-coded 0L sentinel.
  # With the n = 1 case above this pins the payload as the real sum(complete).
  m <- poly_matrix(n = 4L, p = 6L, seed = 50L)
  m[, 1L] <- NA                      # every respondent has a missing cell
  w <- expect_warning(cier_gnormed(m, poly_items(6L)),
                      class = "cier_warning_insufficient_items")
  msg <- flatten_cnd_message(w)
  expect_match(msg, "only 0 complete respondents", fixed = TRUE)
  expect_identical(cier_condition_data(w)$n_used, 0L)
})

test_that("the unattained-extreme abort names complete-case dropping as a cause", {
  # Beyond the unused-extreme contract, the abort adds an i bullet naming the likely
  # real cause: complete-case dropping removed the only rows holding the top/bottom
  # category. Pin the new bullet and re-confirm the screen's skip text (data$reason).
  withr::with_seed(33L, {
    m <- matrix(sample(2:4, 30L * 6L, replace = TRUE), nrow = 30L, ncol = 6L)
  })
  storage.mode(m) <- "double"
  err <- tryCatch(cier_gnormed(m, poly_items(6L)), error = function(e) e)
  msg <- flatten_cnd_message(err)
  expect_match(msg, "any missing cell are dropped before scoring", fixed = TRUE)
  expect_match(msg, "does not count toward the observed span", fixed = TRUE)
  expect_identical(
    cier_condition_data(err)$reason,
    "sample does not attain both scale extremes (every end category must be observed)"
  )
})

test_that("a fractional (non-integer) response is a typed input error", {
  # check_responses() only rejects NaN/Inf; a fractional cell (e.g. averaged or
  # imputed) would otherwise be silently truncated by the integer cast before
  # scoring. The bridge catches it as a typed error -- no silent coercion.
  m <- poly_matrix(n = 30L, p = 6L, seed = 34L)
  m[1L, 1L] <- 2.5
  expect_error(cier_gnormed(m, poly_items(6L)), class = "cier_error_input")
})

# ---- Input validation -------------------------------------------------------

test_that("a non-matrix / non-numeric / non-finite payload is a typed input error", {
  expect_error(cier_gnormed(1:10, poly_items(10L)), class = "cier_error_input")
  expect_error(cier_gnormed(matrix(letters[1:6], nrow = 2L), poly_items(3L)),
               class = "cier_error_input")
  bad <- poly_matrix(n = 6L, p = 6L)
  bad[1L, 1L] <- Inf
  expect_error(cier_gnormed(bad, poly_items(6L)), class = "cier_error_input")
})

test_that("malformed items (max) are typed input errors", {
  m <- poly_matrix(n = 10L, p = 6L)
  # max column absent
  expect_error(cier_gnormed(m, data.frame(reverse_keyed = rep(FALSE, 6L))),
               class = "cier_error_input")
  # NA max
  expect_error(cier_gnormed(m, data.frame(max = c(5, 5, 5, 5, 5, NA))),
               class = "cier_error_input")
  # fewer than two response options (max == min with the default min = 1)
  expect_error(cier_gnormed(m, data.frame(max = rep(1, 6L))),
               class = "cier_error_input")
  # wrong number of item rows
  expect_error(cier_gnormed(m, data.frame(max = rep(5, 3L))),
               class = "cier_error_input")
  # fractional and non-finite max are PLAIN input errors: although both also make
  # spans unequal, per-item validity is classified before homogeneity, so neither
  # may carry the backend-limit subclass (the screen must propagate, not skip).
  for (bad_max in list(c(5, 5, 5, 5, 5, 2.5), c(5, 5, 5, 5, 5, Inf))) {
    err <- tryCatch(cier_gnormed(m, data.frame(max = bad_max)),
                    error = function(e) e)
    expect_s3_class(err, "cier_error_input")
    expect_false(inherits(err, "cier_error_backend_limit"))
  }
})

test_that("a heterogeneous span is a typed backend limit, not a plain input error", {
  # Accurate metadata describing genuinely mixed-format data (here five-option and
  # four-option items together) is not a malformed frame: the single-Ncat contract
  # cannot score it. The abort carries the cier_error_backend_limit subclass (the same
  # line the unattained-scale-extreme case draws) so cier_screen() can skip-with-reason
  # instead of dying. It stays cier_error_input for direct callers.
  m <- poly_matrix(n = 10L, p = 6L)
  het <- data.frame(max = c(5, 5, 5, 5, 5, 4))
  expect_error(cier_gnormed(m, het), class = "cier_error_input")
  expect_error(cier_gnormed(m, het), class = "cier_error_backend_limit")
  # The limit is on the span (number of options), NOT on max itself: equal max with
  # differing min is just as heterogeneous (spans 4 vs 5)...
  het_min <- data.frame(max = 5, min = c(1, 1, 1, 1, 1, 0))
  expect_error(cier_gnormed(m, het_min), class = "cier_error_backend_limit")
  # ...while a malformed per-item max (NA) stays a PLAIN input error: the backend-
  # limit subclass must not swallow genuine metadata defects.
  bad <- tryCatch(cier_gnormed(m, data.frame(max = c(5, 5, 5, 5, 5, NA))),
                  error = function(e) e)
  expect_false(inherits(bad, "cier_error_backend_limit"))
})

test_that("a non-data.frame items or a non-integer min is a typed input error", {
  m <- poly_matrix(n = 10L, p = 6L)
  expect_error(cier_gnormed(m, "not a frame"), class = "cier_error_input")
  # `min` (the scale base) must be a finite whole number on every item.
  bad_min <- data.frame(max = 5L, min = c(1, 1, 1, 1, 1, 1.5))
  expect_error(cier_gnormed(m, bad_min), class = "cier_error_input")
})

test_that("invalid fpr / cutoff / seed values are typed input errors", {
  m <- poly_matrix(n = 10L, p = 6L)
  items <- poly_items(6L)
  expect_error(cier_gnormed(m, items, fpr = 0), class = "cier_error_input")
  expect_error(cier_gnormed(m, items, fpr = 1), class = "cier_error_input")
  expect_error(cier_gnormed(m, items, fpr = c(0.05, 0.1)),
               class = "cier_error_input")
  # A normed Guttman value is in [0, 1]; a threshold outside flags everyone or no
  # one, so it is rejected.
  expect_error(cier_gnormed(m, items, cutoff = -0.1), class = "cier_error_input")
  expect_error(cier_gnormed(m, items, cutoff = 1.5), class = "cier_error_input")
  expect_error(cier_gnormed(m, items, seed = 1.5), class = "cier_error_input")
})

test_that("supplying both fpr and cutoff is a typed input error", {
  m <- poly_matrix(n = 10L, p = 6L)
  expect_error(cier_gnormed(m, poly_items(6L), fpr = 0.1, cutoff = 0.5),
               class = "cier_error_input")
})

test_that("cier_gnormed scores without PerFit (no runtime backend gate)", {
  # Targeted regression guard for the removed runtime gate: Gnormed is pure R, so
  # it scores normally without consulting PerFit.
  m <- poly_matrix(n = 40L, p = 10L, seed = 26L)
  items <- poly_items(10L)
  out <- cier_gnormed(m, items, cutoff = 0.5)
  expect_s3_class(out, "cier_index")
  expect_true(any(is.finite(out$value)))
  ref <- ref_personfit_gnormed_poly(list(responses = m, items = items))
  expect_equal(out$value, ref, tolerance = 1e-12)
})

# ---- print snapshot (locked, design-first; direction = upper) ---------------

test_that("print renders the locked cli summary (upper direction)", {
  # A literal cutoff keeps the printed threshold stable across platforms (the null
  # would vary the number); the snapshot pins the print FORMAT.
  m <- poly_matrix(n = 30L, seed = 11L)
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    expect_snapshot(print(cier_gnormed(m, poly_items(), cutoff = 0.5)))
  })
})

test_that("print reports abstaining respondents as '(no score)'", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    m <- poly_matrix(n = 30L, p = 12L, seed = 11L)
    m[1L, 3L] <- NA                 # one missing-cell respondent -> abstains
    expect_snapshot(print(cier_gnormed(m, poly_items(12L), cutoff = 0.5)))
  })
})

# ---- missing reverse_keyed inform -------------------------------------------

test_that("cier_gnormed informs when items omits reverse_keyed", {
  # Known-good shape (spans 1..5 globally, so it scores rather than tripping the
  # unattained-extreme backend limit); only difference is absent reverse_keyed.
  m <- poly_matrix(n = 60L, p = 12L, seed = 11L)
  it <- data.frame(max = rep(5L, 12L))   # max on every item, but no reverse_keyed
  expect_message(suppressWarnings(cier_gnormed(m, it, seed = 1L)),
                 class = "cier_message_forward_keyed")
})
