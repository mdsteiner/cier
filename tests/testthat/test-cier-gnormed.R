# Tests for cier_gnormed() -- the normed polytomous Guttman-error person-fit
# index (Emons 2008; Molenaar 1991), scored via PerFit::Gnormed.poly with a
# PerFit Monte-Carlo null cutoff.
#
# Trust model: the Gnormed scorer IS PerFit::Gnormed.poly (single-kernel rule),
# so the genuine independent check is the closed-form oracle
# (ref_personfit_gnormed_poly: popularity-rank numerator + max-plus-knapsack
# denominator, from scratch, never calling the production bridge). PerFit rounds
# its scores to 4 decimals and the oracle is exact, so round(oracle, 4) matches
# bytewise (tolerance 0; see tests/reference/TOLERANCES.md). The PerFit
# direct-call parity exists to pin the BRIDGE PREPROCESSING -- zero-basing to
# 0..(Ncat-1), persons-as-rows orientation, reverse-keying, complete-casing --
# not the scorer (which is shared); an n != p fixture makes a missing-transpose
# mutant break. The cutoff is the PerFit Monte-Carlo null (PerFit::cutoff): fast
# and reproducible under a seed, pinned against a same-seed PerFit::cutoff call.

source(test_path("..", "reference", "ref-personfit-niessen-2016.R"))

# A reproducible polytomous response matrix in 1..ncat coding (the wrapper
# zero-bases internally). Uniform sampling gives a spread of Gnormed values
# (needed for the direction test).
poly_matrix <- function(n = 60L, p = 12L, ncat = 5L, seed = 21L) {
  withr::with_seed(seed, {
    m <- matrix(sample.int(ncat, n * p, replace = TRUE), nrow = n, ncol = p)
  })
  storage.mode(m) <- "double"
  m
}

# Item metadata for cier_gnormed(): homogeneous span (the fixtures are 1-based,
# 1..ncat, so `max = ncat`), optional `reverse_keyed`. The same data.frame
# doubles as the oracle's `data$items` (the oracle reads items$reverse_keyed
# and items$max under its own 1-based contract).
poly_items <- function(p = 12L, reverse = FALSE, ncat = 5L) {
  rk <- if (length(reverse) == 1L) rep(reverse, p) else reverse
  data.frame(reverse_keyed = rk, max = ncat)
}

# Recode 1..ncat -> 0..(ncat-1) for a hand-built PerFit call.
zero_base <- function(m) {
  z <- m - 1L
  storage.mode(z) <- "integer"
  z
}

# ---- Schema -----------------------------------------------------------------

test_that("cier_gnormed returns a list-based cier_index with the pinned schema", {
  skip_if_not_installed("PerFit")
  out <- cier_gnormed(poly_matrix(), poly_items(), seed = 1L)
  expect_s3_class(out, "cier_index")
  expect_type(out, "list")
  expect_identical(names(out),
                   c("value", "flagged", "method", "cutoff", "direction"))
  expect_type(out$value, "double")
  expect_type(out$flagged, "logical")
  expect_identical(length(out$value), 60L)
  expect_identical(length(out$flagged), 60L)
  expect_identical(out$method, "cier_gnormed")
  expect_identical(out$direction, "upper")
  expect_type(out$cutoff, "double")
})

test_that("as.data.frame.cier_index returns the tidy per-respondent frame", {
  skip_if_not_installed("PerFit")
  df <- as.data.frame(cier_gnormed(poly_matrix(n = 12L), poly_items(),
                                   cutoff = 0.5))
  expect_s3_class(df, "data.frame")
  expect_identical(names(df), c("value", "flagged"))
  expect_identical(nrow(df), 12L)
})

# ---- Independent oracle parity (round-to-4dp, tolerance 0) -------------------

test_that("cier_gnormed$value equals the oracle on a complete matrix", {
  skip_if_not_installed("PerFit")
  m <- poly_matrix(n = 60L, p = 12L, seed = 21L)
  items <- poly_items(12L)
  got <- cier_gnormed(m, items, cutoff = 0.5)$value
  ref <- round(ref_personfit_gnormed_poly(list(responses = m, items = items)), 4)
  expect_equal(got, ref, tolerance = 0)
})

test_that("cier_gnormed$value equals the oracle when rows carry NAs", {
  skip_if_not_installed("PerFit")
  m <- poly_matrix(n = 60L, p = 12L, seed = 24L)
  m[2L, 3L] <- NA                   # one missing cell -> respondent excluded
  m[10L, ] <- NA                    # all-NA respondent -> excluded
  items <- poly_items(12L)
  got <- cier_gnormed(m, items, cutoff = 0.5)$value
  ref <- round(ref_personfit_gnormed_poly(list(responses = m, items = items)), 4)
  expect_equal(got, ref, tolerance = 0)
  expect_true(is.na(got[[2L]]))
  expect_true(is.na(got[[10L]]))
})

test_that("cier_gnormed reverse-scores keyed items (oracle parity, not a no-op)", {
  skip_if_not_installed("PerFit")
  m <- poly_matrix(n = 60L, p = 12L, seed = 22L)
  rk <- rep(c(FALSE, TRUE), 6L)
  items <- poly_items(12L, reverse = rk)
  got <- cier_gnormed(m, items, cutoff = 0.5)$value
  ref <- round(ref_personfit_gnormed_poly(list(responses = m, items = items)), 4)
  expect_equal(got, ref, tolerance = 0)
  # Reverse-scoring is not a no-op: ignoring the keying gives a different score.
  raw <- ref_personfit_gnormed_poly(list(responses = m, items = poly_items(12L)))
  expect_false(isTRUE(all.equal(got, round(raw, 4))))
})

# ---- Cross-package parity: PerFit preprocessing + dichotomous reduction ------

test_that("cier_gnormed matches a hand-built PerFit::Gnormed.poly call (bytewise)", {
  # Pins the bridge preprocessing: zero-basing to 0..(Ncat-1) and persons-as-rows
  # orientation. n != p so a missing-transpose mutant returns the wrong length /
  # values; a raw-1..k coding mutant shifts every score. The scorer is shared, so
  # this is bytewise (tolerance 0).
  skip_if_not_installed("PerFit")
  m <- poly_matrix(n = 50L, p = 12L, seed = 7L)
  got <- cier_gnormed(m, poly_items(12L), cutoff = 0.5)$value
  fit <- PerFit::Gnormed.poly(matrix = zero_base(m), Ncat = 5L)
  ref <- as.numeric(fit$PFscores$PFscores)
  expect_equal(got, ref, tolerance = 0)
})

test_that("cier_gnormed reduces to the dichotomous PerFit::Gnormed at Ncat = 2", {
  skip_if_not_installed("PerFit")
  b <- poly_matrix(n = 80L, p = 12L, ncat = 2L, seed = 11L)
  got <- cier_gnormed(b, poly_items(12L, ncat = 2L), cutoff = 0.5)$value
  dich <- as.numeric(PerFit::Gnormed(zero_base(b))$PFscores$PFscores)
  ok <- is.finite(got) & is.finite(dich)
  expect_gt(sum(ok), 10L)
  expect_lt(max(abs(got[ok] - dich[ok])), 1e-9)
})

# ---- Cutoff: PerFit Monte-Carlo null -----------------------------------------

test_that("the default cutoff is the seeded PerFit Monte-Carlo null", {
  skip_if_not_installed("PerFit")
  m <- poly_matrix(n = 60L, p = 12L, seed = 5L)
  fit <- PerFit::Gnormed.poly(matrix = zero_base(m), Ncat = 5L)
  # Same seed immediately before PerFit::cutoff on both sides: Gnormed.poly is
  # deterministic, so the simulated null draws from an identical RNG stream.
  set.seed(123L)
  expected <- as.numeric(PerFit::cutoff(fit, Blvl = 0.05)$Cutoff)
  got <- cier_gnormed(m, poly_items(12L), seed = 123L)$cutoff
  expect_equal(got, expected, tolerance = 0)
})

test_that("the fpr argument moves the Monte-Carlo nominal level (Blvl)", {
  skip_if_not_installed("PerFit")
  m <- poly_matrix(n = 60L, p = 12L, seed = 5L)
  fit <- PerFit::Gnormed.poly(matrix = zero_base(m), Ncat = 5L)
  set.seed(7L)
  expected <- as.numeric(PerFit::cutoff(fit, Blvl = 0.10)$Cutoff)
  got <- cier_gnormed(m, poly_items(12L), fpr = 0.10, seed = 7L)$cutoff
  expect_equal(got, expected, tolerance = 0)
})

test_that("a seeded null cutoff is reproducible and leaves the ambient RNG intact", {
  skip_if_not_installed("PerFit")
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
  # The NULL-seed path must NOT save/restore: it consumes the caller's RNG, the
  # complement of the seeded path's restore contract. A mutant that always
  # save/restores (or ignores the RNG and returns a fixed cutoff) is caught by
  # the advanced-stream assertion.
  skip_if_not_installed("PerFit")
  m <- poly_matrix(n = 50L, p = 10L, seed = 3L)
  items <- poly_items(10L)
  set.seed(1L)
  before <- get(".Random.seed", envir = globalenv())
  a <- cier_gnormed(m, items)$cutoff
  after <- get(".Random.seed", envir = globalenv())
  expect_true(is.finite(a))
  expect_false(identical(before, after))   # ambient stream advanced (no restore)
})

# ---- Direction (upper) ------------------------------------------------------

test_that("direction is upper: a high-Guttman-error row flags, a low one does not", {
  # Deterministic comparator check with a literal mid-cutoff (independent of the
  # simulated null): the highest score flags, the lowest does not, and the flag
  # rule is value >= cutoff. A flag-lower mutant inverts both.
  skip_if_not_installed("PerFit")
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

test_that("fewer than three items abstains (PerFit denominator needs >= 3)", {
  skip_if_not_installed("PerFit")
  m <- poly_matrix(n = 40L, p = 2L, seed = 26L)
  expect_warning(out <- cier_gnormed(m, poly_items(2L)),
                 class = "cier_warning_insufficient_items")
  expect_true(all(is.na(out$value)))
  expect_true(is.na(out$cutoff))
  expect_true(all(is.na(out$flagged)))
  # Length pins the full-length abstain vectors -- `all(is.na(numeric(0)))` is
  # vacuously TRUE, so a zero-length / scalar result would slip past the asserts
  # above; the index must stay one row per respondent even when it abstains.
  expect_identical(length(out$value), 40L)
  expect_identical(length(out$flagged), 40L)
})

test_that("exactly three items scores (the items knife-edge, scoring side)", {
  # Brackets the >= 3 boundary from below: p = 2 abstains (above), p = 3 must
  # score, so an off-by-one (abstain at p <= 3) is caught. Parity (not just
  # finiteness) pins value correctness at the minimum item count, where a small-p
  # preprocessing/orientation bug could otherwise hide.
  skip_if_not_installed("PerFit")
  m <- poly_matrix(n = 40L, p = 3L, seed = 28L)
  items <- poly_items(3L)
  out <- cier_gnormed(m, items, cutoff = 0.5)
  expect_true(any(is.finite(out$value)))
  ref <- round(ref_personfit_gnormed_poly(list(responses = m, items = items)), 4)
  expect_equal(out$value, ref, tolerance = 0)
})

test_that("fewer than two complete respondents abstains", {
  skip_if_not_installed("PerFit")
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
  # 2 complete must score, so a `< 3 complete` off-by-one is caught. The two
  # complete rows are forced to span the 1..5 scale (PerFit needs both extreme
  # categories present); the count threshold -- not the span -- is what this pins.
  skip_if_not_installed("PerFit")
  m <- poly_matrix(n = 6L, p = 6L, seed = 29L)
  m[1L, 1L] <- 1                    # row 1 alone spans the 1..5 scale (min..max)
  m[1L, 2L] <- 5
  m[3:6, 1L] <- NA                  # only respondents 1 and 2 are complete
  items <- poly_items(6L)
  out <- cier_gnormed(m, items, cutoff = 0.5)
  expect_true(is.finite(out$value[[1L]]))
  expect_true(is.finite(out$value[[2L]]))
  expect_true(all(is.na(out$value[3:6])))
  # Parity pins the two scored values (not just their finiteness) and the NA
  # placement against the oracle at the 2-complete boundary.
  ref <- round(ref_personfit_gnormed_poly(list(responses = m, items = items)), 4)
  expect_equal(out$value, ref, tolerance = 0)
})

test_that("a single-respondent matrix abstains", {
  skip_if_not_installed("PerFit")
  m <- poly_matrix(n = 1L, p = 5L, seed = 1L)
  expect_warning(out <- cier_gnormed(m, poly_items(5L)),
                 class = "cier_warning_insufficient_items")
  expect_true(all(is.na(out$value)))
  expect_identical(length(out$value), 1L)
  expect_identical(length(out$flagged), 1L)
})

test_that("a straightliner gets a finite (low) score, not an abstention", {
  # Documented blind spot: unlike Ht, the Guttman-error index scores a
  # straightliner; straightlining evades detection rather than abstaining.
  skip_if_not_installed("PerFit")
  m <- poly_matrix(n = 40L, p = 10L, seed = 23L)
  m[1L, ] <- 3                      # midpoint straightline (complete row)
  out <- cier_gnormed(m, poly_items(10L), cutoff = 0.5)
  expect_true(is.finite(out$value[[1L]]))
  expect_true(all(is.finite(out$value)))
})

test_that("a respondent with one missing cell abstains; rows stay aligned", {
  skip_if_not_installed("PerFit")
  m <- poly_matrix(n = 40L, p = 10L, seed = 24L)
  m[2L, 3L] <- NA
  out <- cier_gnormed(m, poly_items(10L), cutoff = 0.5)
  expect_true(is.na(out$value[[2L]]))
  expect_true(is.na(out$flagged[[2L]]))
  expect_gt(sum(is.finite(out$value)), 2L)
})

test_that("a fully-constant item column matches the oracle (popularity-rank edge)", {
  # A zero-variance item ties its popularity rank; the popularity-rank numerator
  # and the knapsack denominator must still agree with the oracle. Distinct from
  # the straightliner ROW above -- this is a constant COLUMN.
  skip_if_not_installed("PerFit")
  m <- poly_matrix(n = 50L, p = 10L, seed = 30L)
  m[, 4L] <- 3                      # one constant (straightlined) item column
  items <- poly_items(10L)
  got <- cier_gnormed(m, items, cutoff = 0.5)$value
  ref <- round(ref_personfit_gnormed_poly(list(responses = m, items = items)), 4)
  expect_equal(got, ref, tolerance = 0)
})

test_that("the scale base (items$min) is honoured for keying and zero-basing", {
  # Min-invariance: the same responses in 1..ncat coding (min = 1, default) and
  # in 0..(ncat-1) coding (min = 0) must score identically -- with a reverse-keyed
  # item present so the (min + max) - x reflection is exercised, not just the
  # zero-base. A bridge that hardcodes `- 1` zero-basing or `(max + 1) - x`
  # reflection (ignoring min), or that derives Ncat as `max` instead of
  # `max - min + 1` (the 0-based declaration has max = 4 but five categories),
  # diverges or errors on the 0-based scale. The ported oracle hardcodes a base
  # of 1, so this is asserted as an internal invariance.
  skip_if_not_installed("PerFit")
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
  # PerFit's single-Ncat contract is about the NUMBER of response options
  # (max - min + 1), not about identical min/max pairs: a battery mixing 1..5
  # and 0..4 items (both five options) is valid and must equal the all-1..5
  # scoring of the same underlying responses. A homogeneity check on `max`
  # alone (or a global rather than per-item zero-base) would reject or
  # mis-score it.
  skip_if_not_installed("PerFit")
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
  # The smallest valid scale through the personfit validator: max >= min + 1 is
  # the bound, so a validator demanding max >= 2 regardless of min would
  # wrongly reject this valid declaration. Equality with the same data in 1..2
  # coding pins the full min/max plumbing at Ncat = 2.
  skip_if_not_installed("PerFit")
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
  # The bridge zero-bases to 0..(Ncat - 1) and range-checks; a code above
  # `max` (or below `min`) is a data/contract violation, not silently
  # coerced. check_responses only rejects NaN/Inf, so this pins the bridge's own
  # range guard. Crucially it is a PLAIN cier_error_input, NOT a backend limit:
  # an out-of-range value is a genuine data defect that must keep propagating
  # through cier_screen(), so the backend-limit subclass (reserved for
  # otherwise-valid data the backend cannot score) must NOT be attached here --
  # a mutant that tags both zero-base branches alike would wrongly let the screen
  # swallow corrupt data.
  skip_if_not_installed("PerFit")
  m <- poly_matrix(n = 10L, p = 6L, seed = 32L)
  m[1L, 1L] <- 6                    # exceeds max = 5
  err <- tryCatch(cier_gnormed(m, poly_items(6L)), error = function(e) e)
  expect_s3_class(err, "cier_error_input")
  expect_false(inherits(err, "cier_error_backend_limit"))
})

test_that("data that does not span the declared scale is a screen-survivable backend limit", {
  # PerFit's item-step popularities are undefined when a declared extreme category
  # never occurs (min == 0, max == Ncat - 1 after zero-basing). Responses confined
  # to 2..4 of a declared 1..5 scale never reach the extremes, so the bridge
  # surfaces a typed error instead of PerFit's terse abort. This is OTHERWISE-VALID
  # data the backend cannot score (sample-dependent, exactly like the
  # heterogeneous-span and mokken-ceiling cases), NOT a metadata defect -- so the
  # abort carries the cier_error_backend_limit subclass with a compact data$reason,
  # and cier_screen() records Gnormed as skipped-with-reason instead of aborting
  # the whole battery. It remains a cier_error_input for direct callers.
  skip_if_not_installed("PerFit")
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
    "sample does not attain both scale extremes (PerFit needs every end category observed)"
  )
})

test_that("a fractional (non-integer) response is a typed input error", {
  # check_responses() only rejects NaN/Inf; a fractional cell (e.g. an averaged
  # or imputed value) would otherwise be silently truncated by the integer cast
  # before PerFit. The bridge catches it as a typed error -- no silent coercion.
  skip_if_not_installed("PerFit")
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
  # fractional and non-finite max are PLAIN input errors -- although both also
  # make the spans unequal, per-item validity is classified before homogeneity,
  # so neither may carry the backend-limit subclass (the screen must propagate
  # them, not skip).
  for (bad_max in list(c(5, 5, 5, 5, 5, 2.5), c(5, 5, 5, 5, 5, Inf))) {
    err <- tryCatch(cier_gnormed(m, data.frame(max = bad_max)),
                    error = function(e) e)
    expect_s3_class(err, "cier_error_input")
    expect_false(inherits(err, "cier_error_backend_limit"))
  }
})

test_that("a heterogeneous span is a typed backend limit, not a plain input error", {
  # Accurate metadata describing genuinely mixed-format data (here five-option
  # and four-option items together) is not a malformed frame: it is PerFit's
  # single-Ncat contract that cannot score it. The abort therefore carries the
  # cier_error_backend_limit subclass -- the same line the mokken 10-category
  # ceiling draws -- so cier_screen() can skip-with-reason instead of dying.
  # It remains a cier_error_input for direct callers.
  m <- poly_matrix(n = 10L, p = 6L)
  het <- data.frame(max = c(5, 5, 5, 5, 5, 4))
  expect_error(cier_gnormed(m, het), class = "cier_error_input")
  expect_error(cier_gnormed(m, het), class = "cier_error_backend_limit")
  # The limit is on the span (number of options), NOT on max itself: equal max
  # with differing min is just as heterogeneous (spans 4 vs 5)...
  het_min <- data.frame(max = 5, min = c(1, 1, 1, 1, 1, 0))
  expect_error(cier_gnormed(m, het_min), class = "cier_error_backend_limit")
  # ...while a malformed per-item max (NA) stays a PLAIN input error: the
  # backend-limit subclass must not swallow genuine metadata defects.
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
  # A normed Guttman value is in [0, 1]; a threshold outside that flags everyone
  # or no one, so it is rejected.
  expect_error(cier_gnormed(m, items, cutoff = -0.1), class = "cier_error_input")
  expect_error(cier_gnormed(m, items, cutoff = 1.5), class = "cier_error_input")
  expect_error(cier_gnormed(m, items, seed = 1.5), class = "cier_error_input")
})

test_that("supplying both fpr and cutoff is a typed input error", {
  m <- poly_matrix(n = 10L, p = 6L)
  expect_error(cier_gnormed(m, poly_items(6L), fpr = 0.1, cutoff = 0.5),
               class = "cier_error_input")
})

test_that("cier_gnormed aborts with a typed error when PerFit is absent", {
  # The standalone index requires PerFit to score; the bridge gates through the
  # mockable cier_namespace_present() so absence surfaces a typed input error.
  testthat::local_mocked_bindings(cier_namespace_present = function(...) FALSE)
  expect_error(cier_gnormed(poly_matrix(n = 10L, p = 6L), poly_items(6L)),
               class = "cier_error_input")
})

# ---- print snapshot (locked, design-first; direction = upper) ---------------

test_that("print renders the locked cli summary (upper direction)", {
  skip_if_not_installed("PerFit")
  # A literal cutoff keeps the printed threshold stable across platforms (the
  # simulated null would vary the number); the snapshot pins the print FORMAT.
  m <- poly_matrix(n = 30L, seed = 11L)
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    expect_snapshot(print(cier_gnormed(m, poly_items(), cutoff = 0.5)))
  })
})

test_that("print reports abstaining respondents as '(no score)'", {
  skip_if_not_installed("PerFit")
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    m <- poly_matrix(n = 30L, p = 12L, seed = 11L)
    m[1L, 3L] <- NA                 # one missing-cell respondent -> abstains
    expect_snapshot(print(cier_gnormed(m, poly_items(12L), cutoff = 0.5)))
  })
})

test_that("the abstaining-row noun is '(no score)' and not '(no responses)'", {
  # Pinned independently of the generated snapshot: a wrapper that forgets to
  # extend abstention_noun() for cier_gnormed would bless the wrong noun on first
  # generation. Gnormed abstains on a fully-answered-but-one-missing-cell row, so
  # '(no score)' is the honest wording.
  skip_if_not_installed("PerFit")
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    m <- poly_matrix(n = 30L, p = 12L, seed = 11L)
    m[1L, 3L] <- NA
    out <- capture.output(print(cier_gnormed(m, poly_items(12L), cutoff = 0.5)))
  })
  expect_true(any(grepl("(no score)", out, fixed = TRUE)))
  expect_false(any(grepl("(no responses)", out, fixed = TRUE)))
})
