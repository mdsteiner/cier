# Tests for cier_ht() -- the polytomous person-scalability Ht person-fit index
# (Sijtsma & Meijer 1992; Molenaar 1991), scored via mokken::coefH on the
# transposed scale with an empirical-percentile cutoff (lower tail).
#
# Trust model: the Ht scorer IS mokken::coefH(t(z))$Hi (single-kernel rule), so
# the genuine independent check is the closed-form oracle (ref_personfit_ht_poly:
# the Frechet / rearrangement collapse, re-derived from scratch, never calling
# the production kernel). coefH returns full precision (unlike PerFit's 4-dp
# rounding), so the oracle holds to ~1e-14 (tolerance 1e-12; see
# tests/reference/TOLERANCES.md). The mokken direct-call parity exists to pin the
# BRIDGE PREPROCESSING -- the transpose, global zero-basing, reverse-keying, and
# complete-casing -- not the scorer (which is shared); an n != p fixture makes a
# missing-transpose mutant break. The dichotomous reduction to PerFit::Ht pins
# that the polytomous coefH path reduces to the classic Ht at Ncat = 2 (1e-4,
# PerFit's 4-dp rounding). The cutoff is the empirical lower-tail percentile (a
# ranking convention; no PerFit Monte-Carlo null exists for the mokken-backed
# polytomous Ht).

source(test_path("..", "reference", "ref-personfit-niessen-2016.R"))

# A reproducible polytomous response matrix in 1..ncat coding (the wrapper
# zero-bases internally). Uniform sampling gives a spread of Ht values (needed
# for the direction and percentile tests).
poly_matrix <- function(n = 60L, p = 12L, ncat = 5L, seed = 21L) {
  withr::with_seed(seed, {
    m <- matrix(sample.int(ncat, n * p, replace = TRUE), nrow = n, ncol = p)
  })
  storage.mode(m) <- "double"
  m
}

# Item metadata for cier_ht(): 1-based fixtures (1..ncat, so `max = ncat`),
# optional `reverse_keyed`. The same data.frame doubles as the oracle's
# `data$items` (the oracle reads items$reverse_keyed and items$max under its
# own 1-based contract).
poly_items <- function(p = 12L, reverse = FALSE, ncat = 5L) {
  rk <- if (length(reverse) == 1L) rep(reverse, p) else reverse
  data.frame(reverse_keyed = rk, max = ncat)
}

# Production-equivalent hand call of the mokken scorer on a complete,
# non-keyed block: global zero-base (z - min(z)) then coefH on the TRANSPOSE,
# with the same NaN -> NA reduction the kernel applies. Used by the bridge-
# preprocessing parity check.
mokken_ht_ref <- function(m) {
  z <- m - min(m)
  storage.mode(z) <- "integer"
  invisible(utils::capture.output(
    res <- suppressWarnings(mokken::coefH(t(z), se = FALSE))
  ))
  hi <- res$Hi
  v <- as.numeric(if (is.null(dim(hi))) hi else hi[, 1L])
  v[!is.finite(v)] <- NA_real_
  v
}

# ---- Schema -----------------------------------------------------------------

test_that("cier_ht returns a list-based cier_index with the pinned schema", {
  skip_if_not_installed("mokken")
  out <- cier_ht(poly_matrix(), poly_items())
  expect_s3_class(out, "cier_index")
  expect_type(out, "list")
  expect_identical(names(out),
                   c("value", "flagged", "method", "cutoff", "direction"))
  expect_type(out$value, "double")
  expect_type(out$flagged, "logical")
  expect_identical(length(out$value), 60L)
  expect_identical(length(out$flagged), 60L)
  expect_identical(out$method, "cier_ht")
  expect_identical(out$direction, "lower")
  expect_type(out$cutoff, "double")
})

test_that("as.data.frame.cier_index returns the tidy per-respondent frame", {
  skip_if_not_installed("mokken")
  df <- as.data.frame(cier_ht(poly_matrix(n = 12L), poly_items(), cutoff = 0))
  expect_s3_class(df, "data.frame")
  expect_identical(names(df), c("value", "flagged"))
  expect_identical(nrow(df), 12L)
})

# ---- Independent oracle parity (1e-12) --------------------------------------

test_that("cier_ht$value equals the independent oracle on a complete matrix", {
  skip_if_not_installed("mokken")
  m <- poly_matrix(n = 60L, p = 12L, seed = 21L)
  items <- poly_items(12L)
  got <- cier_ht(m, items, cutoff = 0)$value
  ref <- ref_personfit_ht_poly(list(responses = m, items = items))
  expect_equal(got, ref, tolerance = 1e-12)
})

test_that("cier_ht$value equals the oracle when rows carry NAs", {
  skip_if_not_installed("mokken")
  m <- poly_matrix(n = 60L, p = 12L, seed = 24L)
  m[2L, 3L] <- NA                   # one missing cell -> respondent excluded
  m[10L, ] <- NA                    # all-NA respondent -> excluded
  items <- poly_items(12L)
  got <- cier_ht(m, items, cutoff = 0)$value
  ref <- ref_personfit_ht_poly(list(responses = m, items = items))
  expect_equal(got, ref, tolerance = 1e-12)
  expect_true(is.na(got[[2L]]))
  expect_true(is.na(got[[10L]]))
})

test_that("cier_ht reverse-scores keyed items (oracle parity, not a no-op)", {
  skip_if_not_installed("mokken")
  m <- poly_matrix(n = 60L, p = 12L, seed = 22L)
  rk <- rep(c(FALSE, TRUE), 6L)
  items <- poly_items(12L, reverse = rk)
  got <- cier_ht(m, items, cutoff = 0)$value
  ref <- ref_personfit_ht_poly(list(responses = m, items = items))
  expect_equal(got, ref, tolerance = 1e-12)
  # Reverse-scoring is not a no-op: ignoring the keying gives a different score
  # (a misaligned scale collapses person scalability).
  raw <- ref_personfit_ht_poly(list(responses = m, items = poly_items(12L)))
  expect_false(isTRUE(all.equal(got, raw)))
})

# ---- Cross-package parity: mokken preprocessing + dichotomous reduction ------

test_that("cier_ht matches a hand-built mokken::coefH(t(.)) call (bytewise)", {
  # Pins the bridge preprocessing's persons-as-rows orientation (the TRANSPOSE):
  # n != p so a missing-transpose mutant (coefH(z), item scalability) returns the
  # wrong length and errors / mismatches. The scorer is shared, so this is
  # bytewise (tolerance 0). (Global zero-basing is translation-invariant for the
  # covariance-based Ht, so it is the oracle parity above -- not this check --
  # that exercises the rest of the preprocessing.)
  skip_if_not_installed("mokken")
  m <- poly_matrix(n = 50L, p = 12L, seed = 7L)
  got <- cier_ht(m, poly_items(12L), cutoff = 0)$value
  ref <- mokken_ht_ref(m)
  expect_equal(got, ref, tolerance = 0)
})

test_that("cier_ht reduces to the dichotomous PerFit::Ht at Ncat = 2", {
  skip_if_not_installed("mokken")
  skip_if_not_installed("PerFit")
  b <- poly_matrix(n = 80L, p = 12L, ncat = 2L, seed = 11L)
  got <- cier_ht(b, poly_items(12L, ncat = 2L), cutoff = 0)$value
  invisible(utils::capture.output(pf <- PerFit::Ht(b - 1L)))
  perf <- as.numeric(pf$PFscores$PFscores)
  ok <- is.finite(got) & is.finite(perf)
  expect_gt(sum(ok), 10L)
  expect_lt(max(abs(got[ok] - perf[ok])), 1e-4)
})

test_that("the scale base (items$min) is honoured for keying and zero-basing", {
  # Min-invariance: the same responses in 1..ncat coding (min = 1, default) and
  # in 0..(ncat-1) coding (min = 0) must score identically -- with a reverse-keyed
  # item present so the (min + max) - x reflection is exercised, not just the
  # zero-base. A wrapper that hardcodes `(max + 1) - x` reflection
  # (ignoring min) diverges on the 0-based scale. The ported oracle hardcodes a
  # base of 1, so this is asserted as an internal invariance.
  skip_if_not_installed("mokken")
  rk <- rep(c(FALSE, TRUE), 5L)
  m1 <- poly_matrix(n = 50L, p = 10L, ncat = 5L, seed = 31L)   # 1..5 coding
  items1 <- poly_items(10L, reverse = rk)                       # min defaults 1
  m0 <- m1 - 1L                                                 # 0..4 coding
  items0 <- data.frame(reverse_keyed = rk, max = 4L, min = 0L)
  expect_equal(cier_ht(m1, items1, cutoff = 0)$value,
               cier_ht(m0, items0, cutoff = 0)$value,
               tolerance = 0)
})

# ---- Cutoff: empirical lower-tail percentile --------------------------------

test_that("the default cutoff is the lower-tail empirical percentile (fpr 0.05)", {
  # The single direction flip: a LOWER index takes the `fpr` quantile (not the
  # 1 - fpr quantile). A mutant flipping to the upper tail would differ.
  skip_if_not_installed("mokken")
  m <- poly_matrix(n = 60L, p = 12L, seed = 21L)
  out <- cier_ht(m, poly_items(12L))
  fin <- out$value[is.finite(out$value)]
  expect_equal(out$cutoff,
               as.numeric(stats::quantile(fin, 0.05, names = FALSE, type = 7L)),
               tolerance = 0)
})

test_that("the fpr argument moves the percentile tail mass", {
  skip_if_not_installed("mokken")
  m <- poly_matrix(n = 60L, p = 12L, seed = 21L)
  out <- cier_ht(m, poly_items(12L), fpr = 0.10)
  fin <- out$value[is.finite(out$value)]
  expect_equal(out$cutoff,
               as.numeric(stats::quantile(fin, 0.10, names = FALSE, type = 7L)),
               tolerance = 0)
})

# ---- Direction (lower) ------------------------------------------------------

test_that("direction is lower: a low-Ht row flags, a high one does not", {
  # Deterministic comparator check with a literal mid-cutoff (independent of the
  # percentile): the lowest score flags, the highest does not, and the flag rule
  # is value <= cutoff. A flag-upper mutant inverts both.
  skip_if_not_installed("mokken")
  m <- poly_matrix(n = 60L, p = 12L, seed = 11L)
  items <- poly_items(12L)
  v <- cier_ht(m, items, cutoff = 0)$value
  mid <- stats::median(v, na.rm = TRUE)
  out <- cier_ht(m, items, cutoff = mid)
  expect_true(out$flagged[[which.min(out$value)]])
  expect_false(out$flagged[[which.max(out$value)]])
  expect_identical(out$flagged, !is.na(out$value) & out$value <= out$cutoff)
})

# ---- Edge cases -------------------------------------------------------------

test_that("a straightliner abstains (NA), unlike Gnormed which scores it", {
  # Headline difference from Gnormed: a complete zero-variance row is
  # structurally unscorable, so coefH returns NaN -> NA. The OTHER respondents
  # still score, and the kept rows match the oracle (which abstains on the same
  # row). A mutant that scores the straightliner (or drops it from alignment)
  # breaks both the NA placement and the oracle parity.
  skip_if_not_installed("mokken")
  m <- poly_matrix(n = 40L, p = 10L, seed = 40L)
  m[1L, ] <- 3                      # midpoint straightline (complete row)
  items <- poly_items(10L)
  out <- cier_ht(m, items, cutoff = 0)
  expect_true(is.na(out$value[[1L]]))
  expect_identical(sum(is.finite(out$value)), 39L)
  ref <- ref_personfit_ht_poly(list(responses = m, items = items))
  expect_equal(out$value, ref, tolerance = 1e-12)
})

test_that("a fully-constant item column matches the oracle (variance edge)", {
  # A zero-variance item COLUMN (distinct from the straightliner ROW above): the
  # covariance numerator and the rearrangement denominator must still agree with
  # the oracle.
  skip_if_not_installed("mokken")
  m <- poly_matrix(n = 50L, p = 10L, seed = 30L)
  m[, 4L] <- 3                      # one constant (straightlined) item column
  items <- poly_items(10L)
  got <- cier_ht(m, items, cutoff = 0)$value
  ref <- ref_personfit_ht_poly(list(responses = m, items = items))
  expect_equal(got, ref, tolerance = 1e-12)
})

test_that("an all-constant block abstains (fewer than two non-constant rows)", {
  # Every respondent straightlines -> coefH would error / return all NaN; the
  # kernel's < 2 non-constant guard returns all NA, and the percentile cutoff
  # abstains with the shared insufficient-items warning.
  skip_if_not_installed("mokken")
  m <- matrix(((0:7) %% 5L) + 1L, nrow = 8L, ncol = 5L)  # every row constant
  storage.mode(m) <- "double"
  expect_warning(out <- cier_ht(m, poly_items(5L)),
                 class = "cier_warning_insufficient_items")
  expect_true(all(is.na(out$value)))
  expect_true(is.na(out$cutoff))
  expect_true(all(is.na(out$flagged)))
  expect_identical(length(out$value), 8L)
  expect_identical(length(out$flagged), 8L)
})

test_that("fewer than two complete respondents abstains", {
  skip_if_not_installed("mokken")
  m <- poly_matrix(n = 5L, p = 6L, seed = 27L)
  m[1:4, 1L] <- NA                  # only respondent 5 is complete
  expect_warning(out <- cier_ht(m, poly_items(6L)),
                 class = "cier_warning_insufficient_items")
  expect_true(all(is.na(out$value)))
  expect_identical(length(out$value), 5L)
  expect_identical(length(out$flagged), 5L)
})

test_that("a single-respondent matrix abstains", {
  skip_if_not_installed("mokken")
  m <- poly_matrix(n = 1L, p = 5L, seed = 1L)
  expect_warning(out <- cier_ht(m, poly_items(5L)),
                 class = "cier_warning_insufficient_items")
  expect_true(all(is.na(out$value)))
  expect_identical(length(out$value), 1L)
  expect_identical(length(out$flagged), 1L)
})

test_that("a single-item (p = 1) matrix abstains cleanly, not an untyped crash", {
  # Person scalability is undefined on one item; with a single column the per-row
  # variance is NA, which (without an ncol guard) would turn the non-constant-row
  # guard into `if (NA)` -- an untyped base error. The kernel must abstain (all
  # NA) gracefully, the same way cier_gnormed does on too-few items.
  skip_if_not_installed("mokken")
  m <- matrix(c(1, 2, 3, 4, 5, 2, 3, 1, 4, 5), ncol = 1L)
  expect_warning(out <- cier_ht(m, poly_items(1L)),
                 class = "cier_warning_insufficient_items")
  expect_true(all(is.na(out$value)))
  expect_true(is.na(out$cutoff))
  expect_identical(length(out$value), 10L)
  expect_identical(length(out$flagged), 10L)
})

test_that("a respondent with one missing cell abstains; rows stay aligned", {
  skip_if_not_installed("mokken")
  m <- poly_matrix(n = 40L, p = 10L, seed = 24L)
  m[2L, 3L] <- NA
  out <- cier_ht(m, poly_items(10L), cutoff = 0)
  expect_true(is.na(out$value[[2L]]))
  expect_true(is.na(out$flagged[[2L]]))
  expect_gt(sum(is.finite(out$value)), 2L)
})

test_that("a fractional (non-integer) response is a typed input error", {
  # check_responses() only rejects NaN/Inf; a fractional cell (e.g. an averaged
  # or imputed value) would otherwise be silently truncated by the integer cast
  # before coefH. The kernel catches it as a typed error -- no silent coercion.
  skip_if_not_installed("mokken")
  m <- poly_matrix(n = 30L, p = 6L, seed = 34L)
  m[1L, 1L] <- 2.5
  expect_error(cier_ht(m, poly_items(6L)), class = "cier_error_input")
})

test_that("a fractional response errors even when the matrix would abstain", {
  # The whole-number contract is checked BEFORE the abstention guards, so a
  # fractional cell surfaces a typed error rather than being silently swallowed
  # by an all-constant / too-small block that abstains. (Otherwise the user is
  # never told their data is non-integer.)
  skip_if_not_installed("mokken")
  m <- matrix(2, nrow = 8L, ncol = 5L)   # all-constant -> would abstain (all NA)
  m[1L, 1L] <- 2.5                        # but a fractional cell must error first
  expect_error(cier_ht(m, poly_items(5L)), class = "cier_error_input")
})

# ---- Input validation -------------------------------------------------------

test_that("a non-matrix / non-numeric / non-finite payload is a typed input error", {
  expect_error(cier_ht(1:10, poly_items(10L)), class = "cier_error_input")
  expect_error(cier_ht(matrix(letters[1:6], nrow = 2L), poly_items(3L)),
               class = "cier_error_input")
  bad <- poly_matrix(n = 6L, p = 6L)
  bad[1L, 1L] <- Inf
  expect_error(cier_ht(bad, poly_items(6L)), class = "cier_error_input")
})

test_that("Ht requires `max` only on reverse-keyed items", {
  # mokken::coefH accepts a mix of category counts and the kernel never reads
  # `max`, so -- unlike Gnormed's single-Ncat contract -- Ht needs it only
  # to reverse-score keyed items. A reverse-keyed item with an absent / NA /
  # at-or-below-min `max` is a typed error (it cannot be reverse-scored).
  m <- poly_matrix(n = 10L, p = 6L)
  rev_last <- c(rep(FALSE, 5L), TRUE)
  no_max <- data.frame(reverse_keyed = rev_last)
  na_max <- data.frame(reverse_keyed = rev_last, max = c(5, 5, 5, 5, 5, NA))
  low_max <- data.frame(reverse_keyed = rev_last, max = c(5, 5, 5, 5, 5, 1))
  frac_max <- data.frame(reverse_keyed = rev_last, max = c(5, 5, 5, 5, 5, 2.5))
  inf_max <- data.frame(reverse_keyed = rev_last, max = c(5, 5, 5, 5, 5, Inf))
  expect_error(cier_ht(m, no_max), class = "cier_error_input")
  expect_error(cier_ht(m, na_max), class = "cier_error_input")
  expect_error(cier_ht(m, low_max), class = "cier_error_input")
  expect_error(cier_ht(m, frac_max), class = "cier_error_input")
  expect_error(cier_ht(m, inf_max), class = "cier_error_input")
})

test_that("Ht allows NA max / min on forward items in a mixed battery", {
  # max / min are read only for the reflection of keyed items, so NA on a
  # FORWARD item must not abort a battery that does contain a reverse item
  # (a require-on-every-item validator would be over-strict).
  skip_if_not_installed("mokken")
  m <- poly_matrix(n = 30L, p = 6L)
  it <- data.frame(reverse_keyed = c(rep(FALSE, 5L), TRUE),
                   max = c(NA, 5, 5, 5, 5, 5), min = c(NA, 1, 1, 1, 1, 1))
  expect_s3_class(cier_ht(m, it, cutoff = 0), "cier_index")
})

test_that("a reverse-keyed two-option item (min = 0, max = 1) is accepted", {
  # max >= min + 1 is the bound: a 0/1 reverse item reflects by (0 + 1) - x.
  # A validator demanding max >= 2 regardless of min would reject it.
  skip_if_not_installed("mokken")
  m <- poly_matrix(n = 30L, p = 6L)              # items 1-5: 1..5 coding
  m[, 6L] <- m[, 6L] %% 2                        # item 6: 0/1 coding
  it <- data.frame(reverse_keyed = c(rep(FALSE, 5L), TRUE),
                   max = c(rep(5L, 5L), 1L), min = c(rep(1L, 5L), 0L))
  expect_s3_class(cier_ht(m, it, cutoff = 0), "cier_index")
})

test_that("Ht accepts forward-only items without (or with heterogeneous) max", {
  # No reverse keys -> `max` is never used for scoring, so an absent or
  # heterogeneous `max` column must be accepted (mokken handles mixed
  # scales); the homogeneity Gnormed needs is not Ht's contract.
  skip_if_not_installed("mokken")
  m <- poly_matrix(n = 30L, p = 6L)
  expect_s3_class(cier_ht(m, data.frame(reverse_keyed = rep(FALSE, 6L)),
                          cutoff = 0), "cier_index")
  expect_s3_class(cier_ht(m, data.frame(max = c(5, 5, 5, 5, 5, 4)),
                          cutoff = 0), "cier_index")
})

test_that("a scale wider than mokken's 10-category ceiling is a typed backend-limit error", {
  # mokken::coefH's check.data hard-stops with a raw simpleError when the global
  # zero-based range exceeds 9 ("mokken cannot ... handle [more than] 10
  # categories"). The kernel converts that backend ceiling into a typed
  # cier_error_input BEFORE coefH runs, with the cier_error_backend_limit
  # subclass cier_screen() catches to skip-with-reason. An 11-point scale
  # (spanning 1..11) trips it; so does a heterogeneous battery where a single
  # item spans more than 10 points (the ceiling is the GLOBAL range, not a
  # per-item property).
  skip_if_not_installed("mokken")
  m11 <- poly_matrix(n = 30L, p = 4L, ncat = 11L)
  m11[1L, 1L] <- 1   # force both extremes so the global range is exactly 0..10
  m11[2L, 1L] <- 11
  items11 <- data.frame(reverse_keyed = rep(FALSE, 4L))
  expect_error(cier_ht(m11, items11), class = "cier_error_backend_limit")
  expect_error(cier_ht(m11, items11), class = "cier_error_input")
  # Heterogeneous battery: 5-point items plus one wide item spanning 1..11.
  m_mix <- poly_matrix(n = 30L, p = 5L, ncat = 5L)
  wide <- poly_matrix(n = 30L, p = 1L, ncat = 11L, seed = 9L)
  wide[1L, 1L] <- 1
  wide[2L, 1L] <- 11
  m_mix <- cbind(m_mix, wide)
  expect_error(cier_ht(m_mix, data.frame(reverse_keyed = rep(FALSE, 6L))),
               class = "cier_error_backend_limit")
})

test_that("a 10-point scale (the ceiling boundary) still scores", {
  # 1..10 zero-bases to 0..9, exactly at mokken's limit: the guard must use a
  # strict > 9 comparison, not >= 9 (an off-by-one mutant would reject valid
  # 10-point data).
  skip_if_not_installed("mokken")
  m10 <- poly_matrix(n = 30L, p = 6L, ncat = 10L)
  m10[1L, 1L] <- 1
  m10[2L, 1L] <- 10
  out <- cier_ht(m10, data.frame(reverse_keyed = rep(FALSE, 6L)))
  expect_s3_class(out, "cier_index")
  expect_true(any(is.finite(out$value)))
})

test_that("a wrong number of item rows is a typed input error", {
  m <- poly_matrix(n = 10L, p = 6L)
  expect_error(cier_ht(m, data.frame(max = rep(5, 3L))),
               class = "cier_error_input")
})

test_that("a non-data.frame items or a non-integer min (reverse item) is a typed input error", {
  m <- poly_matrix(n = 10L, p = 6L)
  expect_error(cier_ht(m, "not a frame"), class = "cier_error_input")
  # `min` (the scale base) is validated on reverse-keyed items, where it sets the
  # reflection base; a fractional min on a reverse item is a typed error.
  bad_min <- data.frame(reverse_keyed = c(rep(FALSE, 5L), TRUE),
                        max = 5L, min = c(1, 1, 1, 1, 1, 1.5))
  expect_error(cier_ht(m, bad_min), class = "cier_error_input")
})

test_that("invalid fpr / cutoff values are typed input errors", {
  m <- poly_matrix(n = 10L, p = 6L)
  items <- poly_items(6L)
  expect_error(cier_ht(m, items, fpr = 0), class = "cier_error_input")
  expect_error(cier_ht(m, items, fpr = 1), class = "cier_error_input")
  expect_error(cier_ht(m, items, fpr = c(0.05, 0.1)), class = "cier_error_input")
  # An Ht value is in [-1, 1]; a threshold outside that flags everyone or no one.
  expect_error(cier_ht(m, items, cutoff = -1.5), class = "cier_error_input")
  expect_error(cier_ht(m, items, cutoff = 1.5), class = "cier_error_input")
})

test_that("supplying both fpr and cutoff is a typed input error", {
  m <- poly_matrix(n = 10L, p = 6L)
  expect_error(cier_ht(m, poly_items(6L), fpr = 0.1, cutoff = 0),
               class = "cier_error_input")
})

test_that("cier_ht aborts with a typed error when mokken is absent", {
  # The standalone index requires mokken to score; the wrapper gates through the
  # mockable cier_namespace_present() so absence surfaces a typed input error.
  testthat::local_mocked_bindings(cier_namespace_present = function(...) FALSE)
  expect_error(cier_ht(poly_matrix(n = 10L, p = 6L), poly_items(6L)),
               class = "cier_error_input")
})

# ---- print snapshot (locked, design-first; direction = lower) ---------------

test_that("print renders the locked cli summary (lower direction)", {
  skip_if_not_installed("mokken")
  m <- poly_matrix(n = 30L, seed = 11L)
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    expect_snapshot(print(cier_ht(m, poly_items())))
  })
})

test_that("print reports abstaining respondents as '(no score)'", {
  skip_if_not_installed("mokken")
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    m <- poly_matrix(n = 30L, p = 12L, seed = 11L)
    m[1L, ] <- 3                    # a straightliner -> abstains on a full row
    expect_snapshot(print(cier_ht(m, poly_items(12L))))
  })
})

test_that("the abstaining-row noun is '(no score)' and not '(no responses)'", {
  # Pinned independently of the generated snapshot: a wrapper that forgets to
  # extend abstention_noun() for cier_ht would bless the wrong noun on first
  # generation. Ht abstains on a fully-answered straightliner, so '(no score)'
  # is the honest wording.
  skip_if_not_installed("mokken")
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    m <- poly_matrix(n = 30L, p = 12L, seed = 11L)
    m[1L, ] <- 3
    out <- capture.output(print(cier_ht(m, poly_items(12L))))
  })
  expect_true(any(grepl("(no score)", out, fixed = TRUE)))
  expect_false(any(grepl("(no responses)", out, fixed = TRUE)))
})
