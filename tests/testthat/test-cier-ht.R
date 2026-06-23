# Tests for cier_ht() -- polytomous person-scalability Ht person-fit (Sijtsma &
# Meijer 1992; Molenaar 1991), computed in-package by the closed-form (Frechet /
# rearrangement) kernel with an empirical-percentile cutoff (lower tail).
#
# Trust model (INVERTED -- cier is now the implementation, not a bridge): the scorer
# is the in-package closed form, so the primary check is the INDEPENDENT closed-form
# oracle (ref_personfit_ht_poly: the same Frechet / rearrangement collapse derived
# afresh, never calling the kernel) at 1e-12. The mokken direct-call parity is the
# CROSS-PACKAGE check that our kernel reproduces mokken::coefH(t(z))$Hi (the original
# backend, now the oracle) to machine precision (~1e-12, no longer bytewise -- the
# kernel and mokken are independent code paths). The dichotomous reduction to
# PerFit::Ht pins that the polytomous form collapses to the classic Ht at Ncat = 2
# (1e-4, PerFit's 4-dp rounding). The DEFAULT cutoff is the empirical lower-tail
# percentile (a ranking convention); an opt-in model-conforming Monte-Carlo null
# (method = "mc_null") is selectable, reusing the shared WI-2 null engine with the Ht
# kernel, the lower tail, and perfect-vector exclusion -- its mechanism (seeded
# reproducibility, tail/perfect handling, a positive attentive cutoff) is pinned
# directly, plus a check that the default percentile path is unchanged. Because Ht no
# longer needs a backend, the oracle-parity and edge tests run with or without mokken;
# only the cross-package and dichotomous checks skip when absent.

source(test_path("..", "reference", "ref-personfit-niessen-2016.R"))

# poly_matrix() / poly_items() fixtures are shared with test-cier-gnormed.R via
# helper-personfit.R.

# Hand call of the mokken scorer on a complete, non-keyed block: global zero-base
# (z - min(z)) then coefH on the TRANSPOSE, with the same NaN -> NA reduction the
# kernel applies. The cross-package parity oracle -- the kernel must reproduce it to
# machine precision (mokken is no longer the runtime scorer).
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

test_that("cier_ht returns the pinned cier_index schema", {
  expect_cier_index_schema(cier_ht(poly_matrix(), poly_items()),
                           "cier_ht", "lower", 60L)
})

# ---- Independent oracle parity (1e-12) --------------------------------------
# The PRIMARY check now that cier owns the implementation: the kernel value must
# equal the independent closed-form oracle. mokken is not needed for these.

test_that("cier_ht$value equals the independent oracle on a complete matrix", {
  m <- poly_matrix(n = 60L, p = 12L, seed = 21L)
  items <- poly_items(12L)
  got <- cier_ht(m, items, cutoff = 0)$value
  ref <- ref_personfit_ht_poly(list(responses = m, items = items))
  expect_equal(got, ref, tolerance = 1e-12)
})

test_that("cier_ht$value equals the oracle when rows carry NAs", {
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

test_that("cier_ht matches mokken::coefH(t(.)) to machine precision", {
  # Cross-package parity (the inverted trust model): our in-package closed form must
  # reproduce mokken -- the original backend, now the oracle -- to machine precision.
  # n != p, so a mutant that dropped the transpose-equivalent collapse (item rather
  # than person scalability) would diverge here. Tolerance is 1e-12, NOT bytewise:
  # the kernel and mokken are independent code paths agreeing only to ~1e-15.
  skip_if_not_installed("mokken")
  m <- poly_matrix(n = 50L, p = 12L, seed = 7L)
  got <- cier_ht(m, poly_items(12L), cutoff = 0)$value
  ref <- mokken_ht_ref(m)
  expect_equal(got, ref, tolerance = 1e-12)
})

test_that("cier_ht reduces to the dichotomous PerFit::Ht at Ncat = 2", {
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
  # Min-invariance: the same responses in 1..ncat coding (min = 1, default) and in
  # 0..(ncat-1) coding (min = 0) must score identically, with a reverse-keyed item
  # present so the (min + max) - x reflection is exercised, not just the zero-base.
  # A wrapper hardcoding `(max + 1) - x` reflection (ignoring min) diverges on the
  # 0-based scale. The oracle hardcodes base 1, so this is an internal invariance.
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
  # 1 - fpr quantile). A mutant flipping to the upper tail differs.
  m <- poly_matrix(n = 60L, p = 12L, seed = 21L)
  out <- cier_ht(m, poly_items(12L))
  fin <- out$value[is.finite(out$value)]
  expect_equal(out$cutoff,
               as.numeric(stats::quantile(fin, 0.05, names = FALSE, type = 7L)),
               tolerance = 0)
})

test_that("the fpr argument moves the percentile tail mass", {
  m <- poly_matrix(n = 60L, p = 12L, seed = 21L)
  out <- cier_ht(m, poly_items(12L), fpr = 0.10)
  fin <- out$value[is.finite(out$value)]
  expect_equal(out$cutoff,
               as.numeric(stats::quantile(fin, 0.10, names = FALSE, type = 7L)),
               tolerance = 0)
})

# ---- Cutoff: opt-in model-conforming Monte-Carlo null (method = "mc_null") ---
# WI-3: the polytomous Ht null reuses the shared WI-2 engine (personfit_null_cutoff)
# with the Ht kernel ht_scores, tail = "lower", perfect = "excluded". Percentile stays
# the default; the null is opt-in and reproducible via `seed`.

test_that("method = 'mc_null' selects the Monte-Carlo null cutoff (provenance)", {
  m <- poly_matrix(n = 60L, p = 12L, seed = 21L)
  out <- cier_ht(m, poly_items(12L), method = "mc_null", seed = 1L)
  expect_s3_class(out, "cier_index")
  expect_identical(out$cutoff_method, "mc_null")   # routed through the null, not pctile
  expect_identical(out$cutoff_rate, 0.05)          # default nominal level recorded
  expect_identical(out$direction, "lower")
  expect_true(is.finite(out$cutoff) && out$cutoff >= -1 && out$cutoff <= 1)
})

test_that("a seeded mc_null cutoff is reproducible and leaves the ambient RNG intact", {
  # local_preserve_seed: this test pins the global stream to set.seed(999L); restore
  # the caller's ambient state on exit so the deterministic seed does not leak to the
  # next test (the seeded cier_ht calls save/restore around themselves, so the global
  # seed would otherwise exit pinned at the 999 state).
  withr::local_preserve_seed()
  m <- poly_matrix(n = 50L, p = 10L, seed = 3L)
  items <- poly_items(10L)
  set.seed(999L)
  before <- get(".Random.seed", envir = globalenv())
  c1 <- cier_ht(m, items, method = "mc_null", seed = 42L)$cutoff
  after <- get(".Random.seed", envir = globalenv())
  c2 <- cier_ht(m, items, method = "mc_null", seed = 42L)$cutoff
  expect_identical(c1, c2)             # same seed -> identical cutoff
  expect_identical(after, before)      # seeded call restored the caller's RNG
})

test_that("a NULL-seed mc_null cutoff draws from (and advances) the ambient stream", {
  # The complement of the seeded restore contract: a NULL seed consumes the caller's
  # RNG. A mutant that always save/restores is caught by the stream advance; a mutant
  # returning a fixed cutoff is caught by two consecutive draws differing.
  # local_preserve_seed stops this test's deliberate stream advance from leaking to the
  # next test.
  withr::local_preserve_seed()
  m <- poly_matrix(n = 50L, p = 10L, seed = 3L)
  items <- poly_items(10L)
  set.seed(1L)
  before <- get(".Random.seed", envir = globalenv())
  a <- cier_ht(m, items, method = "mc_null")$cutoff
  b <- cier_ht(m, items, method = "mc_null")$cutoff
  after <- get(".Random.seed", envir = globalenv())
  expect_true(is.finite(a) && is.finite(b))
  expect_false(identical(before, after))   # ambient stream advanced (no restore)
  expect_false(identical(a, b))            # consecutive NULL-seed draws differ (not fixed)
})

test_that("a larger fpr raises the lower-tail mc_null cutoff (same seed)", {
  # Ht is LOWER-tail, so the null cutoff is the `fpr` (not 1 - fpr) quantile. With the
  # null generation and bootstrap seed-fixed and independent of fpr, a larger fpr can
  # only RAISE the cutoff. An upper-tail mutant (blvl.use = 1 - fpr) would invert this.
  m <- poly_matrix(n = 60L, p = 12L, seed = 5L)
  items <- poly_items(12L)
  o05 <- cier_ht(m, items, method = "mc_null", fpr = 0.05, seed = 7L)
  o10 <- cier_ht(m, items, method = "mc_null", fpr = 0.10, seed = 7L)
  expect_gt(o10$cutoff, o05$cutoff)
  # The user fpr is recorded as the null's nominal rate (not hardcoded to the 0.05
  # default): a mutant fixing provenance$rate = 0.05 passes every default-fpr test but
  # dies here.
  expect_identical(o05$cutoff_rate, 0.05)
  expect_identical(o10$cutoff_rate, 0.10)
})

test_that("the Ht null engine is reproducible, lower-tailed, and perfect-excluded", {
  # Direct mechanism pin on the shared engine driven by the Ht kernel: same seed ->
  # identical cutoff; the lower-tail cutoff sits below the upper-tail one on the same
  # null; and perfect = "excluded" runs on polytomous data without over-pruning to an
  # empty matrix (the constant-row drop is matrixStats::rowVars > 0, not PerFit's
  # dichotomous heuristic). The cutoff is a finite Ht value in [-1, 1].
  m <- poly_matrix(n = 80L, p = 12L, seed = 5L)
  z <- m - min(m)
  storage.mode(z) <- "integer"
  ncat <- max(z) + 1L
  score <- function(null) ht_scores(null)
  run <- function(tl) {
    personfit_null_cutoff(z, ncat, score, tail = tl, perfect = "excluded",
                          fpr = 0.05)
  }
  c1 <- with_local_seed(9L, function() run("lower"))
  c2 <- with_local_seed(9L, function() run("lower"))
  expect_identical(c1, c2)
  expect_true(is.finite(c1) && c1 >= -1 && c1 <= 1)
  up <- with_local_seed(9L, function() run("upper"))
  expect_lt(c1, up)
})

test_that("the null engine returns NA (not a crash) when exclusion empties the null", {
  # perfect = "excluded" drops every constant null row. On a straightliner-dominated
  # block the sum-score-conditional resample can only produce constant rows, so the
  # null prunes below two non-constant vectors. The engine must return NA, NOT call
  # ht_scores on a 0-row matrix -- max(numeric(0)) = -Inf makes personfit_sorted_rows
  # evaluate 0:-Inf and abort with an untyped 'result would be a too long vector'.
  # Deterministic: an all-constant z yields an all-constant null for any RNG draw.
  z <- rbind(matrix(0L, 50L, 2L), matrix(4L, 50L, 2L))   # every row constant
  cut <- with_local_seed(1L, function() {
    personfit_null_cutoff(z, ncat = 5L, score = function(null) ht_scores(null),
                          tail = "lower", perfect = "excluded", fpr = 0.05)
  })
  expect_true(is.na(cut))                  # clean NA, no error
})

test_that("resolve_ht_cutoff warns (not silently NA) when the null degenerates", {
  # A SCORED block (res$z non-NULL) whose Monte-Carlo null degenerates must surface the
  # typed abstention, not a silent NA cutoff that flags no one. A hand-built all-constant
  # block (kernel_ht would never pass one, but resolve_ht_cutoff trusts its caller)
  # deterministically prunes the null to empty -> NA -> a degenerate-null warning.
  z <- rbind(matrix(0L, 50L, 2L), matrix(4L, 50L, 2L))
  res <- list(value = rep(NA_real_, 100L), z = z, ncat = 5L, abstain = NULL)
  expect_warning(
    cut <- resolve_ht_cutoff(res, fpr = 0.05, seed = 1L),
    class = "cier_warning_insufficient_items"
  )
  expect_true(is.na(cut))
})

test_that("the mc_null cutoff is a sensible positive Ht on attentive data", {
  # The model-conforming null's point: on a genuinely scalable (attentive) sample the
  # null Ht values are positive, so the lower-tail cutoff is positive -- the attentive
  # cutoff the percentile ranking convention cannot give. Attentive data comes from the
  # package's own GRM generator (cier_simulate, prevalence = 0). Person scalability
  # needs an item-difficulty ORDERING, so the items carry difficulty-shifted thresholds
  # (a common 4-cut shape slid across the 20 items); without that spread every item has
  # the same difficulty and Ht is ~0 for everyone (attentive included).
  p <- 20L
  base <- c(-1.5, -0.5, 0.5, 1.5)
  thresholds <- lapply(seq(-1.2, 1.2, length.out = p), function(o) base + o)
  it <- data.frame(scale = rep("E", p), max = 5)
  sim <- cier_simulate(n = 400L, it, prevalence = 0, thresholds = thresholds,
                       seed = 4L)
  items <- data.frame(reverse_keyed = rep(FALSE, p), max = 5L)
  # Sanity: the attentive sample really is scalable (a clearly positive observed Ht),
  # so a positive null cutoff is a property of the null, not of degenerate data.
  expect_gt(stats::median(cier_ht(sim$responses, items)$value, na.rm = TRUE), 0.2)
  out <- cier_ht(sim$responses, items, method = "mc_null", seed = 7L)
  expect_identical(out$cutoff_method, "mc_null")
  expect_gt(out$cutoff, 0)              # attentive null Ht cutoff is positive
  expect_lt(out$cutoff, 1)
})

test_that("method defaults to percentile -- the default path is unchanged", {
  # The headline no-silent-behaviour-change guarantee: leaving `method` unset is
  # byte-identical to method = "percentile", resolves the exact empirical lower-tail
  # percentile, and -- crucially -- builds no Monte-Carlo null, so it consumes NO RNG.
  m <- poly_matrix(n = 60L, p = 12L, seed = 21L)
  items <- poly_items(12L)
  default <- cier_ht(m, items)
  explicit <- cier_ht(m, items, method = "percentile")
  expect_identical(default$cutoff, explicit$cutoff)
  expect_identical(default$cutoff_method, "percentile")
  fin <- default$value[is.finite(default$value)]
  expect_equal(default$cutoff,
               as.numeric(stats::quantile(fin, 0.05, names = FALSE, type = 7L)),
               tolerance = 0)
  set.seed(123L)
  before <- get(".Random.seed", envir = globalenv())
  invisible(cier_ht(m, items))                       # default (percentile) path
  after <- get(".Random.seed", envir = globalenv())
  expect_identical(before, after)                    # the default path draws no RNG
})

test_that("method = 'mc_null' abstains (NA cutoff + warning) on an unscalable block", {
  # When the block cannot be scored (here every row is constant -> < 2 varying rows)
  # there is no null to build: the cutoff is NA, nobody is flagged, and the typed
  # insufficient-items warning names the cause -- matching the percentile path's
  # abstention class. The mc_null provenance is still recorded.
  m <- matrix(((0:7) %% 5L) + 1L, nrow = 8L, ncol = 5L)  # every row constant
  storage.mode(m) <- "double"
  expect_warning(out <- cier_ht(m, poly_items(5L), method = "mc_null"),
                 class = "cier_warning_insufficient_items")
  expect_true(is.na(out$cutoff))
  expect_true(all(is.na(out$flagged)))
  expect_identical(out$cutoff_method, "mc_null")
})

test_that("a literal cutoff overrides method = 'mc_null' (provenance is literal)", {
  # A literal threshold wins over the selected null (and `seed` is then unused): the
  # cutoff passes through verbatim, provenance is "literal", and the lower-tail flag
  # rule applies. A mutant letting method win would build a (stochastic) null instead.
  m <- poly_matrix(n = 60L, p = 12L, seed = 21L)
  items <- poly_items(12L)
  out <- cier_ht(m, items, method = "mc_null", cutoff = 0, seed = 1L)
  expect_identical(out$cutoff, 0)
  expect_identical(out$cutoff_method, "literal")
  expect_identical(out$flagged, !is.na(out$value) & out$value <= 0)
})

test_that("an invalid method, or a fractional mc_null seed, is a typed input error", {
  m <- poly_matrix(n = 10L, p = 6L)
  items <- poly_items(6L)
  expect_error(cier_ht(m, items, method = "bogus"), class = "cier_error_input")
  expect_error(cier_ht(m, items, method = c("percentile", "mc_null")),
               class = "cier_error_input")
  expect_error(cier_ht(m, items, method = "mc_null", seed = 1.5),
               class = "cier_error_input")
})

# ---- Direction (lower) ------------------------------------------------------

test_that("direction is lower: a low-Ht row flags, a high one does not", {
  # Deterministic comparator with a literal mid-cutoff (independent of the
  # percentile): the lowest score flags, the highest does not, flag rule is
  # value <= cutoff. A flag-upper mutant inverts both.
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
  # Headline difference from Gnormed: a complete zero-variance row is structurally
  # unscorable (its rearrangement denominator collapses to 0), so the kernel returns
  # NA. The OTHER respondents still score, and the kept rows match the oracle (which
  # abstains on the same row). A mutant that scores the straightliner (or drops it
  # from alignment) breaks both the NA placement and the oracle parity.
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
  # covariance numerator and the rearrangement denominator must still agree with the
  # oracle.
  m <- poly_matrix(n = 50L, p = 10L, seed = 30L)
  m[, 4L] <- 3                      # one constant (straightlined) item column
  items <- poly_items(10L)
  got <- cier_ht(m, items, cutoff = 0)$value
  ref <- ref_personfit_ht_poly(list(responses = m, items = items))
  expect_equal(got, ref, tolerance = 1e-12)
})

test_that("tied responses match the oracle (the row sort handles ties)", {
  # Heavy within-row ties (a 3-point scale, plus a deliberate two-value row and a
  # near-straightliner) exercise the per-respondent sort and the comonotonic
  # rearrangement denominator. The kernel and the oracle sort identically, so they
  # agree exactly even where ranks tie.
  m <- poly_matrix(n = 50L, p = 12L, ncat = 3L, seed = 71L)
  m[1L, ] <- c(rep(1L, 6L), rep(3L, 6L))   # a two-value tie pattern
  m[2L, ] <- c(rep(2L, 11L), 1L)           # a near-straightliner with one break
  items <- poly_items(12L, ncat = 3L)
  got <- cier_ht(m, items, cutoff = 0)$value
  ref <- ref_personfit_ht_poly(list(responses = m, items = items))
  expect_equal(got, ref, tolerance = 1e-12)
})

test_that("an all-constant block abstains (fewer than two non-constant rows)", {
  # Every respondent straightlines; the kernel's < 2 non-constant guard returns all
  # NA, and the percentile cutoff abstains with the shared insufficient-items warning.
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
  m <- poly_matrix(n = 5L, p = 6L, seed = 27L)
  m[1:4, 1L] <- NA                  # only respondent 5 is complete
  expect_warning(out <- cier_ht(m, poly_items(6L)),
                 class = "cier_warning_insufficient_items")
  expect_true(all(is.na(out$value)))
  expect_identical(length(out$value), 5L)
  expect_identical(length(out$flagged), 5L)
})

test_that("a single-respondent matrix abstains", {
  m <- poly_matrix(n = 1L, p = 5L, seed = 1L)
  expect_warning(out <- cier_ht(m, poly_items(5L)),
                 class = "cier_warning_insufficient_items")
  expect_true(all(is.na(out$value)))
  expect_identical(length(out$value), 1L)
  expect_identical(length(out$flagged), 1L)
})

test_that("a single-item (p = 1) matrix abstains cleanly, not an untyped crash", {
  # Person scalability is undefined on one item; with a single column the per-row
  # variance is NA, which (without an ncol guard) turns the non-constant-row guard
  # into `if (NA)` -- an untyped base error. The kernel must abstain (all NA)
  # gracefully, like cier_gnormed on too-few items.
  m <- matrix(c(1, 2, 3, 4, 5, 2, 3, 1, 4, 5), ncol = 1L)
  expect_warning(out <- cier_ht(m, poly_items(1L)),
                 class = "cier_warning_insufficient_items")
  expect_true(all(is.na(out$value)))
  expect_true(is.na(out$cutoff))
  expect_identical(length(out$value), 10L)
  expect_identical(length(out$flagged), 10L)
})

test_that("a respondent with one missing cell abstains; rows stay aligned", {
  m <- poly_matrix(n = 40L, p = 10L, seed = 24L)
  m[2L, 3L] <- NA
  out <- cier_ht(m, poly_items(10L), cutoff = 0)
  expect_true(is.na(out$value[[2L]]))
  expect_true(is.na(out$flagged[[2L]]))
  expect_gt(sum(is.finite(out$value)), 2L)
})

test_that("a fractional (non-integer) response is a typed input error", {
  # check_responses() only rejects NaN/Inf; a fractional cell (e.g. averaged or
  # imputed) is not a valid category code. assert_integer_responses() in the kernel
  # catches it as a typed error -- no silent coercion into a wrong score.
  m <- poly_matrix(n = 30L, p = 6L, seed = 34L)
  m[1L, 1L] <- 2.5
  expect_error(cier_ht(m, poly_items(6L)), class = "cier_error_input")
})

test_that("a fractional response errors even when the matrix would abstain", {
  # The whole-number contract is checked BEFORE the abstention guards, so a
  # fractional cell surfaces a typed error rather than being silently swallowed by
  # an all-constant / too-small block that abstains (otherwise the user is never
  # told their data is non-integer).
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
  # The kernel accepts a mix of category counts and never reads `max` for scoring,
  # so -- unlike Gnormed's single-Ncat contract -- Ht needs it only to reverse-score
  # keyed items. A reverse-keyed item with an absent / NA / at-or-below-min `max` is
  # a typed error (it cannot be reverse-scored).
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
  # max / min are read only for the reflection of keyed items, so NA on a FORWARD
  # item must not abort a battery that does contain a reverse item (a
  # require-on-every-item validator would be over-strict).
  m <- poly_matrix(n = 30L, p = 6L)
  it <- data.frame(reverse_keyed = c(rep(FALSE, 5L), TRUE),
                   max = c(NA, 5, 5, 5, 5, 5), min = c(NA, 1, 1, 1, 1, 1))
  expect_s3_class(cier_ht(m, it, cutoff = 0), "cier_index")
})

test_that("a reverse-keyed two-option item (min = 0, max = 1) is accepted", {
  # max >= min + 1 is the bound: a 0/1 reverse item reflects by (0 + 1) - x.
  # A validator demanding max >= 2 regardless of min would reject it.
  m <- poly_matrix(n = 30L, p = 6L)              # items 1-5: 1..5 coding
  m[, 6L] <- m[, 6L] %% 2                        # item 6: 0/1 coding
  it <- data.frame(reverse_keyed = c(rep(FALSE, 5L), TRUE),
                   max = c(rep(5L, 5L), 1L), min = c(rep(1L, 5L), 0L))
  expect_s3_class(cier_ht(m, it, cutoff = 0), "cier_index")
})

test_that("Ht accepts forward-only items without (or with heterogeneous) max", {
  # No reverse keys -> `max` is never used for scoring, so an absent or heterogeneous
  # `max` column must be accepted (the kernel handles mixed scales); the homogeneity
  # Gnormed needs is not Ht's contract. The forward out-of-range cross-check is
  # therefore scoped OFF for Ht's forward items -- declaring max = 4 on a column
  # carrying 5s is harmless metadata here (forward_range = FALSE), unlike the
  # split-half / Gnormed path.
  m <- poly_matrix(n = 30L, p = 6L)
  expect_s3_class(cier_ht(m, data.frame(reverse_keyed = rep(FALSE, 6L)),
                          cutoff = 0), "cier_index")
  expect_s3_class(cier_ht(m, data.frame(max = c(5, 5, 5, 5, 5, 4)),
                          cutoff = 0), "cier_index")
})

test_that("Ht still range-checks REVERSE items against the declared max", {
  # forward_range = FALSE drops only the FORWARD check; a reverse-keyed item whose
  # observed responses exceed its declared max would reflect off-scale, so it must
  # still abort (the reverse range is a contract for every index).
  m <- poly_matrix(n = 30L, p = 6L, seed = 21L)   # 1..5
  m[1L, 1L] <- 6L                                  # item 1 out of declared max 5
  it <- poly_items(6L, reverse = c(TRUE, rep(FALSE, 5L)))   # item 1 reverse, max 5
  expect_error(cier_ht(m, it, cutoff = 0), class = "cier_error_input")
})

test_that("wide scales past the former 10-category ceiling now score (oracle parity)", {
  # The mokken backend hard-stopped above 10 categories; the in-package closed form
  # has no such ceiling. An 11-point scale (1..11) and a heterogeneous battery whose
  # global range exceeds ten points -- both previously a typed backend-limit error --
  # now score, and the values match the independent oracle.
  m11 <- poly_matrix(n = 30L, p = 4L, ncat = 11L, seed = 51L)
  m11[1L, 1L] <- 1    # force both extremes so the global range is the full 1..11
  m11[2L, 1L] <- 11
  items11 <- data.frame(reverse_keyed = rep(FALSE, 4L))
  out11 <- cier_ht(m11, items11, cutoff = 0)
  expect_s3_class(out11, "cier_index")
  expect_true(any(is.finite(out11$value)))
  ref11 <- ref_personfit_ht_poly(list(responses = m11, items = items11))
  expect_equal(out11$value, ref11, tolerance = 1e-12)
  # Heterogeneous battery: 5-point items plus one wide item spanning 1..11 (the
  # global range, 0..10 zero-based, is what tripped the old ceiling).
  m_mix <- poly_matrix(n = 30L, p = 5L, ncat = 5L, seed = 52L)
  wide <- poly_matrix(n = 30L, p = 1L, ncat = 11L, seed = 9L)
  wide[1L, 1L] <- 1
  wide[2L, 1L] <- 11
  m_mix <- cbind(m_mix, wide)
  items_mix <- data.frame(reverse_keyed = rep(FALSE, 6L))
  got_mix <- cier_ht(m_mix, items_mix, cutoff = 0)$value
  ref_mix <- ref_personfit_ht_poly(list(responses = m_mix, items = items_mix))
  expect_equal(got_mix, ref_mix, tolerance = 1e-12)
  # A genuinely wide 0..100 scale (101 categories, 1..101 coding) -- the case the
  # docs explicitly promise -- scores and matches the oracle (mokken cannot run it).
  m101 <- poly_matrix(n = 40L, p = 6L, ncat = 101L, seed = 54L)
  m101[1L, 1L] <- 1
  m101[2L, 1L] <- 101
  items101 <- data.frame(reverse_keyed = rep(FALSE, 6L))
  got101 <- cier_ht(m101, items101, cutoff = 0)$value
  expect_true(any(is.finite(got101)))
  ref101 <- ref_personfit_ht_poly(list(responses = m101, items = items101))
  expect_equal(got101, ref101, tolerance = 1e-12)
})

# A per-row ascending sort as an unnamed double matrix -- the independent reference
# the counting sort must reproduce (NOT the fallback's own body, so the check is a
# genuine cross-implementation comparison, not a tautology).
base_sorted <- function(z) {
  s <- t(apply(z, 1L, sort))
  dimnames(s) <- NULL
  storage.mode(s) <- "double"
  s
}

test_that("the counting-sort branch equals a per-row sort across its regime", {
  # kernel_ht sorts via personfit_sorted_rows()'s difference-array counting sort for
  # max(z) <= 100. Validate that novel algorithm against base sort across small AND
  # many-pass scales (including the 11..100 middle the kernel tests otherwise skip),
  # on both double and integer storage (check_responses() preserves integer input).
  for (mx in c(5L, 18L, 50L, 100L)) {
    z <- poly_matrix(n = 50L, p = 9L, ncat = mx + 1L, seed = 70L + mx) - 1
    z[1L, 1L] <- 0                            # pin the full 0..mx span (both extremes)
    z[2L, 1L] <- mx
    expect_lte(max(z), 100)                   # counting-sort branch
    expect_identical(personfit_sorted_rows(z), base_sorted(z))
    zi <- z                                   # integer input must give the same result
    storage.mode(zi) <- "integer"
    expect_identical(personfit_sorted_rows(zi), base_sorted(z))
  }
})

test_that("the wide-scale fallback returns the right shape and matches the oracle", {
  # max(z) > 100 takes the comparison-sort fallback, which IS base R's sort, so a
  # value comparison to base_sorted would be a tautology (same code path). Pin only
  # what the fallback adds beyond base sort -- an unnamed double n x p matrix -- and
  # validate the covariance VALUES end-to-end through the kernel vs the independent
  # oracle (the oracle runs its own covariance formula, so that check is real).
  z <- poly_matrix(n = 40L, p = 8L, ncat = 130L, seed = 91L) - 1
  z[1L, 1L] <- 0
  z[2L, 1L] <- 129
  expect_gt(max(z), 100)                      # comparison-sort fallback
  out <- personfit_sorted_rows(z)
  expect_true(is.double(out))
  expect_null(dimnames(out))
  expect_identical(dim(out), c(40L, 8L))
  # Pin the global span on the matrix actually handed to cier_ht, so the kernel
  # deterministically zero-bases to 0..129 (m = 129 > 100) and takes the fallback --
  # not left to the random draw.
  m130 <- poly_matrix(n = 40L, p = 8L, ncat = 130L, seed = 92L)
  m130[1L, 1L] <- 1
  m130[2L, 1L] <- 130
  items130 <- data.frame(reverse_keyed = rep(FALSE, 8L))
  got <- cier_ht(m130, items130, cutoff = 0)$value
  ref <- ref_personfit_ht_poly(list(responses = m130, items = items130))
  expect_equal(got, ref, tolerance = 1e-12)
})

test_that("cier_ht scores an integer responses matrix (preserved storage path)", {
  # check_responses() keeps integer storage, so an integer matrix flows through the
  # kernel as integer (the path the covariance's `^2` and matrix-vector product must
  # handle); it must score identically to the oracle. Every other fixture is double
  # (poly_matrix casts), so this is the only end-to-end exercise of the integer path.
  mi <- poly_matrix(n = 40L, p = 10L, seed = 26L)
  storage.mode(mi) <- "integer"
  items <- poly_items(10L)
  out <- cier_ht(mi, items, cutoff = 0)
  expect_s3_class(out, "cier_index")
  expect_true(any(is.finite(out$value)))
  ref <- ref_personfit_ht_poly(list(responses = mi, items = items))
  expect_equal(out$value, ref, tolerance = 1e-12)
})

test_that("a 10-point scale scores (no off-by-one at the former boundary)", {
  # 1..10 was the old limit; it still scores and matches the oracle, so the boundary
  # carries no regression now that the ceiling is gone.
  m10 <- poly_matrix(n = 30L, p = 6L, ncat = 10L, seed = 53L)
  m10[1L, 1L] <- 1
  m10[2L, 1L] <- 10
  items10 <- data.frame(reverse_keyed = rep(FALSE, 6L))
  out <- cier_ht(m10, items10, cutoff = 0)
  expect_s3_class(out, "cier_index")
  expect_true(any(is.finite(out$value)))
  ref10 <- ref_personfit_ht_poly(list(responses = m10, items = items10))
  expect_equal(out$value, ref10, tolerance = 1e-12)
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

test_that("cier_ht scores without mokken (no runtime backend gate)", {
  # Targeted regression guard for the removed runtime gate: Ht is pure R, so it
  # scores normally without consulting mokken.
  m <- poly_matrix(n = 30L, p = 8L, seed = 61L)
  items <- poly_items(8L)
  out <- cier_ht(m, items, cutoff = 0)
  expect_s3_class(out, "cier_index")
  expect_true(any(is.finite(out$value)))
  ref <- ref_personfit_ht_poly(list(responses = m, items = items))
  expect_equal(out$value, ref, tolerance = 1e-12)
})

# ---- print snapshot (locked, design-first; direction = lower) ---------------

test_that("print renders the locked cli summary (lower direction)", {
  m <- poly_matrix(n = 30L, seed = 11L)
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    expect_snapshot(print(cier_ht(m, poly_items())))
  })
})

test_that("print reports abstaining respondents as '(no score)'", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    m <- poly_matrix(n = 30L, p = 12L, seed = 11L)
    m[1L, ] <- 3                    # a straightliner -> abstains on a full row
    expect_snapshot(print(cier_ht(m, poly_items(12L))))
  })
})

# ---- missing reverse_keyed inform -------------------------------------------

test_that("cier_ht informs when items omits reverse_keyed", {
  m <- poly_matrix(n = 60L, p = 12L)     # default seed 21 -- known-good shape
  it <- data.frame(max = rep(5L, 12L))   # no reverse_keyed column
  expect_message(suppressWarnings(cier_ht(m, it)),
                 class = "cier_message_forward_keyed")
})
