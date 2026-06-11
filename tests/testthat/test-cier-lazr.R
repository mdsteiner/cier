# Tests for cier_lazr() (Laz.R, Biemann et al. 2025).
#
# Trust model: the independent oracle (ref_lazr) re-derives the per-respondent
# Laz.R with an explicit transition double-loop and NEVER calls the production
# kernel. Laz.R has no CRAN parity partner (verified 2026-06-10), so the oracle
# plus the paper's published worked examples (John = 33/49; the one-liner
# Laz.R(c(1,2,3,4,5,4,3,2,1,2)) = 2/3) are the parity checks, like PR / RPR.
# Conventions under test (signed off this slice): drop-NA transitions, abstain
# below two valid transitions, matrix-only (anchor-count-invariant) scoring,
# integer-coded responses, percentile / upper cutoff.

source(test_path("..", "reference", "ref-lazr-biemann2025.R"))

# ---- Fixtures ---------------------------------------------------------------

lazr_fixture <- function(n = 20L, p = 30L, seed = 2026L) {
  withr::with_seed(seed, {
    m <- matrix(sample.int(5L, n * p, replace = TRUE), nrow = n, ncol = p)
  })
  storage.mode(m) <- "double"
  m
}

john_matrix <- function(extra_rows = 0L) {
  s <- ref_lazr_john_sequence()
  p <- length(s)
  m <- matrix(s, nrow = 1L, ncol = p)
  if (extra_rows > 0L) {
    m <- rbind(m, matrix(rep(s, extra_rows), nrow = extra_rows, byrow = TRUE))
  }
  storage.mode(m) <- "double"
  m
}

# ---- Schema -----------------------------------------------------------------

test_that("cier_lazr returns a list-based cier_index schema", {
  out <- cier_lazr(lazr_fixture(n = 12L))
  expect_s3_class(out, "cier_index")
  expect_type(out, "list")
  expect_identical(names(out),
                   c("value", "flagged", "method", "cutoff", "direction"))
  expect_type(out$value, "double")
  expect_type(out$flagged, "logical")
  expect_identical(length(out$value), 12L)
  expect_identical(length(out$flagged), 12L)
  expect_identical(out$method, "cier_lazr")
  expect_identical(out$direction, "upper")
  expect_type(out$cutoff, "double")
})

test_that("as.data.frame.cier_index returns the tidy per-respondent frame", {
  df <- as.data.frame(cier_lazr(lazr_fixture(n = 8L)))
  expect_s3_class(df, "data.frame")
  expect_identical(names(df), c("value", "flagged"))
  expect_identical(nrow(df), 8L)
})

# ---- Published worked examples (the spec) -----------------------------------

test_that("cier_lazr reproduces John's worked example (Eq. 4): 33/49", {
  out <- cier_lazr(john_matrix())
  expect_equal(out$value[[1L]], 33 / 49, tolerance = 1e-12)
})

test_that("footnote-2 one-liner: Laz.R(c(1,2,3,4,5,4,3,2,1,2)) = 2/3", {
  m <- matrix(c(1, 2, 3, 4, 5, 4, 3, 2, 1, 2), nrow = 1L)
  storage.mode(m) <- "double"
  expect_equal(cier_lazr(m)$value[[1L]], 2 / 3, tolerance = 1e-12)
})

test_that("the oracle reproduces the paper's John transition matrix", {
  # Guards the oracle itself (independent of production, which exposes no T):
  # T_12 = 9 and T_21 = T_23 = T_32 = T_34 = T_43 = 8 over anchors 1..4.
  tmat <- ref_lazr_row(ref_lazr_john_sequence())$transition_matrix
  expect_equal(tmat,
               matrix(c(0L, 9L, 0L, 0L,
                        8L, 0L, 8L, 0L,
                        0L, 8L, 0L, 8L,
                        0L, 0L, 8L, 0L),
                      nrow = 4L, byrow = TRUE))
})

# ---- Independent oracle parity (1e-12) --------------------------------------

test_that("cier_lazr$value equals the oracle on random integer data", {
  x <- lazr_fixture(n = 25L, p = 30L)
  expect_equal(cier_lazr(x)$value, ref_lazr(x), tolerance = 1e-12)
})

test_that("NA transitions drop from the denominator (drop-NA, not N-1)", {
  # 10-item row 1,2,NA,4,3,2,1,2,3,4 -> the (2,NA) and (NA,4) pairs drop, so 7
  # valid transitions remain. sum(P*T) = 5, so Laz.R = 5/7 (NOT 5/9 = 5/(N-1)).
  row <- c(1, 2, NA, 4, 3, 2, 1, 2, 3, 4)
  m <- matrix(row, nrow = 1L)
  storage.mode(m) <- "double"
  out <- cier_lazr(m)
  expect_equal(out$value[[1L]], 5 / 7, tolerance = 1e-12)        # NOT 5/9
  expect_equal(out$value[[1L]], ref_lazr_row(row)$value, tolerance = 1e-12)
})

# ---- Convention pins --------------------------------------------------------

test_that("a straightliner scores exactly 1 and is flagged, never NA", {
  out <- cier_lazr(matrix(rep(3, 10L), nrow = 1L))
  expect_identical(out$value, 1)
  expect_false(is.na(out$value))
  expect_true(out$flagged[[1L]])
})

test_that("a diagonal-liner (1..5 repeated) scores 1 (deterministic chain)", {
  m <- matrix(rep_len(1:5, 20L), nrow = 1L)
  storage.mode(m) <- "double"
  expect_equal(cier_lazr(m)$value[[1L]], 1, tolerance = 1e-12)
})

test_that("values are bounded in (0, 1]", {
  v <- cier_lazr(lazr_fixture(n = 40L, p = 30L))$value
  expect_true(all(v > 0 & v <= 1 + 1e-12))
})

test_that("value is invariant to a constant integer shift (anchor-count-free)", {
  # The s-invariance that justifies the matrix-only design: predictability
  # depends on the transition structure, not the absolute anchor labels. A
  # mutant that hard-codes anchors 1..max (ignoring the base) fails the 0-based
  # recoding here.
  x <- lazr_fixture(n = 15L, p = 20L)
  expect_equal(cier_lazr(x)$value, cier_lazr(x + 10L)$value, tolerance = 1e-12)
  expect_equal(cier_lazr(x)$value, cier_lazr(x - 1L)$value, tolerance = 1e-12)
})

test_that("a stray large integer does not inflate the bin space (overflow guard)", {
  # An un-recoded numeric missing code / sentinel is a whole number, so it passes
  # validation. The kernel ranks DISTINCT observed anchors, so a single 1e5 in
  # 1..5 data leaves the bin space bounded (s = 6 here) instead of sizing it by
  # the value span (s = 1e5), which would overflow `n * s * s` and crash. The
  # value matches the oracle (which also ranks distinct anchors), and rows
  # without the sentinel are unchanged by it.
  x <- lazr_fixture(n = 10L, p = 20L)
  with_sentinel <- x
  with_sentinel[1L, 5L] <- 1e5
  out <- expect_no_error(cier_lazr(with_sentinel))
  expect_equal(out$value, ref_lazr(with_sentinel), tolerance = 1e-12)
  expect_equal(out$value[-1L], cier_lazr(x)$value[-1L], tolerance = 1e-12)
})

# ---- Abstention / NA edges --------------------------------------------------

test_that("an all-NA row abstains and keeps rows aligned", {
  m <- rbind(
    c(1, 2, 3, 4, 3, 2, 1, 2, 3, 4),
    rep(NA_real_, 10L),                    # abstains, in the middle
    c(5, 4, 3, 2, 1, 2, 3, 4, 5, 4)
  )
  storage.mode(m) <- "double"
  out <- cier_lazr(m)
  expect_true(is.na(out$value[[2L]]))
  expect_true(is.na(out$flagged[[2L]]))
  expect_false(is.na(out$value[[1L]]))
  expect_false(is.na(out$value[[3L]]))
})

test_that("a single valid transition abstains (the < 2 rule)", {
  # Row 1 has exactly one valid transition (2,5); it must abstain rather than
  # report the degenerate 1.0. Row 2 scores so the matrix does not wholly
  # abstain (no percentile warning here).
  m <- rbind(
    c(2, 5, rep(NA_real_, 8L)),
    c(1, 2, 3, 4, 5, 4, 3, 2, 1, 2)
  )
  storage.mode(m) <- "double"
  out <- cier_lazr(m)
  expect_true(is.na(out$value[[1L]]))
  expect_true(is.na(out$flagged[[1L]]))
  expect_false(is.na(out$value[[2L]]))
})

test_that("exactly two valid transitions scores (the < 2 boundary, positive edge)", {
  # Two valid transitions (2,5),(5,1): above the threshold, so the row scores.
  # Pins the positive edge so a regression to an `n >= 1` rule is visible.
  m <- matrix(c(2, 5, 1, rep(NA_real_, 7L)), nrow = 1L)
  storage.mode(m) <- "double"
  out <- cier_lazr(m)
  expect_false(is.na(out$value[[1L]]))
  expect_true(is.finite(out$value[[1L]]))
})

test_that("a single-column matrix yields no transitions: all abstain", {
  expect_warning(
    out <- cier_lazr(matrix(c(1, 2, 3), ncol = 1L)),
    class = "cier_warning_insufficient_items"
  )
  expect_true(all(is.na(out$value)))
})

test_that("a wholly abstaining matrix warns and flags nobody", {
  expect_warning(
    out <- cier_lazr(matrix(NA_real_, nrow = 3L, ncol = 10L)),
    class = "cier_warning_insufficient_items"
  )
  expect_true(all(is.na(out$value)))
  expect_true(is.na(out$cutoff))
  expect_true(all(is.na(out$flagged)))
})

# ---- Input validation -------------------------------------------------------

test_that("a non-matrix / non-numeric payload is a typed input error", {
  expect_error(cier_lazr(1:10), class = "cier_error_input")
  expect_error(cier_lazr(matrix(letters[1:8], nrow = 2L)),
               class = "cier_error_input")
  expect_error(cier_lazr(NULL), class = "cier_error_input")
})

test_that("NaN / infinite values are typed input errors", {
  expect_error(cier_lazr(matrix(c(1, 2, 3, Inf, 1, 2), nrow = 2L)),
               class = "cier_error_input")
  expect_error(cier_lazr(matrix(c(1, 2, 3, NaN, 1, 2), nrow = 2L)),
               class = "cier_error_input")
})

test_that("non-integer responses are a typed input error", {
  expect_error(cier_lazr(matrix(c(1, 2, 2.5, 4, 3, 2), nrow = 2L)),
               class = "cier_error_input")
})

test_that("bad fpr values are typed input errors", {
  x <- lazr_fixture(n = 6L)
  expect_error(cier_lazr(x, fpr = 0), class = "cier_error_input")
  expect_error(cier_lazr(x, fpr = 1), class = "cier_error_input")
  expect_error(cier_lazr(x, fpr = -0.1), class = "cier_error_input")
  expect_error(cier_lazr(x, fpr = NA_real_), class = "cier_error_input")
  expect_error(cier_lazr(x, fpr = c(0.05, 0.1)), class = "cier_error_input")
  expect_error(cier_lazr(x, fpr = "x"), class = "cier_error_input")
})

test_that("invalid literal cutoff values are typed input errors", {
  x <- lazr_fixture(n = 6L)
  expect_error(cier_lazr(x, cutoff = 1.5), class = "cier_error_input")
  expect_error(cier_lazr(x, cutoff = -0.1), class = "cier_error_input")
  expect_error(cier_lazr(x, cutoff = NA_real_), class = "cier_error_input")
  expect_error(cier_lazr(x, cutoff = c(0.1, 0.2)), class = "cier_error_input")
  expect_error(cier_lazr(x, cutoff = "x"), class = "cier_error_input")
})

test_that("supplying both fpr and cutoff is a typed input error", {
  expect_error(cier_lazr(lazr_fixture(n = 6L), fpr = 0.1, cutoff = 0.5),
               class = "cier_error_input")
})

# ---- Cutoff: default, fpr override, NO double-flip, direction ---------------

test_that("default cutoff is the upper-tail 95th percentile (NO double-flip)", {
  x <- lazr_fixture(n = 60L, p = 20L)
  out <- cier_lazr(x)
  expect_equal(out$cutoff,
               as.numeric(stats::quantile(out$value, 0.95,
                                          names = FALSE, type = 7L)),
               tolerance = 1e-12)
})

test_that("the fpr argument moves the percentile target", {
  x <- lazr_fixture(n = 60L, p = 20L)
  out <- cier_lazr(x, fpr = 0.10)
  expect_equal(out$cutoff,
               as.numeric(stats::quantile(out$value, 0.90,
                                          names = FALSE, type = 7L)),
               tolerance = 1e-12)
})

test_that("direction is upper: high-predictability rows flag, low ones do not", {
  x <- rbind(rep(3, 20L), lazr_fixture(n = 20L, p = 20L))   # constant -> value 1
  storage.mode(x) <- "double"
  out <- cier_lazr(x)
  expect_true(out$flagged[[1L]])                            # constant, value 1
  expect_false(out$flagged[[which.min(out$value)]])         # least predictable
  expect_identical(out$flagged,
                   !is.na(out$value) & out$value >= out$cutoff)
})

test_that("a literal cutoff passes through and flags via the upper direction", {
  x <- rbind(rep(3, 20L), lazr_fixture(n = 10L, p = 20L))
  storage.mode(x) <- "double"
  out <- cier_lazr(x, cutoff = 0.95)
  expect_identical(out$cutoff, 0.95)
  expect_true(out$flagged[[1L]])                            # value 1 >= 0.95
  expect_identical(out$flagged, !is.na(out$value) & out$value >= 0.95)
})

# ---- Cutoff: the paper-faithful Kneedle elbow (kneedle = TRUE) --------------
# Biemann et al.'s own studies flag the top 5% (= the fpr default); the Kneedle
# elbow (Satopaa et al. 2011) is the sample-specific cutoff they offer in their
# companion app. cier ships it as an opt-in `kneedle = TRUE`, oracle-only trust.

source(test_path("..", "reference", "ref-kneedle-satopaa-2011.R"))

test_that("kneedle = TRUE sets the cutoff to the Satopaa elbow (oracle, tol 0)", {
  x <- lazr_fixture(n = 70L, p = 18L)
  out <- cier_lazr(x, kneedle = TRUE)
  # End-to-end paper-faithfulness: the cutoff IS the convex/increasing elbow of
  # the observed scores. The oracle drops the NA-abstaining rows, as the resolver
  # does. This also proves kneedle is not aliased to the percentile default.
  expect_identical(out$cutoff,
                   ref_kneedle(out$value, "convex", "increasing")$value)
  expect_identical(out$direction, "upper")
  expect_s3_class(out, "cier_index")
})

test_that("kneedle = TRUE flags the predictable high tail, not the bulk", {
  # Five straightliners (Laz.R = 1) on top of a careful bulk: the elbow sits at
  # the bulk/spike boundary, so every straightliner flags while the least
  # predictable respondent does not.
  # Five constant rows (each a straightliner -> Laz.R = 1) over a careful bulk.
  x <- rbind(matrix(c(2, 4, 1, 5, 3), nrow = 5L, ncol = 18L),
             lazr_fixture(n = 35L, p = 18L))
  storage.mode(x) <- "double"
  out <- cier_lazr(x, kneedle = TRUE)
  expect_true(all(out$flagged[1:5]))                       # the spike flags
  expect_false(out$flagged[[which.min(out$value)]])        # the bulk floor does not
  expect_false(all(out$flagged))                           # a tail, not everyone
})

test_that("kneedle = TRUE is deterministic (no RNG)", {
  x <- lazr_fixture(n = 50L, p = 16L)
  expect_identical(cier_lazr(x, kneedle = TRUE)$cutoff,
                   cier_lazr(x, kneedle = TRUE)$cutoff)
})

test_that("kneedle is mutually exclusive with fpr and with cutoff", {
  x <- lazr_fixture(n = 12L)
  expect_error(cier_lazr(x, kneedle = TRUE, fpr = 0.1),
               class = "cier_error_input")
  expect_error(cier_lazr(x, kneedle = TRUE, cutoff = 0.5),
               class = "cier_error_input")
})

test_that("a non-flag kneedle argument is a typed input error", {
  x <- lazr_fixture(n = 12L)
  expect_error(cier_lazr(x, kneedle = "yes"), class = "cier_error_input")
  expect_error(cier_lazr(x, kneedle = NA), class = "cier_error_input")
  expect_error(cier_lazr(x, kneedle = c(TRUE, FALSE)), class = "cier_error_input")
})

test_that("kneedle abstains (NA cutoff + warning, flags nobody) below three scores", {
  # Two scoring rows plus all-NA rows: only two finite Laz.R values remain, so
  # the elbow is undefined and the cutoff abstains -- the percentile path's
  # contract, reached through the kneedle resolver.
  x <- rbind(c(1, 2, 3, 4, 5, 4, 3, 2, 1, 2),
             c(5, 4, 3, 2, 1, 2, 3, 4, 5, 4),
             matrix(NA_real_, nrow = 3L, ncol = 10L))
  storage.mode(x) <- "double"
  expect_warning(out <- cier_lazr(x, kneedle = TRUE),
                 class = "cier_warning_insufficient_items")
  expect_true(is.na(out$cutoff))
  expect_false(any(out$flagged, na.rm = TRUE))     # an NA cutoff flags nobody
  expect_false(out$flagged[[1L]])                  # a scored row -> FALSE, not NA
  expect_true(all(is.na(out$flagged[3:5])))        # the all-NA rows abstain
})

# ---- print snapshot (locked; reuses the shared cier_index print) ------------

test_that("print renders the locked cli summary (upper direction)", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    expect_snapshot(print(cier_lazr(lazr_fixture(n = 30L, p = 12L))))
  })
})

test_that("print reports abstaining respondents on their own line", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    x <- lazr_fixture(n = 29L, p = 12L)
    x <- rbind(x, rep(NA_real_, 12L))     # one abstaining respondent
    expect_snapshot(print(cier_lazr(x)))
  })
})
