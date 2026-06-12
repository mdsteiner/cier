# Tests for cier_irv().
#
# Trust model: the independent oracle (ref_irv) re-derives the per-row sample SD
# with a hand-rolled two-pass sum-of-squares and never calls the production
# kernel; the cross-package check pins parity with careless::irv at the 1e-10
# tolerance recorded in tests/reference/TOLERANCES.md (matrixStats::rowSds vs
# stats::sd differ only at ulp level); the property / mutant-killer block targets
# each mutant named in dev/restart/index-specs.md card 2 (population vs sample
# SD, missing na.rm, flagging upper instead of lower).

source(test_path("..", "reference", "ref-irv-marjanovic-2015.R"))

# Analytic fixture of present-data rows with hand-checkable SDs:
#   row 1: 1..5            -> sd = sqrt(2.5)
#   row 2: constant 3      -> sd = 0
#   row 3: 1,5,1,5,1       -> sd = sqrt(4.8)
#   row 4: 2,2,2,2,NA      -> sd of c(2,2,2,2) = 0 (NA dropped, 4 present)
present_fixture <- function() {
  matrix(
    c(
      1, 2, 3, 4, 5,
      3, 3, 3, 3, 3,
      1, 5, 1, 5, 1,
      2, 2, 2, 2, NA
    ),
    nrow = 4L, byrow = TRUE
  )
}

# ---- Schema -----------------------------------------------------------------

test_that("cier_irv returns a list-based cier_index with the pinned schema", {
  # WP3: small/saturated fixtures trip the percentile-cutoff degeneracy guard
  # (D2/D7); these value/oracle tests assert the score, not the flag, so the
  # (correct) warning is muffled.
  out <- suppressWarnings(cier_irv(present_fixture()))
  expect_s3_class(out, "cier_index")
  expect_type(out, "list")
  expect_identical(names(out),
                   c("value", "flagged", "method", "cutoff", "direction"))
  expect_type(out$value, "double")
  expect_type(out$flagged, "logical")
  expect_identical(length(out$value), 4L)
  expect_identical(length(out$flagged), 4L)
  expect_identical(out$method, "cier_irv")
  expect_identical(out$direction, "lower")
  expect_type(out$cutoff, "double")
})

test_that("as.data.frame.cier_index returns the tidy per-respondent frame", {
  df <- as.data.frame(suppressWarnings(cier_irv(present_fixture())))
  expect_s3_class(df, "data.frame")
  expect_identical(names(df), c("value", "flagged"))
  expect_identical(nrow(df), 4L)
  expect_equal(df$value,
               c(sqrt(2.5), 0, sqrt(4.8), 0), tolerance = 1e-10)
})

# ---- Independent oracle parity (1e-10) --------------------------------------

test_that("cier_irv$value equals the hand-rolled oracle on the fixture", {
  x <- present_fixture()
  expect_equal(suppressWarnings(cier_irv(x))$value, ref_irv(x), tolerance = 1e-10)
})

test_that("cier_irv$value equals the oracle on a random complete matrix", {
  withr::with_seed(2026, {
    x <- matrix(sample.int(5L, 40L * 15L, replace = TRUE), nrow = 40L)
  })
  storage.mode(x) <- "double"
  expect_equal(cier_irv(x)$value, ref_irv(x), tolerance = 1e-10)
})

test_that("cier_irv$value equals the oracle when rows carry NAs (na.rm path)", {
  withr::with_seed(99, {
    x <- matrix(sample.int(5L, 40L * 15L, replace = TRUE), nrow = 40L)
  })
  storage.mode(x) <- "double"
  x[3L, c(1L, 5L)] <- NA      # 13 present
  x[10L, 1L:13L] <- NA        # 2 present (still scored)
  expect_equal(cier_irv(x)$value, ref_irv(x), tolerance = 1e-10)
})

# ---- Cross-package parity: careless::irv (1e-10) ----------------------------

test_that("cier_irv matches careless::irv on a random matrix at 1e-10", {
  skip_if_not_installed("careless")
  withr::with_seed(2026, {
    x <- matrix(sample.int(5L, 30L * 20L, replace = TRUE), nrow = 30L)
  })
  storage.mode(x) <- "double"
  expect_equal(cier_irv(x)$value, as.numeric(careless::irv(x)),
               tolerance = 1e-10)
})

test_that("cier_irv matches careless::irv on careless_dataset at 1e-10", {
  skip_if_not_installed("careless")
  x <- unname(as.matrix(careless::careless_dataset))
  storage.mode(x) <- "double"
  expect_equal(cier_irv(x)$value, as.numeric(careless::irv(x)),
               tolerance = 1e-10)
})

test_that("cier_irv matches careless::irv on NA-bearing rows incl. abstention", {
  # Substantiates the no-explicit-guard design: rowSds(na.rm = TRUE) must agree
  # with careless::irv on scored NA rows AND return the same NA where careless
  # abstains (n < 2). careless::irv defaults to na.rm = TRUE, so the only
  # divergence allowed is ulp-level matrixStats-vs-stats::sd reordering.
  skip_if_not_installed("careless")
  withr::with_seed(99, {
    x <- matrix(sample.int(5L, 40L * 15L, replace = TRUE), nrow = 40L)
  })
  storage.mode(x) <- "double"
  x[3L, c(1L, 5L)] <- NA      # 13 present (scored)
  x[10L, 1L:13L] <- NA        # 2 present  (scored, n >= 2 boundary)
  x[20L, 1L:14L] <- NA        # 1 present  (abstains -> NA in both)
  ours <- cier_irv(x)$value
  expect_equal(ours, as.numeric(careless::irv(x)), tolerance = 1e-10)
  expect_true(is.na(ours[[20L]]))   # abstention agrees with careless
})

# ---- Property / invariant + mutant-killers ----------------------------------

test_that("a constant row has IRV exactly 0, never NA (constant->NA mutant)", {
  expect_identical(suppressWarnings(cier_irv(matrix(rep(3, 10L), nrow = 1L)))$value,
                   0)
})

test_that("IRV is the SAMPLE sd (denominator n-1), not the population sd", {
  # c(1, 5): mean 3, ss = 8. sample sd = sqrt(8/1) = sqrt(8); a population-sd
  # mutant (divide by n) would return sqrt(8/2) = 2. sqrt(8) ~ 2.828 vs 2.
  expect_equal(suppressWarnings(cier_irv(matrix(c(1, 5), nrow = 1L)))$value, sqrt(8),
               tolerance = 1e-10)
})

test_that("na.rm is honoured: an NA-bearing row scores over present items", {
  # c(1, 5, NA): present (1, 5) -> sqrt(8). A mutant dropping na.rm = TRUE would
  # see an NA in the row and return NA for the whole respondent.
  out <- suppressWarnings(cier_irv(matrix(c(1, 5, NA), nrow = 1L)))$value
  expect_false(is.na(out))
  expect_equal(out, sqrt(8), tolerance = 1e-10)
})

test_that("scored IRV values are finite and non-negative", {
  withr::with_seed(7, {
    x <- matrix(sample.int(5L, 50L * 30L, replace = TRUE), nrow = 50L)
  })
  storage.mode(x) <- "double"
  v <- cier_irv(x)$value
  expect_true(all(is.finite(v)))
  expect_true(all(v >= 0))
})

test_that("direction is lower: low-variability rows flag, high ones do not", {
  # A constant row (sd 0, the global minimum) plus 20 random rows. At the
  # default 5th-percentile cutoff the constant row is below the cutoff and the
  # most-variable row is above it. A flag-upper mutant would invert both.
  withr::with_seed(11, {
    rnd <- matrix(sample.int(5L, 20L * 12L, replace = TRUE), nrow = 20L)
  })
  x <- rbind(rep(3, 12L), rnd)
  storage.mode(x) <- "double"
  out <- cier_irv(x)
  expect_true(out$flagged[[1L]])                       # constant row, sd 0
  expect_false(out$flagged[[which.max(out$value)]])    # most-variable row
  # The flag is exactly the lower comparator against the cutoff (NA-safe).
  expect_identical(out$flagged, !is.na(out$value) & out$value <= out$cutoff)
})

# ---- Edge cases -------------------------------------------------------------

test_that("an all-NA row abstains and an abstainer keeps rows aligned", {
  # Guards a row-indexing mutant: the abstaining row is in the middle, so value
  # and flagged must stay aligned to their respondents.
  x <- matrix(
    c(1, 2, 3, 4, 5,         # row 1: scored, sd sqrt(2.5)
      NA, NA, NA, NA, NA,    # row 2: abstains
      3, 3, 3, 3, 3),        # row 3: constant, sd 0
    nrow = 3L, byrow = TRUE
  )
  out <- suppressWarnings(cier_irv(x))
  expect_true(is.na(out$value[[2L]]))
  expect_true(is.na(out$flagged[[2L]]))
  expect_false(is.na(out$value[[1L]]))
  expect_identical(out$value[[3L]], 0)
})

test_that("a row with a single present value abstains (value NA, flagged NA)", {
  # Two scored rows give the percentile cutoff finite values to work with, so
  # this isolates the n = 1 abstention from the all-abstain warning path (which
  # the single-column test below covers).
  x <- matrix(
    c(1, 2, 3, 4, 5,          # scored
      5, 4, 3, 2, 1,          # scored
      3, NA, NA, NA, NA),     # one present value -> abstains
    nrow = 3L, byrow = TRUE
  )
  out <- suppressWarnings(cier_irv(x))
  expect_true(is.na(out$value[[3L]]))
  expect_true(is.na(out$flagged[[3L]]))
  expect_false(is.na(out$value[[1L]]))
  expect_false(is.na(out$value[[2L]]))
})

test_that("a single-column matrix abstains for every row and flags nobody", {
  # Each row has one item -> fewer than two present -> all abstain -> the
  # percentile cutoff has no finite values to work with: it warns and returns
  # NA, and an NA cutoff flags no one.
  expect_warning(
    out <- cier_irv(matrix(c(3, 4, 5), ncol = 1L)),
    class = "cier_warning_insufficient_items"
  )
  expect_true(all(is.na(out$value)))
  expect_true(is.na(out$cutoff))
  expect_true(all(is.na(out$flagged)))
})

test_that("a constant matrix abstains rather than flagging 100% (F02/F18)", {
  # The reported bug: a constant score distribution flagged EVERY respondent at
  # any fpr. With 25 rows the small-sample guard (D2) is satisfied, so this
  # isolates D1: the cutoff would flag everyone, so it abstains (NA + warning) and
  # flags nobody. Values are still the exact IRV scores (0); only the cutoff moves.
  expect_warning(out <- cier_irv(matrix(rep(2, 25L * 6L), nrow = 25L)),
                 class = "cier_warning_insufficient_items")
  expect_true(all(out$value == 0))
  expect_identical(out$cutoff, NA_real_)
  expect_false(any(out$flagged))                 # NA cutoff -> nobody flagged
})

test_that("a single respondent abstains rather than flagging itself (F37)", {
  # The reported bug: n = 1 flagged the only respondent even on a perfect score.
  # A 5%-tail cutoff needs >= 20 scored respondents (D2), so a lone respondent
  # abstains. Its value is still computed (the sd of 1..5); only the cutoff moves.
  expect_warning(out <- cier_irv(matrix(c(1, 2, 3, 4, 5), nrow = 1L)),
                 class = "cier_warning_insufficient_items")
  expect_false(is.na(out$value))
  expect_identical(out$cutoff, NA_real_)
  expect_false(out$flagged[[1L]])
})

test_that("on healthy continuous data the realised rate tracks fpr (no regression)", {
  # The guards must not perturb clean data: on continuous scores the percentile
  # cutoff still flags about fpr by construction (no tie mass, no degeneracy).
  withr::with_seed(123L,
                   x <- matrix(stats::runif(1000L * 12L, 1, 5), nrow = 1000L))
  for (f in c(0.01, 0.05, 0.10)) {
    fl <- cier_irv(x, fpr = f)$flagged
    expect_lt(abs(mean(fl) - f), 0.015)
  }
})

# ---- Cutoff: default, fpr override, NO-FLIP direction ------------------------

test_that("default cutoff is the lower-tail 5th percentile (NO-FLIP)", {
  withr::with_seed(5, {
    x <- matrix(stats::runif(60L * 12L, 1, 5), nrow = 60L)
  })
  out <- cier_irv(x)
  # Lower direction takes the fpr quantile directly (NOT 1 - fpr): the registry
  # stores the literal directional quantile and the kernel must not re-flip.
  expect_equal(out$cutoff,
               as.numeric(stats::quantile(out$value, 0.05,
                                          names = FALSE, type = 7L)),
               tolerance = 1e-12)
})

test_that("the fpr argument moves the percentile target", {
  withr::with_seed(5, {
    x <- matrix(stats::runif(60L * 12L, 1, 5), nrow = 60L)
  })
  out <- cier_irv(x, fpr = 0.10)
  expect_equal(out$cutoff,
               as.numeric(stats::quantile(out$value, 0.10,
                                          names = FALSE, type = 7L)),
               tolerance = 1e-12)
})

test_that("invalid fpr values are typed input errors", {
  x <- present_fixture()
  expect_error(cier_irv(x, fpr = 0), class = "cier_error_input")
  expect_error(cier_irv(x, fpr = 1), class = "cier_error_input")
  expect_error(cier_irv(x, fpr = -0.1), class = "cier_error_input")
  expect_error(cier_irv(x, fpr = NA_real_), class = "cier_error_input")
  expect_error(cier_irv(x, fpr = c(0.05, 0.1)), class = "cier_error_input")
  expect_error(cier_irv(x, fpr = "x"), class = "cier_error_input")
})

test_that("an absolute cutoff overrides the percentile and flags via direction", {
  # A literal IRV threshold instead of a percentile. Lower direction => values
  # at or below the cutoff flag. present_fixture sds are ~1.58, 0, ~2.19, 0.
  out <- cier_irv(present_fixture(), cutoff = 1)
  expect_identical(out$cutoff, 1)
  expect_identical(out$flagged, c(FALSE, TRUE, FALSE, TRUE))  # the two sd-0 rows
})

test_that("supplying both fpr and cutoff is a typed input error", {
  expect_error(cier_irv(present_fixture(), fpr = 0.1, cutoff = 0.5),
               class = "cier_error_input")
})

test_that("an invalid absolute cutoff is a typed input error", {
  x <- present_fixture()
  expect_error(cier_irv(x, cutoff = c(0.5, 1)), class = "cier_error_input")
  expect_error(cier_irv(x, cutoff = NA_real_), class = "cier_error_input")
  expect_error(cier_irv(x, cutoff = "x"), class = "cier_error_input")
  # IRV is a sample SD (>= 0); a negative threshold is out of domain and would
  # silently flag nobody, so it is rejected rather than accepted.
  expect_error(cier_irv(x, cutoff = -0.5), class = "cier_error_input")
})

test_that("cutoff = 0 is accepted and flags only exact straightliners", {
  # 0 is a legitimate lower-tail IRV threshold: value <= 0 flags sd-0 rows only.
  out <- cier_irv(present_fixture(), cutoff = 0)   # sds ~1.58, 0, ~2.19, 0
  expect_identical(out$cutoff, 0)
  expect_identical(out$flagged, c(FALSE, TRUE, FALSE, TRUE))
})

test_that("a non-matrix / non-numeric payload is a typed input error", {
  expect_error(cier_irv(1:10), class = "cier_error_input")
  expect_error(cier_irv(matrix(letters[1:6], nrow = 2L)),
               class = "cier_error_input")
})

# ---- print snapshot (locked, design-first; direction = lower) ---------------

test_that("print renders the locked cli summary (lower direction)", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    withr::with_seed(11, {
      x <- matrix(sample.int(5L, 30L * 6L, replace = TRUE), nrow = 30L)
    })
    storage.mode(x) <- "double"
    expect_snapshot(print(cier_irv(x)))
  })
})

test_that("print reports abstaining respondents on their own line", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    withr::with_seed(11, {
      x <- matrix(sample.int(5L, 29L * 6L, replace = TRUE), nrow = 29L)
    })
    storage.mode(x) <- "double"
    x <- rbind(x, rep(NA_real_, 6L))    # one abstaining respondent
    expect_snapshot(print(cier_irv(x)))
  })
})
