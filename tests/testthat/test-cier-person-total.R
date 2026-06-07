# Tests for cier_person_total() -- person-total correlation (r_pbis).
#
# Trust model: the independent oracle (ref_person_total) re-derives each
# respondent's correlation with the whole-sample item means via a per-row
# stats::cor() loop and never calls the production kernel (which evaluates the
# whole battery with vectorised masked sums); a separate hand-computed fixture
# pins exact +/-1 values that share no masked-sum machinery. The cross-package
# check pins parity with PerFit::r.pbis at the 1e-4 tolerance recorded in
# tests/reference/TOLERANCES.md (PerFit rounds its output to 4 dp; the
# correlation itself is exact). The property / mutant-killer block targets each
# mutant named in dev/restart/index-specs.md card 4 (item-rest instead of
# item-total, reverse-keying applied, scale-level means, flag upper vs lower).

source(test_path("..", "reference", "ref-person-total-donlon-fischer-1968.R"))

# Hand-computed fixture. The column means are (1.75, 2.25, 2.75, 3.25) -- a
# perfectly linear increasing sequence -- so an increasing row correlates +1
# with them and a decreasing row -1. Values are c(1, 1, -1, 1), worked by hand
# and confirmed against the oracle, with no masked-sum machinery in the
# expectation: a gross formula error cannot survive this.
hand_fixture <- function() {
  m <- rbind(c(1, 2, 3, 4), c(1, 2, 3, 4), c(4, 3, 2, 1), c(1, 2, 3, 4))
  storage.mode(m) <- "double"
  m
}

# A reproducible discrete fixture (n respondents x p items).
rand_matrix <- function(n = 30L, p = 8L, seed = 7L) {
  withr::with_seed(seed, {
    x <- matrix(sample.int(5L, n * p, replace = TRUE), nrow = n)
  })
  storage.mode(x) <- "double"
  x
}

# Independent item-REST (leave-one-out) reference, used only to prove the
# production statistic is item-TOTAL and not item-rest. Hand-rolled; never calls
# the kernel.
ref_item_rest <- function(x) {
  vapply(seq_len(nrow(x)), function(i) {
    rest_means <- colMeans(x[-i, , drop = FALSE], na.rm = TRUE)
    row <- x[i, ]
    ok <- !is.na(row)
    if (sum(ok) < 3L) {
      return(NA_real_)
    }
    suppressWarnings(stats::cor(as.numeric(row[ok]), rest_means[ok]))
  }, numeric(1L))
}

# ---- Schema -----------------------------------------------------------------

test_that("cier_person_total returns a list-based cier_index with the pinned schema", {
  out <- cier_person_total(hand_fixture())
  expect_s3_class(out, "cier_index")
  expect_type(out, "list")
  expect_identical(names(out),
                   c("value", "flagged", "method", "cutoff", "direction"))
  expect_type(out$value, "double")
  expect_type(out$flagged, "logical")
  expect_identical(length(out$value), 4L)
  expect_identical(length(out$flagged), 4L)
  expect_identical(out$method, "cier_person_total")
  expect_identical(out$direction, "lower")
  expect_type(out$cutoff, "double")
})

test_that("as.data.frame.cier_index returns the tidy per-respondent frame", {
  df <- as.data.frame(cier_person_total(hand_fixture()))
  expect_s3_class(df, "data.frame")
  expect_identical(names(df), c("value", "flagged"))
  expect_identical(nrow(df), 4L)
  expect_equal(df$value, c(1, 1, -1, 1), tolerance = 1e-12)
})

# ---- Independent oracle parity (1e-12) --------------------------------------

test_that("cier_person_total$value equals the hand-computed fixture exactly", {
  expect_equal(cier_person_total(hand_fixture())$value,
               c(1, 1, -1, 1), tolerance = 1e-12)
})

test_that("cier_person_total$value equals the oracle on a random complete matrix", {
  x <- rand_matrix(n = 40L, p = 15L, seed = 2026L)
  expect_equal(cier_person_total(x)$value, ref_person_total(x),
               tolerance = 1e-12)
})

test_that("cier_person_total$value equals the oracle when rows carry NAs", {
  x <- rand_matrix(n = 40L, p = 15L, seed = 99L)
  x[3L, c(1L, 5L)] <- NA       # 13 answered
  x[10L, 1L:12L] <- NA         # 3 answered (boundary, still scored)
  x[20L, 1L:13L] <- NA         # 2 answered (abstains)
  expect_equal(cier_person_total(x)$value, ref_person_total(x),
               tolerance = 1e-12)
})

# ---- Cross-package parity: PerFit::r.pbis (1e-4) ----------------------------

test_that("cier_person_total matches PerFit::r.pbis (cross-package, 4-dp)", {
  # GENUINE cross-package parity: PerFit::r.pbis is an independent point-biserial
  # implementation. The person-total correlation is invariant to an additive
  # shift, so our raw 1/2 coding agrees with PerFit's zero-based 0/1 coding; the
  # only divergence is PerFit's 4-dp output rounding (observed ~4.6e-5). See
  # tests/reference/TOLERANCES.md.
  skip_if_not_installed("PerFit")
  withr::with_seed(101L, {
    m <- matrix(sample.int(2L, 40L * 10L, replace = TRUE), 40L, 10L)
  })
  storage.mode(m) <- "double"
  ours <- cier_person_total(m)$value
  ref_matrix <- m - 1L
  storage.mode(ref_matrix) <- "integer"
  ref <- NULL
  invisible(utils::capture.output(ref <- PerFit::r.pbis(matrix = ref_matrix)))
  expect_equal(ours, as.numeric(ref$PFscores$PFscores), tolerance = 1e-4)
})

# ---- Property / invariant + mutant-killers ----------------------------------

test_that("the statistic is item-TOTAL, not item-rest (leave-one-out)", {
  # On a small sample the item-total and item-rest correlations diverge sharply
  # (here by ~0.66). We must match the item-total oracle to 1e-12 and be far from
  # the item-rest values; an item-rest mutant would invert both assertions.
  x <- rand_matrix(n = 5L, p = 6L, seed = 13L)
  expect_equal(cier_person_total(x)$value, ref_person_total(x),
               tolerance = 1e-12)
  expect_gt(max(abs(cier_person_total(x)$value - ref_item_rest(x))), 1e-6)
})

test_that("the wrapper is matrix-only: no item-metadata channel exists", {
  # Pins keying-insensitivity and whole-sample (not scale-level) means
  # structurally: there is no items / scale / reverse_keyed / categories
  # argument, so neither a reverse-keying nor a scale-means mutant can be wired
  # in. The whole-sample colMeans definition is enforced by the oracle parity.
  fmls <- names(formals(cier_person_total))
  expect_setequal(fmls, c("responses", "fpr", "cutoff"))
  expect_false(any(c("items", "scale", "reverse_keyed",
                     "categories") %in% fmls))
})

test_that("scored person-total values are finite and within [-1, 1]", {
  v <- cier_person_total(rand_matrix(n = 50L, p = 30L, seed = 7L))$value
  expect_true(all(is.finite(v)))
  expect_true(all(v >= -1 & v <= 1))
})

test_that("direction is lower: a careless (low-correlation) row flags, a consistent one does not", {
  # 20 rows that track an increasing item profile plus one decreasing
  # (anti-correlated) row. The anti row is the global minimum correlation -- the
  # most careless -- so it flags at the lower-tail cutoff and the most consistent
  # row does not. A flag-upper mutant would invert both.
  withr::with_seed(11L, {
    consistent <- t(vapply(seq_len(20L), function(i) {
      (1:12) + sample(-1:1, 12L, replace = TRUE)
    }, numeric(12L)))
  })
  x <- rbind(rev(1:12), consistent)        # row 1: decreasing -> negative cor
  storage.mode(x) <- "double"
  out <- cier_person_total(x)
  expect_true(out$flagged[[which.min(out$value)]])    # most careless row flags
  expect_false(out$flagged[[which.max(out$value)]])   # most consistent does not
  # The flag is exactly the lower comparator against the cutoff (NA-safe).
  expect_identical(out$flagged, !is.na(out$value) & out$value <= out$cutoff)
})

# ---- Edge cases -------------------------------------------------------------

test_that("a respondent with exactly two answered items abstains; three scores", {
  # cor() of two points is always +/-1, so without the k < 3 guard a 2-item row
  # would score a spurious +/-1. The guard forces NA at k = 2 and scores at k = 3
  # (the boundary): a k < 2 mutant scores row 1, a k < 4 mutant abstains row 2.
  x <- rand_matrix(n = 6L, p = 10L, seed = 5L)
  x[1L, 3L:10L] <- NA        # 2 answered -> abstains
  x[2L, 4L:10L] <- NA        # 3 answered -> scored
  out <- cier_person_total(x)
  expect_true(is.na(out$value[[1L]]))
  expect_false(is.na(out$value[[2L]]))
})

test_that("a constant (straightliner) row abstains (zero variance -> NA)", {
  x <- rand_matrix(n = 5L, p = 8L, seed = 4L)
  x[2L, ] <- 3               # constant straightliner -> zero variance
  out <- cier_person_total(x)
  expect_true(is.na(out$value[[2L]]))
  expect_false(is.na(out$value[[1L]]))
})

test_that("a flat item-mean profile abstains (zero variance on the item side)", {
  # Distinct from the k < 3 path: this circulant answers every item (k = 4) and
  # every row is non-constant, but the column means are all equal (2.5), so the
  # correlation has zero variance on the item-mean side and must be NA. A mutant
  # that only guards the respondent side (constant rows) would return NaN here.
  circ <- rbind(c(1, 2, 3, 4), c(2, 3, 4, 1), c(3, 4, 1, 2), c(4, 1, 2, 3))
  storage.mode(circ) <- "double"
  expect_warning(out <- cier_person_total(circ),
                 class = "cier_warning_insufficient_items")
  expect_true(all(is.na(out$value)))
})

test_that("a column answered by nobody yields finite scores (no NaN) and matches the oracle", {
  # An all-NA column gives a 0/0 = NaN item mean; the kernel forces it finite (0)
  # so it drops out of the masked products instead of poisoning every row. A
  # mutant that skips that guard turns every value NaN.
  x <- rand_matrix(n = 30L, p = 8L, seed = 7L)
  x[, 4L] <- NA              # item 4 answered by nobody
  out <- cier_person_total(x)
  expect_true(all(is.finite(out$value)))
  expect_equal(out$value, ref_person_total(x), tolerance = 1e-12)
})

test_that("an all-NA row abstains and keeps the remaining rows aligned", {
  # The abstaining row sits in the middle, so value/flagged must stay aligned to
  # their respondents (row-indexing mutant guard).
  x <- rand_matrix(n = 10L, p = 8L, seed = 4L)
  x[5L, ] <- NA
  out <- cier_person_total(x)
  expect_true(is.na(out$value[[5L]]))
  expect_true(is.na(out$flagged[[5L]]))
  expect_false(is.na(out$value[[1L]]))
  expect_false(is.na(out$value[[10L]]))
  expect_equal(out$value, ref_person_total(x), tolerance = 1e-12)
})

test_that("a two-column matrix abstains for every row and flags nobody", {
  # Each row has at most two answered items -> all abstain -> the percentile
  # cutoff has no finite values: it warns and returns NA, and an NA cutoff flags
  # no one.
  expect_warning(
    out <- cier_person_total(matrix(stats::runif(6L, 1, 5), ncol = 2L)),
    class = "cier_warning_insufficient_items"
  )
  expect_true(all(is.na(out$value)))
  expect_true(is.na(out$cutoff))
  expect_true(all(is.na(out$flagged)))
})

test_that("a non-matrix / non-numeric payload is a typed input error", {
  expect_error(cier_person_total(1:10), class = "cier_error_input")
  expect_error(cier_person_total(matrix(letters[1:6], nrow = 2L)),
               class = "cier_error_input")
})

# ---- Cutoff: default, fpr override, NO-FLIP direction ------------------------

test_that("default cutoff is the lower-tail 5th percentile (NO-FLIP)", {
  withr::with_seed(5L, {
    x <- matrix(stats::runif(60L * 12L, 1, 5), nrow = 60L)
  })
  out <- cier_person_total(x)
  # Lower direction takes the fpr quantile directly (NOT 1 - fpr): the registry
  # stores the literal directional quantile and the kernel must not re-flip.
  expect_equal(out$cutoff,
               as.numeric(stats::quantile(out$value, 0.05,
                                          names = FALSE, type = 7L)),
               tolerance = 1e-12)
})

test_that("the fpr argument moves the percentile target", {
  withr::with_seed(5L, {
    x <- matrix(stats::runif(60L * 12L, 1, 5), nrow = 60L)
  })
  out <- cier_person_total(x, fpr = 0.10)
  expect_equal(out$cutoff,
               as.numeric(stats::quantile(out$value, 0.10,
                                          names = FALSE, type = 7L)),
               tolerance = 1e-12)
})

test_that("invalid fpr values are typed input errors", {
  x <- hand_fixture()
  expect_error(cier_person_total(x, fpr = 0), class = "cier_error_input")
  expect_error(cier_person_total(x, fpr = 1), class = "cier_error_input")
  expect_error(cier_person_total(x, fpr = -0.1), class = "cier_error_input")
  expect_error(cier_person_total(x, fpr = NA_real_), class = "cier_error_input")
  expect_error(cier_person_total(x, fpr = c(0.05, 0.1)),
               class = "cier_error_input")
  expect_error(cier_person_total(x, fpr = "x"), class = "cier_error_input")
})

test_that("an absolute cutoff overrides the percentile and flags via direction", {
  # hand_fixture values are c(1, 1, -1, 1); lower direction => value <= 0 flags,
  # so only the anti-correlated (-1) respondent is flagged.
  out <- cier_person_total(hand_fixture(), cutoff = 0)
  expect_identical(out$cutoff, 0)
  expect_identical(out$flagged, c(FALSE, FALSE, TRUE, FALSE))
})

test_that("a respondent exactly at the cutoff is flagged (<= ties, not <)", {
  # apply_flag uses <= for the lower tail, so a correlation exactly equal to the
  # cutoff flags. Set the literal cutoff to one respondent's own score; a
  # strict-less (<) mutant would leave it unflagged.
  x <- rand_matrix(n = 30L, p = 8L, seed = 7L)
  v <- cier_person_total(x)$value
  k <- which.min(v)
  out <- cier_person_total(x, cutoff = v[[k]])
  expect_true(out$flagged[[k]])
})

test_that("supplying both fpr and cutoff is a typed input error", {
  expect_error(cier_person_total(hand_fixture(), fpr = 0.1, cutoff = 0.5),
               class = "cier_error_input")
})

test_that("an invalid absolute cutoff is a typed input error", {
  x <- hand_fixture()
  expect_error(cier_person_total(x, cutoff = c(0.5, 1)),
               class = "cier_error_input")
  expect_error(cier_person_total(x, cutoff = NA_real_),
               class = "cier_error_input")
  expect_error(cier_person_total(x, cutoff = "x"), class = "cier_error_input")
  # A person-total value is a correlation in [-1, 1]; a threshold outside that
  # range is degenerate (flags everyone or no one), so it is rejected.
  expect_error(cier_person_total(x, cutoff = -1.5), class = "cier_error_input")
  expect_error(cier_person_total(x, cutoff = 1.5), class = "cier_error_input")
})

# ---- print snapshot (locked, design-first; direction = lower) ---------------

test_that("print renders the locked cli summary (lower direction)", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    expect_snapshot(print(cier_person_total(rand_matrix(n = 30L, p = 12L,
                                                        seed = 11L))))
  })
})

test_that("print reports abstaining respondents on their own line", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    x <- rand_matrix(n = 29L, p = 12L, seed = 11L)
    x <- rbind(x, rep(NA_real_, 12L))      # one abstaining respondent
    expect_snapshot(print(cier_person_total(x)))
  })
})
