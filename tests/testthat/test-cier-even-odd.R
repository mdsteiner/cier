# Tests for cier_even_odd() -- even-odd consistency (Curran 2016; Jackson 1976).
#
# Trust model: the independent oracle (ref_even_odd) re-derives each respondent's
# -SB(r) step by step and never calls the production kernel; a separate set of
# hand-computed analytic fixtures pins exact -1 / +1 values that share no kernel
# machinery. The cross-package check pins BYTEWISE parity with careless::evenodd
# on no-reverse-key data (tolerance 0; see tests/reference/TOLERANCES.md); the
# NA / reverse-keying paths are covered by the 1e-12 oracle instead. The
# property / mutant-killer block targets each mutant named in
# dev/restart/index-specs.md card 5 (first/second-half split, missing
# Spearman-Brown, not negating, reverse-keying applied to all items or twice).

source(test_path("..", "reference", "ref-evenodd-curran-2016.R"))

# Scale-blocked `items`: `n_scales` scales of `per_scale` items each.
# `reverse_keyed` defaults to an alternating pattern; pass FALSE for analytic /
# parity fixtures that must isolate the correlation math from reverse-scoring.
blocked_items <- function(n_scales = 3L, per_scale = 4L, categories = 5L,
                          reverse_keyed = NULL) {
  scale <- rep(LETTERS[seq_len(n_scales)], each = per_scale)
  if (is.null(reverse_keyed)) {
    reverse_keyed <- rep(c(FALSE, TRUE), length.out = length(scale))
  }
  data.frame(scale = scale, reverse_keyed = reverse_keyed,
             categories = categories, stringsAsFactors = FALSE)
}

# Build the scale blocks INDEPENDENTLY of production (the oracle needs them but
# must not borrow scale_block_indices()).
blocks_from_scale <- function(items) {
  uniq <- unique(items$scale)
  lapply(uniq, function(s) which(items$scale == s))
}

rand_matrix <- function(n = 30L, p = 12L, seed = 7L) {
  withr::with_seed(seed, {
    x <- matrix(sample.int(5L, n * p, replace = TRUE), nrow = n)
  })
  storage.mode(x) <- "double"
  x
}

# Independently reverse-score (categories + 1) - x on reverse items.
prescore <- function(x, items) {
  rk <- items$reverse_keyed
  x[, rk] <- (items$categories[rk] + 1L) - x[, rk]
  x
}

# ---- Schema -----------------------------------------------------------------

test_that("cier_even_odd returns a list-based cier_index with the pinned schema", {
  out <- cier_even_odd(rand_matrix(20L, 12L, 1L),
                       blocked_items(3L, 4L, reverse_keyed = FALSE))
  expect_s3_class(out, "cier_index")
  expect_type(out, "list")
  expect_identical(names(out),
                   c("value", "flagged", "method", "cutoff", "direction"))
  expect_type(out$value, "double")
  expect_type(out$flagged, "logical")
  expect_identical(length(out$value), 20L)
  expect_identical(length(out$flagged), 20L)
  expect_identical(out$method, "cier_even_odd")
  expect_identical(out$direction, "upper")
  expect_type(out$cutoff, "double")
})

test_that("as.data.frame.cier_index returns the tidy per-respondent frame", {
  df <- as.data.frame(cier_even_odd(rand_matrix(20L, 12L, 2L),
                                    blocked_items(3L, 4L, reverse_keyed = FALSE)))
  expect_s3_class(df, "data.frame")
  expect_identical(names(df), c("value", "flagged"))
  expect_identical(nrow(df), 20L)
})

# ---- Independent oracle parity (1e-12) --------------------------------------

test_that("cier_even_odd$value equals the oracle on a complete matrix", {
  it <- blocked_items(5L, 5L, reverse_keyed = FALSE)
  x <- rand_matrix(30L, 25L, 2026L)
  expect_equal(cier_even_odd(x, it)$value,
               ref_even_odd(x, blocks_from_scale(it)), tolerance = 1e-12)
})

test_that("cier_even_odd$value equals the oracle when rows carry NAs", {
  it <- blocked_items(5L, 5L, reverse_keyed = FALSE)
  x <- rand_matrix(30L, 25L, 99L)
  x[3L, c(1L, 2L)] <- NA       # thins one scale's even/odd side (still a pair)
  x[10L, 1L:23L] <- NA         # only the last scale keeps a pair -> row NA
  expect_equal(cier_even_odd(x, it)$value,
               ref_even_odd(x, blocks_from_scale(it)), tolerance = 1e-12)
})

# ---- Analytic fixtures (hand, 1e-12): SB + negation -------------------------

test_that("a perfectly consistent respondent scores -1", {
  # Two scales x two items, within-scale identical, varying across scales.
  # No reverse, so the value isolates the Spearman-Brown correlation math.
  it <- blocked_items(2L, 2L, reverse_keyed = FALSE)
  x <- matrix(c(1, 1, 5, 5), nrow = 1L)
  storage.mode(x) <- "double"
  expect_equal(cier_even_odd(x, it)$value, -1, tolerance = 1e-12)
})

test_that("a perfectly inversely consistent respondent scores +1", {
  it <- blocked_items(2L, 2L, reverse_keyed = FALSE)
  x <- matrix(c(1, 5, 5, 1), nrow = 1L)
  storage.mode(x) <- "double"
  expect_equal(cier_even_odd(x, it)$value, 1, tolerance = 1e-12)
})

# ---- Split function + block builder ----------------------------------------

test_that("even_odd_split_fn assigns even positions first, odd positions second", {
  # cor() is symmetric, so even-as-first / odd-as-second reproduces the
  # even-vs-odd correlation. A first/second-half mutant would change this map.
  expect_identical(even_odd_split_fn(2L),
                   list(first_idx = 2L, second_idx = 1L))
  expect_identical(even_odd_split_fn(3L),
                   list(first_idx = 2L, second_idx = c(1L, 3L)))
  expect_identical(even_odd_split_fn(4L),
                   list(first_idx = c(2L, 4L), second_idx = c(1L, 3L)))
})

test_that("scale_block_indices groups columns by scale in first-appearance order", {
  it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
  expect_identical(scale_block_indices(it),
                   stats::setNames(list(1:4, 5:8, 9:12), c("A", "B", "C")))
})

test_that("scale_block_indices groups by scale LABEL, not by column position", {
  # Interleaved layout: a contiguous-chunk grouping mutant (ignoring the scale
  # column) would give A = {1,2,3}; the by-label rule gives A = {1,3,5}.
  it <- data.frame(scale = c("A", "B", "A", "B", "A", "B"))
  expect_identical(scale_block_indices(it),
                   stats::setNames(list(c(1L, 3L, 5L), c(2L, 4L, 6L)),
                                   c("A", "B")))
})

test_that("cier_even_odd groups by scale label on an interleaved layout (oracle parity)", {
  # End-to-end: pins the WRAPPER path, not just the helper. A contiguous-chunk
  # mutant (A = cols 1-4) diverges from the by-label oracle (A = cols 1,4,7,10).
  it <- data.frame(scale = rep(c("A", "B", "C"), times = 4L),
                   reverse_keyed = FALSE)
  x <- rand_matrix(20L, 12L, 21L)
  expect_equal(cier_even_odd(x, it)$value,
               ref_even_odd(x, blocks_from_scale(it)), tolerance = 1e-12)
})

# ---- Direction (upper, NO-FLIP) --------------------------------------------

test_that("direction is upper: a careless (high) row flags, a consistent one does not", {
  # A flag-lower mutant inverts both assertions.
  it <- blocked_items(5L, 4L, reverse_keyed = FALSE)
  out <- cier_even_odd(rand_matrix(40L, 20L, 11L), it)
  expect_true(out$flagged[[which.max(out$value)]])
  expect_false(out$flagged[[which.min(out$value)]])
  expect_identical(out$flagged, !is.na(out$value) & out$value >= out$cutoff)
})

test_that("default cutoff is the upper-tail 95th percentile (NO-FLIP)", {
  # Upper takes the 1 - fpr quantile directly; a double-flip mutant uses fpr.
  it <- blocked_items(5L, 4L, reverse_keyed = FALSE)
  out <- cier_even_odd(rand_matrix(60L, 20L, 5L), it)
  expect_equal(out$cutoff,
               as.numeric(stats::quantile(out$value[is.finite(out$value)],
                                          0.95, names = FALSE, type = 7L)),
               tolerance = 1e-12)
})

test_that("the fpr argument moves the percentile target", {
  it <- blocked_items(5L, 4L, reverse_keyed = FALSE)
  out <- cier_even_odd(rand_matrix(60L, 20L, 5L), it, fpr = 0.10)
  expect_equal(out$cutoff,
               as.numeric(stats::quantile(out$value[is.finite(out$value)],
                                          0.90, names = FALSE, type = 7L)),
               tolerance = 1e-12)
})

# ---- Reverse-keying (kills keying / double-reflection mutants) --------------

test_that("honouring reverse_keyed equals independently pre-scoring the items", {
  # Contract: the index reverse-scores raw responses internally, so self-keying
  # raw data must equal pre-scoring it and declaring no reverse items. Applying
  # the reflection twice (a mutant) would return the raw value, breaking this.
  it_rev <- blocked_items(4L, 4L)          # alternating reverse
  x <- rand_matrix(25L, 16L, 101L)
  it_fwd <- it_rev
  it_fwd$reverse_keyed <- FALSE
  honoured  <- cier_even_odd(x, it_rev)$value
  prescored <- cier_even_odd(prescore(x, it_rev), it_fwd)$value
  expect_equal(honoured, prescored, tolerance = 1e-12)
})

test_that("honoured reverse-keying also equals the oracle on pre-scored input", {
  it_rev <- blocked_items(4L, 4L)
  x <- rand_matrix(25L, 16L, 7L)
  expect_equal(cier_even_odd(x, it_rev)$value,
               ref_even_odd(prescore(x, it_rev), blocks_from_scale(it_rev)),
               tolerance = 1e-12)
})

test_that("apply_split_half_keying is a strict no-op without reverse items", {
  it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
  x <- rand_matrix(8L, 12L, 102L)
  expect_identical(apply_split_half_keying(x, it), x)
})

test_that("apply_split_half_keying aborts on a reverse item with NA categories", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, FALSE, TRUE, FALSE),
                   categories = c(5L, 5L, NA, 5L))
  x <- matrix(c(1, 2, 3, 4), nrow = 1L)
  storage.mode(x) <- "double"
  expect_error(apply_split_half_keying(x, it), class = "cier_error_input")
})

test_that("apply_split_half_keying aborts on a reverse item with NA min", {
  # Symmetric with the categories guard: a direct caller (bypassing check_items)
  # that supplies a reverse item with no scale base must error, not silently
  # reflect to NA.
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, FALSE, TRUE, FALSE),
                   categories = 5L, min = c(1L, 1L, NA, 1L))
  x <- matrix(c(1, 2, 3, 4), nrow = 1L)
  storage.mode(x) <- "double"
  expect_error(apply_split_half_keying(x, it), class = "cier_error_input")
})

# ---- Reverse-keying respects the response-scale base (min) ------------------

test_that("apply_split_half_keying reflects with the declared min: (min + max) - x", {
  # 0..4 scale (min=0, categories=5 -> max=4): reverse items reflect by 4 - x,
  # forward items untouched. A reflection that ignores min ((categories+1)-x)
  # would send item 2 to 6 and escape the 0..4 range.
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, TRUE),
                   categories = 5L, min = 0L)
  x <- matrix(c(0, 0, 4, 4), nrow = 1L)
  storage.mode(x) <- "double"
  expect_equal(as.numeric(apply_split_half_keying(x, it)),
               c(0, 4, 4, 0), tolerance = 1e-12)
})

test_that("apply_split_half_keying defaults min to 1 when no min column (backward compat)", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, TRUE), categories = 5L)
  x <- matrix(c(1, 1, 5, 5), nrow = 1L)
  storage.mode(x) <- "double"
  # min defaults to 1: reflect (1 + 5) - x = 6 - x.
  expect_equal(as.numeric(apply_split_half_keying(x, it)),
               c(1, 5, 5, 1), tolerance = 1e-12)
})

test_that("apply_split_half_keying handles per-item categories with a declared min", {
  # Per-item categories + a 0 base: each reverse item reflects by its OWN
  # (min + categories - 1) - x, guarding the vectorisation of max. A scalar
  # recycle of categories[1] would mis-reflect items 2 and 4.
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, TRUE),
                   categories = c(5L, 4L, 6L, 3L), min = 0L)
  x <- matrix(c(1, 2, 3, 1), nrow = 1L)
  storage.mode(x) <- "double"
  # item 2: cats 4, min 0 -> max 3 -> 3 - 2 = 1; item 4: cats 3, min 0 -> 2 - 1 = 1.
  expect_equal(as.numeric(apply_split_half_keying(x, it)),
               c(1, 1, 3, 1), tolerance = 1e-12)
})

test_that("responses outside the declared reverse-keying range are a typed error", {
  # A type-valid but WRONG declaration used to reflect to off-scale values and
  # silently corrupt the consistency score (flipping flags with no signal). The
  # keying step now cross-checks the observed range of every reverse-keyed
  # column against [min, min + categories - 1] and aborts naming the offenders
  # -- the same mistake the person-fit bridges already catch in
  # personfit_zero_base(). The classic trigger: 0-based data (0..4, five
  # options) declared categories = 5 but with the default min = 1, which would
  # reflect 0 -> 6 and 4 -> 2.
  x <- rand_matrix(20L, 12L, 3L) - 1     # 0..4 coding
  it <- blocked_items(3L, 4L)            # categories = 5, min defaults to 1
  expect_error(cier_even_odd(x, it), class = "cier_error_input")
  expect_error(cier_personal_reliability(x, it), class = "cier_error_input")
  # Declaring the true base scores cleanly (the guard keys on the declaration).
  it$min <- 0L
  expect_s3_class(cier_even_odd(x, it), "cier_index")
  # A declared categories SMALLER than the data is the other direction: 1..5
  # data with categories = 3 would reflect 5 -> -1.
  it2 <- blocked_items(3L, 4L, categories = 3L)
  expect_error(cier_even_odd(rand_matrix(20L, 12L, 3L), it2),
               class = "cier_error_input")
})

test_that("the range cross-check ignores forward items and all-NA reverse columns", {
  # Only reverse-keyed columns are reflected, so only they are checked: a
  # forward item may exceed the declared range without error (categories is not
  # read for it), and an all-NA reverse column has no observed range to violate.
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, TRUE),
                   categories = 5L)
  x <- rand_matrix(10L, 4L, 5L)
  x[, 1L] <- 9                           # forward item off the declared range
  x[, 2L] <- NA_real_                    # all-NA reverse column
  expect_s3_class(suppressWarnings(cier_even_odd(x, it)), "cier_index")
})

# A keying pattern that is NOT aligned to the even/odd partition: one reverse
# item per scale at a rotating within-scale position (1, 2, 3, 4). With the
# alternating default, every reverse item lands on an even position, so a
# min-ignoring reflection shifts the whole even-mean vector by a constant and the
# across-scale correlation is unchanged -- the end-to-end test would then pass
# for a buggy impl. The rotation breaks that location-invariance.
rotating_rev_items <- function() {
  rk <- rep(FALSE, 16L)
  rk[c(1L, 6L, 11L, 16L)] <- TRUE
  data.frame(scale = rep(LETTERS[1:4], each = 4L),
             reverse_keyed = rk, categories = 5L, stringsAsFactors = FALSE)
}

test_that("even-odd is invariant to the response-scale base when min is declared", {
  # The SAME respondents coded 1..5 (min=1) and 0..4 (min=0) must score
  # identically -- the reflection must use the declared base. A min-ignoring
  # (categories+1)-x reflection diverges by up to the full [-1, 1] range.
  it1 <- rotating_rev_items()              # min defaults to 1 (1..5 coding)
  x1 <- rand_matrix(30L, 16L, 303L)
  it0 <- it1
  it0$min <- 0L
  x0 <- x1 - 1L                            # SAME information, 0..4 coding
  expect_equal(cier_even_odd(x1, it1)$value,
               cier_even_odd(x0, it0)$value, tolerance = 1e-12)
})

test_that("omitting min reproduces an explicit min = 1 (default base, end-to-end)", {
  # Kills a default-to-0 / default-to-NA mutant: with the rotating keying the
  # default base genuinely changes the score, so no-min must equal min = 1.
  it_no_min <- rotating_rev_items()
  it_min1 <- it_no_min
  it_min1$min <- 1L
  x <- rand_matrix(30L, 16L, 7L)
  expect_equal(cier_even_odd(x, it_no_min)$value,
               cier_even_odd(x, it_min1)$value, tolerance = 1e-12)
})

test_that("off-midpoint straightliner is scored, not abstained, with reverse-keying", {
  # Reverse-scoring reflects a constant raw row into a NON-constant one unless it
  # sits at the scale midpoint, so -- unlike a forward-keyed battery -- the row is
  # scored rather than abstaining. Pins the documented @details behaviour (and
  # guards against the help page reverting to "a constant row always abstains").
  it <- rotating_rev_items()                    # reverse items at rotating positions
  filler <- rand_matrix(8L, 16L, 51L)           # so the percentile cutoff resolves
  x <- rbind(rep(5, 16L), rep(3, 16L), filler)  # row 1 off-midpoint, row 2 midpoint
  storage.mode(x) <- "double"
  out <- cier_even_odd(x, it)
  expect_false(is.na(out$value[[1L]]))          # off-midpoint constant -> finite score
  expect_true(is.na(out$value[[2L]]))           # midpoint constant -> reflects to itself -> NA
})

# ---- Optional metadata defaults --------------------------------------------

test_that("reverse_keyed and categories are optional when nothing is reverse-keyed", {
  # `items` carrying only `scale` must equal an explicit all-FALSE / categories
  # frame -- exercising the conditional-categories contract.
  it_scale_only <- data.frame(scale = rep(LETTERS[1:3], each = 4L))
  it_fwd <- blocked_items(3L, 4L, reverse_keyed = FALSE)
  x <- rand_matrix(20L, 12L, 3L)
  expect_equal(cier_even_odd(x, it_scale_only)$value,
               cier_even_odd(x, it_fwd)$value, tolerance = 1e-12)
})

# ---- Edge cases -------------------------------------------------------------

test_that("a single-item scale is skipped; remaining scales still score", {
  it <- data.frame(scale = c("A", "B", "B", "C", "C"),
                   reverse_keyed = FALSE, categories = 5L)
  x <- matrix(c(3, 1, 5, 2, 4,
                4, 2, 4, 1, 5), nrow = 2L, byrow = TRUE)
  storage.mode(x) <- "double"
  out <- cier_even_odd(x, it)
  expect_false(any(is.na(out$value)))     # scales B and C give two finite pairs
  expect_equal(out$value, ref_even_odd(x, blocks_from_scale(it)),
               tolerance = 1e-12)
})

test_that("a constant (straightliner) row abstains (zero variance -> NA)", {
  it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
  x <- rand_matrix(6L, 12L, 4L)
  x[2L, ] <- 3                            # constant -> identical even/odd means
  out <- cier_even_odd(x, it)
  expect_true(is.na(out$value[[2L]]))
  expect_false(is.na(out$value[[1L]]))
})

test_that("an all-NA row abstains and keeps the remaining rows aligned", {
  it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
  x <- rand_matrix(10L, 12L, 4L)
  x[5L, ] <- NA
  out <- cier_even_odd(x, it)
  expect_true(is.na(out$value[[5L]]))
  expect_true(is.na(out$flagged[[5L]]))
  expect_false(is.na(out$value[[1L]]))
  expect_false(is.na(out$value[[10L]]))
  expect_equal(out$value, ref_even_odd(x, blocks_from_scale(it)),
               tolerance = 1e-12)
})

test_that("when every respondent abstains the cutoff warns and flags nobody", {
  it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
  x <- matrix(3, nrow = 5L, ncol = 12L)  # all straightliners -> all NA
  storage.mode(x) <- "double"
  expect_warning(out <- cier_even_odd(x, it),
                 class = "cier_warning_insufficient_items")
  expect_true(all(is.na(out$value)))
  expect_true(is.na(out$cutoff))
  expect_true(all(is.na(out$flagged)))
})

# ---- Input validation (typed) ----------------------------------------------

test_that("fewer than two distinct scales is a typed input error", {
  it <- data.frame(scale = rep("A", 4L), reverse_keyed = FALSE, categories = 5L)
  x <- matrix(c(1, 2, 3, 4, 5, 4, 3, 2), nrow = 2L, byrow = TRUE)
  storage.mode(x) <- "double"
  expect_error(cier_even_odd(x, it), class = "cier_error_input")
})

test_that("a missing scale column is a typed input error", {
  it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
  it$scale <- NULL
  expect_error(cier_even_odd(rand_matrix(5L, 12L, 1L), it),
               class = "cier_error_input")
})

test_that("an items frame not aligned to the columns is a typed input error", {
  it <- blocked_items(3L, 4L, reverse_keyed = FALSE)   # 12 items
  expect_error(cier_even_odd(rand_matrix(5L, 10L, 1L), it),  # 10 columns
               class = "cier_error_input")
})

test_that("a non-matrix / non-numeric payload is a typed input error", {
  it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
  expect_error(cier_even_odd(1:10, it), class = "cier_error_input")
  expect_error(cier_even_odd(matrix(letters[1:12], nrow = 1L), it),
               class = "cier_error_input")
})

# ---- Cutoff overrides -------------------------------------------------------

test_that("an absolute cutoff overrides the percentile and flags via the upper tail", {
  # Two respondents: consistent (-1) and inverse (+1). cutoff 0 flags only +1.
  it <- blocked_items(2L, 2L, reverse_keyed = FALSE)
  x <- rbind(c(1, 1, 5, 5), c(1, 5, 5, 1))
  storage.mode(x) <- "double"
  out <- cier_even_odd(x, it, cutoff = 0)
  expect_identical(out$cutoff, 0)
  expect_identical(out$flagged, c(FALSE, TRUE))
})

test_that("a respondent exactly at the cutoff is flagged (>= ties, not >)", {
  it <- blocked_items(5L, 4L, reverse_keyed = FALSE)
  x <- rand_matrix(30L, 20L, 7L)
  v <- cier_even_odd(x, it)$value
  k <- which.max(v)
  out <- cier_even_odd(x, it, cutoff = v[[k]])
  expect_true(out$flagged[[k]])
})

test_that("invalid fpr values are typed input errors", {
  it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
  x <- rand_matrix(10L, 12L, 1L)
  expect_error(cier_even_odd(x, it, fpr = 0), class = "cier_error_input")
  expect_error(cier_even_odd(x, it, fpr = 1), class = "cier_error_input")
  expect_error(cier_even_odd(x, it, fpr = -0.1), class = "cier_error_input")
  expect_error(cier_even_odd(x, it, fpr = NA_real_), class = "cier_error_input")
  expect_error(cier_even_odd(x, it, fpr = c(0.05, 0.1)),
               class = "cier_error_input")
  expect_error(cier_even_odd(x, it, fpr = "x"), class = "cier_error_input")
})

test_that("an invalid absolute cutoff is a typed input error", {
  it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
  x <- rand_matrix(10L, 12L, 1L)
  expect_error(cier_even_odd(x, it, cutoff = c(0.5, 1)),
               class = "cier_error_input")
  expect_error(cier_even_odd(x, it, cutoff = NA_real_),
               class = "cier_error_input")
  expect_error(cier_even_odd(x, it, cutoff = "x"), class = "cier_error_input")
  # An even-odd value is a (negated, SB-corrected) correlation in [-1, 1];
  # a threshold outside that range flags everyone or no one.
  expect_error(cier_even_odd(x, it, cutoff = -1.5), class = "cier_error_input")
  expect_error(cier_even_odd(x, it, cutoff = 1.5), class = "cier_error_input")
})

test_that("supplying both fpr and cutoff is a typed input error", {
  it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
  expect_error(cier_even_odd(rand_matrix(10L, 12L, 1L), it,
                             fpr = 0.1, cutoff = 0.5),
               class = "cier_error_input")
})

# ---- Cross-package parity: bytewise vs careless (no-reverse-key data) -------

test_that("cier_even_odd matches careless::evenodd bytewise on no-reverse data", {
  skip_if_not_installed("careless")
  raw <- careless::careless_dataset
  responses <- unname(as.matrix(raw))
  storage.mode(responses) <- "double"
  # 10 contiguous scales of 5 items; no reverse-keyed items (careless::evenodd
  # does not reverse-key), so categories are not even required here.
  it <- data.frame(scale = rep(paste0("s", LETTERS[1:10]), each = 5L),
                   reverse_keyed = FALSE, stringsAsFactors = FALSE)
  ours <- cier_even_odd(responses, it)$value
  theirs <- suppressWarnings(careless::evenodd(raw, factors = rep(5L, 10L)))
  expect_equal(ours, as.numeric(theirs), tolerance = 0)
})

# ---- print snapshot (locked; shared upper-direction format) -----------------

test_that("print renders the locked cli summary (upper direction)", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
    expect_snapshot(print(cier_even_odd(rand_matrix(11L, 12L, 11L), it)))
  })
})

test_that("print reports abstaining respondents on their own line", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
    x <- rand_matrix(11L, 12L, 11L)
    x <- rbind(x, rep(NA_real_, 12L))
    expect_snapshot(print(cier_even_odd(x, it)))
  })
})
