# Tests for cier_even_odd() -- even-odd consistency (Curran 2016; Jackson 1976).
#
# Trust model: oracle ref_even_odd re-derives each respondent's -SB(r) step by
# step, never calling the kernel; hand-computed analytic fixtures pin exact
# -1 / +1 values sharing no kernel machinery. The cross-package check pins parity
# with careless::evenodd on no-reverse-key data at 1e-12 (observed 0 on this
# integer fixture -- masked-sum half-means and careless's base mean() round
# identically -- but the contract carries a forward-looking ulp margin; the
# NA / reverse-keying paths are covered by the 1e-12 oracle instead). The
# property / regression block targets each documented mutant (first/second-half
# split, missing Spearman-Brown, not negating, reverse-keying applied to all
# items or twice).

source(test_path("..", "reference", "ref-evenodd-curran-2016.R"))

# blocked_items() / blocks_from_scale() / prescore() fixtures are shared with
# test-cier-personal-reliability.R via helper-split-half.R. rand_matrix() stays local
# (the name is reused with different defaults across other index tests).
rand_matrix <- function(n = 30L, p = 12L, seed = 7L) {
  withr::with_seed(seed, {
    x <- matrix(sample.int(5L, n * p, replace = TRUE), nrow = n)
  })
  storage.mode(x) <- "double"
  x
}

# ---- Schema -----------------------------------------------------------------

test_that("cier_even_odd returns the pinned cier_index schema", {
  # Small/saturated fixtures trip the percentile-cutoff degeneracy guard; assert the
  # score, not the flag, so the (correct) warning is muffled.
  it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
  out <- suppressWarnings(cier_even_odd(rand_matrix(20L, 12L, 1L), it))
  expect_cier_index_schema(out, "cier_even_odd", "upper", 20L)
})

# ---- Independent oracle parity (1e-12) --------------------------------------

test_that("cier_even_odd$value equals the oracle on a complete matrix", {
  it <- blocked_items(5L, 5L, reverse_keyed = FALSE)
  x <- rand_matrix(30L, 25L, 2026L)
  expect_equal(suppressWarnings(cier_even_odd(x, it))$value,
               ref_even_odd(x, blocks_from_scale(it)), tolerance = 1e-12)
})

test_that("cier_even_odd$value equals the oracle when rows carry NAs", {
  it <- blocked_items(5L, 5L, reverse_keyed = FALSE)
  x <- rand_matrix(30L, 25L, 99L)
  x[3L, c(1L, 2L)] <- NA       # thins one scale's even/odd side (still a pair)
  x[10L, 1L:23L] <- NA         # only the last scale keeps a pair -> row NA
  expect_equal(suppressWarnings(cier_even_odd(x, it))$value,
               ref_even_odd(x, blocks_from_scale(it)), tolerance = 1e-12)
})

# ---- Analytic fixtures (hand, 1e-12): SB + negation -------------------------

test_that("a perfectly consistent respondent scores -1", {
  # Two scales x two items, within-scale identical, varying across scales. No
  # reverse, so the value isolates the Spearman-Brown correlation math.
  it <- blocked_items(2L, 2L, reverse_keyed = FALSE)
  x <- matrix(c(1, 1, 5, 5), nrow = 1L)
  storage.mode(x) <- "double"
  expect_equal(suppressWarnings(cier_even_odd(x, it))$value, -1, tolerance = 1e-12)
})

test_that("a perfectly inversely consistent respondent scores +1", {
  it <- blocked_items(2L, 2L, reverse_keyed = FALSE)
  x <- matrix(c(1, 5, 5, 1), nrow = 1L)
  storage.mode(x) <- "double"
  expect_equal(suppressWarnings(cier_even_odd(x, it))$value, 1, tolerance = 1e-12)
})

# ---- Split function + block builder ----------------------------------------

test_that("even_odd_split_fn assigns even positions first, odd positions second", {
  # cor() is symmetric, so even-as-first / odd-as-second reproduces the even-vs-odd
  # correlation. A first/second-half mutant would change this map.
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
  # End-to-end pins the WRAPPER path, not just the helper. A contiguous-chunk
  # mutant (A = cols 1-4) diverges from the by-label oracle (A = cols 1,4,7,10).
  it <- data.frame(scale = rep(c("A", "B", "C"), times = 4L),
                   reverse_keyed = FALSE)
  x <- rand_matrix(20L, 12L, 21L)
  expect_equal(suppressWarnings(cier_even_odd(x, it))$value,
               ref_even_odd(x, blocks_from_scale(it)), tolerance = 1e-12)
})

# ---- Direction (upper, NO-FLIP) --------------------------------------------

test_that("direction is upper: a careless (high) row flags, a consistent one does not", {
  # A flag-lower mutant inverts both assertions.
  it <- blocked_items(5L, 4L, reverse_keyed = FALSE)
  out <- suppressWarnings(cier_even_odd(rand_matrix(40L, 20L, 11L), it))
  expect_true(out$flagged[[which.max(out$value)]])
  expect_false(out$flagged[[which.min(out$value)]])
  expect_identical(out$flagged, !is.na(out$value) & out$value >= out$cutoff)
})

test_that("default cutoff is the upper-tail 95th percentile (NO-FLIP)", {
  # Upper takes the 1 - fpr quantile directly; a double-flip mutant uses fpr.
  it <- blocked_items(5L, 4L, reverse_keyed = FALSE)
  out <- suppressWarnings(cier_even_odd(rand_matrix(60L, 20L, 5L), it))
  expect_equal(out$cutoff,
               as.numeric(stats::quantile(out$value[is.finite(out$value)],
                                          0.95, names = FALSE, type = 7L)),
               tolerance = 1e-12)
})

test_that("the fpr argument moves the percentile target", {
  it <- blocked_items(5L, 4L, reverse_keyed = FALSE)
  out <- suppressWarnings(cier_even_odd(rand_matrix(60L, 20L, 5L), it, fpr = 0.10))
  expect_equal(out$cutoff,
               as.numeric(stats::quantile(out$value[is.finite(out$value)],
                                          0.90, names = FALSE, type = 7L)),
               tolerance = 1e-12)
})

# ---- Reverse-keying (kills keying / double-reflection mutants) --------------

test_that("honouring reverse_keyed equals independently pre-scoring the items", {
  # The index reverse-scores raw responses internally, so self-keying raw data must
  # equal pre-scoring it and declaring no reverse items. Reflecting twice (a mutant)
  # returns the raw value, breaking this.
  it_rev <- blocked_items(4L, 4L)          # alternating reverse
  x <- rand_matrix(25L, 16L, 101L)
  it_fwd <- it_rev
  it_fwd$reverse_keyed <- FALSE
  honoured  <- suppressWarnings(cier_even_odd(x, it_rev))$value
  prescored <- suppressWarnings(cier_even_odd(prescore(x, it_rev), it_fwd))$value
  expect_equal(honoured, prescored, tolerance = 1e-12)
})

test_that("honoured reverse-keying also equals the oracle on pre-scored input", {
  it_rev <- blocked_items(4L, 4L)
  x <- rand_matrix(25L, 16L, 7L)
  expect_equal(suppressWarnings(cier_even_odd(x, it_rev))$value,
               ref_even_odd(prescore(x, it_rev), blocks_from_scale(it_rev)),
               tolerance = 1e-12)
})

test_that("apply_split_half_keying is a strict no-op without reverse items", {
  it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
  x <- rand_matrix(8L, 12L, 102L)
  expect_identical(apply_split_half_keying(x, it), x)
})

test_that("apply_split_half_keying aborts on a reverse item with NA max", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, FALSE, TRUE, FALSE),
                   max = c(5L, 5L, NA, 5L))
  x <- matrix(c(1, 2, 3, 4), nrow = 1L)
  storage.mode(x) <- "double"
  expect_error(apply_split_half_keying(x, it), class = "cier_error_input")
})

test_that("apply_split_half_keying aborts on a reverse item with NA min", {
  # Symmetric with the max guard: a direct caller (bypassing check_items) supplying
  # a reverse item with no scale base must error, not silently reflect to NA.
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, FALSE, TRUE, FALSE),
                   max = 5L, min = c(1L, 1L, NA, 1L))
  x <- matrix(c(1, 2, 3, 4), nrow = 1L)
  storage.mode(x) <- "double"
  expect_error(apply_split_half_keying(x, it), class = "cier_error_input")
})

# ---- Reverse-keying respects the response-scale base (min) ------------------

test_that("apply_split_half_keying reflects with the declared min: (min + max) - x", {
  # 0..4 scale (min=0, max=4): reverse items reflect by 4 - x, forward untouched.
  # A min-ignoring ((max+1)-x) reflection sends item 2 to 5 (off the 0..4 range);
  # treating max as a category COUNT (min + (min + max - 1) - x = 3 - x) sends item
  # 2 to 3, not 4. The exact-value pin kills both.
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, TRUE),
                   max = 4L, min = 0L)
  x <- matrix(c(0, 0, 4, 4), nrow = 1L)
  storage.mode(x) <- "double"
  expect_equal(as.numeric(apply_split_half_keying(x, it)),
               c(0, 4, 4, 0), tolerance = 1e-12)
})

test_that("apply_split_half_keying defaults min to 1 when no min column (backward compat)", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, TRUE), max = 5L)
  x <- matrix(c(1, 1, 5, 5), nrow = 1L)
  storage.mode(x) <- "double"
  # min defaults to 1: reflect (1 + 5) - x = 6 - x.
  expect_equal(as.numeric(apply_split_half_keying(x, it)),
               c(1, 5, 5, 1), tolerance = 1e-12)
})

test_that("apply_split_half_keying handles per-item max with a declared min", {
  # Per-item max + a 0 base: each reverse item reflects by its OWN (min + max) - x,
  # guarding the vectorisation. A scalar recycle of max[1] mis-reflects items 2, 4.
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, TRUE),
                   max = c(4L, 3L, 5L, 2L), min = 0L)
  x <- matrix(c(1, 2, 3, 1), nrow = 1L)
  storage.mode(x) <- "double"
  # item 2: max 3, min 0 -> 3 - 2 = 1; item 4: max 2, min 0 -> 2 - 1 = 1.
  expect_equal(as.numeric(apply_split_half_keying(x, it)),
               c(1, 1, 3, 1), tolerance = 1e-12)
})

test_that("responses outside the declared reverse-keying range are a typed error", {
  # A type-valid but WRONG declaration used to reflect to off-scale values and
  # silently corrupt the consistency score (flipping flags with no signal). The
  # keying step now cross-checks every reverse-keyed column's observed range
  # against [min, max] and aborts naming the offenders -- the same mistake
  # personfit_zero_base() already catches. Classic trigger: 0-based data (0..4)
  # declared max = 5 with default min = 1, reflecting 0 -> 6 and 4 -> 2.
  x <- rand_matrix(20L, 12L, 3L) - 1     # 0..4 coding
  it <- blocked_items(3L, 4L)            # max = 5, min defaults to 1
  expect_error(cier_even_odd(x, it), class = "cier_error_input")
  expect_error(cier_personal_reliability(x, it), class = "cier_error_input")
  # Declaring the true range scores cleanly (the guard keys on the declaration).
  it$min <- 0L
  it$max <- 4L
  expect_s3_class(suppressWarnings(cier_even_odd(x, it)), "cier_index")
  # A declared max SMALLER than the data is the other direction: 1..5 data
  # with max = 3 would reflect 5 -> -1.
  it2 <- blocked_items(3L, 4L, max = 3L)
  expect_error(cier_even_odd(rand_matrix(20L, 12L, 3L), it2),
               class = "cier_error_input")
})

test_that("out-of-range responses are caught on forward AND reverse items", {
  # The range cross-check now covers EVERY item declaring a min/max, not only
  # reverse-keyed ones -- a 99-style sentinel must be caught in any column. A
  # forward item could once exceed its declared max and silently corrupt the
  # consistency score. This replaces the old test that pinned forward items as
  # exempt.
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, TRUE),
                   max = 5L)
  x <- rand_matrix(10L, 4L, 5L)
  x[, 1L] <- 9                           # FORWARD item ABOVE the declared max = 5
  expect_error(cier_even_odd(x, it), class = "cier_error_input")
  expect_error(cier_personal_reliability(x, it), class = "cier_error_input")
  # The other direction: a FORWARD item BELOW the declared min (default 1) must
  # also abort -- the check covers both bounds, not just max.
  x_lo <- rand_matrix(10L, 4L, 5L)
  x_lo[, 1L] <- 0                        # below min = 1 on a forward item
  expect_error(cier_even_odd(x_lo, it), class = "cier_error_input")
})

test_that("the range cross-check still ignores all-NA columns and undeclared-max items", {
  # Only items declaring a max are checked: a NA (undeclared) `max` is exempt, and
  # an all-NA column has no observed range to violate. (`min` defaults to 1, so a
  # non-NA max is the gate.)
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, TRUE),
                   max = c(NA, 5L, NA, 5L))   # forward items 1, 3 declare no max
  x <- rand_matrix(10L, 4L, 5L)
  x[, 1L] <- 9                           # forward item with no declared max -> ok
  x[, 2L] <- NA_real_                    # all-NA reverse column -> ok
  expect_s3_class(suppressWarnings(cier_even_odd(x, it)), "cier_index")
})

test_that("a forward item declaring max but an explicit NA min is still upper-checked", {
  # The validator permits an explicit NA `min` on a forward item (a `min` column
  # declared for a 0-based reverse item, left NA on forward items) while keeping a
  # declared `max`. The upper bound must still catch a sentinel: a 99 must not slip
  # through just because that item's `min` is NA.
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, TRUE),
                   max = 5L, min = c(NA, 1L, NA, 1L))
  x <- rand_matrix(10L, 4L, 5L)
  x[, 1L] <- 99                          # forward item: declares max 5, min NA
  expect_error(cier_even_odd(x, it), class = "cier_error_input")
})

test_that("cier_even_odd informs once when items omits reverse_keyed", {
  x <- rand_matrix(10L, 4L, 5L)
  it <- data.frame(scale = rep(c("A", "B"), each = 2L), max = 5L)   # no reverse_keyed
  expect_message(suppressWarnings(cier_even_odd(x, it)),
                 class = "cier_message_forward_keyed")
  # An explicit reverse_keyed = FALSE is the deliberate forward-keyed case: silent.
  it2 <- data.frame(scale = rep(c("A", "B"), each = 2L),
                    reverse_keyed = FALSE, max = 5L)
  expect_no_message(suppressWarnings(cier_even_odd(x, it2)))
})

test_that("reordering items relative to named responses is caught, not silently scored", {
  # A reordered items frame silently reshuffles scale / reverse-keying. With
  # matching identifiers it scores; once the `item` ids no longer line up with the
  # response columns it is a typed error.
  x <- rand_matrix(10L, 4L, 5L)
  colnames(x) <- c("q1", "q2", "q3", "q4")
  it <- data.frame(item = c("q1", "q2", "q3", "q4"),
                   scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, TRUE), max = 5L)
  expect_s3_class(suppressWarnings(cier_even_odd(x, it)), "cier_index")
  it_reordered <- it[c(1L, 3L, 2L, 4L), ]
  expect_error(cier_even_odd(x, it_reordered), class = "cier_error_input")
})

# A keying pattern NOT aligned to the even/odd partition: one reverse item per
# scale at a rotating within-scale position (1, 2, 3, 4). Under the alternating
# default every reverse item lands on an even position, so a min-ignoring
# reflection shifts the whole even-mean vector by a constant, leaving the
# across-scale correlation unchanged -- the end-to-end test would then pass for a
# buggy impl. The rotation breaks that location-invariance.
rotating_rev_items <- function() {
  rk <- rep(FALSE, 16L)
  rk[c(1L, 6L, 11L, 16L)] <- TRUE
  data.frame(scale = rep(LETTERS[1:4], each = 4L),
             reverse_keyed = rk, max = 5L, stringsAsFactors = FALSE)
}

test_that("even-odd is invariant to the response-scale base when min is declared", {
  # The SAME respondents coded 1..5 (min=1, max=5) and 0..4 (min=0, max=4) must
  # score identically -- the reflection must use the declared base. A min-ignoring
  # (max+1)-x reflection diverges by up to the full [-1, 1] range.
  it1 <- rotating_rev_items()              # min defaults to 1 (1..5 coding)
  x1 <- rand_matrix(30L, 16L, 303L)
  it0 <- it1
  it0$min <- 0L
  it0$max <- 4L
  x0 <- x1 - 1L                            # SAME information, 0..4 coding
  expect_equal(suppressWarnings(cier_even_odd(x1, it1))$value,
               suppressWarnings(cier_even_odd(x0, it0))$value, tolerance = 1e-12)
})

test_that("omitting min reproduces an explicit min = 1 (default base, end-to-end)", {
  # Kills a default-to-0 / default-to-NA mutant: with the rotating keying the
  # default base genuinely changes the score, so no-min must equal min = 1.
  it_no_min <- rotating_rev_items()
  it_min1 <- it_no_min
  it_min1$min <- 1L
  x <- rand_matrix(30L, 16L, 7L)
  expect_equal(suppressWarnings(cier_even_odd(x, it_no_min))$value,
               suppressWarnings(cier_even_odd(x, it_min1))$value, tolerance = 1e-12)
})

test_that("off-midpoint straightliner is scored, not abstained, with reverse-keying", {
  # Reverse-scoring reflects a constant raw row into a NON-constant one unless it
  # sits at the scale midpoint, so -- unlike a forward-keyed battery -- the row is
  # scored, not abstained. Pins the documented @details behaviour (guards the help
  # page reverting to "a constant row always abstains").
  it <- rotating_rev_items()                    # reverse items at rotating positions
  filler <- rand_matrix(8L, 16L, 51L)           # so the percentile cutoff resolves
  x <- rbind(rep(5, 16L), rep(3, 16L), filler)  # row 1 off-midpoint, row 2 midpoint
  storage.mode(x) <- "double"
  out <- suppressWarnings(cier_even_odd(x, it))
  expect_false(is.na(out$value[[1L]]))          # off-midpoint constant -> finite score
  expect_true(is.na(out$value[[2L]]))           # midpoint constant -> reflects to itself -> NA
})

# ---- Optional metadata defaults --------------------------------------------

test_that("reverse_keyed and max are optional when nothing is reverse-keyed", {
  # `items` carrying only `scale` must equal an explicit all-FALSE / max frame --
  # exercising the conditional-max contract.
  it_scale_only <- data.frame(scale = rep(LETTERS[1:3], each = 4L))
  it_fwd <- blocked_items(3L, 4L, reverse_keyed = FALSE)
  x <- rand_matrix(20L, 12L, 3L)
  expect_equal(suppressWarnings(cier_even_odd(x, it_scale_only))$value,
               suppressWarnings(cier_even_odd(x, it_fwd))$value, tolerance = 1e-12)
})

# ---- Edge cases -------------------------------------------------------------

test_that("a single-item scale is skipped; remaining scales still score", {
  it <- data.frame(scale = c("A", "B", "B", "C", "C"),
                   reverse_keyed = FALSE, max = 5L)
  x <- matrix(c(3, 1, 5, 2, 4,
                4, 2, 4, 1, 5), nrow = 2L, byrow = TRUE)
  storage.mode(x) <- "double"
  out <- suppressWarnings(cier_even_odd(x, it))
  expect_false(any(is.na(out$value)))     # scales B and C give two finite pairs
  expect_equal(out$value, ref_even_odd(x, blocks_from_scale(it)),
               tolerance = 1e-12)
})

test_that("a constant (straightliner) row abstains (zero variance -> NA)", {
  it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
  x <- rand_matrix(6L, 12L, 4L)
  x[2L, ] <- 3                            # constant -> identical even/odd means
  out <- suppressWarnings(cier_even_odd(x, it))
  expect_true(is.na(out$value[[2L]]))
  expect_false(is.na(out$value[[1L]]))
})

test_that("an all-NA row abstains and keeps the remaining rows aligned", {
  it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
  x <- rand_matrix(10L, 12L, 4L)
  x[5L, ] <- NA
  out <- suppressWarnings(cier_even_odd(x, it))
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

# ---- Two-scale degeneracy ---------------------------------------------------
# With exactly two scorable scale blocks (>= 2 items each) the across-block
# consistency correlation is over two points, so it is +/-1 by construction -- a
# degenerate score min_scales = 2 used to admit silently. The wrapper now warns
# (typed cier_warning_two_scale_consistency) from a purely structural trigger
# (block sizes), without changing the statistic.

# A 2-scale x 2-item battery: `n_consistent` perfectly consistent rows (-1) and
# `n_inverse` perfectly inverse rows (+1), so every finite score is exactly +/-1.
two_scale_pm1 <- function(n_consistent, n_inverse) {
  x <- rbind(
    matrix(rep(c(1, 1, 5, 5), n_consistent), ncol = 4L, byrow = TRUE),
    matrix(rep(c(1, 5, 5, 1), n_inverse),    ncol = 4L, byrow = TRUE)
  )
  storage.mode(x) <- "double"
  x
}

test_that("two scorable scale blocks warn and score the degenerate +/-1 (even-odd)", {
  # Regression: a 2-scale battery flags on a binary +/-1 score. The silent
  # degenerate flagging now carries a typed warning; cutoff = 0 isolates this
  # warning from the percentile path.
  it <- blocked_items(2L, 2L, reverse_keyed = FALSE)
  x <- two_scale_pm1(12L, 8L)                    # 12 consistent (-1), 8 inverse (+1)
  expect_warning(out <- cier_even_odd(x, it, cutoff = 0),
                 class = "cier_warning_two_scale_consistency",
                 regexp = "scorable")
  expect_true(all(out$value %in% c(-1, 1)))      # a +/-1 point mass, no gradation
  expect_identical(sum(out$flagged), 8L)         # the 8 inverse rows flag at cutoff 0
})

test_that("the two-scale warning keys on scorable blocks, not scale labels (even-odd)", {
  # A one-item scale is skipped, so three labels with a singleton scale leave
  # exactly two scorable blocks and must warn; three full scales must not (kills a
  # label-count / >= 2 / length(blocks) trigger).
  it_two <- data.frame(scale = c("A", "B", "B", "C", "C"), reverse_keyed = FALSE)
  expect_warning(cier_even_odd(rand_matrix(20L, 5L, 11L), it_two, cutoff = 0),
                 class = "cier_warning_two_scale_consistency")
  it_three <- blocked_items(3L, 4L, reverse_keyed = FALSE)    # 3 scorable blocks
  expect_no_warning(cier_even_odd(rand_matrix(20L, 12L, 7L), it_three, cutoff = 0),
                    class = "cier_warning_two_scale_consistency")
})

test_that("a single scorable scale block abstains without the two-scale warning", {
  # Two labels but one is a singleton -> one scorable block -> every row abstains
  # (< 2 finite half-mean pairs). This is the insufficient-blocks case, NOT the
  # two-scale degeneracy, so the warning stays silent (separates == 2 from <= 2).
  # cutoff = 0 keeps the percentile abstention warning out of the way.
  it_one <- data.frame(scale = c("A", "B", "B"), reverse_keyed = FALSE)
  out <- expect_no_warning(
    cier_even_odd(rand_matrix(20L, 3L, 9L), it_one, cutoff = 0),
    class = "cier_warning_two_scale_consistency"
  )
  expect_true(all(is.na(out$value)))
})

test_that("the two-scale warning fires once per call, not once per respondent", {
  # The helper is called once at wrapper level, so a 20-respondent battery warns
  # exactly once. A per-respondent row-loop mutant would signal 20 times while
  # still passing expect_warning.
  it <- blocked_items(2L, 2L, reverse_keyed = FALSE)
  ww <- capture_warnings(cier_even_odd(two_scale_pm1(12L, 8L), it, cutoff = 0))
  expect_length(grep("scorable", ww), 1L)
})

test_that("warn_two_scale_consistency keys on scorable blocks and carries n_scorable", {
  # Direct helper coverage: the trigger is purely structural (block sizes), so it
  # runs on hand-built block lists with NO responses -- localizing a count mutant
  # to the helper and pinning the n_scorable payload. Index values are irrelevant;
  # only the per-block lengths drive the trigger.
  cond <- expect_warning(
    warn_two_scale_consistency(list(A = 1:2, B = 3:4)),     # 2 scorable -> warn
    class = "cier_warning_two_scale_consistency"
  )
  expect_identical(cier_condition_data(cond)$n_scorable, 2L)
  expect_no_warning(                                        # 3 scorable -> silent
    warn_two_scale_consistency(list(A = 1:2, B = 3:4, C = 5:6)),
    class = "cier_warning_two_scale_consistency"
  )
  expect_no_warning(                                        # 1 scorable -> silent
    warn_two_scale_consistency(list(A = 1L, B = 2:3)),
    class = "cier_warning_two_scale_consistency"
  )
})

# ---- Input validation (typed) ----------------------------------------------

test_that("bad scale structure / payload / fpr / cutoff is a typed input error", {
  # Shared with cier_personal_reliability via helper-split-half.R: identical scale /
  # items / payload validation and the check_percentile_overrides contract (an even-odd
  # value is a correlation in [-1, 1]).
  expect_split_half_input_rejected(cier_even_odd)
})

# ---- Cutoff overrides -------------------------------------------------------

test_that("an absolute cutoff overrides the percentile and flags via the upper tail", {
  # Two respondents: consistent (-1) and inverse (+1). cutoff 0 flags only +1.
  it <- blocked_items(2L, 2L, reverse_keyed = FALSE)
  x <- rbind(c(1, 1, 5, 5), c(1, 5, 5, 1))
  storage.mode(x) <- "double"
  # 2-scale fixture: the two-scale warning fires alongside the cutoff override under
  # test, so it is muffled here (covered by the degeneracy block above).
  out <- suppressWarnings(cier_even_odd(x, it, cutoff = 0))
  expect_identical(out$cutoff, 0)
  expect_identical(out$flagged, c(FALSE, TRUE))
})

test_that("a respondent exactly at the cutoff is flagged (>= ties, not >)", {
  it <- blocked_items(5L, 4L, reverse_keyed = FALSE)
  x <- rand_matrix(30L, 20L, 7L)
  v <- suppressWarnings(cier_even_odd(x, it))$value
  k <- which.max(v)
  out <- cier_even_odd(x, it, cutoff = v[[k]])
  expect_true(out$flagged[[k]])
})

# ---- Cross-package parity: vs careless (no-reverse-key data, 1e-12) ---------

test_that("cier_even_odd matches careless::evenodd on no-reverse data (1e-12)", {
  # careless::evenodd forms half-means with base mean(); cier's split-half kernel
  # uses masked-sum sum / count. On this integer Likert fixture they round
  # identically (observed difference 0), but the contract is 1e-12, not bytewise:
  # masked-sum can differ from base mean()'s long-double accumulation by <= 1 ulp on
  # adversarial inputs (mirrors the psychsyn parity).
  skip_if_not_installed("careless")
  raw <- careless::careless_dataset
  responses <- unname(as.matrix(raw))
  storage.mode(responses) <- "double"
  # 10 contiguous scales of 5 items; no reverse-keyed items (careless::evenodd
  # does not reverse-key), so max is not even required here.
  it <- data.frame(scale = rep(paste0("s", LETTERS[1:10]), each = 5L),
                   reverse_keyed = FALSE, stringsAsFactors = FALSE)
  ours <- suppressWarnings(cier_even_odd(responses, it))$value
  theirs <- suppressWarnings(careless::evenodd(raw, factors = rep(5L, 10L)))
  expect_equal(ours, as.numeric(theirs), tolerance = 1e-12)
})

# ---- print snapshot (locked; shared upper-direction format) -----------------

test_that("print renders the locked cli summary (upper direction)", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
    # 24 rows keep the percentile cutoff resolvable (>= 20 scored), so the snapshot
    # locks the resolved-cutoff summary format; suppressWarnings keeps the saturation
    # note (even-odd's +1 point mass) out of the captured print.
    out <- suppressWarnings(cier_even_odd(rand_matrix(24L, 12L, 11L), it))
    expect_snapshot(print(out))
  })
})

test_that("print reports abstaining respondents on their own line", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
    x <- rbind(rand_matrix(24L, 12L, 11L), rep(NA_real_, 12L))
    out <- suppressWarnings(cier_even_odd(x, it))
    expect_snapshot(print(out))
  })
})
