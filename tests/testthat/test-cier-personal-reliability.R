# Tests for cier_personal_reliability() -- Personal Reliability (PR; Jackson
# 1976) and its resampled variant (RPR; Goldammer et al. 2024).
#
# Trust model: two INDEPENDENT oracles re-derive the statistic and never call
# the production kernel. `ref_pr` pins the deterministic first/second-half split
# at 1e-12; `ref_rpr` pins the seeded resampled mean at 1e-10 (it additionally
# coordinates its random-draw order with production so a fixed seed matches --
# a deliberate reproducibility constraint, NOT a tautology: the statistic is
# re-derived from scratch). There is no cross-package partner (careless / psych
# / PerFit / mokken implement no PR/RPR), so the paper oracles ARE the parity.
# Analytic hand fixtures pin exact -1 / +1 values sharing no kernel machinery.
# The mutant block targets each mutant in dev/restart/index-specs.md card 6:
# even/odd split instead of first/second-half, missing Spearman-Brown, not
# negating, RPR aggregation level (per-iteration value vs per-iteration r),
# seed not honoured, and a draw desync on one-item scales.

source(test_path("..", "reference", "ref-pr-jackson-1976.R"))
source(test_path("..", "reference", "ref-rpr-goldammer-2024.R"))

# Scale-blocked `items`: `n_scales` scales of `per_scale` items each.
blocked_items <- function(n_scales = 3L, per_scale = 4L, max = 5L,
                          reverse_keyed = NULL) {
  scale <- rep(LETTERS[seq_len(n_scales)], each = per_scale)
  if (is.null(reverse_keyed)) {
    reverse_keyed <- rep(c(FALSE, TRUE), length.out = length(scale))
  }
  data.frame(scale = scale, reverse_keyed = reverse_keyed,
             max = max, stringsAsFactors = FALSE)
}

# Build the scale blocks INDEPENDENTLY of production (the oracles need them but
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

# Independently reverse-score (min + max) - x on reverse items; the fixtures
# here are 1-based, so the reflection is (1 + max) - x.
prescore <- function(x, items) {
  rk <- items$reverse_keyed
  x[, rk] <- (items$max[rk] + 1L) - x[, rk]
  x
}

# ---- Schema -----------------------------------------------------------------

test_that("cier_personal_reliability returns the pinned list-based schema", {
  # WP3: small/saturated fixtures trip the percentile-cutoff degeneracy guard
  # (D2/D7); these value/oracle tests assert the score, not the flag, so the
  # (correct) warning is muffled.
  it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
  out <- suppressWarnings(
    cier_personal_reliability(rand_matrix(20L, 12L, 1L), it, resample = FALSE)
  )
  expect_s3_class(out, "cier_index")
  expect_type(out, "list")
  expect_identical(names(out),
                   c("value", "flagged", "method", "cutoff", "direction"))
  expect_type(out$value, "double")
  expect_type(out$flagged, "logical")
  expect_identical(length(out$value), 20L)
  expect_identical(length(out$flagged), 20L)
  expect_identical(out$method, "cier_personal_reliability")
  expect_identical(out$direction, "upper")
  expect_type(out$cutoff, "double")
})

test_that("the schema is identical for the RPR (default) variant", {
  it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
  out <- suppressWarnings(
    cier_personal_reliability(rand_matrix(20L, 12L, 1L), it, seed = 1L)
  )
  expect_s3_class(out, "cier_index")
  expect_identical(out$method, "cier_personal_reliability")
  expect_identical(out$direction, "upper")
  expect_identical(length(out$value), 20L)
})

test_that("as.data.frame.cier_index returns the tidy per-respondent frame", {
  it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
  out <- suppressWarnings(
    cier_personal_reliability(rand_matrix(20L, 12L, 2L), it, resample = FALSE)
  )
  df <- as.data.frame(out)
  expect_s3_class(df, "data.frame")
  expect_identical(names(df), c("value", "flagged"))
  expect_identical(nrow(df), 20L)
})

# ---- PR oracle parity (1e-12) -----------------------------------------------

test_that("PR $value equals the ref_pr oracle on a complete matrix", {
  it <- blocked_items(5L, 5L, reverse_keyed = FALSE)
  x <- rand_matrix(30L, 25L, 2026L)
  expect_equal(suppressWarnings(cier_personal_reliability(x, it,
                                                          resample = FALSE))$value,
               ref_pr(x, blocks_from_scale(it)), tolerance = 1e-12)
})

test_that("PR $value equals the ref_pr oracle when rows carry NAs", {
  it <- blocked_items(5L, 5L, reverse_keyed = FALSE)
  x <- rand_matrix(30L, 25L, 99L)
  x[3L, c(1L, 2L)] <- NA       # thins one scale's first half (still a pair)
  x[10L, 1L:23L] <- NA         # only the last scale keeps a pair -> row NA
  expect_equal(suppressWarnings(cier_personal_reliability(x, it,
                                                          resample = FALSE))$value,
               ref_pr(x, blocks_from_scale(it)), tolerance = 1e-12)
})

# ---- RPR oracle parity (1e-10, seeded) --------------------------------------

test_that("RPR $value equals the ref_rpr oracle at a fixed seed", {
  it <- blocked_items(5L, 5L, reverse_keyed = FALSE)
  x <- rand_matrix(30L, 25L, 2026L)
  expect_equal(cier_personal_reliability(x, it, n_resamples = 25L, seed = 42L)$value,
               ref_rpr(x, blocks_from_scale(it), 25L, 42L), tolerance = 1e-10)
})

test_that("RPR $value equals the ref_rpr oracle with NA-thinned rows", {
  # Oracle parity under within-row NAs: row 3 keeps a finite pair on every split,
  # row 10 (only one scale carries data) abstains on every iteration and stays NA.
  # (The mixed finite/NA per-iteration averaging is pinned by the test below.)
  it <- blocked_items(5L, 5L, reverse_keyed = FALSE)
  x <- rand_matrix(30L, 25L, 99L)
  x[3L, c(1L, 2L)] <- NA
  x[10L, 1L:23L] <- NA
  expect_equal(cier_personal_reliability(x, it, n_resamples = 25L, seed = 7L)$value,
               ref_rpr(x, blocks_from_scale(it), 25L, 7L), tolerance = 1e-10)
})

test_that("RPR averages only the finite iterations when a row abstains on some", {
  # A row with one full scale (A) and one 2-of-4-present scale (B) abstains on the
  # random splits that put both present B items in the same half (B loses a half
  # -> only one finite pair -> iteration NA), and is finite otherwise. RPR must
  # average ONLY the finite iterations (na.rm); a no-na.rm mean would return NA.
  it <- data.frame(scale = rep(c("A", "B"), each = 4L),
                   reverse_keyed = FALSE, max = 5L)
  x <- rand_matrix(12L, 8L, 77L)
  x[1L, c(7L, 8L)] <- NA          # scale B (cols 5:8) keeps only items 5, 6
  blocks <- blocks_from_scale(it)
  # Independently materialise the per-iteration values for the target row, in the
  # oracle's draw order, to PROVE the fixture mixes finite and NA iterations.
  per_iter <- withr::with_seed(123L, {
    m <- matrix(NA_real_, nrow = 25L, ncol = nrow(x))
    for (b in seq_len(25L)) m[b, ] <- ref_rpr_iter(x, blocks)
    m
  })
  col1 <- per_iter[, 1L]
  expect_true(anyNA(col1) && any(!is.na(col1)))        # genuinely mixed iterations
  out <- suppressWarnings(
    cier_personal_reliability(x, it, n_resamples = 25L, seed = 123L)
  )
  expect_false(is.na(out$value[[1L]]))                 # finite iterations averaged
  expect_equal(out$value[[1L]], mean(col1, na.rm = TRUE), tolerance = 1e-10)
  expect_equal(out$value, ref_rpr(x, blocks, 25L, 123L), tolerance = 1e-10)
})

test_that("RPR with n_resamples = 1 equals a single seeded iteration", {
  it <- blocked_items(4L, 4L, reverse_keyed = FALSE)
  x <- rand_matrix(15L, 16L, 5L)
  expect_equal(suppressWarnings(cier_personal_reliability(x, it, n_resamples = 1L,
                                                          seed = 3L))$value,
               ref_rpr(x, blocks_from_scale(it), 1L, 3L), tolerance = 1e-10)
})

test_that("the default n_resamples is 25 (anchored to the oracle)", {
  # No n_resamples argument -> the documented default of 25 must be used. A
  # mutant defaulting to any other count (e.g. 10 or 200) diverges from
  # ref_rpr(..., 25L, seed) on the same seed.
  it <- blocked_items(5L, 5L, reverse_keyed = FALSE)
  x <- rand_matrix(30L, 25L, 2026L)
  expect_equal(cier_personal_reliability(x, it, seed = 42L)$value,
               ref_rpr(x, blocks_from_scale(it), 25L, 42L), tolerance = 1e-10)
})

# ---- Analytic fixtures (hand, 1e-12): SB + negation -------------------------

test_that("a perfectly consistent respondent scores -1 (PR)", {
  # Two scales x two items; within-scale identical, varying across scales.
  # No reverse, so the value isolates the Spearman-Brown correlation math.
  it <- blocked_items(2L, 2L, reverse_keyed = FALSE)
  x <- matrix(c(1, 1, 5, 5), nrow = 1L)
  storage.mode(x) <- "double"
  expect_equal(suppressWarnings(cier_personal_reliability(x, it,
                                                          resample = FALSE))$value,
               -1, tolerance = 1e-12)
})

test_that("a perfectly inversely consistent respondent scores +1 (PR)", {
  it <- blocked_items(2L, 2L, reverse_keyed = FALSE)
  x <- matrix(c(1, 5, 5, 1), nrow = 1L)
  storage.mode(x) <- "double"
  expect_equal(suppressWarnings(cier_personal_reliability(x, it,
                                                          resample = FALSE))$value,
               1, tolerance = 1e-12)
})

# ---- Split function + first/second-half vs even/odd ------------------------

test_that("split_half_indices halves by first ceil(k/2) / remaining", {
  # A first/second-half split, NOT the even/odd partition. An even/odd mutant
  # would map split_half_indices(4) to list(c(2,4), c(1,3)).
  expect_identical(split_half_indices(2L), list(first_idx = 1L, second_idx = 2L))
  expect_identical(split_half_indices(3L),
                   list(first_idx = c(1L, 2L), second_idx = 3L))
  expect_identical(split_half_indices(4L),
                   list(first_idx = c(1L, 2L), second_idx = c(3L, 4L)))
  expect_identical(split_half_indices(5L),
                   list(first_idx = c(1L, 2L, 3L), second_idx = c(4L, 5L)))
})

test_that("PR uses the first/second-half split, NOT even/odd", {
  # On 4-item scales the two splits give different half-means, so a PR that
  # secretly used the even/odd partition would equal cier_even_odd. Pins PR
  # against ref_pr AND shows it diverges from even-odd.
  it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
  x <- rand_matrix(20L, 12L, 21L)
  pr <- suppressWarnings(cier_personal_reliability(x, it, resample = FALSE))$value
  eo <- suppressWarnings(cier_even_odd(x, it))$value
  expect_equal(pr, ref_pr(x, blocks_from_scale(it)), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(pr, eo)))
})

# ---- RPR aggregation level (per-iteration value, not per-iteration r) -------

test_that("RPR averages the per-iteration -SB values, not the raw correlations", {
  # ref_rpr applies Spearman-Brown + negation PER ITERATION and then averages;
  # a mutant that averaged the across-scale correlations first and SB-corrected
  # once would diverge (SB is non-linear). The seeded oracle parity pins the
  # value-level aggregation; here we additionally bound the result to [-1, 1].
  it <- blocked_items(5L, 5L, reverse_keyed = FALSE)
  x <- rand_matrix(40L, 25L, 13L)
  v <- suppressWarnings(
    cier_personal_reliability(x, it, n_resamples = 25L, seed = 11L)
  )$value
  expect_equal(v, ref_rpr(x, blocks_from_scale(it), 25L, 11L), tolerance = 1e-10)
  finite <- v[is.finite(v)]
  expect_true(all(finite >= -1 - 1e-9 & finite <= 1 + 1e-9))
})

# ---- Seed handling ----------------------------------------------------------

test_that("the same seed yields identical RPR values; a different seed differs", {
  it <- blocked_items(4L, 4L, reverse_keyed = FALSE)
  x <- rand_matrix(25L, 16L, 8L)
  a <- cier_personal_reliability(x, it, seed = 42L)$value
  b <- cier_personal_reliability(x, it, seed = 42L)$value
  cc <- cier_personal_reliability(x, it, seed = 7L)$value
  expect_identical(a, b)
  expect_false(isTRUE(all.equal(a, cc)))
})

test_that("seed = NULL draws from the ambient RNG stream (no internal reseed)", {
  # A NULL seed must use the session RNG: reproducible UNDER an outer with_seed,
  # and dependent on it. A hardcoded internal seed would make r1 == r2; a
  # set.seed(NULL) entropy reseed would make r1 != r1b.
  it <- blocked_items(4L, 4L, reverse_keyed = FALSE)
  x <- rand_matrix(25L, 16L, 8L)
  r1  <- withr::with_seed(1L, cier_personal_reliability(x, it)$value)
  r1b <- withr::with_seed(1L, cier_personal_reliability(x, it)$value)
  r2  <- withr::with_seed(2L, cier_personal_reliability(x, it)$value)
  expect_equal(r1, r1b)
  expect_false(isTRUE(all.equal(r1, r2)))
})

test_that("a seeded RPR call preserves the caller's RNG state", {
  # The seed is applied locally: a seeded call saves and restores the global
  # .Random.seed, so it does not advance/clobber the caller's random stream.
  it <- blocked_items(4L, 4L, reverse_keyed = FALSE)
  x <- rand_matrix(20L, 16L, 8L)
  withr::with_seed(123L, {
    before <- get(".Random.seed", envir = globalenv())
    cier_personal_reliability(x, it, seed = 1L)
    expect_identical(get(".Random.seed", envir = globalenv()), before)
  })
})

# ---- One-item-scale draw synchronisation (precompute-in-callers contract) ---

test_that("RPR draws NO permutation for a one-item scale (seed stays in sync)", {
  # A one-item scale is skipped without consuming a random draw. If production
  # drew a permutation for it, the whole RNG stream would desync from the oracle
  # and every value would diverge. The mixed battery (a 1-item scale plus three
  # >=2-item scales) makes the bytewise ref_rpr parity pin that skip.
  it <- data.frame(scale = c("A", "B", "B", "C", "C", "C", "D", "D"),
                   reverse_keyed = FALSE, max = 5L)
  x <- rand_matrix(20L, 8L, 314L)
  expect_equal(suppressWarnings(cier_personal_reliability(x, it, n_resamples = 25L,
                                                          seed = 5L))$value,
               ref_rpr(x, blocks_from_scale(it), 25L, 5L), tolerance = 1e-10)
})

# ---- Default variant --------------------------------------------------------

test_that("the default variant is RPR, not deterministic PR", {
  it <- blocked_items(4L, 4L, reverse_keyed = FALSE)
  x <- rand_matrix(30L, 16L, 4L)
  default_v <- cier_personal_reliability(x, it, seed = 1L)$value          # default
  rpr_v     <- cier_personal_reliability(x, it, resample = TRUE, seed = 1L)$value
  pr_v      <- suppressWarnings(
    cier_personal_reliability(x, it, resample = FALSE)
  )$value
  expect_identical(default_v, rpr_v)                 # default == explicit RPR
  expect_false(isTRUE(all.equal(default_v, pr_v)))   # default != PR
})

test_that("PR ignores seed and n_resamples (deterministic)", {
  # The deterministic PR path must not consume the RNG or vary with the
  # resampling controls -- its value is identical with or without them.
  it <- blocked_items(4L, 4L, reverse_keyed = FALSE)
  x <- rand_matrix(20L, 16L, 4L)
  base <- suppressWarnings(cier_personal_reliability(x, it, resample = FALSE))$value
  with_seed <- suppressWarnings(
    cier_personal_reliability(x, it, resample = FALSE, seed = 99L)
  )$value
  with_resamples <- suppressWarnings(
    cier_personal_reliability(x, it, resample = FALSE, n_resamples = 3L)
  )$value
  expect_identical(base, with_seed)
  expect_identical(base, with_resamples)
})

# ---- Direction (upper, NO-FLIP) --------------------------------------------

test_that("direction is upper: a careless (high) row flags, a consistent one does not", {
  it <- blocked_items(5L, 4L, reverse_keyed = FALSE)
  out <- suppressWarnings(
    cier_personal_reliability(rand_matrix(40L, 20L, 11L), it, resample = FALSE)
  )
  expect_true(out$flagged[[which.max(out$value)]])
  expect_false(out$flagged[[which.min(out$value)]])
  expect_identical(out$flagged, !is.na(out$value) & out$value >= out$cutoff)
})

test_that("default cutoff is the upper-tail 95th percentile (NO-FLIP)", {
  it <- blocked_items(5L, 4L, reverse_keyed = FALSE)
  out <- suppressWarnings(
    cier_personal_reliability(rand_matrix(60L, 20L, 5L), it, resample = FALSE)
  )
  expect_equal(out$cutoff,
               as.numeric(stats::quantile(out$value[is.finite(out$value)],
                                          0.95, names = FALSE, type = 7L)),
               tolerance = 1e-12)
})

test_that("the cutoff is the 95th percentile for the RPR (default) variant too", {
  it <- blocked_items(5L, 4L, reverse_keyed = FALSE)
  out <- cier_personal_reliability(rand_matrix(60L, 20L, 5L), it, seed = 9L)
  expect_equal(out$cutoff,
               as.numeric(stats::quantile(out$value[is.finite(out$value)],
                                          0.95, names = FALSE, type = 7L)),
               tolerance = 1e-12)
})

test_that("the fpr argument moves the percentile target", {
  it <- blocked_items(5L, 4L, reverse_keyed = FALSE)
  out <- suppressWarnings(
    cier_personal_reliability(rand_matrix(60L, 20L, 5L), it, resample = FALSE, fpr = 0.10)
  )
  expect_equal(out$cutoff,
               as.numeric(stats::quantile(out$value[is.finite(out$value)],
                                          0.90, names = FALSE, type = 7L)),
               tolerance = 1e-12)
})

# ---- Reverse-keying (kills keying / double-reflection mutants) --------------

test_that("honouring reverse_keyed equals independently pre-scoring (PR)", {
  it_rev <- blocked_items(4L, 4L)          # alternating reverse
  x <- rand_matrix(25L, 16L, 101L)
  it_fwd <- it_rev
  it_fwd$reverse_keyed <- FALSE
  honoured  <- suppressWarnings(
    cier_personal_reliability(x, it_rev, resample = FALSE)
  )$value
  prescored <- suppressWarnings(
    cier_personal_reliability(prescore(x, it_rev), it_fwd, resample = FALSE)
  )$value
  expect_equal(honoured, prescored, tolerance = 1e-12)
})

test_that("honoured PR reverse-keying also equals the oracle on pre-scored input", {
  it_rev <- blocked_items(4L, 4L)
  x <- rand_matrix(25L, 16L, 7L)
  expect_equal(suppressWarnings(cier_personal_reliability(x, it_rev,
                                                          resample = FALSE))$value,
               ref_pr(prescore(x, it_rev), blocks_from_scale(it_rev)),
               tolerance = 1e-12)
})

test_that("honouring reverse_keyed equals pre-scoring for RPR (same seed)", {
  # Reverse-keying is applied BEFORE the kernel, so the seeded random draws act
  # on the same reflected data either way -> identical values.
  it_rev <- blocked_items(4L, 4L)
  x <- rand_matrix(25L, 16L, 55L)
  it_fwd <- it_rev
  it_fwd$reverse_keyed <- FALSE
  honoured  <- cier_personal_reliability(x, it_rev, seed = 2L)$value
  prescored <- cier_personal_reliability(prescore(x, it_rev), it_fwd,
                                         seed = 2L)$value
  expect_equal(honoured, prescored, tolerance = 1e-12)
})

test_that("honoured RPR reverse-keying also equals the oracle on pre-scored input", {
  # The PR analogue above pins keying against the independent oracle; mirror it
  # for RPR so the resampled path has its own non-self-referential anchor (the
  # self-consistency check above routes BOTH sides through the production
  # kernel and could not catch a keying bug symmetric across them).
  it_rev <- blocked_items(4L, 4L)
  x <- rand_matrix(25L, 16L, 55L)
  expect_equal(cier_personal_reliability(x, it_rev, seed = 2L)$value,
               ref_rpr(prescore(x, it_rev), blocks_from_scale(it_rev),
                       25L, 2L),
               tolerance = 1e-10)
})

# ---- Edge cases -------------------------------------------------------------

test_that("a single-item scale is skipped; remaining scales still score (PR)", {
  it <- data.frame(scale = c("A", "B", "B", "C", "C"),
                   reverse_keyed = FALSE, max = 5L)
  x <- matrix(c(3, 1, 5, 2, 4,
                4, 2, 4, 1, 5), nrow = 2L, byrow = TRUE)
  storage.mode(x) <- "double"
  out <- suppressWarnings(cier_personal_reliability(x, it, resample = FALSE))
  expect_false(any(is.na(out$value)))     # scales B and C give two finite pairs
  expect_equal(out$value, ref_pr(x, blocks_from_scale(it)), tolerance = 1e-12)
})

test_that("a constant (straightliner) row abstains (zero variance -> NA, PR)", {
  it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
  x <- rand_matrix(6L, 12L, 4L)
  x[2L, ] <- 3                            # constant -> identical half means
  out <- suppressWarnings(cier_personal_reliability(x, it, resample = FALSE))
  expect_true(is.na(out$value[[2L]]))
  expect_false(is.na(out$value[[1L]]))
})

test_that("an all-NA row abstains and keeps the remaining rows aligned (PR)", {
  it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
  x <- rand_matrix(10L, 12L, 4L)
  x[5L, ] <- NA
  out <- suppressWarnings(cier_personal_reliability(x, it, resample = FALSE))
  expect_true(is.na(out$value[[5L]]))
  expect_true(is.na(out$flagged[[5L]]))
  expect_false(is.na(out$value[[1L]]))
  expect_false(is.na(out$value[[10L]]))
  expect_equal(out$value, ref_pr(x, blocks_from_scale(it)), tolerance = 1e-12)
})

test_that("when every respondent abstains the cutoff warns and flags nobody", {
  it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
  x <- matrix(3, nrow = 5L, ncol = 12L)  # all straightliners -> all NA
  storage.mode(x) <- "double"
  expect_warning(out <- cier_personal_reliability(x, it, resample = FALSE),
                 class = "cier_warning_insufficient_items")
  expect_true(all(is.na(out$value)))
  expect_true(is.na(out$cutoff))
  expect_true(all(is.na(out$flagged)))
})

# ---- Input validation (typed) ----------------------------------------------

test_that("fewer than two distinct scales is a typed input error", {
  it <- data.frame(scale = rep("A", 4L), reverse_keyed = FALSE, max = 5L)
  x <- matrix(c(1, 2, 3, 4, 5, 4, 3, 2), nrow = 2L, byrow = TRUE)
  storage.mode(x) <- "double"
  expect_error(cier_personal_reliability(x, it, resample = FALSE),
               class = "cier_error_input")
})

test_that("a missing scale column is a typed input error", {
  it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
  it$scale <- NULL
  expect_error(cier_personal_reliability(rand_matrix(5L, 12L, 1L), it),
               class = "cier_error_input")
})

test_that("an items frame not aligned to the columns is a typed input error", {
  it <- blocked_items(3L, 4L, reverse_keyed = FALSE)   # 12 items
  expect_error(cier_personal_reliability(rand_matrix(5L, 10L, 1L), it),
               class = "cier_error_input")
})

test_that("a non-matrix / non-numeric payload is a typed input error", {
  it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
  expect_error(cier_personal_reliability(1:10, it), class = "cier_error_input")
  expect_error(cier_personal_reliability(matrix(letters[1:12], nrow = 1L), it),
               class = "cier_error_input")
})

# ---- Cutoff overrides -------------------------------------------------------

test_that("an absolute cutoff overrides the percentile and flags via the upper tail", {
  it <- blocked_items(2L, 2L, reverse_keyed = FALSE)
  x <- rbind(c(1, 1, 5, 5), c(1, 5, 5, 1))   # consistent (-1) and inverse (+1)
  storage.mode(x) <- "double"
  out <- cier_personal_reliability(x, it, resample = FALSE, cutoff = 0)
  expect_identical(out$cutoff, 0)
  expect_identical(out$flagged, c(FALSE, TRUE))
})

test_that("a respondent exactly at the cutoff is flagged (>= ties, not >)", {
  it <- blocked_items(5L, 4L, reverse_keyed = FALSE)
  x <- rand_matrix(30L, 20L, 7L)
  v <- suppressWarnings(cier_personal_reliability(x, it, resample = FALSE))$value
  k <- which.max(v)
  out <- cier_personal_reliability(x, it, resample = FALSE, cutoff = v[[k]])
  expect_true(out$flagged[[k]])
})

test_that("invalid fpr values are typed input errors", {
  it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
  x <- rand_matrix(10L, 12L, 1L)
  expect_error(cier_personal_reliability(x, it, resample = FALSE, fpr = 0),
               class = "cier_error_input")
  expect_error(cier_personal_reliability(x, it, resample = FALSE, fpr = 1),
               class = "cier_error_input")
  expect_error(cier_personal_reliability(x, it, resample = FALSE, fpr = -0.1),
               class = "cier_error_input")
  expect_error(cier_personal_reliability(x, it, resample = FALSE,
                                         fpr = NA_real_),
               class = "cier_error_input")
  expect_error(cier_personal_reliability(x, it, resample = FALSE,
                                         fpr = c(0.05, 0.1)),
               class = "cier_error_input")
})

test_that("an invalid absolute cutoff is a typed input error", {
  it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
  x <- rand_matrix(10L, 12L, 1L)
  expect_error(cier_personal_reliability(x, it, resample = FALSE,
                                         cutoff = c(0.5, 1)),
               class = "cier_error_input")
  expect_error(cier_personal_reliability(x, it, resample = FALSE,
                                         cutoff = NA_real_),
               class = "cier_error_input")
  # A PR/RPR value is a (negated, SB-corrected) correlation in [-1, 1]; a
  # threshold outside that range flags everyone or no one.
  expect_error(cier_personal_reliability(x, it, resample = FALSE, cutoff = -1.5),
               class = "cier_error_input")
  expect_error(cier_personal_reliability(x, it, resample = FALSE, cutoff = 1.5),
               class = "cier_error_input")
})

test_that("supplying both fpr and cutoff is a typed input error", {
  it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
  expect_error(cier_personal_reliability(rand_matrix(10L, 12L, 1L), it,
                                         resample = FALSE,
                                         fpr = 0.1, cutoff = 0.5),
               class = "cier_error_input")
})

# ---- RPR argument validation ------------------------------------------------

test_that("an invalid n_resamples is a typed input error", {
  it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
  x <- rand_matrix(10L, 12L, 1L)
  expect_error(cier_personal_reliability(x, it, n_resamples = 0L),
               class = "cier_error_input")
  expect_error(cier_personal_reliability(x, it, n_resamples = -5L),
               class = "cier_error_input")
  expect_error(cier_personal_reliability(x, it, n_resamples = 2.5),
               class = "cier_error_input")
  expect_error(cier_personal_reliability(x, it, n_resamples = c(10L, 20L)),
               class = "cier_error_input")
  expect_error(cier_personal_reliability(x, it, n_resamples = NA_integer_),
               class = "cier_error_input")
})

test_that("an invalid seed is a typed input error", {
  it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
  x <- rand_matrix(10L, 12L, 1L)
  expect_error(cier_personal_reliability(x, it, seed = "x"),
               class = "cier_error_input")
  expect_error(cier_personal_reliability(x, it, seed = c(1L, 2L)),
               class = "cier_error_input")
  expect_error(cier_personal_reliability(x, it, seed = NA_integer_),
               class = "cier_error_input")
  # A non-integer seed would be silently truncated by set.seed(); reject it.
  expect_error(cier_personal_reliability(x, it, seed = 1.9),
               class = "cier_error_input")
})

test_that("an invalid resample flag is a typed input error", {
  it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
  x <- rand_matrix(10L, 12L, 1L)
  expect_error(cier_personal_reliability(x, it, resample = NA),
               class = "cier_error_input")
  expect_error(cier_personal_reliability(x, it, resample = "yes"),
               class = "cier_error_input")
  expect_error(cier_personal_reliability(x, it, resample = c(TRUE, FALSE)),
               class = "cier_error_input")
})

test_that("n_resamples and seed are validated even on the PR path", {
  # The wrapper validates the resampling controls up front (early fail),
  # regardless of `resample`, so a bad value errors even when it is unused.
  it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
  x <- rand_matrix(10L, 12L, 1L)
  expect_error(cier_personal_reliability(x, it, resample = FALSE,
                                         n_resamples = 0L),
               class = "cier_error_input")
  expect_error(cier_personal_reliability(x, it, resample = FALSE, seed = "x"),
               class = "cier_error_input")
})

# ---- print snapshot (locked; shared upper-direction format) -----------------

test_that("print renders the locked cli summary (deterministic PR)", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
    # WP3: 24 rows keep the percentile cutoff resolvable (>= 20 scored, D2), so
    # the snapshot locks the resolved-cutoff summary format; suppressWarnings keeps
    # the saturation note out of the captured print.
    x <- rand_matrix(24L, 12L, 6L)
    out <- suppressWarnings(cier_personal_reliability(x, it, resample = FALSE))
    expect_snapshot(print(out))
  })
})

test_that("print reports abstaining respondents on their own line (PR)", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
    x <- rbind(rand_matrix(24L, 12L, 6L), rep(NA_real_, 12L))
    out <- suppressWarnings(cier_personal_reliability(x, it, resample = FALSE))
    expect_snapshot(print(out))
  })
})

test_that("a complete straightliner abstains as '(no score)', not '(no responses)'", {
  # Regression guard for the reason-neutral abstention noun. A fully-answered
  # straightliner has zero within-person variance, so the split-half correlation
  # is undefined and PR abstains -- with NO missing data. The respondent answered
  # every item, so the old method-keyed '(no responses)' was a false reason; the
  # reason-neutral '(no score)' is the honest wording the print must use.
  it <- blocked_items(3L, 4L, reverse_keyed = FALSE)
  x  <- rbind(rand_matrix(11L, 12L, 6L), rep(3, 12L))  # row 12: complete straightliner
  out <- suppressWarnings(cier_personal_reliability(x, it, resample = FALSE))
  expect_true(is.na(out$value[12L]))                   # abstains (unscalable)
  expect_false(anyNA(x[12L, ]))                         # yet has no missing data
  printed <- withr::with_options(
    list(cli.width = 80, cli.unicode = FALSE),
    capture.output(print(out))
  )
  expect_true(any(grepl("(no score)", printed, fixed = TRUE)))
  expect_false(any(grepl("(no responses)", printed, fixed = TRUE)))
})
