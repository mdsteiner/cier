# Tests for cier_psychant() -- the psychometric-antonyms C/IER index, the mirror
# of cier_psychsyn() over the strongly NEGATIVELY correlated item pairs.
#
# Trust model: the independent oracle (ref_psychant) re-derives each respondent's
# correlation across the antonym pairs with an explicit nested-loop pair search
# (r < -|critval|, i > j) + a per-row stacked-pair stats::cor() (reusing the
# synonyms row function ref_psychsyn_row), and never calls the production kernel
# (which discovers pairs from the column-major lower triangle and scores every
# respondent in one vectorised masked-sum pass). The kernel computes the same
# Pearson correlation re-summed in a different order, so it matches the oracle to
# ~1e-13 and is held at 1e-12. A separate orthogonal-contrast hand fixture pins
# the exact pair SET and exact scored values (-1) with no oracle machinery. The
# cross-package check pins parity with careless::psychsyn(anto = TRUE,
# resample_na = FALSE) at 1e-12 (see tests/reference/TOLERANCES.md): the same
# vectorised masked-sum kernel that psychsyn uses, so it inherits psychsyn's
# 1e-12 (no longer bytewise after the kernel was vectorised for a 4-6x speedup);
# careless::psychant() does not surface resample_na, so the deterministic
# comparison calls the underlying psychsyn(anto = TRUE, resample_na = FALSE)
# directly. The property / mutant-killer block targets each mutant named in
# dev/restart/index-specs.md card 8 (synonym threshold instead of antonym, flag
# lower instead of upper, wrong sign branch, reverse-keying before pairing).

source(test_path("..", "reference", "ref-psychsyn-meade-craig-2012.R"))
source(test_path("..", "reference", "ref-psychant-meade-craig-2012.R"))

# A reproducible matrix with genuine antonym pairs: `k` mutually-orthogonal latent
# factors, each driving ONE positively-loading item (cols 1..k) and ONE
# negatively-loading item (cols k+1..2k) with small noise, so each factor's two
# items correlate well below -0.60 (about -0.89 at noise 0.35) and items across
# factors correlate near 0. The antonym pairs are therefore EXACTLY (k+j, j) for
# j = 1..k; the per-respondent within-pair opposition varies with the noise, so
# the scored values span a range below 0 (not a constant -1) -- needed for the
# direction test.
ant_matrix <- function(n = 60L, k = 4L, seed = 42L, noise = 0.35) {
  withr::with_seed(seed, {
    factors <- lapply(seq_len(k), function(j) stats::rnorm(n))
    pos <- vapply(factors, function(f) f + stats::rnorm(n, 0, noise), numeric(n))
    neg <- vapply(factors, function(f) -f + stats::rnorm(n, 0, noise), numeric(n))
    x <- cbind(pos, neg)
  })
  storage.mode(x) <- "double"
  x
}

# Orthogonal-contrast hand fixture (n = 4 respondents, 6 items in 3 antonym
# pairs). Columns 1/3/5 are the mutually orthogonal Hadamard contrasts of order 4
# (every cross-correlation exactly 0); columns 2/4/6 are their exact negations
# (within-pair r = -1). The qualifying ANTONYM pair set is therefore EXACTLY
# {(2,1), (4,3), (6,5)} -- worked by hand, no oracle -- and NO synonym pair
# exists (no off-diagonal r is positive). For each non-constant respondent the
# stacked larger-item values are the exact negatives of the smaller-item values,
# so the antonym correlation is -1; the all-constant straightliner (row 1) abstains
# (NA). This is the exact mirror of the psychsyn hand fixture (which scored +1).
hand_fixture <- function() {
  a  <- c(4, 4, 2, 2)   # contrast (1, 1, -1, -1)
  na <- c(2, 2, 4, 4)   # -a:  r(a, na) = -1
  b  <- c(4, 2, 4, 2)   # contrast (1, -1, 1, -1)
  nb <- c(2, 4, 2, 4)   # -b
  c3 <- c(4, 2, 2, 4)   # contrast (1, -1, -1, 1)
  nc <- c(2, 4, 4, 2)   # -c3
  m <- cbind(a, na, b, nb, c3, nc)
  dimnames(m) <- NULL
  storage.mode(m) <- "double"
  m
}

# A constructed cross-package fixture with planted antonym structure: items
# (k+1)..2k are the negatives of items 1..k plus noise, mapped to a 1..5 Likert.
# careless_dataset carries NO antonym pairs at r < -0.60, so a genuine careless
# parity check needs constructed negative structure (the empty-antonym
# careless_dataset would be a trivial all-NA match).
ant_careless_fixture <- function(n = 80L, half = 7L, seed = 202L) {
  withr::with_seed(seed, {
    z <- matrix(stats::rnorm(n * (2L * half)), nrow = n, ncol = 2L * half)
    z[, (half + 1L):(2L * half)] <-
      -z[, 1L:half] + matrix(stats::rnorm(n * half, sd = 0.4), nrow = n)
    responses <- round(pmin(pmax(z + 3, 1), 5))
  })
  storage.mode(responses) <- "double"
  responses
}

# ---- Schema -----------------------------------------------------------------

test_that("cier_psychant returns a list-based cier_index with the pinned schema", {
  out <- cier_psychant(ant_matrix())
  expect_s3_class(out, "cier_index")
  expect_type(out, "list")
  expect_identical(names(out),
                   c("value", "flagged", "method", "cutoff", "direction"))
  expect_type(out$value, "double")
  expect_type(out$flagged, "logical")
  expect_identical(length(out$value), 60L)
  expect_identical(length(out$flagged), 60L)
  expect_identical(out$method, "cier_psychant")
  expect_identical(out$direction, "upper")
  expect_type(out$cutoff, "double")
})

test_that("as.data.frame.cier_index returns the tidy per-respondent frame", {
  df <- as.data.frame(cier_psychant(ant_matrix(n = 12L)))
  expect_s3_class(df, "data.frame")
  expect_identical(names(df), c("value", "flagged"))
  expect_identical(nrow(df), 12L)
})

# ---- Independent oracle parity (1e-12) --------------------------------------

test_that("cier_psychant$value equals the oracle on a complete matrix", {
  x <- ant_matrix(n = 50L, seed = 2026L)
  expect_equal(cier_psychant(x)$value, ref_psychant(x, 0.60),
               tolerance = 1e-12)
})

test_that("cier_psychant$value equals the oracle when rows carry NAs", {
  x <- ant_matrix(n = 50L, seed = 99L)
  x[3L, c(1L, 5L)] <- NA           # a couple of cells dropped, still scored
  x[10L, 1L:6L] <- NA              # most pairs gone -> too few -> abstains
  x[20L, ] <- NA                   # all-NA respondent -> abstains
  ours <- cier_psychant(x)$value
  expect_equal(ours, ref_psychant(x, 0.60), tolerance = 1e-12)
  expect_true(is.na(ours[[20L]]))
})

test_that("cier_psychant$value equals the hand-computed fixture exactly", {
  # Pair set is the three orthogonally-separated antonym pairs; every non-constant
  # respondent scores exactly -1, the all-constant row abstains (NA).
  expect_equal(cier_psychant(hand_fixture())$value, c(NA, -1, -1, -1),
               tolerance = 1e-12)
})

# ---- Cross-package parity: careless (anto = TRUE, 1e-12) --------------------

test_that("cier_psychant matches careless antonyms on a planted fixture (1e-12)", {
  # GENUINE cross-package parity: careless's antonym scorer is an independent pair
  # search + stacked-pair correlation. careless::psychant() does not expose
  # resample_na, so the deterministic comparison calls the underlying
  # psychsyn(anto = TRUE, resample_na = FALSE) directly (careless's default
  # resample_na = TRUE is a non-deterministic NA fallback cier does not
  # reproduce). The shared masked-sum kernel computes that same correlation as one
  # vectorised pass, agreeing to ~1e-13 -- held at 1e-12, the same tolerance
  # psychsyn carries post-vectorisation. The fixture carries injected NAs (a few
  # cells + one all-NA row) so this also pins NA agreement. See
  # tests/reference/TOLERANCES.md.
  skip_if_not_installed("careless")
  responses <- ant_careless_fixture()
  responses[5L, c(2L, 9L)] <- NA          # a few dropped cells, still scored
  responses[40L, ] <- NA                  # an all-NA respondent -> NA both sides
  ours <- cier_psychant(responses, critical_r = 0.60)$value
  theirs <- careless::psychsyn(responses, critval = -0.60, anto = TRUE,
                               resample_na = FALSE)
  expect_equal(ours, as.numeric(theirs), tolerance = 1e-12)
})

# ---- Pair discovery / mutant-killers ----------------------------------------

test_that("find_item_pairs selects exactly the antonym pairs (no cross, no dupes)", {
  # The orthogonal hand fixture has EXACTLY three qualifying antonym pairs. A
  # mutant that scans the upper triangle too would return six (each pair twice)
  # and break the oracle / careless parity; the count + set pin it here directly.
  pairs <- find_item_pairs(hand_fixture(), 0.60, "ant")
  expect_identical(nrow(pairs), 3L)
  expect_true(all(pairs[, 1L] > pairs[, 2L]))    # larger item index first
  got <- apply(pairs, 1L, function(r) paste(sort(r), collapse = "-"))
  expect_setequal(got, c("1-2", "3-4", "5-6"))
  # The same orthogonal data has no synonym pairs (no off-diagonal r is positive).
  expect_identical(nrow(find_item_pairs(hand_fixture(), 0.60, "syn")), 0L)
})

test_that("the antonym branch selects on the NEGATIVE tail, not the positive one", {
  # Items 1 and 2 are perfectly NEGATIVELY correlated (r = -1); items 3 and 4 are
  # perfectly POSITIVELY correlated (r = +1); the rest orthogonal. The antonym
  # branch must select ONLY (2, 1); the synonym branch must select ONLY (4, 3). A
  # wrong-sign mutant (> r vs < -r swapped) inverts both -- the core psychant
  # mutant.
  a   <- c(4, 4, 2, 2)
  na  <- c(2, 2, 4, 4)             # -a: r(a, na) = -1  (antonym)
  b   <- c(4, 2, 4, 2)
  b2  <- c(4, 2, 4, 2)             # = b: r(b, b2) = +1 (synonym)
  m <- cbind(a, na, b, b2)
  dimnames(m) <- NULL
  storage.mode(m) <- "double"
  ant <- find_item_pairs(m, 0.60, "ant")
  expect_identical(nrow(ant), 1L)
  expect_identical(as.integer(ant[1L, ]), c(2L, 1L))   # larger item index first
  syn <- find_item_pairs(m, 0.60, "syn")
  expect_identical(nrow(syn), 1L)
  expect_identical(as.integer(syn[1L, ]), c(4L, 3L))
})

test_that("the wrapper is matrix-only: no item-metadata channel exists", {
  # Pins keying-insensitivity structurally: there is no items / scale /
  # reverse_keyed / categories argument, so a reverse-keying-before-pairing
  # mutant cannot be wired in. Pairing uses the raw responses (reverse-keying
  # would collapse genuine antonym correlations toward zero).
  fmls <- names(formals(cier_psychant))
  expect_setequal(fmls, c("responses", "critical_r", "fpr", "cutoff"))
  expect_false(any(c("items", "scale", "reverse_keyed",
                     "categories") %in% fmls))
})

test_that("scored antonym values are finite, within [-1, 1], and vary per respondent", {
  v <- cier_psychant(ant_matrix(n = 60L, seed = 7L))$value
  scored <- v[!is.na(v)]
  expect_gt(length(scored), 0L)
  expect_true(all(is.finite(scored)))
  expect_true(all(scored >= -1 & scored <= 1))
  # The score is ONE stacked correlation per respondent, so it varies row to row.
  # A mutant that returns the (constant) whole-sample inter-item r averaged over
  # the qualifying pairs would make every scored value identical.
  expect_gt(length(unique(round(scored, 8L))), 1L)
})

test_that("a straightliner (zero-variance pair side) abstains, not a resampled value", {
  # A constant row has zero variance on each pair side, so the antonym
  # correlation is NA. A resample_na-fallback mutant would impute a non-NA
  # permutation value instead; we require NA.
  x <- ant_matrix(n = 20L, seed = 5L)
  x[2L, ] <- 0                      # midpoint straightliner (zero variance both sides)
  out <- cier_psychant(x)
  expect_true(is.na(out$value[[2L]]))
  expect_false(is.na(out$value[[1L]]))
})

# ---- Direction (upper) ------------------------------------------------------

test_that("direction is upper: a decoupled row flags, a strongly-opposed one does not", {
  # The least within-pair-opposed respondent has the HIGHEST (least negative)
  # antonym correlation -- the most careless -- so it flags at the upper-tail
  # cutoff; the most strongly opposed (most negative) does not. A flag-lower
  # mutant inverts both.
  out <- cier_psychant(ant_matrix(n = 60L, seed = 11L))
  expect_true(out$flagged[[which.max(out$value)]])
  expect_false(out$flagged[[which.min(out$value)]])
  expect_identical(out$flagged, !is.na(out$value) & out$value >= out$cutoff)
})

# ---- Edge cases -------------------------------------------------------------

test_that("critical_r too high finds no pairs: every row abstains and flags nobody", {
  # No inter-item correlation clears -0.99, so there are no antonym pairs, every
  # value is NA, and the percentile cutoff has no finite values to rank: it warns
  # and returns NA, and an NA cutoff flags no one.
  expect_warning(out <- cier_psychant(ant_matrix(n = 30L), critical_r = 0.99),
                 class = "cier_warning_insufficient_items")
  expect_true(all(is.na(out$value)))
  expect_true(is.na(out$cutoff))
  expect_true(all(is.na(out$flagged)))
})

test_that("a respondent with fewer than three complete pairs abstains; rows stay aligned", {
  x <- ant_matrix(n = 12L, seed = 3L)
  x[5L, c(2L, 3L, 4L, 6L, 7L, 8L)] <- NA   # leaves only pair (5,1) complete -> abstains
  out <- cier_psychant(x)
  expect_true(is.na(out$value[[5L]]))
  expect_true(is.na(out$flagged[[5L]]))
  expect_false(is.na(out$value[[1L]]))
  expect_false(is.na(out$value[[12L]]))
})

test_that("the complete-pair boundary is exactly three: three scores, two/one abstain", {
  # ant_matrix(k = 4) has four antonym pairs (5,1) (6,2) (7,3) (8,4). NA-ing one
  # item of one pair leaves EXACTLY three complete pairs -> scored (n_complete = 3);
  # NA-ing one item of two pairs leaves two complete -> abstains (n_complete = 2);
  # NA-ing one item of three pairs leaves one complete -> abstains. This pins the
  # n_complete <= 2 -> NA / >= 3 -> scored knife-edge from both sides; an
  # off-by-one in either direction flips an assertion.
  x <- ant_matrix(n = 30L, seed = 8L)
  x[3L, 8L] <- NA                  # breaks pair (8,4) -> 3 complete -> scored
  x[7L, c(7L, 8L)] <- NA           # breaks (7,3),(8,4) -> 2 complete -> abstains
  x[9L, c(6L, 7L, 8L)] <- NA       # breaks 3 pairs -> 1 complete -> abstains
  out <- cier_psychant(x)
  expect_false(is.na(out$value[[3L]]))
  expect_true(is.na(out$value[[7L]]))
  expect_true(is.na(out$value[[9L]]))
})

test_that("a constant item column drops out of pairing without leaking a warning", {
  # A zero-variance item makes stats::cor() emit a base-R, locale-dependent
  # "standard deviation is zero" warning; pairing_cor() must suppress it (the
  # package signals only typed cli conditions) and simply not pair the constant
  # item. The remaining antonym pairs still score, so no abstention warning fires
  # either -- the call is warning-free.
  x <- ant_matrix(n = 40L, seed = 6L)
  x[, 5L] <- 3                       # a constant (straightlined) item column
  expect_no_warning(out <- cier_psychant(x))
  expect_gt(sum(!is.na(out$value)), 0L)
})

test_that("a non-matrix / non-numeric / non-finite payload is a typed input error", {
  expect_error(cier_psychant(1:10), class = "cier_error_input")
  expect_error(cier_psychant(matrix(letters[1:6], nrow = 2L)),
               class = "cier_error_input")
  bad <- ant_matrix(n = 6L)
  bad[1L, 1L] <- Inf
  expect_error(cier_psychant(bad), class = "cier_error_input")
})

# ---- Cutoff: default, fpr override, NO-FLIP direction -----------------------

test_that("default cutoff is the upper-tail 95th percentile (NO-FLIP)", {
  out <- cier_psychant(ant_matrix(n = 80L, seed = 5L))
  # Upper direction takes the 1 - fpr quantile (NOT fpr): the registry stores the
  # literal tail mass fpr = 0.05 and resolve_cutoff() applies the single flip.
  expect_equal(out$cutoff,
               as.numeric(stats::quantile(out$value[is.finite(out$value)], 0.95,
                                          names = FALSE, type = 7L)),
               tolerance = 1e-12)
})

test_that("the fpr argument moves the percentile target", {
  x <- ant_matrix(n = 80L, seed = 5L)
  out <- cier_psychant(x, fpr = 0.10)
  expect_equal(out$cutoff,
               as.numeric(stats::quantile(out$value[is.finite(out$value)], 0.90,
                                          names = FALSE, type = 7L)),
               tolerance = 1e-12)
})

test_that("an absolute cutoff overrides the percentile and flags via direction (>= ties)", {
  x <- ant_matrix(n = 40L, seed = 7L)
  v <- cier_psychant(x)$value
  k <- which.max(v)
  out <- cier_psychant(x, cutoff = v[[k]])
  expect_identical(out$cutoff, v[[k]])
  expect_true(out$flagged[[k]])     # exactly at the cutoff flags (>=, not >)
})

test_that("invalid critical_r values are typed input errors (positive magnitude only)", {
  x <- ant_matrix(n = 10L)
  expect_error(cier_psychant(x, critical_r = 0), class = "cier_error_input")
  expect_error(cier_psychant(x, critical_r = 1), class = "cier_error_input")
  # The convention is a POSITIVE magnitude in (0, 1); a negative value (the sign
  # careless's critval uses) is rejected -- pins the documented sign convention.
  expect_error(cier_psychant(x, critical_r = -0.5), class = "cier_error_input")
  expect_error(cier_psychant(x, critical_r = -0.60), class = "cier_error_input")
  expect_error(cier_psychant(x, critical_r = c(0.3, 0.4)),
               class = "cier_error_input")
  expect_error(cier_psychant(x, critical_r = NA_real_),
               class = "cier_error_input")
})

test_that("invalid fpr / cutoff values are typed input errors", {
  x <- ant_matrix(n = 10L)
  expect_error(cier_psychant(x, fpr = 0), class = "cier_error_input")
  expect_error(cier_psychant(x, fpr = 1), class = "cier_error_input")
  expect_error(cier_psychant(x, fpr = c(0.05, 0.1)), class = "cier_error_input")
  # An antonym value is a correlation in [-1, 1]; a threshold outside that range
  # is degenerate (flags everyone or no one), so it is rejected.
  expect_error(cier_psychant(x, cutoff = -1.5), class = "cier_error_input")
  expect_error(cier_psychant(x, cutoff = 1.5), class = "cier_error_input")
  expect_error(cier_psychant(x, cutoff = NA_real_), class = "cier_error_input")
})

test_that("supplying both fpr and cutoff is a typed input error", {
  expect_error(cier_psychant(ant_matrix(n = 10L), fpr = 0.1, cutoff = 0.5),
               class = "cier_error_input")
})

# ---- print snapshot (locked, design-first; direction = upper) ---------------

test_that("print renders the locked cli summary (upper direction)", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    expect_snapshot(print(cier_psychant(ant_matrix(n = 30L, seed = 11L))))
  })
})

test_that("print reports abstaining respondents as '(no score)'", {
  # A psychant respondent can abstain despite answering (too few complete pairs),
  # so the shared abstention line reads '(no score)', not '(no responses)'.
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    x <- ant_matrix(n = 29L, seed = 11L)
    x <- rbind(x, rep(NA_real_, ncol(x)))   # one abstaining respondent
    expect_snapshot(print(cier_psychant(x)))
  })
})

test_that("the abstaining-row noun is '(no score)' and not '(no responses)'", {
  # Pinned independently of the generated snapshot above: the snapshot records
  # whatever the implementation emits, so a wrapper that forgets to extend
  # abstention_noun() for cier_psychant would bless the wrong '(no responses)' on
  # first generation. This asserts the intended noun directly -- psychant can
  # abstain on a fully-answered row (too few antonym pairs), so '(no score)' is
  # the honest wording.
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    x <- ant_matrix(n = 29L, seed = 11L)
    x <- rbind(x, rep(NA_real_, ncol(x)))
    out <- capture.output(print(cier_psychant(x)))
  })
  expect_true(any(grepl("(no score)", out, fixed = TRUE)))
  expect_false(any(grepl("(no responses)", out, fixed = TRUE)))
})
