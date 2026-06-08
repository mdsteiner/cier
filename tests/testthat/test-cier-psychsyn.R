# Tests for cier_psychsyn() -- the psychometric-synonyms C/IER index.
#
# Trust model: the independent oracle (ref_psychsyn) re-derives each
# respondent's correlation across the strongly-positively-correlated item pairs
# with an explicit nested-loop pair search + a per-row stacked-pair stats::cor(),
# and never calls the production kernel (which discovers pairs from the
# column-major lower triangle and scores every respondent in one vectorised
# masked-sum pass). The kernel computes the same Pearson correlation re-summed in
# a different order, so it matches the oracle to ~1e-13 and is held at 1e-12. A
# separate orthogonal-contrast hand fixture pins the exact pair SET and exact
# +/-1 values with no oracle machinery. The cross-package check pins parity with
# careless::psychsyn(resample_na = FALSE) at 1e-12 (see
# tests/reference/TOLERANCES.md): careless::syn_for_one scores a per-row cor()
# over the same stacked pair vectors, which the vectorised masked-sum kernel
# reproduces to ~1e-13 -- no longer bytewise, the kernel having been vectorised
# for a 4-6x speedup (see the ADR entry); the careless default resample_na = TRUE
# is a non-deterministic NA fallback deliberately not reproduced. The property /
# mutant-killer block targets each mutant named in dev/restart/index-specs.md
# card 7 (antonym threshold instead of synonym, upper-triangle pair duplication,
# resample_na fallback instead of NA, reverse-keying before pairing).

source(test_path("..", "reference", "ref-psychsyn-meade-craig-2012.R"))

# A reproducible matrix with genuine synonym pairs: three latent factors each
# driving a cluster of `per` items with small noise, so within-cluster items
# correlate well above 0.60 (about 0.9) and cross-cluster items near 0. The
# per-respondent within-pair agreement varies with the noise, so the scored
# values span a range (not a constant +1) -- needed for the direction test.
syn_matrix <- function(n = 60L, per = 3L, seed = 42L, noise = 0.35) {
  withr::with_seed(seed, {
    mk <- function() {
      f <- stats::rnorm(n)
      vapply(seq_len(per), function(j) f + stats::rnorm(n, 0, noise),
             numeric(n))
    }
    x <- cbind(mk(), mk(), mk())
  })
  storage.mode(x) <- "double"
  x
}

# Orthogonal-contrast hand fixture (n = 4 respondents, 6 items in 3 synonym
# pairs). Items {1,2}, {3,4}, {5,6} are equal columns (within-pair r = 1); the
# three base columns are the mutually orthogonal Hadamard contrasts of order 4,
# so every cross-pair correlation is exactly 0. The qualifying pair set is
# therefore EXACTLY {(2,1), (4,3), (6,5)} -- worked by hand, no oracle. Since
# each pair's two columns are identical, every respondent's stacked pair vectors
# are equal, so the synonym correlation is +1 for a varying row and NA for the
# all-constant straightliner (row 1).
hand_fixture <- function() {
  a <- c(4, 4, 2, 2)   # contrast (1, 1, -1, -1)
  b <- c(4, 2, 4, 2)   # contrast (1, -1, 1, -1)
  c3 <- c(4, 2, 2, 4)  # contrast (1, -1, -1, 1)
  m <- cbind(a, a, b, b, c3, c3)
  dimnames(m) <- NULL
  storage.mode(m) <- "double"
  m
}

# ---- Schema -----------------------------------------------------------------

test_that("cier_psychsyn returns a list-based cier_index with the pinned schema", {
  out <- cier_psychsyn(syn_matrix())
  expect_s3_class(out, "cier_index")
  expect_type(out, "list")
  expect_identical(names(out),
                   c("value", "flagged", "method", "cutoff", "direction"))
  expect_type(out$value, "double")
  expect_type(out$flagged, "logical")
  expect_identical(length(out$value), 60L)
  expect_identical(length(out$flagged), 60L)
  expect_identical(out$method, "cier_psychsyn")
  expect_identical(out$direction, "lower")
  expect_type(out$cutoff, "double")
})

test_that("as.data.frame.cier_index returns the tidy per-respondent frame", {
  df <- as.data.frame(cier_psychsyn(syn_matrix(n = 12L)))
  expect_s3_class(df, "data.frame")
  expect_identical(names(df), c("value", "flagged"))
  expect_identical(nrow(df), 12L)
})

# ---- Independent oracle parity (1e-12) --------------------------------------

test_that("cier_psychsyn$value equals the oracle on a complete matrix", {
  x <- syn_matrix(n = 50L, seed = 2026L)
  expect_equal(cier_psychsyn(x)$value, ref_psychsyn(x, 0.60),
               tolerance = 1e-12)
})

test_that("cier_psychsyn$value equals the oracle when rows carry NAs", {
  x <- syn_matrix(n = 50L, seed = 99L)
  x[3L, c(1L, 4L)] <- NA           # a few cells dropped, still scored
  x[10L, 1L:8L] <- NA              # most pairs gone -> too few -> abstains
  x[20L, ] <- NA                   # all-NA respondent -> abstains
  ours <- cier_psychsyn(x)$value
  expect_equal(ours, ref_psychsyn(x, 0.60), tolerance = 1e-12)
  expect_true(is.na(ours[[20L]]))
})

test_that("cier_psychsyn$value equals the hand-computed fixture exactly", {
  # Pair set is the three orthogonally-separated synonym pairs; every non-constant
  # respondent scores exactly +1, the all-constant row abstains (NA).
  expect_equal(cier_psychsyn(hand_fixture())$value, c(NA, 1, 1, 1),
               tolerance = 1e-12)
})

# ---- Cross-package parity: careless::psychsyn (bytewise) --------------------

test_that("cier_psychsyn matches careless::psychsyn on careless_dataset (1e-12)", {
  # GENUINE cross-package parity: careless::psychsyn is an independent pair
  # search + stacked-pair correlation. careless::syn_for_one scores each row as
  # cor(cbind(a, b), use = "pairwise.complete.obs")[1, 2] over the FULL (not
  # complete-filtered) pair vectors. The kernel computes that same correlation as
  # a vectorised masked-sum (one rowSums pass, not a per-row cor()), which agrees
  # to ~1e-13 -- held at 1e-12, not bytewise, the difference being summation order
  # only (the per-row-cor() kernel was bytewise; it was vectorised for a 4-6x
  # speedup). The claim holds for resample_na = FALSE (the
  # careless default resample_na = TRUE is a non-deterministic NA fallback cier
  # deliberately does not reproduce). The dataset (careless 1.2.2) carries 53
  # synonym pairs at r > 0.60 and two natural NA respondents, so this also pins NA
  # agreement. See tests/reference/TOLERANCES.md.
  skip_if_not_installed("careless")
  raw <- careless::careless_dataset
  responses <- unname(as.matrix(raw))
  storage.mode(responses) <- "double"
  ours <- cier_psychsyn(responses)$value
  theirs <- careless::psychsyn(raw, critval = 0.60, resample_na = FALSE)
  expect_equal(ours, as.numeric(theirs), tolerance = 1e-12)
})

# ---- Pair discovery / mutant-killers ----------------------------------------

test_that("find_item_pairs selects exactly the synonym pairs (no cross, no dupes)", {
  # The orthogonal hand fixture has EXACTLY three qualifying synonym pairs. A
  # mutant that scans the upper triangle too would return six (each pair twice)
  # and break the oracle / bytewise parity; the count + set pin it here directly.
  pairs <- find_item_pairs(hand_fixture(), 0.60, "syn")
  expect_identical(nrow(pairs), 3L)
  expect_true(all(pairs[, 1L] > pairs[, 2L]))    # larger item index first
  got <- apply(pairs, 1L, function(r) paste(sort(r), collapse = "-"))
  expect_setequal(got, c("1-2", "3-4", "5-6"))
  # The same orthogonal data has no antonym pairs (all cross-correlations are 0).
  expect_identical(nrow(find_item_pairs(hand_fixture(), 0.60, "ant")), 0L)
})

test_that("the synonym / antonym branch selects on the correct sign of r", {
  # Items 1 and 2 are perfectly NEGATIVELY correlated (r = -1); the rest are
  # orthogonal. The antonym branch must select (2, 1); the synonym branch must
  # select nothing. A wrong-sign mutant (< -r vs > r swapped) inverts both. This
  # also exercises the shared kernel's "ant" path that Slice 9 (psychant) wires.
  a <- c(4, 4, 2, 2)
  neg <- c(2, 2, 4, 4)             # = -a (shifted): r(a, neg) = -1
  b <- c(4, 2, 4, 2)
  c3 <- c(4, 2, 2, 4)
  m <- cbind(a, neg, b, c3)
  dimnames(m) <- NULL
  storage.mode(m) <- "double"
  ant <- find_item_pairs(m, 0.60, "ant")
  expect_identical(nrow(ant), 1L)
  expect_identical(as.integer(ant[1L, ]), c(2L, 1L))   # larger item index first
  expect_identical(nrow(find_item_pairs(m, 0.60, "syn")), 0L)
})

test_that("the wrapper is matrix-only: no item-metadata channel exists", {
  # Pins keying-insensitivity structurally: there is no items / scale /
  # reverse_keyed / categories argument, so a reverse-keying-before-pairing
  # mutant cannot be wired in. Pairing uses the raw responses.
  fmls <- names(formals(cier_psychsyn))
  expect_setequal(fmls, c("responses", "critical_r", "fpr", "cutoff"))
  expect_false(any(c("items", "scale", "reverse_keyed",
                     "categories") %in% fmls))
})

test_that("scored synonym values are finite, within [-1, 1], and vary per respondent", {
  v <- cier_psychsyn(syn_matrix(n = 60L, seed = 7L))$value
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
  # A constant row has zero variance on each pair side, so the synonym
  # correlation is NA. A resample_na-fallback mutant would impute a non-NA
  # permutation value instead; we require NA.
  x <- syn_matrix(n = 20L, seed = 5L)
  x[2L, ] <- 3                      # constant straightliner
  out <- cier_psychsyn(x)
  expect_true(is.na(out$value[[2L]]))
  expect_false(is.na(out$value[[1L]]))
})

# ---- Direction (lower) ------------------------------------------------------

test_that("direction is lower: a low-consistency row flags, a consistent one does not", {
  # The least within-pair-consistent respondent is the global minimum -- the most
  # careless -- so it flags at the lower-tail cutoff; the most consistent does
  # not. A flag-upper mutant inverts both.
  out <- cier_psychsyn(syn_matrix(n = 60L, seed = 11L))
  expect_true(out$flagged[[which.min(out$value)]])
  expect_false(out$flagged[[which.max(out$value)]])
  expect_identical(out$flagged, !is.na(out$value) & out$value <= out$cutoff)
})

# ---- Edge cases -------------------------------------------------------------

test_that("critical_r too high finds no pairs: every row abstains and flags nobody", {
  # No inter-item correlation clears 0.99, so there are no synonym pairs, every
  # value is NA, and the percentile cutoff has no finite values to rank: it warns
  # and returns NA, and an NA cutoff flags no one. (The card's documented
  # 100%-NA-at-high-r path.)
  expect_warning(out <- cier_psychsyn(syn_matrix(n = 30L), critical_r = 0.99),
                 class = "cier_warning_insufficient_items")
  expect_true(all(is.na(out$value)))
  expect_true(is.na(out$cutoff))
  expect_true(all(is.na(out$flagged)))
})

test_that("a respondent with fewer than three complete pairs abstains; rows stay aligned", {
  x <- syn_matrix(n = 12L, seed = 3L)
  x[5L, 1L:8L] <- NA                # leaves too few complete pairs -> abstains
  out <- cier_psychsyn(x)
  expect_true(is.na(out$value[[5L]]))
  expect_true(is.na(out$flagged[[5L]]))
  expect_false(is.na(out$value[[1L]]))
  expect_false(is.na(out$value[[12L]]))
})

test_that("the complete-pair boundary is exactly three: three scores, two abstains", {
  # syn_matrix has three clusters of three items (cols 1-3, 4-6, 7-9), and all
  # three within-cluster pairs of each cluster qualify (r ~ 0.89). Leaving only
  # cluster 1 present gives EXACTLY three complete pairs -> scored (n_complete = 3);
  # leaving only its first two items gives one complete pair -> abstains. This
  # pins the n_complete <= 2 -> NA / >= 3 -> scored knife-edge from both sides; an
  # off-by-one in either direction flips one of the two assertions.
  x <- syn_matrix(n = 30L, seed = 8L)
  x[3L, 4L:9L] <- NA               # cluster 1 intact -> exactly 3 complete pairs
  x[7L, 3L:9L] <- NA               # only items 1-2 left -> 1 complete pair
  out <- cier_psychsyn(x)
  expect_false(is.na(out$value[[3L]]))
  expect_true(is.na(out$value[[7L]]))
})

test_that("a constant item column drops out of pairing without leaking a warning", {
  # A zero-variance item makes stats::cor() emit a base-R, locale-dependent
  # "standard deviation is zero" warning; pairing_cor() must suppress it (the
  # package signals only typed cli conditions) and simply not pair the constant
  # item. The remaining synonym pairs still score, so no abstention warning fires
  # either -- the call is warning-free.
  x <- syn_matrix(n = 40L, seed = 6L)
  x[, 5L] <- 3                       # a constant (straightlined) item column
  expect_no_warning(out <- cier_psychsyn(x, critical_r = 0.5))
  expect_gt(sum(!is.na(out$value)), 0L)
})

test_that("a non-matrix / non-numeric / non-finite payload is a typed input error", {
  expect_error(cier_psychsyn(1:10), class = "cier_error_input")
  expect_error(cier_psychsyn(matrix(letters[1:6], nrow = 2L)),
               class = "cier_error_input")
  bad <- syn_matrix(n = 6L)
  bad[1L, 1L] <- Inf
  expect_error(cier_psychsyn(bad), class = "cier_error_input")
})

# ---- Cutoff: default, fpr override, NO-FLIP direction -----------------------

test_that("default cutoff is the lower-tail 5th percentile (NO-FLIP)", {
  out <- cier_psychsyn(syn_matrix(n = 80L, seed = 5L))
  # Lower direction takes the fpr quantile directly (NOT 1 - fpr): the registry
  # stores the literal directional quantile and the kernel must not re-flip.
  expect_equal(out$cutoff,
               as.numeric(stats::quantile(out$value[is.finite(out$value)], 0.05,
                                          names = FALSE, type = 7L)),
               tolerance = 1e-12)
})

test_that("the fpr argument moves the percentile target", {
  x <- syn_matrix(n = 80L, seed = 5L)
  out <- cier_psychsyn(x, fpr = 0.10)
  expect_equal(out$cutoff,
               as.numeric(stats::quantile(out$value[is.finite(out$value)], 0.10,
                                          names = FALSE, type = 7L)),
               tolerance = 1e-12)
})

test_that("an absolute cutoff overrides the percentile and flags via direction (<= ties)", {
  x <- syn_matrix(n = 40L, seed = 7L)
  v <- cier_psychsyn(x)$value
  k <- which.min(v)
  out <- cier_psychsyn(x, cutoff = v[[k]])
  expect_identical(out$cutoff, v[[k]])
  expect_true(out$flagged[[k]])     # exactly at the cutoff flags (<=, not <)
})

test_that("invalid critical_r values are typed input errors", {
  x <- syn_matrix(n = 10L)
  expect_error(cier_psychsyn(x, critical_r = 0), class = "cier_error_input")
  expect_error(cier_psychsyn(x, critical_r = 1), class = "cier_error_input")
  expect_error(cier_psychsyn(x, critical_r = -0.5), class = "cier_error_input")
  expect_error(cier_psychsyn(x, critical_r = c(0.3, 0.4)),
               class = "cier_error_input")
  expect_error(cier_psychsyn(x, critical_r = NA_real_),
               class = "cier_error_input")
})

test_that("invalid fpr / cutoff values are typed input errors", {
  x <- syn_matrix(n = 10L)
  expect_error(cier_psychsyn(x, fpr = 0), class = "cier_error_input")
  expect_error(cier_psychsyn(x, fpr = 1), class = "cier_error_input")
  expect_error(cier_psychsyn(x, fpr = c(0.05, 0.1)), class = "cier_error_input")
  # A synonym value is a correlation in [-1, 1]; a threshold outside that range
  # is degenerate (flags everyone or no one), so it is rejected.
  expect_error(cier_psychsyn(x, cutoff = -1.5), class = "cier_error_input")
  expect_error(cier_psychsyn(x, cutoff = 1.5), class = "cier_error_input")
  expect_error(cier_psychsyn(x, cutoff = NA_real_), class = "cier_error_input")
})

test_that("supplying both fpr and cutoff is a typed input error", {
  expect_error(cier_psychsyn(syn_matrix(n = 10L), fpr = 0.1, cutoff = 0.5),
               class = "cier_error_input")
})

# ---- print snapshot (locked, design-first; direction = lower) ---------------

test_that("print renders the locked cli summary (lower direction)", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    expect_snapshot(print(cier_psychsyn(syn_matrix(n = 30L, seed = 11L))))
  })
})

test_that("print reports abstaining respondents as '(no score)'", {
  # A psychsyn respondent can abstain despite answering (too few complete pairs),
  # so the shared abstention line reads '(no score)', not '(no responses)'.
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    x <- syn_matrix(n = 29L, seed = 11L)
    x <- rbind(x, rep(NA_real_, ncol(x)))   # one abstaining respondent
    expect_snapshot(print(cier_psychsyn(x)))
  })
})
