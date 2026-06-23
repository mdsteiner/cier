# Tests for cier_psychant() -- the psychometric-antonyms C/IER index, the mirror
# of cier_psychsyn() over the strongly NEGATIVELY correlated item pairs.
#
# Trust model: the independent oracle (ref_psychant) re-derives each respondent's
# correlation across the antonym pairs (nested-loop pair search r < -|critval|,
# i > j, + per-row stacked-pair stats::cor() reusing ref_psychsyn_row) and never
# calls the production kernel (column-major lower-triangle discovery, one vectorised
# masked-sum pass). The kernel re-sums the same Pearson correlation in a different
# order, matching the oracle to ~1e-13, held at 1e-12. A separate orthogonal-contrast
# hand fixture pins the exact pair SET and scored values (-1) with no oracle
# machinery. The cross-package check pins parity with careless::psychsyn(anto = TRUE,
# resample_na = FALSE) at 1e-12 -- the same shared kernel, so it inherits psychsyn's
# tolerance (no longer bytewise after the 4-6x-speedup vectorisation); careless::
# psychant() does not surface resample_na, so the deterministic comparison calls the
# underlying psychsyn(anto = TRUE, resample_na = FALSE) directly. The property /
# regression block targets each documented mutant (synonym threshold instead of
# antonym, flag lower instead of upper, wrong sign branch, reverse-keying before
# pairing).

source(test_path("..", "reference", "ref-psychsyn-meade-craig-2012.R"))
source(test_path("..", "reference", "ref-psychant-meade-craig-2012.R"))

# A reproducible matrix with genuine antonym pairs: `k` mutually-orthogonal latent
# factors, each driving ONE positively-loading item (cols 1..k) and ONE
# negatively-loading item (cols k+1..2k) with small noise, so each factor's two
# items correlate well below -0.60 (about -0.89 at noise 0.35) and cross-factor
# items near 0. The antonym pairs are therefore EXACTLY (k+j, j) for j = 1..k;
# within-pair opposition varies with the noise, so the scored values span a range
# below 0 (not a constant -1) -- needed for the direction test.
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
# {(2,1), (4,3), (6,5)} -- worked by hand, no oracle -- and NO synonym pair exists
# (no off-diagonal r is positive). For each non-constant respondent the stacked
# larger-item values are the exact negatives of the smaller-item values, so the
# antonym correlation is -1; the all-constant straightliner (row 1) abstains (NA).
# The exact mirror of the psychsyn hand fixture (which scored +1).
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
# careless_dataset has NO antonym pairs at r < -0.60, so genuine parity needs
# constructed negative structure (the empty-antonym dataset is a trivial all-NA
# match).
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

test_that("cier_psychant returns the pinned cier_index schema", {
  expect_cier_index_schema(cier_psychant(ant_matrix()),
                           "cier_psychant", "upper", 60L)
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
  expect_equal(suppressWarnings(cier_psychant(hand_fixture()))$value,
               c(NA, -1, -1, -1), tolerance = 1e-12)
})

# ---- Cross-package parity: careless (anto = TRUE, 1e-12) --------------------

test_that("cier_psychant matches careless antonyms on a planted fixture (1e-12)", {
  # GENUINE cross-package parity: careless's antonym scorer is an independent pair
  # search + stacked-pair correlation. careless::psychant() does not expose
  # resample_na, so the deterministic comparison calls the underlying
  # psychsyn(anto = TRUE, resample_na = FALSE) directly (careless's default
  # resample_na = TRUE is a non-deterministic NA fallback cier does not reproduce).
  # The shared masked-sum kernel computes the same correlation in one vectorised
  # pass, agreeing to ~1e-13, held at 1e-12 (psychsyn's post-vectorisation
  # tolerance). Injected NAs (a few cells + one all-NA row) also pin NA agreement.
  skip_if_not_installed("careless")
  responses <- ant_careless_fixture()
  responses[5L, c(2L, 9L)] <- NA          # a few dropped cells, still scored
  responses[40L, ] <- NA                  # an all-NA respondent -> NA both sides
  ours <- cier_psychant(responses, critical_r = 0.60)$value
  theirs <- careless::psychsyn(responses, critval = -0.60, anto = TRUE,
                               resample_na = FALSE)
  expect_equal(ours, as.numeric(theirs), tolerance = 1e-12)
})

# ---- Pair discovery / regression --------------------------------------------

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
  # Items 1 and 2 are perfectly NEGATIVELY correlated (r = -1); items 3 and 4
  # perfectly POSITIVELY (r = +1); the rest orthogonal. The antonym branch must
  # select ONLY (2, 1); the synonym branch ONLY (4, 3). The core psychant mutant --
  # wrong-sign (> r vs < -r swapped) -- inverts both.
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

test_that("chunked antonym scoring matches the default kernel", {
  x <- ant_matrix(n = 30L, seed = 4L)
  expect_equal(kernel_psychsyn(x, 0.60, "ant", chunk_cells = 60L),
               kernel_psychsyn(x, 0.60, "ant"),
               tolerance = 1e-12)
})

test_that("the wrapper is matrix-only: no item-metadata channel exists", {
  # Pins keying-insensitivity structurally: no items / scale / reverse_keyed / max
  # argument, so a reverse-keying-before-pairing mutant cannot be wired in (pairing
  # uses raw responses; reverse-keying would collapse genuine antonym correlations
  # toward zero). The opt-in `reference` is a pairing-SAMPLE selector (which rows
  # discover pairs), NOT item metadata, so keying-insensitivity is unaffected.
  fmls <- names(formals(cier_psychant))
  expect_setequal(fmls, c("responses", "critical_r", "fpr", "cutoff",
                          "reference"))
  expect_false(any(c("items", "scale", "reverse_keyed",
                     "max") %in% fmls))
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
  out <- suppressWarnings(cier_psychant(x))
  expect_true(is.na(out$value[[2L]]))
  expect_false(is.na(out$value[[1L]]))
})

# ---- Direction (upper) ------------------------------------------------------

test_that("direction is upper: a decoupled row flags, a strongly-opposed one does not", {
  # The least within-pair-opposed respondent has the HIGHEST (least negative)
  # antonym correlation (most careless), so it flags at the upper-tail cutoff; the
  # most strongly opposed (most negative) does not. A flag-lower mutant inverts both.
  out <- cier_psychant(ant_matrix(n = 60L, seed = 11L))
  expect_true(out$flagged[[which.max(out$value)]])
  expect_false(out$flagged[[which.min(out$value)]])
  expect_identical(out$flagged, !is.na(out$value) & out$value >= out$cutoff)
})

# ---- Edge cases -------------------------------------------------------------

test_that("critical_r too high finds no pairs: every row abstains and flags nobody", {
  # No inter-item correlation clears -0.99: no antonym pairs, every value NA, and the
  # percentile cutoff has no finite values to rank -- it warns and returns NA, and
  # an NA cutoff flags no one.
  expect_warning(out <- cier_psychant(ant_matrix(n = 30L), critical_r = 0.99),
                 class = "cier_warning_insufficient_items")
  expect_true(all(is.na(out$value)))
  expect_true(is.na(out$cutoff))
  expect_true(all(is.na(out$flagged)))
})

test_that("the no-pairs warning is the tailored antonym one, raised once", {
  # Mirrors the psychsyn assertion: the antonym wrapper shares the no-pairs tail, so
  # the single warning carries the cier_warning_no_pairs subclass and the antonym
  # noun (not "synonym").
  w <- testthat::capture_warnings(cier_psychant(ant_matrix(n = 30L),
                                                critical_r = 0.99))
  expect_length(w, 1L)
  expect_match(w, "No antonym pairs clear")
  # The sweep hint is tail-aware -- the antonym path points at the antonym
  # sweep, cier_psychsyn_critval(antonym = TRUE), not the bare synonym sweep.
  expect_match(w, "cier_psychsyn_critval(antonym = TRUE)", fixed = TRUE)
  # ...rendered cleanly: a `{.fun ...}` span would append a stray "()" -- use
  # `{.code ...}` so the call is not malformed.
  expect_no_match(w, "TRUE)()", fixed = TRUE)
  cond <- tryCatch(cier_psychant(ant_matrix(n = 30L), critical_r = 0.99),
                   warning = function(w) w)
  expect_s3_class(cond, "cier_warning_no_pairs")
})

test_that("a respondent with fewer than three complete pairs abstains; rows stay aligned", {
  x <- ant_matrix(n = 12L, seed = 3L)
  x[5L, c(2L, 3L, 4L, 6L, 7L, 8L)] <- NA   # leaves only pair (5,1) complete -> abstains
  out <- suppressWarnings(cier_psychant(x))
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

test_that("a bad payload or fpr / cutoff override is a typed input error", {
  # Shared with cier_psychsyn via helper-pairing.R: identical check_responses +
  # check_percentile_overrides contract (a literal cutoff is a correlation in [-1, 1]).
  expect_pairing_overrides_rejected(cier_psychant, ant_matrix(n = 10L))
})

# ---- Clean-reference pairing ------------------------------------------------

# A contaminated antonym battery: a careful arm (30 rows) carrying ant_matrix's
# strong forward/reverse opposition (most-negative inter-item r ~ -0.94) plus a
# careless arm (90 rows) of uniform noise. The careless majority pulls the
# whole-sample most-negative correlation up to ~-0.18, so self-pairing finds NO
# antonym pair and the index abstains for everyone; the careful arm restores them
# -- the antonym mirror of the psychsyn recovery.
contaminated_ant <- function(n_careful = 30L, n_careless = 90L, seed = 14L) {
  careful_mat <- ant_matrix(n = n_careful, k = 4L, seed = 101L, noise = 0.35)
  careless <- withr::with_seed(seed, {
    matrix(stats::runif(n_careless * ncol(careful_mat),
                        min = min(careful_mat), max = max(careful_mat)),
           nrow = n_careless, ncol = ncol(careful_mat))
  })
  x <- rbind(careful_mat, careless)
  storage.mode(x) <- "double"
  list(x = x, careful = c(rep(TRUE, n_careful), rep(FALSE, n_careless)))
}

test_that("reference == the analysis sample reproduces the default path byte-for-byte", {
  x <- ant_matrix(n = 50L, seed = 2026L)
  base <- cier_psychant(x)$value
  expect_identical(cier_psychant(x, reference = NULL)$value, base)
  expect_identical(cier_psychant(x, reference = rep(TRUE, nrow(x)))$value, base)
  expect_identical(cier_psychant(x, reference = seq_len(nrow(x)))$value, base)
  expect_identical(cier_psychant(x, reference = x)$value, base)
})

test_that("reference discovers antonym pairs but scores the full sample (oracle 1e-12)", {
  fx <- contaminated_ant()
  out <- cier_psychant(fx$x, critical_r = 0.50, reference = fx$careful)
  expect_identical(length(out$value), nrow(fx$x))
  expect_equal(out$value,
               ref_psychant_ref(fx$x, fx$x[fx$careful, ], 0.50),
               tolerance = 1e-12)
  ext <- ant_matrix(n = 40L, seed = 77L)               # clean reference sample
  y <- ant_matrix(n = 25L, seed = 5L)
  expect_equal(cier_psychant(y, reference = ext)$value,
               ref_psychant_ref(y, ext, 0.60), tolerance = 1e-12)
})

test_that("a clean reference recovers antonym pairs that contamination loses", {
  # The antonym mirror of the psychsyn recovery: whole-sample pairing finds no
  # antonym pair under contamination (warns, all NA); the careful-arm reference
  # surfaces them, scores all 120, and the careless rows -- whose opposed-item
  # answers decouple toward zero -- fall on the flagged UPPER tail.
  fx <- contaminated_ant()
  cond <- tryCatch(cier_psychant(fx$x, critical_r = 0.50), warning = function(w) w)
  expect_s3_class(cond, "cier_warning_no_pairs")
  default <- suppressWarnings(cier_psychant(fx$x, critical_r = 0.50))
  expect_true(all(is.na(default$value)))
  out <- cier_psychant(fx$x, critical_r = 0.50, reference = fx$careful)
  expect_identical(sum(!is.na(out$value)), nrow(fx$x))
  expect_gt(mean(out$value[!fx$careful]), mean(out$value[fx$careful]))  # direction
  expect_gt(sum(out$flagged[!fx$careful], na.rm = TRUE), 0L)            # careless flag
  expect_identical(sum(out$flagged[fx$careful], na.rm = TRUE), 0L)      # careful do not
})

test_that("the no-pairs warning's strongest_r comes from the reference (most-negative)", {
  # For antonyms the strongest in-tail correlation is the MOST NEGATIVE. When the
  # reference clears no antonym pair, the warning names the reference's most-negative
  # r, not the analysis sample's.
  x <- ant_matrix(n = 40L, seed = 7L)                  # strong opposition
  ref <- withr::with_seed(3L, matrix(stats::runif(30L * ncol(x), 0, 1),
                                     nrow = 30L, ncol = ncol(x)))  # weak noise ref
  cond <- tryCatch(cier_psychant(x, critical_r = 0.60, reference = ref),
                   warning = function(w) w)
  expect_s3_class(cond, "cier_warning_no_pairs")
  strongest_ref <- min(pairing_cor(ref)[lower.tri(diag(ncol(x)))], na.rm = TRUE)
  expect_equal(cier_condition_data(cond)$strongest_r, strongest_ref,
               tolerance = 1e-12)
  expect_gt(cier_condition_data(cond)$strongest_r, -0.60)  # not the strong analysis r
})

test_that("invalid reference selectors are typed input errors", {
  # Shared selector validation with cier_psychsyn via helper-pairing.R.
  expect_pairing_reference_rejected(cier_psychant, ant_matrix(n = 20L, seed = 4L))
})

test_that("a reference of exactly three rows is accepted (the >= 3 boundary, accept side)", {
  # Accept-side boundary mirror: a mutant demanding >= 4 reference rows would abort
  # here. Both the subset-index and external-matrix forms select exactly three rows.
  x <- ant_matrix(n = 25L, seed = 4L)
  expect_equal(suppressWarnings(cier_psychant(x, reference = 1:3L)$value),
               ref_psychant_ref(x, x[1:3, ], 0.60), tolerance = 1e-12)
  expect_equal(suppressWarnings(cier_psychant(x, reference = x[1:3, ])$value),
               ref_psychant_ref(x, x[1:3, ], 0.60), tolerance = 1e-12)
})

test_that("a double (non-integer-typed) index vector is a subset selector, not a sample", {
  fx <- contaminated_ant()
  idx <- as.double(which(fx$careful))
  expect_equal(cier_psychant(fx$x, critical_r = 0.50, reference = idx)$value,
               ref_psychant_ref(fx$x, fx$x[fx$careful, ], 0.50), tolerance = 1e-12)
  expect_error(cier_psychant(fx$x, reference = c(1, 2, 999)),
               class = "cier_error_input")
})

# ---- print snapshot (locked; direction = upper) -----------------------------

test_that("print renders the locked cli summary (upper direction)", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    expect_snapshot(print(cier_psychant(ant_matrix(n = 30L, seed = 11L))))
  })
})

test_that("print reports abstaining respondents as '(no score)'", {
  # A psychant respondent can abstain despite answering (too few complete pairs),
  # so the abstention line reads '(no score)', not '(no responses)'.
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    x <- ant_matrix(n = 29L, seed = 11L)
    x <- rbind(x, rep(NA_real_, ncol(x)))   # one abstaining respondent
    expect_snapshot(print(cier_psychant(x)))
  })
})
