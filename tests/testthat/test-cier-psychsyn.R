# Tests for cier_psychsyn() -- the psychometric-synonyms C/IER index.
#
# Trust model: the independent oracle (ref_psychsyn) re-derives each respondent's
# correlation across the strongly-positively-correlated item pairs (nested-loop
# pair search + per-row stacked-pair stats::cor()) and never calls the production
# kernel (which discovers pairs from the column-major lower triangle and scores
# every respondent in one vectorised masked-sum pass). The kernel re-sums the same
# Pearson correlation in a different order, so it matches the oracle to ~1e-13,
# held at 1e-12. A separate orthogonal-contrast hand fixture pins the exact pair
# SET and exact +/-1 values with no oracle machinery. The cross-package check pins
# parity with careless::psychsyn(resample_na = FALSE) at 1e-12 (careless's per-row
# cor() over the same stacked pair vectors; no longer bytewise after the kernel was
# vectorised for a 4-6x speedup; the careless default resample_na = TRUE is a
# non-deterministic NA fallback deliberately not reproduced). The property /
# regression block targets each documented mutant (antonym threshold instead of
# synonym, upper-triangle pair duplication, resample_na fallback instead of NA,
# reverse-keying before pairing).

source(test_path("..", "reference", "ref-psychsyn-meade-craig-2012.R"))

# A reproducible matrix with genuine synonym pairs: three latent factors each
# driving a cluster of `per` items with small noise, so within-cluster items
# correlate well above 0.60 (about 0.9) and cross-cluster items near 0. Within-pair
# agreement varies with the noise, so the scored values span a range (not a
# constant +1) -- needed for the direction test.
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
# three base columns are the mutually orthogonal Hadamard contrasts of order 4
# (every cross-pair correlation exactly 0). The qualifying pair set is therefore
# EXACTLY {(2,1), (4,3), (6,5)} -- worked by hand, no oracle. Each pair's columns
# are identical, so every respondent's stacked pair vectors are equal: synonym
# correlation +1 for a varying row, NA for the all-constant straightliner (row 1).
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

test_that("cier_psychsyn returns the pinned cier_index schema", {
  expect_cier_index_schema(cier_psychsyn(syn_matrix()),
                           "cier_psychsyn", "lower", 60L)
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
  expect_equal(suppressWarnings(cier_psychsyn(hand_fixture()))$value,
               c(NA, 1, 1, 1), tolerance = 1e-12)
})

# ---- Cross-package parity: careless::psychsyn (bytewise) --------------------

test_that("cier_psychsyn matches careless::psychsyn on careless_dataset (1e-12)", {
  # GENUINE cross-package parity: careless::syn_for_one scores each row as
  # cor(cbind(a, b), use = "pairwise.complete.obs")[1, 2] over the FULL (not
  # complete-filtered) pair vectors. The vectorised masked-sum kernel (one rowSums
  # pass, not a per-row cor()) agrees to ~1e-13, held at 1e-12 -- summation order
  # only, not bytewise after the 4-6x-speedup vectorisation. Holds for
  # resample_na = FALSE (careless's default resample_na = TRUE is a
  # non-deterministic NA fallback cier does not reproduce). The dataset (careless
  # 1.2.2) carries 53 synonym pairs at r > 0.60 and two natural NA respondents, so
  # this also pins NA agreement.
  skip_if_not_installed("careless")
  raw <- careless::careless_dataset
  responses <- unname(as.matrix(raw))
  storage.mode(responses) <- "double"
  ours <- cier_psychsyn(responses)$value
  theirs <- careless::psychsyn(raw, critval = 0.60, resample_na = FALSE)
  expect_equal(ours, as.numeric(theirs), tolerance = 1e-12)
})

# ---- Pair discovery / regression --------------------------------------------

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

test_that("an injected cor_mat reproduces the default path byte-identically", {
  # cor_mat lets a multi-threshold caller (the critical_r sweep, the wrappers'
  # shared discovery + scoring) build the p x p pairing correlation once; injecting
  # pairing_cor(x) must change nothing.
  x <- syn_matrix(n = 30L, seed = 4L)
  cm <- pairing_cor(x)
  expect_identical(find_item_pairs(x, 0.60, "syn", cor_mat = cm),
                   find_item_pairs(x, 0.60, "syn"))
  expect_identical(kernel_psychsyn(x, 0.60, "syn", cor_mat = cm),
                   kernel_psychsyn(x, 0.60, "syn"))
})

test_that("chunked synonym scoring matches the default kernel", {
  x <- syn_matrix(n = 30L, seed = 4L)
  expect_equal(kernel_psychsyn(x, 0.60, "syn", chunk_cells = 60L),
               kernel_psychsyn(x, 0.60, "syn"),
               tolerance = 1e-12)
})

test_that("the synonym / antonym branch selects on the correct sign of r", {
  # Items 1 and 2 are perfectly NEGATIVELY correlated (r = -1); the rest are
  # orthogonal. The antonym branch must select (2, 1); the synonym branch must
  # select nothing. A wrong-sign mutant (< -r vs > r swapped) inverts both. This
  # also exercises the shared kernel's "ant" path that psychant wires.
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
  # Pins keying-insensitivity structurally: no items / scale / reverse_keyed / max
  # argument, so a reverse-keying-before-pairing mutant cannot be wired in (pairing
  # uses raw responses). The opt-in `reference` is a pairing-SAMPLE selector (which
  # rows discover pairs), NOT item metadata -- the only knob on the matrix-only
  # surface, leaving keying-insensitivity unaffected.
  fmls <- names(formals(cier_psychsyn))
  expect_setequal(fmls, c("responses", "critical_r", "fpr", "cutoff",
                          "reference"))
  expect_false(any(c("items", "scale", "reverse_keyed",
                     "max") %in% fmls))
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
  out <- suppressWarnings(cier_psychsyn(x))
  expect_true(is.na(out$value[[2L]]))
  expect_false(is.na(out$value[[1L]]))
})

test_that("a NON-INTEGER straightliner abstains exactly, with no leaked warning or score", {
  # An integer constant cancels the pair-side deviation sum-of-squares to exactly
  # 0, but a decimal constant lands it a few ulp on EITHER side of zero. The
  # pmax(, 0) clamp catches only the tiny-NEGATIVE side (den 0 -> NaN -> NA);
  # tiny-POSITIVE used to leak a spurious finite ~1.0 -- a perfect straightliner
  # score -- into the percentile pool and flag count. The kernel now detects a
  # constant pair side exactly (masked min == max), so every constant abstains
  # regardless of which way the cancellation fell.
  for (const in c(2.597092, 1.663422, 0.3, 0.05)) {
    # 25 rows keep the cutoff resolvable (>= 20 scored), so expect_no_warning
    # isolates the kernel's silence on the constant pair side from the cutoff layer.
    x <- syn_matrix(n = 25L, seed = 5L)
    x[2L, ] <- const                 # constant DECIMAL straightliner
    expect_no_warning(out <- cier_psychsyn(x))
    expect_true(is.na(out$value[[2L]]))
    expect_true(is.na(out$flagged[[2L]]))
    expect_false(is.na(out$value[[1L]]))
  }
})

# ---- Direction (lower) ------------------------------------------------------

test_that("direction is lower: a low-consistency row flags, a consistent one does not", {
  # The least within-pair-consistent respondent is the global minimum (most
  # careless), so it flags at the lower-tail cutoff; the most consistent does not.
  # A flag-upper mutant inverts both.
  out <- cier_psychsyn(syn_matrix(n = 60L, seed = 11L))
  expect_true(out$flagged[[which.min(out$value)]])
  expect_false(out$flagged[[which.max(out$value)]])
  expect_identical(out$flagged, !is.na(out$value) & out$value <= out$cutoff)
})

# ---- Edge cases -------------------------------------------------------------

test_that("critical_r too high finds no pairs: every row abstains and flags nobody", {
  # No inter-item correlation clears 0.99: no synonym pairs, every value NA, and the
  # percentile cutoff has no finite values to rank -- it warns and returns NA, and
  # an NA cutoff flags no one (the documented 100%-NA-at-high-r path).
  expect_warning(out <- cier_psychsyn(syn_matrix(n = 30L), critical_r = 0.99),
                 class = "cier_warning_insufficient_items")
  expect_true(all(is.na(out$value)))
  expect_true(is.na(out$cutoff))
  expect_true(all(is.na(out$flagged)))
})

test_that("the no-pairs warning is the tailored one: cause, remedy, and ONE warning only", {
  # The generic percentile abstention names neither cause nor fix. The no-pairs
  # path replaces it with a cier_warning_no_pairs naming critical_r, the strongest
  # in-tail r, and the cier_psychsyn_critval() sweep, muffling the redundant generic
  # warning so exactly one reaches the user. The subclass still carries
  # cier_warning_insufficient_items so cier_screen()'s targeted muffler covers it.
  w <- testthat::capture_warnings(out <- cier_psychsyn(syn_matrix(n = 30L),
                                                       critical_r = 0.99))
  expect_length(w, 1L)
  expect_match(w, "No synonym pairs clear")
  expect_match(w, "cier_psychsyn_critval")
  # The synonym sweep hint stays bare -- only the antonym path adds
  # `antonym = TRUE` (asserted in test-cier-psychant.R).
  expect_false(any(grepl("antonym", w, fixed = TRUE)))
  cond <- tryCatch(cier_psychsyn(syn_matrix(n = 30L), critical_r = 0.99),
                   warning = function(w) w)
  expect_s3_class(cond, "cier_warning_no_pairs")
  expect_identical(cier_condition_data(cond)$critical_r, 0.99)
  # A respondent-level abstention with pairs PRESENT keeps the generic warning (the
  # tailored message would mis-state the cause): only cluster 1's first two items
  # survive, so ONE pair qualifies but every respondent has a single complete pair
  # (k = 1 <= 2) and abstains.
  x <- syn_matrix(n = 30L, seed = 3L)
  x[, 3L:9L] <- NA_real_
  cond_na <- tryCatch(cier_psychsyn(x), warning = function(w) w)
  expect_false(inherits(cond_na, "cier_warning_no_pairs"))
  expect_s3_class(cond_na, "cier_warning_insufficient_items")
})

test_that("a respondent with fewer than three complete pairs abstains; rows stay aligned", {
  x <- syn_matrix(n = 12L, seed = 3L)
  x[5L, 1L:8L] <- NA                # leaves too few complete pairs -> abstains
  out <- suppressWarnings(cier_psychsyn(x))
  expect_true(is.na(out$value[[5L]]))
  expect_true(is.na(out$flagged[[5L]]))
  expect_false(is.na(out$value[[1L]]))
  expect_false(is.na(out$value[[12L]]))
})

test_that("the complete-pair boundary is exactly three: three scores, two abstains", {
  # syn_matrix has three clusters of three items (cols 1-3, 4-6, 7-9); all three
  # within-cluster pairs qualify (r ~ 0.89). Leaving only cluster 1 gives EXACTLY
  # three complete pairs -> scored (n_complete = 3); leaving only its first two
  # items gives one complete pair -> abstains. Pins the n_complete <= 2 -> NA /
  # >= 3 -> scored knife-edge from both sides; an off-by-one flips one assertion.
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
  # package signals only typed cli conditions) and just not pair the constant item.
  # The remaining synonym pairs still score, so no abstention warning fires either
  # -- the call is warning-free.
  x <- syn_matrix(n = 40L, seed = 6L)
  x[, 5L] <- 3                       # a constant (straightlined) item column
  expect_no_warning(out <- cier_psychsyn(x, critical_r = 0.5))
  expect_gt(sum(!is.na(out$value)), 0L)
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

test_that("a bad payload or fpr / cutoff override is a typed input error", {
  # Shared with cier_psychant via helper-pairing.R: identical check_responses +
  # check_percentile_overrides contract (a literal cutoff is a correlation in [-1, 1]).
  expect_pairing_overrides_rejected(cier_psychsyn, syn_matrix(n = 10L))
})

# ---- Clean-reference pairing ------------------------------------------------

# A contaminated battery: a careful arm (30 rows) carrying syn_matrix's three
# strong synonym clusters (strongest inter-item r ~ 0.94) plus a careless arm
# (90 rows) of uniform noise over the same value range. The careless majority
# dilutes the WHOLE-sample inter-item correlations to ~0.21 (below any sane
# critical_r), so whole-sample self-pairing finds ZERO pairs and the index abstains
# for everyone -- the documented contamination failure. Discovering pairs on the
# careful arm (reference) restores them; `careful` is the study's leakage-protocol
# mask.
contaminated_syn <- function(n_careful = 30L, n_careless = 90L, seed = 13L) {
  careful_mat <- syn_matrix(n = n_careful, per = 3L, seed = 101L, noise = 0.35)
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
  # Self-reference must change NOTHING: NULL, a full TRUE mask, all row indices, and
  # the analysis matrix itself all discover pairs on the same data, so values are
  # byte-identical to the no-reference call. Pins that `reference` only moves pair
  # DISCOVERY, never scoring (always the full sample) -- default path unchanged.
  x <- syn_matrix(n = 50L, seed = 2026L)
  base <- cier_psychsyn(x)$value
  expect_identical(cier_psychsyn(x, reference = NULL)$value, base)
  expect_identical(cier_psychsyn(x, reference = rep(TRUE, nrow(x)))$value, base)
  expect_identical(cier_psychsyn(x, reference = seq_len(nrow(x)))$value, base)
  expect_identical(cier_psychsyn(x, reference = x)$value, base)
})

test_that("reference discovers pairs on the subset but scores the full sample (oracle, 1e-12)", {
  # Pairs come from the careful arm; all 120 respondents are scored on those pairs.
  # The independent oracle re-derives the construction (pairs from cor(reference),
  # per-row stacked-pair correlation over the full sample) and NEVER calls the
  # production resolver. A mutant scoring on the reference subset would return 30
  # values, not 120, and miss the oracle.
  fx <- contaminated_syn()
  out <- cier_psychsyn(fx$x, critical_r = 0.50, reference = fx$careful)
  expect_identical(length(out$value), nrow(fx$x))
  expect_equal(out$value,
               ref_psychsyn_ref(fx$x, fx$x[fx$careful, ], 0.50),
               tolerance = 1e-12)
  # An external clean SAMPLE (a separate matrix, not a subset) is the same contract.
  ext <- syn_matrix(n = 40L, seed = 77L)               # clean reference sample
  y <- syn_matrix(n = 25L, seed = 5L)
  expect_equal(cier_psychsyn(y, reference = ext)$value,
               ref_psychsyn_ref(y, ext, 0.60), tolerance = 1e-12)
})

test_that("a clean reference recovers pairs that whole-sample pairing loses to contamination", {
  # The documented contamination finding: at the relaxed critical_r = 0.50 the
  # whole-sample inter-item correlations (~0.21) clear no pair, so the default call
  # abstains for EVERYONE and warns; the careful-arm reference surfaces the pairs,
  # scores all 120, and the careless rows -- whose synonym answers decouple -- fall
  # on the flagged lower tail.
  fx <- contaminated_syn()
  cond <- tryCatch(cier_psychsyn(fx$x, critical_r = 0.50), warning = function(w) w)
  expect_s3_class(cond, "cier_warning_no_pairs")
  default <- suppressWarnings(cier_psychsyn(fx$x, critical_r = 0.50))
  expect_true(all(is.na(default$value)))                # whole-sample: nobody scored
  out <- cier_psychsyn(fx$x, critical_r = 0.50, reference = fx$careful)
  expect_identical(sum(!is.na(out$value)), nrow(fx$x))  # reference: everyone scored
  expect_lt(mean(out$value[!fx$careful]), mean(out$value[fx$careful]))  # direction
  expect_gt(sum(out$flagged[!fx$careful], na.rm = TRUE), 0L)            # careless flag
  expect_identical(sum(out$flagged[fx$careful], na.rm = TRUE), 0L)      # careful do not
})

test_that("the no-pairs warning's strongest_r comes from the reference, not the analysis data", {
  # When the reference clears no pair, the warning must name the reference's
  # strongest in-tail correlation (where discovery happened), not the analysis
  # sample's. Here the analysis data carries strong synonym structure but the
  # reference (noise) does not, so the fired warning's strongest_r reflects the
  # weak reference.
  x <- syn_matrix(n = 40L, seed = 7L)                  # strong structure
  ref <- withr::with_seed(3L, matrix(stats::runif(30L * ncol(x), 0, 1),
                                     nrow = 30L, ncol = ncol(x)))  # weak noise ref
  cond <- tryCatch(cier_psychsyn(x, critical_r = 0.60, reference = ref),
                   warning = function(w) w)
  expect_s3_class(cond, "cier_warning_no_pairs")
  strongest_ref <- max(pairing_cor(ref)[lower.tri(diag(ncol(x)))], na.rm = TRUE)
  expect_equal(cier_condition_data(cond)$strongest_r, strongest_ref,
               tolerance = 1e-12)
  expect_lt(cier_condition_data(cond)$strongest_r, 0.60)  # not the strong analysis r
})

test_that("invalid reference selectors are typed input errors", {
  x <- syn_matrix(n = 20L, seed = 4L)
  expect_pairing_reference_rejected(cier_psychsyn, x)   # shared selector validation
  # A malformed external sample must blame `reference`, not the valid `responses`:
  # the typed condition's $arg names the offending argument.
  cond <- tryCatch(cier_psychsyn(x, reference = matrix("a", 5L, ncol(x))),
                   error = function(e) e)
  expect_identical(cier_condition_data(cond)$arg, "reference")
})

test_that("a reference of exactly three rows is accepted (the >= 3 boundary, accept side)", {
  # The rejected side is pinned above (2 rows -> error); this pins the ACCEPT side so
  # a mutant demanding >= 4 rows cannot hide. Both the subset-index and
  # external-matrix forms select exactly three reference rows and must score (or
  # abstain) WITHOUT a too-few-rows abort -- under a >= 4 mutant the wrapper would
  # raise cier_error_input here and the value comparison would error.
  x <- syn_matrix(n = 25L, seed = 4L)
  expect_equal(suppressWarnings(cier_psychsyn(x, reference = 1:3L)$value),
               ref_psychsyn_ref(x, x[1:3, ], 0.60), tolerance = 1e-12)
  expect_equal(suppressWarnings(cier_psychsyn(x, reference = x[1:3, ])$value),
               ref_psychsyn_ref(x, x[1:3, ], 0.60), tolerance = 1e-12)
})

test_that("a double (non-integer-typed) index vector is a subset selector, not a sample", {
  # Real callers pass c(1, 2, 3) (double), not seq_len()'s integer. A double vector
  # has no dim, so it is a SUBSET selector (row indices), never an external sample.
  # A mutant accepting only integer-typed indices would route the double to the
  # external-sample branch -> check_responses on a bare vector -> error.
  fx <- contaminated_syn()
  idx <- as.double(which(fx$careful))
  expect_equal(cier_psychsyn(fx$x, critical_r = 0.50, reference = idx)$value,
               ref_psychsyn_ref(fx$x, fx$x[fx$careful, ], 0.50), tolerance = 1e-12)
  expect_error(cier_psychsyn(fx$x, reference = c(1, 2, 999)),   # double out-of-range
               class = "cier_error_input")
})

test_that("a column-reordered external reference is aligned by name, not silently mispaired", {
  # Pairing is positional, so a same-items sample in a different column order would
  # discover the WRONG pairs. When both carry names, the sample is reordered to the
  # analysis column order: a reshuffled sample must give the SAME result as the
  # in-order one (and match the oracle on the in-order sample).
  x <- syn_matrix(n = 30L, seed = 4L)
  colnames(x) <- paste0("i", seq_len(ncol(x)))
  ext <- syn_matrix(n = 40L, seed = 8L)
  colnames(ext) <- colnames(x)
  perm <- c(3L, 1L, 2L, 6L, 4L, 5L, 9L, 7L, 8L)     # a within/cross-cluster shuffle
  ext_shuffled <- ext[, perm]
  expect_equal(cier_psychsyn(x, reference = ext_shuffled)$value,
               cier_psychsyn(x, reference = ext)$value, tolerance = 1e-12)
  expect_equal(cier_psychsyn(x, reference = ext_shuffled)$value,
               ref_psychsyn_ref(x, ext, 0.60), tolerance = 1e-12)
})

test_that("a named external reference measuring different items is a typed error (same ncol)", {
  # An ncol-only check would pass this (same item COUNT); name alignment catches the
  # differing item SETS that positional pairing would silently mispair.
  x <- syn_matrix(n = 30L, seed = 4L)
  colnames(x) <- paste0("i", seq_len(ncol(x)))
  ext <- syn_matrix(n = 40L, seed = 8L)
  colnames(ext) <- paste0("j", seq_len(ncol(ext)))
  expect_error(cier_psychsyn(x, reference = ext), class = "cier_error_input")
})

# ---- print snapshot (locked; direction = lower) -----------------------------

test_that("print renders the locked cli summary (lower direction)", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    expect_snapshot(print(cier_psychsyn(syn_matrix(n = 30L, seed = 11L))))
  })
})

test_that("print reports abstaining respondents as '(no score)'", {
  # A psychsyn respondent can abstain despite answering (too few complete pairs),
  # so the abstention line reads '(no score)', not '(no responses)'.
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    x <- syn_matrix(n = 29L, seed = 11L)
    x <- rbind(x, rep(NA_real_, ncol(x)))   # one abstaining respondent
    expect_snapshot(print(cier_psychsyn(x)))
  })
})
