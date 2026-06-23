# Tests for cier_mahalanobis().
#
# Trust model: oracle ref_mahalanobis re-derives each row's squared Mahalanobis
# distance with an explicit per-respondent quadratic-form loop, never calling the
# kernel (which evaluates the whole battery with one BLAS matrix multiply); a
# hand-computed 2-column fixture pins exact values sharing no cov/solve machinery.
# The cross-package check pins parity with careless::mahad AND psych::outlier at
# 1e-10 tolerance.

source(test_path("..", "reference", "ref-mahalanobis-curran-2016.R"))

# A reproducible discrete fixture (30 respondents x 6 items, df = 6 != n).
rand_matrix <- function(n = 30L, p = 6L, seed = 7L) {
  withr::with_seed(seed, {
    x <- matrix(sample.int(5L, n * p, replace = TRUE), nrow = n)
  })
  storage.mode(x) <- "double"
  x
}

# Hand-computed fixture (worked by hand in the oracle header's algebra):
#   x = (0,0),(1,2),(2,2),(1,4); mu = (1, 2); Sigma = [[2/3, 2/3],[2/3, 8/3]];
#   Sigma^-1 = [[2, -1/2],[-1/2, 1/2]]  => D^2 = c(2, 0, 2, 2), sum = (n-1)*p = 6.
hand_fixture <- function() {
  rbind(c(0, 0), c(1, 2), c(2, 2), c(1, 4))
}

# Heavy-missingness fixture whose pairwise covariance is invertible but NOT
# positive definite. Likert 1..5, 30 x 5, 40% MCAR; seed 11 is the first seed in
# 1:500 that, under the default RNG kinds (fresh session, as in testthat / CI),
# gives: no all-NA row, solve() succeeds with a finite inverse, chol() FAILS (min
# eigenvalue ~ -0.12), and the unguarded bilinear form yields a negative D^2 on 16
# rows.
indefinite_fixture <- function() {
  x <- withr::with_seed(11L, {
    m <- matrix(sample.int(5L, 30L * 5L, replace = TRUE), nrow = 30L)
    m[sample.int(length(m), 60L)] <- NA # 60 / 150 = 40% missing
    m
  })
  storage.mode(x) <- "double"
  x
}

# A reproducible partial-row fixture: the canonical seed-7 matrix with three rows
# carrying some-but-not-all NA (rows 3, 10, 20), each SCORED (mean-imputed) not
# abstaining. n_partial = 3, n_scored = 30 (no all-NA row). Same NA layout as the
# zero-fill oracle-parity test, so it is known to score.
partial_fixture <- function() {
  x <- rand_matrix(n = 30L, p = 6L, seed = 7L)
  x[3L, c(1L, 5L)] <- NA
  x[10L, 2L] <- NA
  x[20L, c(2L, 3L, 4L)] <- NA
  x
}

# Collect the most-specific class of every warning a call signals, muffling each,
# so a test can assert one class fires while another does NOT (and count repeats)
# -- expect_warning() only inspects the first match. An environment accumulator
# avoids the package-forbidden `<<-`.
collect_warning_classes <- function(expr) {
  acc <- new.env(parent = emptyenv())
  acc$cls <- character()
  withCallingHandlers(
    expr,
    warning = function(w) {
      acc$cls <- c(acc$cls, class(w)[[1L]])
      invokeRestart("muffleWarning")
    }
  )
  acc$cls
}

# ---- Schema -----------------------------------------------------------------

test_that("cier_mahalanobis returns the pinned cier_index schema", {
  expect_cier_index_schema(cier_mahalanobis(rand_matrix()),
                           "cier_mahalanobis", "upper", 30L)
})

# ---- Independent oracle parity (1e-10) --------------------------------------

test_that("cier_mahalanobis$value equals the hand-computed fixture exactly", {
  # Pins exact hand-derived D^2 values (no shared cov/solve machinery in the
  # expectation): a gross formula error cannot survive this.
  expect_equal(cier_mahalanobis(hand_fixture())$value,
               c(2, 0, 2, 2), tolerance = 1e-12)
})

test_that("cier_mahalanobis$value equals the oracle on a random complete matrix", {
  x <- rand_matrix()
  expect_equal(cier_mahalanobis(x)$value, ref_mahalanobis(x), tolerance = 1e-10)
})

test_that("cier_mahalanobis$value equals the oracle when rows carry NAs (zero-fill path)", {
  x <- rand_matrix()
  x[3L, c(1L, 5L)] <- NA
  x[10L, 2L] <- NA
  x[20L, c(2L, 3L, 4L)] <- NA
  # A na.rm-on-the-bilinear-form mutant collapses these rows to 0 and diverges; a
  # complete-case-cov mutant drops the NA rows from cov() and diverges too.
  # (Scoring partial rows now warns; this asserts the value, not the warning.)
  expect_equal(suppressWarnings(cier_mahalanobis(x))$value,
               ref_mahalanobis(x), tolerance = 1e-10)
})

# ---- Cross-package parity: careless::mahad + psych::outlier (1e-10) ----------

test_that("cier_mahalanobis matches careless::mahad on a random matrix at 1e-10", {
  skip_if_not_installed("careless")
  x <- rand_matrix(n = 40L, p = 8L, seed = 1L)
  expect_equal(cier_mahalanobis(x)$value,
               as.numeric(careless::mahad(x, flag = FALSE, plot = FALSE)),
               tolerance = 1e-10)
})

test_that("cier_mahalanobis matches careless::mahad on careless_dataset at 1e-10", {
  skip_if_not_installed("careless")
  x <- unname(as.matrix(careless::careless_dataset))
  storage.mode(x) <- "double"
  expect_equal(cier_mahalanobis(x)$value,
               as.numeric(careless::mahad(x, flag = FALSE, plot = FALSE)),
               tolerance = 1e-10)
})

test_that("cier_mahalanobis matches psych::outlier (D^2, centre + pairwise cov) at 1e-10", {
  skip_if_not_installed("psych")
  x <- rand_matrix(n = 40L, p = 8L, seed = 1L)
  expect_equal(cier_mahalanobis(x)$value,
               as.numeric(psych::outlier(x, plot = FALSE, bad = 0, na.rm = TRUE)),
               tolerance = 1e-10)
})

test_that("both partners agree on NA-bearing rows (zero-fill, not complete-case)", {
  skip_if_not_installed("careless")
  skip_if_not_installed("psych")
  x <- rand_matrix(n = 40L, p = 8L, seed = 1L)
  x[3L, c(1L, 5L)] <- NA
  x[10L, 2L] <- NA
  ours <- suppressWarnings(cier_mahalanobis(x))$value   # partial rows now warn
  expect_equal(ours, as.numeric(careless::mahad(x, flag = FALSE, plot = FALSE)),
               tolerance = 1e-10)
  expect_equal(ours, as.numeric(psych::outlier(x, plot = FALSE, bad = 0,
                                               na.rm = TRUE)),
               tolerance = 1e-10)
})

# ---- Property / invariant + regression --------------------------------------

test_that("the squared distances obey the trace identity sum(D^2) == (n-1)*p", {
  # A standard result for D^2 from the sample mean with sample covariance; guards
  # the (n-1) denominator (a population-cov mutant dividing by n inflates the sum by
  # n/(n-1)). The centre-only-vs-scale=TRUE mutation (standardising the centred
  # matrix while keeping the raw covariance) is a mismatched quadratic form that
  # oracle parity on random data catches; a fully self-consistent standardisation
  # (correlation matrix) is the SAME squared-distance statistic, nothing to catch.
  x <- rand_matrix(n = 25L, p = 6L, seed = 3L)
  expect_equal(sum(cier_mahalanobis(x)$value), (25 - 1) * 6, tolerance = 1e-8)
})

test_that("scores stay finite through the zero-fill path and non-negative on complete data", {
  # Complete data: pairwise cov == complete-case cov is positive semi-definite, so
  # the quadratic form is genuinely >= 0 here (not a vacuous claim).
  v <- cier_mahalanobis(rand_matrix(n = 50L, p = 6L, seed = 9L))$value
  expect_true(all(is.finite(v)))
  expect_true(all(v >= 0))
  # Partial rows exercise the na.rm / zero-fill path: it must not leak NaN. A mutant
  # that na.rm-collapses the bilinear form (instead of zero-filling the centred
  # coordinates) would turn these rows non-finite.
  x <- rand_matrix(n = 50L, p = 6L, seed = 9L)
  x[c(2L, 7L, 30L), c(1L, 4L)] <- NA
  expect_true(all(is.finite(suppressWarnings(cier_mahalanobis(x))$value)))
})

test_that("direction is upper: the most extreme row flags, a central row does not", {
  # A strong multivariate outlier (alternating extremes) against 24 mid rows. The
  # relaxed tail (alpha = 0.2) yields a genuine flag: the highest-D^2 row clears the
  # cutoff and the most central (lowest-D^2) row does not. A flag-lower mutant would
  # invert both.
  rnd <- rand_matrix(n = 24L, p = 6L, seed = 21L)
  x <- rbind(c(5, 1, 5, 1, 5, 1), rnd)
  out <- cier_mahalanobis(x, alpha = 0.2)
  expect_true(out$flagged[[which.max(out$value)]])  # the most extreme row flags
  expect_false(out$flagged[[which.min(out$value)]]) # the most central row does not
  # The flag is exactly the upper comparator against the cutoff (NA-safe).
  expect_identical(out$flagged, !is.na(out$value) & out$value >= out$cutoff)
})

# ---- Edge cases -------------------------------------------------------------

test_that("an all-NA row abstains and keeps the remaining rows aligned", {
  # The abstaining row sits in the middle, so value/flagged must stay aligned to
  # their respondents (row-indexing mutant guard).
  x <- rand_matrix(n = 10L, p = 6L, seed = 4L)
  x[5L, ] <- NA
  # A single all-NA respondent is silent per-row abstention, NOT the wholesale
  # degenerate-covariance case: it must NOT emit the singular-covariance warning (an
  # over-eager mutant that warns on every all-NA row dies here).
  expect_no_warning(out <- cier_mahalanobis(x))
  expect_true(is.na(out$value[[5L]]))
  expect_true(is.na(out$flagged[[5L]]))
  expect_false(is.na(out$value[[1L]]))
  expect_false(is.na(out$value[[10L]]))
  # The scored rows still match the oracle (which also abstains on the all-NA row).
  expect_equal(out$value, ref_mahalanobis(x), tolerance = 1e-10)
})

test_that("fewer than two respondents with data abstains with a typed warning", {
  x <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 1L)   # one respondent
  expect_warning(out <- cier_mahalanobis(x),
                 class = "cier_warning_singular_covariance")
  expect_true(all(is.na(out$value)))
  expect_true(all(is.na(out$flagged)))
})

test_that("a singular covariance abstains with a typed warning and flags nobody", {
  # A perfectly collinear (duplicated) column makes Sigma singular: solve() fails.
  x <- rand_matrix(n = 30L, p = 5L, seed = 6L)
  x <- cbind(x, x[, 1L])                          # duplicate column 1
  expect_warning(out <- cier_mahalanobis(x),
                 class = "cier_warning_singular_covariance")
  expect_true(all(is.na(out$value)))
  expect_true(all(is.na(out$flagged)))
})

test_that("an NA-bearing pairwise covariance also abstains with the typed warning", {
  # Distinct path from the duplicated-column case above: when two items are never
  # co-answered, cov(use = "pairwise") has an NA cell and solve() returns an all-NA
  # inverse WITHOUT erroring -- so the guard must also check anyNA(), not only catch
  # a solve() error. (An item answered by nobody is the same failure.) Every
  # respondent abstains; no one is flagged.
  x <- matrix(c(
    1, NA,
    2, NA,
    NA, 4,
    NA, 5,
    3, NA,
    NA, 2
  ), ncol = 2L, byrow = TRUE)                     # items 1 and 2 never co-answered
  expect_warning(out <- cier_mahalanobis(x),
                 class = "cier_warning_singular_covariance")
  w <- tryCatch(cier_mahalanobis(x), warning = function(w) w)
  expect_identical(cier_condition_data(w)$reason, "singular_covariance")
  expect_true(all(is.na(out$value)))
  expect_true(all(is.na(out$flagged)))
})

test_that("an indefinite pairwise covariance abstains wholesale with the typed warning", {
  # Heavy MCAR missingness can leave cov(use = "pairwise") INVERTIBLE yet not
  # positive definite (a negative eigenvalue): solve() succeeds with a finite
  # inverse, so neither earlier guard fires, but the bilinear form is then signed --
  # a respondent can score a NEGATIVE "squared distance" the upper-tail chi-square
  # flag can never reach, and the ranking among positive rows is distorted too. The
  # distance is invalid for EVERY row, so the kernel abstains wholesale.
  # careless::mahad reproduces the signed values on this input, so cross-package
  # parity cannot guard this path -- this regression test is the only guard (fixture
  # verified: solve() succeeds, chol() fails, the unguarded form yields negative
  # D^2).
  x <- indefinite_fixture()
  expect_warning(out <- cier_mahalanobis(x),
                 class = "cier_warning_singular_covariance")
  w <- tryCatch(cier_mahalanobis(x), warning = function(w) w)
  expect_identical(cier_condition_data(w)$reason, "indefinite_covariance")
  expect_true(all(is.na(out$value)))
  expect_true(all(is.na(out$flagged)))
})

test_that("scored Mahalanobis values are never negative (PD-guarded quadratic form)", {
  # With the positive-definiteness guard, ANY input either abstains (typed warning,
  # all NA) or scores D^2 from a positive-definite quadratic form, >= 0 for every
  # row -- including missing-item rows on the zero-fill path. Sweeping missingness
  # fixtures pins the property the complete-data non-negativity test above cannot
  # reach (its pairwise cov is PSD by construction).
  for (seed in 1:8) {
    x <- withr::with_seed(seed, {
      m <- matrix(sample.int(5L, 40L * 6L, replace = TRUE), nrow = 40L)
      m[sample.int(length(m), 72L)] <- NA # 72 / 240 = 30% missing
      m
    })
    storage.mode(x) <- "double"
    v <- suppressWarnings(cier_mahalanobis(x))$value
    expect_true(all(v >= 0, na.rm = TRUE))
  }
})

test_that("the warning distinguishes the singular vs insufficient-respondents cause", {
  # The class is shared, but the structured payload names which case occurred, so
  # the two messages differ (asserted without matching message text).
  w_few <- tryCatch(cier_mahalanobis(matrix(c(1, 2, 3), nrow = 1L)),
                    warning = function(w) w)
  x_sing <- rand_matrix(n = 20L, p = 4L, seed = 8L)
  x_sing <- cbind(x_sing, x_sing[, 1L])
  w_sing <- tryCatch(cier_mahalanobis(x_sing), warning = function(w) w)
  expect_identical(cier_condition_data(w_few)$reason, "insufficient_responses")
  expect_identical(cier_condition_data(w_sing)$reason, "singular_covariance")
})

test_that("the degenerate-cause message names each cause", {
  # The shared-class warning relies on mahalanobis_abstain_message() to render a
  # cause-specific first line. Testing the helper directly catches a mutant that
  # hardcodes one cause (collapsing the if/else) -- the reason payload alone would
  # not, as it is sourced from the status, not the message.
  sing <- mahalanobis_abstain_message("singular_covariance")
  few <- mahalanobis_abstain_message("insufficient_responses")
  indef <- mahalanobis_abstain_message("indefinite_covariance")
  expect_match(sing[[1L]], "singular")
  expect_match(few[[1L]], "fewer than two")
  expect_match(indef[[1L]], "not positive definite")
  expect_identical(anyDuplicated(c(sing[[1L]], few[[1L]], indef[[1L]])), 0L)
})

# ---- Partial-row warning (scored rows carry imputed missing items) ----------

test_that("scoring partial rows warns once, naming the scored-partial count", {
  x <- partial_fixture()
  # One warning per call, not one per partial row: a per-row-loop mutant would
  # signal three times on this three-partial-row fixture.
  cls <- collect_warning_classes(cier_mahalanobis(x))
  expect_identical(sum(cls == "cier_warning_partial_rows"), 1L)
  w <- tryCatch(cier_mahalanobis(x), warning = function(w) w)
  # Three partial rows (3, 10, 20) of 30 scored respondents -- the count is SCORED
  # rows carrying an NA, not the total NA cell count or nrow.
  expect_identical(cier_condition_data(w)$n_partial, 3L)
  expect_identical(cier_condition_data(w)$n_scored, 30L)
  # The warning is purely additive: the statistic (oracle parity) is unchanged.
  expect_equal(suppressWarnings(cier_mahalanobis(x))$value,
               ref_mahalanobis(x), tolerance = 1e-10)
})

test_that("the partial-row count excludes all-NA abstaining rows", {
  # One partial (scored) row plus one all-NA (abstaining) row appended: scored
  # partial rows = 1, and n_scored excludes the abstainer (30 of 31), so a mutant
  # counting nrow() or raw NA-bearing rows reports one too many.
  x <- rand_matrix(n = 30L, p = 6L, seed = 7L)
  x[3L, 1L] <- NA
  x <- rbind(x, rep(NA_real_, 6L))
  w <- tryCatch(cier_mahalanobis(x), warning = function(w) w)
  expect_identical(cier_condition_data(w)$n_partial, 1L)
  expect_identical(cier_condition_data(w)$n_scored, 30L)
})

test_that("the partial-rows warning names the count and the sensitivity check", {
  # Wording is locked here (the typed-warning convention is an expect_match on the
  # rendered message, not a snapshot). cnd_message renders the info bullets too;
  # conditionMessage alone would drop them. Whitespace is collapsed so a cli line
  # wrap inside a matched phrase does not break the match.
  w <- tryCatch(cier_mahalanobis(partial_fixture()), warning = function(w) w)
  msg <- gsub("\\s+", " ", cli::ansi_strip(rlang::cnd_message(w)))
  expect_match(msg, "3 of 30 scored respondents")
  expect_match(msg, "mean-imputed")
  expect_match(msg, "inverse covariance")   # the over-flag mechanism bullet
  expect_match(msg, "complete cases")       # the sensitivity-check action
})

test_that("complete data scores silently (no partial-rows warning)", {
  expect_no_warning(cier_mahalanobis(rand_matrix(n = 30L, p = 6L, seed = 7L)))
})

test_that("a wholesale-abstaining input does not also warn about partial rows", {
  # The indefinite-covariance fixture is 40% MCAR -- nearly every row is partial --
  # but the kernel abstains wholesale (every value NA), so NO row is scored and the
  # partial-rows warning stays silent; only the singular-covariance warning fires. A
  # mutant reading rowSums(is.na(responses)) off the raw matrix, without the
  # scored-value mask, would wrongly fire here.
  cls <- collect_warning_classes(cier_mahalanobis(indefinite_fixture()))
  expect_true("cier_warning_singular_covariance" %in% cls)
  expect_false("cier_warning_partial_rows" %in% cls)
})

test_that("the partial-rows warning is independent of the cutoff override", {
  # The over-flagging is in the reference, not the threshold: partial rows warn at
  # the chi-square default, a moved alpha, or a literal D^2 cutoff alike.
  expect_warning(cier_mahalanobis(partial_fixture(), cutoff = 10),
                 class = "cier_warning_partial_rows")
  expect_warning(cier_mahalanobis(partial_fixture(), alpha = 0.05),
                 class = "cier_warning_partial_rows")
})

test_that("a non-matrix / non-numeric payload is a typed input error", {
  expect_error(cier_mahalanobis(1:10), class = "cier_error_input")
  expect_error(cier_mahalanobis(matrix(letters[1:6], nrow = 2L)),
               class = "cier_error_input")
})

# ---- Cutoff: chi-square default, alpha override, df = item count -------------

test_that("default cutoff is qchisq(1 - 0.001, df = p) with df the ITEM count", {
  x <- rand_matrix(n = 30L, p = 6L, seed = 7L)   # n = 30 != p = 6
  out <- cier_mahalanobis(x)
  # df is the item count (ncol), NOT the respondent count (nrow): a df = n mutant
  # would give qchisq(0.999, 30) ~ 59.7 instead of ~22.46.
  expect_identical(out$cutoff, as.numeric(stats::qchisq(1 - 0.001, df = 6)))
  expect_false(isTRUE(all.equal(out$cutoff,
                                as.numeric(stats::qchisq(1 - 0.001, df = 30)))))
})

test_that("the alpha argument moves the chi-square tail", {
  x <- rand_matrix(n = 30L, p = 6L, seed = 7L)
  out <- cier_mahalanobis(x, alpha = 0.01)
  expect_identical(out$cutoff, as.numeric(stats::qchisq(1 - 0.01, df = 6)))
})

test_that("an absolute cutoff overrides the chi-square default and flags via direction", {
  x <- rand_matrix(n = 30L, p = 6L, seed = 7L)
  out <- cier_mahalanobis(x, cutoff = 10)
  expect_identical(out$cutoff, 10)
  expect_identical(out$flagged, !is.na(out$value) & out$value >= 10)
})

test_that("a respondent exactly at the cutoff is flagged (>= ties, not >)", {
  # apply_flag uses >= for the upper tail, so a D^2 equal to the cutoff flags. Set
  # the literal cutoff to one respondent's own score; a strict-greater (>) mutant
  # would leave it unflagged.
  x <- rand_matrix(n = 30L, p = 6L, seed = 7L)
  v <- cier_mahalanobis(x)$value
  k <- which.max(v)
  out <- cier_mahalanobis(x, cutoff = v[[k]])
  expect_true(out$flagged[[k]])
})

test_that("supplying both alpha and cutoff is a typed input error", {
  expect_error(cier_mahalanobis(rand_matrix(), alpha = 0.01, cutoff = 10),
               class = "cier_error_input")
})

test_that("an invalid alpha is a typed input error, rejected at the open (0, 1) bounds", {
  x <- rand_matrix()
  expect_error(cier_mahalanobis(x, alpha = -0.1), class = "cier_error_input")
  expect_error(cier_mahalanobis(x, alpha = 2), class = "cier_error_input")
  expect_error(cier_mahalanobis(x, alpha = "x"), class = "cier_error_input")
  expect_error(cier_mahalanobis(x, alpha = c(0.01, 0.05)),
               class = "cier_error_input")
  # Boundary values rejected (open interval): alpha = 1 yields a zero cutoff that
  # flags everyone, alpha = 0 an infinite cutoff that flags no one. The documented
  # domain is (0, 1), so both abort.
  expect_error(cier_mahalanobis(x, alpha = 1), class = "cier_error_input")
  expect_error(cier_mahalanobis(x, alpha = 0), class = "cier_error_input")
})

test_that("an out-of-range alpha fails BEFORE the kernel runs (earliest point)", {
  # Validation is in the wrapper, so a bad alpha aborts even on input that would make
  # the kernel raise a (different) condition -- the alpha error wins.
  singular <- cbind(rand_matrix(n = 8L, p = 3L, seed = 2L),
                    rand_matrix(n = 8L, p = 3L, seed = 2L)[, 1L])
  expect_error(cier_mahalanobis(singular, alpha = 1.5),
               class = "cier_error_input")
})

test_that("an invalid absolute cutoff is a typed input error", {
  x <- rand_matrix()
  expect_error(cier_mahalanobis(x, cutoff = c(1, 2)), class = "cier_error_input")
  expect_error(cier_mahalanobis(x, cutoff = NA_real_), class = "cier_error_input")
  expect_error(cier_mahalanobis(x, cutoff = "x"), class = "cier_error_input")
  # D^2 is non-negative; a negative threshold is out of domain (would flag everyone),
  # so it is rejected rather than accepted silently.
  expect_error(cier_mahalanobis(x, cutoff = -1), class = "cier_error_input")
})

# ---- print snapshot (locked, design-first; direction = upper) ---------------

test_that("print renders the locked cli summary (upper direction)", {
  # A relaxed tail (alpha = 0.2) surfaces flags so the flagged-line render (the `>=`
  # comparator and a nonzero count) is what gets locked.
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    rnd <- rand_matrix(n = 29L, p = 6L, seed = 21L)
    x <- rbind(c(5, 1, 5, 1, 5, 1), rnd)
    expect_snapshot(print(cier_mahalanobis(x, alpha = 0.2)))
  })
})

test_that("print reports abstaining respondents on their own line", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    rnd <- rand_matrix(n = 28L, p = 6L, seed = 21L)
    x <- rbind(c(5, 1, 5, 1, 5, 1), rnd, rep(NA_real_, 6L))
    expect_snapshot(print(cier_mahalanobis(x, alpha = 0.2)))
  })
})
