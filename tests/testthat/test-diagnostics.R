# Cross-index agreement diagnostic (flag_agreement). The Poisson-binomial
# independence baseline is cross-checked against an independent 2^m-enumeration
# oracle; see tests/reference/ref-poisson-binomial-enumeration.R.

source(test_path("..", "reference", "ref-poisson-binomial-enumeration.R"))

# A logical n x m flag matrix whose columns fire at distinct marginal rates, so the
# Poisson-binomial baseline (unequal p) is genuinely exercised.
make_flags <- function(n, m, seed) {
  withr::with_seed(seed, {
    p <- stats::runif(m, 0.05, 0.5)
    do.call(cbind, lapply(p, function(pi) stats::runif(n) < pi))
  })
}

test_that("agreement$expected matches the independent enumeration oracle (1e-12)", {
  for (seed in c(1L, 7L, 42L)) {
    for (m in c(2L, 5L, 9L)) {
      flags <- make_flags(200L, m, seed + m)
      res <- flag_agreement(flags)
      expect_identical(res$agreement$k, seq_len(m))
      expect_equal(res$agreement$expected,
                   ref_poisson_binomial_tail(colMeans(flags)),
                   tolerance = 1e-12)
    }
  }
})

test_that("agreement$observed equals the hand-counted share flagged by >= k votes", {
  flags <- make_flags(150L, 6L, 3L)
  counts <- rowSums(flags)
  obs <- vapply(seq_len(6L), function(k) mean(counts >= k), numeric(1L))
  expect_equal(flag_agreement(flags)$agreement$observed, obs)
})

test_that("per_vote reports observed rates; null_rate drives excess + informative", {
  flags <- make_flags(100L, 3L, 9L)
  p <- as.numeric(colMeans(flags))

  # No null supplied: every column is tautological (an empirical percentile flags
  # fpr by construction), so nothing is informative.
  pv0 <- flag_agreement(flags)$per_vote
  expect_equal(pv0$observed, p)
  expect_identical(pv0$null, rep(NA_real_, 3L))
  expect_identical(pv0$excess, rep(NA_real_, 3L))
  expect_identical(pv0$informative, rep(FALSE, 3L))

  # The vote identifier falls back to vote1..voteM when the matrix is unnamed, and
  # carries the column names through when present.
  expect_identical(pv0$vote, paste0("vote", 1:3))
  named <- flags
  colnames(named) <- c("cier_irv", "cier_even_odd", "cier_mahalanobis")
  expect_identical(flag_agreement(named)$per_vote$vote, colnames(named))

  # A null nominal on the third (null-referenced) vote only: its excess over the
  # calibrated null is informative; the percentile columns stay NA.
  nr <- c(NA_real_, NA_real_, 0.001)
  pv1 <- flag_agreement(flags, null_rate = nr)$per_vote
  expect_identical(pv1$informative, c(FALSE, FALSE, TRUE))
  expect_equal(pv1$excess[[3L]], p[[3L]] - 0.001)
  expect_identical(pv1$excess[1:2], rep(NA_real_, 2L))

  # `informative` means "a calibrated null was SUPPLIED", not "observed exceeds it":
  # a vote observed BELOW its null stays informative with a negative excess. Pins the
  # documented contract so a sign-gating refactor (informative = excess > 0) cannot
  # pass silently.
  nr_high <- c(NA_real_, NA_real_, 0.9)
  pv2 <- flag_agreement(flags, null_rate = nr_high)$per_vote
  expect_identical(pv2$informative, c(FALSE, FALSE, TRUE))
  expect_lt(pv2$excess[[3L]], 0)
})

test_that("null_rate of the wrong length is a typed input error", {
  flags <- make_flags(40L, 4L, 5L)
  expect_error(flag_agreement(flags, null_rate = c(0.001, 0.05)),
               class = "cier_error_input")
})

test_that("agreement tails are monotone non-increasing in k", {
  res <- flag_agreement(make_flags(120L, 7L, 11L))$agreement
  expect_true(all(diff(res$observed) <= 0))
  expect_true(all(diff(res$expected) <= 0))
})

test_that("degenerate flag matrices behave sensibly", {
  n <- 50L
  m <- 4L

  res_none <- flag_agreement(matrix(FALSE, n, m))$agreement
  expect_equal(res_none$observed, rep(0, m))
  expect_equal(res_none$expected, rep(0, m))

  res_all <- flag_agreement(matrix(TRUE, n, m))$agreement
  expect_equal(res_all$observed, rep(1, m))
  expect_equal(res_all$expected, rep(1, m))

  # Single vote: no co-occurrence is possible, so observed == expected == p.
  one <- matrix(c(rep(TRUE, 10L), rep(FALSE, 40L)), ncol = 1L)
  r1 <- flag_agreement(one)$agreement
  expect_identical(r1$k, 1L)
  expect_equal(r1$observed, 0.2)
  expect_equal(r1$expected, 0.2)
})

test_that("flag_agreement rejects non-logical, NA, or empty flags", {
  expect_error(flag_agreement(matrix(c(1, 0, 1, 1), 2L, 2L)),
               class = "cier_error_input")
  expect_error(flag_agreement(matrix(c(TRUE, NA, FALSE, TRUE), 2L, 2L)),
               class = "cier_error_input")
  expect_error(flag_agreement(matrix(logical(0), 0L, 0L)),
               class = "cier_error_input")
  expect_error(flag_agreement(matrix(logical(0), 5L, 0L)),
               class = "cier_error_input")
})
