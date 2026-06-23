# Tests for the attention-check injection -- the direct half of cier_simulate().
# Checks are generated OUTSIDE the response matrix, shaped exactly for the shipped
# cier_attention(checks, pass) contract; failure is probabilistic on BOTH sides
# (careless rows fail at p_fail_careless = 0.75, attentive at p_fail_attentive =
# 0.05 -- documented conventions, not estimates), so neither group separates
# perfectly.
#
# Trust model: oracle-only. The oracle (ref-sim-times.R) re-derives the injection by
# hand with a coordinated draw order; the deterministic p_fail edges are pinned
# THROUGH the shipped cier_attention.

source(test_path("..", "reference", "ref-sim-times.R"))

# =============================================================================
# FAST TIER
# =============================================================================

# ---- Shape + membership contract --------------------------------------------

test_that("checks and pass are shaped for cier_attention, values in 1..5", {
  careless <- rep(c(TRUE, FALSE), c(10L, 30L))
  out <- withr::with_seed(8L, sim_direct_checks(40L, careless, 3L))
  expect_true(is.matrix(out$checks))
  expect_identical(dim(out$checks), c(40L, 3L))
  expect_true(is.integer(out$checks))
  expect_true(all(out$checks >= 1L & out$checks <= 5L))
  expect_true(is.list(out$pass))
  expect_length(out$pass, 3L)
  for (j in seq_len(3L)) {
    expect_true(is.integer(out$pass[[j]]) && length(out$pass[[j]]) == 1L)
    expect_true(out$pass[[j]] >= 1L && out$pass[[j]] <= 5L)
  }
  # consumable by the shipped wrapper as-is.
  att <- cier_attention(out$checks, out$pass)
  expect_length(att$value, 40L)
})

test_that("n_checks = 0 returns NULL slots and consumes no RNG", {
  out <- sim_direct_checks(5L, rep(FALSE, 5L), 0L)
  # the full named-list shape, not a bare NULL (whose $checks is also NULL).
  expect_identical(out, list(checks = NULL, pass = NULL))
  a <- withr::with_seed(1L, stats::runif(1L))
  b <- withr::with_seed(1L, {
    sim_direct_checks(5L, rep(FALSE, 5L), 0L)
    stats::runif(1L)
  })
  expect_identical(a, b)
})

# ---- Oracle parity (coordinated draw order) ---------------------------------

test_that("the injection matches the independent hand re-derivation", {
  careless <- c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE)
  prod <- withr::with_seed(20260612L, sim_direct_checks(8L, careless, 4L))
  ref <- ref_direct_checks(20260612L, 8L, careless, 4L)
  expect_identical(prod$checks, ref$checks)
  expect_identical(prod$pass, ref$pass)
  prod2 <- withr::with_seed(20260613L, sim_direct_checks(8L, careless, 4L))
  expect_false(identical(prod$checks, prod2$checks))      # RNG genuinely used
})

test_that("oracle parity holds under overridden failure probabilities", {
  careless <- rep(c(TRUE, FALSE), 6L)
  prod <- withr::with_seed(
    55L, sim_direct_checks(12L, careless, 2L,
                           p_fail_careless = 0.4, p_fail_attentive = 0.2)
  )
  ref <- ref_direct_checks(55L, 12L, careless, 2L,
                           p_fail_careless = 0.4, p_fail_attentive = 0.2)
  expect_identical(prod$checks, ref$checks)
  expect_identical(prod$pass, ref$pass)
})

# ---- Deterministic p_fail edges, pinned through cier_attention --------------

test_that("p_fail 1 / 0 splits the groups exactly through cier_attention", {
  careless <- rep(c(TRUE, FALSE), c(7L, 13L))
  out <- withr::with_seed(
    21L, sim_direct_checks(20L, careless, 3L,
                           p_fail_careless = 1, p_fail_attentive = 0)
  )
  value <- cier_attention(out$checks, out$pass)$value
  expect_true(all(value[careless] == 3))                  # every check failed
  expect_true(all(value[!careless] == 0))                 # every check passed
  # failure semantics: a failing response is never the key, a passing one is.
  for (j in seq_len(3L)) {
    expect_true(all(out$checks[careless, j] != out$pass[[j]]))
    expect_true(all(out$checks[!careless, j] == out$pass[[j]]))
  }
})

test_that("flipped edge: attentive can be forced to fail and careless to pass", {
  careless <- rep(c(TRUE, FALSE), c(4L, 6L))
  out <- withr::with_seed(
    22L, sim_direct_checks(10L, careless, 2L,
                           p_fail_careless = 0, p_fail_attentive = 1)
  )
  value <- cier_attention(out$checks, out$pass)$value
  expect_true(all(value[careless] == 0))
  expect_true(all(value[!careless] == 2))
})

# ---- Input validation --------------------------------------------------------

test_that("the injection rejects malformed inputs", {
  ok_careless <- rep(FALSE, 5L)
  expect_error(sim_direct_checks(5L, ok_careless, -1L), class = "cier_error_input")
  expect_error(sim_direct_checks(5L, ok_careless, 1.5), class = "cier_error_input")
  expect_error(sim_direct_checks(5L, ok_careless, NA), class = "cier_error_input")
  expect_error(sim_direct_checks(5L, c(1, 0, 1, 0, 1), 2L),  # not logical
               class = "cier_error_input")
  expect_error(sim_direct_checks(5L, rep(FALSE, 4L), 2L),    # wrong length
               class = "cier_error_input")
  expect_error(sim_direct_checks(5L, c(TRUE, NA, FALSE, FALSE, FALSE), 2L),
               class = "cier_error_input")                   # NA carelessness
  expect_error(sim_direct_checks(5L, ok_careless, 2L, p_fail_careless = 1.2),
               class = "cier_error_input")
  expect_error(sim_direct_checks(5L, ok_careless, 2L, p_fail_attentive = -0.1),
               class = "cier_error_input")
})

# =============================================================================
# SLOW TIER -- large-n failure-rate contract
# =============================================================================

test_that("large-n failure rates land on the documented conventions", {
  skip_on_cran()
  skip_if_slow()
  careless <- rep(c(TRUE, FALSE), each = 2000L)
  out <- withr::with_seed(20260612L, sim_direct_checks(4000L, careless, 4L))
  key_mat <- matrix(unlist(out$pass), nrow = 4000L, ncol = 4L, byrow = TRUE)
  failed <- out$checks != key_mat
  rate_careless <- mean(failed[careless, ])
  rate_attentive <- mean(failed[!careless, ])
  expect_lt(abs(rate_careless - 0.75), 0.03)
  expect_lt(abs(rate_attentive - 0.05), 0.02)
})
