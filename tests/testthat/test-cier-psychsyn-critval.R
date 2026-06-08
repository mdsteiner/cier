# Tests for cier_psychsyn_critval() -- the critical_r sweep diagnostic that shows,
# per candidate critical_r, how many synonym (or antonym) pairs qualify and how
# many respondents become scorable, so users can weigh the cutoff against coverage.
#
# Trust model: n_pairs is pinned against the production find_item_pairs() and
# n_scored against the INDEPENDENT oracle ref_psychsyn (a per-row cor loop, not the
# production kernel); the strongest_r and abstention arithmetic are recomputed by
# hand. The reproduction of the documented BFI numbers is an end-to-end anchor.

source(test_path("..", "reference", "ref-psychsyn-meade-craig-2012.R"))

syn_matrix_c <- function(n = 60L, per = 3L, seed = 42L, noise = 0.35) {
  withr::with_seed(seed, {
    mk <- function() {
      f <- stats::rnorm(n)
      vapply(seq_len(per), function(j) f + stats::rnorm(n, 0, noise), numeric(n))
    }
    x <- cbind(mk(), mk(), mk())
  })
  storage.mode(x) <- "double"
  x
}

# Several antonym pairs of differing strength (cols 1-3 on f, cols 4-5 on -f).
ant_matrix_c <- function(n = 60L, seed = 5L) {
  withr::with_seed(seed, {
    f <- stats::rnorm(n)
    pos <- vapply(seq_len(3L), function(j) f + stats::rnorm(n, 0, 0.30), numeric(n))
    neg <- vapply(seq_len(2L), function(j) -f + stats::rnorm(n, 0, 0.30), numeric(n))
    x <- cbind(pos, neg)
  })
  storage.mode(x) <- "double"
  x
}

# ---- Schema -----------------------------------------------------------------

test_that("cier_psychsyn_critval returns the pinned sweep schema", {
  out <- cier_psychsyn_critval(syn_matrix_c(), critical_r = c(0.7, 0.5, 0.3))
  expect_s3_class(out, "data.frame")
  expect_identical(names(out),
                   c("critical_r", "n_pairs", "n_scored", "n_abstain",
                     "abstention_rate", "strongest_r"))
  expect_identical(nrow(out), 3L)
  expect_identical(out$critical_r, c(0.7, 0.5, 0.3))
  expect_type(out$n_pairs, "integer")
  expect_type(out$n_scored, "integer")
})

# ---- Pinned against production pairing + independent oracle counts -----------

test_that("n_pairs matches find_item_pairs and n_scored matches the oracle", {
  x <- syn_matrix_c(n = 50L, seed = 99L)
  grid <- c(0.8, 0.6, 0.4)
  out <- cier_psychsyn_critval(x, critical_r = grid)
  for (g in seq_along(grid)) {
    cr <- grid[[g]]
    expect_identical(out$n_pairs[[g]], nrow(find_item_pairs(x, cr, "syn")))
    # Independent oracle: a per-row cor loop, NOT the production kernel.
    expect_identical(out$n_scored[[g]], sum(!is.na(ref_psychsyn(x, cr))))
  }
  expect_identical(out$n_abstain, nrow(x) - out$n_scored)
  expect_equal(out$abstention_rate, round(out$n_abstain / nrow(x), 3),
               tolerance = 1e-9)
})

test_that("strongest_r is the strongest in-tail inter-item correlation, rounded", {
  x <- syn_matrix_c(n = 40L, seed = 7L)
  cm <- stats::cor(x, use = "pairwise.complete.obs")
  cm[upper.tri(cm, diag = TRUE)] <- NA
  out <- cier_psychsyn_critval(x, critical_r = c(0.6, 0.4))
  expect_true(all(out$strongest_r == round(max(cm, na.rm = TRUE), 3)))
})

test_that("n_pairs is non-decreasing as critical_r falls (more pairs at a looser cutoff)", {
  # The grid runs from the strictest (0.8, fewest pairs) to the loosest (0.3, most),
  # so n_pairs along the rows is non-decreasing.
  out <- cier_psychsyn_critval(syn_matrix_c(n = 50L),
                               critical_r = c(0.8, 0.7, 0.6, 0.5, 0.4, 0.3))
  expect_false(is.unsorted(out$n_pairs))
})

# ---- End-to-end BFI anchor (the user's reported case) -----------------------

test_that("the sweep reproduces the documented BFI coverage numbers", {
  bfi <- as.matrix(bfi_careless[, 1:44])
  storage.mode(bfi) <- "double"
  out <- cier_psychsyn_critval(bfi, critical_r = c(0.60, 0.50))
  # The BFI tops out at r = 0.59, so the default 0.60 finds NO pairs (the user's
  # bug); lowering to 0.50 surfaces 23 pairs and scores 370 of 394 respondents.
  expect_identical(out$n_pairs, c(0L, 23L))
  expect_identical(out$n_scored, c(0L, 370L))
  expect_true(all(out$strongest_r == 0.590))
})

# ---- Antonym tail -----------------------------------------------------------

test_that("antonym = TRUE sweeps the negative tail", {
  a <- c(4, 4, 2, 2)
  neg <- c(2, 2, 4, 4)
  b <- c(4, 2, 4, 2)
  c3 <- c(4, 2, 2, 4)
  x <- cbind(a, neg, b, c3)
  dimnames(x) <- NULL
  storage.mode(x) <- "double"
  out <- cier_psychsyn_critval(x, critical_r = 0.6, antonym = TRUE)
  expect_identical(out$n_pairs, nrow(find_item_pairs(x, 0.6, "ant")))
  cm <- stats::cor(x, use = "pairwise.complete.obs")
  cm[upper.tri(cm, diag = TRUE)] <- NA
  expect_true(all(out$strongest_r == round(min(cm, na.rm = TRUE), 3)))
})

test_that("antonym strongest_r is the MOST NEGATIVE when several antonym pairs exist", {
  # The single-pair Hadamard fixture cannot distinguish min from max on the
  # negative tail; this one has several. strongest_r must be the global minimum.
  x <- ant_matrix_c()
  out <- cier_psychsyn_critval(x, critical_r = c(0.7, 0.6), antonym = TRUE)
  cm <- stats::cor(x, use = "pairwise.complete.obs")
  cm[upper.tri(cm, diag = TRUE)] <- NA
  expect_true(all(out$strongest_r == round(min(cm, na.rm = TRUE), 3)))
  expect_gt(out$n_pairs[[2L]], 1L)
})

test_that("an all-constant battery yields NA strongest_r and no pairs", {
  # No finite inter-item correlation exists, so there is nothing to be strongest:
  # strongest_r is NA (not -Inf from an empty max) and every threshold finds 0.
  x <- matrix(3, nrow = 20L, ncol = 5L)
  storage.mode(x) <- "double"
  out <- suppressMessages(cier_psychsyn_critval(x, critical_r = c(0.5, 0.4)))
  expect_true(all(out$n_pairs == 0L))
  expect_true(all(is.na(out$strongest_r)))
})

test_that("a constant column does not poison strongest_r (na.rm)", {
  # A constant item gives NA correlations; max()/min() without na.rm would make
  # strongest_r NA. It must reflect the strongest FINITE inter-item correlation.
  x <- syn_matrix_c(n = 40L, seed = 6L)
  x[, 5L] <- 3
  out <- cier_psychsyn_critval(x, critical_r = 0.5)
  cm <- suppressWarnings(stats::cor(x, use = "pairwise.complete.obs"))
  cm[upper.tri(cm, diag = TRUE)] <- NA
  expect_true(all(is.finite(out$strongest_r)))
  expect_true(all(out$strongest_r == round(max(cm, na.rm = TRUE), 3)))
})

# ---- Guidance message -------------------------------------------------------

test_that("a typed message fires when the highest grid value finds no pairs", {
  bfi <- as.matrix(bfi_careless[, 1:44])
  storage.mode(bfi) <- "double"
  expect_message(cier_psychsyn_critval(bfi, critical_r = c(0.60, 0.50)),
                 class = "cier_message_no_pairs")
})

test_that("the message keys on max(critical_r), not on grid position", {
  # Unsorted grid: the largest threshold (0.60) sits in the MIDDLE and finds 0
  # pairs. A mutant keying on the first/last row (0.50 / 0.45, both find pairs)
  # would not fire; the correct condition (the strictest threshold tested) does.
  bfi <- as.matrix(bfi_careless[, 1:44])
  storage.mode(bfi) <- "double"
  expect_message(cier_psychsyn_critval(bfi, critical_r = c(0.50, 0.60, 0.45)),
                 class = "cier_message_no_pairs")
})

test_that("no message fires when every grid value finds pairs", {
  x <- syn_matrix_c(n = 40L)
  expect_no_message(cier_psychsyn_critval(x, critical_r = c(0.5, 0.4)))
})

test_that("the guidance message renders the locked text", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    bfi <- as.matrix(bfi_careless[, 1:44])
    storage.mode(bfi) <- "double"
    expect_snapshot(invisible(cier_psychsyn_critval(bfi, critical_r = c(0.60, 0.50))))
  })
})

# ---- Validation -------------------------------------------------------------

test_that("invalid inputs are typed errors", {
  x <- syn_matrix_c(n = 10L)
  expect_error(cier_psychsyn_critval(1:10), class = "cier_error_input")
  expect_error(cier_psychsyn_critval(x, critical_r = 0), class = "cier_error_input")
  expect_error(cier_psychsyn_critval(x, critical_r = c(0.5, 1.2)),
               class = "cier_error_input")
  expect_error(cier_psychsyn_critval(x, critical_r = numeric(0L)),
               class = "cier_error_input")
  expect_error(cier_psychsyn_critval(x, critical_r = "x"),
               class = "cier_error_input")
  expect_error(cier_psychsyn_critval(x, antonym = NA), class = "cier_error_input")
  bad <- x
  bad[1L, 1L] <- Inf                       # check_responses must reject NaN/Inf
  expect_error(cier_psychsyn_critval(bad), class = "cier_error_input")
})
