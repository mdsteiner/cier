# Tests for cier_synonym_pairs() -- the pair-listing diagnostic that surfaces the
# whole-sample inter-item correlations behind psychsyn / psychant, so users can
# choose a critical_r.
#
# Trust model: an independent nested-loop reference (ref_pairs) re-derives the
# lower-triangle pair table + correlations by hand and never calls the production
# helper; the qualifying-pair SET is additionally pinned against the production
# find_item_pairs() (the same pairing the index scores) and, cross-package,
# against careless:::psychsyn_critval().

# Independent lower-triangle pair table (item_i > item_j), optionally filtered to
# a critical_r tail, with the same sort as the helper. Hand-rolled; no call to the
# production code.
ref_pairs <- function(x, critical_r = NULL, antonym = FALSE) {
  cm <- stats::cor(x, use = "pairwise.complete.obs")
  p <- ncol(cm)
  rows <- list()
  for (j in seq_len(p - 1L)) {
    for (i in seq.int(j + 1L, p)) {
      rr <- cm[i, j]
      if (is.na(rr)) next
      if (!is.null(critical_r)) {
        keep <- if (antonym) rr < -critical_r else rr > critical_r
        if (!keep) next
      }
      rows[[length(rows) + 1L]] <- c(i = i, j = j, r = rr)
    }
  }
  if (!length(rows)) {
    return(data.frame(item_i = integer(0L), item_j = integer(0L), r = numeric(0L)))
  }
  m <- do.call(rbind, rows)
  data.frame(item_i = as.integer(m[, "i"]), item_j = as.integer(m[, "j"]),
             r = m[, "r"])
}

pair_keys <- function(df) paste(df$item_i, df$item_j, sep = "-")

# Reuse the latent-factor fixture from the psychsyn tests (synonym pairs exist).
syn_matrix_d <- function(n = 60L, per = 3L, seed = 42L, noise = 0.35) {
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

# A fixture with MULTIPLE antonym pairs of differing strength: cols 1-3 load on a
# latent factor f (mutually positive = synonyms), cols 4-5 load on -f (each
# negatively correlated with 1-3 = antonyms, with several distinct r values), so
# the antonym tail has > 1 pair to pin the ascending sort and the min strongest_r.
ant_matrix_d <- function(n = 60L, seed = 5L) {
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

test_that("cier_synonym_pairs returns the pinned data.frame schema", {
  out <- cier_synonym_pairs(syn_matrix_d(), critical_r = 0.5)
  expect_s3_class(out, "data.frame")
  expect_identical(names(out), c("item_i", "item_j", "name_i", "name_j", "r"))
  expect_type(out$item_i, "integer")
  expect_type(out$item_j, "integer")
  expect_type(out$r, "double")
  expect_true(all(out$item_i > out$item_j))   # larger item index first
})

# ---- Independent reference + production-pairing consistency ------------------

test_that("the qualifying pairs and r values match an independent recomputation", {
  x <- syn_matrix_d(n = 50L, seed = 2026L)
  out <- cier_synonym_pairs(x, critical_r = 0.5)
  ref <- ref_pairs(x, critical_r = 0.5)
  expect_setequal(pair_keys(out), pair_keys(ref))
  # r matches the hand correlation for each pair (matched by key, not row order).
  ref_r <- setNames(ref$r, pair_keys(ref))
  expect_equal(out$r, unname(ref_r[pair_keys(out)]), tolerance = 1e-12)
  expect_true(all(out$r > 0.5))
})

test_that("the filtered pair set equals the production find_item_pairs() set", {
  # The lister must surface exactly the pairs the index scores at the same
  # critical_r -- otherwise the diagnostic would mislead the cutoff choice.
  x <- syn_matrix_d(n = 40L, seed = 7L)
  out <- cier_synonym_pairs(x, critical_r = 0.5)
  fip <- find_item_pairs(x, 0.5, "syn")
  expect_setequal(pair_keys(out),
                  paste(fip[, 1L], fip[, 2L], sep = "-"))
})

test_that("critical_r = NULL lists every non-NA lower-triangle pair", {
  x <- syn_matrix_d(n = 30L)
  out <- cier_synonym_pairs(x, critical_r = NULL)
  p <- ncol(x)
  expect_identical(nrow(out), as.integer(p * (p - 1L) / 2L))
  expect_setequal(pair_keys(out), pair_keys(ref_pairs(x, critical_r = NULL)))
})

test_that("the critical_r threshold is strict (> not >=) at an exact pair r", {
  # Continuous random fixtures never land a pair exactly on the cutoff, so an
  # off-by-one (>= vs >) survives expect_setequal. Set critical_r to the exact
  # correlation of a real pair: a strict > excludes that pair; a >= mutant keeps
  # it. Pinned against the production find_item_pairs at the same boundary.
  x <- syn_matrix_d(n = 40L, seed = 3L)
  cm <- stats::cor(x, use = "pairwise.complete.obs")
  cm[upper.tri(cm, diag = TRUE)] <- NA
  cut <- sort(cm[!is.na(cm)], decreasing = TRUE)[[3L]]   # exact r of the 3rd pair
  out <- cier_synonym_pairs(x, critical_r = cut)
  expect_true(all(out$r > cut))
  expect_identical(nrow(out), sum(cm > cut, na.rm = TRUE))
  fip <- find_item_pairs(x, cut, "syn")
  expect_setequal(pair_keys(out), paste(fip[, 1L], fip[, 2L], sep = "-"))
})

test_that("a constant (zero-variance) column forms no pair and leaves r finite", {
  # A constant item gives NA correlations; it must be dropped from the listing and
  # never leak an NA r (and, in the sweep, never poison strongest_r).
  x <- syn_matrix_d(n = 40L, seed = 6L)
  x[, 5L] <- 3
  out <- cier_synonym_pairs(x, critical_r = NULL)
  expect_false(any(out$item_i == 5L | out$item_j == 5L))
  expect_true(all(is.finite(out$r)))
})

# ---- Sort + direction -------------------------------------------------------

test_that("synonym output is sorted by descending r", {
  out <- cier_synonym_pairs(syn_matrix_d(n = 40L), critical_r = NULL)
  expect_false(is.unsorted(rev(out$r)))     # non-increasing
})

test_that("antonym = TRUE surfaces the negative tail, sorted ascending", {
  # Two strongly negatively correlated items among orthogonal others.
  a <- c(4, 4, 2, 2)
  neg <- c(2, 2, 4, 4)
  b <- c(4, 2, 4, 2)
  c3 <- c(4, 2, 2, 4)
  x <- cbind(a, neg, b, c3)
  dimnames(x) <- NULL
  storage.mode(x) <- "double"
  out <- cier_synonym_pairs(x, critical_r = 0.6, antonym = TRUE)
  expect_identical(nrow(out), 1L)
  expect_identical(c(out$item_i, out$item_j), c(2L, 1L))
  expect_true(out$r < -0.6)
  # Consistency with the production antonym pairing.
  fip <- find_item_pairs(x, 0.6, "ant")
  expect_setequal(pair_keys(out), paste(fip[, 1L], fip[, 2L], sep = "-"))
})

test_that("antonym output with several pairs is sorted by ASCending r (most negative first)", {
  # The single-pair Hadamard fixture cannot test sort direction; this one has
  # multiple antonym pairs of differing strength. A descending / no-sort mutant on
  # the negative tail survives a 1-row frame but fails here.
  x <- ant_matrix_d()
  out <- cier_synonym_pairs(x, critical_r = 0.6, antonym = TRUE)
  expect_gt(nrow(out), 1L)
  expect_false(is.unsorted(out$r))            # non-decreasing == ascending
  expect_true(all(out$r < -0.6))
  ref <- ref_pairs(x, critical_r = 0.6, antonym = TRUE)
  expect_setequal(pair_keys(out), pair_keys(ref))
})

test_that("critical_r = NULL with antonym = TRUE lists all pairs, ascending", {
  x <- ant_matrix_d()
  out <- cier_synonym_pairs(x, critical_r = NULL, antonym = TRUE)
  p <- ncol(x)
  expect_identical(nrow(out), as.integer(p * (p - 1L) / 2L))
  expect_false(is.unsorted(out$r))
})

# ---- Names ------------------------------------------------------------------

test_that("item names are carried when present and NA when absent", {
  x <- syn_matrix_d(n = 30L)
  colnames(x) <- paste0("Q", seq_len(ncol(x)))
  out <- cier_synonym_pairs(x, critical_r = 0.5)
  expect_identical(out$name_i, paste0("Q", out$item_i))
  expect_identical(out$name_j, paste0("Q", out$item_j))
  x2 <- syn_matrix_d(n = 30L)                # no colnames
  out2 <- cier_synonym_pairs(x2, critical_r = 0.5)
  expect_true(all(is.na(out2$name_i)))
  expect_true(all(is.na(out2$name_j)))
})

# ---- Cross-package: careless::psychsyn_critval ------------------------------

# Map careless's name-keyed pairs (var1, var2 factors) to "i-j" keys on the
# column indices, larger index first, to compare against cier's item_i / item_j.
careless_pair_keys <- function(cv, item_names) {
  i <- match(as.character(cv$var1), item_names)
  j <- match(as.character(cv$var2), item_names)
  paste(pmax(i, j), pmin(i, j), sep = "-")
}

test_that("cier_synonym_pairs matches careless::psychsyn_critval (synonyms)", {
  # GENUINE cross-package parity: careless::psychsyn_critval() returns the full
  # p x p long frame (var1, var2, cor) with the upper triangle + diagonal NA; its
  # non-NA rows are the lower-triangle pairs sorted by descending correlation --
  # exactly what cier_synonym_pairs(critical_r = NULL) produces. Same pairs, same
  # correlations, same order.
  skip_if_not_installed("careless")
  raw <- careless::careless_dataset
  responses <- unname(as.matrix(raw))
  storage.mode(responses) <- "double"
  cv <- careless::psychsyn_critval(raw, anto = FALSE)
  cv <- cv[!is.na(cv$cor), ]
  ours <- cier_synonym_pairs(responses, critical_r = NULL)
  expect_identical(nrow(ours), nrow(cv))
  expect_equal(ours$r, cv$cor, tolerance = 1e-12)          # values, in careless's order
  expect_setequal(pair_keys(ours), careless_pair_keys(cv, colnames(raw)))
})

test_that("cier_synonym_pairs matches careless::psychsyn_critval (antonyms)", {
  # The antonym tail: careless sorts anto = TRUE ascending (most negative first),
  # which cier reproduces with antonym = TRUE.
  skip_if_not_installed("careless")
  raw <- careless::careless_dataset
  responses <- unname(as.matrix(raw))
  storage.mode(responses) <- "double"
  cv <- careless::psychsyn_critval(raw, anto = TRUE)
  cv <- cv[!is.na(cv$cor), ]
  ours <- cier_synonym_pairs(responses, critical_r = NULL, antonym = TRUE)
  expect_identical(nrow(ours), nrow(cv))
  expect_equal(ours$r, cv$cor, tolerance = 1e-12)
  expect_setequal(pair_keys(ours), careless_pair_keys(cv, colnames(raw)))
})

# ---- Edges + validation -----------------------------------------------------

test_that("a critical_r that no pair clears yields a zero-row frame, not an error", {
  # careless::get_item_pairs stop()s here; the cier lister returns an empty frame.
  out <- cier_synonym_pairs(syn_matrix_d(n = 30L), critical_r = 0.99)
  expect_s3_class(out, "data.frame")
  expect_identical(nrow(out), 0L)
  expect_identical(names(out), c("item_i", "item_j", "name_i", "name_j", "r"))
})

test_that("invalid inputs are typed errors", {
  x <- syn_matrix_d(n = 10L)
  expect_error(cier_synonym_pairs(1:10), class = "cier_error_input")
  expect_error(cier_synonym_pairs(x, critical_r = 0), class = "cier_error_input")
  expect_error(cier_synonym_pairs(x, critical_r = 1), class = "cier_error_input")
  expect_error(cier_synonym_pairs(x, critical_r = c(0.4, 0.6)),
               class = "cier_error_input")
  expect_error(cier_synonym_pairs(x, antonym = NA), class = "cier_error_input")
  bad <- x
  bad[1L, 1L] <- Inf                       # check_responses must reject NaN/Inf
  expect_error(cier_synonym_pairs(bad), class = "cier_error_input")
})
