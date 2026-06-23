# Tests for cier_flagged_cases() -- the accessor turning a cier_index or cier_screen
# into the row positions of the flagged respondents.
#
# A pure extractor (no statistics), so the contract is: the cier_index method
# returns which(flagged) (NA-safe, 1-based, ascending); the cier_screen method
# thresholds the COLLAPSED votes (rowSums(votes) >= min_votes), so it agrees exactly
# with the agreement diagnostic (length == the ">= k votes" count) and never
# double-counts the consistency construct. The decisive test is the collapse one: a
# "counts raw per-index flags" implementation flags a respondent hit by both even-odd
# and PR at min_votes = 2, whereas the collapsed count is 1.

flag_matrix <- function(n = 60L, p = 20L, seed = 1L) {
  withr::with_seed(seed, {
    m <- matrix(sample.int(5L, n * p, replace = TRUE), nrow = n, ncol = p)
  })
  storage.mode(m) <- "double"
  m
}

flag_items <- function(n_scales = 10L, per_scale = 3L, max = 5L) {
  scale <- rep(paste0("s", seq_len(n_scales)), each = per_scale)
  data.frame(scale = scale,
             reverse_keyed = rep(c(FALSE, TRUE), length.out = length(scale)),
             max = max, stringsAsFactors = FALSE)
}

base8 <- function() {
  c("cier_longstring", "cier_irv", "cier_even_odd", "cier_personal_reliability",
    "cier_psychsyn", "cier_psychant", "cier_mahalanobis", "cier_person_total")
}

# ---- cier_index method ------------------------------------------------------

test_that("the cier_index method returns flagged positions, NA-safe", {
  # value[2] = NA forces flagged[2] = NA (new_cier_index's abstention rule), so an
  # abstaining respondent is never returned; FALSE rows drop too.
  idx <- new_cier_index(value = c(1, NA, 3, 4),
                        flagged = c(TRUE, NA, FALSE, TRUE),
                        method = "cier_irv", cutoff = 2, direction = "lower")
  out <- cier_flagged_cases(idx)
  expect_identical(out, c(1L, 4L))
  expect_type(out, "integer")
  # min_votes is a cier_screen-only argument; passing it to a single index is a
  # likely mistake, so it WARNS (typed) and is ignored -- the index returns every
  # flagged respondent regardless.
  expect_warning(out2 <- cier_flagged_cases(idx, min_votes = 99),
                 class = "cier_warning_ignored_min_votes")
  expect_identical(out2, c(1L, 4L))
})

test_that("none flagged -> integer(0); all flagged -> every position", {
  none <- new_cier_index(c(1, 2, 3), c(FALSE, FALSE, FALSE),
                         "cier_irv", 0, "lower")
  expect_identical(cier_flagged_cases(none), integer(0))
  all_f <- new_cier_index(c(1, 2, 3), c(TRUE, TRUE, TRUE),
                          "cier_irv", 5, "lower")
  expect_identical(cier_flagged_cases(all_f), c(1L, 2L, 3L))
})

test_that("the cier_index method equals which(out$flagged) end-to-end", {
  x <- flag_matrix(60L, 20L, 2L)
  x[5L, ] <- NA                       # an abstaining row (irv -> NA)
  out <- cier_irv(x)
  fc <- cier_flagged_cases(out)
  expect_identical(fc, which(out$flagged))
  expect_false(5L %in% fc)
})

# ---- cier_screen method: agreement consistency ------------------------------

test_that("flagged-case count matches the agreement table, on real multi-vote data", {
  # bfi_careless is contaminated real data, so its careless subgroup is flagged by
  # several constructs -- giving respondents with 0, 1, 2 and 3 votes. That spread
  # makes the per-level pins bite: which(>= 1), which(>= 2), which(>= 3) are all
  # different non-empty sets, so an "ignores min_votes" mutant (returns the >= 1 set
  # for every threshold) and an "always empty above 1" mutant both fail.
  nm <- names(bfi_careless)[1:44]
  items <- data.frame(scale = gsub("^v_BFI_|[0-9_R]+$", "", nm),
                      reverse_keyed = grepl("_R$", nm), max = 5L)
  sc <- cier_screen(bfi_careless[, 1:44], items, methods = base8(),
                    control = list(cier_personal_reliability = list(seed = 1)))
  n <- sc$n_respondents
  ag <- sc$agreement$agreement
  for (k in ag$k) {
    fc <- cier_flagged_cases(sc, min_votes = k)
    # Ties the extractor to the diagnostic: both count the COLLAPSED votes.
    expect_identical(length(fc), as.integer(round(ag$observed[[k]] * n)))
    # And independently: the exact rowSums(votes) >= k set (membership + order).
    expect_identical(fc, which(rowSums(as.matrix(sc$votes)) >= k))
  }
  expect_type(cier_flagged_cases(sc), "integer")
  # Positive control + strict shrinkage: respondents DO reach >= 2 votes, and the
  # set strictly shrinks as the threshold rises. The decisive guard against an
  # "always integer(0) above 1" or "min_votes-ignoring" implementation.
  expect_gt(length(cier_flagged_cases(sc, min_votes = 2)), 0L)
  expect_gt(length(cier_flagged_cases(sc, min_votes = 1)),
            length(cier_flagged_cases(sc, min_votes = 2)))
  expect_gt(length(cier_flagged_cases(sc, min_votes = 2)),
            length(cier_flagged_cases(sc, min_votes = 3)))
})

# ---- cier_screen method: the collapse is respected (votes, not raw flags) ----

test_that("min_votes counts collapsed votes, not raw per-index flags", {
  # Force BOTH consistency members to flag every scored respondent (cutoff = -1).
  # Each then carries TWO raw flags but only ONE construct vote, so min_votes = 2
  # must return nobody. A raw-flag implementation returns them all.
  x <- flag_matrix(80L, 20L, 23L)
  it <- flag_items(4L, 5L)
  sc <- cier_screen(
    x, it, methods = c("cier_even_odd", "cier_personal_reliability"),
    control = list(cier_even_odd = list(cutoff = -1),
                   cier_personal_reliability = list(cutoff = -1, seed = 1))
  )
  scored <- !is.na(sc$indices$cier_even_odd$value) &
    !is.na(sc$indices$cier_personal_reliability$value)
  expect_gt(sum(scored), 0L)
  expect_true(all(which(scored) %in% cier_flagged_cases(sc, min_votes = 1)))
  expect_identical(cier_flagged_cases(sc, min_votes = 2), integer(0))
})

# ---- cier_screen method: default + monotonicity + bounds --------------------

test_that("min_votes defaults to 1, and a threshold beyond the votes is empty", {
  x <- flag_matrix(120L, 30L, 14L)
  it <- flag_items(10L, 3L)
  sc <- cier_screen(x, it, methods = base8(),
                    control = list(cier_personal_reliability = list(seed = 1)))
  expect_identical(cier_flagged_cases(sc), cier_flagged_cases(sc, min_votes = 1))
  # A threshold beyond the number of votes is reachable by no one.
  expect_identical(cier_flagged_cases(sc, min_votes = ncol(sc$votes) + 1L),
                   integer(0))
})

test_that("an empty screen (no index ran) yields no flagged cases", {
  # items = NULL -> both metadata indices skip -> 0 indices, 0 votes.
  sc <- cier_screen(flag_matrix(30L, 20L, 20L),
                    methods = c("cier_even_odd", "cier_gnormed"))
  expect_identical(length(sc$indices), 0L)
  expect_identical(cier_flagged_cases(sc), integer(0))
  expect_identical(cier_flagged_cases(sc, min_votes = 3), integer(0))
})

# ---- validation + dispatch --------------------------------------------------

test_that("an invalid min_votes is a typed input error", {
  sc <- cier_screen(flag_matrix(40L, 20L, 15L))   # matrix-only, six indices
  for (bad in list(0, -1, 1.5, NA_integer_, c(1L, 2L), "x")) {
    expect_error(cier_flagged_cases(sc, min_votes = bad),
                 class = "cier_error_input")
  }
})

test_that("a non-cier object is a typed input error (default method)", {
  expect_error(cier_flagged_cases(42), class = "cier_error_input")
  expect_error(cier_flagged_cases(list(a = 1)), class = "cier_error_input")
  expect_error(cier_flagged_cases(data.frame(x = 1)), class = "cier_error_input")
})

test_that("output is unnamed positional integers even with row-named input", {
  # A data.frame input with rownames makes an index's `flagged` named; the accessor
  # must still return BARE positional indices (the documented contract), and the
  # index and screen methods must agree.
  x <- flag_matrix(40L, 20L, 7L)
  df <- as.data.frame(x)
  rownames(df) <- paste0("R", seq_len(40L))
  idx <- cier_irv(df)
  fc <- cier_flagged_cases(idx)
  expect_null(names(fc))
  expect_type(fc, "integer")
  it <- flag_items(4L, 5L)
  sc <- cier_screen(df, it, methods = c("cier_irv", "cier_even_odd"))
  expect_null(names(cier_flagged_cases(sc, min_votes = 1)))
})

# ---- new_cier_index invariants ----------------------------------------------

test_that("new_cier_index enforces value / flagged length + type invariants", {
  # `flagged` shorter than `value`: the abstention rule flagged[is.na(value)] <- NA
  # would silently GROW the short vector -- guard it (the footgun) before that line.
  expect_error(
    new_cier_index(c(1, 2, 3), c(TRUE, FALSE), "cier_irv", 2, "lower"),
    class = "cier_error_state"
  )
  # The decisive ordering case: `value` has an NA at the LAST position, so
  # flagged[is.na(value)] <- NA would GROW c(TRUE, FALSE) to length 3 and MASK the
  # mismatch if the length check ran AFTER it. The guard must run BEFORE.
  expect_error(
    new_cier_index(c(1, 2, NA), c(TRUE, FALSE), "cier_irv", 2, "lower"),
    class = "cier_error_state"
  )
  # `value` not numeric.
  expect_error(
    new_cier_index(c("a", "b", "c"), c(TRUE, FALSE, TRUE), "cier_irv", 2, "lower"),
    class = "cier_error_state"
  )
  # `flagged` not logical.
  expect_error(
    new_cier_index(c(1, 2, 3), c(1, 0, 1), "cier_irv", 2, "lower"),
    class = "cier_error_state"
  )
  # `cutoff` not numeric.
  expect_error(
    new_cier_index(c(1, 2, 3), c(TRUE, FALSE, TRUE), "cier_irv", "x", "lower"),
    class = "cier_error_state"
  )
  # `direction` outside {upper, lower}.
  expect_error(
    new_cier_index(c(1, 2, 3), c(TRUE, FALSE, TRUE), "cier_irv", 2, "sideways"),
    class = "cier_error_state"
  )
  # `cutoff` must be a length-1 numeric (NA allowed).
  expect_error(
    new_cier_index(c(1, 2, 3), c(TRUE, FALSE, TRUE), "cier_irv", c(2, 3), "lower"),
    class = "cier_error_state"
  )
  # A valid construction (NA cutoff, NA-aligned flagged) passes unchanged.
  ok <- new_cier_index(c(1, NA, 3), c(TRUE, NA, FALSE), "cier_irv",
                       NA_real_, "lower")
  expect_s3_class(ok, "cier_index")
  expect_identical(ok$flagged, c(TRUE, NA, FALSE))
})
