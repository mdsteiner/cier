# The package-wide missing-data / abstention contract.
#
# The n = 1 / p = 1 / all-missing boundary is handled differently by each index (a
# friendly abstention, a typed error, or a silently-computed score), but ONE safety
# property must hold uniformly across every entry point: missing data must never
# manufacture a flag. This file pins that contract once, across all 15 index
# wrappers and cier_screen, so a future index or refactor cannot regress it:
#
#   1. An entirely-missing respondent scores NA, its flag is NA (never TRUE), and it
#      is never returned by cier_flagged_cases (the new_cier_index abstention rule,
#      end to end through every wrapper).
#   2. In cier_screen that respondent contributes a FALSE vote (abstention is a
#      no-vote, not a missing vote), gets zero collapsed votes, and is excluded.
#   3. The "flagged / scored" denominator is sum(!is.na(value)) in BOTH the
#      cier_index print and the cier_screen per-index line -- one convention.
#   4. The n = 1 / p = 1 boundary always resolves to either a valid cier_index or a
#      TYPED cier_error_input, never an opaque base-R error.
#
# (1) reconciles with the all-NA-matrix acceptance fix in test-checks.R: matrix(NA)
# abstains rather than erroring, and that abstention is never a flag.

# ---- fixtures ---------------------------------------------------------------

# Real BFI responses (genuine within-scale correlations, so psychsyn / psychant find
# qualifying pairs -- the contract is then non-vacuous for the pair indices too) with
# one interior respondent set entirely missing.
abst_responses <- function(n = 80L, na_row = 5L) {
  resp <- as.matrix(bfi_careless[seq_len(n), 1:44])
  storage.mode(resp) <- "double"
  resp[na_row, ] <- NA_real_
  resp
}

abst_items <- function() {
  nm <- names(bfi_careless)[1:44]
  data.frame(scale = sub("^v_BFI_([A-Za-z]+)[0-9].*$", "\\1", nm),
             reverse_keyed = grepl("_R$", nm), max = 5L,
             stringsAsFactors = FALSE)
}

# Deterministic timing / direct-check inputs (BFI has no per-page times), each with
# the same interior respondent (row / element 5) entirely missing.
abst_seconds <- function() {
  c(420, 412, 505, 388, NA, 401, 440, 415, 380, 460, 405, 390,
    430, 365, 480, 445, 398, 455, 425, 372, 95, 47, 88, 110)
}

abst_page_seconds <- function() {
  withr::with_seed(7, {
    m <- matrix(stats::runif(24L * 3L, 5, 60), nrow = 24L, ncol = 3L)
  })
  m[5L, ] <- NA_real_
  m
}

abst_checks <- function() {
  withr::with_seed(8, {
    m <- matrix(sample(0:7, 24L * 3L, replace = TRUE), nrow = 24L, ncol = 3L)
  })
  m[5L, ] <- NA_real_
  m
}

# A copy of the BFI fixture with one ITEM (column) entirely missing -- a dead item,
# the all-NA-column companion to the all-NA row above.
abst_responses_deadcol <- function(n = 80L, dead_col = 3L) {
  resp <- as.matrix(bfi_careless[seq_len(n), 1:44])
  storage.mode(resp) <- "double"
  resp[, dead_col] <- NA_real_
  resp
}

# The twelve response-matrix / items index thunks for given responses + items, each
# tagged with the optional Suggests backend it needs. Shared by the missing-row
# scenario (abst_runners) and the dead-column scenario, so both reach the same index
# set.
abst_matrix_runners <- function(resp, items) {
  list(
    list(id = "cier_longstring", fn = function() cier_longstring(resp)),
    list(id = "cier_irv", fn = function() cier_irv(resp)),
    list(id = "cier_psychsyn", fn = function() cier_psychsyn(resp)),
    list(id = "cier_psychant", fn = function() cier_psychant(resp)),
    list(id = "cier_mahalanobis", fn = function() cier_mahalanobis(resp)),
    list(id = "cier_person_total", fn = function() cier_person_total(resp)),
    list(id = "cier_autocorrelation", fn = function() cier_autocorrelation(resp)),
    list(id = "cier_lazr", fn = function() cier_lazr(resp)),
    list(id = "cier_even_odd", fn = function() cier_even_odd(resp, items)),
    list(id = "cier_personal_reliability",
         fn = function() cier_personal_reliability(resp, items, seed = 1)),
    # Both person-fit indices are pure R now (no backend), so neither carries a
    # `needs` -- they always run.
    list(id = "cier_gnormed", fn = function() cier_gnormed(resp, items, seed = 1)),
    list(id = "cier_ht", fn = function() cier_ht(resp, items))
  )
}

# Every index entry point as a thunk returning its cier_index, tagged with the row
# that is entirely missing (na) and the optional Suggests backend it needs. The
# matrix / items indices share one all-missing-row matrix; the timing and direct
# families bring their own. Backend indices are filtered out (not skipped) when their
# package is absent, so the rest of the contract still runs.
abst_runners <- function() {
  resp  <- abst_responses()
  items <- abst_items()
  secs  <- abst_seconds()
  pgs   <- abst_page_seconds()
  ipp   <- c(5L, 6L, 4L)
  chk   <- abst_checks()
  pass  <- list(c(1, 2), 0, c(3, 4))
  mat <- lapply(abst_matrix_runners(resp, items),
                function(r) c(r, list(na = 5L)))
  c(mat, list(
    list(id = "cier_total_time", na = 5L, fn = function() cier_total_time(secs)),
    list(id = "cier_page_time", na = 5L, fn = function() cier_page_time(pgs, ipp)),
    list(id = "cier_attention", na = 5L, fn = function() cier_attention(chk, pass))
  ))
}

# Drop runners whose Suggests backend is not installed (keep the rest running).
abst_available <- function(runners) {
  Filter(function(r) {
    is.null(r$needs) || requireNamespace(r$needs, quietly = TRUE)
  }, runners)
}

# Run a thunk expected to abstain on a missing row, muffling the degenerate-input
# warnings every index legitimately raises here (saturated cutoff, insufficient
# items, no pairs, singular covariance, short battery) and the forward-keyed message.
run_quiet <- function(fn) {
  suppressMessages(suppressWarnings(fn()))
}

# ---- 1. an entirely-missing respondent is never flagged, everywhere ----------

# Indices that legitimately flag NOBODY on this clean-ish BFI fixture, so the
# "flags >= 1 respondent" liveness control below does not apply to them:
# cier_psychant finds no antonym pairs at the default critical_r (abstains
# whole-sample on BFI's positively-correlated items), and cier_mahalanobis at the
# chi-square p = 0.001 cutoff flags no outlier here. Every OTHER index flags at least
# one respondent, which makes "the NA row is excluded from a NON-EMPTY flagged set" a
# non-vacuous claim and kills a "flag nobody, ever" regression (which would otherwise
# pass assertions 1-3 vacuously).
abst_no_flag_on_fixture <- function() c("cier_psychant", "cier_mahalanobis")

test_that("a fully-missing respondent abstains and is never flagged (every index, C18)", {
  for (r in abst_available(abst_runners())) {
    out <- run_quiet(r$fn)
    expect_s3_class(out, "cier_index")
    # The all-missing respondent scores NA ...
    expect_true(is.na(out$value[r$na]), info = paste(r$id, "value[na] NA"))
    # ... its flag is forced NA (never TRUE) ...
    expect_true(is.na(out$flagged[r$na]), info = paste(r$id, "flagged[na] NA"))
    # ... and it is never returned as a flagged case.
    fc <- cier_flagged_cases(out)
    expect_false(r$na %in% fc, info = paste(r$id, "not in flagged_cases"))
    # General invariant: no abstaining (NA-valued) respondent is ever flagged TRUE
    # (%in% TRUE treats NA as not-TRUE), so missing data manufactures no flag.
    abst <- is.na(out$value)
    expect_false(any(out$flagged[abst] %in% TRUE), info = paste(r$id, "no abstainer flagged"))
    # Liveness: for every index that flags on this fixture, the flagged set is
    # non-empty AND scores most respondents -- so the NA-row exclusion above is
    # tested against a real, populated flag set, not an all-abstain no-op.
    if (!(r$id %in% abst_no_flag_on_fixture())) {
      expect_gt(length(fc), 0L, label = paste(r$id, "flagged count"))
      expect_gt(sum(!is.na(out$value)), 1L, label = paste(r$id, "scored count"))
    }
  }
})

# ---- 1b. a dead item (all-missing column) manufactures no flag --------------

test_that("a dead item (all-missing column) is handled gracefully, no flag", {
  # The all-NA-column companion: one item is entirely missing. Every response index
  # must still return a typed cier_index (no crash, no opaque error) -- several
  # whole-abstain because the dead column breaks their complete-case / covariance
  # requirement, which is exactly the graceful path -- and crucially no NA-valued
  # respondent is ever flagged.
  resp <- abst_responses_deadcol()
  items <- abst_items()
  for (r in abst_available(abst_matrix_runners(resp, items))) {
    out <- run_quiet(r$fn)
    expect_s3_class(out, "cier_index")
    expect_false(any(out$flagged[is.na(out$value)] %in% TRUE),
                 info = paste(r$id, "dead column manufactures no flag"))
  }
})

# ---- 2. abstention is a no-vote in cier_screen ------------------------------

test_that("a fully-missing respondent gets zero votes and no flag in cier_screen", {
  na_row <- 5L
  resp <- abst_responses(na_row = na_row)
  items <- abst_items()
  sc <- run_quiet(function() {
    cier_screen(resp, items,
                control = list(cier_personal_reliability = list(seed = 1),
                               cier_gnormed = list(seed = 1)))
  })
  # Every index abstained on the missing row, so its whole flag row is NA ...
  expect_true(all(is.na(as.matrix(sc$flags)[na_row, ])))
  # ... the collapsed votes carry no NA (abstain -> FALSE) and that respondent
  # reaches zero votes ...
  expect_false(anyNA(as.matrix(sc$votes)))
  expect_identical(rowSums(as.matrix(sc$votes))[[na_row]], 0)
  # ... so it is excluded from the flagged cases at every threshold.
  expect_false(na_row %in% cier_flagged_cases(sc, min_votes = 1L))
  # Positive control: the screen is otherwise live (some respondents ARE flagged),
  # so the exclusion above is not vacuous.
  expect_gt(length(cier_flagged_cases(sc, min_votes = 1L)), 0L)
})

# ---- 3. one scored-denominator convention across index and screen -----------

test_that("the flagged/scored denominator is sum(!is.na(value)) in index AND screen", {
  na_row <- 5L
  resp <- abst_responses(na_row = na_row)
  items <- abst_items()
  withr::with_options(list(cli.width = 200, cli.unicode = FALSE), {
    idx <- run_quiet(function() cier_irv(resp))
    scored <- sum(!is.na(idx$value))
    expect_true(scored < length(idx$value))        # an abstainer is present
    # The single index reports "Flagged: F of {scored} scored".
    idx_txt <- paste(capture.output(print(idx)), collapse = " ")
    expect_match(idx_txt, sprintf("Flagged: [0-9]+ of %d scored", scored))

    sc <- run_quiet(function() cier_screen(resp, items, methods = "cier_irv"))
    # The same index reported via the screen uses the SAME scored denominator on
    # its per-index line "cier_irv  F / {scored} (..%)".
    sc_txt <- paste(capture.output(print(sc)), collapse = " ")
    expect_match(sc_txt, sprintf("cier_irv +[0-9]+ / %d ", scored))
    # And the screen's own count of scored respondents for that index equals it.
    expect_identical(sum(!is.na(sc$indices$cier_irv$value)), scored)
  })
})

# ---- 4. the n=1 / p=1 boundary is always a typed condition ------------------

test_that("n=1 returns a cier_index for every index", {
  # The single-respondent boundary is UNIFORM: every entry point returns a valid
  # cier_index (abstaining or scoring), never an error. A mutant that errors out at
  # n = 1 -- or returns something other than a cier_index -- fails here.
  items <- abst_items()
  one_resp <- as.matrix(bfi_careless[1L, 1:44, drop = FALSE])
  storage.mode(one_resp) <- "double"
  expect_index <- function(thunk, info) {
    out <- suppressMessages(suppressWarnings(thunk()))
    expect_s3_class(out, "cier_index")
  }
  for (m in c("cier_longstring", "cier_irv", "cier_psychsyn", "cier_psychant",
              "cier_mahalanobis", "cier_person_total", "cier_autocorrelation",
              "cier_lazr")) {
    f <- match.fun(m)
    expect_index(function() f(one_resp), info = paste(m, "n=1"))
  }
  expect_index(function() cier_even_odd(one_resp, items), "even_odd n=1")
  expect_index(function() cier_personal_reliability(one_resp, items, seed = 1),
               "PR n=1")
  expect_index(function() cier_gnormed(one_resp, items, seed = 1), "gnormed n=1")
  expect_index(function() cier_ht(one_resp, items), "ht n=1")  # pure R, always runs
  expect_index(function() cier_total_time(420), "total_time n=1")
  expect_index(function() cier_page_time(matrix(30, 1L, 2L), c(5L, 6L)),
               "page_time n=1")
  expect_index(function() cier_attention(matrix(c(1, 0), 1L, 2L), list(1, 0)),
               "attention n=1")
})

test_that("p=1 pins the documented error-vs-compute split", {
  # The single-ITEM boundary is the heterogeneous one, and the split is a contract,
  # not an accident: the three indices whose statistic is undefined on one item raise
  # a TYPED cier_error_input (autocorrelation needs >= 4 items; even-odd and personal
  # reliability need >= 2 scales, so one item cannot supply them), while every other
  # index returns a (typically abstaining) cier_index. Pinning the split catches a
  # mutant that, say, silently computes a degenerate even-odd score at p = 1 instead
  # of erroring.
  one_item_resp <- as.matrix(bfi_careless[seq_len(40L), 1, drop = FALSE])
  storage.mode(one_item_resp) <- "double"
  one_item <- abst_items()[1L, , drop = FALSE]
  expect_index <- function(thunk, info) {
    out <- suppressMessages(suppressWarnings(thunk()))
    expect_s3_class(out, "cier_index")
  }
  expect_typed_error <- function(thunk, info) {
    expect_error(suppressMessages(suppressWarnings(thunk())),
                 class = "cier_error_input")
  }
  # must error: statistic undefined on a single item
  expect_typed_error(function() cier_autocorrelation(one_item_resp), "autocorr p=1")
  expect_typed_error(function() cier_even_odd(one_item_resp, one_item), "even_odd p=1")
  expect_typed_error(function() cier_personal_reliability(one_item_resp, one_item, seed = 1),
                     "PR p=1")
  # must return a cier_index (abstain or compute)
  for (m in c("cier_longstring", "cier_irv", "cier_psychsyn", "cier_psychant",
              "cier_mahalanobis", "cier_person_total", "cier_lazr")) {
    f <- match.fun(m)
    expect_index(function() f(one_item_resp), info = paste(m, "p=1"))
  }
  expect_index(function() cier_gnormed(one_item_resp, one_item, seed = 1), "gnormed p=1")
  expect_index(function() cier_ht(one_item_resp, one_item), "ht p=1")  # pure R
  expect_index(function() cier_page_time(matrix(c(30, 40, 25), 3L, 1L), 5L),
               "page_time p=1")
  expect_index(function() cier_attention(matrix(c(1, 3, 2), 3L, 1L), list(1)),
               "attention p=1")
})
