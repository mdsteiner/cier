# Tests for cier_attention() (direct family; Meade & Craig 2012; Goldammer et al.
# 2024).
#
# Trust model: attention is a counting rule -- per respondent, the number of
# attention checks FAILED among those ANSWERED. A check is failed when its
# answered (non-NA) response is NOT in that check's pass-set. The independent
# oracle (ref-attention-meade-craig-2012.R) re-derives that count by a
# per-respondent ROW loop with %in% and never calls the production kernel; there
# is no CRAN parity partner (no package implements this attention-check counting
# rule as a C/IER index; verified 2026-06-10), so the oracle plus base-R
# primitives are the parity checks -- oracle-only trust, like page_time /
# total_time / PR / RPR. Conventions under test (signed off this slice): one
# wrapper with an explicit per-check pass-set (replacing the archive's three
# rule-based functions instructed/bogus/infrequency); NA = no evidence (a missing
# check is neither passed nor failed; an all-NA respondent abstains) -- a
# DELIBERATE deviation from the archive's cier_instructed, where a missing
# instructed response counted as a failure; the cited fixed cutoff of 1 (any
# failed check flags) with a literal `cutoff` count override (no `frac`); raw
# observed-response coding (no reverse-keying, no items metadata).

source(test_path("..", "reference", "ref-attention-meade-craig-2012.R"))

# ---- Fixtures ---------------------------------------------------------------

# Four respondents, two checks: a bogus-like check (pass = below-midpoint {1, 2})
# and an instructed-like check (pass = the directed option {0}). Hand-checkable.
small_checks <- function() {
  matrix(c(1, 0,    # bogus 1 in {1,2} pass; instr 0 in {0} pass    -> 0 failed
           3, 0,    # bogus 3 not in {1,2} FAIL; instr pass         -> 1 failed
           2, 5,    # bogus 2 in {1,2} pass; instr 5 not in {0} FAIL -> 1 failed
           4, 7),   # bogus FAIL; instr FAIL                         -> 2 failed
         nrow = 4L, byrow = TRUE)
}
small_pass <- function() list(c(1, 2), 0)

# An ASYMMETRIC per-column fixture: column 1 passes only on {1}, column 2 only on
# {3}. The exact values c(0, 2, 1) discriminate two wrong implementations a
# symmetric fixture would let survive: a "global pass-set" mutant that pools the
# pass values across columns ({1, 3} for both) scores row 2 as 1 not 2, and a
# "swap / wrong-column" mutant (pass applied to the other column) scores
# c(2, 1, 2). Both diverge from c(0, 2, 1).
asym_checks <- function() {
  matrix(c(1, 3,    # col1 1 in {1} pass;  col2 3 in {3} pass  -> 0
           3, 5,    # col1 3 not in {1};   col2 5 not in {3}   -> 2
           1, 5),   # col1 1 in {1} pass;  col2 5 not in {3}   -> 1
         nrow = 3L, byrow = TRUE)
}
asym_pass <- function() list(1, 3)

# A realistic spread for oracle parity: k checks with different scales and
# pass-sets, scattered NAs, and n != k so a transposed / wrong-axis membership
# test is caught.
spread_checks <- function(n = 30L, k = 4L, seed = 11L) {
  m <- withr::with_seed(seed,
                        matrix(sample(0:6, n * k, replace = TRUE),
                               nrow = n, ncol = k))
  m[withr::with_seed(seed + 1L, sample(seq_len(n * k), 10L))] <- NA
  list(checks = m, pass = list(c(1, 2), 0, c(5, 6), 3))
}

# ---- Schema -----------------------------------------------------------------

test_that("cier_attention returns a list-based cier_index with the schema", {
  out <- cier_attention(small_checks(), small_pass())
  expect_s3_class(out, "cier_index")
  expect_type(out, "list")
  expect_identical(names(out),
                   c("value", "flagged", "method", "cutoff", "direction"))
  expect_type(out$value, "double")
  expect_type(out$flagged, "logical")
  expect_identical(length(out$value), 4L)
  expect_identical(length(out$flagged), 4L)
  expect_identical(out$method, "cier_attention")
  expect_identical(out$direction, "upper")
  expect_type(out$cutoff, "double")
})

test_that("as.data.frame.cier_index returns the tidy per-respondent frame", {
  df <- as.data.frame(cier_attention(small_checks(), small_pass()))
  expect_s3_class(df, "data.frame")
  expect_identical(names(df), c("value", "flagged"))
  expect_identical(nrow(df), 4L)
  expect_identical(df$value, c(0, 1, 1, 2))
})

# ---- Oracle parity (tolerance 0, integer counts) ----------------------------

test_that("the value equals the independent counting-rule oracle exactly", {
  fx <- spread_checks()
  out <- cier_attention(fx$checks, fx$pass)
  expect_identical(out$value, ref_attention(fx$checks, fx$pass))
})

test_that("attention counts failed checks among answered ones", {
  out <- cier_attention(small_checks(), small_pass())
  expect_equal(out$value, c(0, 1, 1, 2))
  expect_identical(out$value, ref_attention(small_checks(), small_pass()))
})

# ---- Pass-set semantics: per-column, not pooled or swapped ------------------

test_that("the pass-set is applied per column, not pooled or swapped", {
  out <- cier_attention(asym_checks(), asym_pass())
  # c(0, 2, 1): a global-pool mutant would score c(0, 1, 1); a swap mutant
  # c(2, 1, 2). Only the correct per-column membership gives c(0, 2, 1).
  expect_identical(out$value, c(0, 2, 1))
  expect_identical(out$value, ref_attention(asym_checks(), asym_pass()))
})

test_that("failure is non-membership: a passing value is never a failure", {
  # Every response is the (sole) passing value of its check -> 0 failures.
  ch <- matrix(c(1, 3,
                 1, 3),
               nrow = 2L, byrow = TRUE)
  out <- cier_attention(ch, list(1, 3))
  expect_equal(out$value, c(0, 0))
  # Inverting the membership (counting PASSES as failures) would score c(2, 2).
})

test_that("a check with multiple passing values passes on any of them", {
  # Instructed-style single pass {0}; bogus-style pass {1, 2}; infrequency-style
  # pass = everything but the infrequent value 1, i.e. {2, 3, 4, 5}.
  ch <- matrix(c(0, 2, 1,    # instr pass; bogus pass; infreq 1 -> FAIL  -> 1
                 0, 4, 5),   # instr pass; bogus 4 -> FAIL; infreq pass  -> 1
               nrow = 2L, byrow = TRUE)
  out <- cier_attention(ch, list(0, c(1, 2), c(2, 3, 4, 5)))
  expect_equal(out$value, c(1, 1))
})

# ---- Bruhlmann (2020) reproduction from the bundled columns -----------------

test_that("reproduces Bruhlmann's 92 bogus / 96 instructed / 136 any-failure", {
  # The two bundled attention checks: a bogus item (1-5; the authors flag a
  # response of 3 or above, so the pass-set is {1, 2}) and an instructed-response
  # item v_IRI (0-7; the directed option is 0, so the pass-set is {0}).
  bogus <- cier_attention(bfi_careless[, "v_Bogus_Item", drop = FALSE],
                          pass = list(c(1, 2)))
  expect_identical(sum(bogus$value >= 1), 92L)

  instr <- cier_attention(bfi_careless[, "v_IRI", drop = FALSE],
                          pass = list(0))
  expect_identical(sum(instr$value >= 1), 96L)

  both <- cier_attention(bfi_careless[, c("v_Bogus_Item", "v_IRI")],
                         pass = list(c(1, 2), 0))
  expect_identical(sum(both$value >= 1), 136L)   # any failure flags
  expect_identical(sum(both$value), 188)          # 92 + 96 total failures
  expect_identical(sum(both$value == 2), 52L)     # both checks failed
  expect_identical(sum(both$flagged), 136L)       # cutoff 1, any failure flags
})

# ---- NA = no evidence; all-NA abstains --------------------------------------

test_that("an NA check contributes no evidence; all-NA respondents abstain", {
  # 3 respondents, 2 checks; pass = list({1}, {3}).
  # Row 1: c(NA, 5)  -> col1 no evidence; col2 5 not in {3} FAIL   -> 1
  # Row 2: c(NA, NA) -> no answered check                          -> NA (abstain)
  # Row 3: c(1, NA)  -> col1 1 in {1} pass; col2 no evidence       -> 0
  ch <- matrix(c(NA, 5,
                 NA, NA,
                 1, NA),
               nrow = 3L, byrow = TRUE)
  # attention uses a FIXED cutoff, so an abstaining respondent must NOT route the
  # cutoff through the percentile abstention (no cier_warning_insufficient_items).
  out <- expect_no_warning(cier_attention(ch, list(1, 3)))
  expect_equal(out$value[c(1L, 3L)], c(1, 0))
  expect_true(is.na(out$value[[2L]]))
  expect_true(is.na(out$flagged[[2L]]))          # abstainer's flag is NA, aligned
  expect_identical(out$value, ref_attention(ch, list(1, 3)))
  # The fixed cutoff stays 1 despite the abstaining row, and the scored rows still
  # flag against it (kills a percentile-misroute mutant that abstains the cutoff).
  expect_identical(out$cutoff, 1)
  expect_true(out$flagged[[1L]])
  expect_false(out$flagged[[3L]])
})

test_that("NA is no evidence, NOT a failure (deviation from cier_instructed)", {
  # Row c(NA, 5), pass {1}/{3}: a "NA = failure" mutant scores 2 (NA + the genuine
  # fail); the no-evidence rule scores 1.
  out <- cier_attention(matrix(c(NA, 5), nrow = 1L), list(1, 3))
  expect_equal(out$value, 1)
  # Row c(1, NA): a "NA = failure" mutant scores 1; no-evidence scores 0.
  out2 <- cier_attention(matrix(c(1, NA), nrow = 1L), list(1, 3))
  expect_equal(out2$value, 0)
})

test_that("a single answered check is enough to score (only all-NA abstains)", {
  # 1 respondent, 3 checks: one answered (failing), two NA. Kills a mutant that
  # abstains whenever ANY check is NA.
  out <- cier_attention(matrix(c(5, NA, NA), nrow = 1L), list(1, 1, 1))
  expect_equal(out$value, 1)
  expect_false(is.na(out$value))
})

# ---- Fixed-cutoff contract (any failed check) -------------------------------

test_that("default cutoff is fixed 1: any failed check flags (>= pin)", {
  out <- cier_attention(small_checks(), small_pass())
  expect_identical(out$cutoff, 1)
  # value c(0, 1, 1, 2) -> flag where count >= 1. A count of exactly 1 MUST flag
  # (kills a `>` mutant).
  expect_identical(out$flagged, c(FALSE, TRUE, TRUE, TRUE))
  expect_identical(out$flagged, ref_attention_flags(out$value, out$cutoff))
})

test_that("a literal cutoff override flags via value >= cutoff", {
  out <- cier_attention(small_checks(), small_pass(), cutoff = 2)
  expect_identical(out$cutoff, 2)
  # value c(0, 1, 1, 2) -> only the both-failed respondent clears 2.
  expect_identical(out$flagged, c(FALSE, FALSE, FALSE, TRUE))
  expect_identical(out$flagged, ref_attention_flags(out$value, out$cutoff))
})

# ---- Direction --------------------------------------------------------------

test_that("direction is upper: more failed checks flag", {
  out <- cier_attention(small_checks(), small_pass())
  expect_identical(out$direction, "upper")
  expect_true(out$flagged[[which.max(out$value)]])   # most failures
  expect_false(out$flagged[[which.min(out$value)]])  # fewest failures
})

test_that("registry flag_direction matches the wrapper direction", {
  declared <- cier_method_row("cier_attention")$flag_direction
  out <- cier_attention(small_checks(), small_pass())
  expect_identical(out$direction, declared)
})

# ---- checks input validation ------------------------------------------------

test_that("a non-matrix / non-numeric checks is a typed input error", {
  expect_error(cier_attention(c(1, 2, 3), list(1)), class = "cier_error_input")
  expect_error(cier_attention(list(1, 2), list(1)), class = "cier_error_input")
  expect_error(cier_attention(matrix("a", 1L, 1L), list(1)),
               class = "cier_error_input")
  expect_error(cier_attention(NULL, list(1)), class = "cier_error_input")
})

test_that("a data.frame checks is accepted (coerced)", {
  df <- as.data.frame(small_checks())
  expect_no_error(cier_attention(df, small_pass()))
})

test_that("an empty checks is a typed input error", {
  expect_error(
    cier_attention(matrix(numeric(0), nrow = 0L, ncol = 0L), list()),
    class = "cier_error_input"
  )
})

test_that("NaN / infinite checks are typed input errors", {
  expect_error(cier_attention(matrix(c(1, NaN), nrow = 1L), list(1, 1)),
               class = "cier_error_input")
  expect_error(cier_attention(matrix(c(1, Inf), nrow = 1L), list(1, 1)),
               class = "cier_error_input")
})

# ---- pass-spec validation ---------------------------------------------------

test_that("a non-list pass is a typed input error", {
  expect_error(cier_attention(small_checks(), c(1, 2)),
               class = "cier_error_input")
  expect_error(cier_attention(small_checks(), 1),
               class = "cier_error_input")
})

test_that("a pass of the wrong length is a typed input error", {
  # small_checks has 2 columns; pass must have length 2.
  expect_error(cier_attention(small_checks(), list(1)),
               class = "cier_error_input")
  expect_error(cier_attention(small_checks(), list(1, 2, 3)),
               class = "cier_error_input")
})

test_that("a non-numeric / empty / NA-bearing pass element is an error", {
  expect_error(cier_attention(small_checks(), list(1, "a")),
               class = "cier_error_input")
  expect_error(cier_attention(small_checks(), list(1, numeric(0))),
               class = "cier_error_input")
  expect_error(cier_attention(small_checks(), list(1, c(2, NA))),
               class = "cier_error_input")
  expect_error(cier_attention(small_checks(), list(1, c(2, NaN))),
               class = "cier_error_input")
  expect_error(cier_attention(small_checks(), list(1, c(2, Inf))),
               class = "cier_error_input")
})

# ---- cutoff argument validation ---------------------------------------------

test_that("invalid literal cutoff values are typed input errors", {
  k <- ncol(small_checks())   # 2 checks -> cutoff in [1, 2]
  expect_error(cier_attention(small_checks(), small_pass(), cutoff = 0),
               class = "cier_error_input")
  expect_error(cier_attention(small_checks(), small_pass(), cutoff = k + 1L),
               class = "cier_error_input")
  expect_error(cier_attention(small_checks(), small_pass(), cutoff = NA_real_),
               class = "cier_error_input")
  expect_error(cier_attention(small_checks(), small_pass(), cutoff = c(1, 2)),
               class = "cier_error_input")
  expect_error(cier_attention(small_checks(), small_pass(), cutoff = "a"),
               class = "cier_error_input")
})

# ---- print snapshots (shared cier_index print; upper direction) -------------

test_that("print renders the locked cli summary (upper direction)", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    expect_snapshot(print(cier_attention(small_checks(), small_pass())))
  })
})

test_that("print reports abstaining respondents on their own line", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    # Row 2 is all-NA -> abstains; the others score.
    ch <- matrix(c(3, 0,
                   NA, NA,
                   1, 0),
                 nrow = 3L, byrow = TRUE)
    expect_snapshot(print(cier_attention(ch, list(c(1, 2), 0))))
  })
})
