# Tests for cier_attention() (direct family; Meade & Craig 2012; Goldammer et al.
# 2024).
#
# Trust model: attention counts, per respondent, the attention checks FAILED among
# those ANSWERED. A check fails when its answered (non-NA) response is NOT in that
# check's pass-set. The oracle (ref-attention-meade-craig-2012.R) re-derives the
# count by a per-respondent ROW loop with %in% and never calls the kernel; no CRAN
# parity partner (no package implements this attention-check counting rule as a
# C/IER index; verified 2026-06-10), so oracle + base-R primitives are the parity
# checks -- oracle-only trust, like page_time / total_time / PR / RPR. Conventions
# under test: one wrapper with an explicit per-check pass-set (replacing the
# archive's three rule-based functions instructed/bogus/infrequency); NA = no
# evidence (a missing check is neither passed nor failed; an all-NA respondent
# abstains) -- a DELIBERATE deviation from the archive's cier_instructed, where a
# missing instructed response counted as a failure; the cited fixed cutoff of 1
# (any failed check flags) with a literal `cutoff` count override (no `frac`); raw
# observed-response coding (no reverse-keying, no items metadata).

source(test_path("..", "reference", "ref-attention-meade-craig-2012.R"))

# ---- Fixtures ---------------------------------------------------------------

# Four respondents, two checks: a bogus-like check (pass = below-midpoint {1, 2})
# and an instructed-like check (pass = directed option {0}). Hand-checkable.
small_checks <- function() {
  matrix(c(1, 0,    # bogus 1 in {1,2} pass; instr 0 in {0} pass    -> 0 failed
           3, 0,    # bogus 3 not in {1,2} FAIL; instr pass         -> 1 failed
           2, 5,    # bogus 2 in {1,2} pass; instr 5 not in {0} FAIL -> 1 failed
           4, 7),   # bogus FAIL; instr FAIL                         -> 2 failed
         nrow = 4L, byrow = TRUE)
}
small_pass <- function() list(c(1, 2), 0)

# An ASYMMETRIC per-column fixture: column 1 passes only on {1}, column 2 only on
# {3}. The values c(0, 2, 1) discriminate two mutants a symmetric fixture would let
# survive: a "global pass-set" mutant pooling pass values across columns ({1, 3}
# for both) scores row 2 as 1 not 2, and a "swap / wrong-column" mutant (pass on
# the other column) scores c(2, 1, 2). Both diverge from c(0, 2, 1).
asym_checks <- function() {
  matrix(c(1, 3,    # col1 1 in {1} pass;  col2 3 in {3} pass  -> 0
           3, 5,    # col1 3 not in {1};   col2 5 not in {3}   -> 2
           1, 5),   # col1 1 in {1} pass;  col2 5 not in {3}   -> 1
         nrow = 3L, byrow = TRUE)
}
asym_pass <- function() list(1, 3)

# A two-column fixture with NAMED columns, for the positional name-binding guard.
# Column "A" all 1 (passes {1}); column "B" all 3 (passes {3}). The correct
# positional pass list(A = 1, B = 3) scores c(0, 0); a name-reversed
# list(B = 3, A = 1) binds POSITIONALLY (pass[[1]] = 3 lands on column A,
# pass[[2]] = 1 on column B) and would silently score c(2, 2) -- the trap the
# guard rejects.
named_checks <- function() {
  matrix(c(1, 3,
           1, 3),
         nrow = 2L, byrow = TRUE,
         dimnames = list(NULL, c("A", "B")))
}

# A realistic spread for oracle parity: k checks with different scales and
# pass-sets, scattered NAs, n != k (catches a transposed / wrong-axis membership
# test).
spread_checks <- function(n = 30L, k = 4L, seed = 11L) {
  m <- withr::with_seed(seed,
                        matrix(sample(0:6, n * k, replace = TRUE),
                               nrow = n, ncol = k))
  m[withr::with_seed(seed + 1L, sample(seq_len(n * k), 10L))] <- NA
  list(checks = m, pass = list(c(1, 2), 0, c(5, 6), 3))
}

# ---- Schema -----------------------------------------------------------------

test_that("cier_attention returns the pinned cier_index schema", {
  expect_cier_index_schema(cier_attention(small_checks(), small_pass()),
                           "cier_attention", "upper", 4L)
})

# ---- Oracle parity (tolerance 0, integer counts) ----------------------------

test_that("the value equals the independent counting-rule oracle exactly", {
  fx <- spread_checks()
  # This spread flags every scored respondent (4 checks scattered over 0:6),
  # tripping the all-flagged diagnostic; suppress the additive warning -- the
  # oracle-parity value is under test here.
  out <- suppressWarnings(cier_attention(fx$checks, fx$pass))
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
  # value c(1, 1): every scored respondent fails a check -> 100% flagged, emitting
  # the all-flagged diagnostic; suppress it -- the multi-value pass semantics, not
  # the warning, are under test (the warning has its own tests below).
  out <- suppressWarnings(cier_attention(ch, list(0, c(1, 2), c(2, 3, 4, 5))))
  expect_equal(out$value, c(1, 1))
})

# ---- Bruhlmann (2020) reproduction from the bundled columns -----------------

test_that("reproduces Bruhlmann's 92 bogus / 96 instructed / 136 any-failure", {
  # The two bundled attention checks: a bogus item (1-5; the authors flag a response
  # of 3 or above, so pass-set {1, 2}) and an instructed-response item v_IRI (0-7;
  # directed option 0, so pass-set {0}).
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
  # Row 1: c(NA, 5) -> col1 no evidence; col2 5 not in {3} FAIL    -> 1
  # Row 2: c(NA, NA) -> no answered check                          -> NA (abstain)
  # Row 3: c(1, 3) -> col1 1 in {1} pass; col2 3 in {3} pass       -> 0
  # Row 3 answers col2 with a PASSING 3, so col2's pass-set {3} is not disjoint from
  # its observed values -- keeping this fixture clear of the disjoint-pass warning,
  # so the expect_no_warning below holds.
  ch <- matrix(c(NA, 5,
                 NA, NA,
                 1, 3),
               nrow = 3L, byrow = TRUE)
  # attention uses a FIXED cutoff, so an abstaining respondent must NOT route the
  # cutoff through the percentile abstention (no cier_warning_insufficient_items).
  out <- expect_no_warning(cier_attention(ch, list(1, 3)))
  expect_equal(out$value[c(1L, 3L)], c(1, 0))
  expect_true(is.na(out$value[[2L]]))
  expect_true(is.na(out$flagged[[2L]]))          # abstainer's flag is NA, aligned
  expect_identical(out$value, ref_attention(ch, list(1, 3)))
  # The fixed cutoff stays 1 despite the abstaining row, and the scored rows still
  # flag against it (kills a percentile-misroute mutant abstaining the cutoff).
  expect_identical(out$cutoff, 1)
  expect_true(out$flagged[[1L]])
  expect_false(out$flagged[[3L]])
})

test_that("NA is no evidence, NOT a failure (deviation from cier_instructed)", {
  # Row c(NA, 5), pass {1}/{3}: a "NA = failure" mutant scores 2 (NA + the genuine
  # fail); the no-evidence rule scores 1. The lone scored respondent fails (100%
  # flagged) and col2's {3} never appears (disjoint), so both diagnostics fire;
  # suppress -- the NA-as-no-evidence value, not the warnings, is under test.
  out <- suppressWarnings(cier_attention(matrix(c(NA, 5), nrow = 1L), list(1, 3)))
  expect_equal(out$value, 1)
  # Row c(1, NA): a "NA = failure" mutant scores 1; no-evidence scores 0.
  out2 <- cier_attention(matrix(c(1, NA), nrow = 1L), list(1, 3))
  expect_equal(out2$value, 0)
})

test_that("a single answered check is enough to score (only all-NA abstains)", {
  # 1 respondent, 3 checks: one answered (failing), two NA. Kills a mutant abstaining
  # whenever ANY check is NA. The lone failing respondent trips all-flagged +
  # disjoint; suppress them -- the single-answered-check scoring is under test.
  out <- suppressWarnings(cier_attention(matrix(c(5, NA, NA), nrow = 1L), list(1, 1, 1)))
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

test_that("method spec flag_direction matches the wrapper direction", {
  declared <- cier_method_spec("cier_attention")$flag_direction
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

test_that("a fractional literal cutoff is a typed input error", {
  # A failed-check cutoff is a whole count; 1.5 is meaningless.
  expect_error(cier_attention(small_checks(), small_pass(), cutoff = 1.5),
               class = "cier_error_input")
  expect_s3_class(cier_attention(small_checks(), small_pass(), cutoff = 2),
                  "cier_index")
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

# ---- positional name-binding guard ------------------------------------------
# `pass` binds to the columns of `checks` BY POSITION. When the user names the pass
# list AND a name matches some column, that is a name-binding attempt, so the names
# must equal the columns in order -- else the user is silently scoring the wrong
# column. Names matching no column are decoration and ignored.

test_that("a name-reversed pass against named columns is a typed input error", {
  # names c("B", "A") are the columns reversed; positional binding would score
  # c(2, 2) not c(0, 0). Reject rather than mis-score (the verified trap).
  expect_error(
    cier_attention(named_checks(), list(B = 3, A = 1)),
    class = "cier_error_input"
  )
})

test_that("a pass named in column order is accepted and scores positionally", {
  # Identical names in the same order -> no abort; identical result to unnamed.
  named <- cier_attention(named_checks(), list(A = 1, B = 3))
  plain <- cier_attention(named_checks(), list(1, 3))
  expect_identical(named$value, c(0, 0))
  expect_identical(named$value, plain$value)
})

test_that("decorative pass names that match no column are ignored", {
  # Names overlap NO column name -> treated as labels, scored positionally. Kills
  # a strict mutant that aborts on ANY named pass.
  out <- cier_attention(named_checks(), list(bogus = 1, instructed = 3))
  expect_identical(out$value, c(0, 0))
})

test_that("partial naming that reorders by an overlapping name aborts", {
  # names c("B", "") : "B" overlaps a column -> a binding attempt -> demand an
  # exact order. pass[[1]] ("B") lands on column A and would mis-score. Kills a
  # set-only mutant that lets a non-permutation partial name through.
  expect_error(
    cier_attention(named_checks(), stats::setNames(list(3, 1), c("B", ""))),
    class = "cier_error_input"
  )
})

test_that("a named pass with UNNAMED columns is accepted (no cross-check)", {
  # No column names -> nothing to cross-check -> positional binding (as documented).
  # Kills a mutant aborting whenever pass is named, ignoring the columns. The
  # positional binding here is disjoint/all-flagged, hence suppressed.
  plain_cols <- matrix(c(1, 3, 1, 3), nrow = 2L, byrow = TRUE)   # no dimnames
  out <- suppressWarnings(cier_attention(plain_cols, list(B = 3, A = 1)))
  expect_identical(out$value, c(2, 2))
})

test_that("all-empty pass names are treated as unnamed", {
  # setNames(..., c("", "")) has non-NULL but blank names -> not a binding
  # attempt. Kills a mutant triggering on `!is.null(names(pass))` alone.
  out <- cier_attention(named_checks(), stats::setNames(list(1, 3), c("", "")))
  expect_identical(out$value, c(0, 0))
})

test_that("a partial naming with a blank at a named column aborts", {
  # "A" overlaps column A -> a binding attempt -> the positional rule applies, and
  # the blank slot 2 against named column "B" is a mismatch. Kills a mutant comparing
  # only the non-blank names (ignoring blanks at named columns).
  err <- expect_error(
    cier_attention(named_checks(), stats::setNames(list(1, 3), c("A", ""))),
    class = "cier_error_input"
  )
  expect_identical(cier_condition_data(err)$observed, 2L)   # the blank position
})

test_that("the name-mismatch error reports the positions and binds by position", {
  err <- expect_error(
    cier_attention(named_checks(), list(B = 3, A = 1)),
    class = "cier_error_input"
  )
  # Both positions disagree -> payload names them (kills a generic-message mutant
  # aborting without reporting WHICH positions clash).
  expect_identical(cier_condition_data(err)$observed, c(1L, 2L))
  msg <- gsub("\\s+", " ", rlang::cnd_message(err))
  expect_match(msg, "do not match", fixed = TRUE)
  expect_match(msg, "by position", fixed = TRUE)
})

# ---- disjoint pass-set warning ----------------------------------------------
# A pass-set sharing NO value with a column's observed (non-NA) responses flags
# every answered respondent on that column -- almost always a mis-specified pass-set
# (wrong values / outside the column's coding). Warn (not error: a never-endorsed
# value is legal), naming the column. An all-NA column has nothing to share and is
# exempt; the statistic is unchanged.

test_that("a disjoint pass-set warns once, names the column, leaves scoring intact", {
  # Column "A" is answered only 5, but its pass-set is {0}: nothing it endorses ever
  # appears. Column "B" overlaps its pass-set, and the col-B-only respondent (value
  # 0) keeps the sample off the all-flagged path, isolating this warning.
  checks <- matrix(c(5, 1,
                     5, 2,
                     NA, 1),
                   nrow = 3L, byrow = TRUE,
                   dimnames = list(NULL, c("A", "B")))
  cond <- expect_warning(out <- cier_attention(checks, list(0, c(1, 2))),
                         class = "cier_warning_disjoint_pass")
  expect_identical(out$value, c(1, 1, 0))               # diagnostic only
  expect_identical(cier_condition_data(cond)$observed, "A")
  w <- testthat::capture_warnings(cier_attention(checks, list(0, c(1, 2))))
  expect_length(w, 1L)                                  # only the disjoint warn
  expect_match(w, "matches none", fixed = TRUE)
  expect_match(w, "pass = list(0)", fixed = TRUE)       # the approved hint
})

test_that("a disjoint warning on unnamed columns names the column position", {
  checks <- matrix(c(5, 1,
                     5, 2,
                     NA, 1),
                   nrow = 3L, byrow = TRUE)             # no colnames
  cond <- expect_warning(cier_attention(checks, list(0, c(1, 2))),
                         class = "cier_warning_disjoint_pass")
  expect_identical(cier_condition_data(cond)$observed, "column 1")
})

test_that("two disjoint columns warn once, listing both", {
  # Both columns disjoint; cutoff = 2 leaves a one-failure respondent unflagged so
  # the all-flagged warning stays out and disjoint is the only one.
  checks <- matrix(c(5, 5,
                     5, NA),
                   nrow = 2L, byrow = TRUE,
                   dimnames = list(NULL, c("A", "B")))
  cond <- expect_warning(cier_attention(checks, list(0, 0), cutoff = 2),
                         class = "cier_warning_disjoint_pass")
  expect_identical(cier_condition_data(cond)$observed, c("A", "B"))
  w <- testthat::capture_warnings(cier_attention(checks, list(0, 0), cutoff = 2))
  expect_length(w, 1L)                                  # ONE warning, not per-col
  expect_match(w, "pass-sets", fixed = TRUE)            # plural rendered
})

test_that("a pass-set overlapping the observed coding does not warn", {
  # small_pass overlaps every column's observed values and the sample is not
  # all-flagged -> no warning (the guards stay silent on clean data).
  expect_no_warning(cier_attention(small_checks(), small_pass()))
})

test_that("an all-NA column is not treated as disjoint", {
  # Column "A" is entirely NA (no observed value to share) -> not a misspecified
  # pass-set. Kills a mutant missing the observed-non-empty guard (an empty
  # intersection is vacuously disjoint).
  checks <- matrix(c(NA, 1,
                     NA, 2),
                   nrow = 2L, byrow = TRUE,
                   dimnames = list(NULL, c("A", "B")))
  expect_no_warning(cier_attention(checks, list(0, c(1, 2))))
})

test_that("a pass-set sharing at least one observed value does not warn", {
  # Observed {2, 3, 4}; pass-set {1, 2} shares 2 -> not disjoint, even though 3 and 4
  # fail. Kills a mutant warning unless EVERY observed value passes.
  checks <- matrix(c(2, 3, 4), ncol = 1L)
  expect_no_warning(out <- cier_attention(checks, list(c(1, 2))))
  expect_identical(out$value, c(0, 1, 1))
})

# ---- all-flagged warning ----------------------------------------------------
# When every SCORED respondent is flagged, a mis-specified pass-set is the usual
# cause. Mirrors page_time's no-minimum-n all-flagged diagnostic, over scored rows
# only (an abstaining row neither blocks nor triggers it).

test_that("every scored respondent flagged emits the all-flagged warning", {
  # Each column's pass-set overlaps its observed values (no disjoint warning), yet
  # every respondent fails a check -> 100% flagged. Isolates this warning.
  checks <- matrix(c(5, 3,
                     1, 5,
                     5, 5),
                   nrow = 3L, byrow = TRUE,
                   dimnames = list(NULL, c("A", "B")))
  cond <- expect_warning(out <- cier_attention(checks, list(1, 3)),
                         class = "cier_warning_all_flagged")
  expect_identical(out$value, c(1, 1, 2))
  expect_true(all(out$flagged))
  expect_identical(cier_condition_data(cond)$n_scored, 3L)
  w <- testthat::capture_warnings(cier_attention(checks, list(1, 3)))
  expect_length(w, 1L)
  expect_match(w, "Every scored respondent was flagged", fixed = TRUE)
  expect_match(w, "mis-specified", fixed = TRUE)
})

test_that("all-flagged warns at small n with no minimum threshold", {
  # n = 2, both flagged, neither column disjoint -> only the all-flagged warning,
  # at the smallest non-trivial sample. Kills any `n_scored >= k` mutant.
  checks <- matrix(c(1, 5,
                     5, 3),
                   nrow = 2L, byrow = TRUE)
  w <- testthat::capture_warnings(cier_attention(checks, list(1, 3)))
  expect_length(w, 1L)
  expect_match(w, "Every scored respondent was flagged", fixed = TRUE)
})

test_that("all-flagged warns when fewer rows are scored than there are checks", {
  # 1 scored respondent, 3 checks: n_scored (1) < k (3). Kills an `n_scored >= k`
  # mutant the n = k = 2 fixture above cannot (2 >= 2 is still TRUE). A single-row
  # all-fail is necessarily disjoint, so match by text among the several warnings
  # rather than asserting a single one.
  w <- testthat::capture_warnings(cier_attention(matrix(c(5, 3, 5), nrow = 1L), list(1, 1, 1)))
  expect_true(any(grepl("Every scored respondent was flagged", w, fixed = TRUE)))
})

test_that("a mixed flag set does not emit the all-flagged warning", {
  # small flags c(F, T, T, T): not all -> no warning of that class at all.
  cond <- rlang::catch_cnd(cier_attention(small_checks(), small_pass()),
                           classes = "cier_warning_all_flagged")
  expect_null(cond)
})

test_that("an all-abstaining sample does not emit the all-flagged warning", {
  # Every respondent all-NA -> no scored row -> the predicate is vacuously false.
  # Kills an `all(flagged)` mutant ignoring the scored mask (NA -> not TRUE).
  checks <- matrix(NA_real_, nrow = 3L, ncol = 2L)
  expect_no_warning(out <- cier_attention(checks, list(1, 3)))
  expect_true(all(is.na(out$flagged)))
})

test_that("all-flagged warns on the scored set despite an abstaining row", {
  # Two scored rows both flagged + one all-NA abstainer; the abstainer neither
  # blocks nor triggers the warning.
  checks <- matrix(c(1, 5,
                     5, 3,
                     NA, NA),
                   nrow = 3L, byrow = TRUE)
  expect_warning(cier_attention(checks, list(1, 3)),
                 class = "cier_warning_all_flagged")
})

test_that("the cutoff provenance survives the explicit assemble", {
  # The all-flagged warning forces the explicit apply_flag/assemble form (not
  # flag_and_assemble); the default cutoff stays fixed_count and a literal cutoff
  # "literal", each with no rate.
  def <- cier_attention(small_checks(), small_pass())
  expect_identical(def$cutoff_method, "fixed_count")
  expect_true(is.na(def$cutoff_rate))
  lit <- cier_attention(small_checks(), small_pass(), cutoff = 2)
  expect_identical(lit$cutoff_method, "literal")
  expect_true(is.na(lit$cutoff_rate))
})
