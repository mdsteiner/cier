# Tests for cier_lazr() (Laz.R, Biemann et al. 2025).
#
# Trust model: the oracle (ref_lazr) re-derives the per-respondent Laz.R with an
# explicit transition double-loop and NEVER calls the kernel. Laz.R has no CRAN
# parity partner (verified 2026-06-10), so the oracle plus the paper's worked
# examples (John = 33/49; the one-liner Laz.R(c(1,2,3,4,5,4,3,2,1,2)) = 2/3) are
# the parity checks, like PR / RPR. Conventions under test: drop-NA transitions,
# abstain below two valid transitions, matrix-only (anchor-count-invariant)
# scoring, integer-coded responses, percentile / upper cutoff.

source(test_path("..", "reference", "ref-lazr-biemann2025.R"))

# ---- Fixtures ---------------------------------------------------------------

lazr_fixture <- function(n = 20L, p = 30L, seed = 2026L) {
  withr::with_seed(seed, {
    m <- matrix(sample.int(5L, n * p, replace = TRUE), nrow = n, ncol = p)
  })
  storage.mode(m) <- "double"
  m
}

john_matrix <- function(extra_rows = 0L) {
  s <- ref_lazr_john_sequence()
  p <- length(s)
  m <- matrix(s, nrow = 1L, ncol = p)
  if (extra_rows > 0L) {
    m <- rbind(m, matrix(rep(s, extra_rows), nrow = extra_rows, byrow = TRUE))
  }
  storage.mode(m) <- "double"
  m
}

# ---- Schema -----------------------------------------------------------------

test_that("cier_lazr returns the pinned cier_index schema", {
  out <- suppressWarnings(cier_lazr(lazr_fixture(n = 12L)))
  expect_cier_index_schema(out, "cier_lazr", "upper", 12L)
})

test_that("cier_lazr strips respondent row names from $value and $flagged", {
  # n_trans = rowSums(valid) carries the input row names into value; strip them so
  # $value/$flagged are bare positional vectors like the direct/timing kernels (the
  # contract cier_flagged_cases() and as.data.frame.cier_screen() rely on). A
  # >= 20-row fixture RESOLVES the percentile cutoff, so $flagged is a real named
  # TRUE/FALSE vector without the fix -- not the unnamed all-FALSE an abstaining (NA)
  # cutoff would hand back for free, which would let a "strip $value only" mutant
  # survive.
  x <- lazr_fixture(n = 30L)
  df <- as.data.frame(x)
  rownames(df) <- paste0("R", seq_len(nrow(df)))
  named <- suppressWarnings(cier_lazr(df))
  expect_false(is.na(named$cutoff))    # cutoff resolved -> $flagged is real
  expect_null(names(named$value))
  expect_null(names(named$flagged))
  # unname() strips names ONLY: the scores match the unnamed-matrix input exactly.
  expect_identical(named$value, suppressWarnings(cier_lazr(x))$value)
})

# ---- Published worked examples (the spec) -----------------------------------

test_that("cier_lazr reproduces John's worked example (Eq. 4): 33/49", {
  out <- suppressWarnings(cier_lazr(john_matrix()))
  expect_equal(out$value[[1L]], 33 / 49, tolerance = 1e-12)
})

test_that("footnote-2 one-liner: Laz.R(c(1,2,3,4,5,4,3,2,1,2)) = 2/3", {
  m <- matrix(c(1, 2, 3, 4, 5, 4, 3, 2, 1, 2), nrow = 1L)
  storage.mode(m) <- "double"
  expect_equal(suppressWarnings(cier_lazr(m))$value[[1L]], 2 / 3, tolerance = 1e-12)
})

test_that("the oracle reproduces the paper's John transition matrix", {
  # Guards the oracle itself (production exposes no T): T_12 = 9 and
  # T_21 = T_23 = T_32 = T_34 = T_43 = 8 over anchors 1..4.
  tmat <- ref_lazr_row(ref_lazr_john_sequence())$transition_matrix
  expect_equal(tmat,
               matrix(c(0L, 9L, 0L, 0L,
                        8L, 0L, 8L, 0L,
                        0L, 8L, 0L, 8L,
                        0L, 0L, 8L, 0L),
                      nrow = 4L, byrow = TRUE))
})

# ---- Independent oracle parity (1e-12) --------------------------------------

test_that("cier_lazr$value equals the oracle on random integer data", {
  x <- lazr_fixture(n = 25L, p = 30L)
  expect_equal(cier_lazr(x)$value, ref_lazr(x), tolerance = 1e-12)
})

test_that("NA transitions drop from the denominator (drop-NA, not N-1)", {
  # 10-item row 1,2,NA,4,3,2,1,2,3,4 -> the (2,NA) and (NA,4) pairs drop, leaving 7
  # valid transitions. sum(P*T) = 5, so Laz.R = 5/7 (NOT 5/9 = 5/(N-1)).
  row <- c(1, 2, NA, 4, 3, 2, 1, 2, 3, 4)
  m <- matrix(row, nrow = 1L)
  storage.mode(m) <- "double"
  out <- suppressWarnings(cier_lazr(m))
  expect_equal(out$value[[1L]], 5 / 7, tolerance = 1e-12)        # NOT 5/9
  expect_equal(out$value[[1L]], ref_lazr_row(row)$value, tolerance = 1e-12)
})

# ---- Convention pins --------------------------------------------------------

test_that("a straightliner scores exactly 1 (zero-variance chain -> 1, never NA)", {
  # The value convention: a zero-variance row has maximal predictability 1, not NA.
  # A lone respondent cannot resolve a percentile cutoff, so it is not flagged here;
  # flagging once the cutoff is resolvable is pinned by the direction test below.
  expect_warning(out <- cier_lazr(matrix(rep(3, 10L), nrow = 1L)),
                 class = "cier_warning_insufficient_items")
  expect_identical(out$value, 1)
  expect_false(is.na(out$value))
  expect_false(out$flagged[[1L]])
})

test_that("a diagonal-liner (1..5 repeated) scores 1 (deterministic chain)", {
  m <- matrix(rep_len(1:5, 20L), nrow = 1L)
  storage.mode(m) <- "double"
  expect_equal(suppressWarnings(cier_lazr(m))$value[[1L]], 1, tolerance = 1e-12)
})

test_that("values are bounded in (0, 1]", {
  v <- cier_lazr(lazr_fixture(n = 40L, p = 30L))$value
  expect_true(all(v > 0 & v <= 1 + 1e-12))
})

test_that("value is invariant to a constant integer shift (anchor-count-free)", {
  # The s-invariance that justifies the matrix-only design: predictability depends
  # on the transition structure, not the absolute anchor labels. A mutant
  # hard-coding anchors 1..max (ignoring the base) fails the 0-based recoding here.
  x <- lazr_fixture(n = 15L, p = 20L)
  expect_equal(suppressWarnings(cier_lazr(x))$value,
               suppressWarnings(cier_lazr(x + 10L))$value, tolerance = 1e-12)
  expect_equal(suppressWarnings(cier_lazr(x))$value,
               suppressWarnings(cier_lazr(x - 1L))$value, tolerance = 1e-12)
})

test_that("a stray large integer does not inflate the bin space (overflow guard)", {
  # An un-recoded numeric missing code / sentinel is a whole number, so it passes
  # validation. The kernel ranks DISTINCT observed anchors, so a single 1e5 in 1..5
  # data leaves the bin space bounded (s = 6 here), not sized by the value span
  # (s = 1e5). This covers ONE stray value (s grows by one); the section below
  # covers MANY distinct values (a slider, or an ID column with ~n distinct anchors)
  # which the chunked tabulate / ceiling abort handle. The value matches the oracle
  # (which also ranks distinct anchors), and rows without the sentinel are unchanged.
  x <- lazr_fixture(n = 10L, p = 20L)
  with_sentinel <- x
  with_sentinel[1L, 5L] <- 1e5
  out <- expect_no_error(suppressWarnings(cier_lazr(with_sentinel)))
  expect_equal(out$value, ref_lazr(with_sentinel), tolerance = 1e-12)
  expect_equal(out$value[-1L], suppressWarnings(cier_lazr(x))$value[-1L],
               tolerance = 1e-12)
})

# ---- anchor explosion -- chunked tabulation + ceiling abort -----------------
# Laz.R tabulates an s x s transition table per respondent over the s pooled
# DISTINCT anchors. A single stray value adds one anchor (above); MANY distinct
# values -- a 0-100 slider (s ~ 101), or a stray unique-integer ID / timestamp /
# un-recoded-missing column (s ~ n) -- blow up the n * s * s tabulate: ~2 GB on a
# slider at n = 10k, and integer overflow past 2^31 + an untyped crash on an ID
# column. Chunk the tabulate over respondent rows so peak allocation is bounded by
# chunk_n * s^2 (byte-identical -- per-row counts are independent), and abort with a
# typed stray-column hint when s exceeds the survey-plausible ceiling
# floor(sqrt(cell_budget)) = 2048 at the default 2^22 budget. The internal
# `cell_budget` arg is NOT exposed by the wrapper (it mirrors psychsyn's `cor_mat`
# seam), so the chunk boundary and ceiling are testable on tiny fixtures.

test_that("chunking is byte-identical across the chunk boundary (oracle parity)", {
  # Per-row independence: the chunked tabulate must reproduce the single-chunk result
  # -- and the oracle -- exactly, including the NA / abstaining rows whose
  # row-subsetting per chunk must stay aligned. With s = 5 (values 1..5),
  # cell_budget = 25 gives chunk_n = 1 (one respondent per chunk, many chunks)
  # without tripping the ceiling (5 not > floor(sqrt(25)) = 5); a huge budget is one
  # chunk.
  x <- rbind(
    lazr_fixture(n = 9L, p = 14L),
    rep(NA_real_, 14L),                          # all-NA: abstains
    c(2, 5, rep(NA_real_, 12L)),                 # one valid transition: abstains
    lazr_fixture(n = 4L, p = 14L, seed = 7L)
  )
  storage.mode(x) <- "double"
  many_chunks  <- kernel_lazr(x, cell_budget = 25)         # chunk_n = 1
  # cell_budget = 50 => chunk_n = 2, so a chunk mixes a valid + abstaining row.
  paired_chunk <- kernel_lazr(x, cell_budget = 50)
  single_chunk <- kernel_lazr(x, cell_budget = 1e9)        # one chunk
  expect_identical(many_chunks, single_chunk)              # chunk size changes nothing
  expect_identical(paired_chunk, single_chunk)             # ... incl. mixed chunks
  expect_equal(many_chunks, ref_lazr(x), tolerance = 1e-12)
  expect_identical(is.na(many_chunks), is.na(ref_lazr(x))) # abstentions aligned
})

test_that("lazr_row_chunks splits slider-scale data into multiple tiling chunks", {
  # The partition is the load-bearing artefact a "never chunks" implementation (one
  # un-chunked tabulate of n * s^2 -- the ~2 GB allocation) would skip. Pinning it
  # directly is cheaper than forcing the un-chunked path to fail (which needs
  # n * s^2 > 2^31). At the shipped default budget a slider battery (s ~ 101) MUST
  # split into > 1 chunk, and the chunks must tile 1..n exactly.
  ch <- lazr_row_chunks(500L, 101^2, 4194304)              # default budget 2^22
  expect_gt(length(ch), 1L)                                # chunking is active
  starts <- vapply(ch, `[`, integer(1L), 1L)
  stops  <- vapply(ch, `[`, integer(1L), 2L)
  expect_identical(unlist(Map(seq.int, starts, stops)), seq_len(500L)) # tiles, no gap/overlap
  expect_identical(ch[[1L]], c(1L, 411L))                  # chunk_n = floor(2^22 / 101^2) = 411
})

test_that("lazr_row_chunks is a single chunk when the matrix fits, chunk_n >= 1 at the edge", {
  expect_identical(lazr_row_chunks(30L, 25, 1e9), list(c(1L, 30L)))  # one chunk
  # ss == budget is the reachable boundary (s = 2048 at the default) -> chunk_n = 1.
  expect_identical(lazr_row_chunks(3L, 25, 25),
                   list(c(1L, 1L), c(2L, 2L), c(3L, 3L)))
})

test_that("a multi-chunk slider battery scores correctly at the default budget", {
  # 0-100 slider data: s ~ 101, so the shipped default (cell_budget 2^22,
  # chunk_n = floor(2^22 / 101^2) = 411) genuinely splits n = 500 into > 1 chunk --
  # the case that allocated ~2 GB before chunking. It must complete within the
  # bounded per-chunk allocation and match the oracle to 1e-12.
  withr::with_seed(99L, {
    x <- matrix(sample.int(101L, 500L * 16L, replace = TRUE) - 1L,   # 0..100
                nrow = 500L, ncol = 16L)
  })
  storage.mode(x) <- "double"
  # s >= 92 => chunk_n = floor(2^22 / s^2) <= 455 < 500, i.e. >= 2 chunks.
  expect_gte(length(unique(as.vector(x))), 92L)
  # Sparse anchors make most rows near-maximally predictable, so the percentile
  # cutoff saturates (a correct diagnostic); this test pins value/oracle parity, not
  # the flag, so the warning is muffled.
  expect_no_error(out <- suppressWarnings(cier_lazr(x)))
  expect_equal(out$value, ref_lazr(x), tolerance = 1e-12)
})

test_that("the ceiling abort fires exactly when s exceeds floor(sqrt(cell_budget))", {
  # The boundary, pinned via the seam on a tiny s = 5 fixture: budget 24 gives ceiling
  # floor(sqrt(24)) = 4 (5 > 4 -> abort); budget 25 gives ceiling 5 (5 not > 5 ->
  # scores). A mutant using >= or s^2 >= budget would flip one side.
  x <- lazr_fixture(n = 6L, p = 12L)                  # values 1..5 -> s = 5
  expect_error(kernel_lazr(x, cell_budget = 24), class = "cier_error_input")
  expect_no_error(suppressWarnings(kernel_lazr(x, cell_budget = 25)))
})

test_that("the ceiling abort names the stray non-item column and reports the count", {
  # The hint must point at the actionable cause (a non-item column) and report the
  # distinct-anchor count, not surface a raw allocation / tabulate error. The remedy
  # lives in an `i` bullet, so match the full cnd_message (conditionMessage carries
  # only the header) with whitespace collapsed (wrap-independent).
  x <- lazr_fixture(n = 6L, p = 12L)                  # values 1..5 -> s = 5
  err <- tryCatch(kernel_lazr(x, cell_budget = 1),
                  cier_error_input = function(e) e)
  expect_s3_class(err, "cier_error_input")
  msg <- gsub("\\s+", " ", cli::ansi_strip(rlang::cnd_message(err)))
  expect_match(msg, "distinct", fixed = TRUE)          # reports the anchor pile-up
  expect_match(msg, "item columns", fixed = TRUE)      # the actionable remedy
  expect_identical(err$data$observed, 5L)              # the reported count IS s
})

test_that("a stray unique-integer ID column aborts typed at the default budget", {
  # The repro through the PUBLIC wrapper: an ID column contributes one distinct
  # anchor per respondent, so s ~ n. At n = 2100 (> 2048) the old kernel overflowed
  # n * s * s past 2^31 and died with an untyped tabulate error; now it is a typed
  # cier_error_input. The abort is computed from s alone, before any large
  # allocation, so the 2100-row fixture stays cheap.
  ids <- seq_len(2100L) + 1000L
  x <- cbind(ids,
             item1 = rep_len(1:5, 2100L),
             item2 = rep_len(5:1, 2100L))
  storage.mode(x) <- "double"
  expect_error(cier_lazr(x), class = "cier_error_input")
})

test_that("the shipped default ceiling is 2048 distinct anchors (both sides)", {
  # One row carrying exactly K distinct anchors (1..K then back to 1). The default
  # ceiling is floor(sqrt(2^22)) = 2048: K = 2049 aborts, K = 2048 scores. n = 1
  # keeps both sides cheap (abort is pre-allocation; no-abort is a single bounded
  # chunk). The lone respondent trips the small-sample cutoff warning, muffled here
  # -- only the abort / no-abort is under test. Pins the shipped default budget,
  # complementing the seam boundary test above.
  row_k <- function(k) {
    m <- matrix(c(seq_len(k), 1), nrow = 1L)
    storage.mode(m) <- "double"
    m
  }
  expect_error(suppressWarnings(cier_lazr(row_k(2049L))), class = "cier_error_input")
  expect_no_error(out <- suppressWarnings(cier_lazr(row_k(2048L))))
  expect_true(is.finite(out$value[[1L]]))
})

# ---- Abstention / NA edges --------------------------------------------------

test_that("an all-NA row abstains and keeps rows aligned", {
  m <- rbind(
    c(1, 2, 3, 4, 3, 2, 1, 2, 3, 4),
    rep(NA_real_, 10L),                    # abstains, in the middle
    c(5, 4, 3, 2, 1, 2, 3, 4, 5, 4)
  )
  storage.mode(m) <- "double"
  out <- suppressWarnings(cier_lazr(m))
  expect_true(is.na(out$value[[2L]]))
  expect_true(is.na(out$flagged[[2L]]))
  expect_false(is.na(out$value[[1L]]))
  expect_false(is.na(out$value[[3L]]))
})

test_that("a single valid transition abstains (the < 2 rule)", {
  # Row 1 has exactly one valid transition (2,5); it must abstain rather than report
  # the degenerate 1.0. Row 2 scores so the matrix does not wholly abstain (no
  # percentile warning here).
  m <- rbind(
    c(2, 5, rep(NA_real_, 8L)),
    c(1, 2, 3, 4, 5, 4, 3, 2, 1, 2)
  )
  storage.mode(m) <- "double"
  out <- suppressWarnings(cier_lazr(m))
  expect_true(is.na(out$value[[1L]]))
  expect_true(is.na(out$flagged[[1L]]))
  expect_false(is.na(out$value[[2L]]))
})

test_that("exactly two valid transitions scores (the < 2 boundary, positive edge)", {
  # Two valid transitions (2,5),(5,1): above the threshold, so the row scores. Pins
  # the positive edge so a regression to an `n >= 1` rule is visible.
  m <- matrix(c(2, 5, 1, rep(NA_real_, 7L)), nrow = 1L)
  storage.mode(m) <- "double"
  out <- suppressWarnings(cier_lazr(m))
  expect_false(is.na(out$value[[1L]]))
  expect_true(is.finite(out$value[[1L]]))
})

test_that("a single-column matrix yields no transitions: all abstain", {
  expect_warning(
    out <- cier_lazr(matrix(c(1, 2, 3), ncol = 1L)),
    class = "cier_warning_insufficient_items"
  )
  expect_true(all(is.na(out$value)))
})

test_that("a wholly abstaining matrix warns and flags nobody", {
  expect_warning(
    out <- cier_lazr(matrix(NA_real_, nrow = 3L, ncol = 10L)),
    class = "cier_warning_insufficient_items"
  )
  expect_true(all(is.na(out$value)))
  expect_true(is.na(out$cutoff))
  expect_true(all(is.na(out$flagged)))
})

# ---- Input validation -------------------------------------------------------

test_that("a non-matrix / non-numeric payload is a typed input error", {
  expect_error(cier_lazr(1:10), class = "cier_error_input")
  expect_error(cier_lazr(matrix(letters[1:8], nrow = 2L)),
               class = "cier_error_input")
  expect_error(cier_lazr(NULL), class = "cier_error_input")
})

test_that("NaN / infinite values are typed input errors", {
  expect_error(cier_lazr(matrix(c(1, 2, 3, Inf, 1, 2), nrow = 2L)),
               class = "cier_error_input")
  expect_error(cier_lazr(matrix(c(1, 2, 3, NaN, 1, 2), nrow = 2L)),
               class = "cier_error_input")
})

test_that("non-integer responses are a typed input error", {
  expect_error(cier_lazr(matrix(c(1, 2, 2.5, 4, 3, 2), nrow = 2L)),
               class = "cier_error_input")
})

test_that("bad fpr values are typed input errors", {
  x <- lazr_fixture(n = 6L)
  expect_error(cier_lazr(x, fpr = 0), class = "cier_error_input")
  expect_error(cier_lazr(x, fpr = 1), class = "cier_error_input")
  expect_error(cier_lazr(x, fpr = -0.1), class = "cier_error_input")
  expect_error(cier_lazr(x, fpr = NA_real_), class = "cier_error_input")
  expect_error(cier_lazr(x, fpr = c(0.05, 0.1)), class = "cier_error_input")
  expect_error(cier_lazr(x, fpr = "x"), class = "cier_error_input")
})

test_that("invalid literal cutoff values are typed input errors", {
  x <- lazr_fixture(n = 6L)
  expect_error(cier_lazr(x, cutoff = 1.5), class = "cier_error_input")
  expect_error(cier_lazr(x, cutoff = -0.1), class = "cier_error_input")
  expect_error(cier_lazr(x, cutoff = NA_real_), class = "cier_error_input")
  expect_error(cier_lazr(x, cutoff = c(0.1, 0.2)), class = "cier_error_input")
  expect_error(cier_lazr(x, cutoff = "x"), class = "cier_error_input")
})

test_that("supplying both fpr and cutoff is a typed input error", {
  expect_error(cier_lazr(lazr_fixture(n = 6L), fpr = 0.1, cutoff = 0.5),
               class = "cier_error_input")
})

# ---- Cutoff: default, fpr override, NO double-flip, direction ---------------

test_that("default cutoff is the upper-tail 95th percentile (NO double-flip)", {
  x <- lazr_fixture(n = 60L, p = 20L)
  out <- cier_lazr(x)
  expect_equal(out$cutoff,
               as.numeric(stats::quantile(out$value, 0.95,
                                          names = FALSE, type = 7L)),
               tolerance = 1e-12)
})

test_that("the fpr argument moves the percentile target", {
  x <- lazr_fixture(n = 60L, p = 20L)
  out <- cier_lazr(x, fpr = 0.10)
  expect_equal(out$cutoff,
               as.numeric(stats::quantile(out$value, 0.90,
                                          names = FALSE, type = 7L)),
               tolerance = 1e-12)
})

test_that("direction is upper: high-predictability rows flag, low ones do not", {
  x <- rbind(rep(3, 20L), lazr_fixture(n = 20L, p = 20L))   # constant -> value 1
  storage.mode(x) <- "double"
  out <- cier_lazr(x)
  expect_true(out$flagged[[1L]])                            # constant, value 1
  expect_false(out$flagged[[which.min(out$value)]])         # least predictable
  expect_identical(out$flagged,
                   !is.na(out$value) & out$value >= out$cutoff)
})

test_that("a literal cutoff passes through and flags via the upper direction", {
  x <- rbind(rep(3, 20L), lazr_fixture(n = 10L, p = 20L))
  storage.mode(x) <- "double"
  out <- cier_lazr(x, cutoff = 0.95)
  expect_identical(out$cutoff, 0.95)
  expect_true(out$flagged[[1L]])                            # value 1 >= 0.95
  expect_identical(out$flagged, !is.na(out$value) & out$value >= 0.95)
})

# ---- Cutoff: the paper-faithful Kneedle elbow (kneedle = TRUE) --------------
# Biemann et al.'s studies flag the top 5% (= the fpr default); the Kneedle elbow
# (Satopaa et al. 2011) is the sample-specific cutoff from their companion app.
# cier ships it as an opt-in `kneedle = TRUE`, oracle-only trust.

source(test_path("..", "reference", "ref-kneedle-satopaa-2011.R"))

test_that("kneedle = TRUE sets the cutoff to the Satopaa elbow (oracle, tol 0)", {
  x <- lazr_fixture(n = 70L, p = 18L)
  out <- cier_lazr(x, kneedle = TRUE)
  # End-to-end paper-faithfulness: the cutoff IS the convex/increasing elbow of the
  # observed scores. The oracle drops the NA-abstaining rows, as the resolver does.
  # Also proves kneedle is not aliased to the percentile default.
  expect_identical(out$cutoff,
                   ref_kneedle(out$value, "convex", "increasing")$value)
  expect_identical(out$direction, "upper")
  expect_s3_class(out, "cier_index")
})

test_that("kneedle = TRUE flags the predictable high tail, not the bulk", {
  # Five constant rows (each a straightliner -> Laz.R = 1) over a careful bulk: the
  # elbow sits at the bulk/spike boundary, so every straightliner flags while the
  # least predictable respondent does not.
  x <- rbind(matrix(c(2, 4, 1, 5, 3), nrow = 5L, ncol = 18L),
             lazr_fixture(n = 35L, p = 18L))
  storage.mode(x) <- "double"
  out <- cier_lazr(x, kneedle = TRUE)
  expect_true(all(out$flagged[1:5]))                       # the spike flags
  expect_false(out$flagged[[which.min(out$value)]])        # the bulk floor does not
  expect_false(all(out$flagged))                           # a tail, not everyone
})

test_that("kneedle = TRUE is deterministic (no RNG)", {
  x <- lazr_fixture(n = 50L, p = 16L)
  expect_identical(cier_lazr(x, kneedle = TRUE)$cutoff,
                   cier_lazr(x, kneedle = TRUE)$cutoff)
})

test_that("kneedle is mutually exclusive with fpr and with cutoff", {
  x <- lazr_fixture(n = 12L)
  expect_error(cier_lazr(x, kneedle = TRUE, fpr = 0.1),
               class = "cier_error_input")
  expect_error(cier_lazr(x, kneedle = TRUE, cutoff = 0.5),
               class = "cier_error_input")
})

test_that("a non-flag kneedle argument is a typed input error", {
  x <- lazr_fixture(n = 12L)
  expect_error(cier_lazr(x, kneedle = "yes"), class = "cier_error_input")
  expect_error(cier_lazr(x, kneedle = NA), class = "cier_error_input")
  expect_error(cier_lazr(x, kneedle = c(TRUE, FALSE)), class = "cier_error_input")
})

test_that("kneedle abstains (NA cutoff + warning, flags nobody) below three scores", {
  # Two scoring rows plus all-NA rows: only two finite Laz.R values remain, so the
  # elbow is undefined and the cutoff abstains -- the percentile path's contract,
  # reached through the kneedle resolver.
  x <- rbind(c(1, 2, 3, 4, 5, 4, 3, 2, 1, 2),
             c(5, 4, 3, 2, 1, 2, 3, 4, 5, 4),
             matrix(NA_real_, nrow = 3L, ncol = 10L))
  storage.mode(x) <- "double"
  expect_warning(out <- cier_lazr(x, kneedle = TRUE),
                 class = "cier_warning_insufficient_items")
  expect_true(is.na(out$cutoff))
  expect_false(any(out$flagged, na.rm = TRUE))     # an NA cutoff flags nobody
  expect_false(out$flagged[[1L]])                  # a scored row -> FALSE, not NA
  expect_true(all(is.na(out$flagged[3:5])))        # the all-NA rows abstain
})

test_that("kneedle = TRUE warns but still resolves when the elbow flags a majority", {
  # A low-quality panel: a careless MAJORITY (25 straightliners -> Laz.R = 1) over a
  # small careful minority. The convex elbow sits at/below the careful bulk, so the
  # Kneedle cutoff flags a majority. End-to-end this must surface the saturation
  # diagnostic (not silently flag 100%) AND still resolve the cutoff -- a
  # true-majority panel is not an abstention.
  straight <- matrix(rep(rep_len(1:5, 25L), 18L), nrow = 25L)   # 25 constant rows
  x <- rbind(straight, lazr_fixture(n = 8L, p = 18L))
  storage.mode(x) <- "double"
  expect_warning(out <- cier_lazr(x, kneedle = TRUE),
                 class = "cier_warning_saturated_cutoff")
  expect_false(is.na(out$cutoff))                  # resolved, NOT abstained
  expect_gt(mean(out$flagged, na.rm = TRUE), 0.5)  # a majority flagged
})

test_that("kneedle = TRUE does not warn on a healthy careless minority", {
  # The companion to the majority case: a small careless spike (5 straightliners over
  # a careful bulk, ~15% flagged) flags a minority, so no saturation warning fires.
  x <- rbind(matrix(c(2, 4, 1, 5, 3), nrow = 5L, ncol = 18L),
             lazr_fixture(n = 35L, p = 18L))
  storage.mode(x) <- "double"
  expect_no_warning(cier_lazr(x, kneedle = TRUE),
                    class = "cier_warning_saturated_cutoff")
})

# ---- print snapshot (locked; reuses the shared cier_index print) ------------

test_that("print renders the locked cli summary (upper direction)", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    expect_snapshot(print(cier_lazr(lazr_fixture(n = 30L, p = 12L))))
  })
})

test_that("print reports abstaining respondents on their own line", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    x <- lazr_fixture(n = 29L, p = 12L)
    x <- rbind(x, rep(NA_real_, 12L))     # one abstaining respondent
    expect_snapshot(print(cier_lazr(x)))
  })
})
