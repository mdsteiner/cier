# Tests for the shared cier_index object: cutoff-provenance metadata, the cutoff
# descriptor + ordinal rendering helpers, the new_cier_index schema invariants, and
# the summary.cier_index method.
#
# Trust model: the provenance label + rate every wrapper records is re-derived
# INDEPENDENTLY here from the provenance spec (cutoff_method vocabulary + its rate
# source), never read back from the resolver -- a mutant tagging the wrong method
# or rate must fail. The descriptor strings are pinned by a hand-written oracle;
# the print / summary snapshots then lock only the format.

# A reproducible continuous matrix (no tie mass) so the percentile cutoff resolves
# cleanly and the small-sample guard (>= 20 scored) is satisfied.
prov_matrix <- function(n = 30L, p = 6L, seed = 1L) {
  withr::with_seed(seed, matrix(stats::runif(n * p, 1, 5), nrow = n, ncol = p))
}

# 3-scale item metadata (>= 3 scales avoids the two-scale degeneracy warning) for
# the upper-direction percentile index even-odd.
prov_items <- function() {
  data.frame(scale = rep(c("a", "b", "c"), each = 3L),
             reverse_keyed = rep(c(FALSE, TRUE, FALSE), each = 3L),
             max = 5L, stringsAsFactors = FALSE)
}

# ---- ordinal() oracle -------------------------------------------------------

test_that("ordinal() forms English ordinals incl. the 11-13 exception", {
  expect_identical(ordinal(1L), "1st")
  expect_identical(ordinal(2L), "2nd")
  expect_identical(ordinal(3L), "3rd")
  expect_identical(ordinal(4L), "4th")
  expect_identical(ordinal(5L), "5th")
  expect_identical(ordinal(10L), "10th")
  expect_identical(ordinal(11L), "11th")
  expect_identical(ordinal(12L), "12th")
  expect_identical(ordinal(13L), "13th")
  expect_identical(ordinal(21L), "21st")
  expect_identical(ordinal(95L), "95th")
  expect_identical(ordinal(99L), "99th")
})

# ---- cutoff_descriptor() oracle (position + rate phrase per kind) -----------

test_that("cutoff_descriptor re-derives each provenance kind's phrase", {
  # percentile: the quantile position flips on direction (lower = fpr tail, upper =
  # 1 - fpr tail); the rate phrase always names fpr.
  expect_identical(cutoff_descriptor("percentile", 0.05, "lower"),
                   list(position = "5th sample percentile",
                        rate_phrase = "fpr = 0.05"))
  expect_identical(cutoff_descriptor("percentile", 0.05, "upper"),
                   list(position = "95th sample percentile",
                        rate_phrase = "fpr = 0.05"))
  expect_identical(cutoff_descriptor("percentile", 0.1, "lower"),
                   list(position = "10th sample percentile",
                        rate_phrase = "fpr = 0.1"))
  # A non-integer percentile position must ROUND, not floor/truncate: fpr = 0.125 is
  # exactly representable, so the upper tail is 100*(1-0.125) = 87.5 -> 88 (a floor
  # mutant gives 87). Pins the rounding rule away from the boundary.
  expect_identical(cutoff_descriptor("percentile", 0.125, "upper")$position,
                   "88th sample percentile")
  expect_identical(cutoff_descriptor("chisq", 0.001, "upper"),
                   list(position = "chi-square tail",
                        rate_phrase = "alpha = 0.001"))
  expect_identical(cutoff_descriptor("mc_null", 0.05, "upper"),
                   list(position = "Monte-Carlo null",
                        rate_phrase = "nominal = 0.05"))
  expect_identical(cutoff_descriptor("fixed_fraction", 0.5, "upper"),
                   list(position = "fixed fraction 0.5 of the item count",
                        rate_phrase = ""))
  expect_identical(cutoff_descriptor("median_relative", 0.5, "lower"),
                   list(position = "0.5 x the sample median",
                        rate_phrase = ""))
  expect_identical(cutoff_descriptor("fixed_count", NA_real_, "upper"),
                   list(position = "fixed count", rate_phrase = ""))
  expect_identical(cutoff_descriptor("kneedle", NA_real_, "upper"),
                   list(position = "Kneedle elbow (parameter-free)",
                        rate_phrase = ""))
  expect_identical(cutoff_descriptor("literal", NA_real_, "lower"),
                   list(position = "user-supplied threshold", rate_phrase = ""))
})

# ---- new_cier_index schema + invariants -------------------------------------

test_that("new_cier_index carries the provenance fields in the pinned schema", {
  ix <- new_cier_index(c(1, 2, 3), c(TRUE, FALSE, FALSE), "cier_irv", 2, "lower",
                       cutoff_method = "percentile", cutoff_rate = 0.05)
  expect_identical(names(ix),
                   c("value", "flagged", "method", "cutoff", "direction",
                     "cutoff_method", "cutoff_rate"))
  expect_identical(ix$cutoff_method, "percentile")
  expect_identical(ix$cutoff_rate, 0.05)
})

test_that("new_cier_index defaults the provenance fields to NA (back-compat)", {
  ix <- new_cier_index(c(1, 2, 3), c(TRUE, FALSE, FALSE), "cier_irv", 2, "lower")
  expect_identical(ix$cutoff_method, NA_character_)
  expect_identical(ix$cutoff_rate, NA_real_)
})

test_that("new_cier_index rejects an out-of-vocabulary cutoff_method", {
  expect_error(
    new_cier_index(c(1, 2), c(TRUE, FALSE), "cier_irv", 2, "lower",
                   cutoff_method = "made_up", cutoff_rate = NA_real_),
    class = "cier_error_state"
  )
})

test_that("new_cier_index rejects a non-scalar cutoff_rate", {
  expect_error(
    new_cier_index(c(1, 2), c(TRUE, FALSE), "cier_irv", 2, "lower",
                   cutoff_method = "percentile", cutoff_rate = c(0.05, 0.1)),
    class = "cier_error_state"
  )
})

test_that("new_cier_index rejects a non-string, non-NA cutoff_method", {
  expect_error(
    new_cier_index(c(1, 2), c(TRUE, FALSE), "cier_irv", 2, "lower",
                   cutoff_method = 5L, cutoff_rate = NA_real_),
    class = "cier_error_state"
  )
})

test_that("new_cier_index's contract error names the provenance fields", {
  # The guard validates cutoff_method / cutoff_rate, so its message must name them
  # (the old stopped at direction) and attach the offending values. The remedy lives
  # in an `i` bullet, so match the whitespace-collapsed cnd_message (conditionMessage
  # carries only the header).
  err <- tryCatch(
    new_cier_index(c(1, 2), c(TRUE, FALSE), "cier_irv", 2, "lower",
                   cutoff_method = "made_up", cutoff_rate = NA_real_),
    cier_error_state = function(e) e
  )
  expect_s3_class(err, "cier_error_state")
  msg <- gsub("\\s+", " ", cli::ansi_strip(rlang::cnd_message(err)))
  expect_match(msg, "cutoff_method", fixed = TRUE)
  expect_match(msg, "cutoff_rate", fixed = TRUE)
  expect_identical(err$data$cutoff_method, "made_up")
  expect_identical(err$data$cutoff_rate, NA_real_)
})

test_that("as.data.frame.cier_index drops respondent names from the row index", {
  # Even when a value vector arrives NAMED, the frame must carry the plain integer
  # row index so the positional contract holds (cier_flagged_cases() /
  # as.data.frame.cier_screen()). Build the object with a NAMED value directly,
  # bypassing kernels, to pin the defensive rownames(out) <- NULL on its own.
  v <- c(a = 1, b = 2, c = 3)
  f <- c(a = TRUE, b = FALSE, c = FALSE)
  ix <- new_cier_index(v, f, "cier_irv", 2, "lower")
  df <- as.data.frame(ix)
  # The tidy frame is exactly the two per-respondent columns (this is the single,
  # index-independent pin of the as.data.frame.cier_index column contract; the
  # per-index schema tests assert only the cier_index object, not its frame).
  expect_s3_class(df, "data.frame")
  expect_identical(names(df), c("value", "flagged"))
  # The row index is THE discriminator: without the fix data.frame() adopts the
  # value vector's names (a, b, c) as row labels; with it the index is 1..n.
  expect_identical(rownames(df), as.character(seq_len(3L)))
  # data.frame() strips the column's own names regardless, so these confirm only
  # that rownames(out) <- NULL leaves the column values/order intact (not the
  # row-name strip itself).
  expect_identical(df$value, unname(v))
  expect_identical(df$flagged, unname(f))
})

# ---- Provenance oracle: every wrapper tags the right (method, rate) ----------

test_that("percentile indices record method = percentile and the fpr rate", {
  x <- prov_matrix()
  irv <- cier_irv(x)                                   # lower tail
  expect_identical(irv$cutoff_method, "percentile")
  expect_identical(irv$cutoff_rate, 0.05)
  expect_identical(cier_irv(x, fpr = 0.10)$cutoff_rate, 0.10)
  x9 <- prov_matrix(p = 9L)
  eo <- suppressWarnings(cier_even_odd(x9, prov_items()))  # upper tail
  expect_identical(eo$cutoff_method, "percentile")
  expect_identical(eo$cutoff_rate, 0.05)
})

test_that("a literal cutoff override records method = literal, rate NA", {
  x <- prov_matrix()
  lit <- cier_irv(x, cutoff = 0.5)
  expect_identical(lit$cutoff_method, "literal")
  expect_identical(lit$cutoff_rate, NA_real_)
})

test_that("provenance is the resolution STRATEGY, recorded even when it abstains", {
  # A single-column matrix abstains for every row, so the percentile cutoff is NA.
  # The provenance still records the attempted strategy (percentile, the default
  # fpr) -- a "no value resolved -> no provenance" mutant (cutoff_method dropped to
  # NA) must fail.
  out <- suppressWarnings(cier_irv(matrix(c(3, 4, 5), ncol = 1L)))
  expect_identical(out$cutoff, NA_real_)
  expect_identical(out$cutoff_method, "percentile")
  expect_identical(out$cutoff_rate, 0.05)
})

test_that("longstring records a fixed fraction of the item count", {
  x <- prov_matrix()
  ls_default <- cier_longstring(x)
  expect_identical(ls_default$cutoff_method, "fixed_fraction")
  expect_identical(ls_default$cutoff_rate, 0.5)
  expect_identical(cier_longstring(x, frac = 0.4)$cutoff_rate, 0.4)
  expect_identical(cier_longstring(x, cutoff = 5)$cutoff_method, "literal")
})

test_that("mahalanobis records the chi-square tail with the alpha rate", {
  x <- prov_matrix()
  mh <- cier_mahalanobis(x)
  expect_identical(mh$cutoff_method, "chisq")
  expect_identical(mh$cutoff_rate, 0.001)
  expect_identical(cier_mahalanobis(x, alpha = 0.01)$cutoff_rate, 0.01)
  # The literal-override branch of index_provenance() must win here too.
  expect_identical(cier_mahalanobis(x, cutoff = 60)$cutoff_method, "literal")
})

test_that("the timing family records its own provenance kinds", {
  secs <- prov_matrix(n = 30L, p = 1L)[, 1]
  tt <- cier_total_time(secs)
  expect_identical(tt$cutoff_method, "percentile")
  expect_identical(tt$cutoff_rate, 0.05)
  med <- cier_total_time(secs, frac_median = 0.5)
  expect_identical(med$cutoff_method, "median_relative")
  expect_identical(med$cutoff_rate, 0.5)
  expect_identical(cier_total_time(secs, cutoff = 100)$cutoff_method, "literal")
  expect_identical(cier_total_time(secs, fpr = 0.10)$cutoff_rate, 0.10)
  # page time: registry default is a literal count of 1 (fixed_count); frac is a
  # fraction of the page count (fixed_fraction); a literal cutoff is "literal".
  pages <- matrix(c(10, 2, 4, 2), nrow = 2L)
  pt <- cier_page_time(pages, c(2L, 2L))
  expect_identical(pt$cutoff_method, "fixed_count")
  expect_identical(pt$cutoff_rate, NA_real_)
  pt_frac <- cier_page_time(pages, c(2L, 2L), frac = 0.5)
  expect_identical(pt_frac$cutoff_method, "fixed_fraction")
  expect_identical(pt_frac$cutoff_rate, 0.5)
  expect_identical(cier_page_time(pages, c(2L, 2L), cutoff = 1)$cutoff_method,
                   "literal")
})

test_that("attention records a fixed count, lazr records percentile/kneedle", {
  checks <- matrix(c(1, 0, 3, 0, 2, 5, 4, 7), nrow = 4L, byrow = TRUE)
  att <- cier_attention(checks, pass = list(c(1, 2), 0))
  expect_identical(att$cutoff_method, "fixed_count")
  expect_identical(att$cutoff_rate, NA_real_)
  expect_identical(cier_attention(checks, pass = list(c(1, 2), 0),
                                  cutoff = 2)$cutoff_method, "literal")
  x <- matrix(withr::with_seed(3L, sample.int(5L, 30L * 8L, replace = TRUE)),
              nrow = 30L)
  storage.mode(x) <- "double"
  lz <- cier_lazr(x)
  expect_identical(lz$cutoff_method, "percentile")
  expect_identical(lz$cutoff_rate, 0.05)
  kn <- suppressWarnings(cier_lazr(x, kneedle = TRUE))
  expect_identical(kn$cutoff_method, "kneedle")
  expect_identical(kn$cutoff_rate, NA_real_)
})

test_that("gnormed records the Monte-Carlo null nominal", {
  x <- prov_matrix(n = 40L, p = 8L)
  x <- round(x)                       # whole-number categories spanning 1..5
  x[1L, 1L] <- 1
  x[2L, 1L] <- 5
  it <- data.frame(scale = rep(c("a", "b"), each = 4L), reverse_keyed = FALSE,
                   max = 5L, stringsAsFactors = FALSE)
  gn <- suppressWarnings(cier_gnormed(x, it, seed = 1))
  expect_identical(gn$cutoff_method, "mc_null")
  expect_identical(gn$cutoff_rate, 0.05)
})

# ---- print: the Cutoff method line (values pinned; snapshot locks format) ----

test_that("print adds a Cutoff method line naming the provenance", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    out <- capture.output(print(cier_irv(prov_matrix())))
  })
  line <- out[grepl("Cutoff method:", out, fixed = TRUE)]
  expect_length(line, 1L)
  expect_match(line, "5th sample percentile (fpr = 0.05)", fixed = TRUE)
})

test_that("print omits the Cutoff method line when provenance is unknown (NA)", {
  ix <- new_cier_index(c(1, 2, 3), c(TRUE, FALSE, FALSE), "cier_irv", 2, "lower")
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    out <- capture.output(print(ix))
  })
  expect_false(any(grepl("Cutoff method:", out, fixed = TRUE)))
})

test_that("print renders the locked cli summary with the Cutoff method line", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    expect_snapshot(print(cier_irv(prov_matrix())))
  })
})

# ---- summary.cier_index -----------------------------------------------------

test_that("summary.cier_index returns the object invisibly", {
  out <- cier_irv(prov_matrix())
  expect_identical(withVisible(summary(out))$visible, FALSE)
  expect_identical(suppressMessages(summary(out)), out)
})

test_that("summary.cier_index adds quartiles + position that print does NOT show", {
  x <- prov_matrix()
  out <- cier_irv(x)
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    s <- capture.output(summary(out))
    p <- capture.output(print(out))
  })
  # The most plausible wrong implementation just delegates to print(); summary MUST
  # differ (quartiles + cutoff position are extra).
  expect_false(identical(s, p))
  joined <- paste(s, collapse = " ")
  expect_match(joined, "Score quartiles")
  expect_match(joined, "Scored 30 of 30")
  # The five-number summary is the quantile of the SCORED values (independent
  # re-derivation), formatted as the method does (signif 3), so a mutant summarising
  # the raw/flagged vector would diverge at min/max.
  fivenum_ref <- stats::quantile(out$value[!is.na(out$value)],
                                 c(0, .25, .5, .75, 1), names = FALSE)
  expect_match(joined, format(signif(fivenum_ref[[1L]], 3)), fixed = TRUE)
  expect_match(joined, format(signif(fivenum_ref[[5L]], 3)), fixed = TRUE)
  # The cutoff-position line is direction-aware (lower -> "at or below") and
  # value-pinned to the scored flag count, not a bare noun phrase.
  n_flagged <- sum(out$flagged, na.rm = TRUE)
  expect_match(joined, "at or below the cutoff", ignore.case = TRUE)
  expect_match(joined, sprintf("%d of 30", n_flagged), fixed = TRUE)
})

test_that("summary.cier_index says 'at or above' for an upper-direction index", {
  x9 <- prov_matrix(p = 9L)
  out <- suppressWarnings(cier_even_odd(x9, prov_items()))   # upper tail
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    joined <- paste(capture.output(summary(out)), collapse = " ")
  })
  expect_match(joined, "at or above the cutoff", ignore.case = TRUE)
  expect_false(grepl("at or below the cutoff", joined, ignore.case = TRUE))
})

test_that("summary.cier_index handles an all-abstaining / unresolved cutoff", {
  # A single-column matrix abstains for every row; the percentile cutoff is NA. No
  # scored values -> no quartiles and no position line, but the scored/abstain count
  # is still reported (something print's wording does not supply verbatim).
  out <- suppressWarnings(cier_irv(matrix(c(3, 4, 5), ncol = 1L)))
  expect_identical(out$cutoff, NA_real_)
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    txt <- capture.output(suppressMessages(summary(out)))
  })
  joined <- paste(txt, collapse = " ")
  expect_match(joined, "Scored 0 of 3", fixed = TRUE)
  # The cutoff-position line is the only place "at or below/above the cutoff"
  # appears; with no scored values it must be absent (don't match "the cutoff"
  # broadly, which can appear in unrelated prose).
  expect_false(grepl("at or below the cutoff", joined, ignore.case = TRUE))
  expect_false(grepl("at or above the cutoff", joined, ignore.case = TRUE))
  expect_false(grepl("Score quartiles", joined, fixed = TRUE))
})

test_that("summary.cier_index renders the locked snapshot (scored + abstain)", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    expect_snapshot(summary(cier_irv(prov_matrix())))
    withr::with_seed(11, {
      x <- matrix(sample.int(5L, 29L * 6L, replace = TRUE), nrow = 29L)
    })
    storage.mode(x) <- "double"
    x <- rbind(x, rep(NA_real_, 6L))
    expect_snapshot(summary(suppressWarnings(cier_irv(x))))
  })
})
