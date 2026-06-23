# check_responses(): coerces a data.frame / tibble to a numeric matrix and rejects anything not
# numeric-and-finite-or-NA. All failures are typed cier_error_input (asserted by class, not text).

test_that("a numeric matrix passes through unchanged in value", {
  m <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2L)
  out <- check_responses(m)
  expect_true(is.matrix(out) && is.numeric(out))
  expect_equal(unname(out), matrix(c(1, 2, 3, 4, 5, 6), nrow = 2L))
})

test_that("a data.frame of numeric columns is coerced, values preserved", {
  df <- data.frame(a = c(1, 2), b = c(3, 4), c = c(5, 6))
  out <- check_responses(df)
  expect_true(is.matrix(out) && is.numeric(out))
  expect_equal(unname(out), matrix(c(1, 2, 3, 4, 5, 6), nrow = 2L))
})

test_that("NA is an allowed missing marker", {
  m <- matrix(c(1, NA, 3, 4), nrow = 2L)
  expect_no_error(check_responses(m))
})

test_that("NaN is rejected", {
  m <- matrix(c(1, NaN, 3, 4), nrow = 2L)
  expect_error(check_responses(m), class = "cier_error_input")
})

test_that("Inf and -Inf are rejected", {
  expect_error(check_responses(matrix(c(1, Inf, 3, 4), nrow = 2L)),
               class = "cier_error_input")
  expect_error(check_responses(matrix(c(1, -Inf, 3, 4), nrow = 2L)),
               class = "cier_error_input")
})

test_that("a non-numeric (character) column is rejected", {
  df <- data.frame(a = c("x", "y"), b = c(1, 2), stringsAsFactors = FALSE)
  expect_error(check_responses(df), class = "cier_error_input")
})

test_that("a factor column is rejected", {
  df <- data.frame(a = factor(c("x", "y")), b = c(1, 2))
  expect_error(check_responses(df), class = "cier_error_input")
})

test_that("a logical matrix is rejected (not numeric)", {
  expect_error(check_responses(matrix(c(TRUE, FALSE, TRUE, FALSE), nrow = 2L)),
               class = "cier_error_input")
})

test_that("zero-row and zero-column inputs are rejected", {
  expect_error(check_responses(matrix(numeric(0), nrow = 0L, ncol = 3L)),
               class = "cier_error_input")
  expect_error(check_responses(matrix(numeric(0), nrow = 3L, ncol = 0L)),
               class = "cier_error_input")
})

test_that("a bare numeric vector is rejected, not silently shaped into a column", {
  # as.matrix(c(3,3,3,3)) would yield a 4x1 matrix (4 'respondents' of 1 item); a vector is
  # ambiguous, so require an explicit 2-D matrix / data.frame.
  expect_error(check_responses(c(3, 3, 3, 3)), class = "cier_error_input")
})

test_that("a higher-dimensional array is rejected", {
  expect_error(check_responses(array(1:24, dim = c(2L, 3L, 4L))),
               class = "cier_error_input")
})

# check_items(): validates the per-item `items` frame the split-half family uses: `scale`
# (>= min_scales distinct), optional logical `reverse_keyed` (defaults all-FALSE), and `max` (the
# largest response option) -- required (integer >= min + 1, non-NA) only on reverse-keyed items.
# All failures are typed cier_error_input.

test_that("check_items returns normalized scale / reverse_keyed / max", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, TRUE),
                   max = 5L)
  out <- check_items(it, n_items = 4L)
  expect_identical(out$scale, c("A", "A", "B", "B"))
  expect_identical(out$reverse_keyed, c(FALSE, TRUE, FALSE, TRUE))
  expect_identical(out$max, rep(5L, 4L))
})

test_that("check_items defaults reverse_keyed to all-FALSE", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L))
  out <- check_items(it, n_items = 4L)
  expect_identical(out$reverse_keyed, rep(FALSE, 4L))
})

test_that("check_items coerces a factor scale column to character", {
  it <- data.frame(scale = factor(rep(c("A", "B"), each = 2L)),
                   reverse_keyed = FALSE)
  out <- check_items(it, n_items = 4L)
  expect_type(out$scale, "character")
  expect_identical(out$scale, c("A", "A", "B", "B"))
})

test_that("check_items requires the items row-count to equal n_items", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L))
  expect_error(check_items(it, n_items = 5L), class = "cier_error_input")
})

test_that("check_items requires at least min_scales distinct scales", {
  it <- data.frame(scale = rep("A", 4L))
  expect_error(check_items(it, n_items = 4L, min_scales = 2L),
               class = "cier_error_input")
})

test_that("check_items requires a scale column", {
  it <- data.frame(reverse_keyed = rep(FALSE, 4L))
  expect_error(check_items(it, n_items = 4L), class = "cier_error_input")
})

test_that("check_items rejects a non-data-frame items argument", {
  expect_error(check_items(list(scale = c("A", "B")), n_items = 2L),
               class = "cier_error_input")
})

test_that("check_items rejects a non-logical reverse_keyed column", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(0L, 1L, 0L, 1L))
  expect_error(check_items(it, n_items = 4L), class = "cier_error_input")
})

test_that("check_items rejects NA in reverse_keyed", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, NA, FALSE, TRUE))
  expect_error(check_items(it, n_items = 4L), class = "cier_error_input")
})

test_that("check_items requires max when an item is reverse-keyed", {
  # No max column but a reverse item -> cannot reverse-score -> error.
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, FALSE))
  expect_error(check_items(it, n_items = 4L), class = "cier_error_input")
})

test_that("check_items rejects NA max on a reverse item", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, FALSE),
                   max = c(5L, NA, 5L, 5L))
  expect_error(check_items(it, n_items = 4L), class = "cier_error_input")
})

test_that("check_items allows NA max on non-reverse items", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, FALSE),
                   max = c(NA, 5L, NA, 3L))
  expect_no_error(check_items(it, n_items = 4L))
})

test_that("check_items rejects max at or below min on a reverse item", {
  # max == min is a one-option "scale" (nothing to reflect); with the default min = 1, max = 1 is
  # the malformed declaration the old categories >= 2 bound caught.
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, FALSE),
                   max = c(5L, 1L, 5L, 5L))
  expect_error(check_items(it, n_items = 4L), class = "cier_error_input")
  # Same violation with an explicit min ABOVE the absolute-2 line: max == min == 2 must be rejected
  # too, so a validator testing `max >= 2` (ignoring min) is wrong.
  it2 <- data.frame(scale = rep(c("A", "B"), each = 2L),
                    reverse_keyed = c(FALSE, TRUE, FALSE, FALSE),
                    max = c(5L, 2L, 5L, 5L), min = c(1L, 2L, 1L, 1L))
  expect_error(check_items(it2, n_items = 4L), class = "cier_error_input")
})

test_that("check_items accepts a two-option scale (max == min + 1)", {
  # The smallest valid scale: two response options. With min = 0 this is a 0/1 item; a validator
  # demanding max >= min + 2 (or max >= 2 regardless of min) would wrongly reject it.
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, TRUE),
                   max = 1L, min = 0L)
  expect_no_error(check_items(it, n_items = 4L))
})

test_that("check_items rejects non-integer max on a reverse item", {
  # The reflection (min + max) - x assumes whole-number response options; a
  # fractional value is a malformed item definition.
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, FALSE),
                   max = c(5, 2.5, 5, 5))
  expect_error(check_items(it, n_items = 4L), class = "cier_error_input")
})

test_that("check_items rejects a non-finite (Inf) max on a reverse item", {
  # Inf passes is.numeric / >= min + 1 / == round(); only the is.finite guard rejects it. Without
  # it the item would reflect to (min + Inf) - x = Inf and silently poison the reverse columns.
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, FALSE),
                   max = c(5, Inf, 5, 5))
  expect_error(check_items(it, n_items = 4L), class = "cier_error_input")
})

test_that("check_items does not require max on a forward-keyed battery", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L), reverse_keyed = FALSE)
  expect_no_error(check_items(it, n_items = 4L))
})

# `min` -- the response-scale base. Optional; defaults to 1 (1..max coding). When supplied it
# generalises the reverse-keying reflection to (min + max) - x so a 0-based (or bipolar) scale
# reflects onto itself.

test_that("check_items defaults min to all-1 when the column is absent", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, TRUE), max = 5L)
  out <- check_items(it, n_items = 4L)
  expect_identical(out$min, rep(1L, 4L))
})

test_that("check_items returns a supplied min column", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, TRUE),
                   max = 5L, min = 0L)
  out <- check_items(it, n_items = 4L)
  expect_identical(out$min, rep(0L, 4L))
})

test_that("check_items allows a zero, negative, or bipolar min base on reverse items", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, TRUE),
                   max = 5L, min = c(0L, -2L, 0L, -2L))
  expect_no_error(check_items(it, n_items = 4L))
})

test_that("check_items rejects a non-finite min on a reverse item", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, FALSE),
                   max = 5L, min = c(1, Inf, 1, 1))
  expect_error(check_items(it, n_items = 4L), class = "cier_error_input")
})

test_that("check_items rejects a non-integer min on a reverse item", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, FALSE),
                   max = 5L, min = c(1, 0.5, 1, 1))
  expect_error(check_items(it, n_items = 4L), class = "cier_error_input")
})

test_that("check_items rejects NA min on a reverse item", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, FALSE),
                   max = 5L, min = c(1L, NA, 1L, 1L))
  expect_error(check_items(it, n_items = 4L), class = "cier_error_input")
})

test_that("check_items allows NA min on non-reverse items", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L),
                   reverse_keyed = c(FALSE, TRUE, FALSE, FALSE),
                   max = 5L, min = c(NA, 0L, NA, 3L))
  expect_no_error(check_items(it, n_items = 4L))
})

# ---- input-validation hardening ---------------------------------------------

# A non-numeric `responses` names the offending column(s) and adds the ID / free-text fix-it hint.
# A bullet's text lives in the rlang condition message, not the abort header, so collapse
# whitespace over cnd_message().
cnd_text <- function(err) gsub("\\s+", " ", rlang::cnd_message(err))

test_that("a non-numeric data.frame column is named in the responses error", {
  df <- data.frame(respondent_id = c("a", "b"), q1 = c(1, 2), q2 = c(3, 4),
                   stringsAsFactors = FALSE)
  err <- tryCatch(check_responses(df), error = function(e) e)
  expect_s3_class(err, "cier_error_input")
  msg <- cnd_text(err)
  expect_match(msg, "respondent_id")          # the offending column is named
  expect_match(msg, "data columns")           # the ID / free-text fix-it hint
})

test_that("a multi-column non-numeric data.frame names every offender", {
  df <- data.frame(id = c("a", "b"), q1 = c(1, 2),
                   note = c("x", "y"), stringsAsFactors = FALSE)
  err <- tryCatch(check_responses(df), error = function(e) e)
  msg <- cnd_text(err)
  expect_match(msg, "id")
  expect_match(msg, "note")
})

test_that("a non-numeric character matrix still carries the ID hint", {
  # A matrix has no per-column types to name, but the ID / free-text hint must still appear (the
  # most likely cause of a non-numeric matrix).
  err <- tryCatch(check_responses(matrix("a", 2L, 2L)), error = function(e) e)
  expect_s3_class(err, "cier_error_input")
  expect_match(cnd_text(err), "data columns")
})

# A whole-number guard on `check_number(whole = TRUE)` for the count cutoffs.
test_that("check_number(whole = TRUE) rejects a fractional value", {
  expect_error(check_number(2.5, "cutoff", lower = 1, upper = 5, whole = TRUE),
               class = "cier_error_input")
  expect_no_error(check_number(2, "cutoff", lower = 1, upper = 5, whole = TRUE))
  # The bounds still bite under whole = TRUE (below lower / above upper / NA).
  expect_error(check_number(0, "cutoff", lower = 1, upper = 5, whole = TRUE),
               class = "cier_error_input")
  expect_error(check_number(6, "cutoff", lower = 1, upper = 5, whole = TRUE),
               class = "cier_error_input")
  # Default (whole = FALSE) still accepts a fraction.
  expect_no_error(check_number(2.5, "x", lower = 0))
})

# When `items` carries explicit row identifiers (an `item` column, else character rownames) AND
# `responses` has column names, they must align positionally -- a reordered metadata frame is a
# typed error, not a silent reshuffle of reverse-keying / scale assignment.
test_that("check_items aborts when an item column mismatches the response columns", {
  it <- data.frame(item = c("q1", "q2", "q3", "q4"),
                   scale = rep(c("A", "B"), each = 2L), reverse_keyed = FALSE)
  expect_error(
    check_items(it, n_items = 4L, response_names = c("q1", "q3", "q2", "q4")),
    class = "cier_error_input"
  )
  expect_no_error(
    check_items(it, n_items = 4L, response_names = c("q1", "q2", "q3", "q4"))
  )
})

test_that("check_items aborts on a character-rowname mismatch", {
  it <- data.frame(scale = rep(c("A", "B"), each = 2L), reverse_keyed = FALSE)
  rownames(it) <- c("q1", "q2", "q3", "q4")
  expect_error(
    check_items(it, n_items = 4L, response_names = c("q1", "q2", "X", "q4")),
    class = "cier_error_input"
  )
  expect_no_error(
    check_items(it, n_items = 4L, response_names = c("q1", "q2", "q3", "q4"))
  )
})

test_that("the name cross-check names the mismatched position and labels", {
  # Non-digit labels so the only "2" in the message is the reported position (q1/q2 labels would
  # collide with the position number).
  it <- data.frame(item = c("alpha", "beta", "gamma", "delta"),
                   scale = rep(c("A", "B"), each = 2L), reverse_keyed = FALSE)
  err <- tryCatch(
    check_items(it, n_items = 4L,
                response_names = c("alpha", "OTHER", "gamma", "delta")),
    error = function(e) e
  )
  msg <- cnd_text(err)
  expect_match(msg, "[Pp]osition")            # positions are reported
  expect_match(msg, "\\b2\\b")                # the single offender is position 2
  expect_match(msg, "beta")                   # both mismatched labels are surfaced
  expect_match(msg, "OTHER")
})

test_that("the cross-check flags an NA item id as a mismatch (not silently dropped)", {
  # `which(ids != response_names)` would drop an NA-vs-value comparison (NA), so an unlabelled item
  # at a reordered position could slip through; NA-ness is tested explicitly.
  it <- data.frame(item = c("q1", NA, "q3", "q4"),
                   scale = rep(c("A", "B"), each = 2L), reverse_keyed = FALSE)
  expect_error(
    check_items(it, n_items = 4L, response_names = c("q1", "q2", "q3", "q4")),
    class = "cier_error_input"
  )
})

test_that("the cross-check is a no-op without identifiers or column names", {
  # Default item frames (auto integer rownames, no `item` column) never trigger it; neither does a
  # response matrix with no column names.
  auto <- data.frame(scale = rep(c("A", "B"), each = 2L), reverse_keyed = FALSE)
  expect_no_error(check_items(auto, n_items = 4L,
                              response_names = c("a", "b", "c", "d")))
  named <- data.frame(item = c("q1", "q2", "q3", "q4"),
                      scale = rep(c("A", "B"), each = 2L), reverse_keyed = FALSE)
  expect_no_error(check_items(named, n_items = 4L, response_names = NULL))
})

# A one-time message when `items` omits the reverse_keyed column entirely (treating every item as
# forward-keyed); an explicit reverse_keyed = FALSE is silent.
test_that("inform_if_unkeyed messages only when reverse_keyed is absent", {
  absent <- data.frame(scale = rep(c("A", "B"), each = 2L))
  expect_message(inform_if_unkeyed(absent, 4L),
                 class = "cier_message_forward_keyed")
  explicit <- data.frame(scale = rep(c("A", "B"), each = 2L),
                         reverse_keyed = FALSE)
  expect_no_message(inform_if_unkeyed(explicit, 4L))
  # A non-data-frame is left to the real validator -- no message here.
  expect_no_message(inform_if_unkeyed(list(scale = "A"), 1L))
})

test_that("the shared check_items_reverse stays silent (the inform lives in wrappers)", {
  # check_items_reverse defaults the column to all-FALSE for EVERY caller, including cier_simulate's
  # items path -- so the inform must NOT come from here (else it fires on every simulation). It
  # lives in inform_if_unkeyed, called only by the four index wrappers.
  expect_no_message(
    check_items_reverse(data.frame(scale = "A"), 1L, "items", rlang::current_env()),
    class = "cier_message_forward_keyed"
  )
})

# ---- arg-aware hints, all-NA acceptance, the single logical rule -------------

# The validator is shared by every payload argument (responses, checks, page_seconds, reference).
# Three foot-guns closed together: (a) the fix-it hint names the failing argument, not a hardcoded
# `responses[, item_cols]`; (b) an all-NA payload abstains rather than being rejected as "not
# numeric"; (c) a logical column is rejected like a logical matrix, not coerced to 0/1. The single
# logical rule: logical content with an OBSERVED (non-NA) value is non-numeric and rejected; an
# all-NA payload (of any type) is accepted and abstains.

# (b) all-NA acceptance --------------------------------------------------------

test_that("an all-NA logical matrix is accepted and coerced to double", {
  # matrix(NA, n, p) is logical-typed; it was once rejected as "not numeric". It carries no observed
  # responses, so it now abstains.
  out <- check_responses(matrix(NA, nrow = 3L, ncol = 4L))
  expect_true(is.matrix(out) && is.double(out))
  expect_identical(dim(out), c(3L, 4L))
  expect_true(all(is.na(out)))
})

test_that("an all-NA logical matrix matches the all-NA double matrix exactly", {
  # The single all-NA outcome: matrix(NA) is made identical to matrix(NA_real_), which already
  # passed and abstains -- not a separate typed error.
  expect_identical(check_responses(matrix(NA, 3L, 4L)),
                   check_responses(matrix(NA_real_, 3L, 4L)))
})

test_that("an all-NA logical data.frame is accepted and abstains", {
  df <- data.frame(a = c(NA, NA), b = c(NA, NA))   # both columns logical, all-NA
  out <- check_responses(df)
  expect_true(is.double(out) && all(is.na(out)))
})

test_that("all-NA logical and double matrices reach the same index abstention", {
  # Reconciles with the cross-index all-NA invariant: the coerced all-NA matrix flows to the SAME
  # abstention as a double all-NA matrix (no score, flags none) -- an error here would mean
  # matrix(NA) and matrix(NA_real_) diverge.
  a <- suppressWarnings(cier_irv(matrix(NA, nrow = 4L, ncol = 6L)))
  b <- suppressWarnings(cier_irv(matrix(NA_real_, nrow = 4L, ncol = 6L)))
  expect_identical(a$value, b$value)
  expect_identical(a$flagged, b$flagged)
  expect_true(all(is.na(a$value)) && all(is.na(a$flagged)))
})

test_that("an empty all-NA input is still rejected", {
  # all(is.na()) is vacuously TRUE on an empty matrix; the coercion is guarded by a length check
  # (NOT by gate ordering), so it never accepts a 0-row / 0-column payload -- an empty input is
  # still rejected on every axis and storage type.
  expect_error(check_responses(matrix(NA, nrow = 0L, ncol = 3L)),
               class = "cier_error_input")
  expect_error(check_responses(matrix(NA, nrow = 3L, ncol = 0L)),
               class = "cier_error_input")
  expect_error(check_responses(matrix(numeric(0), nrow = 0L, ncol = 3L)),
               class = "cier_error_input")
})

# (c) the single logical rule --------------------------------------------------

test_that("a logical data.frame column is rejected and named", {
  # This once coerced silently to active -> 0/1 and passed; a logical *matrix* of the same content
  # was rejected. Close the asymmetry by rejecting.
  df <- data.frame(active = c(TRUE, FALSE), q1 = c(1, 2), q2 = c(3, 4))
  err <- tryCatch(check_responses(df), error = function(e) e)
  expect_s3_class(err, "cier_error_input")
  expect_match(cnd_text(err), "active")            # the logical column is named
  expect_identical(err$data$observed, "active")    # ... and recorded in the payload
})

test_that("an all-logical data.frame with observed values is rejected", {
  df <- data.frame(a = c(TRUE, FALSE), b = c(FALSE, TRUE))
  err <- tryCatch(check_responses(df), error = function(e) e)
  expect_s3_class(err, "cier_error_input")
  expect_setequal(err$data$observed, c("a", "b"))  # both logical columns named
})

test_that("an all-NA logical column beside a numeric column is NOT rejected", {
  # The crux of the single rule: a logical column with no observed value carries no evidence, so it
  # abstains -- it must NOT be named or rejected. This separates the fix from a "coerce / reject any
  # logical" variant.
  df <- data.frame(a = c(NA, NA), b = c(1, 2))      # a is logical, all-NA
  expect_no_error(check_responses(df))
})

test_that("a mixed frame names ONLY the observed logical column (per-column)", {
  # The observed filter is per-column, not per-frame: `a` is observed (named), `b` is an all-NA
  # logical (abstains, NOT named) -- pins against a filter that names every logical column once any
  # logical column is observed.
  df <- data.frame(a = c(TRUE, FALSE), b = c(NA, NA), c = c(1, 2))
  err <- tryCatch(check_responses(df), error = function(e) e)
  expect_s3_class(err, "cier_error_input")
  expect_identical(err$data$observed, "a")
})

test_that("a logical column with a single observed value is rejected (any-not-all)", {
  # The predicate is `any(!is.na())`, not `all(!is.na())`: one observed TRUE is enough to make the
  # column non-numeric. An all-observed predicate would wrongly abstain on c(TRUE, NA).
  df <- data.frame(a = c(TRUE, NA), b = c(1, 2))
  err <- tryCatch(check_responses(df), error = function(e) e)
  expect_s3_class(err, "cier_error_input")
  expect_identical(err$data$observed, "a")
})

test_that("a blank- or NA-named logical column is still rejected (decision != display)", {
  # Rejection turns on whether an offender EXISTS, not on whether it has a printable name. A
  # header-less frame (e.g. read.csv(check.names = FALSE)) whose logical column has an empty name
  # once tripped silently: the name filter that trims the printed list also gated the rejection, so
  # the column coerced to 0/1 and slipped through. Reject it; the printed name is just absent.
  df <- data.frame(active = c(TRUE, FALSE), q1 = c(1, 2))
  names(df) <- c("", "")
  expect_error(check_responses(df), class = "cier_error_input")
  na_named <- data.frame(active = c(TRUE, FALSE), q1 = c(1, 2))
  names(na_named) <- c(NA, "q1")
  expect_error(check_responses(na_named), class = "cier_error_input")
})

test_that("an all-NA non-numeric column beside numeric data is rejected AND named", {
  # An all-NA *character* (or factor) column has no observed value but still poisons as.matrix() to
  # character, so the frame is non-numeric and rejected -- and the offending column must be NAMED
  # (the observed filter applies only when the matrix coerced cleanly to numeric, i.e. logical-0/1).
  df <- data.frame(comment = c(NA_character_, NA_character_), q1 = c(1, 2))
  err <- tryCatch(check_responses(df), error = function(e) e)
  expect_s3_class(err, "cier_error_input")
  expect_identical(err$data$observed, "comment")
})

test_that("a zero-row frame names its non-numeric column, not just its dimensions", {
  # The non-numeric gate runs before the dimension check, so a malformed empty frame with a stray
  # non-numeric column points at the column (the actionable cause) rather than only the empty dims.
  df <- data.frame(id = character(0), q1 = numeric(0))
  err <- tryCatch(check_responses(df), error = function(e) e)
  expect_s3_class(err, "cier_error_input")
  expect_identical(err$data$observed, "id")
})

test_that("a logical matrix with observed values is still rejected (locked)", {
  # The asymmetry is closed by rejecting the logical column, NOT by accepting the logical matrix: a
  # "coerce any logical -> double" rule would wrongly pass this.
  expect_error(check_responses(matrix(c(TRUE, FALSE, TRUE, FALSE), nrow = 2L)),
               class = "cier_error_input")
})

# (a) arg-aware, softened fix-it hint ------------------------------------------

test_that("the non-numeric hint names the failing argument, not responses", {
  # check_page_seconds / cier_attention / the reference= path all route through check_responses with
  # a non-default arg; the hint must reference {arg}.
  df <- data.frame(label = c("x", "y"), q1 = c(1, 2), stringsAsFactors = FALSE)
  msg <- cnd_text(tryCatch(check_responses(df, arg = "page_seconds"),
                           error = function(e) e))
  expect_match(msg, "page_seconds")                # the actual arg is referenced
  expect_match(msg, "data columns")                # softened from "item columns"
  expect_false(grepl("item columns", msg))         # the responses-only noun is gone
  expect_false(grepl("responses\\[", msg))         # no hardcoded responses[...]
  err <- tryCatch(check_responses(df, arg = "page_seconds"), error = function(e) e)
  expect_identical(err$data$arg, "page_seconds")   # the payload arg is the real one too
})

test_that("the default-arg non-numeric hint still references responses", {
  df <- data.frame(label = c("x", "y"), q1 = c(1, 2), stringsAsFactors = FALSE)
  msg <- cnd_text(tryCatch(check_responses(df), error = function(e) e))
  expect_match(msg, "responses")
  expect_match(msg, "data columns")
})
