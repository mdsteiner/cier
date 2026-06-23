# cli / glue interpolation safety.
#
# User-supplied names (response columns, item ids, pass labels, screen method ids/control keys)
# flow into cli error messages and may carry literal braces or glue-like tokens -- "score_{q1}",
# "x{1+1}", even "z{options(...)}". cier interpolates them through the nested {.val {var}} /
# {.code {var}} form, so cli treats the value as LITERAL DATA, never re-parsed as a glue template.
# These guards pin that contract end to end: a clean typed cier_error_input, the token verbatim,
# and decisively no execution of an injected expression. A regression to direct interpolation
# (paste()-ing a user name into a template) would surface as a glue eval error or executed side
# effect instead of the expected validation error.

cnd_text <- function(e) gsub("[[:space:]]+", " ", rlang::cnd_message(e))

# Force `expr` lazily inside tryCatch: expect a typed input error whose message reproduces `token`
# verbatim -- proof the braces were data, not an evaluated template.
expect_literal_token <- function(expr, token, info = token) {
  err <- tryCatch(suppressMessages(suppressWarnings(expr)), error = function(e) e)
  expect_s3_class(err, "cier_error_input")
  seen <- if (inherits(err, "condition")) cnd_text(err) else "(no error raised)"
  expect_true(grepl(token, seen, fixed = TRUE), info = paste0(info, " :: ", seen))
}

# ---- payload validators (check_responses / items / pass) --------------------

test_that("a non-numeric column name with glue braces is rendered literally", {
  df <- data.frame(a = 1:5, b = 1:5, z = letters[1:5], stringsAsFactors = FALSE)
  names(df)[3] <- "score_{q1}"
  expect_literal_token(check_responses(df), "score_{q1}")
  # The literal name is also preserved in the structured payload, unmangled.
  err <- tryCatch(check_responses(df), error = function(e) e)
  expect_identical(err$data$observed, "score_{q1}")
})

test_that("a brace name that would be a valid glue expression is not evaluated", {
  # `{1+1}` is a syntactically valid glue expression; if cier re-parsed the value
  # it would interpolate to "2" (or error). It must appear verbatim instead.
  df <- data.frame(q1 = 1:5, bad = letters[1:5], stringsAsFactors = FALSE)
  names(df)[2] <- "x{1+1}"
  expect_literal_token(check_responses(df), "x{1+1}")
})

test_that("an item-id glue token that would run code is inert", {
  # The decisive injection guard, by OBSERVABLE side effect: an item id whose braces, if evaluated
  # as a glue template, would call options() and leave a trace. The token must surface verbatim in
  # the ordinary alignment error, and the side effect must NOT have run.
  withr::defer(options(cier.brace.injected = NULL))
  options(cier.brace.injected = NULL)
  token <- "z{options(cier.brace.injected='breached')}"
  r <- matrix(seq_len(30L), nrow = 5L, ncol = 6L)
  colnames(r) <- paste0("c", 1:6)
  it <- data.frame(scale = rep(c("A", "B"), each = 3L), reverse_keyed = FALSE,
                   max = 5L,
                   item = c("c1", "c2", "c3", "c4", "c5", token),
                   stringsAsFactors = FALSE)
  err <- tryCatch(suppressMessages(suppressWarnings(cier_even_odd(r, it))),
                  error = function(e) e)
  expect_s3_class(err, "cier_error_input")        # a clean validation error ...
  expect_true(grepl(token, cnd_text(err), fixed = TRUE))   # ... id verbatim ...
  # ... and decisively: the injected options() call never executed.
  expect_null(getOption("cier.brace.injected"))
})

test_that("an attention pass-set label with braces is rendered literally", {
  chk <- matrix(c(1, 0, 3, 0, 2, 5), ncol = 2L, byrow = TRUE)
  colnames(chk) <- c("q1", "q2")
  # A named pass that overlaps the columns triggers the positional cross-check;
  # the mismatching label carries a brace token.
  expect_literal_token(
    cier_attention(chk, pass = stats::setNames(list(c(1, 2), 0), c("q1", "q2{x}"))),
    "q2{x}"
  )
})

# ---- cier_screen control-surface errors -------------------------------------

test_that("cier_screen method and control errors render glue-token names literally", {
  r <- matrix(sample(1:5, 30L, replace = TRUE), nrow = 5L, ncol = 6L)
  # Unknown method id with a glue token.
  expect_literal_token(cier_screen(r, methods = c("cier_irv", "boom{1+1}")),
                       "boom{1+1}")
  # Control keyed by a not-selected method whose name carries braces.
  expect_literal_token(
    cier_screen(r, control = stats::setNames(list(list(fpr = 0.1)), "evil{q}")),
    "evil{q}"
  )
  # Control argument name with braces for a real, selected method.
  expect_literal_token(
    cier_screen(r, methods = "cier_irv",
                control = list(cier_irv = stats::setNames(list(1), "zz{1}"))),
    "zz{1}"
  )
})
