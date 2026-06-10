# Typed condition system (R/conditions.R). Conditions are asserted by
# class, never by message text.

test_that("cier_abort raises an error of the requested sub-class", {
  err <- expect_error(
    cier_abort("cier_error_input", "bad input"),
    class = "cier_error_input"
  )
  expect_s3_class(err, "cier_error_input")
  expect_s3_class(err, "cier_error")
  expect_s3_class(err, "error")
  expect_s3_class(err, "condition")
})

test_that("each v0 error sub-class is class-stable", {
  subs <- c("cier_error_input", "cier_error_data", "cier_error_state")
  for (s in subs) {
    cond <- tryCatch(cier_abort(s, "x"), error = function(e) e)
    expect_s3_class(cond, s)
    expect_s3_class(cond, "cier_error")
    expect_s3_class(cond, "error")
  }
})

test_that("cier_warn raises a typed warning catchable via the parent class", {
  caught <- NULL
  withCallingHandlers(
    cier_warn("cier_warning_insufficient_items", "soft issue"),
    cier_warning = function(w) {
      caught <<- w
      invokeRestart("muffleWarning")
    }
  )
  expect_s3_class(caught, "cier_warning_insufficient_items")
  expect_s3_class(caught, "cier_warning")
  expect_s3_class(caught, "warning")
})

test_that("cier_inform raises a typed message catchable via the parent", {
  caught <- NULL
  withCallingHandlers(
    cier_inform("cier_message_progress", "noted"),
    cier_message = function(m) {
      caught <<- m
      invokeRestart("muffleMessage")
    }
  )
  expect_s3_class(caught, "cier_message_progress")
  expect_s3_class(caught, "cier_message")
  expect_s3_class(caught, "message")
})

test_that("the minimal class chain is exactly what is documented", {
  cond_err <- tryCatch(cier_abort("cier_error_input", "x"), error = function(e) e)
  expect_identical(
    class(cond_err)[1:4],
    c("cier_error_input", "cier_error", "error", "condition")
  )
  cond_warn <- NULL
  withCallingHandlers(
    cier_warn("cier_warning_insufficient_items", "x"),
    cier_warning = function(w) {
      cond_warn <<- w
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(
    class(cond_warn)[1:4],
    c("cier_warning_insufficient_items", "cier_warning", "warning", "condition")
  )
})

test_that("data payload is preserved on the condition object", {
  payload <- list(arg = "scale", respondent = 42L, value = 7)
  cond <- tryCatch(
    cier_abort("cier_error_input", "x", data = payload),
    error = function(e) e
  )
  expect_identical(cond$data, payload)
  extracted <- cier_condition_data(cond)
  expect_identical(extracted, payload)
  expect_identical(extracted$respondent, 42L)
})

# The `class` and `data` argument guards inside cier_abort/warn/inform are
# deliberately base-R `stop()`s (a typed cier condition cannot validate the
# very inputs used to build a typed cier condition), so these assert a plain
# error, not a `cier_error_*` class.
test_that("data must be a (possibly empty) named list", {
  expect_error(cier_abort("cier_error_input", "x", data = 1:3))
  expect_error(cier_abort("cier_error_input", "x", data = list(1, 2)))
  expect_error(cier_warn("cier_warning_insufficient_items", "x", data = 1:3))
  expect_error(cier_inform("cier_message_progress", "x", data = list(1, 2)))
})

test_that("an empty named list payload is accepted", {
  expect_error(cier_abort("cier_error_input", "x", data = list()),
               class = "cier_error_input")
})

test_that("class argument must be a non-empty character vector with no NAs", {
  expect_error(cier_abort(character(0), "x"))
  expect_error(cier_abort(NULL, "x"))
  expect_error(cier_warn(123L, "x"))
  expect_error(cier_inform(NA_character_, "x"))
})

test_that("cli inline markup interpolates in the message body", {
  field <- "scale"
  err <- tryCatch(
    cier_abort("cier_error_input", "Bad {.field {field}}"),
    error = function(e) e
  )
  expect_s3_class(err, "cier_error_input")
  expect_true(grepl("scale", conditionMessage(err), fixed = TRUE))
})

test_that("cier_condition_data returns empty list when payload absent", {
  raw <- structure(
    list(message = "x", call = NULL),
    class = c("cier_error_input", "cier_error", "error", "condition")
  )
  expect_identical(cier_condition_data(raw), list())
})

test_that("cier_warn and cier_inform return invisibly", {
  expect_null(suppressWarnings(cier_warn("cier_warning_insufficient_items", "x")))
  expect_null(suppressMessages(cier_inform("cier_message_progress", "x")))
})
