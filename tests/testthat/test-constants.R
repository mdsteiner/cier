# Controlled-vocabulary accessors (R/constants.R), trimmed to what the
# registry loader/validator pin against.

test_that("family vocabulary is the v0.2 set", {
  expect_identical(cier_family_levels(),
                   c("indirect", "personfit", "timing", "direct"))
})

test_that("cutoff-method vocabulary is the v0 set", {
  expect_identical(cier_cutoff_methods(),
                   c("percentile", "fixed", "chisq", "perfit_null"))
})

test_that("flag-direction vocabulary is upper/lower", {
  expect_identical(cier_flag_directions(), c("upper", "lower"))
})

test_that("accessors are pure (return a fresh copy each call)", {
  a <- cier_cutoff_methods()
  a[1] <- "tampered"
  expect_identical(cier_cutoff_methods()[1], "percentile")
})
