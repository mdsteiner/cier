# Controlled-vocabulary accessors in R/constants.R.

test_that("cutoff-method vocabulary is the v0 set", {
  expect_identical(cier_cutoff_methods(),
                   c("percentile", "fixed", "chisq", "mc_null"))
})

test_that("flag-direction vocabulary is upper/lower", {
  expect_identical(cier_flag_directions(), c("upper", "lower"))
})

test_that("accessors are pure (return a fresh copy each call)", {
  a <- cier_cutoff_methods()
  a[1] <- "tampered"
  expect_identical(cier_cutoff_methods()[1], "percentile")
})
