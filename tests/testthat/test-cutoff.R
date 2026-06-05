# Cutoff resolution (resolve_cutoff) and flag application (apply_flag).
# The reference for the percentile / chi-square branches is base R itself
# (stats::quantile / stats::qchisq), asserted inline. See
# tests/reference/TOLERANCES.md.

# ---- resolve_cutoff: percentile + the NO-FLIP guard --------------------------

test_that("percentile cutoff applies exactly one direction flip (NO-FLIP)", {
  x <- as.numeric(1:100)
  q_lo <- as.numeric(stats::quantile(x, 0.05, names = FALSE, type = 7L))
  q_hi <- as.numeric(stats::quantile(x, 0.95, names = FALSE, type = 7L))

  # lower flags the low tail -> probs = fpr; upper flags the high tail ->
  # probs = 1 - fpr (the 0.95 analogue). The registry stores these literal
  # directional quantiles, so the kernel must flip ONCE here, never twice.
  expect_identical(resolve_cutoff(x, "lower", fpr = 0.05), q_lo)
  expect_identical(resolve_cutoff(x, "upper", fpr = 0.05), q_hi)

  # Exactly one flip: upper != lower at the same fpr. A missing-flip mutant
  # returns q_lo for both; a double-flip mutant returns q_lo for upper. Either
  # collapses this inequality.
  expect_false(identical(resolve_cutoff(x, "upper", fpr = 0.05),
                         resolve_cutoff(x, "lower", fpr = 0.05)))

  # The flip lives only in the value; a non-default fpr still maps cleanly.
  expect_identical(resolve_cutoff(x, "upper", fpr = 0.10),
                   as.numeric(stats::quantile(x, 0.90, names = FALSE, type = 7L)))
  expect_identical(resolve_cutoff(x, "lower", fpr = 0.01),
                   as.numeric(stats::quantile(x, 0.01, names = FALSE, type = 7L)))
})

test_that("the registry percentile defaults equal the single-flip output at fpr = 0.05", {
  reg <- cier_methods()
  pct <- reg[reg$default_cutoff_method == "percentile", , drop = FALSE]
  expect_gt(nrow(pct), 0L)
  for (i in seq_len(nrow(pct))) {
    expected_p <- if (identical(pct$flag_direction[i], "upper")) 1 - 0.05 else 0.05
    expect_equal(pct$default_cutoff_value[i], expected_p, info = pct$method[i])
  }
})

# ---- resolve_cutoff: chi-square ---------------------------------------------

test_that("chisq cutoff equals qchisq(1 - alpha, df)", {
  expect_identical(resolve_cutoff(method = "chisq", df = 44),
                   as.numeric(stats::qchisq(1 - 0.001, df = 44)))
  expect_identical(resolve_cutoff(method = "chisq", df = 10, alpha = 0.01),
                   as.numeric(stats::qchisq(1 - 0.01, df = 10)))
})

test_that("chisq without df is a typed input error", {
  expect_error(resolve_cutoff(method = "chisq"), class = "cier_error_input")
})

test_that("chisq with a non-positive or non-numeric df is a typed input error", {
  expect_error(resolve_cutoff(method = "chisq", df = 0), class = "cier_error_input")
  expect_error(resolve_cutoff(method = "chisq", df = "x"), class = "cier_error_input")
})

# ---- resolve_cutoff: fixed ---------------------------------------------------

test_that("fixed cutoff returns the supplied value verbatim", {
  # The longstring wrapper resolves ceiling(0.5 * ncol) and passes the integer;
  # resolve_cutoff just routes it through the single cutoff path.
  expect_identical(resolve_cutoff(method = "fixed", value = 22), 22)
  expect_identical(resolve_cutoff(method = "fixed", value = 3.5), 3.5)
})

test_that("fixed without value is a typed input error", {
  expect_error(resolve_cutoff(method = "fixed"), class = "cier_error_input")
})

# ---- resolve_cutoff: abstention + edges -------------------------------------

test_that("percentile abstains (NA + typed warning) when no finite values", {
  expect_warning(
    res <- resolve_cutoff(rep(NA_real_, 5L), "lower", fpr = 0.05),
    class = "cier_warning_insufficient_items"
  )
  expect_identical(res, NA_real_)
  expect_warning(resolve_cutoff(numeric(0), "upper"),
                 class = "cier_warning_insufficient_items")
})

test_that("percentile drops NaN/Inf before the quantile", {
  x <- c(1, 2, 3, 4, Inf, NaN, NA)
  finite <- c(1, 2, 3, 4)
  expect_identical(
    resolve_cutoff(x, "lower", fpr = 0.25),
    as.numeric(stats::quantile(finite, 0.25, names = FALSE, type = 7L))
  )
})

test_that("percentile handles a single finite value and constant input", {
  expect_identical(resolve_cutoff(c(5, NA), "lower", fpr = 0.05), 5)
  expect_identical(resolve_cutoff(rep(7, 10L), "upper", fpr = 0.05), 7)
})

test_that("invalid method or direction is a typed input error", {
  expect_error(resolve_cutoff(1:10, method = "kneedle"), class = "cier_error_input")
  expect_error(resolve_cutoff(1:10, direction = "sideways"), class = "cier_error_input")
})

test_that("fpr outside the open interval (0, 1) is a typed input error", {
  # A target false-positive rate of 0% or 100% is meaningless, as are values
  # outside [0, 1]; reject the whole closed boundary and beyond.
  expect_error(resolve_cutoff(1:10, fpr = 1.5), class = "cier_error_input")
  expect_error(resolve_cutoff(1:10, fpr = 1), class = "cier_error_input")
  expect_error(resolve_cutoff(1:10, fpr = 0), class = "cier_error_input")
  expect_error(resolve_cutoff(1:10, fpr = -0.1), class = "cier_error_input")
})

# ---- apply_flag --------------------------------------------------------------

test_that("apply_flag flags the correct tail and includes ties at the cutoff", {
  v <- c(1, 5, 10, NA)
  expect_identical(apply_flag(v, cutoff = 5, direction = "upper"),
                   c(FALSE, TRUE, TRUE, FALSE))
  expect_identical(apply_flag(v, cutoff = 5, direction = "lower"),
                   c(TRUE, TRUE, FALSE, FALSE))
})

test_that("apply_flag never flags NA values and abstains on an NA cutoff", {
  v <- c(1, 5, NA, 10)
  expect_identical(apply_flag(v, cutoff = NA_real_, direction = "upper"),
                   rep(FALSE, 4L))
  expect_false(apply_flag(v, cutoff = 5, direction = "upper")[[3L]])
})

test_that("apply_flag rejects an invalid direction (silent tail inversion guard)", {
  expect_error(apply_flag(c(1, 2), cutoff = 1, direction = "sideways"),
               class = "cier_error_input")
})
