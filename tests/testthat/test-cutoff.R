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

test_that("every percentile registry default stores the fpr (0.05), direction-agnostic", {
  reg <- cier_methods()
  pct <- reg[reg$default_cutoff_method == "percentile", , drop = FALSE]
  expect_gt(nrow(pct), 0L)
  # The registry stores the tail mass (fpr), NOT the post-flip quantile, so a
  # lower-tail index and an upper-tail index both store 0.05; resolve_cutoff
  # applies the single direction flip.
  expect_equal(pct$default_cutoff_value, rep(0.05, nrow(pct)))

  # Fed that one fpr with the row's direction, the cutoff lands on the correct
  # tail: lower -> q05, upper -> q95. (Guards the inversion an upper index would
  # suffer if its registry value were passed through the flip as if it were the
  # post-flip 0.95 quantile.)
  x <- as.numeric(1:100)
  q05 <- as.numeric(stats::quantile(x, 0.05, names = FALSE, type = 7L))
  q95 <- as.numeric(stats::quantile(x, 0.95, names = FALSE, type = 7L))
  for (i in seq_len(nrow(pct))) {
    got <- resolve_cutoff(x, pct$flag_direction[i],
                          fpr = pct$default_cutoff_value[i])
    want <- if (identical(pct$flag_direction[i], "upper")) q95 else q05
    expect_identical(got, want, info = pct$method[i])
  }
})

# ---- resolve_cutoff: chi-square ---------------------------------------------

test_that("chisq cutoff equals qchisq(1 - alpha, df)", {
  expect_identical(resolve_cutoff(method = "chisq", df = 44),
                   as.numeric(stats::qchisq(1 - 0.001, df = 44)))
  expect_identical(resolve_cutoff(method = "chisq", df = 10, alpha = 0.01),
                   as.numeric(stats::qchisq(1 - 0.01, df = 10)))
})

# Note: resolve_cutoff() is an internal resolver that trusts its inputs (df,
# alpha, method, direction, fpr). The public wrappers validate the user-supplied
# rate / cutoff before calling -- e.g. cier_irv() rejects fpr outside (0, 1) and
# cier_mahalanobis() rejects alpha outside (0, 1) -- so those input-error cases
# are tested at the wrapper level, not here.

# ---- resolve_cutoff: fixed ---------------------------------------------------

test_that("fixed cutoff returns the supplied value verbatim (no n_items)", {
  # A literal threshold (e.g. an absolute score cutoff) passes through unchanged.
  expect_identical(resolve_cutoff(method = "fixed", value = 22), 22)
  expect_identical(resolve_cutoff(method = "fixed", value = 3.5), 3.5)
})

test_that("fixed cutoff scales by the item count when n_items is given (fraction)", {
  # The fraction path (longstring's `frac`): ceiling(value * n_items).
  expect_identical(resolve_cutoff(method = "fixed", value = 0.5, n_items = 10), 5)
  expect_identical(resolve_cutoff(method = "fixed", value = 0.25, n_items = 8), 2)
  expect_identical(resolve_cutoff(method = "fixed", value = 1, n_items = 8), 8)
})

test_that("the fraction path is robust to floating-point error (no ceiling bump)", {
  # 0.28 * 25 == 7.0000000000000009 in IEEE-754; ceiling must still give 7, not
  # 8. Likewise 0.56 * 25 == 14.000...; both are silent off-by-one footguns if
  # the product is not rounded before ceiling().
  expect_identical(resolve_cutoff(method = "fixed", value = 0.28, n_items = 25), 7)
  expect_identical(resolve_cutoff(method = "fixed", value = 0.56, n_items = 25), 14)
  # A genuine fraction still rounds up: 0.21 * 25 = 5.25 -> 6.
  expect_identical(resolve_cutoff(method = "fixed", value = 0.21, n_items = 25), 6)
})

test_that("an unknown cutoff method is a typed state error, not a silent percentile", {
  # resolve_cutoff() handles only the value-only methods (percentile / fixed /
  # chisq). A model-referenced method like "perfit_null" (Gnormed) is resolved at
  # its bridge from the fitted object and must NOT silently fall through to the
  # percentile branch -- that would mis-resolve a future registry-driven caller.
  expect_error(resolve_cutoff(values = 1:10, method = "perfit_null"),
               class = "cier_error_state")
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

# ---- kneedle: the paper-faithful Laz.R cutoff (Satopaa et al. 2011) ----------
# The independent oracle (ref_kneedle) re-derives the parameter-free elbow and
# never calls the production kernel; oracle-only trust, tolerance 0 (kernel and
# oracle do byte-identical float ops). cier ships only the convex + increasing
# case (high Laz.R = careless = upper tail). See tests/reference/TOLERANCES.md.

source(test_path("..", "reference", "ref-kneedle-satopaa-2011.R"))

test_that("kneedle_knee equals a hand-computed elbow (independent of the oracle)", {
  # Sorted v = c(0, 0, 0, 0, 1, 2, 5, 9), n = 8. x = (0:7)/7; y = v / 9. The
  # deviations y - x are 0, -.143, -.286, -.429, -.460, -.492, -.302, 0, so the
  # most-negative (the convex knee) is index 6 -> value 2. Computed by hand, so a
  # mutant that drifts the normalisation, the argmin, or the sort fails here even
  # if it still happens to match the oracle.
  expect_identical(kneedle_knee(c(9, 0, 5, 0, 2, 0, 1, 0)), 2)   # unsorted input
})

test_that("kneedle_knee matches the oracle byte-for-byte on a right-skewed sample", {
  vals <- withr::with_seed(20260611L, sort(rexp(80L, rate = 4)))
  expect_identical(kneedle_knee(vals),
                   ref_kneedle(vals, "convex", "increasing")$value)
})

test_that("kneedle_knee matches the oracle on a smooth convex curve (x^2)", {
  vals <- (seq(0, 1, length.out = 500L))^2
  expect_identical(kneedle_knee(vals),
                   ref_kneedle(vals, "convex", "increasing")$value)
  # Analytic elbow of y = x^2 on [0, 1] sits at x = 0.5 -> value 0.25.
  expect_lt(abs(kneedle_knee(vals) - 0.25), 0.01)
})

test_that("kneedle_knee lands at the bulk/tail boundary, not in the bulk or tail", {
  # A gentle moderate bulk (40 values 0.30..0.34) then a sharp careless spike
  # (10 values 0.90..1.00). The convex elbow sits at the top of the bulk / start
  # of the gap -- at or above every bulk value and strictly below the spike -- so
  # flagging at/above it isolates the careless cluster, never the whole bulk and
  # never nobody. The knee equals max(bulk) exactly here (verified), so pin
  # `>=` not `>`; an argmax / inverted-elbow mutant lands at min(bulk) and fails.
  bulk <- seq(0.30, 0.34, length.out = 40L)
  tail <- seq(0.90, 1.00, length.out = 10L)
  knee <- kneedle_knee(c(bulk, tail))
  expect_gte(knee, max(bulk))                            # at/above the bulk
  expect_lt(knee, min(tail))                             # strictly below the spike
})

test_that("kneedle_knee pins the (i-1)/(n-1) normalisation (off-by-one guard)", {
  # A convex increasing curve with an INTERIOR elbow at index 4 -> value 4.591
  # (verified in R). The (i-1)/n off-by-one moves the argmin to the left endpoint
  # (2.287) and the concave/argmax mutant gives 3.484, so this fixture kills the
  # normalisation off-by-one that every smooth / right-skewed parity fixture
  # happens to be degenerate against. (The i/(n-1) variant is a uniform shift of
  # the x-axis: it can never change the argmin, so it is an equivalent mutant
  # that needs no fixture.)
  v <- c(2.287, 3.484, 4.178, 4.591, 5.561, 6.509, 7.257)
  expect_identical(kneedle_knee(v), ref_kneedle(v, "convex", "increasing")$value)
  expect_identical(kneedle_knee(v), 4.591)
})

test_that("resolve_kneedle_cutoff drops non-finite values, then matches the kernel", {
  value <- c(NA_real_, 0.1, 0.2, NaN, 0.3, 0.4, Inf, 0.5, 0.95)
  finite <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.95)
  expect_identical(resolve_kneedle_cutoff(value), kneedle_knee(finite))
})

test_that("resolve_kneedle_cutoff abstains (NA + warning) below three finite values", {
  expect_warning(res <- resolve_kneedle_cutoff(c(0.4, 0.9, NA, Inf)),
                 class = "cier_warning_insufficient_items")
  expect_identical(res, NA_real_)
  expect_warning(resolve_kneedle_cutoff(rep(NA_real_, 5L)),
                 class = "cier_warning_insufficient_items")
})

test_that("resolve_kneedle_cutoff abstains on a constant distribution (no knee)", {
  # The signed-off deviation from the archive (which returned the constant value,
  # flagging everyone): a knee is undefined on a degenerate distribution, so the
  # cutoff abstains and flags nobody.
  expect_warning(res <- resolve_kneedle_cutoff(rep(0.42, 50L)),
                 class = "cier_warning_insufficient_items")
  expect_identical(res, NA_real_)
})

test_that("resolve_kneedle_cutoff is deterministic (no RNG, no jitter)", {
  vals <- withr::with_seed(11L, sort(rexp(60L)))
  expect_identical(resolve_kneedle_cutoff(vals), resolve_kneedle_cutoff(vals))
})

test_that("the kneedle elbow differs from the fpr=0.05 percentile (a distinct rule)", {
  # A 40-point bulk at 0.30 then a graded 20-point ramp 0.50..1.00. The elbow
  # sits at the bulk/ramp boundary (~0.3-0.5); the 95th percentile sits high in
  # the ramp (~0.92). They are different cutoff rules, so a mutant that aliases
  # kneedle to the percentile cannot survive: the elbow is strictly the lower.
  v <- c(rep(0.30, 40L), seq(0.50, 1.00, length.out = 20L))
  expect_lt(resolve_kneedle_cutoff(v),
            resolve_cutoff(v, "upper", fpr = 0.05))
})
