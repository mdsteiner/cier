# Cutoff resolution (resolve_cutoff) and flag application (apply_flag). The percentile / chi-square
# reference is base R itself (stats::quantile / stats::qchisq), asserted inline.

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

test_that("every percentile method default stores the fpr (0.05), direction-agnostic", {
  specs <- cier_method_specs()
  pct <- specs[specs$default_cutoff_method == "percentile", , drop = FALSE]
  expect_gt(nrow(pct), 0L)
  # The method spec stores the tail mass (fpr), NOT the post-flip quantile, so a
  # lower-tail index and an upper-tail index both store 0.05; resolve_cutoff
  # applies the single direction flip.
  expect_equal(pct$default_cutoff_value, rep(0.05, nrow(pct)))

  # Fed that one fpr with the row's direction, the cutoff lands on the correct
  # tail: lower -> q05, upper -> q95. (Guards the inversion an upper index would
  # suffer if its method default were passed through the flip as if it were the
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

# resolve_cutoff() is an internal resolver that trusts its inputs (df, alpha, method, direction,
# fpr). The public wrappers validate the user rate / cutoff first -- cier_irv() rejects fpr outside
# (0, 1), cier_mahalanobis() rejects alpha outside (0, 1) -- so input-error cases are tested there.

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
  # 0.28 * 25 == 7.0000000000000009 in IEEE-754; ceiling must still give 7, not 8. Likewise
  # 0.56 * 25 == 14.000...; both are off-by-one footguns if the product is not rounded first.
  expect_identical(resolve_cutoff(method = "fixed", value = 0.28, n_items = 25), 7)
  expect_identical(resolve_cutoff(method = "fixed", value = 0.56, n_items = 25), 14)
  # A genuine fraction still rounds up: 0.21 * 25 = 5.25 -> 6.
  expect_identical(resolve_cutoff(method = "fixed", value = 0.21, n_items = 25), 6)
})

test_that("an unknown cutoff method is a typed state error, not a silent percentile", {
  # resolve_cutoff() handles only the value-only methods (percentile / fixed / chisq). A
  # null-referenced method like "mc_null" (Gnormed) is resolved at its bridge from the
  # scored response block and must NOT fall through to percentile -- that mis-resolves
  # a registry-driven caller.
  expect_error(resolve_cutoff(values = 1:10, method = "mc_null"),
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

# ---- resolve_cutoff: degeneracy guards (constant, small-n, saturation) -------
# A percentile cutoff is meaningless on a degenerate distribution. Three guards in the order
# small-n -> constant -> saturation (mutually exclusive). The kneedle resolver above abstains on
# the same degeneracy; these align the percentile path. Expectations are re-derived by hand (the
# threshold ceiling(1/fpr) and quantile order statistics), never by calling the resolver.

test_that("percentile abstains below ceiling(1/fpr) finite scores", {
  # The 5%-tail cutoff is undefined until at least ceiling(1/0.05) = 20 scores exist; 19 abstains,
  # 20 resolves. (A single finite value flagging a perfect score -- the old behaviour -- is wrong.)
  expect_warning(res <- resolve_cutoff(as.numeric(1:19), "upper", fpr = 0.05),
                 class = "cier_warning_insufficient_items")
  expect_identical(res, NA_real_)
  # n = 20 clears the threshold and resolves to the hand-computed 95th percentile.
  expect_identical(resolve_cutoff(as.numeric(1:20), "upper", fpr = 0.05),
                   as.numeric(stats::quantile(1:20, 0.95, names = FALSE,
                                              type = 7L)))      # = 19.05
  # ...silently: the interior quantile must not trip the saturation warning.
  expect_no_warning(resolve_cutoff(as.numeric(1:20), "upper", fpr = 0.05))
  # The threshold scales with fpr: 1/0.10 -> 10, 1/0.01 -> 100. Both directions.
  expect_warning(resolve_cutoff(as.numeric(1:9), "upper", fpr = 0.10),
                 class = "cier_warning_insufficient_items")
  expect_identical(resolve_cutoff(as.numeric(1:10), "lower", fpr = 0.10),
                   as.numeric(stats::quantile(1:10, 0.10, names = FALSE,
                                              type = 7L)))
  expect_warning(resolve_cutoff(as.numeric(1:99), "lower", fpr = 0.01),
                 class = "cier_warning_insufficient_items")
  expect_identical(resolve_cutoff(as.numeric(1:100), "lower", fpr = 0.01),
                   as.numeric(stats::quantile(1:100, 0.01, names = FALSE,
                                              type = 7L)))
})

test_that("the threshold is ceiling(1/fpr), not floor (non-integer reciprocal)", {
  # fpr = 0.03 -> 1/fpr = 33.33..; ceiling = 34. n = 33 must abstain (a floor=33 mutant resolves
  # it); n = 34 resolves. The {0.05, 0.10, 0.01} fixtures all have an integer reciprocal, so this
  # is the only test that separates the two.
  expect_warning(res <- resolve_cutoff(as.numeric(1:33), "upper", fpr = 0.03),
                 class = "cier_warning_insufficient_items")
  expect_identical(res, NA_real_)
  expect_identical(resolve_cutoff(as.numeric(1:34), "upper", fpr = 0.03),
                   as.numeric(stats::quantile(1:34, 0.97, names = FALSE,
                                              type = 7L)))
})

test_that("small-n guard counts finite scores, not total length (NA padding does not lift it)", {
  # 19 finite values padded with 10 NAs: total length 29 (>= 20) but finite count 19 (< 20). A
  # mutant comparing length(values) to min_n would resolve; the spec counts finite scores only.
  padded <- c(as.numeric(1:19), rep(NA_real_, 10L))
  expect_warning(res <- resolve_cutoff(padded, "upper", fpr = 0.05),
                 class = "cier_warning_insufficient_items")
  expect_identical(res, NA_real_)
})

test_that("a single finite value (was 5/7) now abstains, both directions", {
  # Pinned change from the earlier behaviour (deliberately replaces the old `-> 5` / `-> 7` pin):
  # one finite value, or ten, is below the 20-score floor at fpr = 0.05.
  expect_warning(res1 <- resolve_cutoff(c(5, NA), "lower", fpr = 0.05),
                 class = "cier_warning_insufficient_items")
  expect_identical(res1, NA_real_)
  expect_warning(res2 <- resolve_cutoff(rep(7, 10L), "upper", fpr = 0.05),
                 class = "cier_warning_insufficient_items")
  expect_identical(res2, NA_real_)
})

test_that("a constant distribution abstains even above the small-n floor", {
  # 25 finite values clears the small-n floor (>= 20), isolating the constant guard: a high/low
  # quantile of a constant equals that constant, so the cutoff flags EVERY respondent (upper:
  # cutoff <= min; lower: cutoff >= max). Abstain, do not flag 100%. Both directions.
  expect_warning(up <- resolve_cutoff(rep(7, 25L), "upper", fpr = 0.05),
                 class = "cier_warning_insufficient_items")
  expect_identical(up, NA_real_)
  expect_warning(lo <- resolve_cutoff(rep(7, 25L), "lower", fpr = 0.05),
                 class = "cier_warning_insufficient_items")
  expect_identical(lo, NA_real_)
})

test_that("the constant guard fires as insufficient_items, never the saturation warning", {
  # Ordering guard: a wholly-constant input is the all-tie case (constant), not a partial tie mass
  # (saturation). It must raise insufficient_items and abstain, NOT saturated_cutoff. (A
  # guard-ordering mutant would resolve + warn saturated.)
  seen <- new.env()
  seen$classes <- character(0)
  withCallingHandlers(
    res <- resolve_cutoff(rep(2, 30L), "upper", fpr = 0.05),
    warning = function(w) {
      seen$classes <- c(seen$classes, class(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(res, NA_real_)
  expect_true("cier_warning_insufficient_items" %in% seen$classes)
  expect_false("cier_warning_saturated_cutoff" %in% seen$classes)
})

test_that("a partial tie mass at the extreme warns but resolves the cutoff", {
  # Upper: 5 of 40 values tie at the maximum (100), > 5%, so the 95th percentile lands ON the max.
  # The cutoff IS resolvable (finite), but flagging >= it flags the whole 12.5% tie mass -- more
  # than fpr. Warn, do not abstain.
  up_vals <- c(as.numeric(1:35), rep(100, 5L))                # n = 40
  expect_warning(cut_up <- resolve_cutoff(up_vals, "upper", fpr = 0.05),
                 class = "cier_warning_saturated_cutoff")
  expect_identical(cut_up, 100)                               # == max, finite
  expect_identical(cut_up,
                   as.numeric(stats::quantile(up_vals, 0.95, names = FALSE,
                                              type = 7L)))
  # Lower mirror: 5 of 40 tie at the minimum (0); the 5th percentile lands on min.
  lo_vals <- c(rep(0, 5L), as.numeric(1:35))                  # n = 40
  expect_warning(cut_lo <- resolve_cutoff(lo_vals, "lower", fpr = 0.05),
                 class = "cier_warning_saturated_cutoff")
  expect_identical(cut_lo, 0)                                 # == min, finite
})

test_that("healthy continuous data resolves silently (no saturation warning)", {
  # Continuous scores with no tie at the extreme: the quantile interpolates strictly inside the
  # range, so neither the constant nor the saturation guard fires.
  vals <- withr::with_seed(20260612L, stats::runif(60L))
  expect_no_warning(out <- resolve_cutoff(vals, "upper", fpr = 0.05))
  expect_false(is.na(out))
  expect_gt(out, min(vals))
  expect_lt(out, max(vals))
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
# The independent oracle (ref_kneedle) re-derives the parameter-free elbow and never calls the
# production kernel; oracle-only trust, tolerance 0 (kernel and oracle do byte-identical float
# ops). cier ships only the convex + increasing case (high Laz.R = careless = upper tail).

source(test_path("..", "reference", "ref-kneedle-satopaa-2011.R"))

test_that("kneedle_knee equals a hand-computed elbow (independent of the oracle)", {
  # Sorted v = c(0, 0, 0, 0, 1, 2, 5, 9), n = 8. x = (0:7)/7; y = v / 9. The deviations y - x are
  # 0, -.143, -.286, -.429, -.460, -.492, -.302, 0, so the most-negative (the convex knee) is
  # index 6 -> value 2. Computed by hand, so a mutant drifting the normalisation, argmin, or sort
  # fails here even if it still matches the oracle.
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
  # A moderate bulk (40 values 0.30..0.34) then a sharp careless spike (10 values 0.90..1.00). The
  # convex elbow sits at the top of the bulk / start of the gap -- at or above every bulk value and
  # strictly below the spike -- so flagging at/above it isolates the careless cluster, never the
  # whole bulk and never nobody. The knee equals max(bulk) exactly here (verified), so pin `>=` not
  # `>`; an argmax / inverted-elbow mutant lands at min(bulk) and fails.
  bulk <- seq(0.30, 0.34, length.out = 40L)
  tail <- seq(0.90, 1.00, length.out = 10L)
  knee <- kneedle_knee(c(bulk, tail))
  expect_gte(knee, max(bulk))                            # at/above the bulk
  expect_lt(knee, min(tail))                             # strictly below the spike
})

test_that("kneedle_knee pins the (i-1)/(n-1) normalisation (off-by-one guard)", {
  # A convex increasing curve with an INTERIOR elbow at index 4 -> value 4.591 (verified in R). The
  # (i-1)/n off-by-one moves the argmin to the left endpoint (2.287) and the concave/argmax mutant
  # gives 3.484, so this fixture kills the normalisation off-by-one that every smooth / right-skewed
  # parity fixture is degenerate against. (The i/(n-1) variant uniformly shifts the x-axis: it can
  # never change the argmin, so it is an equivalent mutant needing no fixture.)
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
  # A deliberate deviation from the predecessor (which returned the constant value, flagging
  # everyone): a knee is undefined on a degenerate distribution, so the cutoff abstains.
  expect_warning(res <- resolve_kneedle_cutoff(rep(0.42, 50L)),
                 class = "cier_warning_insufficient_items")
  expect_identical(res, NA_real_)
})

test_that("resolve_kneedle_cutoff is deterministic (no RNG, no jitter)", {
  vals <- withr::with_seed(11L, sort(rexp(60L)))
  expect_identical(resolve_kneedle_cutoff(vals), resolve_kneedle_cutoff(vals))
})

test_that("the kneedle elbow differs from the fpr=0.05 percentile (a distinct rule)", {
  # A 40-point bulk at 0.30 then a graded 20-point ramp 0.50..1.00. The elbow sits at the bulk/ramp
  # boundary (~0.3-0.5); the 95th percentile sits high in the ramp (~0.92). Different cutoff rules,
  # so a mutant aliasing kneedle to the percentile cannot survive: the elbow is strictly the lower.
  # The elbow lands on the bulk (0.30), flagging 100% -- a saturation warning (the diagnostic under
  # test below); muffled here, where only the cutoff is asserted.
  v <- c(rep(0.30, 40L), seq(0.50, 1.00, length.out = 20L))
  expect_lt(suppressWarnings(resolve_kneedle_cutoff(v)),
            resolve_cutoff(v, "upper", fpr = 0.05))
})

# ---- kneedle saturation diagnostic (the elbow flags a majority) --------------
# The percentile path surfaces an over-target realised rate (the saturation guard); the kneedle
# path did not. When the resolved elbow would flag > 50% of scored respondents -- a strongly
# bimodal sample whose knee sits at/below the careful bulk -- resolve_kneedle_cutoff() now WARNS
# (reusing cier_warning_saturated_cutoff) but still RESOLVES the elbow. It does NOT abstain: a
# true-majority low-quality panel should still be flagged. The realised rate is over the finite
# scores only, as mean(finite >= knee) (apply_flag's `>=`).

test_that("resolve_kneedle_cutoff warns but resolves when the elbow flags a majority", {
  # A bottom plateau at 0.30 places the convex elbow ON the plateau (knee = 0.30 = min), so flagging
  # >= it flags 100%. The flags-everyone case must still RESOLVE (warn, do not abstain), unlike the
  # percentile constant-guard abstain. mean(v >= 0.30) = 1.0 (warn); mean(v > 0.30) = 0.4 (< 0.5),
  # so a strict-`>` rate mutant would not warn here and dies on the expect_warning.
  v <- c(rep(0.30, 6L), 0.7, 0.8, 0.9, 1.0)
  expect_warning(res <- resolve_kneedle_cutoff(v),
                 class = "cier_warning_saturated_cutoff")
  expect_identical(res, 0.30)                                          # resolves to the elbow
  expect_identical(res, ref_kneedle(v, "convex", "increasing")$value)  # independent oracle
  expect_false(is.na(res))                                             # NOT an abstention
  expect_gt(mean(v >= res), 0.5)                                       # the realised majority
})

test_that("resolve_kneedle_cutoff does not warn at exactly 50% flagged (strict > 0.5)", {
  # A convex curve whose elbow flags exactly 5 / 10 = 50%. The strict `> 0.5` rule must NOT warn
  # here; a `>= 0.5` mutant would. (A symmetric split is not saturation -- the kneedle reading is
  # doing exactly its job.)
  v <- c(0, 1, 2, 3, 4, 5, 20, 40, 70, 110)
  knee <- ref_kneedle(v, "convex", "increasing")$value
  expect_identical(knee, 5)
  expect_identical(mean(v >= knee), 0.5)                              # exactly half
  expect_no_warning(res <- resolve_kneedle_cutoff(v),
                    class = "cier_warning_saturated_cutoff")
  expect_identical(res, 5)                                            # still resolves
})

test_that("resolve_kneedle_cutoff warns when the elbow flags just over half (60%)", {
  # The threshold is at 50%, not some higher share: an elbow flagging 6 / 10 = 60% warns. With the
  # boundary test above this pins the cut at (50%, 60%], killing a too-high threshold mutant
  # (e.g. `> 0.65`).
  v <- c(0, 1, 2, 3, 8, 30, 55, 85, 120, 160)
  knee <- ref_kneedle(v, "convex", "increasing")$value
  expect_identical(knee, 8)
  expect_identical(mean(v >= knee), 0.6)
  expect_warning(res <- resolve_kneedle_cutoff(v),
                 class = "cier_warning_saturated_cutoff")
  expect_identical(res, 8)
  # Pin the count / percent arithmetic on a NON-100% share: a hardcoded "100%" or a mis-rounded
  # percent would survive the all-flagged plateau message test alone.
  w <- tryCatch(resolve_kneedle_cutoff(v),
                cier_warning_saturated_cutoff = function(w) w)
  msg <- gsub("\\s+", " ", cli::ansi_strip(rlang::cnd_message(w)))
  expect_match(msg, "(8)", fixed = TRUE)             # the elbow value reported
  expect_match(msg, "6 of 10 scored", fixed = TRUE)  # count over the finite scores
  expect_match(msg, "60%", fixed = TRUE)             # percent of scored respondents
})

test_that("kneedle saturation reports finite-score counts and points to fpr", {
  # The realised rate / count is over the FINITE scores only: non-finite values are dropped before
  # both the elbow and the rate. With two non-finite values prepended to the 10-value plateau, the
  # message must read "10 of 10" (scored count), not "10 of 12" (vector length). The wording reports
  # count + percent and steers to the percentile `fpr` cutoff.
  v <- c(NA_real_, Inf, rep(0.30, 6L), 0.7, 0.8, 0.9, 1.0)
  w <- tryCatch(resolve_kneedle_cutoff(v),
                cier_warning_saturated_cutoff = function(w) w)
  msg <- gsub("\\s+", " ", cli::ansi_strip(rlang::cnd_message(w)))
  expect_match(msg, "10 of 10 scored", fixed = TRUE)   # finite-only denominator
  expect_match(msg, "100%", fixed = TRUE)              # the realised majority
  expect_match(msg, "more than half", fixed = TRUE)
  expect_match(msg, "fpr", fixed = TRUE)               # the remedy
})

test_that("saturation_pct never reports 100% unless every scored respondent is flagged", {
  # round(100 * n_flagged / n_scored) would show "100%" for 999 / 1000 (99.9%), contradicting the
  # "999 of 1000" count beside it; the display reaches 100 ONLY when n_flagged == n_scored. Exact
  # mid-shares round normally.
  expect_identical(saturation_pct(10L, 10L), 100)      # everyone -> 100%
  expect_identical(saturation_pct(6L, 10L), 60)        # exact half-plus -> 60%
  expect_identical(saturation_pct(999L, 1000L), 99)    # 99.9% -> 99%, never 100%
  expect_identical(saturation_pct(199L, 200L), 99)     # 99.5% rounds to 100 -> capped
})
