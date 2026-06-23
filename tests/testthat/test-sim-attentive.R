# Tests for the attentive GRM generator -- the careful-respondent layer of
# cier_simulate(). Pure internal kernels, called directly in the package namespace.
#
# Trust model: oracle-only (a generator, not an index, no CRAN partner). The oracle
# (ref-sim-attentive-grm.R) re-derives the preset pmf, categorisation, and implied
# marginal from scratch, never calling production. Fast tier: closed-form marginals
# (pnorm of resolved thresholds == target pmf at 1e-12, no RNG), categorisation tol
# 0, API / validator / shape edges, deterministic structure, seed reproducibility.
# Slow tier (skip_if_slow + skip_on_cran): large-n marginal-frequency oracle
# (n = 10000, band 0.03), factor-structure recovery, non-normal trait properties.

source(test_path("..", "reference", "ref-sim-attentive-grm.R"))

# ---- Local builders ---------------------------------------------------------

# Raw items data.frame in the lean min/max/scale schema (what sim_attentive()
# validates internally). `max` / `min` / `reverse_keyed` recycle.
sim_items <- function(scale, max, min = NULL, reverse_keyed = NULL) {
  df <- data.frame(scale = as.character(scale), max = max,
                   stringsAsFactors = FALSE)
  if (!is.null(min)) df$min <- min
  if (!is.null(reverse_keyed)) df$reverse_keyed <- reverse_keyed
  df
}

# The validated items list the lower kernels consume (scale, reverse_keyed, min,
# max, categories), via the production validator.
vit <- function(scale, max, min = NULL, reverse_keyed = NULL) {
  it <- sim_items(scale, max, min, reverse_keyed)
  check_items_simulate(it, nrow(it))
}

# Population moments (skewness, Pearson kurtosis = 3 for normal) for the
# non-normal trait property checks.
pop_skew <- function(x) {
  d <- x - mean(x)
  mean(d^3) / mean(d^2)^1.5
}
pop_kurt <- function(x) {
  d <- x - mean(x)
  mean(d^4) / mean(d^2)^2
}

# The full preset set in one place (production switch and oracle enumerate it
# independently).
sim_presets <- c("uniform", "peaked", "skewed_left", "skewed_right", "bimodal")

# =============================================================================
# FAST TIER
# =============================================================================

# ---- check_items_simulate ---------------------------------------------------

test_that("check_items_simulate returns the normalised list with categories", {
  it <- sim_items(c("E", "E", "A"), max = c(5, 5, 7), min = c(1, 1, 1))
  out <- check_items_simulate(it, 3L)
  expect_identical(out$scale, c("E", "E", "A"))
  expect_identical(out$min, c(1, 1, 1))
  expect_identical(out$max, c(5, 5, 7))
  expect_identical(out$categories, c(5L, 5L, 7L))   # max - min + 1
  expect_identical(out$reverse_keyed, c(FALSE, FALSE, FALSE))
})

test_that("check_items_simulate defaults min to 1 when the column is absent", {
  out <- check_items_simulate(sim_items(rep("F", 3L), max = 5), 3L)
  expect_identical(out$min, rep(1L, 3L))            # NOT required; defaults to 1
  expect_identical(out$categories, rep(5L, 3L))
})

test_that("check_items_simulate allows a single scale (min_scales = 1)", {
  out <- check_items_simulate(sim_items(rep("F", 4L), max = 5), 4L)
  expect_identical(length(unique(out$scale)), 1L)
})

test_that("check_items_simulate requires max on every item", {
  expect_error(check_items_simulate(sim_items(rep("F", 3L), max = c(5, 5, NA)), 3L),
               class = "cier_error_input")
  # max column absent entirely is also an error (unlike the index wrappers).
  bare <- data.frame(scale = rep("F", 3L), stringsAsFactors = FALSE)
  expect_error(check_items_simulate(bare, 3L), class = "cier_error_input")
})

test_that("check_items_simulate rejects max < min + 1 (need two options)", {
  expect_error(check_items_simulate(sim_items("F", max = 1, min = 1), 1L),
               class = "cier_error_input")
})

test_that("check_items_simulate rejects a range too large to enumerate", {
  # max - min + 1 overflowing the integer type must be a typed error, not an NA
  # categories count that crashes downstream untyped.
  expect_error(check_items_simulate(sim_items("F", max = 3e9, min = 0), 1L),
               class = "cier_error_input")
})

test_that("check_items_simulate accepts 0-based and bipolar bases", {
  zero <- check_items_simulate(sim_items(rep("F", 2L), max = 4, min = 0), 2L)
  expect_identical(zero$categories, rep(5L, 2L))    # 0..4 = five options
  bip <- check_items_simulate(sim_items(rep("F", 2L), max = 3, min = -3), 2L)
  expect_identical(bip$categories, rep(7L, 2L))     # -3..3
})

test_that("check_items_simulate allows heterogeneous ranges across items", {
  it <- sim_items(c("F", "F", "F"), max = c(5, 7, 4), min = c(1, 1, 0))
  out <- check_items_simulate(it, 3L)
  expect_identical(out$categories, c(5L, 7L, 5L))
})

test_that("check_items_simulate rejects non-whole / NA / missing min and scale", {
  expect_error(check_items_simulate(sim_items(rep("F", 2L), max = 5, min = 1.5), 2L),
               class = "cier_error_input")
  expect_error(check_items_simulate(sim_items(rep("F", 2L), max = 5, min = NA), 2L),
               class = "cier_error_input")
  no_scale <- data.frame(max = c(5, 5))
  expect_error(check_items_simulate(no_scale, 2L), class = "cier_error_input")
})

test_that("check_items_simulate rejects a non-frame or mis-sized items", {
  expect_error(check_items_simulate(list(scale = "F", max = 5), 1L),
               class = "cier_error_input")
  expect_error(check_items_simulate(sim_items(rep("F", 2L), max = 5), 3L),
               class = "cier_error_input")
})

# ---- Preset pmf (sim_preset_pmf) --------------------------------------------

test_that("sim_preset_pmf matches the oracle for every preset and K", {
  for (preset in sim_presets) {
    for (k in c(2L, 3L, 5L, 7L)) {
      expect_equal(sim_preset_pmf(preset, k), ref_grm_pmf(preset, k),
                   tolerance = 1e-12,
                   info = paste(preset, k))
    }
  }
})

test_that("sim_preset_pmf reproduces the agreed reference vectors (K = 2, 4, 5)", {
  # K = 5.
  expect_equal(sim_preset_pmf("peaked", 5L), c(1, 2, 3, 2, 1) / 9, tolerance = 1e-12)
  expect_equal(sim_preset_pmf("bimodal", 5L), c(3, 2, 1, 2, 3) / 11, tolerance = 1e-12)
  expect_equal(sim_preset_pmf("uniform", 5L), rep(0.2, 5L), tolerance = 1e-12)
  # Geometric skew, default strength 1.5: skewed_right is the exact mirror.
  sl <- sim_preset_pmf("skewed_left", 5L)
  expect_equal(sl, 1.5^(0:4) / sum(1.5^(0:4)), tolerance = 1e-12)
  expect_equal(sim_preset_pmf("skewed_right", 5L), rev(sl), tolerance = 1e-12)
  expect_gt(which.max(sl), 4L)                      # mode on the agreement end
  # Independent EVEN-K literals: the oracle shares production's closed form, so
  # these hand-written vectors are the only independent pin on the even-K shape.
  expect_equal(sim_preset_pmf("peaked", 4L), c(1, 2, 2, 1) / 6, tolerance = 1e-12)
  expect_equal(sim_preset_pmf("bimodal", 4L), c(2, 1, 1, 2) / 6, tolerance = 1e-12)
  expect_equal(sim_preset_pmf("peaked", 2L), c(0.5, 0.5), tolerance = 1e-12)
  expect_equal(sim_preset_pmf("bimodal", 2L), c(0.5, 0.5), tolerance = 1e-12)
  expect_equal(sim_preset_pmf("skewed_left", 2L), c(0.4, 0.6), tolerance = 1e-12)
})

test_that("sim_preset_pmf strength steepens the skew (per-spec overridable)", {
  expect_equal(sim_preset_pmf("skewed_left", 5L, strength = 2),
               c(1, 2, 4, 8, 16) / 31, tolerance = 1e-12)
  expect_equal(sim_preset_pmf("skewed_left", 5L, strength = 1),
               rep(0.2, 5L), tolerance = 1e-12)     # s = 1 collapses to uniform
})

test_that("every preset pmf is strictly positive and sums to one", {
  for (preset in sim_presets) {
    p <- sim_preset_pmf(preset, 6L)
    expect_true(all(p > 0))
    expect_equal(sum(p), 1, tolerance = 1e-12)
  }
})

test_that("sim_preset_pmf rejects an unknown preset name", {
  expect_error(sim_preset_pmf("agreement", 5L), class = "cier_error_input")
})

# ---- Marginals -> thresholds (exact closed form, no RNG) ---------------------

test_that("resolved thresholds invert exactly to the target pmf (the oracle)", {
  # tau_k = qnorm(cumsum(p)); feeding tau back through Phi must return p.
  for (preset in sim_presets) {
    for (k in c(2L, 4L, 5L)) {
      p <- ref_grm_pmf(preset, k)
      thr <- sim_marginals_to_thresholds(list(p))[[1L]]
      expect_length(thr, k - 1L)                    # drop the qnorm(1) = Inf cut
      expect_true(all(is.finite(thr)))
      expect_true(all(diff(thr) > 0))               # strictly increasing
      expect_equal(ref_marginal_implied(thr), p, tolerance = 1e-12,
                   info = paste(preset, k))
    }
  }
})

# ---- Marginals resolver (sim_resolve_marginals) -----------------------------

test_that("a NULL marginals spec defaults to peaked on every item", {
  it <- vit(rep("F", 3L), max = 5)
  out <- sim_resolve_marginals(NULL, it)
  expect_length(out, 3L)
  for (p in out) expect_equal(p, ref_grm_pmf("peaked", 5L), tolerance = 1e-12)
})

test_that("a single preset string applies to every item", {
  it <- vit(rep("F", 3L), max = 5)
  out <- sim_resolve_marginals("skewed_left", it)
  for (p in out) expect_equal(p, ref_grm_pmf("skewed_left", 5L), tolerance = 1e-12)
})

test_that("a single list spec carries per-spec strength to every item", {
  it <- vit(rep("F", 2L), max = 5)
  out <- sim_resolve_marginals(list(preset = "skewed_left", strength = 2), it)
  for (p in out) expect_equal(p, ref_grm_pmf("skewed_left", 5L, 2), tolerance = 1e-12)
})

test_that("a preset-keyed list is a single spec, not a per-scale map, across scales", {
  # Disambiguation pin: a list with a `preset` element is ONE spec applied to every
  # item, even on a MULTI-scale battery (a per-scale misread would demand E/N keys
  # and error). A bare list(preset=) with no strength uses the default.
  it <- vit(c("E", "N"), max = 5)
  out <- sim_resolve_marginals(list(preset = "skewed_left"), it)
  expect_equal(out[[1L]], ref_grm_pmf("skewed_left", 5L), tolerance = 1e-12)
  expect_equal(out[[2L]], ref_grm_pmf("skewed_left", 5L), tolerance = 1e-12)
})

test_that("a single numeric pmf vector applies to a homogeneous battery", {
  it <- vit(rep("F", 3L), max = 5)
  v <- c(0.1, 0.2, 0.4, 0.2, 0.1)
  out <- sim_resolve_marginals(v, it)
  for (p in out) expect_equal(p, v, tolerance = 1e-12)
})

test_that("a scale-named list resolves per scale, allowing different K", {
  it <- vit(c("E", "E", "N"), max = c(5, 5, 7))
  out <- sim_resolve_marginals(list(E = "peaked", N = "skewed_left"), it)
  expect_equal(out[[1L]], ref_grm_pmf("peaked", 5L), tolerance = 1e-12)
  expect_equal(out[[2L]], ref_grm_pmf("peaked", 5L), tolerance = 1e-12)
  expect_equal(out[[3L]], ref_grm_pmf("skewed_left", 7L), tolerance = 1e-12)
})

test_that("an unnamed length-p list resolves per item", {
  it <- vit(c("F", "F", "F"), max = 5)
  out <- sim_resolve_marginals(list("peaked", "uniform", "bimodal"), it)
  expect_equal(out[[1L]], ref_grm_pmf("peaked", 5L), tolerance = 1e-12)
  expect_equal(out[[2L]], ref_grm_pmf("uniform", 5L), tolerance = 1e-12)
  expect_equal(out[[3L]], ref_grm_pmf("bimodal", 5L), tolerance = 1e-12)
})

test_that("per-item specs may mix preset strings and explicit pmf vectors", {
  it <- vit(c("F", "F"), max = 5)
  out <- sim_resolve_marginals(list("peaked", c(0.4, 0.3, 0.1, 0.1, 0.1)), it)
  expect_equal(out[[1L]], ref_grm_pmf("peaked", 5L), tolerance = 1e-12)
  expect_equal(out[[2L]], c(0.4, 0.3, 0.1, 0.1, 0.1), tolerance = 1e-12)
})

# ---- Marginals validation edges ---------------------------------------------

test_that("a pmf vector of the wrong length is a typed input error", {
  it <- vit(rep("F", 2L), max = 5)
  expect_error(sim_resolve_marginals(c(0.5, 0.5), it), class = "cier_error_input")
})

test_that("a single pmf vector on a heterogeneous battery errors on the mismatch", {
  it <- vit(c("F", "F"), max = c(5, 7))
  expect_error(sim_resolve_marginals(c(0.1, 0.2, 0.4, 0.2, 0.1), it),
               class = "cier_error_input")
})

test_that("a negative, zero, or non-summing pmf is a typed input error", {
  it <- vit("F", max = 5)
  expect_error(sim_resolve_marginals(list(c(-0.1, 0.3, 0.3, 0.3, 0.2)), it),
               class = "cier_error_input")
  expect_error(sim_resolve_marginals(list(c(0, 0.25, 0.25, 0.25, 0.25)), it),
               class = "cier_error_input")            # zero -> infinite threshold
  expect_error(sim_resolve_marginals(list(c(0.2, 0.2, 0.2, 0.2, 0.1)), it),
               class = "cier_error_input")            # sums to 0.9
})

test_that("an unknown preset string or bad strength is a typed input error", {
  it <- vit("F", max = 5)
  expect_error(sim_resolve_marginals("nope", it), class = "cier_error_input")
  expect_error(sim_resolve_marginals(list(preset = "skewed_left", strength = 0.5), it),
               class = "cier_error_input")             # strength < 1
  expect_error(sim_resolve_marginals(list(preset = "skewed_left", strength = NA), it),
               class = "cier_error_input")
})

test_that("a scale-named list must cover exactly the scales", {
  it <- vit(c("E", "N"), max = 5)
  expect_error(sim_resolve_marginals(list(E = "peaked"), it),
               class = "cier_error_input")             # N missing
  expect_error(sim_resolve_marginals(list(E = "peaked", X = "uniform"), it),
               class = "cier_error_input")             # X is not a scale
})

test_that("an unnamed list of the wrong length is a typed input error", {
  it <- vit(rep("F", 3L), max = 5)
  expect_error(sim_resolve_marginals(list("peaked", "uniform"), it),
               class = "cier_error_input")             # length 2 != 3 items
})

test_that("a partially-named list is rejected, never read positionally", {
  # A forgotten name must NOT silently degrade to positional mapping.
  it <- vit(c("E", "N"), max = 5)
  expect_error(sim_resolve_marginals(list(E = "peaked", "uniform"), it),
               class = "cier_error_input")
})

test_that("a scale-named list with a duplicated scale is a typed input error", {
  # setequal() ignores duplicates; an explicit guard stops the second spec being
  # silently dropped (first-wins).
  it <- vit(c("E", "N"), max = 5)
  spec <- list(E = "peaked", E = "bimodal", N = "uniform")
  expect_error(sim_resolve_marginals(spec, it), class = "cier_error_input")
})

test_that("an unknown key in a preset-keyed spec is a typed input error", {
  # A typo'd strength (or any extra key) must not silently fall back to default.
  it <- vit("F", max = 5)
  spec <- list(preset = "skewed_left", strenght = 2)
  expect_error(sim_resolve_marginals(spec, it), class = "cier_error_input")
})

test_that("a pmf whose sum is a few ulp above 1 still yields finite thresholds", {
  # The validator's 1e-8 sum tolerance admits c(0.25, 0.75, 1e-9); normalisation in
  # the threshold transform keeps the cuts finite and non-erroring (no qnorm(>1) =
  # NaN crash, no qnorm(1) = Inf unreachable category).
  thr <- sim_marginals_to_thresholds(list(c(0.25, 0.75, 1e-9)))[[1L]]
  expect_length(thr, 2L)
  expect_true(all(is.finite(thr)))
  expect_true(all(diff(thr) > 0))
})

# ---- Threshold dispatch & mutual exclusion ----------------------------------

test_that("sim_resolve_thresholds defaults to peaked when nothing is supplied", {
  it <- vit(rep("F", 2L), max = 5)
  thr <- sim_resolve_thresholds(it, marginals = NULL, thresholds = NULL)
  expected <- sim_marginals_to_thresholds(sim_resolve_marginals(NULL, it))
  expect_equal(thr, expected, tolerance = 1e-12)
})

test_that("explicit thresholds pass through after validation", {
  it <- vit("F", max = 5)
  user <- list(c(-2, -0.5, 0.5, 2))
  expect_identical(sim_resolve_thresholds(it, NULL, user), user)
})

test_that("supplying both marginals and thresholds is a typed input error", {
  it <- vit("F", max = 5)
  user <- list(c(-2, -0.5, 0.5, 2))
  expect_error(sim_resolve_thresholds(it, marginals = "uniform", thresholds = user),
               class = "cier_error_input")
  expect_error(sim_attentive(5L, sim_items("F", max = 5),
                             marginals = "uniform", thresholds = user),
               class = "cier_error_input")
})

test_that("marginals and thresholds resolve end-to-end through sim_attentive", {
  # Routes a NON-NULL marginals= AND a non-NULL thresholds= through the full chain.
  # A marginals/thresholds argument swap at the dispatch would either error here or
  # land the wrong thresholds in $metadata.
  it <- sim_items(c("F", "F"), max = 5)
  meta <- sim_attentive_with_metadata(10L, it, marginals = "uniform")$metadata
  expect_equal(meta$thresholds[[1L]],
               sim_marginals_to_thresholds(list(ref_grm_pmf("uniform", 5L)))[[1L]],
               tolerance = 1e-12)
  user <- list(c(-2, -0.5, 0.5, 2), c(-1.5, -0.5, 0.5, 1.5))
  m2 <- sim_attentive_with_metadata(10L, it, thresholds = user)$metadata
  expect_identical(m2$thresholds, user)
})

test_that("user loadings change the generated factor structure (fast)", {
  # Pins that loadings flow through: a mutant substituting the 0.7 defaults could
  # not separate these. Higher loadings -> stronger within-scale correlation, same
  # seed.
  it <- sim_items(c("F", "F"), max = 5)
  hi <- withr::with_seed(5L, sim_attentive(400L, it, loadings = matrix(0.9, 2L, 1L)))
  lo <- withr::with_seed(5L, sim_attentive(400L, it, loadings = matrix(0.2, 2L, 1L)))
  expect_gt(stats::cor(hi[, 1L], hi[, 2L]), stats::cor(lo[, 1L], lo[, 2L]))
})

# ---- GRM categorisation (sim_grm_categorise) --------------------------------

test_that("sim_grm_categorise matches the independent findInterval oracle", {
  eta <- matrix(c(-1.5, 0.5, 1.5, -0.2, 0.0, 0.7), nrow = 3L)
  thr <- list(c(-1, 0, 1), c(-0.5, 0.5))
  out <- sim_grm_categorise(eta, thr)
  expect_identical(out, ref_grm_categorise(eta, thr))
  expect_identical(out[, 1L], c(1L, 3L, 4L))
  expect_identical(out[, 2L], c(2L, 2L, 3L))
})

test_that("a latent value exactly on a threshold falls in the higher category", {
  # eta == tau -> category counts the threshold (left-closed), so cat = 3 here.
  out <- sim_grm_categorise(matrix(0, nrow = 1L), list(c(-1, 0, 1)))
  expect_identical(out[1L, 1L], 3L)
  expect_identical(out, ref_grm_categorise(matrix(0, nrow = 1L), list(c(-1, 0, 1))))
})

# ---- Loadings / factor_cor / thresholds validators --------------------------

test_that("sim_default_loadings builds one 0.7-loaded factor per scale", {
  it <- vit(c("E", "E", "E", "A", "A", "A"), max = 5)
  L <- sim_default_loadings(it)
  expect_identical(dim(L), c(6L, 2L))
  expect_identical(colnames(L), c("A", "E"))        # sorted scale labels
  expect_equal(unname(L[1L, "E"]), 0.7)
  expect_equal(unname(L[1L, "A"]), 0)
  expect_equal(unname(L[4L, "A"]), 0.7)
})

test_that("loadings validation rejects malformed matrices", {
  it <- sim_items(rep("F", 5L), max = 5)
  expect_error(sim_attentive(10L, it, loadings = matrix(0.7, 3L, 1L)),
               class = "cier_error_input")            # wrong row count
  expect_error(sim_attentive(10L, it, loadings = matrix("x", 5L, 1L)),
               class = "cier_error_input")            # non-numeric
  expect_error(sim_attentive(10L, it,
                             loadings = matrix(c(0.7, NA, rep(0.7, 3L)), 5L, 1L)),
               class = "cier_error_input")            # non-finite
  expect_error(sim_attentive(10L, it, loadings = matrix(numeric(0L), 5L, 0L)),
               class = "cier_error_input")            # zero columns
  expect_error(sim_attentive(10L, it, loadings = "x"),
               class = "cier_error_input")            # non-matrix
})

test_that("factor_cor validation rejects non-correlation matrices", {
  it <- sim_items(c("E", "E", "A", "A"), max = 5)
  expect_error(sim_attentive(10L, it, factor_cor = matrix(c(1, 0.4, 0.5, 1), 2L)),
               class = "cier_error_input")            # non-symmetric
  expect_error(sim_attentive(10L, it, factor_cor = matrix(c(1, 0.4, 0.4, 0.9), 2L)),
               class = "cier_error_input")            # non-unit diagonal
  expect_error(sim_attentive(10L, it, factor_cor = matrix(1, 2L, 2L)),
               class = "cier_error_input")            # non-positive-definite
  bad_dim <- diag(3L)
  expect_error(sim_attentive(10L, it, factor_cor = bad_dim),
               class = "cier_error_input")            # wrong dimension vs m = 2
  expect_error(sim_attentive(10L, it, factor_cor = "x"),
               class = "cier_error_input")            # non-matrix
  expect_error(sim_attentive(10L, it,
                             factor_cor = matrix(c(1, NA, NA, 1), 2L)),
               class = "cier_error_input")            # non-finite
})

test_that("a valid factor_cor carrying scale-label dimnames is accepted", {
  # diag() of a named matrix returns a NAMED vector; attribute-blind comparison must
  # not falsely reject it (naming factor_cor by scale is natural usage).
  it <- sim_items(c("A", "A", "E", "E"), max = 5)
  fc <- matrix(c(1, 0.4, 0.4, 1), 2L,
               dimnames = list(c("A", "E"), c("A", "E")))
  expect_no_error(withr::with_seed(1L, sim_attentive(20L, it, factor_cor = fc)))
})

test_that("loadings implying a communality above 1 are a typed input error", {
  # pmax(0, 1 - factor_var) would silently floor the residual at 0 and break the
  # unit-variance / marginals contract; reject instead.
  it <- sim_items("F", max = 5)
  expect_error(sim_attentive(10L, it, loadings = matrix(1.5, 1L, 1L)),
               class = "cier_error_input")
  # Two 0.8 cross-loadings on factors correlated 0.5: factor_var = 1.92 > 1.
  it2 <- sim_items(c("E", "A"), max = 5)
  fc <- matrix(c(1, 0.5, 0.5, 1), 2L)
  cross <- matrix(0.8, 2L, 2L)
  expect_error(sim_attentive(10L, it2, loadings = cross, factor_cor = fc),
               class = "cier_error_input")
})

test_that("user loadings with all-negative reverse rows are a typed keying error", {
  # reverse-keyed rows are negated internally, so CFA-signed (negative) reverse
  # loadings would double-negate to positive and silently generate FORWARD-keyed
  # data. Reject when every nonzero loading on the reverse rows is negative.
  it <- sim_items(rep("F", 4L), max = 5,
                  reverse_keyed = c(FALSE, FALSE, TRUE, TRUE))
  bad <- matrix(c(0.7, 0.7, -0.7, -0.7), 4L, 1L)
  expect_error(sim_attentive(20L, it, loadings = bad), class = "cier_error_input")
})

test_that("positive-magnitude reverse loadings pass (the documented contract)", {
  it <- sim_items(rep("F", 4L), max = 5,
                  reverse_keyed = c(FALSE, FALSE, TRUE, TRUE))
  expect_no_error(withr::with_seed(1L,
                                   sim_attentive(20L, it, loadings = matrix(0.7, 4L, 1L))))
})

test_that("a reverse row with a positive main and a negative cross-loading passes", {
  # mixed signs are not the double-negation trap; only ALL-negative reverse rows are.
  it <- sim_items(c("E", "E", "A", "A"), max = 5,
                  reverse_keyed = c(FALSE, FALSE, TRUE, FALSE))
  mixed <- rbind(c(0.7, 0), c(0.7, 0), c(-0.2, 0.7), c(0, 0.7))
  expect_no_error(withr::with_seed(1L, sim_attentive(20L, it, loadings = mixed)))
})

test_that("the guard counts NONZERO reverse loadings, sparing zero entries", {
  # reverse rows whose only nonzero loading is negative (zero main, negative cross)
  # ARE the trap; pins 'every NONZERO reverse loading < 0' against a cruder 'all
  # entries < 0', which a zero entry would spare.
  it <- sim_items(c("E", "E", "A", "A"), max = 5,
                  reverse_keyed = c(FALSE, FALSE, TRUE, TRUE))
  bad <- rbind(c(0.7, 0), c(0.7, 0), c(0, -0.7), c(0, -0.7))
  expect_error(sim_attentive(20L, it, loadings = bad), class = "cier_error_input")
})

test_that("a positive reverse row alongside a negative one passes (pooled reading)", {
  # the guard fires only when the WHOLE reverse block is negative; a single positive
  # reverse loading means the input is not the clean double-negation trap, so it
  # must NOT error (pooled, not per-row-any).
  it <- sim_items(rep("A", 4L), max = 5,
                  reverse_keyed = c(FALSE, FALSE, TRUE, TRUE))
  mixed <- matrix(c(0.7, 0.7, 0.7, -0.7), 4L, 1L)            # reverse rows: +0.7, -0.7
  expect_no_error(withr::with_seed(1L, sim_attentive(20L, it, loadings = mixed)))
})

test_that("the keying guard is inert without reverse-keyed items", {
  it <- sim_items(rep("F", 4L), max = 5)                     # no reverse_keyed column
  expect_no_error(withr::with_seed(1L,
                                   sim_attentive(20L, it, loadings = matrix(-0.7, 4L, 1L))))
})

test_that("thresholds validation rejects malformed lists", {
  it <- sim_items(rep("F", 2L), max = 5)
  expect_error(sim_attentive(10L, it, thresholds = "x"),
               class = "cier_error_input")            # non-list
  expect_error(sim_attentive(10L, it, thresholds = list(c(-1, 0, 1, 2))),
               class = "cier_error_input")            # wrong length (1 != 2)
  expect_error(sim_attentive(10L, it,
                             thresholds = list(c(-1, 0, 1), c(-1, 0, 1, 2))),
               class = "cier_error_input")            # wrong inner length
  expect_error(sim_attentive(10L, it,
                             thresholds = list(c(0, 0, 1, 2), c(-1, 0, 1, 2))),
               class = "cier_error_input")            # non-increasing
  expect_error(sim_attentive(10L, it,
                             thresholds = list(c(-1, 0, NA, 1), c(-1, 0, 1, 2))),
               class = "cier_error_input")            # NA entry
})

test_that("an unsupported trait distribution is rejected early", {
  it <- sim_items(rep("F", 3L), max = 5)
  expect_error(sim_attentive(10L, it, trait_distribution = "weird"),
               class = "cier_error_input")
})

test_that("unknown or misplaced trait_params keys are a typed input error", {
  it <- sim_items(rep("F", 3L), max = 5)
  # A typo'd key must not silently fall back to the default df = 5.
  expect_error(sim_attentive(10L, it, trait_distribution = "t",
                             trait_params = list(dfs = 20)),
               class = "cier_error_input")
  # A parameter for the wrong distribution (alpha belongs to skew_normal).
  expect_error(sim_attentive(10L, it, trait_distribution = "t",
                             trait_params = list(alpha = 8)),
               class = "cier_error_input")
  # normal takes no parameters.
  expect_error(sim_attentive(10L, it, trait_distribution = "normal",
                             trait_params = list(sep = 2)),
               class = "cier_error_input")
})

test_that("a t df in the unstable region (below 3) is a typed input error", {
  it <- sim_items(rep("F", 3L), max = 5)
  expect_error(sim_attentive(10L, it, trait_distribution = "t",
                             trait_params = list(df = 2.1)),
               class = "cier_error_input")
})

test_that("a non-positive n is a typed input error", {
  it <- sim_items(rep("F", 3L), max = 5)
  expect_error(sim_attentive(0L, it), class = "cier_error_input")
  expect_error(sim_attentive(2.5, it), class = "cier_error_input")
})

# ---- Shape, range, offset to min..max ---------------------------------------

test_that("sim_attentive returns an integer matrix in 1..max for 1-based items", {
  it <- sim_items(rep("F", 4L), max = 5)
  x <- withr::with_seed(11L, sim_attentive(50L, it))
  expect_true(is.matrix(x))
  expect_identical(storage.mode(x), "integer")
  expect_identical(dim(x), c(50L, 4L))
  expect_true(all(x >= 1L & x <= 5L))
})

test_that("sim_attentive offsets a 0-based item to min..max", {
  it <- sim_items(rep("F", 3L), max = 4, min = 0)
  x <- withr::with_seed(12L, sim_attentive(200L, it))
  expect_true(all(x >= 0L & x <= 4L))               # 0..4, not 1..5
  expect_true(any(x == 0L))                          # the base is reachable
})

test_that("sim_attentive honours heterogeneous and binary ranges", {
  it <- sim_items(c("F", "F", "F"), max = c(5, 3, 7), min = c(1, 1, 0))
  x <- withr::with_seed(13L, sim_attentive(300L, it))
  expect_true(all(x[, 1L] %in% 1:5))
  expect_true(all(x[, 2L] %in% 1:3))
  expect_true(all(x[, 3L] %in% 0:7))
  bin <- sim_items(rep("F", 2L), max = 2)            # 1/2 binary
  xb <- withr::with_seed(14L, sim_attentive(50L, bin))
  expect_true(all(xb %in% c(1L, 2L)))
})

test_that("a reverse-keyed item on a 0-based scale stays raw in min..max", {
  # Raw-orientation contract on a non-1 base: output is as-clicked (the reverse item
  # loads negatively but is NOT reverse-SCORED), offset is min - 1 + category. A
  # no-offset / min+category mutant leaves 1..5; a reflect-the-output mutant is
  # caught by the negative-correlation slow test.
  it <- sim_items(c("F", "F"), max = 4, min = 0,
                  reverse_keyed = c(FALSE, TRUE))
  x <- withr::with_seed(17L, sim_attentive(300L, it))
  expect_true(all(x >= 0L & x <= 4L))                # 0..4, both columns
  expect_true(any(x[, 2L] == 0L) && any(x[, 2L] == 4L))  # base + top reachable
})

test_that("sim_attentive supports a single item and a single respondent", {
  one_item <- withr::with_seed(15L, sim_attentive(20L, sim_items("F", max = 5)))
  expect_identical(dim(one_item), c(20L, 1L))
  one_resp <- withr::with_seed(16L, sim_attentive(1L, sim_items(rep("F", 4L), max = 5)))
  expect_identical(dim(one_resp), c(1L, 4L))
})

# ---- Deterministic structure (no RNG) ---------------------------------------

test_that("sim_effective_loadings flips the sign on reverse-keyed rows only", {
  it <- vit(c("F", "F", "F"), max = 5, reverse_keyed = c(FALSE, TRUE, FALSE))
  L <- matrix(c(0.7, 0.6, 0.5), nrow = 3L, ncol = 1L)
  eff <- sim_effective_loadings(L, it)
  expect_identical(eff[, 1L], c(0.7, -0.6, 0.5))     # only row 2 negated, tol 0
})

test_that("sim_residual_variance makes the implied latent variance one", {
  it <- vit(c("E", "E", "A", "A"), max = 5)
  L <- sim_default_loadings(it)
  fc <- diag(2L)
  eff <- sim_effective_loadings(L, it)
  psi <- sim_residual_variance(eff, fc)
  expect_equal(psi, pmax(0, 1 - diag(eff %*% fc %*% t(eff))), tolerance = 1e-12)
  total <- diag(eff %*% fc %*% t(eff)) + psi
  expect_equal(total, rep(1, 4L), tolerance = 1e-12) # unit-variance latent
})

# ---- Seed reproducibility ---------------------------------------------------

test_that("the same seed gives bytewise-identical responses; a different seed does not", {
  it <- sim_items(rep("F", 5L), max = 5)
  a <- withr::with_seed(101L, sim_attentive(40L, it))
  b <- withr::with_seed(101L, sim_attentive(40L, it))
  expect_identical(a, b)
  c <- withr::with_seed(102L, sim_attentive(40L, it))
  expect_false(identical(a, c))
})

test_that("reproducibility holds across every trait distribution", {
  it <- sim_items(rep("F", 4L), max = 5)
  for (dist in sim_trait_distributions()) {
    a <- withr::with_seed(7L, sim_attentive(30L, it, trait_distribution = dist))
    b <- withr::with_seed(7L, sim_attentive(30L, it, trait_distribution = dist))
    expect_identical(a, b, info = dist)
  }
})

test_that("the bimodal trait draw is a distinct, standardised distribution (fast guard)", {
  # bimodal is brand-new (no archive precedent); deep moment / two-camp checks are
  # slow tier. This cheap fast guard kills the cheap mutants: bimodal aliasing to
  # normal, and skipping the unit-variance standardisation. n = 2000 with a 0.06
  # band makes the kill structural, not a sampling accident -- the un-standardised
  # mixture has population sd sqrt(1.25) = 1.118 (outside 0.06) and the
  # divide-by-variance mutant has sd 1/sqrt(1.25) = 0.894 (also outside), while the
  # correct draw lands ~1 +/- 0.02. (t / skew_normal moments stay slow tier --
  # non-normal trait shape is treated as a property.)
  zb <- withr::with_seed(3L, sim_draw_traits(2000L, diag(1L), "bimodal"))
  zn <- withr::with_seed(3L, sim_draw_traits(2000L, diag(1L), "normal"))
  expect_false(identical(as.numeric(zb), as.numeric(zn)))  # not aliased to normal
  expect_equal(stats::sd(as.numeric(zb)), 1, tolerance = 0.06)   # standardised
})

# =============================================================================
# SLOW TIER (large-n distributional checks)
# =============================================================================

# ---- Large-n marginal-frequency oracle (n = 10000, band 0.03) ---------------

test_that("attentive marginals match every preset's target pmf at large n", {
  skip_on_cran()
  skip_if_slow()
  n <- 10000L
  it <- sim_items(rep("F", 4L), max = 5)             # single scale, default 0.7 loadings
  for (preset in sim_presets) {
    x <- withr::with_seed(20260611L, sim_attentive(n, it, marginals = preset))
    target <- ref_grm_pmf(preset, 5L)
    for (j in seq_len(ncol(x))) {
      observed <- as.numeric(table(factor(x[, j], levels = 1:5))) / n
      # The DOCUMENTED contract is an absolute per-category band, so assert max abs
      # deviation directly -- not expect_equal(), whose waldo tolerance is a (much
      # tighter) relative metric on a near-uniform pmf.
      expect_lt(max(abs(observed - target)), 0.03)
    }
  }
})

test_that("a heterogeneous-K battery recovers each item's target marginal", {
  skip_on_cran()
  skip_if_slow()
  n <- 10000L
  it <- sim_items(c("F", "F"), max = c(5, 7))
  x <- withr::with_seed(20260611L, sim_attentive(n, it, marginals = "peaked"))
  obs5 <- as.numeric(table(factor(x[, 1L], levels = 1:5))) / n
  obs7 <- as.numeric(table(factor(x[, 2L], levels = 1:7))) / n
  expect_lt(max(abs(obs5 - ref_grm_pmf("peaked", 5L))), 0.03)
  expect_lt(max(abs(obs7 - ref_grm_pmf("peaked", 7L))), 0.03)
})

# ---- Factor-structure recovery ----------------------------------------------

test_that("within-scale correlation exceeds between-scale correlation", {
  skip_on_cran()
  skip_if_slow()
  it <- sim_items(c(rep("E", 4L), rep("A", 4L)), max = 5)   # two orthogonal scales
  x <- withr::with_seed(20260611L, sim_attentive(8000L, it))
  cm <- stats::cor(x)
  within <- mean(c(cm[1:4, 1:4][lower.tri(cm[1:4, 1:4])],
                   cm[5:8, 5:8][lower.tri(cm[5:8, 5:8])]))
  between <- mean(cm[1:4, 5:8])
  expect_gt(within, 0.2)
  expect_gt(within, between + 0.12)
  expect_lt(abs(between), 0.08)                      # identity factor_cor -> ~0
})

test_that("a reverse-keyed item correlates negatively with its scale-mates", {
  skip_on_cran()
  skip_if_slow()
  # 0-based scale (min = 0): pins the raw-orientation sign flip on a non-1 base,
  # where a reflect-the-output mutant would invert the sign. Loadings 0.8 (latent
  # r = -0.64) clears polychoric attenuation with margin; the point is the sign, not
  # the magnitude.
  it <- sim_items(c("F", "F"), max = 4, min = 0, reverse_keyed = c(FALSE, TRUE))
  x <- withr::with_seed(20260611L, {
    sim_attentive(4000L, it, loadings = matrix(0.8, nrow = 2L, ncol = 1L))
  })
  expect_lt(stats::cor(x[, 1L], x[, 2L]), -0.3)      # raw orientation: flips
})

test_that("correlated factors recover a positive scale-score correlation", {
  skip_on_cran()
  skip_if_slow()
  it <- sim_items(c(rep("E", 6L), rep("A", 6L)), max = 5)
  fc <- matrix(c(1, 0.5, 0.5, 1), nrow = 2L)
  x <- withr::with_seed(20260611L, sim_attentive(4000L, it, factor_cor = fc))
  e_score <- rowMeans(x[, 1:6])
  a_score <- rowMeans(x[, 7:12])
  expect_gt(stats::cor(e_score, a_score), 0.2)
})

# ---- Non-normal trait properties --------------------------------------------

test_that("skew_normal traits are standardised and right-skewed", {
  skip_on_cran()
  skip_if_slow()
  z <- withr::with_seed(20260611L, sim_draw_traits(6000L, diag(1L), "skew_normal"))
  expect_equal(mean(z), 0, tolerance = 0.08)
  expect_equal(stats::sd(as.numeric(z)), 1, tolerance = 0.08)
  expect_gt(pop_skew(as.numeric(z)), 0.3)            # alpha = 5 -> positive skew
})

test_that("t traits are standardised with heavier-than-normal tails", {
  skip_on_cran()
  skip_if_slow()
  zt <- withr::with_seed(20260611L, sim_draw_traits(6000L, diag(1L), "t"))
  zn <- withr::with_seed(20260611L, sim_draw_traits(6000L, diag(1L), "normal"))
  expect_equal(stats::sd(as.numeric(zt)), 1, tolerance = 0.08)
  expect_gt(pop_kurt(as.numeric(zt)), 3.5)           # leptokurtic (normal = 3)
  expect_gt(mean(abs(zt) > 3), mean(abs(zn) > 3))
})

test_that("bimodal traits are standardised, platykurtic, and two-camp", {
  skip_on_cran()
  skip_if_slow()
  zb <- withr::with_seed(20260611L, sim_draw_traits(6000L, diag(1L), "bimodal"))
  zn <- withr::with_seed(20260611L, sim_draw_traits(6000L, diag(1L), "normal"))
  expect_equal(mean(zb), 0, tolerance = 0.08)
  expect_equal(stats::sd(as.numeric(zb)), 1, tolerance = 0.08)
  expect_lt(pop_kurt(as.numeric(zb)), 3)             # bimodal -> platykurtic
  # Two camps: a thinned centre relative to a unimodal normal of equal variance.
  expect_lt(mean(abs(zb) < 0.3), mean(abs(zn) < 0.3))
})

test_that("bimodal sep widens the gap between the camps", {
  skip_on_cran()
  skip_if_slow()
  narrow <- withr::with_seed(1L, {
    sim_draw_traits(6000L, diag(1L), "bimodal", params = list(sep = 0.5))
  })
  wide <- withr::with_seed(1L, {
    sim_draw_traits(6000L, diag(1L), "bimodal", params = list(sep = 2))
  })
  # Both standardised to unit variance, but the wider raw separation leaves a
  # deeper central trough once standardised.
  expect_lt(mean(abs(wide) < 0.3), mean(abs(narrow) < 0.3))
})
