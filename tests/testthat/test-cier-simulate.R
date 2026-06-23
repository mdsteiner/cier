# Tests for the public cier_simulate() orchestrator + the light cier_sim object.
# The orchestrator wires the attentive GRM layer, the careless plan/content engine,
# the span-coupled response times, and the optional attention checks into one
# validated object, under a frozen RNG draw order: attentive -> plan -> content ->
# times -> checks. Attentive-first makes two end-to-end pins possible: the attentive
# layer is invariant to `prevalence` at a fixed seed, and a speeder-only run's
# $responses is byte-identical to a prevalence-0 run (speeders are careless in TIME
# only).
#
# Trust model: oracle-only; engine behaviours are pinned end-to-end THROUGH the
# shipped indices (longstring, lazr, autocorrelation, mahalanobis, total_time,
# page_time, attention) plus bytewise seed reproducibility and a frozen-seed
# RNG-stream digest.

source(test_path("..", "reference", "ref-sim-patterns.R"))

# ---- Local builders ---------------------------------------------------------

# Raw items data.frame in the lean min/max/scale schema (what users pass and what
# $items must hand back verbatim).
sim_items_df <- function(scale, max = 5, min = NULL) {
  df <- data.frame(scale = as.character(scale), max = max,
                   stringsAsFactors = FALSE)
  if (!is.null(min)) df$min <- min
  df
}

it10 <- sim_items_df(rep(c("E", "A"), each = 5L))           # 10 items, 2 scales
it20 <- sim_items_df(rep(c("E", "A"), each = 10L))          # the recovery battery

# Text-serialised digest of the simulator output (the frozen-RNG-stream fixture):
# %g at 7 significant digits, locale-free, no platform EOL translation. 7 (not 17)
# is deliberate: the continuous response times come from exp() of the RNG draws, and
# exp()'s last bits differ across CPU architectures (Apple aarch64 vs x86 libm), so a
# 17-digit digest is not portable. 7 digits sits ~8 orders of magnitude above that
# last-bit noise yet far below a genuinely different draw (which moves the value
# wholesale), so it stays portable AND still catches any added/removed/reordered draw.
sim_digest <- function(sim) {
  fmt <- function(x) formatC(as.vector(x), format = "g", digits = 7)
  key <- c(fmt(sim$responses), fmt(sim$seconds), fmt(sim$page_seconds),
           as.character(sim$items_per_page),
           if (!is.null(sim$checks)) fmt(sim$checks),
           if (!is.null(sim$pass)) fmt(unlist(sim$pass)),
           as.character(sim$truth$careless), sim$truth$pattern,
           sim$truth$extent, as.character(sim$truth$onset_item),
           as.character(sim$truth$offset_item),
           as.character(sim$truth$speeded))
  payload <- paste0(paste(key, collapse = "\n"), "\n")
  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)
  writeBin(charToRaw(payload), tmp)
  unname(tools::md5sum(tmp))
}

# =============================================================================
# FAST TIER
# =============================================================================

# ---- Default pattern mixture (Welz & Alfons equal allocation) ---------------

test_that("the default pattern mixture is equal weights over the eight names", {
  expect_identical(
    sim_default_patterns(),
    stats::setNames(rep(1 / 8, 8L),
                    c("random", "straightline", "midpoint", "extreme",
                      "diagonal", "alternating", "markov", "speeder"))
  )
})

# ---- Return schema ----------------------------------------------------------

test_that("cier_simulate returns the full cier_sim schema", {
  s <- cier_simulate(40L, it10, prevalence = 0.25, n_checks = 2L, seed = 7L)
  expect_s3_class(s, "cier_sim")
  expect_identical(names(s),
                   c("responses", "items", "seconds", "page_seconds",
                     "items_per_page", "checks", "pass", "truth"))
  expect_true(is.matrix(s$responses) && is.integer(s$responses))
  expect_identical(dim(s$responses), c(40L, 10L))
  expect_identical(s$items, it10)                          # verbatim passthrough
  expect_true(is.numeric(s$seconds) && length(s$seconds) == 40L)
  expect_true(all(s$seconds > 0))
  expect_true(is.matrix(s$page_seconds))
  expect_identical(dim(s$page_seconds), c(40L, 2L))        # one page per scale
  expect_identical(s$items_per_page, c(5L, 5L))
  expect_identical(dim(s$checks), c(40L, 2L))
  expect_length(s$pass, 2L)
  expect_identical(names(s$truth),
                   c("careless", "pattern", "extent", "onset_item",
                     "offset_item", "speeded", "params"))
  expect_identical(nrow(s$truth), 40L)
  expect_true(is.logical(s$truth$careless) && !anyNA(s$truth$careless))
  expect_true(is.logical(s$truth$speeded))
  expect_true(is.list(s$truth$params))
  # seconds is the row total of the page times (same generated cells).
  expect_equal(s$seconds, rowSums(s$page_seconds), tolerance = 1e-9)
})

test_that("the truth bookkeeping is internally consistent", {
  s <- cier_simulate(80L, it10, prevalence = 0.5, prop_partial = 0.5,
                     prop_temporary = 0.25, seed = 13L)
  truth <- s$truth
  expect_identical(sum(truth$careless), 40L)               # round(n * prevalence)
  expect_identical(truth$speeded, truth$careless)          # every careless row is
  # exact proportional extent allocation: 20 partial, 10 temporary, 10 full.
  expect_identical(sum(truth$extent == "partial"), 20L)
  expect_identical(sum(truth$extent == "temporary"), 10L)
  expect_identical(sum(truth$extent == "full"), 10L)
  # spans exist exactly on careless rows; attentive rows are all-clean labels.
  expect_true(all(is.na(truth$onset_item[!truth$careless])))
  expect_true(all(is.na(truth$offset_item[!truth$careless])))
  expect_true(all(!is.na(truth$onset_item[truth$careless])))
  expect_true(all(!is.na(truth$offset_item[truth$careless])))
  expect_true(all(truth$pattern[!truth$careless] == "attentive"))
  expect_true(all(truth$extent[!truth$careless] == "none"))
  # the default equal-weight mixture is APPLIED, not just defined: all eight patterns
  # appear among 40 careless draws at this seed (a single-pattern or truncated
  # default would fail; re-pick the seed deliberately if the draw order ever changes).
  expect_setequal(unique(truth$pattern[truth$careless]),
                  names(sim_default_patterns()))
  expect_true(all(vapply(truth$params[!truth$careless],
                         identical, logical(1L), list())))
})

test_that("the careless count is round(n * prevalence), not floor or ceiling", {
  s_up <- cier_simulate(30L, it10, prevalence = 0.29, seed = 1L)
  expect_identical(sum(s_up$truth$careless), 9L)           # round(8.7); floor = 8
  s_down <- cier_simulate(30L, it10, prevalence = 0.21, seed = 1L)
  expect_identical(sum(s_down$truth$careless), 6L)         # round(6.3); ceiling = 7
})

test_that("checks default to NULL and the generator metadata rides along", {
  s <- cier_simulate(20L, it10, seed = 3L,
                     timing = list(mu_car = log(1.2)))
  expect_null(s$checks)
  expect_null(s$pass)
  gen <- attr(s, "generator")
  expect_true(all(c("loadings", "effective_loadings", "factor_cor",
                    "thresholds", "residual_variance", "trait_distribution",
                    "trait_params", "timing", "p_fail_careless",
                    "p_fail_attentive", "seed") %in% names(gen)))
  expect_identical(gen$seed, 3L)
  expect_identical(gen$timing$mu_car, log(1.2))            # the RESOLVED knobs
  expect_identical(gen$timing$mu_att, log(8))
  s2 <- cier_simulate(20L, it10)
  expect_null(attr(s2, "generator")$seed)
})

test_that("validate_cier_sim returns its input and rejects a tampered object", {
  s <- cier_simulate(12L, it10, seed = 1L)
  expect_identical(validate_cier_sim(s), s)
  bad_truth <- s
  bad_truth$truth <- bad_truth$truth[-1L, ]
  expect_error(validate_cier_sim(bad_truth), class = "cier_error_input")
  bad_seconds <- s
  bad_seconds$seconds <- bad_seconds$seconds[-1L]
  expect_error(validate_cier_sim(bad_seconds), class = "cier_error_input")
  bad_pages <- s
  bad_pages$items_per_page <- c(4L, 5L)                    # no longer sums to p
  expect_error(validate_cier_sim(bad_pages), class = "cier_error_input")
})

# ---- Raw range invariant on a heterogeneous battery -------------------------

test_that("responses stay in each item's raw min..max range", {
  it_het <- sim_items_df(rep(c("E", "A"), each = 5L),
                         min = rep(c(0L, 1L), each = 5L),
                         max = rep(c(4L, 7L), each = 5L))
  s <- cier_simulate(60L, it_het, prevalence = 0.5, seed = 9L)
  for (j in seq_len(10L)) {
    expect_true(all(s$responses[, j] >= it_het$min[[j]]))
    expect_true(all(s$responses[, j] <= it_het$max[[j]]))
  }
})

# ---- RNG draw order: attentive-first contracts ------------------------------

test_that("a speeder-only run keeps the attentive content byte-identical", {
  s0 <- cier_simulate(60L, it10, prevalence = 0, seed = 7L)
  s1 <- cier_simulate(60L, it10, prevalence = 0.5, patterns = c(speeder = 1),
                      timing = list(sigma = 0, respondent_sd = 0), seed = 7L)
  # speeders are careless in TIME only: same seed, same attentive layer, no content
  # draw touches it -- the responses match cell for cell.
  expect_identical(s1$responses, s0$responses)
  careless <- s1$truth$careless
  # zero-noise times: speeder rows answer at exactly 1.5 s/item, attentive at 8.
  expect_equal(s1$seconds[careless], rep(10 * 1.5, sum(careless)),
               tolerance = 1e-9)
  expect_equal(s1$seconds[!careless], rep(10 * 8, sum(!careless)),
               tolerance = 1e-9)
  expect_true(all(s1$truth$speeded[careless]))
})

test_that("content patterns rewrite only the careless rows", {
  s0 <- cier_simulate(60L, it10, prevalence = 0, seed = 7L)
  s2 <- cier_simulate(60L, it10, prevalence = 0.5, patterns = c(random = 1),
                      seed = 7L)
  attentive_rows <- !s2$truth$careless
  expect_identical(s2$responses[attentive_rows, ], s0$responses[attentive_rows, ])
  expect_false(identical(s2$responses, s0$responses))
})

# ---- prop_speeded: slow-careless and the byte-identical default --------------

test_that("the default prop_speeded keeps speeded == careless (bytewise vs explicit 1)", {
  s_def <- cier_simulate(40L, it10, prevalence = 0.5, seed = 8L)
  s_one <- cier_simulate(40L, it10, prevalence = 0.5, prop_speeded = 1, seed = 8L)
  expect_identical(s_def, s_one)                             # default == explicit 1
  expect_identical(s_def$truth$speeded, s_def$truth$careless)
})

test_that("prop_speeded < 1 plants slow-careless rows (careless content, attentive pace)", {
  # Regression: today every careless row speeds, so a careless-content row at normal
  # pace is ungeneratable. With prop_speeded = 0 and zero-noise timing, straightline
  # rows keep their planted content but answer at the attentive 8 s/item -- and
  # NOBODY speeds.
  s <- cier_simulate(40L, it10, prevalence = 0.5, patterns = c(straightline = 1),
                     prop_speeded = 0,
                     pattern_params = list(straightline = list(anchor = "value",
                                                               value = 5L)),
                     timing = list(sigma = 0, respondent_sd = 0), seed = 11L)
  careless <- s$truth$careless
  expect_true(all(s$responses[careless, ] == 5L))           # content IS careless
  expect_false(any(s$truth$speeded))                        # ... but nobody speeds
  expect_equal(s$seconds[careless], rep(10 * 8, sum(careless)), tolerance = 1e-9)
  expect_false(identical(s$truth$speeded, s$truth$careless))
})

test_that("the speeded flag and the careless times are the SAME rows (no desync)", {
  # an intermediate prop_speeded with zero-noise timing: the rows marked
  # truth$speeded must be EXACTLY the rows shifted to careless pace. A mutant that
  # computes speeded from one permutation but time-shifts a different content subset
  # would split the two sets and fail here.
  s <- cier_simulate(80L, it10, prevalence = 0.5, patterns = c(straightline = 1),
                     prop_speeded = 0.5,
                     timing = list(sigma = 0, respondent_sd = 0), seed = 21L)
  careless <- s$truth$careless
  speeded  <- s$truth$speeded
  expect_identical(sum(speeded), as.integer(round(sum(careless) * 0.5)))  # half
  expect_equal(s$seconds[speeded], rep(10 * 1.5, sum(speeded)), tolerance = 1e-9)
  expect_equal(s$seconds[careless & !speeded],
               rep(10 * 8, sum(careless & !speeded)), tolerance = 1e-9)
})

test_that("prop_speeded splits content rows by count; speeders always speed", {
  s <- cier_simulate(80L, it10, prevalence = 0.5,
                     patterns = c(straightline = 0.5, speeder = 0.5),
                     prop_speeded = 0, seed = 12L)
  is_speeder <- s$truth$pattern == "speeder"
  is_content <- s$truth$careless & !is_speeder
  expect_true(all(s$truth$speeded[is_speeder]))             # speeders speed regardless
  expect_false(any(s$truth$speeded[is_content]))            # content rows do not (p = 0)
  expect_false(any(s$truth$speeded[!s$truth$careless]))     # attentive never speeds
})

test_that("cier_simulate rejects an out-of-range prop_speeded", {
  expect_error(cier_simulate(10L, it10, prop_speeded = 1.2), class = "cier_error_input")
  expect_error(cier_simulate(10L, it10, prop_speeded = -0.1), class = "cier_error_input")
})

# ---- Provenance: truth$params is what the engine applied --------------------

test_that("recorded pattern_params drive the planted content", {
  pp <- list(straightline = list(anchor = "value", value = 5L))
  s <- cier_simulate(40L, it10, prevalence = 0.5,
                     patterns = c(straightline = 1), pattern_params = pp,
                     seed = 4L)
  careless <- s$truth$careless
  expect_true(all(s$responses[careless, ] == 5L))          # the knob is applied
  expect_true(all(vapply(s$truth$params[careless],        # ... and recorded
                         identical, logical(1L), pp$straightline)))
  expect_true(all(vapply(s$truth$params[!careless],
                         identical, logical(1L), list())))
})

# ---- Partial rows: the time shift starts at the onset ------------------------

test_that("partial rows speed up exactly from the onset page on", {
  it4 <- sim_items_df(rep(c("E", "A", "C", "N"), each = 5L))
  s <- cier_simulate(40L, it4, prevalence = 0.5, patterns = c(random = 1),
                     prop_partial = 1, onset_window = c(11L, 11L),
                     timing = list(sigma = 0, respondent_sd = 0), seed = 11L)
  expect_identical(s$items_per_page, rep(5L, 4L))          # one page per scale
  careless <- s$truth$careless
  expect_true(all(s$truth$onset_item[careless] == 11L))
  expect_true(all(s$truth$extent[careless] == "partial"))
  n_car <- sum(careless)
  # pages 1-2 (items 1-10) are PRE-onset: attentive pace 5 * 8 = 40 s.
  expect_equal(unname(s$page_seconds[careless, 1:2]),
               matrix(40, n_car, 2L), tolerance = 1e-9)
  # pages 3-4 (items 11-20) are in the careless span: 5 * 1.5 = 7.5 s.
  expect_equal(unname(s$page_seconds[careless, 3:4]),
               matrix(7.5, n_car, 2L), tolerance = 1e-9)
  expect_equal(unname(s$page_seconds[!careless, ]),
               matrix(40, sum(!careless), 4L), tolerance = 1e-9)
})

test_that("an explicit items_per_page drives the aggregation, not just validation", {
  # zero-noise, zero-careless: every cell is exactly 8 s, so each page total is
  # items_per_page * 8 -- a resolver that validates the argument but aggregates with
  # the scale-run default produces the wrong page count / totals.
  s_vec <- cier_simulate(6L, it10, prevalence = 0, items_per_page = c(4L, 6L),
                         timing = list(sigma = 0, respondent_sd = 0), seed = 5L)
  expect_identical(s_vec$items_per_page, c(4L, 6L))
  expect_identical(dim(s_vec$page_seconds), c(6L, 2L))
  expect_equal(unname(s_vec$page_seconds),
               matrix(rep(c(4, 6) * 8, each = 6L), 6L, 2L), tolerance = 1e-9)
  s_scl <- cier_simulate(6L, it10, prevalence = 0, items_per_page = 3L,
                         timing = list(sigma = 0, respondent_sd = 0), seed = 5L)
  expect_identical(s_scl$items_per_page, c(3L, 3L, 3L, 1L))  # remainder tail
  expect_equal(unname(s_scl$page_seconds),
               matrix(rep(c(3, 3, 3, 1) * 8, each = 6L), 6L, 4L),
               tolerance = 1e-9)
})

# ---- Attention checks end-to-end --------------------------------------------

test_that("requested checks thread p_fail through to cier_attention exactly", {
  s <- cier_simulate(40L, it10, prevalence = 0.5, n_checks = 3L,
                     p_fail_careless = 1, p_fail_attentive = 0, seed = 2L)
  value <- cier_attention(s$checks, s$pass)$value
  careless <- s$truth$careless
  expect_true(all(value[careless] == 3))
  expect_true(all(value[!careless] == 0))
  expect_identical(attr(s, "generator")$p_fail_careless, 1)
  expect_identical(attr(s, "generator")$p_fail_attentive, 0)
})

# ---- Seed contract ----------------------------------------------------------

test_that("the same seed reproduces the whole object bytewise", {
  a <- cier_simulate(30L, it10, prevalence = 0.3, prop_partial = 0.5,
                     n_checks = 2L, seed = 42L)
  b <- cier_simulate(30L, it10, prevalence = 0.3, prop_partial = 0.5,
                     n_checks = 2L, seed = 42L)
  expect_identical(a, b)
  d <- cier_simulate(30L, it10, prevalence = 0.3, prop_partial = 0.5,
                     n_checks = 2L, seed = 43L)
  expect_false(identical(a$responses, d$responses))
})

test_that("a seeded call leaves the ambient RNG stream untouched", {
  set.seed(99L)
  x1 <- stats::runif(1L)
  set.seed(99L)
  invisible(cier_simulate(10L, it10, seed = 1L))
  x2 <- stats::runif(1L)
  expect_identical(x1, x2)
})

test_that("a NULL seed draws from the ambient stream", {
  set.seed(5L)
  d1 <- cier_simulate(10L, it10)
  set.seed(5L)
  d2 <- cier_simulate(10L, it10)
  expect_identical(d1$responses, d2$responses)             # ambient-reproducible
  d3 <- cier_simulate(10L, it10)                           # stream has advanced
  expect_false(identical(d1$responses, d3$responses))
})

test_that("the frozen-seed RNG-stream digest is stable", {
  # Any new, removed, or reordered draw inside cier_simulate() changes this digest
  # for the recorded seed. That is sometimes intended (a deliberate generator
  # change) -- re-record the hash and say so in the hand-off; it must never drift
  # silently.
  s <- cier_simulate(12L, sim_items_df(rep(c("E", "A"), each = 3L)),
                     prevalence = 0.5, prop_partial = 0.5,
                     prop_temporary = 0.25, n_checks = 2L,
                     seed = 20260612L)
  expect_identical(sim_digest(s), "a911927b19826e9ede79880757a022e3")
})

# ---- Marginals / trait threading --------------------------------------------

test_that("marginals and trait_distribution thread through to the GRM layer", {
  base <- cier_simulate(25L, it10, seed = 6L)
  s <- cier_simulate(25L, it10, marginals = "skewed_right", seed = 6L)
  expect_true(all(s$responses >= 1L & s$responses <= 5L))
  # APPLIED, not just recorded: at the same seed a different marginal shape / trait
  # distribution must move the attentive responses (a knob that only lands in the
  # metadata leaves them identical).
  expect_false(identical(s$responses, base$responses))
  s2 <- cier_simulate(25L, it10, trait_distribution = "bimodal", seed = 6L)
  expect_true(all(s2$responses >= 1L & s2$responses <= 5L))
  expect_false(identical(s2$responses, base$responses))
  expect_identical(attr(s2, "generator")$trait_distribution, "bimodal")
  expect_error(
    cier_simulate(25L, it10, factor_cor = matrix(c(1, 2, 2, 1), 2L), seed = 1L),
    class = "cier_error_input"
  )
})

test_that("the orchestrator composes the kernels in the documented draw order", {
  # White-box stage-order pin (the RPR-style coordination): re-run the four
  # already-oracle-tested kernels in the documented order -- attentive -> plan ->
  # content -> times -> checks -- under one seed and demand the public object
  # bytewise. A reordered, split, or duplicated draw in ANY stage diverges; the
  # frozen digest below then catches cross-version drift on top.
  it6 <- sim_items_df(rep(c("E", "A"), each = 3L))
  s <- cier_simulate(14L, it6, prevalence = 0.5, prop_partial = 0.5,
                     n_checks = 2L, seed = 314L)
  ref <- withr::with_seed(314L, {
    it <- check_items_simulate(it6, 6L)
    att <- sim_attentive_with_metadata(14L, it6)
    truth <- sim_build_plan(14L, 6L, 0.5, sim_default_patterns(),
                            prop_partial = 0.5)
    resp <- sim_apply_patterns(att$responses, it, truth)
    cells <- sim_times_lognormal(14L, 6L, truth$careless, truth$onset_item,
                                 truth$offset_item,
                                 timing = sim_timing_defaults())
    checks <- sim_direct_checks(14L, truth$careless, 2L)
    list(resp = resp, cells = cells, checks = checks, truth = truth)
  })
  expect_identical(s$responses, ref$resp)
  expect_identical(s$truth$careless, ref$truth$careless)
  expect_identical(s$truth$pattern, ref$truth$pattern)
  ipp <- sim_resolve_pages(NULL, it6$scale)
  expect_equal(s$page_seconds, sim_page_totals(ref$cells, ipp),
               tolerance = 1e-12)
  expect_equal(s$seconds, rowSums(ref$cells), tolerance = 1e-9)
  expect_identical(s$checks, ref$checks$checks)
  expect_identical(s$pass, ref$checks$pass)
})

# ---- Input validation --------------------------------------------------------

test_that("cier_simulate rejects malformed inputs with typed errors", {
  expect_error(cier_simulate(0L, it10), class = "cier_error_input")
  expect_error(cier_simulate(10L, data.frame(max = 5)),    # no scale column
               class = "cier_error_input")
  expect_error(cier_simulate(10L, data.frame(scale = "E")),  # no max column
               class = "cier_error_input")
  expect_error(cier_simulate(10L, it10, prevalence = 1.5),
               class = "cier_error_input")
  expect_error(cier_simulate(10L, it10, prevalence = -0.1),
               class = "cier_error_input")
  expect_error(cier_simulate(10L, it10, patterns = c(bogus = 1)),
               class = "cier_error_input")
  expect_error(cier_simulate(10L, it10, patterns = c(random = 0.4)),
               class = "cier_error_input")                 # weights not sum 1
  expect_error(cier_simulate(10L, it10, items_per_page = c(3L, 3L)),
               class = "cier_error_input")                 # does not sum to p
  expect_error(cier_simulate(10L, it10, items_per_page = 25L),
               class = "cier_error_input")
  expect_error(cier_simulate(10L, it10, timing = list(bogus = 1)),
               class = "cier_error_input")
  expect_error(cier_simulate(10L, it10, n_checks = -1L),
               class = "cier_error_input")
  expect_error(cier_simulate(10L, it10, n_checks = 1.5),
               class = "cier_error_input")
  expect_error(cier_simulate(10L, it10, p_fail_careless = 1.2),
               class = "cier_error_input")
  expect_error(cier_simulate(10L, it10, seed = 1.5), class = "cier_error_input")
  expect_error(cier_simulate(10L, it10, seed = "a"), class = "cier_error_input")
  expect_error(cier_simulate(10L, sim_items_df(c("E", "E")), prop_temporary = 1),
               class = "cier_error_input")                 # temporary needs p >= 3
})

test_that("a malformed cheap scalar aborts up front, before any RNG draw", {
  # prop_partial / prop_temporary (and their sum), onset_window, n_checks, and
  # p_fail_* are validated at the boundary, BEFORE the GRM + plan + content + timing
  # generation. With a NULL seed the generator draws from the ambient stream, so a
  # knob aborting only inside a later stage (sim_build_plan / sim_direct_checks)
  # would already have advanced .Random.seed by the time it errored; validated up
  # front it leaves the stream byte-identical. (A non-NULL seed save/restores the
  # stream regardless, so this is the only seed under which the contract is
  # observable.) A mutant that left a guard in its original stage fails only here.
  bad <- list(
    list(prop_partial = 1.5),
    list(prop_temporary = -0.1),
    list(prop_partial = 0.7, prop_temporary = 0.7),    # the <= 1 sum guard
    list(onset_window = c(50L, 60L)),                  # above [1, p]
    list(onset_window = c(0L, 5L)),                    # below the lower bound
    list(onset_window = c(8L, 3L)),                    # lo > hi
    list(prop_temporary = 0.5, onset_window = c(9L, 10L)),  # no recovery column
    list(n_checks = -1L),
    list(n_checks = 1.5),
    list(p_fail_careless = 1.2),
    list(p_fail_attentive = -0.1)
  )
  for (knob in bad) {
    set.seed(404L)
    x1 <- stats::runif(1L)
    set.seed(404L)
    expect_error(do.call(cier_simulate, c(list(60L, it10), knob)),
                 class = "cier_error_input")
    expect_identical(stats::runif(1L), x1)             # the failed call drew nothing
  }
})

test_that("n_checks = 0 is accepted (non-negative, not the positive-only count)", {
  # The boundary hoists n_checks with checkmate::check_count() (non-negative),
  # mirroring sim_direct_checks -- NOT the positive-only check_count() wrapper -- so
  # the default n_checks = 0 must still run rather than abort at the boundary.
  expect_no_error(cier_simulate(10L, it10, n_checks = 0L, seed = 1L))
})

test_that("the hoisted guards keep the inner kernels' verbatim messages", {
  # The boundary copies each guard's wording verbatim from its original stage, so
  # moving the check up front leaves the user-visible diagnostic byte-identical. Pin
  # the three most likely to drift (regexp matches the message header). n_checks is
  # the sharpest: its non-negative wording differs from the positive-only
  # check_count() wrapper, the exact mutant the boundary must avoid.
  expect_error(cier_simulate(60L, it10, n_checks = -1L),
               class = "cier_error_input",
               regexp = "non-negative whole number")
  expect_error(cier_simulate(60L, it10, prop_partial = 0.7, prop_temporary = 0.7),
               class = "cier_error_input", regexp = "must not exceed 1")
  expect_error(cier_simulate(60L, it10, onset_window = c(50L, 60L)),
               class = "cier_error_input", regexp = "two whole numbers")
})

test_that("pattern_params must name known patterns and be lists", {
  expect_error(
    cier_simulate(10L, it10, pattern_params = list(bogus = list())),
    class = "cier_error_input"
  )
  expect_error(
    cier_simulate(10L, it10, patterns = c(random = 1),
                  pattern_params = list(straightline = list(value = 3L))),
    class = "cier_error_input"                             # not in `patterns`
  )
  expect_error(
    cier_simulate(10L, it10, patterns = c(random = 1),
                  pattern_params = list(random = "x")),    # entry not a list
    class = "cier_error_input"
  )
  expect_error(
    cier_simulate(10L, it10, pattern_params = list(list())),  # unnamed entry
    class = "cier_error_input"
  )
})

test_that("pattern_params rejects unknown per-pattern knobs", {
  # A typo'd knob on a real pattern is silently ignored today, yet recorded in
  # truth$params as applied; now a typed error.
  expect_error(
    cier_simulate(10L, it10, patterns = c(straightline = 1),
                  pattern_params = list(straightline = list(ancho = "value"))),
    class = "cier_error_input"
  )
  # a knob belonging to a different pattern (p_high is extreme's, not midpoint's).
  expect_error(
    cier_simulate(10L, it10, patterns = c(midpoint = 1),
                  pattern_params = list(midpoint = list(p_high = 0.5))),
    class = "cier_error_input"
  )
  # a knob real for ANOTHER pattern (period is alternating's) is still unknown for
  # straightline -- distinguishes a per-pattern allowlist from a global union of all
  # knob names.
  expect_error(
    cier_simulate(10L, it10, patterns = c(straightline = 1),
                  pattern_params = list(straightline = list(period = 3L))),
    class = "cier_error_input"
  )
  # random / speeder accept NO knobs.
  expect_error(
    cier_simulate(10L, it10, patterns = c(random = 1),
                  pattern_params = list(random = list(anything = 1))),
    class = "cier_error_input"
  )
  expect_error(
    cier_simulate(10L, it10, patterns = c(speeder = 1),
                  pattern_params = list(speeder = list(p_high = 1))),
    class = "cier_error_input"
  )
  # an unnamed knob entry cannot be matched to a knob name.
  expect_error(
    cier_simulate(10L, it10, patterns = c(straightline = 1),
                  pattern_params = list(straightline = list(5))),
    class = "cier_error_input"
  )
})

test_that("every documented per-pattern knob is accepted (allowlist not too narrow)", {
  # one valid value per knob name, all six knob-bearing patterns in one call: a
  # too-narrow allowlist (a real knob missing) raises "Unknown knob" here.
  pp <- list(
    straightline = list(anchor = "position", q = 0.3, value = 5L, switch_prob = 0.2),
    midpoint     = list(jitter = 1L),
    extreme      = list(p_high = 0.7),
    diagonal     = list(step = 2L, start = 1L, bounce = TRUE),
    alternating  = list(period = 3L, values = c(1L, 3L, 5L), start_offset = 1L),
    markov       = list(transition = matrix(0.2, 5L, 5L), initial = rep(0.2, 5L))
  )
  w <- stats::setNames(rep(1 / 6, 6L),
                       c("straightline", "midpoint", "extreme", "diagonal",
                         "alternating", "markov"))
  expect_no_error(cier_simulate(60L, it10, prevalence = 0.5, patterns = w,
                                pattern_params = pp, seed = 1L))
})

# ---- Recovery smokes: each pattern caught by its matched index --------------

test_that("each planted pattern is recovered by its matched index end-to-end", {
  rec <- function(seed, pattern, params = list()) {
    cier_simulate(200L, it20, prevalence = 0.25,
                  patterns = stats::setNames(1, pattern),
                  pattern_params = params, seed = seed)
  }
  expect_gt(local({
    s <- rec(11L, "straightline",
             params = list(straightline = list(anchor = "value")))
    ref_rank_auc(cier_longstring(s$responses)$value, s$truth$careless)
  }), 0.85)
  expect_gt(local({
    s <- rec(12L, "diagonal")
    # Simulated fixtures can trip the percentile-cutoff degeneracy guard (small-n /
    # saturation); these recovery tests assert rank recovery, not flags, so the
    # (correct) warning is muffled.
    ref_rank_auc(suppressWarnings(cier_lazr(s$responses))$value, s$truth$careless)
  }), 0.85)
  expect_gt(local({
    # supplementary coarse floor; the deterministic honoured-transition test in
    # test-sim-patterns.R is the real markov guard. Do not tighten.
    s <- rec(13L, "markov")
    ref_rank_auc(cier_lazr(s$responses)$value, s$truth$careless)
  }), 0.75)
  expect_gt(local({
    s <- rec(14L, "alternating")
    ref_rank_auc(suppressWarnings(cier_autocorrelation(s$responses, max_lag = 8L))$value,
                 s$truth$careless)
  }), 0.85)
  expect_gt(local({
    s <- rec(15L, "random")
    ref_rank_auc(cier_mahalanobis(s$responses)$value, s$truth$careless)
  }), 0.75)
})

test_that("speeders are caught by timing and missed by content indices", {
  s <- cier_simulate(200L, it20, prevalence = 0.25, patterns = c(speeder = 1),
                     seed = 16L)
  careless <- s$truth$careless
  # timing must catch: low totals rank the speeders high.
  expect_gt(ref_rank_auc(-cier_total_time(s$seconds)$value, careless), 0.70)
  # content must miss: speeder responses are attentive, so a content index ranks
  # them no better than chance (band, not a point).
  expect_lt(ref_rank_auc(cier_longstring(s$responses)$value, careless), 0.65)
})

# ---- Print (snapshots lock the approved mock-up) -----------------------------

test_that("print renders the locked summary (checks present)", {
  it12 <- sim_items_df(rep(c("E", "A", "C"), each = 4L))
  s <- cier_simulate(200L, it12, prevalence = 0.2, prop_partial = 0.2,
                     prop_temporary = 0.1, n_checks = 2L, seed = 2026L)
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    expect_snapshot(print(s))
  })
})

test_that("print renders the no-checks, all-full variant", {
  s <- cier_simulate(50L, it10, prevalence = 0.2, seed = 2027L)
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    expect_snapshot(print(s))
  })
})

test_that("print renders the zero-careless variant", {
  s <- cier_simulate(30L, it10, prevalence = 0, seed = 2028L)
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    expect_snapshot(print(s))
  })
})

# ---- summary delegates to print ----------------------------------------------

test_that("summary.cier_sim prints the same body as print and returns invisibly", {
  s <- cier_simulate(50L, it10, prevalence = 0.2, seed = 2027L)
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    p <- capture.output(print(s))
    out <- capture.output(summary(s))
  })
  expect_identical(out, p)
  expect_identical(withVisible(summary(s))$visible, FALSE)
})
