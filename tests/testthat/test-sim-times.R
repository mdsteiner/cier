# Tests for the response-time generator + page aggregation (Slice 26) -- the
# timing half of cier_simulate(). The kernel couples the careless lognormal
# shift to the planted careless span, so partial / temporary rows speed up
# only from their onset (the within-survey burst cier_page_time() is built to
# see).
#
# Trust model: oracle-only (a generator, not an index, with no CRAN partner).
# The independent oracle (ref-sim-times.R) re-derives the cell times by hand
# with a coordinated draw order, and the calibration of the signed-off timing
# defaults (mu_att = log(8), mu_car = log(1.5), sigma = 0.5, respondent_sd =
# 1.2) is HELD by the slow-tier acceptance test below, not asserted: default
# careless rows must be flaggable by the cited 2 s/item page rule AND the
# total-time rank-AUC must land inside the published .66-.92 band, never 1.0.
# See tests/reference/TOLERANCES.md.

source(test_path("..", "reference", "ref-sim-times.R"))
source(test_path("..", "reference", "ref-sim-patterns.R"))

# =============================================================================
# FAST TIER
# =============================================================================

# ---- Timing defaults + resolver ---------------------------------------------

test_that("sim_timing_defaults pins the signed-off calibration", {
  expect_identical(
    sim_timing_defaults(),
    list(mu_att = log(8), mu_car = log(1.5), sigma = 0.5,
         respondent_sd = 1.2, min_seconds = 0.1)
  )
})

test_that("sim_resolve_timing merges partial overrides over the defaults", {
  expect_identical(sim_resolve_timing(list()), sim_timing_defaults())
  tm <- sim_resolve_timing(list(mu_car = log(2), sigma = 0))
  expect_identical(tm$mu_car, log(2))
  expect_identical(tm$sigma, 0)
  expect_identical(tm$mu_att, log(8))          # untouched defaults survive
  expect_identical(tm$respondent_sd, 1.2)
  expect_identical(tm$min_seconds, 0.1)
  expect_identical(names(tm), names(sim_timing_defaults()))
})

test_that("sim_resolve_timing rejects unknown keys and out-of-range values", {
  expect_error(sim_resolve_timing(list(bogus = 1)), class = "cier_error_input")
  expect_error(sim_resolve_timing(list(mu_atte = log(8))),    # typo must not
               class = "cier_error_input")                    # silently default
  expect_error(sim_resolve_timing("fast"), class = "cier_error_input")
  expect_error(sim_resolve_timing(list(sigma = -0.1)), class = "cier_error_input")
  expect_error(sim_resolve_timing(list(respondent_sd = -1)),
               class = "cier_error_input")
  expect_error(sim_resolve_timing(list(min_seconds = 0)),     # must be > 0 so
               class = "cier_error_input")                    # check_seconds passes
  expect_error(sim_resolve_timing(list(min_seconds = -1)), class = "cier_error_input")
  expect_error(sim_resolve_timing(list(mu_att = Inf)), class = "cier_error_input")
  expect_error(sim_resolve_timing(list(mu_car = NA_real_)), class = "cier_error_input")
  expect_error(sim_resolve_timing(list(sigma = c(0.5, 0.5))),
               class = "cier_error_input")
  expect_error(sim_resolve_timing(list(sigma = "fast")), class = "cier_error_input")
})

# ---- Kernel: oracle parity (coordinated draw order) -------------------------

test_that("the cell times match the independent hand re-derivation", {
  n <- 7L
  p <- 9L
  careless <- c(FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE)
  onset <- c(NA, 1L, NA, NA, 4L, NA, 3L)     # full / partial / temporary spans
  offset <- c(NA, 9L, NA, NA, 9L, NA, 6L)
  tm <- sim_timing_defaults()
  prod <- withr::with_seed(
    20260612L,
    sim_times_lognormal(n, p, careless, onset, offset, timing = tm)
  )
  ref <- ref_times_lognormal(20260612L, n, p, careless, onset, offset,
                             mu_att = tm$mu_att, mu_car = tm$mu_car,
                             sigma = tm$sigma, respondent_sd = tm$respondent_sd,
                             min_seconds = tm$min_seconds)
  expect_equal(prod, ref, tolerance = 1e-12)
  # the kernel genuinely consumes RNG in the documented order: a different seed
  # must move every cell.
  prod2 <- withr::with_seed(
    20260613L,
    sim_times_lognormal(n, p, careless, onset, offset, timing = tm)
  )
  expect_false(identical(prod, prod2))
})

test_that("oracle parity holds under overridden timing parameters", {
  n <- 5L
  p <- 6L
  careless <- c(TRUE, FALSE, TRUE, FALSE, FALSE)
  onset <- c(2L, NA, 1L, NA, NA)
  offset <- c(5L, NA, 6L, NA, NA)
  tm <- sim_resolve_timing(list(mu_att = log(6), mu_car = log(1),
                                sigma = 0.3, respondent_sd = 0.7,
                                min_seconds = 0.5))
  prod <- withr::with_seed(
    77L, sim_times_lognormal(n, p, careless, onset, offset, timing = tm)
  )
  ref <- ref_times_lognormal(77L, n, p, careless, onset, offset,
                             mu_att = log(6), mu_car = log(1), sigma = 0.3,
                             respondent_sd = 0.7, min_seconds = 0.5)
  expect_equal(prod, ref, tolerance = 1e-12)
  # LOAD-BEARING: the floor must bind on some cells here (with a non-zero pace
  # intercept), or this fixture stops distinguishing floor-after-pace from the
  # floor-before-pace mutant, whose difference is sub-tolerance when no cell is
  # floored. Re-pick the seed deliberately if this ever fails.
  expect_true(any(prod == 0.5))
})

# ---- Kernel: zero-noise closed form (span off-by-ones pinned exactly) -------

test_that("with zero noise each cell is exactly exp(mu) on its own span side", {
  n <- 4L
  p <- 10L
  careless <- c(FALSE, TRUE, TRUE, TRUE)
  onset <- c(NA, 1L, 4L, 3L)                 # full [1,10], partial [4,10],
  offset <- c(NA, 10L, 10L, 6L)              # temporary [3,6]
  tm <- sim_resolve_timing(list(sigma = 0, respondent_sd = 0))
  times <- sim_times_lognormal(n, p, careless, onset, offset, timing = tm)
  # attentive row: 8 s/item everywhere.
  expect_equal(times[1L, ], rep(8, p), tolerance = 1e-12)
  # full careless row: 1.5 s/item everywhere.
  expect_equal(times[2L, ], rep(1.5, p), tolerance = 1e-12)
  # partial row [4, 10]: the attentive prefix 1:3 stays 8; the span INCLUDES
  # the onset column (an off-by-one mutant shifts this boundary).
  expect_equal(times[3L, 1:3], rep(8, 3L), tolerance = 1e-12)
  expect_equal(times[3L, 4:10], rep(1.5, 7L), tolerance = 1e-12)
  # temporary row [3, 6]: prefix AND recovery tail attentive; both boundary
  # columns (3 and 6) are inside the careless span.
  expect_equal(times[4L, c(1:2, 7:10)], rep(8, 6L), tolerance = 1e-12)
  expect_equal(times[4L, 3:6], rep(1.5, 4L), tolerance = 1e-12)
})

test_that("the min_seconds floor binds when exp(mu) falls below it", {
  tm <- sim_resolve_timing(list(sigma = 0, respondent_sd = 0,
                                mu_car = log(0.05), min_seconds = 0.1))
  times <- sim_times_lognormal(2L, 3L, c(TRUE, FALSE), c(1L, NA), c(3L, NA),
                               timing = tm)
  expect_identical(times[1L, ], rep(0.1, 3L))            # floored exactly
  expect_equal(times[2L, ], rep(8, 3L), tolerance = 1e-12)
})

# ---- Kernel: input validation -----------------------------------------------

test_that("the kernel rejects malformed careless / onset / offset inputs", {
  tm <- sim_timing_defaults()
  expect_error(                                           # careless not logical
    sim_times_lognormal(3L, 4L, c(1, 0, 1), c(1L, NA, 1L), c(4L, NA, 4L), tm),
    class = "cier_error_input"
  )
  expect_error(                                           # wrong length
    sim_times_lognormal(3L, 4L, c(TRUE, FALSE), c(1L, NA), c(4L, NA), tm),
    class = "cier_error_input"
  )
  expect_error(                                           # NA onset on careless
    sim_times_lognormal(2L, 4L, c(TRUE, FALSE), c(NA, NA), c(4L, NA), tm),
    class = "cier_error_input"
  )
  expect_error(                                           # onset > offset
    sim_times_lognormal(2L, 4L, c(TRUE, FALSE), c(3L, NA), c(2L, NA), tm),
    class = "cier_error_input"
  )
  expect_error(                                           # onset < 1
    sim_times_lognormal(2L, 4L, c(TRUE, FALSE), c(0L, NA), c(4L, NA), tm),
    class = "cier_error_input"
  )
  expect_error(                                           # offset > p
    sim_times_lognormal(2L, 4L, c(TRUE, FALSE), c(1L, NA), c(5L, NA), tm),
    class = "cier_error_input"
  )
})

# ---- Pages: resolver ---------------------------------------------------------

test_that("NULL items_per_page resolves to one page per contiguous scale run", {
  expect_identical(sim_resolve_pages(NULL, rep(c("E", "A", "C"), each = 4L)),
                   rep(4L, 3L))
  expect_identical(sim_resolve_pages(NULL, c("E", "E", "A", "A", "A")),
                   c(2L, 3L))
  # non-contiguous scale labels: one page per RUN, not per unique label.
  expect_identical(sim_resolve_pages(NULL, c("E", "E", "A", "A", "E")),
                   c(2L, 2L, 1L))
  expect_identical(sim_resolve_pages(NULL, rep("F", 6L)), 6L)
})

test_that("a scalar items_per_page chunks uniform pages with a remainder tail", {
  scale <- rep("F", 10L)
  expect_identical(sim_resolve_pages(4L, scale), c(4L, 4L, 2L))
  expect_identical(sim_resolve_pages(5L, scale), c(5L, 5L))
  expect_identical(sim_resolve_pages(10L, scale), 10L)    # one page == p
  expect_identical(sim_resolve_pages(1L, scale), rep(1L, 10L))
  expect_identical(sim_resolve_pages(3, scale), c(3L, 3L, 3L, 1L))  # numeric ok
})

test_that("an explicit items_per_page vector must be whole, positive, sum p", {
  scale <- rep("F", 8L)
  expect_identical(sim_resolve_pages(c(3L, 5L), scale), c(3L, 5L))
  expect_identical(sim_resolve_pages(c(3, 4, 1), scale), c(3L, 4L, 1L))
  expect_error(sim_resolve_pages(c(3L, 4L), scale), class = "cier_error_input")
  expect_error(sim_resolve_pages(c(0L, 8L), scale), class = "cier_error_input")
  expect_error(sim_resolve_pages(c(-2L, 10L), scale), class = "cier_error_input")
  expect_error(sim_resolve_pages(c(2.5, 5.5), scale), class = "cier_error_input")
  expect_error(sim_resolve_pages(c(4L, NA), scale), class = "cier_error_input")
  expect_error(sim_resolve_pages(9L, scale), class = "cier_error_input")  # > p
  expect_error(sim_resolve_pages(0L, scale), class = "cier_error_input")
  expect_error(sim_resolve_pages(2.5, scale), class = "cier_error_input")
  expect_error(sim_resolve_pages("4", scale), class = "cier_error_input")
})

# ---- Pages: totals vs the hand loop -----------------------------------------

test_that("page totals match the independent double-loop re-derivation", {
  cell <- withr::with_seed(3L, matrix(stats::runif(5L * 9L, 1, 20), 5L, 9L))
  ipp <- c(2L, 4L, 3L)                                    # uneven pages
  expect_equal(sim_page_totals(cell, ipp), ref_page_totals(cell, ipp),
               tolerance = 1e-12)
  expect_equal(sim_page_totals(cell, 9L), ref_page_totals(cell, 9L),
               tolerance = 1e-12)                         # single page
  one_row <- cell[1L, , drop = FALSE]                     # n = 1 keeps shape
  expect_identical(dim(sim_page_totals(one_row, ipp)), c(1L, 3L))
  expect_equal(sim_page_totals(one_row, ipp), ref_page_totals(one_row, ipp),
               tolerance = 1e-12)
})

test_that("page totals partition the row total exactly", {
  cell <- withr::with_seed(4L, matrix(stats::runif(6L * 7L, 0.5, 30), 6L, 7L))
  ipp <- c(3L, 3L, 1L)
  expect_equal(rowSums(sim_page_totals(cell, ipp)), rowSums(cell),
               tolerance = 1e-9)
})

# =============================================================================
# SLOW TIER -- the timing-calibration acceptance test (signed off 2026-06-12)
# =============================================================================

test_that("default timing is page-flaggable and lands in the published AUC band", {
  skip_on_cran()
  skip_if_slow()
  # End-to-end at the signed-off defaults: 4,000 respondents, half fully
  # careless, 20 five-point items on 4 scales (so the default one-page-per-scale
  # rule gives 4 five-item pages). Acceptance (recorded in TOLERANCES.md):
  #   (a) total-time rank-AUC inside the published [.66, .92] band -- never 1.0
  #       (overlap is by construction; a perfectly separable simulator is wrong
  #       by spec);
  #   (b) the cited Bowling et al. (2023) 2 s/item page rule -- the shipped
  #       cier_page_time() default -- flags at least half the careless rows
  #       (the archive's mu_car = log(4) flagged essentially none);
  #   (c) attentive rows keep a real but bounded false-flag rate (<= 0.25).
  items <- data.frame(scale = rep(c("E", "A", "C", "N"), each = 5L), max = 5)
  sim <- cier_simulate(4000L, items, prevalence = 0.5,
                       patterns = c(random = 1), seed = 20260612L)
  careless <- sim$truth$careless
  total <- cier_total_time(sim$seconds)
  auc <- ref_rank_auc(-total$value, careless)             # lower totals flag
  expect_gte(auc, 0.66)
  expect_lte(auc, 0.92)
  page <- cier_page_time(sim$page_seconds, sim$items_per_page)
  sensitivity <- mean(page$flagged[careless])
  false_flag <- mean(page$flagged[!careless])
  expect_gte(sensitivity, 0.5)
  expect_lte(false_flag, 0.25)
  expect_gt(sensitivity, false_flag)                      # the rule is informative
})

test_that("partial careless rows are slower than full ones in total time", {
  skip_on_cran()
  skip_if_slow()
  # The shift applies only from the onset on, so a partial row's expected total
  # sits between the attentive and the fully-careless totals -- the coupling
  # cier_page_time's burst detection relies on, aggregated.
  items <- data.frame(scale = rep(c("E", "A"), each = 10L), max = 5)
  sim_full <- cier_simulate(1500L, items, prevalence = 0.5,
                            patterns = c(random = 1), seed = 31L)
  sim_part <- cier_simulate(1500L, items, prevalence = 0.5,
                            patterns = c(random = 1), prop_partial = 1,
                            seed = 31L)
  full_med <- stats::median(sim_full$seconds[sim_full$truth$careless])
  part_med <- stats::median(sim_part$seconds[sim_part$truth$careless])
  att_med <- stats::median(sim_full$seconds[!sim_full$truth$careless])
  expect_gt(part_med, full_med)
  expect_lt(part_med, att_med)
})
