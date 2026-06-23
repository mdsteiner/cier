# Tests for the careless pattern mutators + extent/onset + truth assembly -- the
# careless half of cier_simulate(). These call internal kernels directly in the
# package namespace; the orchestrator and cier_sim object are pinned in
# test-cier-simulate.R. The content engine reads each careless row's knobs from
# truth$params (the single provenance source sim_build_plan records); there is no
# separate pattern_params argument.
#
# Trust model: oracle-only (a generator, not an index, no CRAN partner). The oracle
# (ref-sim-patterns.R) re-derives every deterministic mutator from scratch, never
# calling production. Numerically delicate behaviours are additionally pinned
# THROUGH the shipped indices: the Biemann footnote-2 bounce and cyclic diagonal
# via cier_lazr (2/3 and 1), the value straightliner and positional scale-boundary
# break via cier_longstring, the alternating period via cier_autocorrelation. Fast
# tier: deterministic pins (tol 0), index pins, range invariants, extent
# bookkeeping, truth-schema contract, validation edges, seed reproducibility, and
# small-n recovery floors.

source(test_path("..", "reference", "ref-sim-patterns.R"))

# ---- Local builders ---------------------------------------------------------

# Raw items data.frame in the lean min/max/scale schema; `max` / `min` recycle.
sim_items <- function(scale, max, min = NULL) {
  df <- data.frame(scale = as.character(scale), max = max,
                   stringsAsFactors = FALSE)
  if (!is.null(min)) df$min <- min
  df
}

# The validated items list the engine consumes (scale, reverse_keyed, min, max,
# categories), via the production validator.
vit <- function(scale, max, min = NULL) {
  it <- sim_items(scale, max, min)
  check_items_simulate(it, nrow(it))
}

# A one-row careless-high score through an index, on a single response row -- one
# finite value resolves the percentile cutoff, so no abstention warning fires.
row_value <- function(index_fun, row, ...) {
  # simulated fixtures can trip the percentile-cutoff degeneracy guard; these
  # recovery tests assert rank recovery, not flags, so the (correct) warning is
  # muffled.
  suppressWarnings(index_fun(matrix(as.integer(row), nrow = 1L), ...))$value
}

# Build a contaminated dataset end-to-end (attentive layer + careless plant)
# for the recovery smoke tests. One pattern, all-full, two orthogonal scales.
make_contaminated <- function(seed, pattern, n = 200L, p = 20L,
                              prevalence = 0.25, pattern_params = list()) {
  it_df <- sim_items(rep(c("E", "A"), each = p / 2L), max = 5)
  items <- check_items_simulate(it_df, p)
  withr::with_seed(seed, {
    truth <- sim_build_plan(n, p, prevalence,
                            stats::setNames(1, pattern),
                            pattern_params = pattern_params)
    att <- sim_attentive(n, it_df)
    # the engine reads each row's knobs from truth$params (the single source
    # sim_build_plan recorded); no separate pattern_params argument.
    resp <- sim_apply_patterns(att, items, truth)
    list(responses = resp, truth = truth)
  })
}

# =============================================================================
# FAST TIER
# =============================================================================

# ---- Rounding convention ----------------------------------------------------

test_that("sim_round_half_up rounds .5 boundaries UP, not to even", {
  expect_identical(sim_round_half_up(c(2.5, 3.5, 2.4, 2.6, -0.5, 0.5)),
                   c(3, 4, 2, 3, 0, 1))
  # base-R round() (half to even) gives 2 and 4 for 2.5 / 3.5; the archive midpoint
  # ceiling((K+1)/2) and this helper agree on half-up.
  expect_false(identical(sim_round_half_up(2.5), round(2.5)))
})

# ---- Registry / allowlist ---------------------------------------------------

test_that("the allowlist is the eight patterns; the registry is the seven producers", {
  expect_setequal(
    sim_pattern_names(),
    c("random", "straightline", "midpoint", "extreme", "diagonal",
      "alternating", "markov", "speeder")
  )
  # speeder is a content no-op (it acts on timing, not content), so it is NOT a producer.
  expect_setequal(names(sim_pattern_registry()),
                  setdiff(sim_pattern_names(), "speeder"))
  expect_true(is.function(sim_block_fun("random")))
  expect_null(sim_block_fun("speeder"))
})

# ---- Deterministic mutator pins (tol 0, exact integer) ----------------------

test_that("straightline position emits the round-half-up position per item", {
  mins <- c(1, 1, 1)
  maxs <- c(5, 5, 5)
  blk <- sim_block_straightline(1L, mins, maxs,
                                list(anchor = "position", q = 0.25))
  # round_half_up(1 + 0.25 * 4) = round_half_up(2) = 2, constant across items.
  expect_identical(as.integer(blk), as.integer(ref_position_value(0.25, mins, maxs)))
  expect_identical(storage.mode(blk), "integer")
})

test_that("straightline value emits the constant value across items", {
  blk <- sim_block_straightline(1L, rep(1, 6), rep(5, 6),
                                list(anchor = "value", value = 4))
  expect_identical(as.integer(blk), rep(4L, 6))
})

test_that("the default value draws ONE constant per row from the common range", {
  # No explicit value: each row holds a single integer from the achievable common
  # range [max(min), min(max)] = [2, 5] here. A per-item-min mutant (-> c(1,2,1),
  # non-constant) or a per-cell draw both break the row-constant contract.
  mins <- c(1, 2, 1)
  maxs <- c(5, 5, 5)
  blk <- withr::with_seed(1L, {
    sim_block_straightline(40L, mins, maxs, list(anchor = "value"))
  })
  expect_true(all(apply(blk, 1L, function(r) length(unique(r)) == 1L)))  # row-constant
  expect_true(all(blk >= 2L & blk <= 5L))                                # common range
  # On a homogeneous battery the default value is a true straightliner -> run = p.
  homo <- withr::with_seed(2L, {
    sim_block_straightline(1L, rep(1, 10), rep(5, 10), list(anchor = "value"))
  })
  expect_equal(row_value(cier_longstring, as.integer(homo)), 10, tolerance = 0)
})

test_that("a value outside ANY item's range is a typed error", {
  # value 4 is in range for item 1 (1..5) but out for item 2 (1..3): must error, not
  # silently clip. Pins the "outside any item's range" boundary on a heterogeneous
  # battery, not just a value above every item's max.
  expect_error(
    sim_block_straightline(1L, c(1, 1), c(5, 3), list(anchor = "value", value = 4)),
    class = "cier_error_input"
  )
})

test_that("the switch toggle alternates with switch_prob = 1 (deterministic)", {
  # switch_prob = 1 toggles every column regardless of RNG: position q then 1 - q
  # then q ... For q = 0 (-> min 1) and 1 - q = 1 (-> max 5): 1, 5, 1, 5, 1, 5.
  blk <- sim_block_straightline(1L, rep(1, 6), rep(5, 6),
                                list(anchor = "position", q = 0, switch_prob = 1))
  expect_identical(as.integer(blk), c(1L, 5L, 1L, 5L, 1L, 5L))
  # Value anchor reflects (min + max) - v: v = 2 toggles to 4.
  blk2 <- sim_block_straightline(1L, rep(1, 4), rep(5, 4),
                                 list(anchor = "value", value = 2, switch_prob = 1))
  expect_identical(as.integer(blk2), c(2L, 4L, 2L, 4L))
})

test_that("switch_prob = 0 is a pure straightliner (no toggle)", {
  blk <- sim_block_straightline(1L, rep(1, 8), rep(5, 8),
                                list(anchor = "position", q = 0.5, switch_prob = 0))
  expect_identical(as.integer(blk), rep(3L, 8))           # constant midpoint
})

test_that("midpoint is the upper middle on even-K, the centre on odd-K", {
  # K = 4 (min 1, max 4): round_half_up((1 + 4)/2) = round_half_up(2.5) = 3.
  # K = 5 (min 1, max 5): 3. Heterogeneous battery mixes both.
  blk <- sim_block_midpoint(1L, c(1, 1), c(4, 5), list())
  expect_identical(as.integer(blk), c(3L, 3L))
  blk0 <- sim_block_midpoint(1L, c(0, -3), c(3, 3), list())   # 0..3 (K4), -3..3 (K7)
  expect_identical(as.integer(blk0), c(2L, 0L))              # round_half_up(1.5)=2, mid 0
})

test_that("extreme is min or max only (p_high in {0, 1} pins it)", {
  hi <- sim_block_extreme(3L, c(1, 1, 0), c(5, 7, 4), list(p_high = 1))
  expect_identical(hi, matrix(c(5L, 5L, 5L, 7L, 7L, 7L, 4L, 4L, 4L), nrow = 3L))
  lo <- sim_block_extreme(3L, c(1, 1, 0), c(5, 7, 4), list(p_high = 0))
  expect_identical(lo, matrix(c(1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L), nrow = 3L))
})

test_that("diagonal cyclic wraps; bounce zigzags (Biemann footnote-2)", {
  cyc <- sim_block_diagonal(1L, rep(1, 10), rep(5, 10),
                            list(start = 1, step = 1, bounce = FALSE))
  expect_identical(as.integer(cyc), ref_diagonal_cyclic(1L, 1L, 5L, 10L))
  expect_identical(as.integer(cyc), c(1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L))
  bnc <- sim_block_diagonal(1L, rep(1, 10), rep(5, 10),
                            list(start = 1, step = 1, bounce = TRUE))
  expect_identical(as.integer(bnc), ref_diagonal_bounce(1L, 1L, 5L, 10L))
  expect_identical(as.integer(bnc), ref_biemann_bounce)     # 1,2,3,4,5,4,3,2,1,2
})

test_that("alternating cycles through its values at the given period", {
  blk <- sim_block_alternating(1L, rep(1, 9), rep(5, 9),
                               list(period = 3, values = c(1, 3, 5),
                                    start_offset = 0))
  expect_identical(as.integer(blk), ref_alternating(c(1L, 3L, 5L), 0L, 3L, 9L))
  expect_identical(as.integer(blk), c(1L, 3L, 5L, 1L, 3L, 5L, 1L, 3L, 5L))
})

test_that("alternating with period above K_min cycles 1..K_min in range", {
  # period 7 > K_min 5: the default positions wrap (rep_len of 1..K_min), still
  # >= 2 distinct and in range.
  blk <- sim_block_alternating(1L, rep(1, 9), rep(5, 9),
                               list(period = 7, start_offset = 0))
  expect_true(all(blk >= 1L & blk <= 5L))
})

# ---- Heterogeneous offset (sequential patterns offset by min_j) -------------

test_that("diagonal positions offset to each item's min on a bipolar base", {
  # K_min = 5 over the span; positions 1..5 -> values min_j + pos - 1.
  blk <- sim_block_diagonal(1L, c(0, -2), c(4, 2),
                            list(start = 1, step = 1, bounce = FALSE))
  # positions 1, 2 -> col1 (min 0): 0, 1 ; col2 (min -2): -2, -1.
  expect_identical(as.integer(blk), c(0L, -1L))
})

test_that("sequential patterns cycle on K_min when the span mixes category counts", {
  # span K = 5, 3, 5, 3 -> K_min = 3, so the cyclic walk wraps mod 3 (1,2,3,1), NOT
  # mod 5 (1,2,3,4). A max-cats or first-item-K mutant gives 1,2,3,4 and dies; K_min
  # also keeps every value inside its own (possibly smaller) range.
  blk <- sim_block_diagonal(1L, c(1, 1, 1, 1), c(5, 3, 5, 3),
                            list(start = 1, step = 1, bounce = FALSE))
  expect_identical(as.integer(blk), c(1L, 2L, 3L, 1L))
})

# ---- Stochastic producers: range + seed determinism -------------------------

test_that("random / extreme / markov stay in range and reproduce under a seed", {
  mins <- c(1, 1, 0, -2)
  maxs <- c(5, 7, 4, 2)
  producers <- list(
    random = function() sim_block_random(40L, mins, maxs),
    extreme = function() sim_block_extreme(40L, mins, maxs),
    markov = function() sim_block_markov(40L, mins, maxs)
  )
  for (nm in names(producers)) {
    a <- withr::with_seed(4L, producers[[nm]]())
    b <- withr::with_seed(4L, producers[[nm]]())
    expect_identical(a, b, info = nm)                       # seed reproducible
    for (j in seq_along(mins)) {
      expect_true(all(a[, j] >= mins[j] & a[, j] <= maxs[j]), info = nm)
    }
    expect_identical(storage.mode(a), "integer", info = nm)
  }
})

test_that("markov honours the supplied transition matrix", {
  # Every state transitions to state 1 with probability 1, so from item 2 on the
  # chain is locked to position 1 (value min = 1) regardless of the start draw or
  # seed. A mutant ignoring `transition` (uses the default diag-0.6 chain) produces
  # varied later columns and dies. Deterministic: no seed needed.
  trans <- matrix(0, 5L, 5L)
  trans[, 1L] <- 1
  blk <- sim_block_markov(20L, rep(1, 6), rep(5, 6), list(transition = trans))
  expect_true(all(blk[, 2:6] == 1L))
})

test_that("every deterministic producer keeps values in range under heterogeneity", {
  mins <- c(1, 1, 0, -2)
  maxs <- c(5, 7, 4, 2)
  blocks <- list(
    sim_block_straightline(20L, mins, maxs, list(anchor = "position", q = 0.7)),
    sim_block_straightline(20L, mins, maxs, list(anchor = "position",
                                                 q = 0.3, switch_prob = 0.5)),
    sim_block_midpoint(20L, mins, maxs, list(jitter = 1)),
    sim_block_alternating(20L, mins, maxs, list(period = 2))
  )
  for (blk in blocks) {
    for (j in seq_along(mins)) {
      expect_true(all(blk[, j] >= mins[j] & blk[, j] <= maxs[j]))
    }
  }
})

# ---- Paper-anchored pins THROUGH the shipped indices ------------------------

test_that("the bounce sequence gives cier_lazr = 2/3; cyclic gives 1", {
  bnc <- as.integer(sim_block_diagonal(1L, rep(1, 10), rep(5, 10),
                                       list(start = 1, bounce = TRUE)))
  expect_equal(row_value(cier_lazr, bnc), 2 / 3, tolerance = 1e-12)
  cyc <- as.integer(sim_block_diagonal(1L, rep(1, 10), rep(5, 10),
                                       list(start = 1, bounce = FALSE)))
  expect_equal(row_value(cier_lazr, cyc), 1, tolerance = 1e-12)
})

test_that("a value straightliner gives cier_longstring = p", {
  blk <- as.integer(sim_block_straightline(1L, rep(1, 12), rep(5, 12),
                                           list(anchor = "value", value = 3)))
  expect_equal(row_value(cier_longstring, blk), 12, tolerance = 0)   # exact run = p
})

test_that("a positional straightliner run breaks at the scale boundary", {
  # 5 five-point items then 5 seven-point items; q = 0.5 -> 3s then 4s. The longest
  # run is the within-scale block (5), NOT p (10): longstring breaks at the boundary.
  mins <- rep(1, 10)
  maxs <- c(rep(5, 5), rep(7, 5))
  blk <- sim_block_straightline(1L, mins, maxs, list(anchor = "position", q = 0.5))
  expect_identical(as.integer(blk), c(rep(3L, 5), rep(4L, 5)))
  expect_equal(row_value(cier_longstring, as.integer(blk)), 5, tolerance = 0)
})

test_that("position equals value cell-for-cell on a homogeneous battery", {
  mins <- rep(1, 6)
  maxs <- rep(5, 6)
  q <- 0.25
  pos <- sim_block_straightline(1L, mins, maxs, list(anchor = "position", q = q))
  v <- ref_position_value(q, 1, 5)                          # = 2
  val <- sim_block_straightline(1L, mins, maxs, list(anchor = "value", value = v))
  expect_identical(pos, val)
})

test_that("alternating at period d peaks autocorrelation at lag d", {
  blk <- as.integer(sim_block_alternating(1L, rep(1, 12), rep(5, 12),
                                          list(period = 3, values = c(1, 3, 5),
                                               start_offset = 0)))
  expect_equal(row_value(cier_autocorrelation, blk, min_lag = 3L, max_lag = 3L),
               1, tolerance = 1e-10)
  # the period shows at lag d, not at lag 1 (distinct values -> incomplete cycle).
  expect_lt(row_value(cier_autocorrelation, blk, min_lag = 1L, max_lag = 1L), 1)
})

# ---- Producer validation edges ----------------------------------------------

test_that("producers reject impossible parameters with typed errors", {
  expect_error(                                              # value out of range
    sim_block_straightline(1L, rep(1, 5), rep(5, 5),
                           list(anchor = "value", value = 9)),
    class = "cier_error_input"
  )
  expect_error(                                              # cyclic step a multiple of K
    sim_block_diagonal(1L, rep(1, 5), rep(5, 5), list(step = 5)),
    class = "cier_error_input"
  )
  expect_error(                                              # bounce step = wave period
    sim_block_diagonal(1L, rep(1, 5), rep(5, 5),
                       list(bounce = TRUE, step = 8)),       # 2 * (5 - 1) = 8
    class = "cier_error_input"
  )
  expect_error(                                              # < 2 distinct values
    sim_block_alternating(1L, rep(1, 5), rep(5, 5),
                          list(period = 2, values = c(3, 3))),
    class = "cier_error_input"
  )
  expect_error(                                              # transition wrong dim
    sim_block_markov(1L, rep(1, 5), rep(5, 5),
                     list(transition = matrix(0.5, 2L, 2L))),
    class = "cier_error_input"
  )
})

test_that("producers reject out-of-range per-row params and bad markov specs", {
  expect_error(                                              # start position > K
    sim_block_diagonal(1L, rep(1, 5), rep(5, 5), list(start = 99)),
    class = "cier_error_input"
  )
  expect_error(                                              # position q outside [0,1]
    sim_block_straightline(1L, rep(1, 5), rep(5, 5),
                           list(anchor = "position", q = 2)),
    class = "cier_error_input"
  )
  expect_error(                                              # diagonal step not >= 1
    sim_block_diagonal(1L, rep(1, 5), rep(5, 5), list(step = 0)),
    class = "cier_error_input"
  )
  expect_error(                                              # alternating period < 2
    sim_block_alternating(1L, rep(1, 5), rep(5, 5), list(period = 1)),
    class = "cier_error_input"
  )
  expect_error(                                              # midpoint jitter < 0
    sim_block_midpoint(1L, rep(1, 5), rep(5, 5), list(jitter = -1)),
    class = "cier_error_input"
  )
  expect_error(                                              # no common value on a
    sim_block_straightline(1L, c(1, 5), c(3, 7), list(anchor = "value")),
    class = "cier_error_input"                               # disjoint-range battery
  )
  bad_entry <- diag(5)
  bad_entry[1L, 2L] <- 2                                     # an entry above 1
  expect_error(sim_block_markov(1L, rep(1, 5), rep(5, 5),
                                list(transition = bad_entry)),
               class = "cier_error_input")
  expect_error(sim_block_markov(1L, rep(1, 5), rep(5, 5),    # rows do not sum to 1
                                list(transition = matrix(0.5, 5L, 5L))),
               class = "cier_error_input")
  expect_error(sim_block_markov(1L, rep(1, 5), rep(5, 5),    # initial not a simplex
                                list(initial = rep(0.5, 5L))),
               class = "cier_error_input")
})

test_that("markov accepts an explicit valid initial distribution", {
  # exercises the explicit-initial branch, not just the stationary default.
  blk <- withr::with_seed(7L, {
    sim_block_markov(30L, rep(1, 6), rep(5, 6), list(initial = rep(0.2, 5L)))
  })
  expect_true(all(blk >= 1L & blk <= 5L))
  expect_identical(storage.mode(blk), "integer")
})

# ---- Truth-schema contract --------------------------------------------------

test_that("sim_build_plan returns the truth schema (no speeded)", {
  truth <- withr::with_seed(1L, {
    sim_build_plan(20L, 10L, 0.5, c(straightline = 1))
  })
  expect_s3_class(truth, "data.frame")
  expect_identical(
    names(truth),
    c("careless", "pattern", "extent", "onset_item", "offset_item", "params")
  )
  expect_type(truth$careless, "logical")
  expect_type(truth$pattern, "character")
  expect_type(truth$extent, "character")
  expect_type(truth$onset_item, "integer")
  expect_type(truth$offset_item, "integer")
  expect_type(truth$params, "list")
  expect_false("speeded" %in% names(truth))                 # added by the speeded layer
})

test_that("attentive rows carry the clean labels; careless rows the planted ones", {
  truth <- withr::with_seed(1L, {
    sim_build_plan(20L, 10L, 0.5, c(straightline = 1))
  })
  att <- !truth$careless
  expect_equal(sum(truth$careless), 10L)                    # round(0.5 * 20)
  expect_true(all(truth$pattern[att] == "attentive"))
  expect_true(all(truth$extent[att] == "none"))
  expect_true(all(is.na(truth$onset_item[att])))
  expect_true(all(is.na(truth$offset_item[att])))
  expect_true(all(lengths(truth$params[att]) == 0L))        # empty list()
  cl <- truth$careless
  expect_true(all(truth$pattern[cl] == "straightline"))
  expect_true(all(truth$extent[cl] == "full"))              # default extent
  expect_true(all(truth$onset_item[cl] == 1L))              # full span [1, p]
  expect_true(all(truth$offset_item[cl] == 10L))
})

test_that("the params column records each row's OWN pattern knobs", {
  # Two patterns with distinct knobs: every careless row must record ITS pattern's
  # knobs, so a mutant that keys on the first pattern (or records one list for all
  # rows) is caught; a single-pattern plan could not distinguish that.
  pp <- list(straightline = list(anchor = "value", value = 3),
             diagonal = list(bounce = TRUE))
  truth <- withr::with_seed(1L, {
    sim_build_plan(40L, 6L, 1, c(straightline = 0.5, diagonal = 0.5),
                   pattern_params = pp)
  })
  expect_true(any(truth$pattern == "straightline"))
  expect_true(any(truth$pattern == "diagonal"))
  for (i in which(truth$careless)) {
    expect_identical(truth$params[[i]], pp[[truth$pattern[[i]]]])
  }
})

test_that("partial spans start inside the onset window and run to p", {
  truth <- withr::with_seed(2L, {
    sim_build_plan(60L, 20L, 1, c(random = 1), prop_partial = 1)
  })
  expect_true(all(truth$extent == "partial"))
  # default onset_window = round(c(0.3, 0.8) * 20) = c(6, 16).
  expect_true(all(truth$onset_item >= 6L & truth$onset_item <= 16L))
  expect_true(all(truth$offset_item == 20L))                # partial runs to p
})

test_that("an explicit onset_window is validated and honoured", {
  truth <- withr::with_seed(6L, {
    sim_build_plan(60L, 20L, 1, c(random = 1), prop_partial = 1,
                   onset_window = c(8L, 10L))
  })
  expect_true(all(truth$onset_item >= 8L & truth$onset_item <= 10L))
})

test_that("temporary spans leave a recovery tail (offset <= p - 1)", {
  truth <- withr::with_seed(3L, {
    sim_build_plan(60L, 20L, 1, c(random = 1), prop_temporary = 1)
  })
  expect_true(all(truth$extent == "temporary"))
  expect_true(all(truth$onset_item <= truth$offset_item))
  expect_true(all(truth$offset_item <= 19L))                # at least one recovered col
  expect_true(all(truth$onset_item >= 1L))
})

test_that("extent uses exact proportional counts with the remainder full", {
  # prop_partial = prop_temporary = 0.3 over 100 careless rows -> exactly 30 / 30,
  # remainder 40 full. Exact counts (not a sampled multinomial) so the split is
  # deterministic; kills a mutant that drops the remainder-full bucket or mis-allocates.
  truth <- withr::with_seed(4L, {
    sim_build_plan(100L, 20L, 1, c(random = 1),
                   prop_partial = 0.3, prop_temporary = 0.3)
  })
  tab <- table(factor(truth$extent, levels = c("full", "partial", "temporary")))
  expect_identical(as.integer(tab), c(40L, 30L, 30L))
})

test_that("prop_partial + prop_temporary == 1 is accepted and leaves no full rows", {
  truth <- withr::with_seed(5L, {
    sim_build_plan(80L, 20L, 1, c(random = 1),
                   prop_partial = 0.5, prop_temporary = 0.5)
  })
  expect_false(any(truth$extent == "full"))                 # exact-1 boundary, no error
  expect_identical(sum(truth$extent == "partial"), 40L)
  expect_identical(sum(truth$extent == "temporary"), 40L)
})

test_that("prevalence 0 plants nobody; prevalence 1 plants everybody", {
  t0 <- withr::with_seed(1L, sim_build_plan(10L, 8L, 0))
  expect_false(any(t0$careless))
  expect_true(all(t0$pattern == "attentive"))
  t1 <- withr::with_seed(1L, sim_build_plan(10L, 8L, 1, c(random = 1)))
  expect_true(all(t1$careless))
})

# ---- Speeded assignment (prop_speeded) --------------------------------------

test_that("prop_speeded = 1 returns careless unchanged and draws no RNG", {
  # The default must be byte-identical to a no-knob run: sim_draw_speeded draws
  # NOTHING at prop_speeded = 1, so the downstream timing / check stream is
  # unshifted. A mutant that always permutes (even at 1) consumes RNG and fails the
  # stream-untouched probe below.
  careless <- c(TRUE, FALSE, TRUE, TRUE, FALSE)
  pattern  <- c("straightline", "attentive", "speeder", "random", "attentive")
  expect_identical(sim_draw_speeded(careless, pattern, 1), careless)
  with_call <- withr::with_seed(99L, {
    sim_draw_speeded(careless, pattern, 1)
    stats::runif(1L)
  })
  no_call <- withr::with_seed(99L, stats::runif(1L))
  expect_identical(with_call, no_call)
})

test_that("speeder rows always speed; attentive never; content honours the share", {
  careless <- c(rep(TRUE, 6L), TRUE, TRUE, FALSE, FALSE)
  pattern  <- c(rep("straightline", 6L), "speeder", "speeder",
                "attentive", "attentive")
  sp0 <- withr::with_seed(1L, sim_draw_speeded(careless, pattern, 0))
  expect_identical(sp0[pattern == "speeder"], c(TRUE, TRUE))   # speeders always speed
  expect_false(any(sp0[pattern == "straightline"]))            # content: none at p = 0
  expect_false(any(sp0[!careless]))                            # attentive never speeds
  expect_true(all(sp0 <= careless))                            # speeded subset of careless
})

test_that("prop_speeded uses an exact content-row count (round), not floor / ceiling", {
  # 10 content rows: round(10 * 0.27) = 3 (kills a floor mutant -> 2); round(10 *
  # 0.23) = 2 (kills a ceiling mutant -> 3). No half-way ties.
  careless <- rep(TRUE, 10L)
  pattern  <- rep("random", 10L)
  up <- withr::with_seed(2L, sim_draw_speeded(careless, pattern, 0.27))
  expect_identical(sum(up), 3L)
  dn <- withr::with_seed(2L, sim_draw_speeded(careless, pattern, 0.23))
  expect_identical(sum(dn), 2L)
  # exact count at large n: a per-row Bernoulli hits 50 only ~8% of the time, so an
  # exact 50 forces the round()-count + permute contract (not a Bernoulli).
  big <- withr::with_seed(4L, sim_draw_speeded(rep(TRUE, 100L), rep("random", 100L), 0.5))
  expect_identical(sum(big), 50L)
})

test_that("with no content-careless rows prop_speeded is inert (speeders only)", {
  careless <- c(TRUE, TRUE, FALSE)
  pattern  <- c("speeder", "speeder", "attentive")
  expect_identical(withr::with_seed(3L, sim_draw_speeded(careless, pattern, 0)),
                   careless)                                   # no draw, speeders speed
  # all-attentive degenerate (no careless at all): unchanged, no draw.
  expect_identical(sim_draw_speeded(c(FALSE, FALSE), c("attentive", "attentive"), 0.3),
                   c(FALSE, FALSE))
})

# ---- Apply engine: extent bookkeeping + off-by-one + group-by-row -----------

test_that("the content engine mutates exactly each row's own span", {
  items <- vit(rep("F", 10L), max = 5)
  att <- matrix(2L, nrow = 4L, ncol = 10L)                  # known attentive fill
  truth <- data.frame(
    careless = c(TRUE, TRUE, TRUE, FALSE),
    pattern = c("straightline", "straightline", "straightline", "attentive"),
    extent = c("partial", "temporary", "full", "none"),
    onset_item = c(4L, 3L, 1L, NA),
    offset_item = c(10L, 6L, 10L, NA),
    stringsAsFactors = FALSE
  )
  pp <- list(anchor = "value", value = 5)
  truth$params <- list(pp, pp, pp, list())
  m <- sim_apply_patterns(att, items, truth)
  # row 1 partial [4, 10]: prefix 1:3 untouched, span 4:10 = 5 (off-by-one pinned).
  expect_identical(m[1L, 1:3], rep(2L, 3L))
  expect_identical(m[1L, 4:10], rep(5L, 7L))
  # row 2 temporary [3, 6]: prefix 1:2 AND suffix 7:10 untouched.
  expect_identical(m[2L, c(1:2, 7:10)], rep(2L, 6L))
  expect_identical(m[2L, 3:6], rep(5L, 4L))
  # row 3 full: every column mutated.
  expect_identical(m[3L, ], rep(5L, 10L))
  # row 4 attentive: byte-identical to the input.
  expect_identical(m[4L, ], att[4L, ])
})

test_that("each row mutates over its OWN span, restarting the pattern at the onset", {
  # Two diagonal rows, same fixed start but different spans. A group-by-pattern
  # (shared-span) mutant, or one that indexes the diagonal from global column 1
  # instead of restarting at the onset, gets row 2 wrong. Deterministic via a fixed
  # start, so the exact restarted sequence is pinned.
  items <- vit(rep("F", 6L), max = 5)
  att <- matrix(3L, nrow = 2L, ncol = 6L)
  truth <- data.frame(
    careless = c(TRUE, TRUE),
    pattern = c("diagonal", "diagonal"),
    extent = c("full", "partial"),
    onset_item = c(1L, 3L), offset_item = c(6L, 6L),
    stringsAsFactors = FALSE
  )
  pp <- list(start = 1, step = 1, bounce = FALSE)
  truth$params <- list(pp, pp)
  m <- sim_apply_patterns(att, items, truth)
  expect_identical(m[1L, ], c(1L, 2L, 3L, 4L, 5L, 1L))      # full diagonal from col 1
  expect_identical(m[2L, 1:2], c(3L, 3L))                   # prefix untouched (own span)
  expect_identical(m[2L, 3:6], c(1L, 2L, 3L, 4L))           # restarts at the onset, not 3,4,5,1
})

test_that("speeder rows leave content untouched (speeders act on timing only)", {
  items <- vit(rep("F", 6L), max = 5)
  att <- matrix(3L, nrow = 2L, ncol = 6L)
  truth <- data.frame(
    careless = c(TRUE, FALSE),
    pattern = c("speeder", "attentive"),
    extent = c("full", "none"),
    onset_item = c(1L, NA), offset_item = c(6L, NA),
    stringsAsFactors = FALSE
  )
  truth$params <- list(list(), list())
  expect_identical(sim_apply_patterns(att, items, truth), att)
})

test_that("rows sharing a pattern and span but not params keep their own knobs", {
  # The engine groups producer calls by pattern AND span AND params: a mutant that
  # groups by pattern + span alone would hand row 2 row 1's straightline value.
  # (cier_simulate() always records one params list per pattern, but the engine must
  # honour the truth frame it is actually given.)
  items <- vit(rep("F", 6L), max = 5)
  att <- matrix(3L, nrow = 2L, ncol = 6L)
  truth <- data.frame(
    careless = c(TRUE, TRUE),
    pattern = c("straightline", "straightline"),
    extent = c("full", "full"),
    onset_item = c(1L, 1L), offset_item = c(6L, 6L),
    stringsAsFactors = FALSE
  )
  truth$params <- list(list(anchor = "value", value = 2L),
                       list(anchor = "value", value = 4L))
  m <- sim_apply_patterns(att, items, truth)
  expect_identical(m[1L, ], rep(2L, 6L))
  expect_identical(m[2L, ], rep(4L, 6L))
})

# ---- Plan validation edges --------------------------------------------------

test_that("sim_build_plan rejects malformed weights / proportions / windows", {
  expect_error(sim_build_plan(10L, 8L, 0.5, c(1)),                       # unnamed
               class = "cier_error_input")
  expect_error(sim_build_plan(10L, 8L, 0.5, c(foo = 1)),                 # unknown
               class = "cier_error_input")
  expect_error(sim_build_plan(10L, 8L, 0.5, c(random = -0.5, midpoint = 1.5)),
               class = "cier_error_input")                               # negative
  expect_error(sim_build_plan(10L, 8L, 0.5, c(random = 0.3, midpoint = 0.3)),
               class = "cier_error_input")                               # not sum 1
  expect_error(sim_build_plan(10L, 8L, 1.5, c(random = 1)),              # prevalence
               class = "cier_error_input")
  expect_error(
    sim_build_plan(10L, 8L, 1, c(random = 1), prop_partial = 0.6,
                   prop_temporary = 0.6),                                # props > 1
    class = "cier_error_input"
  )
  expect_error(sim_build_plan(10L, 2L, 1, c(random = 1), prop_temporary = 1),
               class = "cier_error_input")                               # temporary p<3
  expect_error(sim_build_plan(10L, 8L, 0.5, c(random = 1), onset_window = c(5, 3)),
               class = "cier_error_input")                               # lo > hi
  expect_error(sim_build_plan(10L, 8L, 0.5, c(random = 1), onset_window = c(0, 5)),
               class = "cier_error_input")                               # out of [1, p]
})

# ---- Seed reproducibility ---------------------------------------------------

test_that("the plan and the content engine reproduce bytewise under a seed", {
  ta <- withr::with_seed(9L, {
    sim_build_plan(30L, 12L, 0.4, c(straightline = 0.5, diagonal = 0.5),
                   prop_partial = 0.5)
  })
  tb <- withr::with_seed(9L, {
    sim_build_plan(30L, 12L, 0.4, c(straightline = 0.5, diagonal = 0.5),
                   prop_partial = 0.5)
  })
  expect_identical(ta, tb)
  tc <- withr::with_seed(10L, {
    sim_build_plan(30L, 12L, 0.4, c(straightline = 0.5, diagonal = 0.5),
                   prop_partial = 0.5)
  })
  expect_false(identical(ta, tc))
  items <- vit(rep("F", 12L), max = 5)
  att <- matrix(3L, nrow = 30L, ncol = 12L)
  ma <- withr::with_seed(5L, sim_apply_patterns(att, items, ta))
  mb <- withr::with_seed(5L, sim_apply_patterns(att, items, ta))
  expect_identical(ma, mb)
  # the engine genuinely consumes RNG (default position q / diagonal start draws),
  # so a different seed must move the content -- guards against an accidentally
  # deterministic apply.
  mc <- withr::with_seed(6L, sim_apply_patterns(att, items, ta))
  expect_false(identical(ma, mc))
})

# ---- Recovery floors (small n, fast tier) -----------------------------------

test_that("each pattern is recovered by its matched index (rank-AUC floor)", {
  expect_gt(local({
    d <- make_contaminated(11L, "straightline",
                           pattern_params = list(straightline = list(anchor = "value")))
    ref_rank_auc(cier_longstring(d$responses)$value, d$truth$careless)
  }), 0.85)
  expect_gt(local({
    d <- make_contaminated(12L, "diagonal")
    ref_rank_auc(suppressWarnings(cier_lazr(d$responses))$value, d$truth$careless)
  }), 0.85)
  expect_gt(local({
    d <- make_contaminated(13L, "markov")
    ref_rank_auc(cier_lazr(d$responses)$value, d$truth$careless)
  }), 0.75)
  expect_gt(local({
    # max_lag = 8 (the help's recommendation): the default ncol - 3 includes tiny
    # high-lag slices whose |autocorrelation| saturates near 1 for attentive rows
    # too, eroding the seesaw separation.
    d <- make_contaminated(14L, "alternating")
    ref_rank_auc(suppressWarnings(cier_autocorrelation(d$responses, max_lag = 8L))$value,
                 d$truth$careless)
  }), 0.85)
  expect_gt(local({
    # random responders are multivariate outliers -> high Mahalanobis distance (the
    # canonical match; IRV's careless tail is the low-variance straightliner). Coarse
    # smoke floor: at 25% prevalence the careless rows inflate the covariance,
    # deflating their own distances, so ~0.79 here, not near 1.
    d <- make_contaminated(15L, "random")
    ref_rank_auc(cier_mahalanobis(d$responses)$value, d$truth$careless)
  }), 0.75)
})
