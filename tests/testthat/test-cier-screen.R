# Tests for cier_screen() -- the transparent flag-table combiner (slice 12).
#
# cier_screen() is an orchestrator, not a new statistic, so its trust model is
# (1) INTERNAL PARITY: every index it runs is byte-identical to calling that
# index directly (the screen must not alter a single per-respondent value or
# flag), pinned with expect_identical (tolerance 0); and (2) the COMBINER, whose
# only new logic -- collapsing correlated indices to one vote, counting them, and
# feeding flag_agreement() the right inputs -- is checked against the independent
# re-derivation in ref-screen-combiner.R (never the production collapse). There
# is no cross-package partner: careless/psych/PerFit/mokken have no equivalent
# combined screen with this collapse (see tests/reference/TOLERANCES.md).
#
# The headline contract the mutants target: even-odd + personal_reliability are
# ONE construct (they share the registry `vote_group` = "consistency"), so they
# collapse to a single vote -- a respondent flagged by both counts once, an
# abstaining (NA) member contributes FALSE, and the agreement runs on the
# collapsed votes, not the raw per-index flags.

source(test_path("..", "reference", "ref-screen-combiner.R"))

# A reproducible 1..5 Likert matrix (spans the scale, so the person-fit backends
# can zero-base it).
screen_matrix <- function(n = 80L, p = 20L, seed = 1L) {
  withr::with_seed(seed, {
    m <- matrix(sample.int(5L, n * p, replace = TRUE), nrow = n, ncol = p)
  })
  storage.mode(m) <- "double"
  m
}

# Item metadata that serves all four metadata indices at once: >= 2 scales and a
# reverse_keyed column (even-odd, PR) plus a homogeneous categories column
# (Gnormed, Ht). 4 scales x 5 items = 20 columns, alternating reverse keys.
screen_items <- function(n_scales = 4L, per_scale = 5L, categories = 5L,
                         reverse = NULL) {
  scale <- rep(paste0("s", seq_len(n_scales)), each = per_scale)
  if (is.null(reverse)) {
    reverse <- rep(c(FALSE, TRUE), length.out = length(scale))
  }
  data.frame(scale = scale, reverse_keyed = reverse, categories = categories,
             stringsAsFactors = FALSE)
}

# The six matrix-only indices and the eight non-backend indices, in registry
# order (the order cier_screen() must emit).
matrix_methods <- function() {
  c("cier_longstring", "cier_irv", "cier_psychsyn", "cier_psychant",
    "cier_mahalanobis", "cier_person_total")
}
base_methods <- function() {
  c("cier_longstring", "cier_irv", "cier_even_odd",
    "cier_personal_reliability", "cier_psychsyn", "cier_psychant",
    "cier_mahalanobis", "cier_person_total")
}

q <- function(expr) suppressWarnings(expr)

# A clean (uncontaminated) fixture with enough scales that even-odd does not clump
# on the across-scale correlation: 10 scales x 3 items, 300 respondents.
clean_matrix <- function() screen_matrix(300L, 30L, 2024L)
clean_items <- function() screen_items(n_scales = 10L, per_scale = 3L)

# ---- Schema -----------------------------------------------------------------

test_that("cier_screen returns a cier_screen object with the pinned schema", {
  sc <- q(cier_screen(screen_matrix()))
  expect_s3_class(sc, "cier_screen")
  expect_type(sc, "list")
  expect_identical(
    names(sc),
    c("indices", "flags", "vote_group", "votes", "agreement", "skipped",
      "methods", "n_respondents")
  )
  expect_identical(sc$n_respondents, 80L)
  expect_s3_class(sc$flags, "data.frame")
  expect_s3_class(sc$votes, "data.frame")
  expect_s3_class(sc$skipped, "data.frame")
  expect_identical(names(sc$skipped), c("method", "reason"))
  # Every flag column is logical (NA allowed); every vote column is logical with
  # no NA (abstain collapsed to FALSE).
  expect_true(all(vapply(sc$flags, is.logical, logical(1L))))
  expect_true(all(vapply(sc$votes, is.logical, logical(1L))))
  expect_false(anyNA(as.matrix(sc$votes)))
})

test_that("as.data.frame.cier_screen returns a tidy long table", {
  x <- screen_matrix(40L, 20L, 2L)
  it <- screen_items()
  sc <- q(cier_screen(x, it, methods = base_methods(),
                      control = list(cier_personal_reliability = list(seed = 1))))
  df <- as.data.frame(sc)
  expect_s3_class(df, "data.frame")
  expect_identical(names(df),
                   c("respondent", "method", "value", "flagged", "vote_group"))
  expect_identical(nrow(df), length(sc$indices) * 40L)
  expect_setequal(unique(df$method), names(sc$indices))
  # Cell content must match the source index object and the registry vote_group
  # (a mutant filling value from the wrong index, or mis-labelling the construct,
  # survives a names-only check).
  irv_rows <- df[df$method == "cier_irv", ]
  expect_identical(irv_rows$value, sc$indices$cier_irv$value)
  expect_identical(irv_rows$flagged, sc$indices$cier_irv$flagged)
  expect_true(all(irv_rows$vote_group == "cier_irv"))
  expect_true(all(df$vote_group[df$method == "cier_even_odd"] == "consistency"))
  expect_true(all(df$vote_group[df$method == "cier_personal_reliability"] ==
                    "consistency"))
})

# ---- Orchestration parity (the oracle): screen == direct call, tol 0 --------

test_that("each matrix-only index in the screen is identical to a direct call", {
  x <- screen_matrix(80L, 20L, 1L)
  sc <- q(cier_screen(x))                       # items NULL -> the six matrix ones
  expect_identical(names(sc$indices), matrix_methods())
  for (m in matrix_methods()) {
    expect_identical(sc$indices[[m]], q(do.call(m, list(x))))
  }
})

test_that("metadata indices in the screen are identical to direct calls (seeded)", {
  x <- screen_matrix(80L, 20L, 4L)
  it <- screen_items()
  sc <- q(cier_screen(x, it,
                      control = list(cier_personal_reliability = list(seed = 1))))
  expect_identical(sc$indices$cier_even_odd, q(cier_even_odd(x, it)))
  expect_identical(sc$indices$cier_personal_reliability,
                   q(cier_personal_reliability(x, it, seed = 1)))
})

test_that("backend indices in the screen are identical to direct calls (seeded)", {
  skip_if_not_installed("PerFit")
  skip_if_not_installed("mokken")
  x <- screen_matrix(80L, 20L, 6L)
  it <- screen_items()
  sc <- q(cier_screen(x, it, control = list(cier_gnormed = list(seed = 1))))
  expect_identical(sc$indices$cier_gnormed, q(cier_gnormed(x, it, seed = 1)))
  expect_identical(sc$indices$cier_ht, q(cier_ht(x, it)))
})

# ---- Selectability (off-select weaker votes; isolate RPR) -------------------

test_that("methods = a single index runs only that index (one vote)", {
  x <- screen_matrix(60L, 20L, 3L)
  it <- screen_items()
  sc <- q(cier_screen(x, it, methods = "cier_personal_reliability",
                      control = list(cier_personal_reliability = list(seed = 1))))
  expect_identical(names(sc$indices), "cier_personal_reliability")
  expect_identical(ncol(sc$votes), 1L)
  expect_identical(colnames(sc$votes), "consistency")
  expect_identical(sc$agreement$agreement$k, 1L)
})

test_that("a method subset runs only those, in registry order", {
  x <- screen_matrix(40L, 20L, 7L)
  sc <- q(cier_screen(x, methods = c("cier_irv", "cier_longstring")))
  expect_identical(names(sc$indices), c("cier_longstring", "cier_irv"))
  expect_identical(sort(sc$methods), c("cier_irv", "cier_longstring"))
})

test_that("the default selection runs every screenable index", {
  skip_if_not_installed("PerFit")
  skip_if_not_installed("mokken")
  x <- screen_matrix(60L, 20L, 8L)
  it <- screen_items()
  sc <- q(cier_screen(x, it, control = list(cier_gnormed = list(seed = 1))))
  expect_identical(length(sc$indices), 10L)
  expect_identical(nrow(sc$skipped), 0L)
})

test_that("an unknown method id is a typed input error", {
  expect_error(cier_screen(screen_matrix(), methods = "cier_nope"),
               class = "cier_error_input")
})

# ---- Skip-with-reason -------------------------------------------------------

test_that("items = NULL skips the four metadata indices with a reason", {
  x <- screen_matrix(60L, 20L, 9L)
  sc <- q(cier_screen(x))
  expect_identical(names(sc$indices), matrix_methods())
  expect_setequal(sc$skipped$method,
                  c("cier_even_odd", "cier_personal_reliability",
                    "cier_gnormed", "cier_ht"))
  expect_true(all(grepl("item", sc$skipped$reason)))
  # A skipped index never appears in the flag table or the votes.
  expect_false(any(c("cier_even_odd", "cier_gnormed") %in% colnames(sc$flags)))
})

test_that("an absent backend skips Gnormed/Ht with the package reason", {
  testthat::local_mocked_bindings(cier_namespace_present = function(...) FALSE)
  x <- screen_matrix(40L, 20L, 10L)
  it <- screen_items()
  sc <- q(cier_screen(x, it))
  expect_true(all(c("cier_gnormed", "cier_ht") %in% sc$skipped$method))
  expect_match(sc$skipped$reason[sc$skipped$method == "cier_gnormed"], "PerFit")
  expect_match(sc$skipped$reason[sc$skipped$method == "cier_ht"], "mokken")
  # even-odd / PR need no backend, so they still run when items are present.
  expect_false(any(c("cier_even_odd", "cier_personal_reliability")
                   %in% sc$skipped$method))
  expect_true(all(c("cier_even_odd", "cier_personal_reliability")
                  %in% names(sc$indices)))
})

test_that("a typed backend limit records a skip instead of crashing the battery", {
  # cier_ht raises a cier_error_backend_limit on a scale wider than mokken's
  # 10-category ceiling. Valid 11-point data must NOT abort the whole screen:
  # the screen catches exactly that subclass, records ht as skipped with the
  # limit as the reason, and every other selected index still runs and scores.
  skip_if_not_installed("mokken")
  x <- withr::with_seed(31L, {
    m <- matrix(sample.int(11L, 40L * 8L, replace = TRUE), nrow = 40L)
  })
  x[1L, 1L] <- 1   # force both extremes: the global range spans 1..11
  x[2L, 1L] <- 11
  storage.mode(x) <- "double"
  it <- data.frame(scale = rep(c("s1", "s2"), each = 4L),
                   reverse_keyed = FALSE, categories = 11L,
                   stringsAsFactors = FALSE)
  sc <- q(cier_screen(x, it,
                      methods = c("cier_longstring", "cier_irv", "cier_ht")))
  expect_identical(names(sc$indices), c("cier_longstring", "cier_irv"))
  expect_identical(sc$skipped$method, "cier_ht")
  expect_match(sc$skipped$reason, "10-category")
  # The survivors scored: their flag columns are present and the votes collapse.
  expect_identical(colnames(sc$flags), c("cier_longstring", "cier_irv"))
  # A malformed items frame is NOT this path: it still propagates as an error.
  bad_items <- data.frame(scale = "s1", stringsAsFactors = FALSE)
  expect_error(cier_screen(x, bad_items, methods = "cier_ht"),
               class = "cier_error_input")
})

test_that("items = NULL takes precedence over the backend reason", {
  # When BOTH inputs are missing (no items AND no backend), the items reason is
  # reported -- the structural precondition checked first.
  testthat::local_mocked_bindings(cier_namespace_present = function(...) FALSE)
  sc <- q(cier_screen(screen_matrix(20L, 20L, 11L)))
  expect_match(sc$skipped$reason[sc$skipped$method == "cier_gnormed"], "item")
  expect_match(sc$skipped$reason[sc$skipped$method == "cier_ht"], "item")
  # The four metadata indices skip; the six matrix indices still run.
  expect_setequal(sc$skipped$method,
                  c("cier_even_odd", "cier_personal_reliability",
                    "cier_gnormed", "cier_ht"))
  expect_identical(names(sc$indices), matrix_methods())
})

# ---- Redundancy collapse (the headline; enforced + tested) ------------------

test_that("collapse_votes ORs members, maps NA -> FALSE, never double-counts", {
  flags <- data.frame(
    cier_even_odd             = c(TRUE,  FALSE, NA,    FALSE),
    cier_personal_reliability = c(FALSE, TRUE,  NA,    FALSE),
    cier_irv                  = c(FALSE, FALSE, TRUE,  NA),
    stringsAsFactors = FALSE
  )
  map <- c(cier_even_odd = "consistency",
           cier_personal_reliability = "consistency",
           cier_irv = "cier_irv")
  v <- collapse_votes(flags, map)
  # Two votes: the fused consistency construct and the standalone irv.
  expect_identical(colnames(v), c("consistency", "cier_irv"))
  # OR across members; an NA member contributes FALSE (rows 3-4 not flagged).
  expect_identical(v$consistency, c(TRUE, TRUE, FALSE, FALSE))
  expect_identical(v$cier_irv, c(FALSE, FALSE, TRUE, FALSE))
  expect_false(anyNA(as.matrix(v)))
})

test_that("a respondent flagged by BOTH consistency members counts once", {
  flags <- data.frame(cier_even_odd = TRUE, cier_personal_reliability = TRUE)
  map <- c(cier_even_odd = "consistency",
           cier_personal_reliability = "consistency")
  expect_identical(ref_screen_n_flags(collapse_votes(flags, map)), 1L)
})

test_that("vote columns stay in registry order regardless of the methods= order", {
  # collapse_votes derives group order from the flag columns, which follow the
  # ran-methods order; that equals registry order only because
  # screen_resolve_methods re-imposes it upstream. Pin the end-to-end property
  # directly so a future refactor that honours the user's methods= order cannot
  # silently reorder the vote columns (the oracle mirrors the same
  # first-appearance rule and would move in lock-step, not catch it).
  x <- screen_matrix(40L, 20L, 12L)
  sc <- q(cier_screen(x, methods = c("cier_person_total", "cier_irv",
                                     "cier_longstring")))
  expect_identical(names(sc$indices),
                   c("cier_longstring", "cier_irv", "cier_person_total"))
  expect_identical(colnames(sc$votes),
                   c("cier_longstring", "cier_irv", "cier_person_total"))
})

test_that("end-to-end, the screen collapses even-odd + PR to ONE vote", {
  x <- screen_matrix(80L, 20L, 12L)
  it <- screen_items()
  sc <- q(cier_screen(x, it,
                      methods = c("cier_even_odd", "cier_personal_reliability"),
                      control = list(cier_personal_reliability = list(seed = 1))))
  expect_identical(ncol(sc$votes), 1L)              # not two
  eo <- q(cier_even_odd(x, it))$flagged
  pr <- q(cier_personal_reliability(x, it, seed = 1))$flagged
  eo[is.na(eo)] <- FALSE
  pr[is.na(pr)] <- FALSE
  expect_identical(sc$votes$consistency, eo | pr)   # OR, NA -> FALSE
})

test_that("the construct count collapses the pair: votes count < raw flag count", {
  # Force BOTH consistency members to flag every scored respondent (cutoff = -1,
  # the minimum even-odd / PR value, flags the whole upper tail). A respondent
  # then carries TWO raw per-index flags but only ONE construct vote -- pinning
  # that the count derives from the collapsed votes, never rowSums(flags). A
  # raw-count mutant gives 2 where the construct count is 1.
  x <- screen_matrix(80L, 20L, 23L)
  it <- screen_items()
  sc <- q(cier_screen(
    x, it, methods = c("cier_even_odd", "cier_personal_reliability"),
    control = list(cier_even_odd = list(cutoff = -1),
                   cier_personal_reliability = list(cutoff = -1, seed = 1))
  ))
  fl <- as.matrix(sc$flags)
  fl[is.na(fl)] <- FALSE
  raw_count <- rowSums(fl)                          # up to 2 per respondent
  votes_count <- rowSums(sc$votes)                  # collapsed: at most 1 here
  both <- !is.na(sc$indices$cier_even_odd$value) &
    !is.na(sc$indices$cier_personal_reliability$value)
  expect_gt(sum(both), 0L)
  expect_true(all(votes_count[both] == 1L))         # one construct, not two
  expect_true(all(raw_count[both] == 2L))           # two raw flags
  expect_true(all(votes_count <= raw_count))        # collapse never inflates
})

# ---- Agreement wiring (collapsed votes + correct null_rate) ------------------

test_that("agreement runs on the collapsed votes with the right null_rate", {
  skip_if_not_installed("PerFit")
  skip_if_not_installed("mokken")
  x <- screen_matrix(120L, 20L, 13L)
  it <- screen_items()
  ctrl <- list(cier_personal_reliability = list(seed = 1),
               cier_gnormed = list(seed = 1))
  sc <- q(cier_screen(x, it, control = ctrl))
  ran <- names(sc$indices)

  direct <- list(
    cier_longstring = q(cier_longstring(x)),
    cier_irv = q(cier_irv(x)),
    cier_even_odd = q(cier_even_odd(x, it)),
    cier_personal_reliability = q(cier_personal_reliability(x, it, seed = 1)),
    cier_psychsyn = q(cier_psychsyn(x)),
    cier_psychant = q(cier_psychant(x)),
    cier_mahalanobis = q(cier_mahalanobis(x)),
    cier_person_total = q(cier_person_total(x)),
    cier_gnormed = q(cier_gnormed(x, it, seed = 1)),
    cier_ht = q(cier_ht(x, it))
  )[ran]
  flags_direct <- as.data.frame(lapply(direct, function(o) o$flagged),
                                check.names = FALSE)
  reg <- cier_methods()
  vg <- stats::setNames(reg$vote_group[match(ran, reg$method)], ran)
  votes_ref <- ref_collapse_votes(flags_direct, vg)
  null_ref <- vapply(colnames(votes_ref), function(g) {
    members <- ran[vg[ran] == g]
    if (length(members) == 1L) {
      cm <- reg$default_cutoff_method[reg$method == members]
      if (cm %in% c("chisq", "perfit_null")) {
        return(reg$default_cutoff_value[reg$method == members])
      }
    }
    NA_real_
  }, numeric(1L))

  expect_equal(sc$agreement,
               flag_agreement(as.matrix(votes_ref), null_rate = unname(null_ref)))
})

test_that("only the null-referenced votes are marked informative", {
  skip_if_not_installed("PerFit")
  skip_if_not_installed("mokken")
  x <- screen_matrix(120L, 20L, 14L)
  it <- screen_items()
  sc <- q(cier_screen(x, it,
                      control = list(cier_personal_reliability = list(seed = 1),
                                     cier_gnormed = list(seed = 1))))
  pv <- sc$agreement$per_vote
  info <- stats::setNames(pv$informative, pv$vote)
  null <- stats::setNames(pv$null, pv$vote)
  expect_true(info[["cier_mahalanobis"]])
  expect_equal(null[["cier_mahalanobis"]], 0.001)
  expect_true(info[["cier_gnormed"]])
  expect_equal(null[["cier_gnormed"]], 0.05)
  # The percentile votes -- including the all-percentile consistency collapse --
  # are tautological (NA null, not informative).
  expect_false(info[["consistency"]])
  expect_true(is.na(null[["consistency"]]))
  expect_false(info[["cier_irv"]])
})

test_that("a literal cutoff override unhooks a vote's calibrated null", {
  # A null-referenced vote (mahalanobis) flagged on a literal `cutoff` no longer
  # targets a nominal rate, so its calibrated null must read NA / not informative
  # -- not the registry's 0.001, which would be a fiction against an absolute D².
  x <- screen_matrix(60L, 20L, 26L)
  sc <- q(cier_screen(x, methods = c("cier_irv", "cier_mahalanobis"),
                      control = list(cier_mahalanobis = list(cutoff = 60))))
  pv <- sc$agreement$per_vote
  row <- pv[pv$vote == "cier_mahalanobis", ]
  expect_true(is.na(row$null))
  expect_false(row$informative)
  # The default (no override) path still carries the chisq nominal.
  sc2 <- q(cier_screen(x, methods = c("cier_irv", "cier_mahalanobis")))
  pv2 <- sc2$agreement$per_vote
  expect_equal(pv2$null[pv2$vote == "cier_mahalanobis"], 0.001)
})

# ---- Bad metadata propagates (it is NOT a structural skip) -------------------

test_that("malformed items for a selected metadata index propagate the error", {
  x <- screen_matrix(40L, 20L, 15L)
  # Only one scale -> even-odd's typed input error must surface, not a skip.
  one_scale <- data.frame(scale = rep("A", 20L), reverse_keyed = FALSE,
                          categories = 5L)
  expect_error(cier_screen(x, one_scale, methods = "cier_even_odd"),
               class = "cier_error_input")
})

test_that("heterogeneous categories for Gnormed propagate the error", {
  skip_if_not_installed("PerFit")
  x <- screen_matrix(40L, 20L, 16L)
  het <- data.frame(scale = rep(c("A", "B"), each = 10L), reverse_keyed = FALSE,
                    categories = c(rep(5L, 10L), rep(4L, 10L)))
  expect_error(cier_screen(x, het, methods = "cier_gnormed"),
               class = "cier_error_input")
})

# ---- Control forwarding ------------------------------------------------------

test_that("a control override is forwarded to the index verbatim", {
  x <- screen_matrix(60L, 20L, 17L)
  sc <- q(cier_screen(x, methods = "cier_irv",
                      control = list(cier_irv = list(fpr = 0.10))))
  expect_identical(sc$indices$cier_irv, q(cier_irv(x, fpr = 0.10)))
})

test_that("a control seed makes the screen reproducible", {
  x <- screen_matrix(60L, 20L, 18L)
  it <- screen_items()
  args <- list(x, it, methods = "cier_personal_reliability",
               control = list(cier_personal_reliability = list(seed = 1)))
  v1 <- q(do.call(cier_screen, args))$indices$cier_personal_reliability$value
  v2 <- q(do.call(cier_screen, args))$indices$cier_personal_reliability$value
  expect_identical(v1, v2)
})

test_that("a control naming a non-selected method is a typed input error", {
  x <- screen_matrix(40L, 20L, 19L)
  expect_error(
    cier_screen(x, methods = "cier_irv",
                control = list(cier_longstring = list(frac = 0.6))),
    class = "cier_error_input"
  )
  expect_error(q(cier_screen(x, control = list(cier_nope = list()))),
               class = "cier_error_input")
})

# ---- Sane rate on a clean fixture (DoD) -------------------------------------

test_that("a clean sample flags a modest, sane rate", {
  x <- clean_matrix()
  it <- clean_items()
  sc <- q(cier_screen(x, it, methods = base_methods(),
                      control = list(cier_personal_reliability = list(seed = 1))))
  n <- nrow(x)
  # No index flags anywhere near a majority -- an inverted cutoff (flag the wrong
  # tail) would flag ~95%. (colSums, not colMeans: an all-NA abstaining column --
  # psychsyn/psychant find no pairs on random data -- must read as 0, not NaN.)
  n_flagged <- colSums(as.matrix(sc$flags), na.rm = TRUE)
  expect_true(all(n_flagged < 0.40 * n))
  # A clean majority is flagged by NO construct, and almost everyone by at most
  # one -- a screen that flagged everyone (or never collapsed) would fail both.
  nf <- rowSums(sc$votes)
  expect_gt(mean(nf == 0L), 0.50)
  expect_gt(mean(nf <= 1L), 0.90)
  # On independent (clean) data the observed >= 2-vote agreement sits at the
  # independence baseline -- no spurious agreement excess.
  ag <- sc$agreement$agreement
  expect_lt(ag$observed[[2L]], ag$expected[[2L]] + 0.05)
})

test_that("on real contaminated data the screen makes the agreement excess visible", {
  # The complement of the clean-sample test above, on the bundled REAL data:
  # bfi_careless carries a genuine careless subgroup, so flags must CLUSTER on
  # the same respondents -- the observed >= 2-vote share exceeds the
  # Poisson-binomial independence baseline. This pins the screen's headline
  # user-facing behaviour ("contamination is visible") end-to-end; the per-index
  # oracles cannot. Sanity bound: no construct flags an implausible majority.
  nm <- names(bfi_careless)[1:44]
  items <- data.frame(scale = sub("^v_BFI_([A-Za-z]+)[0-9].*$", "\\1", nm),
                      reverse_keyed = grepl("_R$", nm), categories = 5L)
  sc <- q(cier_screen(bfi_careless[, 1:44], items, methods = base_methods(),
                      control = list(cier_personal_reliability = list(seed = 1))))
  ag <- sc$agreement$agreement
  expect_gt(ag$observed[[2L]], ag$expected[[2L]])
  expect_true(all(colMeans(as.matrix(sc$votes)) < 0.5))
})

# ---- Edge cases -------------------------------------------------------------

test_that("when every selected index is skipped the screen is empty, not an error", {
  # Select only metadata indices but pass no items -> all skip; no votes, no
  # agreement, but a well-formed (empty) object.
  sc <- q(cier_screen(screen_matrix(30L, 20L, 20L),
                      methods = c("cier_even_odd", "cier_gnormed")))
  expect_identical(length(sc$indices), 0L)
  expect_identical(nrow(sc$skipped), 2L)
  expect_identical(ncol(sc$votes), 0L)
  expect_null(sc$agreement)
  expect_identical(sc$n_respondents, 30L)
  # as.data.frame of an empty screen is a well-formed zero-row long table.
  df <- as.data.frame(sc)
  expect_identical(nrow(df), 0L)
  expect_identical(names(df),
                   c("respondent", "method", "value", "flagged", "vote_group"))
})

test_that("a single-respondent matrix screens without crashing", {
  sc <- q(cier_screen(screen_matrix(1L, 20L, 21L)))
  expect_s3_class(sc, "cier_screen")
  expect_identical(sc$n_respondents, 1L)
})

test_that("the screen muffles a fully-abstaining index's cutoff warning but records it", {
  # On random data psychsyn / psychant find no qualifying pairs, so they score no
  # respondent. The standalone index warns; the screen instead reports this in
  # $flags / print ("0 / 0"), so NO warning should leak from cier_screen(). The
  # muffle is targeted (cier_warning_insufficient_items only), not suppressWarnings.
  x <- screen_matrix(40L, 20L, 30L)
  expect_no_warning(sc <- cier_screen(x))
  # The abstention is still recorded transparently, not silently dropped.
  expect_true("cier_psychsyn" %in% names(sc$indices))
  expect_true(all(is.na(sc$indices$cier_psychsyn$value)))
  expect_true(all(is.na(sc$flags$cier_psychsyn)))
  # The behaviour change is scoped to the orchestrator: a direct call still warns.
  expect_warning(cier_psychsyn(x), class = "cier_warning_insufficient_items")
})

test_that("a non-matrix payload is a typed input error", {
  expect_error(cier_screen(1:10), class = "cier_error_input")
})

test_that("a malformed control argument is a typed input error", {
  x <- screen_matrix(20L, 20L, 22L)
  expect_error(cier_screen(x, methods = "cier_irv", control = "nope"),
               class = "cier_error_input")
  expect_error(cier_screen(x, methods = "cier_irv",
                           control = list(list(fpr = 0.1))),   # unnamed element
               class = "cier_error_input")
  expect_error(cier_screen(x, methods = "cier_irv",
                           control = list(cier_irv = "nope")), # entry not a list
               class = "cier_error_input")
})

test_that("an unknown or unnamed control argument is a typed input error", {
  # A typo'd / unknown argument name is caught precisely (not left to an opaque
  # do.call failure), and an unnamed inner argument (which would splice
  # positionally) is rejected.
  x <- screen_matrix(20L, 20L, 23L)
  expect_error(cier_screen(x, methods = "cier_irv",
                           control = list(cier_irv = list(fp = 0.1))),  # typo
               class = "cier_error_input")
  expect_error(cier_screen(x, methods = "cier_irv",
                           control = list(cier_irv = list(0.1))),       # unnamed
               class = "cier_error_input")
  # A control entry naming a screen-managed positional (responses / items) must
  # error, not silently override the screen's own responses or collide on items.
  expect_error(cier_screen(x, methods = "cier_irv",
                           control = list(cier_irv = list(responses = x))),
               class = "cier_error_input")
  it <- screen_items()
  expect_error(cier_screen(x, it, methods = "cier_even_odd",
                           control = list(cier_even_odd = list(items = it))),
               class = "cier_error_input")
})

test_that("an NA-named or non-character control / methods argument is a typed error", {
  x <- screen_matrix(20L, 20L, 24L)
  expect_error(cier_screen(x, methods = 123L), class = "cier_error_input")
  # nzchar(NA) is TRUE, so the NA-named entry must be caught explicitly.
  expect_error(
    cier_screen(x, methods = "cier_irv",
                control = stats::setNames(list(list(fpr = 0.1)), NA_character_)),
    class = "cier_error_input"
  )
})

# ---- print values pinned independently (so the snapshots lock only format) --

test_that("print surfaces the skipped count and the skipped method names", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    out <- capture.output(q(print(cier_screen(screen_matrix(60L, 20L, 9L)))))
  })
  expect_true(any(grepl("Skipped: 4", out, fixed = TRUE)))
  for (m in c("cier_even_odd", "cier_personal_reliability",
              "cier_gnormed", "cier_ht")) {
    expect_true(any(grepl(m, out, fixed = TRUE)))
  }
})

test_that("print surfaces per-index counts and the agreement count from the fields", {
  x <- clean_matrix()
  it <- clean_items()
  sc <- q(cier_screen(x, it, methods = base_methods(),
                      control = list(cier_personal_reliability = list(seed = 1))))
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    out <- capture.output(print(sc))
  })
  n <- sc$n_respondents
  # Each index's flagged count (from $flags) appears on the index's own line.
  for (m in names(sc$indices)) {
    n_flagged <- sum(sc$flags[[m]], na.rm = TRUE)
    line <- out[grepl(m, out, fixed = TRUE)]
    expect_true(any(grepl(paste0("\\b", n_flagged, "\\b"), line)))
  }
  # The agreement counts shown are sum(rowSums(votes) >= k) -- derived from the
  # COLLAPSED votes. A print that recomputed from raw flags would differ at the
  # >= 2 level (the consistency pair). Pin the >= 1 and >= 2 counts to $agreement.
  obs_n <- round(sc$agreement$agreement$observed * n)
  expect_true(any(grepl(paste0("\\b", obs_n[[1L]], "\\b"), out)))
  expect_true(any(grepl(paste0("\\b", obs_n[[2L]], "\\b"), out)))
  # Sections are separated by blank lines (cli_verbatim would drop them, so this
  # pins that the separators actually render).
  expect_true(any(out == ""))
})

# ---- print snapshots (design-first; locked after mock-up approval) ----------

test_that("print renders the locked cli summary (default battery)", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    sc <- q(cier_screen(clean_matrix(), clean_items(), methods = base_methods(),
                        control = list(cier_personal_reliability = list(seed = 1))))
    expect_snapshot(print(sc))
  })
})

test_that("print reports skipped methods with their reasons", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    sc <- q(cier_screen(screen_matrix(60L, 20L, 9L)))   # items NULL -> 4 skipped
    expect_snapshot(print(sc))
  })
})

test_that("the excess marker is gated on the binomial tail, not a point comparison", {
  # Under independence the observed >=k count is Binomial(n, expected_k), so it
  # exceeds the expectation about half the time on CLEAN data: a strict
  # observed > expected comparison would advertise sampling noise as
  # contamination (e.g. "10.0% vs expected 9.8% <- excess"). The marker must
  # fire only when the one-sided binomial tail is below 0.05.
  ag <- list(agreement = data.frame(
    k = 1:2,
    observed = c(0.10, 0.05),     # n = 60: counts 6 and 3
    expected = c(0.098, 0.005)
  ))
  lines <- screen_agreement_lines(ag, 60L)
  # k = 1: P(X >= 6 | 60, 0.098) ~ 0.55 -- ordinary noise, no marker.
  expect_false(grepl("excess", lines[[3L]], fixed = TRUE))
  # k = 2: P(X >= 3 | 60, 0.005) ~ 0.0036 -- genuine excess, marker fires.
  expect_true(grepl("excess", lines[[4L]], fixed = TRUE))
  # A zero observed count never marks, even against a zero expectation.
  ag0 <- list(agreement = data.frame(k = 1L, observed = 0, expected = 0))
  expect_false(grepl("excess", screen_agreement_lines(ag0, 60L)[[3L]],
                     fixed = TRUE))
})
