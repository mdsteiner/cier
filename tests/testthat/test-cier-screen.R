# Tests for cier_screen() -- the transparent flag-table combiner.
#
# An orchestrator, not a new statistic, so trust model is (1) INTERNAL PARITY:
# every index it runs is byte-identical to a direct call (no per-respondent value
# or flag altered), pinned expect_identical tol 0; and (2) the COMBINER (collapse
# correlated indices to one vote, count them, feed flag_agreement() right inputs),
# checked against the independent re-derivation in ref-screen-combiner.R. No
# cross-package partner: careless/psych/PerFit/mokken lack this combined screen.
#
# Headline contract the mutants target: even-odd + personal_reliability are ONE
# construct (shared method-spec `vote_group` = "consistency"), so they collapse to a
# single vote -- flagged-by-both counts once, an abstaining (NA) member adds FALSE,
# and agreement runs on the collapsed votes, not the raw per-index flags.

source(test_path("..", "reference", "ref-screen-combiner.R"))

# A reproducible 1..5 Likert matrix (spans the scale so the person-fit backends
# can zero-base it).
screen_matrix <- function(n = 80L, p = 20L, seed = 1L) {
  withr::with_seed(seed, {
    m <- matrix(sample.int(5L, n * p, replace = TRUE), nrow = n, ncol = p)
  })
  storage.mode(m) <- "double"
  m
}

# Item metadata serving all four metadata indices: >= 2 scales + a reverse_keyed
# column (even-odd, PR) + a homogeneous max column (1-based; Gnormed, Ht).
# 4 scales x 5 items = 20 columns, alternating keys.
screen_items <- function(n_scales = 4L, per_scale = 5L, max = 5L,
                         reverse = NULL) {
  scale <- rep(paste0("s", seq_len(n_scales)), each = per_scale)
  if (is.null(reverse)) {
    reverse <- rep(c(FALSE, TRUE), length.out = length(scale))
  }
  data.frame(scale = scale, reverse_keyed = reverse, max = max,
             stringsAsFactors = FALSE)
}

# The six matrix-only and eight non-backend indices, in method-spec order (the
# order cier_screen() must emit).
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

# A clean fixture with enough scales that even-odd does not clump on the
# across-scale correlation: 10 scales x 3 items, 300 respondents.
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
      "notes", "methods", "n_respondents")
  )
  expect_s3_class(sc$notes, "data.frame")
  expect_identical(names(sc$notes), c("method", "note"))
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
  # Cell content must match the source index object and the method-spec vote_group
  # (a wrong-index value or mis-labelled construct survives a names-only check).
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

test_that("a method subset runs only those, in method-spec order", {
  x <- screen_matrix(40L, 20L, 7L)
  sc <- q(cier_screen(x, methods = c("cier_irv", "cier_longstring")))
  expect_identical(names(sc$indices), c("cier_longstring", "cier_irv"))
  expect_identical(sort(sc$methods), c("cier_irv", "cier_longstring"))
})

test_that("the default selection runs every screenable index", {
  skip_if_not_installed("PerFit")
  x <- screen_matrix(60L, 20L, 8L)
  it <- screen_items()
  sc <- q(cier_screen(x, it, control = list(cier_gnormed = list(seed = 1))))
  expect_identical(length(sc$indices), 10L)
  expect_identical(nrow(sc$skipped), 0L)
  # The v0.2 additions ship screenable = FALSE, so they never enter the run set.
  v02 <- c("cier_autocorrelation", "cier_lazr", "cier_total_time",
           "cier_page_time", "cier_attention")
  expect_false(any(v02 %in% names(sc$indices)))
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

test_that("no backend skip: Gnormed and Ht both run in pure R", {
  x <- screen_matrix(40L, 20L, 10L)
  it <- screen_items()
  sc <- q(cier_screen(x, it, control = list(cier_gnormed = list(seed = 1))))
  # Both person-fit indices score in pure R now, so neither has a runtime backend
  # gate and neither is skipped.
  expect_false("cier_gnormed" %in% sc$skipped$method)
  expect_false("cier_ht" %in% sc$skipped$method)
  expect_true(all(c("cier_gnormed", "cier_ht") %in% names(sc$indices)))
  expect_true(any(is.finite(sc$indices$cier_gnormed$value)))
  expect_true(any(is.finite(sc$indices$cier_ht$value)))
  # even-odd / PR need no backend, so they still run when items are present.
  expect_false(any(c("cier_even_odd", "cier_personal_reliability")
                   %in% sc$skipped$method))
  expect_true(all(c("cier_even_odd", "cier_personal_reliability")
                  %in% names(sc$indices)))
})

test_that("a wide (11-point) scale now scores Ht in the screen (ceiling lifted)", {
  # Ht previously raised cier_error_backend_limit above 10 categories and the screen
  # recorded it skipped. The in-package closed form has no ceiling: valid 11-point
  # data now SCORES Ht in the battery -- no skip. (The skip-instead-of-crash path for
  # a genuine backend limit is still covered by the Gnormed tests below.)
  x <- withr::with_seed(31L, {
    m <- matrix(sample.int(11L, 40L * 8L, replace = TRUE), nrow = 40L)
  })
  x[1L, 1L] <- 1   # both extremes present: the global range spans 1..11
  x[2L, 1L] <- 11
  storage.mode(x) <- "double"
  it <- data.frame(scale = rep(c("s1", "s2"), each = 4L),
                   reverse_keyed = FALSE, max = 11L,
                   stringsAsFactors = FALSE)
  sc <- q(cier_screen(x, it,
                      methods = c("cier_longstring", "cier_irv", "cier_ht")))
  expect_identical(names(sc$indices), c("cier_longstring", "cier_irv", "cier_ht"))
  expect_identical(nrow(sc$skipped), 0L)
  expect_true("cier_ht" %in% colnames(sc$flags))
  expect_true(any(is.finite(sc$indices$cier_ht$value)))
  # A malformed items frame is NOT scored away: it still propagates as an error.
  bad_items <- data.frame(scale = "s1", stringsAsFactors = FALSE)
  expect_error(cier_screen(x, bad_items, methods = "cier_ht"),
               class = "cier_error_input")
})

test_that("items = NULL gives the structural skip reason", {
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

test_that("vote columns stay in method-spec order regardless of the methods= order", {
  # collapse_votes derives group order from the flag columns (ran-methods order);
  # that equals method-spec order only because screen_resolve_methods re-imposes it
  # upstream. Pin the end-to-end property directly so a future refactor honouring
  # the user's methods= order cannot silently reorder the vote columns (the oracle
  # mirrors the same first-appearance rule and would move in lock-step, not catch).
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
  # the minimum even-odd / PR value, flags the whole upper tail). A respondent then
  # carries TWO raw flags but ONE construct vote -- pinning the count to collapsed
  # votes, never rowSums(flags). A raw-count mutant gives 2 where it is 1.
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
  specs <- cier_method_specs()
  vg <- stats::setNames(specs$vote_group[match(ran, specs$method)], ran)
  votes_ref <- ref_collapse_votes(flags_direct, vg)
  # Pin the production collapse against the independent oracle directly (not only
  # via the agreement statistic): the two must agree vote-for-vote and on vote-group
  # order. Compared value-wise (unname(as.matrix())) so the pin tracks the collapse
  # spec, not data.frame row-name attributes.
  expect_identical(unname(as.matrix(sc$votes)), unname(as.matrix(votes_ref)))
  expect_identical(colnames(sc$votes), colnames(votes_ref))
  null_ref <- vapply(colnames(votes_ref), function(g) {
    members <- ran[vg[ran] == g]
    if (length(members) == 1L) {
      cm <- specs$default_cutoff_method[specs$method == members]
      if (cm %in% c("chisq", "mc_null")) {
        return(specs$default_cutoff_value[specs$method == members])
      }
    }
    NA_real_
  }, numeric(1L))

  expect_equal(sc$agreement,
               flag_agreement(as.matrix(votes_ref), null_rate = unname(null_ref)))
})

test_that("only the null-referenced votes are marked informative", {
  skip_if_not_installed("PerFit")
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
  # The percentile votes (incl. the all-percentile consistency collapse) are
  # tautological (NA null, not informative).
  expect_false(info[["consistency"]])
  expect_true(is.na(null[["consistency"]]))
  expect_false(info[["cier_irv"]])
})

test_that("a literal cutoff override unhooks a vote's calibrated null", {
  # A null-referenced vote (mahalanobis) flagged on a literal `cutoff` no longer
  # targets a nominal rate, so its calibrated null must read NA / not informative,
  # not the registry's 0.001 (a fiction against an absolute D²).
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

test_that("a method = 'mc_null' control override hooks cier_ht's calibrated null", {
  # WI-3: cier_ht opted into the Monte-Carlo null via control runs a genuinely
  # calibrated null, so its vote must read the nominal / informative -- NOT the
  # registry-default percentile's tautological NA. The null rate keys off the
  # EFFECTIVE method (registry default OR a control `method` override); a regression
  # that keys off the static registry method ("percentile" for cier_ht) marks the
  # deliberately-calibrated vote as non-informative.
  x <- screen_matrix(120L, 20L, 14L)
  it <- screen_items()
  sc <- q(cier_screen(x, it, methods = c("cier_irv", "cier_ht"),
                      control = list(cier_ht = list(method = "mc_null", seed = 1))))
  row <- sc$agreement$per_vote[sc$agreement$per_vote$vote == "cier_ht", ]
  expect_true(row$informative)
  expect_equal(row$null, 0.05)               # the registry-default nominal
  # A non-default fpr on the override is honoured (not the registry 0.05).
  sc2 <- q(cier_screen(x, it, methods = c("cier_irv", "cier_ht"),
                       control = list(cier_ht = list(method = "mc_null",
                                                     fpr = 0.10, seed = 1))))
  pv2 <- sc2$agreement$per_vote
  expect_equal(pv2$null[pv2$vote == "cier_ht"], 0.10)
  # The default (no override) path stays the tautological percentile vote.
  sc3 <- q(cier_screen(x, it, methods = c("cier_irv", "cier_ht")))
  pv3 <- sc3$agreement$per_vote
  expect_true(is.na(pv3$null[pv3$vote == "cier_ht"]))
  expect_false(pv3$informative[pv3$vote == "cier_ht"])
})

# ---- Bad metadata propagates (it is NOT a structural skip) -------------------

test_that("malformed items for a selected metadata index propagate the error", {
  x <- screen_matrix(40L, 20L, 15L)
  # Only one scale -> even-odd's typed input error must surface, not a skip.
  one_scale <- data.frame(scale = rep("A", 20L), reverse_keyed = FALSE,
                          max = 5L)
  expect_error(cier_screen(x, one_scale, methods = "cier_even_odd"),
               class = "cier_error_input")
})

test_that("a heterogeneous span skips Gnormed with a reason; the battery survives", {
  # Accurate metadata for mixed-format data (five- and four-option items together)
  # hits the single-Ncat contract -- a backend limit, not a malformed frame --
  # so the screen records Gnormed skipped (mirroring the unattained-extreme sibling)
  # instead of aborting the battery. The 5-point block alone would zero-base
  # 1..5 -> 0..4; the data stays valid for every other index.
  skip_if_not_installed("PerFit")
  x <- screen_matrix(40L, 20L, 16L)
  x[, 11:20] <- pmin(x[, 11:20], 4)      # items 11-20 genuinely four-option
  het <- data.frame(scale = rep(c("A", "B"), each = 10L), reverse_keyed = FALSE,
                    max = c(rep(5L, 10L), rep(4L, 10L)))
  sc <- q(cier_screen(x, het, methods = c("cier_irv", "cier_gnormed")))
  expect_identical(names(sc$indices), "cier_irv")
  expect_identical(sc$skipped$method, "cier_gnormed")
  expect_match(sc$skipped$reason, "Ncat", fixed = TRUE)
  expect_identical(colnames(sc$flags), "cier_irv")
  # A MALFORMED gnormed frame (NA max on an item) is not a backend limit: it must
  # still propagate, so the skip path cannot over-broadly catch every metadata error.
  bad <- het
  bad$max[1L] <- NA
  expect_error(cier_screen(x, bad, methods = "cier_gnormed"),
               class = "cier_error_input")
})

test_that("an unattained scale extreme skips Gnormed with a reason; the battery survives", {
  # A valid sample whose responses never reach a declared extreme category (1..4 on
  # a declared 1..5 scale -- nobody picked the top option) leaves PerFit's item-step
  # popularities undefined: otherwise-valid data the backend cannot score
  # (sample-dependent, like the heterogeneous-span case), NOT a
  # metadata defect -- so the screen records Gnormed skipped-with-reason and every
  # OTHER index still runs, not a plain input error aborting the battery.
  # Items are all FORWARD-keyed on purpose: a reverse-keyed item reflects
  # (min + max) - x, mapping 1..4 to 5..2 and reintroducing the top category,
  # masking the unattained extreme this test pins.
  skip_if_not_installed("PerFit")
  x <- withr::with_seed(40L, {
    matrix(sample.int(4L, 60L * 20L, replace = TRUE), nrow = 60L, ncol = 20L)
  })
  x[1L, 1L] <- 1                    # bottom category present...
  x[1L, 2L] <- 4                    # ...top declared option (5) never occurs
  storage.mode(x) <- "double"
  it <- screen_items(reverse = rep(FALSE, 20L))   # declared max = 5; forward only
  sc <- q(cier_screen(x, it,
                      control = list(cier_personal_reliability = list(seed = 1))))
  # Gnormed is the ONLY skip; the other nine indices ran.
  expect_identical(sc$skipped$method, "cier_gnormed")
  expect_match(sc$skipped$reason, "scale extremes", fixed = TRUE)
  expect_identical(length(sc$indices), 9L)
  expect_false("cier_gnormed" %in% names(sc$indices))
  expect_false("cier_gnormed" %in% colnames(sc$flags))
  # The catch path leaves the flag table intact: the nine survivors reach it and Ht
  # actually SCORED -- pins genuine survival, not a mere listing. Ht scores fine on
  # the 1..4 data.
  expect_identical(ncol(sc$flags), 9L)
  expect_true("cier_ht" %in% colnames(sc$flags))
  expect_true(any(is.finite(sc$indices$cier_ht$value)))
  # An OUT-OF-RANGE Gnormed input (a value above the declared max) is a data defect,
  # not a backend limit: it must still PROPAGATE, never be silently skipped.
  x_oor <- x
  x_oor[1L, 1L] <- 99               # forward column 1; exceeds declared max = 5
  expect_error(cier_screen(x_oor, it, methods = "cier_gnormed"),
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
                      reverse_keyed = grepl("_R$", nm), max = 5L)
  sc <- q(cier_screen(bfi_careless[, 1:44], items, methods = base_methods(),
                      control = list(cier_personal_reliability = list(seed = 1))))
  ag <- sc$agreement$agreement
  expect_gt(ag$observed[[2L]], ag$expected[[2L]])
  expect_true(all(colMeans(as.matrix(sc$votes)) < 0.5))
})

# ---- Edge cases -------------------------------------------------------------

test_that("when every selected index is skipped the screen is empty, not an error", {
  # Only metadata indices, no items -> all skip; no votes, no agreement, but a
  # well-formed (empty) object.
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
  # On random data psychsyn / psychant find no qualifying pairs -> score nobody.
  # The standalone index warns; the screen reports it in $flags / print ("0 / 0"),
  # so NO warning leaks from cier_screen(). The muffle is targeted
  # (cier_warning_insufficient_items only), not suppressWarnings.
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

test_that("duplicate control names are a typed input error", {
  x <- screen_matrix(20L, 20L, 22L)
  dup_outer <- list(cier_irv = list(fpr = 0.10),
                    cier_irv = list(fpr = 0.20))
  expect_error(cier_screen(x, methods = "cier_irv", control = dup_outer),
               class = "cier_error_input")

  dup_inner <- list(cier_irv = list(fpr = 0.10, fpr = 0.20))
  expect_error(cier_screen(x, methods = "cier_irv", control = dup_inner),
               class = "cier_error_input")
})

test_that("an unknown or unnamed control argument is a typed input error", {
  # A typo'd / unknown argument name is caught precisely (not an opaque do.call
  # failure), and an unnamed inner argument (would splice positionally) is rejected.
  x <- screen_matrix(20L, 20L, 23L)
  expect_error(cier_screen(x, methods = "cier_irv",
                           control = list(cier_irv = list(fp = 0.1))),  # typo
               class = "cier_error_input")
  expect_error(cier_screen(x, methods = "cier_irv",
                           control = list(cier_irv = list(0.1))),       # unnamed
               class = "cier_error_input")
  # A control entry naming a screen-managed positional (responses / items) must
  # error, not override the screen's own responses or collide on items.
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
  # The agreement counts shown are sum(rowSums(votes) >= k), from the COLLAPSED
  # votes. A print recomputing from raw flags would differ at >= 2 (the consistency
  # pair). Pin the >= 1 and >= 2 counts to $agreement.
  obs_n <- round(sc$agreement$agreement$observed * n)
  expect_true(any(grepl(paste0("\\b", obs_n[[1L]], "\\b"), out)))
  expect_true(any(grepl(paste0("\\b", obs_n[[2L]], "\\b"), out)))
  # Sections are separated by blank lines (cli_verbatim would drop them; this pins
  # the separators actually render).
  expect_true(any(out == ""))
})

# ---- print snapshots (locked) -----------------------------------------------

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
  # observed > expected comparison would advertise sampling noise as contamination
  # (e.g. "10.0% vs expected 9.8% <- excess"). The marker must fire only when the
  # one-sided binomial tail is below 0.05.
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

test_that("a nonzero share that rounds below 0.1% prints as <0.1%, not 0.0%", {
  # The marker is gated on the unrounded tail, so a lone hit in a large sample stays
  # unmarked. If its expectation printed "0.0%" beside an observed "0.3%", the row
  # would self-contradict; "<0.1%" removes that, while a genuine (unreachable) zero
  # must still print "0.0%".
  ag <- list(agreement = data.frame(
    k        = 1:3,
    observed = c(1 / 300, 0, 0),
    expected = c(0.0004166667, 0, 7.2e-9)  # ~0.042%, an exact zero, ~7.2e-7%
  ))
  lines <- screen_agreement_lines(ag, 300L)
  # k = 1: observed rounds to 0.3%, expectation rounds below 0.1% -- the gap that
  # used to read "0.3% vs 0.0%" now reads "0.3% vs <0.1%", and is unmarked.
  expect_match(lines[[3L]], "(0.3%); expected <0.1%", fixed = TRUE)
  expect_false(grepl("expected 0.0%", lines[[3L]], fixed = TRUE))
  expect_false(grepl("excess", lines[[3L]], fixed = TRUE))
  # k = 2: a true zero expectation still prints "0.0%", never "<0.1%".
  expect_match(lines[[4L]], "expected 0.0%", fixed = TRUE)
  expect_false(grepl("<0.1%", lines[[4L]], fixed = TRUE))
  # k = 3: a vanishingly small but nonzero expectation prints "<0.1%".
  expect_match(lines[[5L]], "expected <0.1%", fixed = TRUE)
})

# ---- standalone-vs-unknown method error -------------------------------------

test_that("a standalone (registered, non-screenable) method says 'call directly'", {
  err <- tryCatch(
    cier_screen(screen_matrix(30L, 10L, 1L), methods = "cier_lazr"),
    error = function(e) e
  )
  expect_s3_class(err, "cier_error_input")
  msg <- gsub("\\s+", " ", rlang::cnd_message(err))
  expect_match(msg, "standalone")
  expect_match(msg, "cier_lazr")
})

test_that("a genuinely unknown method keeps the unknown-method wording", {
  err <- tryCatch(
    cier_screen(screen_matrix(30L, 10L, 1L), methods = "cier_not_real"),
    error = function(e) e
  )
  expect_s3_class(err, "cier_error_input")
  msg <- gsub("\\s+", " ", rlang::cnd_message(err))
  expect_match(msg, "[Uu]nknown")
  expect_match(msg, "Screenable methods")
  expect_match(msg, "cier_irv")
})

test_that("a mix of standalone and unknown ids lists both causes", {
  err <- tryCatch(
    cier_screen(screen_matrix(30L, 10L, 1L),
                methods = c("cier_lazr", "cier_not_real")),
    error = function(e) e
  )
  msg <- gsub("\\s+", " ", rlang::cnd_message(err))
  # Both causes are reported SEPARATELY: the standalone clause names cier_lazr, the
  # unknown clause names cier_not_real and lists screenable methods. A mutant lumping
  # both into one generic list would miss one of these distinct hints.
  expect_match(msg, "standalone")
  expect_match(msg, "Screenable methods")
  expect_match(msg, "cier_lazr")
  expect_match(msg, "cier_not_real")
})

# ---- missing reverse_keyed informs once at screen level ---------------------

test_that("cier_screen informs once (not per index) when items omits reverse_keyed", {
  x <- screen_matrix(80L, 20L, 7L)
  it <- data.frame(scale = rep(paste0("s", seq_len(4L)), each = 5L), max = 5L)
  counter <- new.env(parent = emptyenv())
  counter$n <- 0L
  withCallingHandlers(
    suppressWarnings(cier_screen(
      x, it,
      methods = c("cier_even_odd", "cier_personal_reliability"),
      control = list(cier_personal_reliability = list(seed = 1))
    )),
    cier_message_forward_keyed = function(m) {
      counter$n <- counter$n + 1L
      invokeRestart("muffleMessage")
    }
  )
  # Two keying indices run, but the inform fires once (per-index muffled, one
  # screen-level emission).
  expect_identical(counter$n, 1L)
})

test_that("cier_screen stays silent when reverse_keyed is declared", {
  x <- screen_matrix(80L, 20L, 7L)
  it <- screen_items(4L, 5L)            # carries reverse_keyed
  expect_no_message(
    suppressWarnings(cier_screen(
      x, it,
      methods = c("cier_even_odd", "cier_personal_reliability"),
      control = list(cier_personal_reliability = list(seed = 1))
    )),
    class = "cier_message_forward_keyed"
  )
})

# ---- battery-wide fpr -------------------------------------------------------

test_that("fpr forwards to the percentile family only, leaving the rest", {
  # The percentile-default indices pick up the battery fpr in their recorded
  # provenance rate; mahalanobis (chisq alpha), longstring (fixed fraction) and
  # gnormed (mc_null nominal) are on their own scales and hold.
  x <- clean_matrix()
  it <- clean_items()
  sc <- q(cier_screen(x, it, methods = base_methods(), fpr = 0.10,
                      control = list(cier_personal_reliability = list(seed = 1))))
  percentile7 <- c("cier_irv", "cier_even_odd", "cier_personal_reliability",
                   "cier_psychsyn", "cier_psychant", "cier_person_total")
  for (m in intersect(percentile7, names(sc$indices))) {
    expect_identical(sc$indices[[m]]$cutoff_method, "percentile")
    expect_identical(sc$indices[[m]]$cutoff_rate, 0.10)
  }
  # Excluded from fpr forwarding: method and rate unchanged from the default.
  expect_identical(sc$indices$cier_mahalanobis$cutoff_method, "chisq")
  expect_identical(sc$indices$cier_mahalanobis$cutoff_rate, 0.001)
  expect_identical(sc$indices$cier_longstring$cutoff_method, "fixed_fraction")
  expect_identical(sc$indices$cier_longstring$cutoff_rate, 0.5)
})

test_that("fpr leaves cier_ht (percentile) moving but cier_gnormed (mc_null) fixed", {
  x <- screen_matrix(120L, 20L, 13L)
  it <- screen_items()
  sc <- q(cier_screen(x, it, fpr = 0.10,
                      control = list(cier_personal_reliability = list(seed = 1),
                                     cier_gnormed = list(seed = 1))))
  expect_identical(sc$indices$cier_ht$cutoff_rate, 0.10)
  expect_identical(sc$indices$cier_gnormed$cutoff_method, "mc_null")
  expect_identical(sc$indices$cier_gnormed$cutoff_rate, 0.05)   # NOT 0.10
})

test_that("a per-index control fpr wins over the battery fpr", {
  x <- clean_matrix()
  it <- clean_items()
  sc <- q(cier_screen(x, it, methods = c("cier_irv", "cier_person_total"),
                      fpr = 0.10, control = list(cier_irv = list(fpr = 0.20))))
  expect_identical(sc$indices$cier_irv$cutoff_rate, 0.20)         # control wins
  expect_identical(sc$indices$cier_person_total$cutoff_rate, 0.10)  # battery fpr
})

test_that("a per-index control cutoff blocks the battery fpr splice (no error)", {
  # fpr and cutoff are mutually exclusive in the wrapper; splicing the battery fpr
  # onto an index already carrying a literal cutoff would error, so the splice must
  # skip it. The literal override stands.
  x <- clean_matrix()
  sc <- q(cier_screen(x, methods = "cier_irv", fpr = 0.10,
                      control = list(cier_irv = list(cutoff = 0.5))))
  expect_identical(sc$indices$cier_irv$cutoff_method, "literal")
  expect_identical(sc$indices$cier_irv$cutoff, 0.5)
})

test_that("the battery fpr actually moves the realised flag count", {
  x <- clean_matrix()
  lo <- q(cier_screen(x, methods = "cier_irv", fpr = 0.01))
  hi <- q(cier_screen(x, methods = "cier_irv", fpr = 0.10))
  expect_gt(sum(hi$indices$cier_irv$flagged, na.rm = TRUE),
            sum(lo$indices$cier_irv$flagged, na.rm = TRUE))
})

test_that("an invalid battery fpr is a typed input error", {
  x <- screen_matrix(40L, 20L, 1L)
  expect_error(cier_screen(x, fpr = 1.5), class = "cier_error_input")
  expect_error(cier_screen(x, fpr = 0), class = "cier_error_input")
  expect_error(cier_screen(x, fpr = c(0.05, 0.1)), class = "cier_error_input")
  expect_error(cier_screen(x, fpr = "x"), class = "cier_error_input")
})

# ---- no-pairs annotation ----------------------------------------------------

test_that("the screen captures the no-pairs reason in $notes", {
  # On random data psychsyn / psychant find no qualifying pairs at the default
  # critical_r = 0.6 -> score 0 / 0. The muffled actionable reason is captured in
  # $notes (method + a one-line explanation), still warning-free.
  x <- screen_matrix(40L, 20L, 30L)
  expect_no_warning(sc <- cier_screen(x))
  expect_s3_class(sc$notes, "data.frame")
  expect_setequal(sc$notes$method, c("cier_psychsyn", "cier_psychant"))
  syn <- sc$notes$note[sc$notes$method == "cier_psychsyn"]
  ant <- sc$notes$note[sc$notes$method == "cier_psychant"]
  expect_match(syn, "synonym")
  expect_match(syn, "critical_r = 0.6", fixed = TRUE)
  expect_match(syn, "cier_psychsyn_critval", fixed = TRUE)
  expect_match(ant, "antonym")
  expect_match(ant, "antonym = TRUE", fixed = TRUE)   # the tail-aware sweep hint
  # The strongest in-tail r is reported (a mutant omitting it or inserting the wrong
  # value must fail). Re-derive via the same pairing helpers the warning uses; the
  # VALUE travels through the captured payload into the note.
  cm <- pairing_cor(x)
  syn_strong <- format(round(strongest_pairing_cor(cm, antonym = FALSE), 3))
  ant_strong <- format(round(strongest_pairing_cor(cm, antonym = TRUE), 3))
  expect_match(syn, paste0("strongest in-tail r = ", syn_strong), fixed = TRUE)
  expect_match(ant, paste0("strongest in-tail r = ", ant_strong), fixed = TRUE)
})

test_that("$notes is empty when the pair indices find pairs (note-scoped to no-pairs)", {
  # Lower critical_r so synonym/antonym pairs qualify on the bundled BFI items ->
  # psychsyn/psychant score respondents, so NEITHER may appear in $notes. A mutant
  # noting whenever a pair index is present (regardless of pairs) must fail.
  nm <- names(bfi_careless)[1:44]
  resp <- bfi_careless[, 1:44]
  sc <- q(cier_screen(resp, methods = c("cier_psychsyn", "cier_psychant"),
                      control = list(cier_psychsyn = list(critical_r = 0.3),
                                     cier_psychant = list(critical_r = 0.3))))
  expect_false("cier_psychsyn" %in% sc$notes$method)
  expect_false("cier_psychant" %in% sc$notes$method)
})

test_that("$notes is empty when no pair index runs", {
  x <- screen_matrix(60L, 20L, 8L)
  sc <- q(cier_screen(x, methods = c("cier_irv", "cier_longstring")))
  expect_identical(nrow(sc$notes), 0L)
})

test_that("print marks ONLY the noted indices with * and lists the Notes section", {
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    out <- capture.output(q(print(cier_screen(screen_matrix(40L, 20L, 30L)))))
  })
  # The 0 / 0 psychsyn line carries a trailing marker, and a Notes: section
  # appears with the actionable per-index explanation.
  syn_line <- out[grepl("cier_psychsyn", out, fixed = TRUE) &
                    grepl("0 / 0", out, fixed = TRUE)]
  expect_true(any(grepl("*", syn_line, fixed = TRUE)))
  expect_true(any(grepl("Notes:", out, fixed = TRUE)))
  expect_true(any(grepl("synonym", out, fixed = TRUE)))
  # A scored, un-noted index line must NOT carry the marker (a mark-every-line
  # mutant survives the positive check alone).
  irv_line <- out[grepl("cier_irv ", out, fixed = TRUE) &
                    grepl("/", out, fixed = TRUE)]
  expect_false(any(grepl("*", irv_line, fixed = TRUE)))
})

# ---- summary delegates to print ---------------------------------------------

test_that("summary.cier_screen prints the same body as print and returns invisibly", {
  sc <- q(cier_screen(screen_matrix(60L, 20L, 9L)))
  withr::with_options(list(cli.width = 80, cli.unicode = FALSE), {
    p <- capture.output(print(sc))
    s <- capture.output(summary(sc))
  })
  expect_identical(s, p)
  expect_identical(withVisible(summary(sc))$visible, FALSE)
})
