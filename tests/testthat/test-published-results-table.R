# Fast-tier contracts for inst/extdata/published-results.csv -- the committed
# table behind the published-results vignette. No network, no recomputation:
# this file pins (1) the schema, (2) the transcribed paper values (digest +
# verbatim spot checks against the cited source tables), (3) the row
# inventory, (4) the status / empty-cell / explanation contracts, and (5) a
# few coarse direction floors that a sign or orientation flip cannot survive.
# The slow tier (test-published-results-reproduction.R) recomputes the cier
# cells from the fetched data.

pubres_statuses <- c("matches", "close", "differs", "abstains",
                     "not_implemented", "not_evaluated")

# Cells allowed to carry status "differs" -- each is a pre-registered,
# explained divergence (see the vignette prose and TOLERANCES.md). A new
# divergence anywhere else must be diagnosed and added here deliberately,
# with its explanation, before the suite passes. The cap bounds how far even
# an allowed divergence may go (recorded in TOLERANCES.md): it converts the
# allowed-differs pocket from "any value" to "the explained mechanism's
# plausible size", so a sign or orientation flip cannot hide inside it.
pubres_allowed_differs <- function(tab) {
  battery <- tab$paper == "goldammer_battery_2024"
  rpr <- tab$paper == "goldammer_rpr_2024"
  bru <- tab$paper == "bruhlmann_2020"
  out <- rep(NA_real_, nrow(tab))
  pick <- function(out, sel, width) ifelse(sel & is.na(out), width, out)
  out <- pick(out, battery & tab$index == "MD", 0.30)
  out <- pick(out, battery & tab$index == "Ht" &
                tab$study %in% c("4", "5"), 0.15)
  # whole-sample synonym pairing on Studies 4/5 finds only a few attenuated
  # pairs on the careless wave: scored, but well below the paper's
  # careful-derived pairing
  out <- pick(out, battery & tab$index == "Synonyms" &
                tab$study %in% c("4", "5"), 0.30)
  out <- pick(out, battery &
                tab$index %in% c("Antonyms (careful-reference)",
                                 "Synonyms (careful-reference)"), 0.10)
  out <- pick(out, rpr & tab$paper_statistic == "SEN95", 0.15)
  # the Bruhlmann consistency counts inherit the authors' correction of one
  # mis-keyed openness item (their analysis script; the vignette demonstrates
  # that replicating it reproduces their counts exactly)
  out <- pick(out, bru & tab$index %in% c("Odd-even consistency",
                                          "Resampled individual reliability",
                                          "Person-total correlation"), 30)
  out
}

test_that("published-results table has the exact schema", {
  tab <- pubres_table()
  expect_identical(names(tab),
                   c("paper", "study", "design", "level", "index",
                     "paper_statistic", "paper_value", "cier_value", "delta",
                     "status", "n", "note"))
  expect_identical(nrow(tab), 264L)
})

test_that("the paper-side transcription is frozen", {
  tab <- pubres_table()
  expect_identical(pubres_paper_digest(tab), "f08f6e196ca1d5919e08a43894c75881")
})

test_that("transcription spot checks match the cited source tables verbatim", {
  tab <- pubres_table()
  pv <- function(...) {
    row <- pubres_cell(tab, ...)
    expect_identical(nrow(row), 1L)
    row$paper_value
  }
  # Goldammer battery, supplement S5 / S6 / S14 / S15 / S42
  expect_identical(pv(c("study", "1"), c("level", "partial"), c("index", "RPR")), .86)
  expect_identical(pv(c("study", "1"), c("level", "partial"), c("index", "IRV")), .43)
  expect_identical(pv(c("study", "1"), c("level", "full"), c("index", "Time")), .87)
  expect_identical(pv(c("study", "2"), c("level", "partial"), c("index", "lz")), .81)
  expect_identical(pv(c("study", "2"), c("level", "full"), c("index", "r_pbis")), .86)
  expect_identical(pv(c("study", "4"), c("design", "WS"), c("level", "full"),
                      c("index", "Synonyms")), .94)
  expect_identical(pv(c("study", "5"), c("design", "BS"), c("level", "partial"),
                      c("index", "Ht")), .65)
  # Goldammer RPR paper, Table 9
  expect_identical(pv(c("paper", "goldammer_rpr_2024"), c("level", "6 facets"),
                      c("index", "RPR50"), c("paper_statistic", "AUC")), .693)
  expect_identical(pv(c("paper", "goldammer_rpr_2024"), c("level", "16 facets"),
                      c("index", "PR"), c("paper_statistic", "SEN95")), .744)
  expect_identical(pv(c("paper", "goldammer_rpr_2024"), c("level", "23 facets"),
                      c("index", "RPR100"), c("paper_statistic", "AUC")), .964)
  # the paper reports no SEN95 for PR on 6 facets (Table 9, footnote a)
  pr6 <- pubres_cell(tab, c("paper", "goldammer_rpr_2024"),
                     c("level", "6 facets"), c("index", "PR"),
                     c("paper_statistic", "SEN95"))
  expect_true(is.na(pr6$paper_value))
  expect_identical(pr6$status, "not_evaluated")
  # Schroeders et al., Table 3
  expect_identical(pv(c("paper", "schroeders_2022"), c("index", "Mahalanobis"),
                      c("paper_statistic", "sensitivity")), .18)
  expect_identical(pv(c("paper", "schroeders_2022"), c("index", "Zh"),
                      c("paper_statistic", "accuracy")), .67)
  expect_identical(pv(c("paper", "schroeders_2022"),
                      c("index", "GBM (responses + times)"),
                      c("paper_statistic", "balanced accuracy")), .66)
  # Bruhlmann et al., Table 2 / Section 4
  expect_identical(pv(c("paper", "bruhlmann_2020"), c("index", "Longstring")), 25)
  expect_identical(pv(c("paper", "bruhlmann_2020"),
                      c("index", "Person-total correlation")), 74)
  expect_identical(pv(c("paper", "bruhlmann_2020"),
                      c("index", "Latent profile class (careless)")), 181)
})

test_that("the row inventory is complete per paper", {
  tab <- pubres_table()
  papers <- c("bruhlmann_2020", "goldammer_battery_2024",
              "goldammer_rpr_2024", "schroeders_2022")
  expect_identical(as.vector(table(tab$paper)[papers]),
                   c(10L, 185L, 24L, 45L))

  bat <- tab[tab$paper == "goldammer_battery_2024" & tab$study != "3", ]
  cells <- unique(bat[, c("study", "design", "level")])
  expect_identical(nrow(cells), 12L)
  paper_idx <- c("RPR", "MD", "lz", "Synonyms", "Gnormed", "Outfit", "Infit",
                 "r_pbis", "Ht", "INDCHI", "Antonyms", "Time", "Longstring",
                 "IRV")
  for (i in seq_len(nrow(cells))) {
    here <- bat[bat$study == cells$study[i] & bat$design == cells$design[i] &
                  bat$level == cells$level[i], ]
    want <- c(paper_idx, if (cells$study[i] %in% c("4", "5")) {
      c("Synonyms (careful-reference)", "Antonyms (careful-reference)")
    })
    expect_setequal(here$index, want)
    expect_identical(anyDuplicated(here$index), 0L)
  }
  s3 <- tab[tab$paper == "goldammer_battery_2024" & tab$study == "3", ]
  expect_identical(nrow(s3), 1L)
  expect_identical(s3$status, "not_evaluated")

  rpr <- tab[tab$paper == "goldammer_rpr_2024", ]
  expect_setequal(rpr$level, c("6 facets", "16 facets", "23 facets"))
  expect_setequal(rpr$index, c("PR", "RPR25", "RPR50", "RPR100"))
  expect_setequal(rpr$paper_statistic, c("AUC", "SEN95"))
  expect_identical(nrow(unique(rpr[, c("level", "index", "paper_statistic")])), 24L)

  sch <- tab[tab$paper == "schroeders_2022", ]
  expect_setequal(sch$index, c("Mahalanobis", "Antonyms", "Even-odd",
                               "Longstring", "IRV", "Zh", "GBM (responses)",
                               "GBM (response times)",
                               "GBM (responses + times)"))
  expect_setequal(sch$paper_statistic, c("accuracy", "sensitivity",
                                         "specificity", "precision",
                                         "balanced accuracy"))
  expect_identical(nrow(unique(sch[, c("index", "paper_statistic")])), 45L)

  bru <- tab[tab$paper == "bruhlmann_2020", ]
  bru_idx <- c("Self-report aggregate", "Bogus item",
               "Instructed response item", "Response time", "Longstring",
               "Odd-even consistency", "Resampled individual reliability",
               "Person-total correlation", "Open-answer quality",
               "Latent profile class (careless)")
  expect_setequal(bru$index, bru_idx)
})

test_that("the not-implemented set is pinned exactly", {
  tab <- pubres_table()
  ni <- tab[tab$status == "not_implemented", ]
  expect_setequal(unique(ni$index[ni$paper == "goldammer_battery_2024"]),
                  c("lz", "Infit", "Outfit", "INDCHI"))
  expect_identical(sum(ni$paper == "goldammer_battery_2024"), 48L)
  expect_setequal(unique(ni$index[ni$paper == "schroeders_2022"]),
                  c("Zh", "GBM (responses)", "GBM (response times)",
                    "GBM (responses + times)"))
  # the response-time row stays empty for a different reason than the rest:
  # the authors flagged nobody by visual inspection and published no numeric
  # cutoff, so there is no rule to reproduce (cier_total_time itself shipped)
  expect_setequal(ni$index[ni$paper == "bruhlmann_2020"],
                  c("Self-report aggregate", "Response time",
                    "Open-answer quality",
                    "Latent profile class (careless)"))
  expect_identical(sum(ni$paper == "goldammer_rpr_2024"), 0L)
})

test_that("status, empty-cell, delta, and explanation contracts hold", {
  tab <- pubres_table()
  expect_true(all(nzchar(tab$status)))
  expect_true(all(tab$status %in% pubres_statuses))

  compared <- tab$status %in% c("matches", "close", "differs")
  expect_true(all(!is.na(tab$paper_value[compared])))
  expect_true(all(!is.na(tab$cier_value[compared])))

  empty_cier <- tab$status %in% c("abstains", "not_implemented")
  expect_true(all(is.na(tab$cier_value[empty_cier])))
  expect_true(all(!is.na(tab$paper_value[empty_cier])))
  expect_true(all(is.na(tab$paper_value[tab$status == "not_evaluated"])))
  expect_true(all(tab$status[is.na(tab$cier_value)] %in%
                    c("abstains", "not_implemented", "not_evaluated")))

  # the abstention inventory is pinned exactly, from the measured data: on
  # Studies 1/2 the single contaminated wave yields no pairs at
  # critical_r = .60 for either pairing index; on Studies 4/5 the careless
  # second wave yields no ANTONYM pairs (synonyms keep a handful and score,
  # as capped differs cells). Nothing else may abstain -- a silently relaxed
  # critical_r must surface here as a deliberate change.
  abst <- tab[tab$status == "abstains", ]
  bat <- tab[tab$paper == "goldammer_battery_2024" & tab$study != "3", ]
  cells <- unique(bat[, c("study", "design", "level")])
  want12 <- merge(cells[cells$study %in% c("1", "2"), ],
                  data.frame(index = c("Synonyms", "Antonyms")))
  want45 <- merge(cells[cells$study %in% c("4", "5"), ],
                  data.frame(index = "Antonyms"))
  want <- rbind(want12, want45)
  key <- function(d) paste(d$study, d$design, d$level, d$index)
  expect_setequal(key(abst), key(want))

  # delta is exactly round(cier - paper, 3) where both sides exist
  both <- !is.na(tab$paper_value) & !is.na(tab$cier_value)
  expect_equal(tab$delta[both],
               round(tab$cier_value[both] - tab$paper_value[both], 3))
  expect_true(all(is.na(tab$delta[!both])))

  # the status vocabulary is recomputable from the stored values
  counts <- tab$paper == "bruhlmann_2020"
  twodp <- both & !counts &
    round(tab$cier_value, 2) == round(tab$paper_value, 2)
  inband <- both & !counts & abs(tab$delta) <= 0.03 + 1e-9
  expect_identical(tab$status[both & !counts],
                   as.character(ifelse(twodp[both & !counts], "matches",
                                       ifelse(inband[both & !counts], "close",
                                              "differs"))))
  bc <- both & counts
  expect_identical(tab$status[bc],
                   as.character(ifelse(tab$delta[bc] == 0, "matches",
                                       ifelse(abs(tab$delta[bc]) <= 2, "close",
                                              "differs"))))

  # not_evaluated means the paper has no value -- cier shows none either
  # (any cier-side observation lives in the note, not the value column)
  expect_true(all(is.na(tab$cier_value[tab$status == "not_evaluated"])))

  # every divergence and every abstention carries a written explanation,
  # and divergences appear only in the pre-registered allowed list, bounded
  # by each mechanism's cap
  needs_note <- tab$status %in% c("differs", "abstains", "not_implemented",
                                  "not_evaluated")
  expect_true(all(!is.na(tab$note[needs_note]) & nzchar(tab$note[needs_note])))
  caps <- pubres_allowed_differs(tab)
  differs <- which(tab$status == "differs")
  expect_true(all(!is.na(caps[differs])))
  expect_true(all(abs(tab$delta[differs]) <= caps[differs]))
})

test_that("values are in range and the reproduction is not vacuous", {
  tab <- pubres_table()
  metric <- tab$paper != "bruhlmann_2020"
  for (col in c("paper_value", "cier_value")) {
    x <- tab[[col]][metric]
    expect_true(all(x >= 0 & x <= 1, na.rm = TRUE))
    y <- tab[[col]][!metric]
    expect_true(all(y >= 0 & y <= 394 & y == round(y), na.rm = TRUE))
  }
  expect_true(all(tab$n > 0, na.rm = TRUE))

  # coarse direction floors, one per family with an allowed-differs pocket
  # included: an orientation flip, label mis-slicing, or wrong-group
  # construction cannot survive these
  floor_cell <- function(lo, ...) {
    row <- pubres_cell(tab, ...)
    expect_identical(nrow(row), 1L)
    expect_gte(row$cier_value, lo)
  }
  floor_cell(0.75, c("study", "4"), c("design", "BS"), c("level", "full"),
             c("index", "r_pbis"))
  floor_cell(0.75, c("study", "1"), c("design", "BS"), c("level", "full"),
             c("index", "RPR"))
  # total time: low totals are careless, so a missed orientation flip on the
  # lower-direction score cannot clear this floor (paper: .87)
  floor_cell(0.75, c("study", "1"), c("design", "BS"), c("level", "full"),
             c("index", "Time"))
  # facet-level RPR splits: domain-level splits depress this cell to ~.83
  floor_cell(0.84, c("study", "1"), c("design", "BS"), c("level", "partial"),
             c("index", "RPR"))
  floor_cell(0.70, c("study", "1"), c("design", "BS"), c("level", "partial"),
             c("index", "MD"))
  floor_cell(0.75, c("study", "4"), c("design", "BS"), c("level", "full"),
             c("index", "Synonyms (careful-reference)"))
  floor_cell(0.75, c("paper", "goldammer_rpr_2024"), c("level", "23 facets"),
             c("index", "RPR25"), c("paper_statistic", "AUC"))
  floor_cell(0.50, c("paper", "goldammer_rpr_2024"), c("level", "23 facets"),
             c("index", "RPR50"), c("paper_statistic", "SEN95"))
  floor_cell(0.50, c("paper", "schroeders_2022"), c("index", "Mahalanobis"),
             c("paper_statistic", "specificity"))
  for (ix in c("Resampled individual reliability",
               "Person-total correlation")) {
    cnt <- pubres_cell(tab, c("paper", "bruhlmann_2020"), c("index", ix))
    expect_true(cnt$cier_value >= 40 && cnt$cier_value <= 120)
  }

  # each paper block reproduces something at the paper's own precision
  for (pp in unique(tab$paper)) {
    expect_true(any(tab$status[tab$paper == pp] == "matches"),
                label = paste("at least one matches row for", pp))
  }
})

test_that("the shared statistics pin their definitions on fixtures", {
  # rank AUC: orientation, value, ties, and NA handling
  expect_identical(pubres_auc(c(1, 2, 3, 4), c(0L, 0L, 1L, 1L)), 1)
  expect_identical(pubres_auc(c(1, 2, 3, 4), c(1L, 1L, 0L, 0L)), 0)
  expect_identical(pubres_auc(c(1, 4, 2, 3), c(0L, 1L, 1L, 0L)), 0.75)
  expect_identical(pubres_auc(c(1, 1), c(0L, 1L)), 0.5)
  expect_true(is.na(pubres_auc(c(NA, 5), c(1L, 0L))))

  # sen95: the cutoff is the CAREFUL group's 95th percentile, exceeded
  # STRICTLY. careful = 1..20 makes the type-7 cutoff exactly 19.05; a
  # careless score equal to it must not count (>= would give 1 here).
  s <- c(1:20, 19.05, 25)
  l <- c(rep(0L, 20), 1L, 1L)
  expect_identical(pubres_sen95(s, l), 0.5)
  # and the quantile must come from the careful group: with all careless
  # scores far above every careful score, sensitivity is exactly 1 (a
  # careless-group cutoff would drop it to 0.25)
  s2 <- c(1:20, 25, 26, 27, 28)
  l2 <- c(rep(0L, 20), rep(1L, 4))
  expect_identical(pubres_sen95(s2, l2), 1)

  # the Schroeders binarisation boundaries (authors' OSF rules)
  expect_identical(pubres_flag_longstring(c(5, 6, 7)), c(FALSE, TRUE, TRUE))
  expect_identical(pubres_flag_antonyms(c(NA, -0.2, 0.001, 0)),
                   c(FALSE, FALSE, TRUE, FALSE))
  expect_identical(pubres_flag_evenodd(c(NA, -0.1, 0.1, 0)),
                   c(FALSE, FALSE, TRUE, FALSE))
  # IRV flags the LOWER tail (straightliners), AT the draw's 10th percentile
  # (<=): with two 1s and 2..10 the type-7 cutoff is exactly 1, so both 1s
  # flag (strict < would flag none; the upper-tail mutant flags neither)
  v <- c(1, 1, 2:10)
  expect_identical(pubres_flag_irv(v), c(TRUE, TRUE, rep(FALSE, 9)))
})

test_that("the shipped vignette displays the committed values (anti-drift)", {
  tab <- pubres_table()
  rmd <- file.path(testthat::test_path("..", ".."), "vignettes",
                   "published-results.Rmd")
  skip_if(!file.exists(rmd), "vignette source not present")
  text <- paste(readLines(rmd, warn = FALSE), collapse = "\n")
  sentinel <- function(...) {
    row <- pubres_cell(tab, ...)
    expect_identical(nrow(row), 1L)
    skip_if(is.na(row$cier_value), "sentinel cell not filled yet")
    expect_match(text, sub("^0", "", sprintf("%.3f", row$cier_value)),
                 fixed = TRUE)
  }
  sentinel(c("study", "1"), c("design", "BS"), c("level", "full"),
           c("index", "RPR"))
  sentinel(c("study", "1"), c("design", "BS"), c("level", "full"),
           c("index", "Time"))
  sentinel(c("paper", "goldammer_rpr_2024"), c("level", "23 facets"),
           c("index", "RPR25"), c("paper_statistic", "AUC"))
  sentinel(c("paper", "schroeders_2022"), c("index", "Mahalanobis"),
           c("paper_statistic", "sensitivity"))
  bru <- pubres_cell(tab, c("paper", "bruhlmann_2020"), c("index", "Longstring"))
  skip_if(is.na(bru$cier_value), "sentinel cell not filled yet")
  expect_match(text, as.character(bru$cier_value), fixed = TRUE)
  bog <- pubres_cell(tab, c("paper", "bruhlmann_2020"), c("index", "Bogus item"))
  skip_if(is.na(bog$cier_value), "sentinel cell not filled yet")
  expect_match(text, as.character(bog$cier_value), fixed = TRUE)
})
