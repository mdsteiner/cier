# Slow-tier reproduction of inst/extdata/published-results.csv: re-fetch the
# external datasets (skip cleanly when offline or when an upstream file has
# changed -- see the integrity pins), recompute every cier cell with the same
# conventions and seeds the vignette uses, and hold the committed values to
# 5e-4 (3-dp storage; all resampling is seeded). Opt in with
# CIER_SLOW_TESTS=true. The Bruhlmann block needs no network (bundled data).

# Compare recomputed cells against the committed table: equal values at the
# recorded tolerance and an identical NA (abstention) pattern.
pubres_compare <- function(committed, recomputed, by, value_col = "cier_value") {
  joined <- merge(committed, recomputed, by = by, suffixes = c("_csv", "_new"))
  expect_identical(nrow(joined), nrow(committed))
  csv <- joined[[paste0(value_col, "_csv")]]
  new <- joined[[paste0(value_col, "_new")]]
  expect_identical(is.na(csv), is.na(new))
  expect_equal(new[!is.na(new)], csv[!is.na(csv)], tolerance = 5e-4)
}

test_that("Bruhlmann flagged counts reproduce from the bundled data", {
  skip_on_cran()
  skip_if_slow()
  tab <- pubres_table()
  committed <- tab[tab$paper == "bruhlmann_2020" &
                     !tab$status %in% c("not_implemented", "not_evaluated"), ]
  expect_identical(nrow(committed), 4L)
  pubres_compare(committed, pubres_bruhlmann_cells(), by = "index")
})

test_that("Goldammer battery AUCs reproduce from the fetched data", {
  skip_on_cran()
  skip_if_slow()
  zip <- pubres_fetch(pubres_urls$battery_zip, "goldammer-battery.zip")
  if (is.null(zip)) skip("Goldammer battery data not fetchable")
  pubres_pin(unname(tools::md5sum(zip)), "5f60c24e1ff48af0b3a802389af5a10b",
             "battery zip md5")
  dir <- pubres_battery_dir()

  s1 <- pubres_battery_read(dir, "Study 1/Study_1.csv")
  pubres_pin(dim(s1), c(357L, 526L), "Study 1 dimensions")
  pubres_pin(as.vector(table(s1$careless_100)), c(119L, 118L),
             "Study 1 careless_100")
  pubres_pin(as.vector(table(s1$careless_33)), c(119L, 120L),
             "Study 1 careless_33")
  item <- grep("^[eacno]_[a-z]+[0-9]+$", names(s1), value = TRUE)
  pubres_pin(length(item), 60L, "Study 1 item count")
  pubres_pin(as.numeric(range(as.matrix(s1[, item]), na.rm = TRUE)), c(1, 6),
             "Study 1 response range")
  twins <- item[vapply(item, function(x) paste0(x, "R") %in% names(s1),
                       logical(1))]
  pubres_pin(length(twins), 30L, "Study 1 reverse-keyed twins")
  pubres_pin(all(s1[[paste0(twins[1], "R")]] == 7 - s1[[twins[1]]],
                 na.rm = TRUE), TRUE, "Study 1 twin coding (7 - x)")
  s2 <- pubres_battery_read(dir, "Study 2/Study_2.csv")
  pubres_pin(dim(s2), c(341L, 480L), "Study 2 dimensions")
  pubres_pin(as.vector(table(s2$careless_100)), c(112L, 115L),
             "Study 2 careless_100")
  s4 <- pubres_battery_read(dir, "Study 4/Study_4_wide.csv")
  pubres_pin(dim(s4), c(481L, 823L), "Study 4 dimensions")
  pubres_pin(as.vector(table(s4$condition)), c(120L, 121L, 122L, 118L),
             "Study 4 conditions")
  s5 <- pubres_battery_read(dir, "Study 5/Study_5_wide.csv")
  pubres_pin(dim(s5), c(481L, 823L), "Study 5 dimensions")
  pubres_pin(as.vector(table(s5$condition)), c(120L, 120L, 121L, 120L),
             "Study 5 conditions")

  tab <- pubres_table()
  committed <- tab[tab$paper == "goldammer_battery_2024" & tab$study != "3" &
                     !tab$status %in% c("not_implemented", "not_evaluated"), ]
  expect_identical(nrow(committed), 124L)  # (7 core + syn + ant) * 12 + 16
  pubres_compare(committed, pubres_battery_cells(dir),
                 by = c("study", "design", "level", "index"))
})

test_that("Goldammer RPR-paper Table 9 reproduces from the fetched data", {
  skip_on_cran()
  skip_if_slow()
  zip <- pubres_fetch(pubres_urls$rpr_zip, "goldammer-rpr.zip")
  if (is.null(zip)) skip("Goldammer RPR data not fetchable")
  pubres_pin(unname(tools::md5sum(zip)), "ad11b6e9a519ad42a5d5d6623249c42c",
             "RPR zip md5")
  d <- pubres_rpr_data()
  pubres_pin(dim(d), c(359L, 393L), "RPR real-data dimensions")
  pubres_pin(as.vector(table(d$carel_grp)), c(121L, 238L), "carel_grp")
  pubres_pin(all(d$oca2_R == 6 - d$oca2), TRUE, "RPR twin coding (6 - x)")
  pubres_pin(sum(is.na(d$PR_F6)), 1L, "stored PR_F6 NA count")
  its <- unlist(pubres_rpr_facets(), use.names = FALSE)
  pubres_pin(all(its %in% names(d)), TRUE, "facet item columns present")
  pubres_pin(anyNA(d[, its]), FALSE, "facet item completeness")

  cells <- pubres_rpr_cells(d)
  tab <- pubres_table()
  committed <- tab[tab$paper == "goldammer_rpr_2024" &
                     tab$status != "not_evaluated", ]
  expect_identical(nrow(committed), 23L)
  long_auc <- data.frame(level = cells$level, index = cells$index,
                         paper_statistic = "AUC", cier_value = cells$AUC,
                         stringsAsFactors = FALSE)
  long_sen <- data.frame(level = cells$level, index = cells$index,
                         paper_statistic = "SEN95", cier_value = cells$SEN95,
                         stringsAsFactors = FALSE)
  long <- rbind(long_auc, long_sen)
  pubres_compare(committed, long, by = c("level", "index", "paper_statistic"))
})

test_that("Schroeders Table 3 metrics reproduce from the fetched data", {
  skip_on_cran()
  skip_if_slow()
  csv <- pubres_fetch(pubres_urls$schroeders_resp,
                      "schroeders-data_mod_resp.csv")
  if (is.null(csv)) skip("Schroeders data not fetchable")
  pubres_pin(unname(tools::md5sum(csv)), "6a58c297ec4814c39525f424e42e03d5",
             "Schroeders data_mod_resp md5")
  d <- pubres_schroeders_data()
  pubres_pin(nrow(d), 605L, "Schroeders row count")
  pubres_pin(as.vector(table(d$Careless)), c(361L, 244L), "Careless label")
  he <- grep("^HE01_", names(d), value = TRUE)
  pubres_pin(length(he), 60L, "HEXACO item count")
  pubres_pin(as.numeric(range(as.matrix(d[, he]))), c(1, 5), "response range")
  pubres_pin(anyNA(d[, he]), FALSE, "item completeness")
  pubres_pin(all(pubres_schroeders_eo_items() %in% he), TRUE,
             "even-odd item naming")

  cells <- pubres_schroeders_cells(d)
  tab <- pubres_table()
  committed <- tab[tab$paper == "schroeders_2022" &
                     tab$status != "not_implemented", ]
  expect_identical(nrow(committed), 25L)
  pubres_compare(committed, cells, by = c("index", "paper_statistic"))
})
