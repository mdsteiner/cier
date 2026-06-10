# Method-properties registry: loader, accessors, and the validator that
# guards the hand-maintained CSV. (cier_methods() is exported: it is the
# accessor the unknown-method error hints point users to.)

registry_schema_columns <- function() {
  c("method", "family", "paper_year", "paper_citation_key", "doi",
    "default_cutoff_method", "default_cutoff_value", "flag_direction",
    "companion_methods", "backend", "screenable", "vote_group", "notes")
}

# The fifteen rows with their behaviourally-binding fields (the spec): the ten
# v0 indices, then the five v0.2 additions (autocorrelation, lazr, total / page
# time, attention). The v0.2 rows ship screenable = FALSE -- the screen wiring is
# deferred until the study has measured redundancy -- so they register metadata
# without entering cier_screen()'s run set.
expected_registry <- function() {
  data.frame(
    method = c(
      "cier_longstring", "cier_irv", "cier_even_odd",
      "cier_personal_reliability", "cier_psychsyn", "cier_psychant",
      "cier_mahalanobis", "cier_person_total", "cier_gnormed", "cier_ht",
      "cier_autocorrelation", "cier_lazr", "cier_total_time",
      "cier_page_time", "cier_attention"
    ),
    family = c(rep("indirect", 8L), "personfit", "personfit",
               "indirect", "indirect", "timing", "timing", "direct"),
    default_cutoff_method = c(
      "fixed", "percentile", "percentile", "percentile", "percentile",
      "percentile", "chisq", "percentile", "perfit_null", "percentile",
      "percentile", "percentile", "percentile", "fixed", "fixed"
    ),
    # Percentile rows store the false-positive tail mass (fpr) uniformly; the
    # flag_direction picks the tail and resolve_cutoff applies the single flip.
    # Gnormed uses the PerFit Monte-Carlo null (perfit_null); its value is the
    # nominal level (Blvl). page_time / attention use a fixed cited count (1).
    default_cutoff_value = c(0.5, 0.05, 0.05, 0.05, 0.05, 0.05, 0.001, 0.05,
                             0.05, 0.05, 0.05, 0.05, 0.05, 1, 1),
    flag_direction = c("upper", "lower", "upper", "upper", "lower", "upper",
                       "upper", "lower", "upper", "lower",
                       "upper", "upper", "lower", "upper", "upper"),
    backend = c(rep(NA_character_, 8L), "PerFit", "mokken",
                rep(NA_character_, 5L)),
    # even-odd and personal_reliability share the `consistency` vote so they
    # collapse to ONE vote in cier_screen(); every other index is its own vote.
    # The v0.2 rows keep own-id votes (the "repetition" grouping is deferred).
    vote_group = c("cier_longstring", "cier_irv", "consistency", "consistency",
                   "cier_psychsyn", "cier_psychant", "cier_mahalanobis",
                   "cier_person_total", "cier_gnormed", "cier_ht",
                   "cier_autocorrelation", "cier_lazr", "cier_total_time",
                   "cier_page_time", "cier_attention"),
    screenable = c(rep(TRUE, 10L), rep(FALSE, 5L)),
    stringsAsFactors = FALSE
  )
}

test_that("the registry loads as a 15-row, 13-column object and caches", {
  reg <- cier_methods()
  expect_s3_class(reg, "cier_method_info")
  expect_identical(nrow(reg), 15L)
  expect_identical(names(reg), registry_schema_columns())
  expect_setequal(reg$method, expected_registry()$method)
  expect_identical(cier_methods(), reg)
})

test_that("cier_methods is exported (the unknown-method error hints point to it)", {
  # Two user-facing aborts say "See cier_methods() for the available set"; the
  # hint must resolve for a user, so the accessor stays exported.
  expect_true("cier_methods" %in% getNamespaceExports("cier"))
})

test_that("each row's behaviourally-binding fields match the spec", {
  reg <- cier_methods()
  exp <- expected_registry()
  reg <- reg[match(exp$method, reg$method), , drop = FALSE]
  expect_identical(reg$family, exp$family)
  expect_identical(reg$default_cutoff_method, exp$default_cutoff_method)
  expect_equal(reg$default_cutoff_value, exp$default_cutoff_value)
  expect_identical(reg$flag_direction, exp$flag_direction)
  expect_identical(reg$backend, exp$backend)
  expect_identical(reg$vote_group, exp$vote_group)
  expect_identical(reg$screenable, exp$screenable)
})

test_that("the v0.2 additions ship screenable = FALSE (screen wiring deferred)", {
  reg <- cier_methods()
  v02 <- c("cier_autocorrelation", "cier_lazr", "cier_total_time",
           "cier_page_time", "cier_attention")
  # The five new indices register metadata but stay OUT of the screen's run set
  # until the study has measured redundancy; the ten v0 indices remain screenable.
  expect_false(any(reg$screenable[reg$method %in% v02]))
  expect_true(all(reg$screenable[!reg$method %in% v02]))
})

test_that("the consistency construct collapses even-odd + personal_reliability", {
  reg <- cier_methods()
  vg <- reg$vote_group[match(c("cier_even_odd", "cier_personal_reliability"),
                             reg$method)]
  expect_identical(vg, c("consistency", "consistency"))
  # Every other index keeps a distinct vote (its own id), so only this pair fuses.
  singletons <- reg$vote_group[!reg$vote_group %in% "consistency"]
  expect_identical(anyDuplicated(singletons), 0L)
})

test_that("person-fit backends are PerFit (Gnormed) and mokken (Ht)", {
  reg <- cier_methods()
  # Only Gnormed (PerFit) and Ht (mokken) carry a backend; the other 13 rows
  # (8 v0 indirect + 5 v0.2 additions) are base-R and have NA backend.
  expect_length(reg$method[is.na(reg$backend)], 13L)
  expect_identical(reg$backend[reg$method == "cier_gnormed"], "PerFit")
  expect_identical(reg$backend[reg$method == "cier_ht"], "mokken")
})

test_that("personal_reliability cites the resampled-PR paper", {
  row <- cier_methods()[cier_methods()$method == "cier_personal_reliability", ]
  expect_identical(row$paper_citation_key, "goldammer2024rpr")
  expect_identical(row$doi, "10.3758/s13428-024-02506-0")
})

test_that("cier_method_row returns one row and rejects bad input", {
  expect_identical(nrow(cier_method_row("cier_longstring")), 1L)
  expect_identical(cier_method_row("cier_gnormed")$family, "personfit")
  expect_error(cier_method_row("not_a_method"), class = "cier_error_input")
  expect_error(cier_method_row(123L), class = "cier_error_input")
})

test_that("the validator catches realistic registry-edit mistakes", {
  wrong_cols <- new_cier_method_info(cbind(as.data.frame(cier_methods()), x = 1))
  expect_error(validate_cier_method_info(wrong_cols), class = "cier_error_data")

  bad_family <- as.data.frame(cier_methods())
  bad_family$family[1] <- "model_based"   # outside the family vocabulary
  expect_error(validate_cier_method_info(new_cier_method_info(bad_family)),
               class = "cier_error_data")

  dup <- as.data.frame(cier_methods())
  dup$method[2] <- dup$method[1]
  expect_error(validate_cier_method_info(new_cier_method_info(dup)),
               class = "cier_error_data")

  bad_comp <- as.data.frame(cier_methods())
  bad_comp$companion_methods[1] <- "cier_not_real"
  expect_error(validate_cier_method_info(new_cier_method_info(bad_comp)),
               class = "cier_error_data")

  na_value <- as.data.frame(cier_methods())
  na_value$default_cutoff_value[1] <- NA_real_
  expect_error(validate_cier_method_info(new_cier_method_info(na_value)),
               class = "cier_error_data")

  # A missing vote_group would silently drop an index from the screen's vote
  # collapse, so it is a registry-data error (not a coercion NA).
  na_group <- as.data.frame(cier_methods())
  na_group$vote_group[1] <- NA_character_
  expect_error(validate_cier_method_info(new_cier_method_info(na_group)),
               class = "cier_error_data")

  blank_group <- as.data.frame(cier_methods())
  blank_group$vote_group[1] <- ""
  expect_error(validate_cier_method_info(new_cier_method_info(blank_group)),
               class = "cier_error_data")
})
