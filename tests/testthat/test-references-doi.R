# Guard: every DOI cited in a man-page \references section must resolve to
# a DOI carried by the method-properties registry. This keeps the
# documented sources and the registry's single source of citation truth in
# lock-step as index wrappers add @references. (DOI-less references, e.g.
# classic pre-DOI papers, are not constrained.)
#
# Data-documentation pages (\docType{data}) are exempt: they cite the source of
# a bundled dataset, which is a data provenance reference, not a method citation,
# and so is recorded in LICENSE.note / the data page rather than the registry.
#
# The simulator's page (cier_simulate.Rd) is a special case: cier_simulate() is a
# data generator, not an index -- by the ADR it has NO registry row (no cutoff, no
# direction, never in cier_screen()), so its citations (the published simulation
# designs and calibration anchors it follows) cannot be carried by the registry.
# Rather than exempt the page from the guard wholesale -- which once let a
# fabricated Goldammer citation through on an otherwise-correct DOI -- we pin its
# DOIs to an explicit, self-contained allow-list and its title/author text to a
# frozen citation core (below): any DOI outside the allow-list, or any drift in the
# pinned citation text, fails.
generator_doc_file <- function() {
  "cier_simulate.Rd"
}

# Exhaustive allow-list of DOIs cier_simulate.Rd may cite (cited-on-page must be a
# subset). Independent of the registry, so it stays correct as the registry
# changes; authors named for the reviewer.
generator_citation_dois <- function() {
  c(
    "10.1177/10944281211056520",   # Bowling et al. (2023)
    "10.1016/j.metip.2020.100022", # Bruhlmann et al. (2020)
    "10.1016/j.jesp.2015.07.006",  # Curran (2016)
    "10.3758/s13428-024-02484-3",  # Goldammer et al. (2024), detection revisited
    "10.1177/00131644211004708",   # Schroeders et al. (2022)
    "10.3758/s13428-023-02074-9",  # Voss (2024)
    "10.48550/arXiv.2303.07167"    # Welz & Alfons (2023)
  )
}

# Frozen author + year + title cores of citations that recur across man pages, so
# a fabricated title or a dropped/reordered author fails without any online lookup.
# The core stops before the journal tail, so per-page volume/DOI formatting does
# not make the fixture brittle. The o-umlaut in "Stockli" is a \u escape so this
# test source stays encoding-independent (man pages are read as UTF-8).
canonical_reference_cores <- function() {
  list(
    goldammer2024_revisited = paste0(
      "Goldammer, P., St\u00f6ckli, P. L., Escher, Y. A., Annen, H., Jonas, K., ",
      "& Antonakis, J. (2024). Careless responding detection revisited: ",
      "Accuracy of direct and indirect measures."
    )
  )
}

# Man page -> the citation-core keys it must contain verbatim.
citation_page_keys <- function() {
  list(
    "cier_screen.Rd" = "goldammer2024_revisited",
    "cier_simulate.Rd" = "goldammer2024_revisited"
  )
}

# TRUE when an Rd file documents a dataset (cites a data source, not a method).
is_data_doc <- function(rd_file) {
  lines <- readLines(rd_file, warn = FALSE, encoding = "UTF-8")
  any(grepl("\\\\docType\\{data\\}", lines))
}

# Extract the lines of every \references{...} block from one Rd file's
# lines, brace-matched so multi-line blocks are captured in full.
extract_references_blocks <- function(lines) {
  starts <- grep("\\\\references\\{", lines)
  out <- character(0)
  for (s in starts) {
    depth <- 0L
    for (i in seq.int(s, length(lines))) {
      ln <- lines[[i]]
      out <- c(out, ln)
      opens  <- lengths(regmatches(ln, gregexpr("\\{", ln)))
      closes <- lengths(regmatches(ln, gregexpr("\\}", ln)))
      depth <- depth + opens - closes
      if (depth <= 0L) break
    }
  }
  out
}

dois_in_references <- function(rd_file) {
  lines <- readLines(rd_file, warn = FALSE, encoding = "UTF-8")
  block <- extract_references_blocks(lines)
  if (length(block) == 0L) return(character(0))
  txt <- paste(block, collapse = "\n")
  m <- regmatches(txt, gregexpr("10\\.\\d{4,9}/[^\\s\"'<>})]+", txt,
                                perl = TRUE))[[1]]
  unique(gsub("[.,;]+$", "", m))
}

# DOIs cited on the generator page that are NOT in the allow-list (should be none).
unlisted_generator_dois <- function(rd_file) {
  setdiff(dois_in_references(rd_file), generator_citation_dois())
}

# A man page's \references block flattened to one whitespace-normalised line, so a
# frozen citation core matches regardless of how roxygen wrapped the text.
normalize_ref_block <- function(rd_file) {
  lines <- readLines(rd_file, warn = FALSE, encoding = "UTF-8")
  block <- extract_references_blocks(lines)
  trimws(gsub("\\s+", " ", paste(block, collapse = " ")))
}

test_that("every @references DOI resolves to a registry citation entry", {
  skip_on_cran()
  skip_if_no_description()
  man_dir <- file.path(cier_pkg_root(), "man")
  skip_if(!dir.exists(man_dir), "man/ not on disk")
  rd <- list.files(man_dir, pattern = "\\.Rd$", full.names = TRUE)
  rd <- rd[basename(rd) != generator_doc_file()]
  rd <- Filter(function(f) !is_data_doc(f), rd)
  skip_if(length(rd) == 0L, "no man pages on disk")

  cited <- unique(unlist(lapply(rd, dois_in_references)))
  reg_dois <- unique(stats::na.omit(cier_methods()$doi))
  unresolved <- setdiff(cited, reg_dois)
  if (length(unresolved) > 0L) {
    cat("\nDOIs cited in man/ but absent from the registry:\n")
    cat(unresolved, sep = "\n")
    cat("\n")
  }
  expect_identical(unresolved, character(0))
})

test_that("the generator page cites only allow-listed DOIs", {
  skip_on_cran()
  skip_if_no_description()
  man <- file.path(cier_pkg_root(), "man", generator_doc_file())
  skip_if(!file.exists(man), "cier_simulate.Rd not on disk")
  unlisted <- unlisted_generator_dois(man)
  if (length(unlisted) > 0L) {
    cat("\nDOIs on the generator page outside the allow-list:\n")
    cat(unlisted, sep = "\n")
    cat("\n")
  }
  expect_identical(unlisted, character(0))
})

test_that("an off-allow-list DOI on the generator page is rejected", {
  # A right-looking citation carrying a DOI outside the allow-list must be caught
  # (the DOI half of the guard: typo'd or uncited new DOIs fail).
  lines <- c("\\references{",
             "Author, A. (2024). A real-looking but uncited paper.",
             "\\doi{10.9999/not-in-the-allow-list}", "}")
  rd <- local({
    f <- tempfile(fileext = ".Rd")
    writeLines(lines, f)
    f
  })
  expect_identical(unlisted_generator_dois(rd), "10.9999/not-in-the-allow-list")
})

test_that("man pages carry the frozen canonical citation cores", {
  skip_on_cran()
  skip_if_no_description()
  man_dir <- file.path(cier_pkg_root(), "man")
  skip_if(!dir.exists(man_dir), "man/ not on disk")
  cores <- canonical_reference_cores()
  pages <- citation_page_keys()
  for (page in names(pages)) {
    rd <- file.path(man_dir, page)
    # A pinned page that has been renamed or removed must FAIL, not silently skip
    # (else the guard goes vacuously green).
    expect_true(file.exists(rd), info = paste(page, "missing from man/"))
    if (!file.exists(rd)) next
    block <- normalize_ref_block(rd)
    for (key in pages[[page]]) {
      expect_true(grepl(cores[[key]], block, fixed = TRUE),
                  info = paste(page, key))
    }
  }
})

test_that("only the generator page cites DOIs beyond the registry", {
  # The generator page is the ONLY non-data man page allowed to cite beyond the
  # registry (its allow-list guards it). If a second generator-like page is added
  # and excluded from the registry guard without an allow-list, this fails --
  # closing the "new page silently escapes both DOI guards" hole.
  skip_on_cran()
  skip_if_no_description()
  man_dir <- file.path(cier_pkg_root(), "man")
  skip_if(!dir.exists(man_dir), "man/ not on disk")
  rd <- list.files(man_dir, pattern = "\\.Rd$", full.names = TRUE)
  rd <- Filter(function(f) !is_data_doc(f), rd)
  skip_if(length(rd) == 0L, "no man pages on disk")
  reg_dois <- unique(stats::na.omit(cier_methods()$doi))
  beyond <- Filter(
    function(f) length(setdiff(dois_in_references(f), reg_dois)) > 0L, rd
  )
  expect_identical(basename(beyond), generator_doc_file())
})

test_that("the frozen citation core discriminates fabrications", {
  # Negative control: pins the matcher's discriminating power independent of the
  # live man pages, so a future re-fabrication is caught even once the pages are
  # green.
  good <- canonical_reference_cores()[["goldammer2024_revisited"]]
  # The exact core, surrounded by other references, matches (containment intended).
  expect_true(grepl(good, paste("A. (2020). X.", good, "B. (2019). Y."),
                    fixed = TRUE))
  # The in-the-wild fabrication (wrong title, truncated authors) does not.
  fabricated <- paste(
    "Goldammer, P., Annen, H., Stockli, P. L., & Jonas, K. (2024).",
    "On the utility of indirect methods for detecting faking and",
    "careless responding."
  )
  expect_false(grepl(good, fabricated, fixed = TRUE))
  # A single dropped author does not.
  dropped <- sub("Escher, Y. A., ", "", good, fixed = TRUE)
  expect_false(grepl(good, dropped, fixed = TRUE))
})

test_that("the DOI extractor reads a DOI from a references block", {
  lines <- c("\\references{", "Author (2024). Title.",
             "\\doi{10.3758/s13428-024-02506-0}", "}")
  expect_identical(
    dois_in_references(local({
      f <- tempfile(fileext = ".Rd")
      writeLines(lines, f)
      f
    })),
    "10.3758/s13428-024-02506-0"
  )
})

test_that("the DOI extractor handles brace-adjacent, prose, and multi DOIs", {
  mk <- function(lines) {
    f <- tempfile(fileext = ".Rd")
    writeLines(lines, f)
    f
  }
  # A DOI flush against the closing brace keeps the brace out of the match.
  expect_identical(
    dois_in_references(mk(c("\\references{", "X. \\doi{10.1000/abc}", "}"))),
    "10.1000/abc"
  )
  # A trailing sentence period is stripped.
  expect_identical(
    dois_in_references(mk(c("\\references{", "See 10.1000/xyz.", "}"))),
    "10.1000/xyz"
  )
  # Two DOIs in one block are both extracted.
  expect_setequal(
    dois_in_references(mk(c("\\references{", "\\doi{10.1000/a}",
                            "\\doi{10.2000/b}", "}"))),
    c("10.1000/a", "10.2000/b")
  )
})
