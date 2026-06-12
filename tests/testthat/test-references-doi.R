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
# The simulator's page (cier_simulate.Rd) is likewise exempt: cier_simulate()
# is a data generator, not an index -- by the ADR it has NO registry row (no
# cutoff, no direction, never in cier_screen()), and its references cite the
# published simulation designs and calibration anchors it follows, not method
# citations the registry could carry.
generator_doc_exemptions <- function() {
  "cier_simulate.Rd"
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

test_that("every @references DOI resolves to a registry citation entry", {
  skip_on_cran()
  skip_if_no_description()
  man_dir <- file.path(cier_pkg_root(), "man")
  skip_if(!dir.exists(man_dir), "man/ not on disk")
  rd <- list.files(man_dir, pattern = "\\.Rd$", full.names = TRUE)
  rd <- rd[!basename(rd) %in% generator_doc_exemptions()]
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
