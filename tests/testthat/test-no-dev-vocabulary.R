# Guard: developer-process vocabulary must never reach the packaged tree.
# Help pages, runtime strings, vignettes, and shipped data are read by end
# users; build-time scaffolding ("Phase 1", "Slice 0", a "ticket", an
# internal planning-doc path) has no meaning to them. This scans the
# packaged surface (R/, man/, vignettes/, inst/); tests/, dev/, and
# archive/ are out of scope (the latter two are .Rbuildignored).
#
# The token set is deliberately narrow and the domain term "decision
# rule" is allowed (a "Decision" hit requires a following number).

dev_vocab_patterns <- function() {
  c(
    "Phase\\s*[0-9]",
    "Step\\s*[0-9]",
    "Slice\\s*[0-9]",
    "\\bCard\\b",
    "Decision\\s*#?\\s*[0-9]",
    "[Tt]icket",
    "[Ii]mplementation plan",
    "DECISIONS",
    "dev/restart",
    "dev-roadmap",
    "\\b(plan|architecture|index-specs|example-data)\\.md\\b"
  )
}

shipped_text_files <- function(root) {
  out <- character(0)
  r_dir <- file.path(root, "R")
  if (dir.exists(r_dir)) {
    out <- c(out, list.files(r_dir, pattern = "\\.R$", full.names = TRUE))
  }
  man_dir <- file.path(root, "man")
  if (dir.exists(man_dir)) {
    out <- c(out, list.files(man_dir, pattern = "\\.Rd$", full.names = TRUE))
  }
  vig_dir <- file.path(root, "vignettes")
  if (dir.exists(vig_dir)) {
    out <- c(out, list.files(vig_dir, pattern = "\\.(Rmd|qmd|Rnw)$",
                             full.names = TRUE))
  }
  inst_dir <- file.path(root, "inst")
  if (dir.exists(inst_dir)) {
    out <- c(out, list.files(inst_dir, recursive = TRUE, full.names = TRUE,
                             pattern = "\\.(csv|md|txt|R|Rmd)$"))
  }
  out
}

scan_for_dev_vocab <- function(files) {
  pats <- dev_vocab_patterns()
  hits <- character(0)
  for (f in files) {
    lines <- readLines(f, warn = FALSE, encoding = "UTF-8")
    for (p in pats) {
      matched <- grep(p, lines, perl = TRUE)
      for (ln in matched) {
        hits <- c(hits, sprintf("%s:%d: /%s/ -> %s",
                                basename(f), ln, p, trimws(lines[[ln]])))
      }
    }
  }
  hits
}

test_that("the packaged tree carries no developer-process vocabulary", {
  skip_on_cran()
  skip_if_no_description()
  files <- shipped_text_files(cier_pkg_root())
  skip_if(length(files) == 0L, "no packaged text files on disk")
  hits <- scan_for_dev_vocab(files)
  if (length(hits) > 0L) {
    cat("\nForbidden developer vocabulary in the packaged tree:\n")
    cat(hits, sep = "\n")
    cat("\n")
  }
  expect_identical(hits, character(0))
})

test_that("the guard actually fires on a planted token", {
  tmp <- tempfile(fileext = ".R")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(c("# implement Phase 2 of the ticket", "x <- 1"), tmp)
  expect_gt(length(scan_for_dev_vocab(tmp)), 0L)
})

test_that("the guard does not fire on the allowed domain term", {
  tmp <- tempfile(fileext = ".R")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(c("# the decision rule combines indices by agreement"), tmp)
  expect_identical(scan_for_dev_vocab(tmp), character(0))
})
