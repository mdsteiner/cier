# Guard: the generated documentation (man/*.Rd) and NAMESPACE must match
# what roxygen2 would produce from the current sources. Runs roxygen on a
# throwaway copy and compares bytewise, so a stale man page or a
# forgotten @export is caught. Skipped when the installed roxygen2 differs
# from the recorded RoxygenNote (cosmetic format drift is not a defect).

test_that("man/ and NAMESPACE are up to date with roxygen2", {
  skip_on_cran()
  skip_if_not_installed("roxygen2")
  skip_if_not_installed("callr")
  skip_if_no_description()
  root <- cier_pkg_root()
  # roxygen2 >= 7.4 records its version in Config/roxygen2/version; older
  # versions used RoxygenNote. Accept whichever is present.
  desc <- read.dcf(file.path(root, "DESCRIPTION"))
  field <- if ("Config/roxygen2/version" %in% colnames(desc)) {
    "Config/roxygen2/version"
  } else {
    "RoxygenNote"
  }
  recorded <- if (field %in% colnames(desc)) desc[1, field] else NA_character_
  installed <- as.character(utils::packageVersion("roxygen2"))
  skip_if(is.na(recorded) || recorded != installed,
          sprintf("recorded roxygen2 %s != installed %s",
                  recorded, installed))

  tmp <- tempfile("cier-roxy-")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)
  # `data` is needed so roxygen can resolve the bundled-dataset documentation.
  for (item in c("DESCRIPTION", "NAMESPACE", "R", "man", "inst", "data")) {
    src <- file.path(root, item)
    if (file.exists(src)) {
      file.copy(src, tmp, recursive = TRUE)
    }
  }
  # Run roxygen in a fresh subprocess: roxygenise() calls load_all(), which
  # would otherwise reload the package namespace mid-suite and corrupt the
  # surrounding tests.
  callr::r(
    function(p) suppressMessages(roxygen2::roxygenise(p)),
    args = list(p = tmp)
  )

  read_norm <- function(path) {
    if (!file.exists(path)) return(character(0))
    readLines(path, warn = FALSE, encoding = "UTF-8")
  }
  expect_identical(read_norm(file.path(tmp, "NAMESPACE")),
                   read_norm(file.path(root, "NAMESPACE")))

  rel <- function(dir) sort(list.files(dir, recursive = TRUE))
  src_man <- file.path(root, "man")
  tmp_man <- file.path(tmp, "man")
  expect_identical(rel(tmp_man), rel(src_man))
  for (f in rel(src_man)) {
    expect_identical(read_norm(file.path(tmp_man, f)),
                     read_norm(file.path(src_man, f)),
                     info = f)
  }
})
