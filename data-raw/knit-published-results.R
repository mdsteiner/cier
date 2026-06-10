# Regenerate the published-results vignette and its committed results table.
#
# Knits vignettes/published-results.Rmd.orig (which fetches the external
# datasets into tools::R_user_dir("cier", "cache") and runs every comparison)
# into the shipped vignettes/published-results.Rmd with frozen output, and
# rewrites inst/extdata/published-results.csv. Needs network access, the
# suggested packages PerFit / mokken / mice, and an installed cier matching
# the working tree. Run from the package root. The knit happens in a fresh
# callr session so a backend crash cannot take the main session down.

stopifnot(file.exists("DESCRIPTION"))
callr::r(
  function() {
    # knit() evaluates chunks in the input file's directory by default;
    # pin the package root so inst/extdata/... resolves correctly
    knitr::opts_knit$set(root.dir = getwd())
    knitr::knit("vignettes/published-results.Rmd.orig",
                output = "vignettes/published-results.Rmd")
  },
  show = TRUE
)
message("wrote vignettes/published-results.Rmd and ",
        "inst/extdata/published-results.csv")
