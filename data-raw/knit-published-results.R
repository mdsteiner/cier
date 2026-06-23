# Regenerate the published-results vignette and its committed results table.
#
# ONE COMMAND (from the package root, on a networked machine):
#   Rscript data-raw/knit-published-results.R
#
# Knits vignettes/published-results.Rmd.orig (which fetches the external
# datasets into tools::R_user_dir("cier", "cache") and runs every comparison)
# into the shipped vignettes/published-results.Rmd with frozen output, and
# rewrites inst/extdata/published-results.csv. Needs network access, the
# suggested package mice (imputation), and an installed cier matching the working
# tree. No person-fit backend is needed -- Gnormed and Ht score in pure R now;
# PerFit / mokken are only the tests' parity oracles. Run from the package root.
# The knit happens in a fresh callr session so a crash cannot take the main
# session down. NOTE: callr is NOT a declared package dependency (it is dev-only
# and was dropped from Suggests), so install it manually if this script errors
# with "no package called 'callr'": install.packages("callr").
#
# WHEN TO RE-RUN: after any change that moves a cier_value cell -- notably a
# change to a randomised index's seeded output (the RPR mechanism switched to a
# vectorised uniform-split generator, so all RPR / RPR25 / RPR50 / RPR100 cells
# and the Goldammer battery RPR AUCs shift; the statistic is unchanged), OR a
# change to a non-random index's exact output (Gnormed moved from PerFit, which
# rounds to 4 dp, to cier's exact in-package closed form -- the extra precision can
# reorder near-tied scores and shift the Gnormed AUC cells). The re-knit recomputes
# every cier_value, delta, and status from the fetched data.
#
# NO TEST EDIT NEEDED afterwards: tests/testthat/test-published-results-table.R
# pins only the PAPER-SIDE transcription (pubres_paper_digest covers paper_value
# / index / n, not the cier columns -- see its comment), so the frozen digest is
# invariant under a re-knit. The Bruhlmann block (bundled data, no network) is
# already committed at the post-switch values (even-odd 88, RPR 87); the re-knit
# reproduces them. If a re-knit trips a coarse floor or an allowed-differs cap in
# that test, the shifted cell is a genuine review point, not a transcription drift.

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
