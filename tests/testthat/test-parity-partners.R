# Fail-loud guard for the cross-package parity layer. The parity checks against
# careless / psych / PerFit / mokken are one of the two trust pillars, but every
# parity test opens with
# skip_if_not_installed(): if a partner silently failed to install on a runner,
# the whole layer would skip and the suite would still report green. CI sets
# CIER_REQUIRE_PARITY=true (the tests-and-guards job), turning a missing partner
# into a FAILURE there; locally the partners stay optional and this test skips.

test_that("every cross-package parity partner is installed when required", {
  if (!identical(tolower(Sys.getenv("CIER_REQUIRE_PARITY")), "true")) {
    skip("CIER_REQUIRE_PARITY not set (parity partners optional locally)")
  }
  for (pkg in c("careless", "psych", "PerFit", "mokken")) {
    expect_true(requireNamespace(pkg, quietly = TRUE),
                label = sprintf("parity partner '%s' is installed", pkg))
  }
})
