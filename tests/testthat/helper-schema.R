# Shared schema assertion for the per-index cier_index contract: every public index
# wrapper returns the pinned 7-field list with $value (double) / $flagged (logical) of
# length n and the registry-declared $method / $direction. Factored so each index test
# pins its own wrapper's schema in one line (it still calls the wrapper with that
# index's own fixture, so a wrong method/direction/shape is caught per index). The
# generic new_cier_index invariants and the as.data.frame.cier_index row-frame contract
# are tested directly, once, in test-cier-index.R.
expect_cier_index_schema <- function(out, method, direction, n) {
  expect_s3_class(out, "cier_index")
  expect_type(out, "list")
  expect_identical(names(out),
                   c("value", "flagged", "method", "cutoff", "direction",
                     "cutoff_method", "cutoff_rate"))
  expect_type(out$value, "double")
  expect_type(out$flagged, "logical")
  expect_identical(length(out$value), n)
  expect_identical(length(out$flagged), n)
  expect_identical(out$method, method)
  expect_identical(out$direction, direction)
  expect_type(out$cutoff, "double")
  invisible(out)
}
