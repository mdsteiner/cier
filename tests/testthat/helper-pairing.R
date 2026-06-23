# Shared input-validation assertions for the synonym / antonym pairing indices
# (cier_psychsyn, cier_psychant). Both share check_responses() + the
# check_percentile_overrides() contract (a literal cutoff is a correlation in [-1, 1])
# and resolve_pairing_cor()'s `reference` validation, so the rejection contract is
# identical -- only the wrapper and its valid fixture differ, passed in by the caller.
# Sign / direction-specific behaviour (critical_r magnitude, cutoff tail, oracle parity)
# stays in each file.

# A bad payload (non-matrix / non-numeric / non-finite) or a malformed fpr / cutoff
# override (out of (0, 1) / outside [-1, 1] / both supplied) is a typed input error.
expect_pairing_overrides_rejected <- function(fn, x) {
  expect_error(fn(1:10), class = "cier_error_input")
  expect_error(fn(matrix(letters[1:6], nrow = 2L)), class = "cier_error_input")
  bad <- x
  bad[1L, 1L] <- Inf
  expect_error(fn(bad), class = "cier_error_input")
  expect_error(fn(x, fpr = 0), class = "cier_error_input")
  expect_error(fn(x, fpr = 1), class = "cier_error_input")
  expect_error(fn(x, fpr = c(0.05, 0.1)), class = "cier_error_input")
  expect_error(fn(x, cutoff = -1.5), class = "cier_error_input")
  expect_error(fn(x, cutoff = 1.5), class = "cier_error_input")
  expect_error(fn(x, cutoff = NA_real_), class = "cier_error_input")
  expect_error(fn(x, fpr = 0.1, cutoff = 0.5), class = "cier_error_input")
  invisible(NULL)
}

# Every malformed `reference` selector is a typed input error: a wrong-length mask, a
# mask selecting nobody, an out-of-range index, fewer than three selected / external
# rows, a wrong column count, and a non-numeric external sample.
expect_pairing_reference_rejected <- function(fn, x) {
  expect_error(fn(x, reference = rep(TRUE, 5L)), class = "cier_error_input")
  expect_error(fn(x, reference = rep(FALSE, nrow(x))), class = "cier_error_input")
  expect_error(fn(x, reference = c(1L, 2L, 999L)), class = "cier_error_input")
  expect_error(fn(x, reference = c(1L, 2L)), class = "cier_error_input")
  expect_error(fn(x, reference = x[1:2, ]), class = "cier_error_input")
  expect_error(fn(x, reference = x[, 1:3]), class = "cier_error_input")
  expect_error(fn(x, reference = matrix("a", 5L, ncol(x))), class = "cier_error_input")
  invisible(NULL)
}
