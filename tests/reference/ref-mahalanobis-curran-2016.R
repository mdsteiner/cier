# Independent reference implementation of the Mahalanobis-distance C/IER index
# per Curran (2016, *Journal of Experimental Social Psychology*, 66, 4-19). This
# reproduces the algorithm `psych::outlier()` uses, which is what
# `careless::mahad()` calls.
#
# Definition:
#   For each respondent r:
#     D^2_r = (x_r - mu)^T Sigma^{-1} (x_r - mu)
#   where Sigma = cov(x, use = "pairwise") and mu = colMeans(x), both over
#   available cases. Missing cells are handled by zero-filling the *centred*
#   coordinate (xc[m] <- 0) before the bilinear form, which drops every term of
#   the quadratic form that involves coordinate m -- exactly what
#   `psych::outlier(na.rm = TRUE)` and `careless::mahad()` compute. This is NOT
#   `sum(xc * (Sigma^{-1} xc), na.rm = TRUE)`: one NA in `xc` makes the whole
#   `Sigma^{-1} xc` product NA, so `na.rm` would collapse the row to 0 instead of
#   dropping only the missing terms.
#
# This oracle re-derives the statistic with an explicit per-respondent
# quadratic-form loop (`vapply`) and NEVER calls the production kernel (which
# evaluates D^2 = diag(X Sigma^{-1} X^T) with a single BLAS matrix multiply), so
# any divergence is attributable to the kernel, not a shared reduction. Rows with
# no data abstain (NA); a singular covariance or fewer than two respondents with
# data yields all-NA.

ref_mahalanobis <- function(x) {
  if (!is.matrix(x)) x <- as.matrix(x)
  all_na <- apply(x, 1L, function(r) all(is.na(r)))
  d <- rep(NA_real_, nrow(x))
  if (sum(!all_na) < 2L) {
    return(d)
  }
  xf <- x[!all_na, , drop = FALSE]
  sx <- stats::cov(xf, use = "pairwise")
  sx_inv <- tryCatch(solve(sx), error = function(e) NULL)
  if (is.null(sx_inv)) {
    return(d)
  }
  centered <- scale(xf, scale = FALSE)
  d2 <- vapply(seq_len(nrow(centered)), function(i) {
    xc <- centered[i, ]
    # Zero-fill the missing centred coordinates so they drop out of the bilinear
    # form. `na.rm = TRUE` on the product would be wrong: one NA in `xc` turns
    # the whole `sx_inv %*% xc` product into NA, collapsing the row to 0 rather
    # than skipping only the missing terms.
    xc[is.na(xc)] <- 0
    sum(xc * (sx_inv %*% xc))
  }, numeric(1L))
  d[!all_na] <- d2
  d
}
