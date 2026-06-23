# Paper-faithful reference implementation of the longstring index per
# Johnson (2005, *J. Res. Pers.*) and Meade & Craig (2012, *Psychol.
# Methods*). Used only by tests/testthat/test-cier-longstring.R as a
# slow but unambiguous cross-check of the production kernel.
#
# Definition (Johnson, 2005, p. 116; Meade & Craig, 2012, p. 442):
#   For each respondent, count the maximum length of consecutive
#   identical responses across the items as administered. Higher values
#   indicate straightlining.
#
# This reference operates row by row using a hand-rolled inner loop so
# any divergence from `cier_longstring()` would be attributable to the
# vectorised kernel, not to the algorithm definition. It never calls the
# production kernel.

ref_longstring_row <- function(row) {
  if (length(row) == 0L) {
    return(c(longest = NA_real_, average = NA_real_))
  }
  lengths <- integer(0L)
  prev <- row[[1L]]
  run <- 1L
  if (length(row) > 1L) {
    for (j in seq.int(2L, length(row))) {
      cur <- row[[j]]
      # Match base::rle(): any NA in either position breaks the run.
      same <- !is.na(prev) && !is.na(cur) && prev == cur
      if (same) {
        run <- run + 1L
      } else {
        lengths <- c(lengths, run)
        prev <- cur
        run <- 1L
      }
    }
  }
  lengths <- c(lengths, run)
  c(longest = max(lengths), average = mean(lengths))
}

ref_longstring <- function(x) {
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  out <- t(apply(x, 1, ref_longstring_row))
  rownames(out) <- NULL
  as.data.frame(out)
}
