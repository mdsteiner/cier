# Independent reference implementation of the IRV index per Marjanovic,
# Holden, Struthers, Cribbie, & Greenglass (2015, *Personality and Individual
# Differences*, 84, 79-83).
#
# Definition (Marjanovic et al., 2015, p. 80):
#   IRV_r = sqrt( sum_j (x_rj - mean(x_r))^2 / (J - 1) )
# i.e. the *sample* standard deviation (denominator n - 1) of respondent r's
# answers across the items they answered (available cases, na.rm).
#
# This oracle re-derives the statistic with a hand-rolled two-pass
# sum-of-squares and NEVER calls the production kernel (which uses
# matrixStats::rowSds), so any divergence is attributable to the kernel, not a
# shared helper. A row with fewer than two answered items has no defined sample
# SD and abstains (NA).

ref_irv_row <- function(row) {
  present <- row[!is.na(row)]
  n <- length(present)
  if (n < 2L) {
    return(NA_real_)
  }
  mean_r <- sum(present) / n
  ss <- sum((present - mean_r)^2)
  sqrt(ss / (n - 1L))
}

ref_irv <- function(x) {
  if (!is.matrix(x)) x <- as.matrix(x)
  unname(apply(x, 1L, ref_irv_row))
}
