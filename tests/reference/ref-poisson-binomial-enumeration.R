# Independent reference for the Poisson-binomial upper tail used by the
# agreement diagnostic `flag_agreement()`.
#
# Setup: each of m votes flags a respondent independently with its own
# probability p[i]. The number of votes flagging a respondent is then a sum
# of m independent Bernoulli(p[i]) variables -- a Poisson-binomial. The
# diagnostic needs the upper tail P(count >= k) for k = 1..m.
#
# This reference brute-enumerates all 2^m flag patterns, weights each pattern
# by its independence probability prod(p^flag * (1 - p)^(1 - flag)), and sums
# the weights over patterns with at least k flags. The production kernel
# `poisson_binomial_pmf()` instead convolves the m two-point pmfs
# [1 - p[i], p[i]]; enumeration vs convolution are genuinely independent
# derivations, cross-checked to 1e-12 in test-diagnostics.R (see
# tests/reference/TOLERANCES.md). Valid only for small m (2^m patterns); the
# tests use m <= 10.

ref_poisson_binomial_tail <- function(p) {
  m <- length(p)
  patterns <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), m)))
  prob <- apply(patterns, 1L, function(flag) prod(ifelse(flag, p, 1 - p)))
  count <- rowSums(patterns)
  vapply(seq_len(m), function(k) sum(prob[count >= k]), numeric(1L))
}
