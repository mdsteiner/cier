# Cross-index agreement diagnostic. A per-index empirical percentile flags its target rate
# by construction, so the per-index rate is tautological; what's informative is whether
# flags cluster on the same respondents across votes (the signature of a real careless
# subgroup). flag_agreement() contrasts the observed share flagged by >= k votes with the
# share expected under independence (exact Poisson-binomial tail of colMeans(flags)). The
# independence baseline is a null of "no shared signal", not a claim the votes are
# independent. A vote is `informative` only where a real null nominal is supplied
# (Mahalanobis chi-square; person-fit Monte-Carlo nulls). `per_vote$excess` = observed -
# null is descriptive: `observed` is over all respondents but the null applies only to
# scored ones, so an abstaining null-referenced index biases the excess low.

# Exact Poisson-binomial pmf for m independent Bernoulli(p_i), by iterative convolution of
# the m two-point pmfs. Returns length m + 1; element j is P(count = j - 1), sums to 1.
poisson_binomial_pmf <- function(p) {
  pmf <- 1
  for (pi in p) {
    pmf <- c(pmf * (1 - pi), 0) + c(0, pmf * pi)
  }
  pmf
}

# Validate + coerce the flag matrix: a non-empty logical matrix with no missing
# values (a data.frame of logical columns is accepted and coerced).
check_flag_matrix <- function(flags, call) {
  flags <- as.matrix(flags)
  if (!is.logical(flags) || nrow(flags) == 0L || ncol(flags) == 0L ||
        anyNA(flags)) {
    cier_abort(
      "cier_error_input",
      c("{.arg flags} must be a non-empty logical matrix with no missing values.",
        "i" = "One TRUE/FALSE column per vote, one row per respondent."),
      data = list(arg = "flags"), call = call
    )
  }
  flags
}

# A calibrated-null nominal per vote: NULL -> all NA (every vote tautological);
# otherwise a numeric vector the width of the flag matrix.
resolve_null_rate <- function(null_rate, m, call) {
  if (is.null(null_rate)) {
    return(rep(NA_real_, m))
  }
  if (!is.numeric(null_rate) || length(null_rate) != m) {
    cier_abort(
      "cier_error_input",
      "{.arg null_rate} must be a numeric vector of length {m} (one per vote).",
      data = list(arg = "null_rate", observed = length(null_rate)), call = call
    )
  }
  as.numeric(null_rate)
}

# `flags`: logical respondents x votes matrix (no NA). `null_rate`: optional per-vote
# nominal null (NA for tautological percentile votes). Returns list(agreement, per_vote):
# agreement$observed[k]/expected[k] are observed and Poisson-binomial share flagged by
# >= k votes; per_vote carries the per-vote rate, null, excess, and informative flag.
flag_agreement <- function(flags, null_rate = NULL,
                           call = rlang::caller_env()) {
  flags <- check_flag_matrix(flags, call)
  m <- ncol(flags)
  ks <- seq_len(m)
  p <- as.numeric(colMeans(flags))
  counts <- rowSums(flags)
  observed <- vapply(ks, function(k) mean(counts >= k), numeric(1L))
  pmf <- poisson_binomial_pmf(p)
  expected <- vapply(ks, function(k) sum(pmf[(k + 1L):(m + 1L)]), numeric(1L))

  null <- resolve_null_rate(null_rate, m, call)
  votes <- colnames(flags)
  if (is.null(votes)) {
    votes <- paste0("vote", ks)
  }
  list(
    agreement = data.frame(k = ks, observed = observed, expected = expected),
    per_vote = data.frame(
      vote        = votes,
      observed    = p,
      null        = null,
      excess      = p - null,
      informative = !is.na(null),
      stringsAsFactors = FALSE
    )
  )
}
