# Purpose: The cross-index agreement diagnostic. A per-index empirical
#          percentile flags exactly its target rate by construction, so that
#          per-index rate is tautological and uninformative. What is
#          informative is whether flags cluster on the SAME respondents across
#          votes -- the signature of a real careless subgroup. `flag_agreement()`
#          contrasts the observed share flagged by >= k votes with the share
#          expected if the votes fired independently (the Poisson-binomial of
#          the per-vote rates), and carries a per-vote table whose excess over a
#          calibrated null is informative for the null-referenced indices.
# Args:    See flag_agreement() below.
# Returns: list(agreement, per_vote) data.frames.
# Invariants:
#   - `expected` is the EXACT Poisson-binomial upper tail of colMeans(flags);
#     it handles unequal per-vote rates (a plain binomial would not).
#   - The independence baseline is a null of "no shared signal", not a claim the
#     votes are independent; observed >> expected makes contamination visible.
#   - A vote is `informative` only where a real null nominal is supplied
#     (Mahalanobis chi-square; the person-fit Monte-Carlo nulls); the empirical
#     percentile votes carry NA and are not informative.

# Exact Poisson-binomial pmf for m independent Bernoulli(p_i), by iterative
# convolution of the m two-point pmfs [1 - p_i, p_i]. Returns a numeric vector
# of length m + 1 whose element j is P(count = j - 1); order-independent, sums
# to 1.
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

# Purpose: Contrast observed cross-vote agreement with the independence baseline.
# Args:
#   flags     - logical matrix, respondents (rows) x votes (columns), no NA.
#   null_rate - optional numeric vector, one nominal null rate per vote (NA for
#               the tautological empirical-percentile votes); length ncol(flags).
#   call      - calling environment for typed conditions.
# Returns: list(agreement, per_vote):
#   $agreement = data.frame(k, observed, expected) for k = 1..ncol(flags), where
#     observed[k] is the share of respondents flagged by >= k votes and
#     expected[k] is the Poisson-binomial P(count >= k) under independence.
#   $per_vote  = data.frame(vote, observed, null, excess, informative).
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
