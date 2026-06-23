# Independent reference for the response-time generator, the page aggregation,
# and the attention-check injection that cier_simulate() couples to the planted
# careless plan.
#
# These oracles re-derive each quantity from scratch and NEVER call the
# production kernels in R/sim-times.R / R/sim-direct.R, so any divergence is
# attributable to the kernel, not a shared helper. Like the RPR oracle, the
# randomised re-derivations coordinate their DRAW ORDER with production.
# The simulator is a generator with no CRAN parity partner (verified
# 2026-06-11), so trust is oracle-only.

# Per-cell lognormal response times, re-derived by hand. Production's frozen
# draw order is: (1) the n*p cell noise, standard-normal, filled column-major;
# (2) the n per-respondent pace intercepts. The careless log-mean mu_car
# applies exactly on the INCLUSIVE careless span [onset_i, offset_i] of a
# careless row; everywhere else the attentive mu_att applies. Each cell is
# exp((sigma * z + mu) + respondent_sd * b_i), floored at min_seconds.
ref_times_lognormal <- function(seed, n, p, careless, onset, offset,
                                mu_att = log(8), mu_car = log(1.5),
                                sigma = 0.5, respondent_sd = 1.2,
                                min_seconds = 0.1) {
  set.seed(seed)
  z <- stats::rnorm(n * p)
  b <- stats::rnorm(n)
  out <- matrix(NA_real_, nrow = n, ncol = p)
  for (i in seq_len(n)) {
    for (j in seq_len(p)) {
      in_span <- careless[[i]] && j >= onset[[i]] && j <= offset[[i]]
      mu <- if (in_span) mu_car else mu_att
      log_rt <- (sigma * z[[(j - 1L) * n + i]] + mu) + respondent_sd * b[[i]]
      out[i, j] <- max(exp(log_rt), min_seconds)
    }
  }
  out
}

# Page totals re-derived by an explicit per-respondent, per-page double loop
# over contiguous item blocks (pages are consecutive: page g covers the
# items_per_page[g] items after page g - 1's last item).
ref_page_totals <- function(cell_times, items_per_page) {
  ends <- cumsum(items_per_page)
  starts <- ends - items_per_page + 1L
  out <- matrix(NA_real_, nrow = nrow(cell_times), ncol = length(items_per_page))
  for (g in seq_along(items_per_page)) {
    for (i in seq_len(nrow(cell_times))) {
      total <- 0
      for (j in starts[[g]]:ends[[g]]) {
        total <- total + cell_times[i, j]
      }
      out[i, g] <- total
    }
  }
  out
}

# Attention-check injection re-derived by hand. Production's frozen draw order,
# per check column j: (1) the passing key, one uniform draw over the five
# response options; (2) one runif(n) failure draw against the per-row failure
# probability (p_fail_careless on careless rows, p_fail_attentive otherwise);
# (3) one uniform draw over the FOUR non-key options per failing row, in row
# order. Passing rows answer the key exactly; failing rows answer a non-key
# value, so failure <=> response != key (the cier_attention counting rule).
ref_direct_checks <- function(seed, n, careless, n_checks,
                              p_fail_careless = 0.75,
                              p_fail_attentive = 0.05) {
  set.seed(seed)
  checks <- matrix(NA_integer_, nrow = n, ncol = n_checks)
  pass <- vector("list", n_checks)
  p_row <- ifelse(careless, p_fail_careless, p_fail_attentive)
  for (j in seq_len(n_checks)) {
    key <- sample.int(5L, 1L)
    fail <- stats::runif(n) < p_row
    wrong <- sample.int(4L, sum(fail), replace = TRUE)
    column <- rep(key, n)
    column[fail] <- (seq_len(5L)[-key])[wrong]
    checks[, j] <- as.integer(column)
    pass[[j]] <- key
  }
  list(checks = checks, pass = pass)
}
