# Attention-check kernel for cier_simulate(). Builds a checks matrix + pass-set list shaped
# for cier_attention(checks, pass): each check is a five-point instructed item whose single
# passing value travels in `pass`. Careless rows fail each check with p_fail_careless,
# attentive rows with p_fail_attentive, so the failed-check count across several checks
# separates the groups. Frozen draw order per check (RNG reproducibility): key, then n
# failure draws, then one non-key value per failing row; n_checks = 0 draws nothing.

sim_direct_checks <- function(n, careless, n_checks,
                              p_fail_careless = 0.75,
                              p_fail_attentive = 0.05,
                              call = rlang::caller_env()) {
  check_count(n, "n", call = call)
  check_count_nonneg(n_checks, "n_checks", call = call)
  if (!is.logical(careless) || length(careless) != n || anyNA(careless)) {
    cier_abort(
      "cier_error_input",
      "{.arg careless} must be a logical vector of length {n} without NA.",
      data = list(arg = "careless"), call = call
    )
  }
  check_number(p_fail_careless, "p_fail_careless",
               lower = 0, upper = 1, call = call)
  check_number(p_fail_attentive, "p_fail_attentive",
               lower = 0, upper = 1, call = call)
  if (n_checks == 0L) {
    return(list(checks = NULL, pass = NULL))
  }
  n_checks <- as.integer(n_checks)
  p_row <- ifelse(careless, p_fail_careless, p_fail_attentive)
  checks <- matrix(NA_integer_, nrow = n, ncol = n_checks)
  pass <- vector("list", n_checks)
  for (j in seq_len(n_checks)) {
    key <- sample.int(5L, 1L)
    fail <- stats::runif(n) < p_row
    wrong <- sample.int(4L, sum(fail), replace = TRUE)
    column <- rep(key, n)
    column[fail] <- (seq_len(5L)[-key])[wrong]
    checks[, j] <- column
    pass[[j]] <- key
  }
  list(checks = checks, pass = pass)
}
