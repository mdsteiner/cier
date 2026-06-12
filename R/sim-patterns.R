# Purpose: Careless-pattern block producers for cier_simulate(). Each producer
#          fills an n_rows x n_cols integer block with one C/IER pattern in the
#          per-item raw coding (min_j..max_j); the content engine (R/sim-plan.R)
#          writes the block into a respondent's careless span. The clean
#          (attentive) layer it overwrites is R/sim-attentive.R. Pure internal
#          kernels: no public cier_simulate() and no S3 object yet.
# Args:    n_rows (block height); per-column mins / maxs already subset to the
#          active span; a pattern-specific `params` list; `call` threaded so a bad
#          parameter is attributed to the public caller.
# Returns: An n_rows x length(mins) integer matrix with values in min_j..max_j.
# Invariants:
#   - One production implementation per pattern; reference re-derivations stay in
#     tests/reference.
#   - Sequential patterns (diagonal, alternating, markov) walk one K_min =
#     min(categories) over the span and offset per item by min_j, so values stay in
#     range on a heterogeneous battery.
#   - speeder is a recognised pattern NAME (attentive content, careless times in a
#     later timing layer) but is NOT a content producer.

# ---- Shared position helpers ------------------------------------------------

# Round half UP: floor(x + 0.5). The agreed position convention (a 4-option
# midpoint is the upper middle, 3), NOT base-R round()'s round-half-to-even.
sim_round_half_up <- function(x) {
  floor(x + 0.5)
}

# The n_rows x p value grid at relative position q (length n_rows): per cell
# round_half_up(min_j + q_i * (max_j - min_j)). q recycles down rows, min / max
# across columns.
sim_position_grid <- function(q, mins, maxs) {
  n <- length(q)
  p <- length(mins)
  qm <- matrix(q, n, p)
  mn <- matrix(mins, n, p, byrow = TRUE)
  mx <- matrix(maxs, n, p, byrow = TRUE)
  sim_round_half_up(mn + qm * (mx - mn))
}

# Resolve a per-row scalar parameter: NULL draws one per row via `draw`; otherwise
# a numeric scalar or a length-n_rows vector in [lo, hi] (whole when `whole`),
# recycled to n_rows. Used for the straightline position q, the diagonal start, and
# the alternating start_offset.
sim_resolve_row_param <- function(value, n_rows, draw, lo, hi, whole, arg, call) {
  if (is.null(value)) {
    return(draw(n_rows))
  }
  ok <- checkmate::test_numeric(value, lower = lo, upper = hi,
                                any.missing = FALSE, finite = TRUE) &&
    (length(value) %in% c(1L, n_rows)) && (!whole || is_finite_whole(value))
  if (!ok) {
    cier_abort(
      "cier_error_input",
      "{.arg {arg}} must be a single value or one per row in [{lo}, {hi}].",
      data = list(arg = arg), call = call
    )
  }
  rep_len(value, n_rows)
}

# The common achievable range [max(min), min(max)] -- the values a single number
# can take on EVERY item. Resolve the straightline `value`: NULL draws one whole
# number per row from that range; an explicit value (scalar or per row) must be a
# whole number in it, so a number outside any item's range is a typed error rather
# than silently clipped.
sim_resolve_value <- function(value, n_rows, mins, maxs, call) {
  lo <- max(mins)
  hi <- min(maxs)
  if (lo > hi) {
    cier_abort(
      "cier_error_input",
      c("No single value lies in every item's range.",
        "x" = "The achievable common range [{lo}, {hi}] is empty.",
        "i" = "Use {.code anchor = \"position\"} on a heterogeneous battery."),
      data = list(arg = "value"), call = call
    )
  }
  if (is.null(value)) {
    return(lo + sample.int(hi - lo + 1L, n_rows, replace = TRUE) - 1L)
  }
  ok <- is_finite_whole(value) && all(value >= lo) && all(value <= hi) &&
    (length(value) %in% c(1L, n_rows))
  if (!ok) {
    cier_abort(
      "cier_error_input",
      c("{.arg value} must be whole number(s) in every item's range [{lo}, {hi}].",
        "i" = "A value above the smallest {.field max} or below the largest \\
               {.field min} cannot be a straightline value for the battery."),
      data = list(arg = "value"), call = call
    )
  }
  as.integer(rep_len(value, n_rows))
}

# The n_rows x p switched-state grid for the straightline toggle: column 1 is
# never switched, and each later column flips when an independent Bernoulli toggle
# fires, so the held position PERSISTS between switches (a cumulative parity, not a
# per-column coin). switch_prob = 1 alternates F, T, F, T ...; 0 never switches.
sim_switch_states <- function(n_rows, p, switch_prob) {
  if (p == 1L || switch_prob == 0) {
    return(matrix(FALSE, n_rows, p))
  }
  toggles <- matrix(stats::rbinom(n_rows * (p - 1L), 1L, switch_prob),
                    n_rows, p - 1L)
  parity <- matrixStats::rowCumsums(toggles) %% 2L == 1L
  cbind(matrix(FALSE, n_rows, 1L), parity)
}

# ---- Block producers --------------------------------------------------------

# random: independent uniform draws on min_j..max_j per column.
sim_block_random <- function(n_rows, mins, maxs, params = list(),
                             call = rlang::caller_env()) {
  p <- length(mins)
  cats <- maxs - mins + 1L
  out <- matrix(0L, n_rows, p)
  for (j in seq_len(p)) {
    out[, j] <- sample.int(cats[[j]], n_rows, replace = TRUE) - 1L + mins[[j]]
  }
  storage.mode(out) <- "integer"
  out
}

# straightline: each row holds one relative position (anchor = "position", a q in
# [0, 1] mapped per item) or one constant value (anchor = "value"). switch_prob
# toggles to the mirrored response (position q -> 1 - q, value -> (min + max) - v)
# and stays there until the next switch, producing longstring runs.
sim_block_straightline <- function(n_rows, mins, maxs, params = list(),
                                   call = rlang::caller_env()) {
  anchor <- params$anchor %||% "position"
  check_choice(anchor, "anchor", c("position", "value"), call = call)
  switch_prob <- params$switch_prob %||% 0
  check_number(switch_prob, "switch_prob", lower = 0, upper = 1, call = call)
  p <- length(mins)
  if (identical(anchor, "position")) {
    q <- sim_resolve_row_param(params$q, n_rows, stats::runif,
                               0, 1, FALSE, "q", call)
    base <- sim_position_grid(q, mins, maxs)
    refl <- sim_position_grid(1 - q, mins, maxs)
  } else {
    v <- sim_resolve_value(params$value, n_rows, mins, maxs, call)
    base <- matrix(v, n_rows, p)
    refl <- matrix(mins + maxs, n_rows, p, byrow = TRUE) - base
  }
  switched <- sim_switch_states(n_rows, p, switch_prob)
  base[switched] <- refl[switched]
  storage.mode(base) <- "integer"
  base
}

# midpoint: the q = 0.5 position per item, optionally jittered by an integer
# uniform on [-jitter, +jitter] and clipped back into min_j..max_j.
sim_block_midpoint <- function(n_rows, mins, maxs, params = list(),
                               call = rlang::caller_env()) {
  jitter <- params$jitter %||% 0L
  if (!is_finite_whole(jitter) || length(jitter) != 1L || jitter < 0) {
    cier_abort("cier_error_input",
               "{.arg jitter} must be a single non-negative whole number.",
               data = list(arg = "jitter"), call = call)
  }
  jitter <- as.integer(jitter)
  p <- length(mins)
  mid <- sim_position_grid(rep(0.5, n_rows), mins, maxs)
  if (jitter > 0L) {
    delta <- matrix(sample.int(2L * jitter + 1L, n_rows * p, replace = TRUE) -
                      jitter - 1L, n_rows, p)
    mn <- matrix(mins, n_rows, p, byrow = TRUE)
    mx <- matrix(maxs, n_rows, p, byrow = TRUE)
    mid <- pmin(pmax(mid + delta, mn), mx)
  }
  storage.mode(mid) <- "integer"
  mid
}

# extreme: each cell is the per-column maximum with probability p_high, else the
# minimum -- independent across cells.
sim_block_extreme <- function(n_rows, mins, maxs, params = list(),
                              call = rlang::caller_env()) {
  p_high <- params$p_high %||% 0.5
  check_number(p_high, "p_high", lower = 0, upper = 1, call = call)
  p <- length(mins)
  hi <- matrix(stats::rbinom(n_rows * p, 1L, p_high), n_rows, p)
  mn <- matrix(mins, n_rows, p, byrow = TRUE)
  mx <- matrix(maxs, n_rows, p, byrow = TRUE)
  out <- hi * mx + (1L - hi) * mn
  storage.mode(out) <- "integer"
  out
}

# diagonal: a cyclic ramp ((start - 1) + step * (j - 1)) mod K + 1, or a bounce
# (triangle-wave) ramp that climbs to K and back. K = K_min over the span; the
# position is offset to min_j + pos - 1. start defaults to a random position.
sim_block_diagonal <- function(n_rows, mins, maxs, params = list(),
                               call = rlang::caller_env()) {
  step <- params$step %||% 1L
  if (!is_finite_whole(step) || length(step) != 1L || step < 1) {
    cier_abort("cier_error_input",
               "{.arg step} must be a single positive whole number.",
               data = list(arg = "step"), call = call)
  }
  step <- as.integer(step)
  bounce <- params$bounce %||% FALSE
  check_flag(bounce, "bounce", call = call)
  k <- min(maxs - mins + 1L)
  if (k < 2L) {
    cier_abort("cier_error_input",
               "A diagonal needs at least two response categories.",
               data = list(arg = "items"), call = call)
  }
  # The walk repeats every `period`: K_min for the cyclic ramp, 2*(K_min - 1) for
  # the bounce triangle wave. A step that is a multiple of the period freezes the
  # position, collapsing the diagonal to a constant column (a mislabelled
  # straightliner), so it is rejected in both modes.
  period <- if (bounce) 2L * (k - 1L) else k
  if (step %% period == 0L) {
    cier_abort(
      "cier_error_input",
      c("{.arg step} must not be a multiple of the diagonal period ({period}).",
        "i" = "Otherwise the diagonal collapses to a constant column."),
      data = list(arg = "step"), call = call
    )
  }
  start <- sim_resolve_row_param(
    params$start, n_rows, function(n) sample.int(k, n, replace = TRUE),
    1L, k, TRUE, "start", call
  )
  p <- length(mins)
  raw <- matrix(start - 1L, n_rows, p) +
    step * matrix(seq_len(p) - 1L, n_rows, p, byrow = TRUE)
  pos <- if (bounce) {
    m <- raw %% period
    ifelse(m < k, m + 1L, 2L * k - 1L - m)
  } else {
    raw %% k + 1L
  }
  out <- matrix(mins, n_rows, p, byrow = TRUE) + pos - 1L
  storage.mode(out) <- "integer"
  out
}

# The default alternating positions: `period` points spread over 1..K_min (the
# classic high/low seesaw c(1, K) at period 2), or a cyclic 1..K_min run when the
# period exceeds K_min.
sim_alternating_default <- function(period, k) {
  if (k >= period) {
    return(as.integer(sim_round_half_up(seq.int(1, k, length.out = period))))
  }
  rep_len(seq_len(k), period)
}

# alternating: a cyclic walk through `values` (positions in 1..K_min) at the given
# period, offset to min_j + pos - 1. period-d distinct values peak the
# autocorrelation at lag d.
sim_block_alternating <- function(n_rows, mins, maxs, params = list(),
                                  call = rlang::caller_env()) {
  period <- params$period %||% 2L
  if (!is_finite_whole(period) || length(period) != 1L || period < 2) {
    cier_abort("cier_error_input",
               "{.arg period} must be a single whole number of at least 2.",
               data = list(arg = "period"), call = call)
  }
  period <- as.integer(period)
  k <- min(maxs - mins + 1L)
  values <- params$values
  if (is.null(values)) {
    values <- sim_alternating_default(period, k)
  } else {
    ok <- is_finite_whole(values) && length(values) == period &&
      all(values >= 1 & values <= k) && length(unique(values)) >= 2L
    if (!ok) {
      cier_abort(
        "cier_error_input",
        c("{.arg values} must be {period} positions in 1..K_min = {k}, with at \\
           least two distinct.",
          "x" = "Got {length(values)} value{?s}."),
        data = list(arg = "values"), call = call
      )
    }
    values <- as.integer(values)
  }
  offset <- sim_resolve_row_param(
    params$start_offset, n_rows,
    function(n) sample.int(period, n, replace = TRUE) - 1L,
    0L, period - 1L, TRUE, "start_offset", call
  )
  p <- length(mins)
  idx <- (matrix(offset, n_rows, p) +
            matrix(seq_len(p) - 1L, n_rows, p, byrow = TRUE)) %% period + 1L
  out <- matrix(mins, n_rows, p, byrow = TRUE) +
    matrix(values[idx], n_rows, p) - 1L
  storage.mode(out) <- "integer"
  out
}

# ---- Markov ------------------------------------------------------------------

# Default transition over K positions: 0.6 of staying, the rest spread uniformly
# -- a predictable but non-degenerate chain (the matched generator for cier_lazr).
sim_markov_default_transition <- function(k) {
  off <- (1 - 0.6) / (k - 1L)
  m <- matrix(off, k, k)
  diag(m) <- 0.6
  m
}

# The stationary distribution of a transition matrix (the left eigenvector for
# eigenvalue 1), falling back to uniform when it is not uniquely real-positive.
sim_markov_stationary <- function(transition) {
  k <- nrow(transition)
  eig <- eigen(t(transition))
  idx <- which(abs(Re(eig$values) - 1) <= 1e-8 & abs(Im(eig$values)) <= 1e-8)
  if (length(idx) != 1L) {
    return(rep(1 / k, k))
  }
  vec <- Re(eig$vectors[, idx])
  if (sum(vec) < 0) vec <- -vec
  if (any(vec < -1e-8) || !is.finite(sum(vec)) ||
        sum(vec) <= .Machine$double.eps) {
    return(rep(1 / k, k))
  }
  vec <- pmax(vec, 0)
  vec / sum(vec)
}

sim_markov_validate_transition <- function(transition, k, call) {
  if (!is.matrix(transition) || !is.numeric(transition) ||
        nrow(transition) != k || ncol(transition) != k) {
    cier_abort("cier_error_input",
               "{.arg transition} must be a numeric {k} x {k} matrix (K_min = {k}).",
               data = list(arg = "transition"), call = call)
  }
  if (any(!is.finite(transition)) || any(transition < -1e-8) ||
        any(transition > 1 + 1e-8)) {
    cier_abort("cier_error_input",
               "{.arg transition} entries must be finite and lie in [0, 1].",
               data = list(arg = "transition"), call = call)
  }
  if (any(abs(rowSums(transition) - 1) > 1e-8)) {
    cier_abort("cier_error_input",
               "{.arg transition} rows must sum to 1.",
               data = list(arg = "transition"), call = call)
  }
  invisible(NULL)
}

sim_markov_validate_initial <- function(initial, k, call) {
  ok <- checkmate::test_numeric(initial, len = k, lower = -1e-8,
                                upper = 1 + 1e-8, any.missing = FALSE,
                                finite = TRUE) &&
    abs(sum(initial) - 1) <= 1e-8
  if (!ok) {
    cier_abort("cier_error_input",
               "{.arg initial} must be {k} probabilities in [0, 1] summing to 1.",
               data = list(arg = "initial"), call = call)
  }
  invisible(NULL)
}

# markov: a first-order chain over positions 1..K_min governed by `transition`
# (the first item drawn from `initial`), offset to min_j + state - 1.
sim_block_markov <- function(n_rows, mins, maxs, params = list(),
                             call = rlang::caller_env()) {
  k <- min(maxs - mins + 1L)
  transition <- params$transition
  if (is.null(transition)) {
    transition <- sim_markov_default_transition(k)
  } else {
    sim_markov_validate_transition(transition, k, call)
  }
  initial <- params$initial
  if (is.null(initial)) {
    initial <- sim_markov_stationary(transition)
  } else {
    sim_markov_validate_initial(initial, k, call)
  }
  p <- length(mins)
  trans_cdf <- t(apply(transition, 1L, cumsum))
  trans_cdf[, k] <- 1
  init_cdf <- cumsum(initial)
  init_cdf[k] <- 1
  states <- matrix(0L, n_rows, p)
  current <- findInterval(stats::runif(n_rows), init_cdf) + 1L
  states[, 1L] <- current
  for (j in seq_len(p)[-1L]) {
    current <- max.col(stats::runif(n_rows) < trans_cdf[current, , drop = FALSE],
                       ties.method = "first")
    states[, j] <- current
  }
  out <- matrix(mins, n_rows, p, byrow = TRUE) + states - 1L
  storage.mode(out) <- "integer"
  out
}

# ---- Registry / allowlist ---------------------------------------------------

# The eight recognised pattern NAMES (the weights allowlist). speeder is content
# attentive (its carelessness is timing only, a later layer), so it is in the
# allowlist but not the content registry below.
sim_pattern_names <- function() {
  c("random", "straightline", "midpoint", "extreme", "diagonal",
    "alternating", "markov", "speeder")
}

# The seven content producers, name -> block function. speeder is absent (its rows
# keep their attentive content), so sim_block_fun("speeder") is NULL.
sim_pattern_registry <- function() {
  list(
    random       = sim_block_random,
    straightline = sim_block_straightline,
    midpoint     = sim_block_midpoint,
    extreme      = sim_block_extreme,
    diagonal     = sim_block_diagonal,
    alternating  = sim_block_alternating,
    markov       = sim_block_markov
  )
}

sim_block_fun <- function(name) {
  sim_pattern_registry()[[name]]
}
