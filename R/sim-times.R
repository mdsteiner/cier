# Internal response-time kernels for cier_simulate(): a per-cell lognormal model (van der
# Linden 2007 family) with a per-respondent pace intercept, plus page resolution/aggregation
# feeding cier_total_time() and cier_page_time(). The careless log-mean applies only inside a
# careless row's inclusive [onset, offset], so partial/temporary rows speed up exactly where
# their content mutates.
#
# Frozen draw order: the n * p cell noise first (column-major), then the n pace intercepts.
# Changing it reorders every downstream draw for a given seed. Output is finite and
# >= min_seconds > 0, so totals pass check_seconds() / check_page_seconds(). The pace
# intercept makes attentive and careless time distributions overlap, by design.

# ---- Timing parameters ------------------------------------------------------

# Timing calibration (single source of truth): attentive median 8 s/item, careless 1.5 s/item
# -- below the 2 s/item per-page floor of Bowling et al. (2023) -- cell-noise SD 0.5, pace SD
# 1.2. Held by the slow-tier acceptance test in test-sim-times.R.
sim_timing_defaults <- function() {
  list(mu_att = log(8), mu_car = log(1.5), sigma = 0.5,
       respondent_sd = 1.2, min_seconds = 0.1)
}

# Merge a partial `timing` override over the defaults and validate every entry. Unknown keys
# are typed errors so a typo cannot silently fall back to a default.
sim_resolve_timing <- function(timing, call = rlang::caller_env()) {
  if (!is.list(timing)) {
    cier_abort("cier_error_input",
               "{.arg timing} must be a named list of timing parameters.",
               data = list(arg = "timing"), call = call)
  }
  defaults <- sim_timing_defaults()
  nms <- names(timing)
  if (length(timing) > 0L) {
    check_unique_names(nms, "timing", call)
  }
  reject_unknown_keys(nms, names(defaults), "timing", "timing parameter", call)
  out <- utils::modifyList(defaults, timing)
  check_number(out$mu_att, "timing$mu_att", call = call)
  check_number(out$mu_car, "timing$mu_car", call = call)
  check_number(out$sigma, "timing$sigma", lower = 0, call = call)
  check_number(out$respondent_sd, "timing$respondent_sd", lower = 0, call = call)
  check_number(out$min_seconds, "timing$min_seconds", call = call)
  # strictly positive, not just non-negative: the floor must keep every duration past
  # check_seconds(), which rejects zero seconds.
  if (out$min_seconds <= 0) {
    cier_abort("cier_error_input",
               "{.arg timing$min_seconds} must be greater than 0.",
               data = list(arg = "timing$min_seconds"), call = call)
  }
  out
}

# ---- Careless-plan validation -----------------------------------------------

# Validate the careless plan: `careless` flags the shifted rows; `onset` / `offset` give each
# careless row's inclusive span (attentive rows carry NA and are never read).
sim_times_validate_plan <- function(careless, onset, offset, n, p, call) {
  if (!is.logical(careless) || length(careless) != n || anyNA(careless)) {
    cier_abort(
      "cier_error_input",
      "{.arg careless} must be a logical vector of length {n} without NA.",
      data = list(arg = "careless"), call = call
    )
  }
  if (length(onset) != n || length(offset) != n) {
    cier_abort(
      "cier_error_input",
      "{.arg onset} and {.arg offset} must each have length {n}.",
      data = list(arg = "onset"), call = call
    )
  }
  sim_times_validate_spans(onset[careless], offset[careless], p, call)
}

# The careless rows' span bounds: whole numbers 1 <= onset <= offset <= p.
sim_times_validate_spans <- function(on_c, off_c, p, call) {
  spans_ok <- length(on_c) == 0L ||
    (is_finite_whole(on_c) && is_finite_whole(off_c) &&
       all(on_c >= 1) && all(off_c <= p) && all(on_c <= off_c))
  if (!spans_ok) {
    cier_abort(
      "cier_error_input",
      c("Careless rows need a whole-number span 1 <= onset <= offset <= {p}.",
        "i" = "{.arg onset} / {.arg offset} may be NA only on attentive rows."),
      data = list(arg = "onset"), call = call
    )
  }
  invisible(NULL)
}

# ---- Lognormal cell times ---------------------------------------------------

# Per-cell response times: exp(noise + mu + pace), floored at min_seconds. noise ~ N(0,
# sigma) per cell; pace ~ N(0, respondent_sd) per respondent, shared across that respondent's
# items so the group distributions overlap; mu is mu_car on a careless row's inclusive span,
# mu_att elsewhere. `timing` is a resolved sim_resolve_timing() list.
sim_times_lognormal <- function(n, p, careless, onset, offset,
                                timing = sim_timing_defaults(),
                                call = rlang::caller_env()) {
  check_count(n, "n", call = call)
  check_count(p, "p", call = call)
  sim_times_validate_plan(careless, onset, offset, n, p, call)
  noise <- matrix(stats::rnorm(n * p, sd = timing$sigma), nrow = n, ncol = p)
  mu <- matrix(timing$mu_att, nrow = n, ncol = p)
  for (i in which(careless)) {
    mu[i, onset[[i]]:offset[[i]]] <- timing$mu_car
  }
  pace <- stats::rnorm(n, sd = timing$respondent_sd)
  # `+ pace` recycles the length-n vector down columns, adding pace[i] to every cell of row i.
  pmax(exp(noise + mu + pace), timing$min_seconds)
}

# ---- Page structure ---------------------------------------------------------

# Resolve per-page item counts. NULL: one page per contiguous run of items$scale. Scalar:
# uniform pages of that size, last takes the remainder. Vector of positive whole counts
# summing to p: used verbatim.
sim_resolve_pages <- function(items_per_page, scale,
                              call = rlang::caller_env()) {
  p <- length(scale)
  if (is.null(items_per_page)) {
    return(as.integer(rle(scale)$lengths))
  }
  ok_shape <- is.numeric(items_per_page) && is.null(dim(items_per_page)) &&
    length(items_per_page) >= 1L && is_finite_whole(items_per_page) &&
    all(items_per_page >= 1)
  if (!ok_shape) {
    cier_abort(
      "cier_error_input",
      "{.arg items_per_page} must be positive whole number(s).",
      data = list(arg = "items_per_page"), call = call
    )
  }
  ipp <- as.integer(items_per_page)
  if (length(ipp) == 1L) {
    if (ipp > p) {
      cier_abort(
        "cier_error_input",
        c("{.arg items_per_page} must be at most the item count ({p}).",
          "x" = "Got {ipp}."),
        data = list(arg = "items_per_page", observed = ipp), call = call
      )
    }
    remainder <- p %% ipp
    return(as.integer(c(rep(ipp, p %/% ipp), if (remainder > 0L) remainder)))
  }
  if (sum(ipp) != p) {
    cier_abort(
      "cier_error_input",
      c("{.arg items_per_page} must sum to the item count ({p}).",
        "x" = "Got {sum(ipp)} over {length(ipp)} page{?s}."),
      data = list(arg = "items_per_page", observed = sum(ipp)), call = call
    )
  }
  ipp
}

# Sum each respondent's cell times into per-page totals over consecutive page blocks (page g
# covers the items_per_page[g] items after page g - 1's last item).
sim_page_totals <- function(cell_times, items_per_page) {
  ends <- cumsum(items_per_page)
  starts <- ends - items_per_page + 1L
  out <- matrix(NA_real_, nrow = nrow(cell_times),
                ncol = length(items_per_page))
  for (g in seq_along(items_per_page)) {
    out[, g] <- matrixStats::rowSums2(cell_times,
                                      cols = starts[[g]]:ends[[g]])
  }
  out
}
