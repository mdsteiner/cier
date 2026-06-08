# Purpose: PerFit-backed kernels for the nonparametric person-fit indices. The
#          normed polytomous Guttman-error index (Gnormed) lives here; the
#          person-scalability index (Ht) will join it. Unlike the pure indirect
#          kernels (R/index-kernels-indirect.R), a bridge kernel here MAY raise a
#          typed condition on a data/contract violation (an out-of-range
#          zero-base) and needs the fitted PerFit object to resolve its
#          Monte-Carlo null cutoff -- the documented exception to the value-only
#          resolve_cutoff() path (see ADR.md, "Gnormed cutoff: PerFit Monte-Carlo
#          null").
# Args:    See per-function documentation.
# Returns: Documented per function.
# Invariants:
#   - The single production implementation of Gnormed is PerFit::Gnormed.poly
#     (single-kernel rule); the closed-form re-derivation lives only in
#     tests/reference/ref-personfit-niessen-2016.R.
#   - Responses are reverse-keyed (in the wrapper) then recoded to PerFit's
#     0..(Ncat-1) contract here; an out-of-range cell is a typed
#     cier_error_input, never silent coercion.
#   - The Monte-Carlo null cutoff is reproducible: a non-NULL seed is applied
#     LOCALLY and the caller's .Random.seed is restored on exit.

# Recode a (reverse-scored) complete-case response block to PerFit's 0..(Ncat-1)
# contract by subtracting the per-item scale base `mins` (default 1) and check the
# two preconditions PerFit::Gnormed.poly enforces (its Sanity.dma.poly): every
# cell must lie in 0..(Ncat-1), AND the block must ATTAIN both extremes (its
# popularity estimates are undefined if the lowest or highest category never
# occurs). The caller passes a complete-case block, so there are no NAs. Either
# violation is a typed cier_error_input -- converting PerFit's terse abort into a
# package condition -- rather than a silent shift or a cryptic upstream stop().
personfit_zero_base <- function(block, mins, ncat, call) {
  m <- sweep(block, 2L, mins)
  # Gnormed scores categorical responses; a fractional cell is a contract
  # violation. check_responses() only rejects NaN/Inf, and the integer cast below
  # would TRUNCATE a fractional value silently (2.5 -> 2) -- and PerFit, which
  # would itself reject a non-integer matrix, never sees the original. Catch it
  # here as a typed error instead of a silent wrong score. (Reverse-keying and
  # zero-basing are integer-exact, so this only fires on genuinely fractional
  # input, never on floating-point noise.)
  if (any(m != round(m))) {
    cier_abort(
      "cier_error_input",
      c("{.arg responses} must be whole-number category codes for Gnormed.",
        "x" = "Found a fractional response value.",
        "i" = "Gnormed scores categorical responses; recode or drop non-integer \\
               cells (it does not accept averaged or imputed fractional scores)."),
      data = list(arg = "responses"), call = call
    )
  }
  rng <- range(m)
  if (rng[[1L]] < 0L || rng[[2L]] > ncat - 1L) {
    cier_abort(
      "cier_error_input",
      c("After zero-basing, every response must lie in {.val 0..(Ncat - 1)}.",
        "x" = "Observed zero-based range: {.val {rng}}.",
        "i" = "Ncat = {ncat}; check {.field categories} / {.field min} against \\
               the data."),
      data = list(arg = "responses", observed = rng, ncat = ncat),
      call = call
    )
  }
  if (rng[[1L]] > 0L || rng[[2L]] < ncat - 1L) {
    cier_abort(
      "cier_error_input",
      c("The responses must use the full declared scale: the lowest and highest \\
         of the {ncat} categories must each occur at least once.",
        "x" = "Observed zero-based range: {.val {rng}} (need {.val {0L}} and \\
               {.val {ncat - 1L}}).",
        "i" = "Check {.field categories} / {.field min}, or that the sample spans \\
               every category (PerFit cannot score a scale with an unused end)."),
      data = list(arg = "responses", observed = rng, ncat = ncat),
      call = call
    )
  }
  storage.mode(m) <- "integer"
  m
}

# Per-respondent Gnormed via PerFit::Gnormed.poly on the complete-case block.
# Returns list(value, fit): `value` is the full-length score vector (NA for any
# respondent with a missing cell, and all-NA when fewer than two respondents are
# complete or fewer than three items are present -- PerFit's recursive
# denominator needs >= 3 items); `fit` is the fitted PerFit object for the
# complete block (NULL when the index abstains wholesale), reused by the
# Monte-Carlo null cutoff so the simulation is not refitted. `responses` has
# already been reverse-keyed by the wrapper. PerFit's printed progress is
# captured so only typed cli conditions reach the user.
kernel_gnormed <- function(responses, categories, mins, call) {
  n <- nrow(responses)
  value <- rep(NA_real_, n)
  if (ncol(responses) < 3L) {
    return(list(value = value, fit = NULL))
  }
  complete <- stats::complete.cases(responses)
  if (sum(complete) < 2L) {
    return(list(value = value, fit = NULL))
  }
  ncat <- as.integer(categories[[1L]])
  z <- personfit_zero_base(responses[complete, , drop = FALSE], mins, ncat, call)
  utils::capture.output(
    fit <- PerFit::Gnormed.poly(matrix = z, Ncat = ncat)
  )
  value[complete] <- as.numeric(fit$PFscores$PFscores)
  list(value = value, fit = fit)
}

# The PerFit Monte-Carlo null cutoff for a fitted person-fit object: simulate
# model-conforming response vectors and take the nominal-FPR (`Blvl`) quantile of
# the statistic (PerFit::cutoff). A non-NULL `seed` makes the simulation
# reproducible and is applied LOCALLY -- the caller's global .Random.seed is saved
# and restored on exit (mirroring kernel_rpr) -- so a seeded call does not disturb
# the session RNG; a NULL seed draws from the ambient stream (the cutoff then
# varies per call). PerFit's printed summary is captured.
resolve_perfit_null_cutoff <- function(fit, fpr, seed) {
  if (!is.null(seed)) {
    saved <- globalenv()[[".Random.seed"]]   # NULL when no RNG drawn yet
    on.exit(restore_random_seed(saved), add = TRUE)
    set.seed(seed)
  }
  utils::capture.output(co <- PerFit::cutoff(fit, Blvl = fpr))
  as.numeric(co$Cutoff)
}

# Resolve the flagging cutoff for cier_gnormed: a literal `cutoff` passes through
# verbatim (already validated by the wrapper); otherwise the PerFit Monte-Carlo
# null at the `fpr` level (NULL `fpr` uses the registry default). When the index
# abstained wholesale (no fitted object) there is no null to build, so the cutoff
# is NA -- flagging no one -- with the shared insufficient-items warning.
resolve_gnormed_cutoff <- function(fit, fpr, cutoff, default_fpr, seed,
                                   call = rlang::caller_env()) {
  if (!is.null(cutoff)) {
    return(cutoff)
  }
  if (is.null(fit)) {
    cier_warn(
      "cier_warning_insufficient_items",
      c("Cannot resolve a cutoff: no respondent could be scored.",
        "i" = "Returning {.val NA} as the cutoff."),
      data = list(n_used = 0L), call = call
    )
    return(NA_real_)
  }
  resolve_perfit_null_cutoff(fit, if (is.null(fpr)) default_fpr else fpr, seed)
}
