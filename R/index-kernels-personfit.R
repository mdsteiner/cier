# Kernels for the nonparametric person-fit indices: the normed polytomous Guttman-error
# index (Gnormed) and the polytomous person-scalability index (Ht). Unlike the pure
# indirect kernels, a kernel here may raise a typed condition on a data/contract violation
# (out-of-range zero-base, fractional response). Gnormed defaults to a Monte-Carlo null
# cutoff (the exception to the value-only resolve_cutoff() path); Ht defaults to the
# percentile path but can opt in to the same Monte-Carlo null (method = "mc_null").
#
# The Monte-Carlo null engine (personfit_null_matrix / bootstrap_tail_cutoff /
# personfit_null_cutoff) is shared by both indices via resolve_personfit_null_cutoff,
# parameterised by scoring kernel, tail direction, and perfect-vector rule: Ht (ht_scores,
# lower tail, perfect vectors excluded) and Gnormed (gnormed_scores, upper tail, perfect
# vectors allowed) draw from the same sum-score-conditional generator.

# Assert that a (reverse-keyed) response matrix carries whole-number category codes. The
# kernels cast to integer before scoring, so a fractional cell (an averaged or imputed
# value) would be silently truncated into a wrong score. Raised on the full matrix (NA cells
# ignored) before any abstention short-circuit; reverse-keying and zero-basing are
# integer-exact, so this fires only on genuinely fractional input. `statistic` names the
# index in the message.
assert_integer_responses <- function(responses, statistic, call) {
  if (any(responses != round(responses), na.rm = TRUE)) {
    cier_abort(
      "cier_error_input",
      c("{.arg responses} must be whole-number category codes for {statistic}.",
        "x" = "Found a fractional response value.",
        "i" = "{statistic} scores categorical responses; recode or drop \\
               non-integer cells (it does not accept averaged or imputed \\
               fractional scores)."),
      data = list(arg = "responses"), call = call
    )
  }
  invisible(responses)
}

# Recode a (reverse-scored) complete-case block to the 0..(Ncat-1) contract by subtracting
# the per-item scale base `mins`, and check two preconditions: every cell in 0..(Ncat-1),
# and the block must attain both extremes (popularity estimates are undefined if the lowest
# or highest category never occurs). Violations are classified differently: an out-of-range
# cell is a data defect (plain cier_error_input); an unused scale extreme is otherwise-valid
# data the closed form cannot score (additionally tagged cier_error_backend_limit, so the
# screen skips-with-reason).
personfit_zero_base <- function(block, mins, ncat, call) {
  m <- sweep(block, 2L, mins)
  rng <- range(m)
  if (rng[[1L]] < 0L || rng[[2L]] > ncat - 1L) {
    cier_abort(
      "cier_error_input",
      c("After zero-basing, every response must lie in {.val 0..(Ncat - 1)}.",
        "x" = "Observed zero-based range: {.val {rng}}.",
        "i" = "Ncat = {ncat}; check {.field max} / {.field min} against \\
               the data."),
      data = list(arg = "responses", observed = rng, ncat = ncat),
      call = call
    )
  }
  if (rng[[1L]] > 0L || rng[[2L]] < ncat - 1L) {
    # An unused extreme is a property of an otherwise-valid sample the closed form cannot
    # score, so it carries the cier_error_backend_limit subclass with a data$reason, letting
    # cier_screen() record Gnormed as skipped-with-reason instead of aborting the battery; a
    # direct cier_gnormed() call still gets the cier_error_input with the remedy.
    cier_abort(
      c("cier_error_backend_limit", "cier_error_input"),
      c("The responses must use the full declared scale: the lowest and highest \\
         of the {ncat} categories must each occur at least once.",
        "x" = "Observed zero-based range: {.val {rng}} (need {.val {0L}} and \\
               {.val {ncat - 1L}}).",
        "i" = "Respondents with any missing cell are dropped before scoring, so \\
               a category present only in incomplete rows does not count toward \\
               the observed span.",
        "i" = "Check {.field max} / {.field min}, or that the sample spans \\
               every category (the item-step popularities are undefined when a \\
               scale end is never used)."),
      data = list(
        arg = "responses", observed = rng, ncat = ncat,
        reason = paste0(
          "sample does not attain both scale extremes ",
          "(every end category must be observed)"
        )
      ),
      call = call
    )
  }
  storage.mode(m) <- "integer"
  m
}

# The Gnormed closed form on a complete, zero-based block `m` (rows = respondents, values
# 0..(ncat - 1); columns = items). Returns the per-respondent normed polytomous
# Guttman-error vector; the single scorer for both the cier_gnormed value and the
# Monte-Carlo null's per-replicate scoring.
#
# Construction (Molenaar 1991; Emons 2008): expand each item into its (ncat - 1) steps and
# give every item-step a popularity position (1 = most popular), ranked over all item-steps.
# For a respondent endorsing NC steps,
#   Gp   = sum(positions of endorsed steps) - NC(NC + 1)/2
# (the second term is the error-free Guttman position-sum; each excess unit is one
# inversion), normed by
#   maxG = max position-sum over allocations of NC steps - NC(NC + 1)/2,
# a max-plus knapsack over items (each contributes a 0..M-step prefix). The knapsack is
# n-independent (depends only on the popularity ordering), so it is computed once and indexed
# by each respondent's NC. Gnormed = Gp / maxG, with 0/0 -> 0 (all-min / all-max rows). `m`
# must be complete with whole-number codes in 0..(ncat - 1); the caller guarantees this.
gnormed_scores <- function(m, ncat) {
  storage.mode(m) <- "double"
  n <- nrow(m)
  p <- ncol(m)
  mm <- ncat - 1L
  # Item-step popularities P(X_j >= s), s = 1..mm, as a p x mm matrix.
  probs <- vapply(seq_len(mm), function(s) colMeans(m >= s), numeric(p))
  if (mm == 1L) probs <- matrix(probs, ncol = 1L)
  # Popularity positions (1 = most popular) over all item-steps, item-major; ties by first
  # appearance.
  pos_vect <- rank(-as.vector(t(probs)), ties.method = "first")
  posmat <- matrix(pos_vect, nrow = p, ncol = mm, byrow = TRUE)
  cs <- apply(posmat, 1L, cumsum)            # mm x p (apply transposes the result)
  if (mm == 1L) cs <- matrix(cs, nrow = 1L)
  cumrank <- cbind(0, t(cs))                 # p x (mm + 1): column k + 1 = steps 0..k
  nc <- rowSums(m)
  # Per-respondent endorsed-position sum: gather cumrank[j, m[r, j] + 1] for every (r, j) in
  # one matrix-index, then row-sum. The gather is item-major, so reshaping n x p recovers
  # contrib[r, j].
  gathered <- cumrank[cbind(rep(seq_len(p), each = n), as.vector(m) + 1L)]
  sumpos <- rowSums(matrix(gathered, nrow = n, ncol = p))
  num <- sumpos - nc * (nc + 1) / 2
  # T[k + 1] = max position-sum over allocations of k steps (max-plus knapsack, each item a
  # 0..mm-step prefix). N-independent; computed once.
  total <- p * mm
  dp <- rep(-Inf, total + 1L)
  dp[[1L]] <- 0
  for (it in seq_len(p)) {
    cr <- cumrank[it, ]
    ndp <- rep(-Inf, total + 1L)
    # Max-plus convolution of dp with this item's cumulative step positions, vectorised over
    # the step count: taking s steps shifts dp by s and adds cr[s + 1], ndp keeps the
    # elementwise max. O(p * mm) vector ops instead of the O(p * total * mm) scalar loop.
    for (s in 0:mm) {
      src <- seq_len(total + 1L - s)
      ndp[src + s] <- pmax(ndp[src + s], dp[src] + cr[[s + 1L]])
    }
    dp <- ndp
  }
  max_g <- dp - (0:total) * (1:(total + 1L)) / 2
  res <- num / max_g[nc + 1L]
  res[is.nan(res)] <- 0
  res
}

# Per-respondent Gnormed on the complete-case block, scored by gnormed_scores. `responses`
# is already reverse-keyed by the wrapper; `ncat` is the validated shared category count,
# `mins` the validated per-item base. Returns list(value, z, ncat, abstain): `value` is the
# full-length score vector (NA per missing-cell respondent, all-NA when fewer than two
# respondents complete or fewer than three items -- the recursive Guttman denominator needs
# >= 3 items); `z` is the validated zero-based complete block (NULL on wholesale abstention),
# passed to the Monte-Carlo null cutoff so it regenerates the null from the same data;
# `abstain` names the wholesale-abstention cause and count (NULL on the scoring path).
kernel_gnormed <- function(responses, ncat, mins, call) {
  n <- nrow(responses)
  value <- rep(NA_real_, n)
  # Check the whole-number contract first, so a fractional cell surfaces a typed error even
  # when the data would otherwise abstain.
  assert_integer_responses(responses, "Gnormed", call)
  if (ncol(responses) < 3L) {
    return(list(value = value, z = NULL, ncat = ncat,
                abstain = list(kind = "few_items", n = ncol(responses))))
  }
  complete <- stats::complete.cases(responses)
  if (sum(complete) < 2L) {
    return(list(value = value, z = NULL, ncat = ncat,
                abstain = list(kind = "few_complete", n = sum(complete))))
  }
  z <- personfit_zero_base(responses[complete, , drop = FALSE], mins, ncat, call)
  value[complete] <- gnormed_scores(z, ncat)
  list(value = value, z = z, ncat = ncat, abstain = NULL)
}

# ---- Shared nonparametric Monte-Carlo null engine ---------------------------
#
# A sum-score-conditional resample scored by a person-fit kernel, summarised by a bootstrap
# tail quantile; parameterised by scoring kernel, tail direction, and perfect-vector rule so
# it serves both Gnormed and Ht's opt-in null.

# Generate an `nreps` x p null matrix conditional on the sum score. `z` is the zero-based
# complete block (codes 0..(ncat - 1)). Resample `nreps` sum-score levels from the observed
# row sums; for each drawn level, build the per-item category-frequency table among observed
# persons at that level and draw each item from it (a category never shown at that level has
# zero frequency, so it is never drawn). Conditioning on the sum score preserves the
# inter-item covariance the person-fit statistics measure, which an unconditional column
# resample would destroy. Each draw is a direct categorical sample.int rather than
# rmultinom + which: distribution-identical but cheaper.
personfit_null_matrix <- function(z, ncat, nreps = 1000L) {
  p <- ncol(z)
  nc <- rowSums(z)
  nc_gen <- sample(nc, size = nreps, replace = TRUE)
  null <- matrix(NA_integer_, nrow = nreps, ncol = p)
  for (lvl in sort(unique(nc_gen))) {
    rows <- which(nc_gen == lvl)
    obs <- z[nc == lvl, , drop = FALSE]
    for (j in seq_len(p)) {
      freq <- tabulate(obs[, j] + 1L, nbins = ncat)
      null[rows, j] <- sample.int(ncat, length(rows), replace = TRUE,
                                  prob = freq) - 1L
    }
  }
  null
}

# Bootstrap tail-quantile cutoff over a null person-fit vector `pfs`. For `breps` resamples
# of `pfs`, take the directional quantile (upper tail: 1 - fpr; lower tail: fpr) and return
# the median of those bootstrap quantiles (quantile type 7, R's default).
#
# Resamples are processed in column chunks to keep the fast matrixStats::colQuantiles path
# without materialising the full n x breps matrix on large nulls. Indexing through sample.int
# (first arg n >= 1) sidesteps sample()'s length-1 convenience trap, so a heavily-pruned null
# leaving a single score still resamples it. `chunk_cells` is internal, so tests can force
# the chunked path on tiny fixtures.
bootstrap_tail_cutoff <- function(pfs, tail = "upper", fpr = 0.05, breps = 1000L,
                                  chunk_cells = 5e6) {
  blvl_use <- if (identical(tail, "upper")) 1 - fpr else fpr
  n <- length(pfs)
  chunk_cols <- max(1L, min(breps, floor(chunk_cells / n)))
  qs <- numeric(breps)
  start <- 1L
  while (start <= breps) {
    end <- min(breps, start + chunk_cols - 1L)
    cols <- end - start + 1L
    resamples <- matrix(pfs[sample.int(n, n * cols, replace = TRUE)], nrow = n)
    qs[start:end] <- as.numeric(
      matrixStats::colQuantiles(resamples, probs = blvl_use, type = 7L)
    )
    start <- end + 1L
  }
  stats::median(qs)
}

# Compose the null generator and bootstrap cutoff into one cutoff value. `score` is the
# kernel applied to the null matrix (a closure over any extra args -- e.g. ncat for Gnormed);
# `tail` is the flag direction; `perfect` is "allowed" (Gnormed scores perfect vectors) or
# "excluded" (drop constant rows, as Ht needs), a perfect vector being a constant row (zero
# variance) tested via matrixStats::rowVars > 0. The RNG is the caller's responsibility; NA
# null scores are dropped before the bootstrap. Returns NA_real_ (a clean abstention) when
# the null is too degenerate to score.
personfit_null_cutoff <- function(z, ncat, score, tail = "upper",
                                  perfect = "allowed", fpr = 0.05,
                                  nreps = 1000L, breps = 1000L) {
  null <- personfit_null_matrix(z, ncat, nreps)
  if (identical(perfect, "excluded")) {
    null <- null[matrixStats::rowVars(null) > 0, , drop = FALSE]
  }
  # Perfect-vector exclusion (Ht) can drop so many constant null rows that fewer than two
  # non-constant vectors remain (a straightliner-dominated sample resamples mostly constant
  # vectors). A whole-sample statistic cannot be scored on such a null (empty matrix errors,
  # a single row's column totals collapse to NaN), so return NA. perfect = "allowed" (Gnormed)
  # never prunes, so this guard is inert.
  if (nrow(null) < 2L) {
    return(NA_real_)
  }
  pfs <- score(null)
  pfs <- pfs[!is.na(pfs)]
  # Backstop: if >= 2 rows survive exclusion but all score NA, nothing is left to resample
  # and the bootstrap would yield NA, so short-circuit.
  if (length(pfs) < 1L) {
    return(NA_real_)
  }
  bootstrap_tail_cutoff(pfs, tail, fpr, breps)
}

# Resolve a person-fit Monte-Carlo null cutoff: build the sum-score-conditional null from the
# fitted block `z`, score it with the same closed form the index uses, and take the `tail`
# bootstrap quantile at level `fpr`. `score`/`tail`/`perfect` are the only per-index choices
# (Gnormed: gnormed_scores / "upper" / "allowed"; Ht: ht_scores / "lower" / "excluded",
# constant vectors being unscalable). A non-NULL `seed` is applied locally via
# with_local_seed() (the caller's global .Random.seed is restored on exit); a NULL seed draws
# from the ambient stream.
resolve_personfit_null_cutoff <- function(z, ncat, fpr, seed, score, tail, perfect) {
  with_local_seed(seed, function() {
    personfit_null_cutoff(z, ncat, score = score, tail = tail,
                          perfect = perfect, fpr = fpr)
  })
}

# Emit the typed wholesale-abstention warning naming why Gnormed could not be scored: too few
# items (Guttman denominator needs >= 3) or too few complete respondents (>= 2). The class
# stays cier_warning_insufficient_items so cier_screen()'s targeted muffler suppresses it.
warn_gnormed_abstention <- function(abstain, call) {
  n <- abstain$n
  if (identical(abstain$kind, "few_items")) {
    cier_warn(
      "cier_warning_insufficient_items",
      c("Cannot resolve a Gnormed cutoff: only {n} item{?s} present.",
        "i" = "The recursive Guttman denominator needs at least three items; \\
               returning {.val NA} as the cutoff."),
      data = list(n_items = n), call = call
    )
  } else {
    cier_warn(
      "cier_warning_insufficient_items",
      c("Cannot resolve a Gnormed cutoff: only {n} complete respondent{?s}.",
        "i" = "Gnormed needs at least two complete respondents; returning \\
               {.val NA} as the cutoff."),
      data = list(n_used = n), call = call
    )
  }
}

# Resolve the flagging cutoff for cier_gnormed from the kernel result `res`: a literal
# `cutoff` passes through verbatim (already validated by the wrapper); otherwise the
# Monte-Carlo null at level `fpr` (NULL `fpr` uses the package default). On wholesale
# abstention (no scored block) there is no null to build, so the cutoff is NA -- flagging no
# one -- and the warning names the cause `res$abstain` carries.
resolve_gnormed_cutoff <- function(res, fpr, cutoff, default_fpr, seed,
                                   call = rlang::caller_env()) {
  if (!is.null(cutoff)) {
    return(cutoff)
  }
  if (is.null(res$z)) {
    warn_gnormed_abstention(res$abstain, call)
    return(NA_real_)
  }
  # Gnormed flags high misfit (upper tail), perfect-ALLOWED (straightliners score).
  resolve_personfit_null_cutoff(
    res$z, res$ncat, if (is.null(fpr)) default_fpr else fpr, seed,
    score = function(null) gnormed_scores(null, res$ncat),
    tail = "upper", perfect = "allowed"
  )
}

# Sort each respondent's responses ascending and return an unnamed double matrix. This
# per-row sort is ~97% of kernel_ht's cost: for small-integer rating-scale codes a counting
# sort avoids the comparison-sort log factor and is ~20x faster than t(apply(., sort)).
#
# PRECONDITION (caller upholds; not re-checked on the hot path): z holds whole numbers in
# 0..m, `m` IS max(z), and z has >= 2 columns. kernel_ht guarantees all three and passes the
# max it already computed.
#
# Scales wider than `max_cats` fall back to the comparison sort. The bound is on memory, not
# speed: the counting sort's cum / index structures are O(n * m), so the fallback caps that
# allocation against pathological input. 100 clears every real rating scale (incl. 0..100).
personfit_sorted_rows <- function(z, m = max(z), max_cats = 100L) {
  if (m > max_cats) {
    sorted <- t(apply(z, 1L, sort))
    dimnames(sorted) <- NULL
    storage.mode(sorted) <- "double"
    return(sorted)
  }
  n <- nrow(z)
  p <- ncol(z)
  # Difference-array counting sort, O(n * p + n * m): per-row category counts give
  # cum[, k] = #(values < k), so value k's run starts at sorted-column cum[, k] + 1 (skip
  # rows where k is absent). Scatter a +1 at each run start, then a row cumsum turns the step
  # markers back into the ascending sorted values.
  zi <- z
  storage.mode(zi) <- "integer"
  cum <- matrixStats::rowCumsums(matrixStats::rowTabulates(zi, values = 0:m))
  starts <- cum[, seq_len(m), drop = FALSE] + 1L
  w <- which(starts <= p)
  lin <- (starts[w] - 1L) * n + (w - 1L) %% n + 1L
  sorted <- matrixStats::rowCumsums(matrix(tabulate(lin, n * p), n, p))
  storage.mode(sorted) <- "double"
  sorted
}

# The Ht closed form on a complete, zero-based block `z` (rows = respondents, codes 0..m;
# columns = items). Returns the per-respondent person-scalability vector -- the Frechet /
# rearrangement collapse of the transposed Mokken scalability; the single scorer for both the
# cier_ht value and the Monte-Carlo null's per-replicate scoring. For each respondent v,
#   Ht_v = cov_p(z_v, Tot - z_v) / cov_p(sort z_v, SortedTot - sort z_v),
# where cov_p is the population covariance across the p items, Tot = colSums(z), and
# SortedTot = colSums of the per-row-sorted block (the comonotonic maximum: sort both, pair).
# This is the bilinearity collapse of sum_{w != v} cov(z_v, z_w) over the off-diagonal, so
# the cost is linear in n (one O(n * p) covariance pass plus the per-respondent counting sort)
# and never materialises the n x n person matrix, with no category ceiling. `m` IS max(z). A
# non-finite ratio -> NA (the denominator collapses to 0 on a zero-variance / straightline
# row); `z` must be complete with whole-number codes 0..m, which the caller guarantees.
ht_scores <- function(z, m = max(z)) {
  p <- ncol(z)
  # The one R-level pass (~97% of the cost): per-row ascending sort (counting sort,
  # comparison-sort fallback for pathologically wide scales).
  sorted <- personfit_sorted_rows(z, m)
  tot <- colSums(z)
  sorted_tot <- colSums(sorted)
  # Per-respondent population covariance across the p items, cov_p(x, T - x), T the per-item
  # column total. Expanded so the only cross term is a matrix-vector product -- no n x p total
  # matrix is materialised:
  #   cov_p(x, T - x) = (x %*% T)/p - E_j[x^2] - E_j[x] * (mean(T) - E_j[x]).
  # `x^2` (not x * x) so an integer x squares in double -- no integer-overflow path.
  pcov_tot <- function(x, col_tot) {
    mx <- rowMeans(x)
    as.numeric(x %*% col_tot) / p - rowMeans(x^2) - mx * (sum(col_tot) / p - mx)
  }
  v <- pcov_tot(z, tot) / pcov_tot(sorted, sorted_tot)
  v[!is.finite(v)] <- NA_real_
  v
}

# Per-respondent polytomous person-scalability Ht on the complete-case block, scored by
# ht_scores. `responses` is already reverse-keyed by the wrapper. Returns list(value, z,
# ncat, abstain): `value` is the full-length score vector (NA per missing-cell respondent; NA
# per complete straightliner, a zero-variance row being structurally unscalable; all NA on
# wholesale abstention -- fewer than two items, complete respondents, or varying respondents).
# `z` is the validated globally zero-based complete block (codes 0..m, NULL on wholesale
# abstention) and `ncat` = m + 1, both handed to the opt-in Monte-Carlo null cutoff;
# `abstain` names the wholesale-abstention cause and count. Ht is a covariance ratio bounded
# in [-1, 1] (the Frechet bound). A fractional cell raises first via assert_integer_responses().
kernel_ht <- function(responses, call = rlang::caller_env()) {
  n <- nrow(responses)
  value <- rep(NA_real_, n)
  assert_integer_responses(responses, "Ht", call)
  # Person scalability needs >= 2 items: with a single column the per-row variance below is
  # NA, turning the guard into `if (NA)` (an untyped error) rather than a clean abstention.
  if (ncol(responses) < 2L) {
    return(list(value = value, z = NULL, ncat = NA_integer_,
                abstain = list(kind = "few_items", n = ncol(responses))))
  }
  complete <- stats::complete.cases(responses)
  if (sum(complete) < 2L) {
    return(list(value = value, z = NULL, ncat = NA_integer_,
                abstain = list(kind = "few_complete", n = sum(complete))))
  }
  z <- responses[complete, , drop = FALSE]
  # rowVars (Bessel n - 1, like stats::var) in one C-vectorised pass instead of a per-row
  # var() dispatch; a row is constant iff its variance is 0.
  row_var <- matrixStats::rowVars(z)
  if (sum(row_var > 0) < 2L) {
    return(list(value = value, z = NULL, ncat = NA_integer_,
                abstain = list(kind = "few_varying", n = sum(row_var > 0))))
  }
  # Global zero-base to codes 0..m. The covariance ratio is translation-invariant, so this
  # does not change the score; it keeps the closed form on small integers, so the
  # min-invariance contract (1..K and 0..(K-1) score identically) holds exactly. z may stay
  # integer: the covariance squares use `^2`, always double.
  z_min <- min(z)
  m <- max(z) - z_min
  z <- z - z_min
  v <- ht_scores(z, m)
  # A straightliner has den == 0 in exact arithmetic but can be floating-point noise, so key
  # its abstention off the row variance; ht_scores already mapped any other non-finite ratio
  # to NA.
  v[row_var == 0] <- NA_real_
  value[complete] <- v
  # A single global ncat = m + 1 is correct even for heterogeneous category counts (an item
  # never reaching the global max has zero frequency at those levels).
  list(value = value, z = z, ncat = m + 1L, abstain = NULL)
}

# Emit the typed warning naming why the Ht Monte-Carlo null could not be resolved: too few
# items (< 2), complete respondents (< 2), non-constant complete respondents (< 2), or a
# "degenerate_null" (a straightliner-dominated sample whose resample leaves < 2 non-constant
# null vectors after perfect-vector exclusion). The class stays
# cier_warning_insufficient_items so cier_screen()'s targeted muffler treats them alike.
warn_ht_abstention <- function(abstain, call) {
  if (identical(abstain$kind, "degenerate_null")) {
    cier_warn(
      "cier_warning_insufficient_items",
      c("Cannot resolve an Ht Monte-Carlo null cutoff: the simulated null \\
         degenerated to fewer than two non-constant response vectors.",
        "i" = "This arises when the scored sample is dominated by straightliners -- \\
               the sum-score-conditional resample then draws mostly constant vectors, \\
               which Ht excludes; returning {.val NA} as the cutoff."),
      data = list(reason = "degenerate Monte-Carlo null"), call = call
    )
    return(invisible(NULL))
  }
  n <- abstain$n
  if (identical(abstain$kind, "few_items")) {
    cier_warn(
      "cier_warning_insufficient_items",
      c("Cannot resolve an Ht Monte-Carlo null cutoff: only {n} item{?s} present.",
        "i" = "Person scalability needs at least two items; returning {.val NA} \\
               as the cutoff."),
      data = list(n_items = n), call = call
    )
  } else if (identical(abstain$kind, "few_complete")) {
    cier_warn(
      "cier_warning_insufficient_items",
      c("Cannot resolve an Ht Monte-Carlo null cutoff: only {n} complete \\
         respondent{?s}.",
        "i" = "Ht needs at least two complete respondents; returning {.val NA} \\
               as the cutoff."),
      data = list(n_used = n), call = call
    )
  } else {
    cier_warn(
      "cier_warning_insufficient_items",
      c("Cannot resolve an Ht Monte-Carlo null cutoff: only {n} complete \\
         respondent{?s} with response variance.",
        "i" = "Ht needs at least two non-constant complete respondents; returning \\
               {.val NA} as the cutoff."),
      data = list(n_used = n), call = call
    )
  }
}

# Resolve the opt-in Ht Monte-Carlo null cutoff from the kernel result `res`. On wholesale
# abstention (no scored block) the cutoff is NA -- flagging no one -- and the warning names
# the cause `res$abstain` carries; otherwise build the null from the scored block at level
# `fpr`. The literal-`cutoff` override is handled in the wrapper. A scored block can still
# yield no cutoff when the null degenerates; surface that as the same typed abstention.
resolve_ht_cutoff <- function(res, fpr, seed, call = rlang::caller_env()) {
  if (is.null(res$z)) {
    warn_ht_abstention(res$abstain, call)
    return(NA_real_)
  }
  # Ht flags low scalability (lower tail), perfect-excluded (constant vectors are
  # unscalable).
  cutoff <- resolve_personfit_null_cutoff(
    res$z, res$ncat, fpr, seed,
    score = ht_scores, tail = "lower", perfect = "excluded"
  )
  if (is.na(cutoff)) {
    warn_ht_abstention(list(kind = "degenerate_null"), call)
  }
  cutoff
}
