# Purpose: Attentive (careful-respondent) data generation for cier_simulate().
#          A graded response model in the factor-analytic parameterisation draws
#          latent trait vectors, forms per-item latent predictions, and
#          categorises them against per-item thresholds. This is the clean layer
#          the careless mutators (later slices) overwrite. Pure internal kernels:
#          there is no public cier_simulate() and no S3 object yet.
# Args:    See per-function signatures.
# Returns: Integer response matrices (raw, in min..max), latent draws, or
#          metadata lists. Validators raise cier_error_input on bad input.
# Invariants:
#   - One production GRM kernel; reference re-derivations stay in tests/reference.
#   - Reverse-keyed items load NEGATIVELY; the returned matrix is raw (as-clicked),
#     the package-wide contract (indices reverse-score internally from `items`).
#   - Trait draws are standardised to mean 0 / unit variance before factor_cor is
#     imposed, so the latent prediction has unit variance and the categorised
#     marginal matches the target pmf under normal traits.

# ---- Items contract ---------------------------------------------------------

# Validate the per-item `items` frame the simulator needs. Stricter than the
# index validators: `scale` is required (it defines the factor structure;
# min_scales = 1, a single-scale simulation is legitimate) and `max` is required
# on EVERY item (the generator needs each item's full response range, not only the
# reverse-keyed ones). `min` defaults to 1 (the package-wide base) and is read on
# every item. Reuses the shared item-frame / scale / reverse / min / max checks.
# Returns a normalised list(scale, reverse_keyed, min, max, categories), where
# categories_j = max_j - min_j + 1 (heterogeneous ranges allowed).
check_items_simulate <- function(items, n_items, arg = "items",
                                 call = rlang::caller_env()) {
  check_items_frame(items, n_items, arg, call)
  scale <- check_items_scale(items, 1L, arg, call)
  reverse_keyed <- check_items_reverse(items, n_items, arg, call)
  minimum <- check_items_min_personfit(items, n_items, arg, call)
  maximum <- check_items_max_personfit(items, minimum, arg, call)
  # max - min + 1 is whole and finite (both validated), but a huge range
  # overflows the integer type to NA; guard it as a typed error rather than let
  # NA categories crash a downstream seq_len() / threshold check untyped.
  int_max <- .Machine$integer.max
  span <- maximum - minimum + 1
  if (any(span > int_max)) {
    cier_abort(
      "cier_error_input",
      c("An item's response range is too large to enumerate as categories.",
        "x" = "max - min + 1 exceeds {int_max} for some item.",
        "i" = "Likert items have a handful of options; check {.field max} / \\
               {.field min} for a miscoded continuous column."),
      data = list(arg = arg), call = call
    )
  }
  list(scale = scale, reverse_keyed = reverse_keyed,
       min = minimum, max = maximum,
       categories = as.integer(span))
}

# ---- Default loadings -------------------------------------------------------

# Default loading matrix: one factor per unique scale label, each item loading
# its own scale factor at 0.7 with zero cross-loadings. `items` is the validated
# check_items_simulate list. Scales are ordered by radix (byte) sort, NOT the
# default locale collation: the column order fixes the trait-draw -> scale
# mapping and the meaning of a user-supplied factor_cor, so it must be identical
# across locales for a seed to reproduce bytewise on every machine.
sim_default_loadings <- function(items) {
  scales <- sort(unique(items$scale), method = "radix")
  m <- length(scales)
  p <- length(items$scale)
  out <- matrix(0, nrow = p, ncol = m)
  for (k in seq_len(m)) {
    out[items$scale == scales[k], k] <- 0.7
  }
  colnames(out) <- scales
  out
}

# ---- Validators -------------------------------------------------------------

sim_validate_loadings <- function(loadings, p, call) {
  if (!is.matrix(loadings) || !is.numeric(loadings)) {
    cier_abort("cier_error_input", "{.arg loadings} must be a numeric matrix.",
               data = list(arg = "loadings"), call = call)
  }
  if (nrow(loadings) != p) {
    cier_abort(
      "cier_error_input",
      c("{.arg loadings} must have {p} rows (one per item).",
        "x" = "Got {nrow(loadings)}."),
      data = list(arg = "loadings", observed = nrow(loadings), expected = p),
      call = call
    )
  }
  if (ncol(loadings) < 1L) {
    cier_abort("cier_error_input",
               "{.arg loadings} must have at least one column.",
               data = list(arg = "loadings"), call = call)
  }
  if (any(!is.finite(loadings))) {
    cier_abort("cier_error_input", "{.arg loadings} must be finite.",
               data = list(arg = "loadings"), call = call)
  }
  invisible(NULL)
}

sim_validate_factor_cor <- function(factor_cor, m, call) {
  if (!is.matrix(factor_cor) || !is.numeric(factor_cor)) {
    cier_abort("cier_error_input", "{.arg factor_cor} must be a numeric matrix.",
               data = list(arg = "factor_cor"), call = call)
  }
  if (!identical(dim(factor_cor), c(m, m))) {
    cier_abort(
      "cier_error_input",
      c("{.arg factor_cor} must be square of dimension {m} x {m}.",
        "x" = "Got {nrow(factor_cor)} x {ncol(factor_cor)}."),
      data = list(arg = "factor_cor"), call = call
    )
  }
  if (any(!is.finite(factor_cor))) {
    cier_abort("cier_error_input", "{.arg factor_cor} must be finite.",
               data = list(arg = "factor_cor"), call = call)
  }
  # check.attributes = FALSE: a valid factor_cor may carry dimnames (scale
  # labels), and diag() of a named matrix returns a NAMED vector, so attribute
  # comparison would falsely reject a perfectly good correlation matrix.
  if (!isTRUE(all.equal(factor_cor, t(factor_cor), tolerance = 1e-8,
                        check.attributes = FALSE))) {
    cier_abort("cier_error_input", "{.arg factor_cor} must be symmetric.",
               data = list(arg = "factor_cor"), call = call)
  }
  if (!isTRUE(all.equal(diag(factor_cor), rep(1, m), tolerance = 1e-8,
                        check.attributes = FALSE))) {
    cier_abort("cier_error_input",
               "{.arg factor_cor} must have unit diagonal (correlation matrix).",
               data = list(arg = "factor_cor"), call = call)
  }
  eigvals <- eigen(factor_cor, symmetric = TRUE, only.values = TRUE)$values
  if (min(eigvals) <= 1e-8) {
    cier_abort("cier_error_input",
               "{.arg factor_cor} must be positive definite.",
               data = list(arg = "factor_cor"), call = call)
  }
  invisible(NULL)
}

sim_validate_thresholds_shape <- function(thresholds, p, call) {
  if (!is.list(thresholds) || length(thresholds) != p) {
    cier_abort(
      "cier_error_input",
      c("{.arg thresholds} must be a list of length {p}.",
        "x" = if (is.list(thresholds)) {
          paste0("Got length ", length(thresholds), ".")
        } else {
          "Got a non-list value."
        }),
      data = list(arg = "thresholds"), call = call
    )
  }
  invisible(NULL)
}

sim_validate_threshold_entry <- function(bj, kj, j, call) {
  if (is.na(kj)) {
    cier_abort(
      "cier_error_input",
      "{.arg thresholds[[{j}]]} cannot be validated: item {j} has no category \\
       count.",
      data = list(arg = "thresholds", item = j), call = call
    )
  }
  shape_ok <- is.numeric(bj) && length(bj) == kj - 1L &&
    !anyNA(bj) && all(is.finite(bj))
  if (!shape_ok) {
    cier_abort(
      "cier_error_input",
      c("{.arg thresholds[[{j}]]} must be a finite numeric vector of \\
         length {kj - 1L}.",
        "x" = "Got length {length(bj)}."),
      data = list(arg = "thresholds", item = j), call = call
    )
  }
  if (length(bj) > 1L && any(diff(bj) <= 0)) {
    cier_abort("cier_error_input",
               "{.arg thresholds[[{j}]]} must be strictly increasing.",
               data = list(arg = "thresholds", item = j), call = call)
  }
  invisible(NULL)
}

# Validate a user-supplied thresholds list against the validated `items`
# (one strictly-increasing cut vector of length categories_j - 1 per item).
sim_validate_thresholds <- function(thresholds, items, call) {
  sim_validate_thresholds_shape(thresholds, length(items$categories), call)
  for (j in seq_along(thresholds)) {
    sim_validate_threshold_entry(thresholds[[j]], items$categories[[j]], j, call)
  }
  invisible(NULL)
}

# ---- Trait draws ------------------------------------------------------------

sim_trait_param_scalar <- function(params, name, default, call,
                                   lower = -Inf, upper = Inf) {
  value <- params[[name]] %||% default
  check_number(value, paste0("trait_params$", name),
               lower = lower, upper = upper, call = call)
  value
}

sim_trait_param_vector <- function(params, name, default, m, call) {
  value <- params[[name]] %||% default
  if (!is.numeric(value) || anyNA(value) || any(!is.finite(value))) {
    cier_abort("cier_error_input",
               "{.arg trait_params${name}} must be a finite numeric vector.",
               data = list(arg = paste0("trait_params$", name)), call = call)
  }
  if (!length(value) %in% c(1L, m)) {
    cier_abort(
      "cier_error_input",
      "{.arg trait_params${name}} must have length 1 or {m}.",
      data = list(arg = paste0("trait_params$", name)), call = call
    )
  }
  rep_len(value, m)
}

# Impose the target factor covariance on standardised marginal draws.
sim_trait_correlate <- function(z, factor_cor) {
  z %*% chol(factor_cor)
}

sim_draw_normal_traits <- function(n, m) {
  matrix(stats::rnorm(n * m), nrow = n, ncol = m)
}

# Multivariate-t traits standardised to unit marginal variance: a shared chi-square
# scaling per respondent gives heavier-than-normal tails (df default 5). The floor
# is 3, not just above 2: a t variance exists for df > 2, but just above 2 the
# 1/chi-square scale is so heavy-tailed that the realised sample variance collapses
# far below 1 in any finite sample, silently breaking the unit-variance contract.
sim_draw_t_traits <- function(n, m, params, call) {
  df <- sim_trait_param_scalar(params, "df", 5, call, lower = 3)
  z <- sim_draw_normal_traits(n, m)
  scale <- sqrt((df - 2) / stats::rchisq(n, df = df))
  z * matrix(scale, nrow = n, ncol = m)
}

# Skew-normal traits standardised to mean 0 / unit variance (alpha default 5;
# positive alpha gives positive skew).
sim_draw_skew_normal_traits <- function(n, m, params, call) {
  alpha <- sim_trait_param_vector(params, "alpha", 5, m, call)
  delta <- alpha / sqrt(1 + alpha^2)
  u <- sim_draw_normal_traits(n, m)
  z <- sim_draw_normal_traits(n, m)
  raw <- sweep(abs(u), 2L, delta, "*") + sweep(z, 2L, sqrt(1 - delta^2), "*")
  location <- delta * sqrt(2 / pi)
  scale <- sqrt(1 - 2 * delta^2 / pi)
  sweep(sweep(raw, 2L, location, "-"), 2L, scale, "/")
}

# Bimodal traits: a two-component normal mixture, component means +/- sep
# (default 1), component SD 0.5, mixing weights (default equal) -- standardised to
# mean 0 / unit variance so `sep` / weights are free and each factor stays
# unit-variance. Larger sep leaves a deeper central trough (two genuine camps).
sim_draw_bimodal_traits <- function(n, m, params, call) {
  sep <- sim_trait_param_scalar(params, "sep", 1, call, lower = 0)
  weights <- params[["weights"]] %||% c(0.5, 0.5)
  ok <- is.numeric(weights) && length(weights) == 2L &&
    all(is.finite(weights)) && all(weights > 0) && abs(sum(weights) - 1) <= 1e-8
  if (!ok) {
    cier_abort(
      "cier_error_input",
      "{.arg trait_params$weights} must be two positive weights summing to 1.",
      data = list(arg = "trait_params$weights"), call = call
    )
  }
  comp_sd <- 0.5
  means <- c(-sep, sep)
  mu <- sum(weights * means)
  v <- sum(weights * (means^2 + comp_sd^2)) - mu^2
  comp <- sample.int(2L, n * m, replace = TRUE, prob = weights)
  raw <- stats::rnorm(n * m, mean = means[comp], sd = comp_sd)
  matrix((raw - mu) / sqrt(v), nrow = n, ncol = m)
}

# The trait_params keys each distribution reads (normal takes none). Used to
# reject misspelled / misplaced parameters so a typo cannot silently fall back to
# a default the user did not intend.
sim_trait_param_keys <- function(distribution) {
  switch(distribution,
    normal      = character(0),
    t           = "df",
    skew_normal = "alpha",
    bimodal     = c("sep", "weights")
  )
}

sim_check_trait_params <- function(distribution, params, call) {
  reject_unknown_keys(names(params), sim_trait_param_keys(distribution),
                      "trait_params", "trait parameter", call)
  invisible(NULL)
}

# Draw n latent trait vectors with target covariance factor_cor. Marginal draws
# are standardised before the correlation is imposed.
sim_draw_traits <- function(n, factor_cor, distribution = "normal",
                            params = list(), call = rlang::caller_env()) {
  if (!is.list(params)) {
    cier_abort("cier_error_input", "{.arg trait_params} must be a list.",
               data = list(arg = "trait_params"), call = call)
  }
  check_choice(distribution, "trait_distribution",
               sim_trait_distributions(), call = call)
  sim_check_trait_params(distribution, params, call)
  m <- nrow(factor_cor)
  z <- switch(distribution,
    normal      = sim_draw_normal_traits(n, m),
    t           = sim_draw_t_traits(n, m, params, call),
    skew_normal = sim_draw_skew_normal_traits(n, m, params, call),
    bimodal     = sim_draw_bimodal_traits(n, m, params, call),
    cier_abort("cier_error_input",
               "Unsupported {.arg trait_distribution}: {.val {distribution}}.",
               data = list(arg = "trait_distribution"), call = call)
  )
  sim_trait_correlate(z, factor_cor)
}

# ---- GRM thresholding -------------------------------------------------------

# Map continuous latent predictions eta to ordered categories 1..K using per-item
# cut-points: category = 1 + #{thresholds at or below eta}. findInterval()'s
# left-closed default puts a latent value exactly on a threshold in the higher
# category. The column loop is unavoidable since thresholds differ per item.
sim_grm_categorise <- function(eta, thresholds) {
  p <- ncol(eta)
  out <- matrix(0L, nrow = nrow(eta), ncol = p)
  for (j in seq_len(p)) {
    out[, j] <- findInterval(eta[, j], thresholds[[j]]) + 1L
  }
  out
}

# ---- Latent structure -------------------------------------------------------

# Negate the loadings of reverse-keyed items so a reverse item correlates
# negatively with its scale-mates (the raw-orientation contract). `items` is the
# validated check_items_simulate list.
sim_effective_loadings <- function(loadings, items) {
  eff <- loadings
  rk <- items$reverse_keyed
  if (any(rk)) {
    eff[rk, ] <- -eff[rk, , drop = FALSE]
  }
  eff
}

# Per-item residual variance that makes the latent prediction unit-variance:
# psi_j = max(0, 1 - factor_var_j), where factor_var_j = diag(L Phi L^T). The
# pmax() floor keeps psi non-negative if a heavy / correlated loading set pushes
# the factor variance above 1.
sim_residual_variance <- function(effective_loadings, factor_cor) {
  factor_var <- diag(effective_loadings %*% factor_cor %*% t(effective_loadings))
  pmax(0, 1 - factor_var)
}

# Reject a loadings / factor_cor pair whose implied per-item factor variance
# (communality) exceeds 1: the residual variance would then be clamped to 0 and
# the latent prediction would have variance above 1, silently breaking the
# closed-form marginals contract. The communality diag(L Phi L^T) is invariant to
# reverse-key sign flips, so the raw loadings are checked. A communality of
# exactly 1 (zero residual) is allowed -- the latent is still unit-variance.
sim_validate_communality <- function(loadings, factor_cor, call) {
  communality <- diag(loadings %*% factor_cor %*% t(loadings))
  if (any(communality > 1 + 1e-8)) {
    cier_abort(
      "cier_error_input",
      c("{.arg loadings} imply an explained variance above 1 for some item.",
        "x" = "Largest communality diag(L Phi L') = {round(max(communality), 3)} \\
               (must be <= 1).",
        "i" = "Lower the loadings or the factor correlations so each item's \\
               factor variance does not exceed 1."),
      data = list(arg = "loadings", observed = max(communality)), call = call
    )
  }
  invisible(NULL)
}

# ---- Assembly ---------------------------------------------------------------

# Validate n / items / distribution and resolve loadings, factor_cor, thresholds
# (defaulting and the marginals-vs-thresholds dispatch). Returns the validated
# items list plus the resolved generator parameters.
sim_attentive_inputs <- function(n, items, loadings, factor_cor, thresholds,
                                 marginals, trait_distribution, call,
                                 validated_items = NULL) {
  check_count(n, "n", call = call)
  # Accept a pre-validated check_items_simulate() list to skip a second walk of
  # the items contract: the orchestrator validates once and threads it through.
  it <- validated_items %||%
    check_items_simulate(items, NROW(items), arg = "items", call = call)
  check_choice(trait_distribution, "trait_distribution",
               sim_trait_distributions(), call = call)
  p <- length(it$scale)
  if (is.null(loadings)) {
    loadings <- sim_default_loadings(it)
  }
  sim_validate_loadings(loadings, p, call)
  m <- ncol(loadings)
  if (is.null(factor_cor)) {
    factor_cor <- diag(m)
  }
  sim_validate_factor_cor(factor_cor, m, call)
  sim_validate_communality(loadings, factor_cor, call)
  thresholds <- sim_resolve_thresholds(it, marginals, thresholds, call)
  list(items = it, loadings = loadings, factor_cor = factor_cor,
       thresholds = thresholds)
}

# Generate the latent predictions eta (n x p) and carry the generator metadata.
# eta_j = traits %*% effective_loadings_j + N(0, psi_j); under normal traits and
# the unit-variance residual it is marginally N(0, 1).
sim_attentive_latent <- function(n, items, loadings = NULL, factor_cor = NULL,
                                 thresholds = NULL, marginals = NULL,
                                 trait_distribution = "normal",
                                 trait_params = list(), validated_items = NULL,
                                 call = rlang::caller_env()) {
  inputs <- sim_attentive_inputs(n, items, loadings, factor_cor, thresholds,
                                 marginals, trait_distribution, call,
                                 validated_items = validated_items)
  it <- inputs$items
  eff <- sim_effective_loadings(inputs$loadings, it)
  traits <- sim_draw_traits(n, inputs$factor_cor,
                            distribution = trait_distribution,
                            params = trait_params, call = call)
  psi <- sim_residual_variance(eff, inputs$factor_cor)
  p <- length(it$scale)
  noise <- matrix(stats::rnorm(n * p) * rep(sqrt(psi), each = n),
                  nrow = n, ncol = p)
  eta <- traits %*% t(eff) + noise
  list(eta = eta, traits = traits, items = it,
       loadings = inputs$loadings, effective_loadings = eff,
       factor_cor = inputs$factor_cor, thresholds = inputs$thresholds,
       residual_variance = psi, trait_distribution = trait_distribution,
       trait_params = trait_params)
}

# Offset categories (1..K_j) to the item's raw coding (min_j..max_j).
sim_offset_to_range <- function(categories, mins) {
  sweep(categories, 2L, mins - 1L, "+")
}

# Generator metadata carried alongside the responses (loadings, thresholds, RT
# parameters are added in later slices).
sim_attentive_metadata <- function(latent) {
  list(loadings = latent$loadings,
       effective_loadings = latent$effective_loadings,
       factor_cor = latent$factor_cor,
       thresholds = latent$thresholds,
       residual_variance = latent$residual_variance,
       trait_distribution = latent$trait_distribution,
       trait_params = latent$trait_params)
}

sim_attentive_with_metadata <- function(n, items, loadings = NULL,
                                        factor_cor = NULL, thresholds = NULL,
                                        marginals = NULL,
                                        trait_distribution = "normal",
                                        trait_params = list(),
                                        validated_items = NULL,
                                        call = rlang::caller_env()) {
  latent <- sim_attentive_latent(n, items, loadings, factor_cor, thresholds,
                                 marginals, trait_distribution, trait_params,
                                 validated_items = validated_items, call = call)
  cats <- sim_grm_categorise(latent$eta, latent$thresholds)
  responses <- sim_offset_to_range(cats, latent$items$min)
  storage.mode(responses) <- "integer"
  list(responses = responses,
       metadata = sim_attentive_metadata(latent), latent = latent)
}

# Generate attentive responses from a graded response model parameterised by a
# factor-analytic structure. Returns a raw integer matrix n x p with values in
# min_j..max_j per item. `call` is captured at this boundary and threaded down so
# validation errors are attributed to the caller, not an internal helper.
sim_attentive <- function(n, items, loadings = NULL, factor_cor = NULL,
                          thresholds = NULL, marginals = NULL,
                          trait_distribution = "normal", trait_params = list()) {
  sim_attentive_with_metadata(n, items, loadings, factor_cor, thresholds,
                              marginals, trait_distribution, trait_params,
                              call = rlang::caller_env())$responses
}
