# Purpose: The careless plan + content engine for cier_simulate(). sim_build_plan()
#          draws which respondents are careless, with which pattern and extent
#          (none / full / partial / temporary) and span, and returns the
#          per-respondent truth frame. sim_apply_patterns() walks that plan and
#          overwrites each careless row's span with its pattern's block
#          (R/sim-patterns.R), leaving every other cell -- and every attentive or
#          speeder row -- byte-identical to the attentive layer (R/sim-attentive.R).
#          Pure internal kernels: no public cier_simulate() and no S3 object yet.
# Args:    See per-function signatures.
# Returns: sim_build_plan() -> a truth data.frame; sim_apply_patterns() -> the
#          mutated raw integer response matrix.
# Invariants:
#   - Extent is orthogonal to pattern: a careless row draws pattern x extent
#     independently. Content mutation is confined to the careless span [onset,
#     offset] (inclusive); the recovery tail of a temporary row stays attentive.
#   - speeded is NOT a truth column here -- it is a timing fact added by a later
#     layer; this slice's truth is careless / pattern / extent / onset_item /
#     offset_item / params.

# ---- Validators -------------------------------------------------------------

# Check the pattern-weight NAMES: a uniquely-named numeric over the pattern
# allowlist (checkmate folds the named-numeric test into one call).
sim_check_pattern_names <- function(patterns, call) {
  named_ok <- checkmate::test_numeric(patterns, min.len = 1L) &&
    checkmate::test_names(names(patterns), type = "unique")
  if (!named_ok) {
    cier_abort("cier_error_input",
               "{.arg patterns} must be a uniquely-named numeric vector of weights.",
               data = list(arg = "patterns"), call = call)
  }
  bad <- setdiff(names(patterns), sim_pattern_names())
  if (length(bad) > 0L) {
    cier_abort("cier_error_input",
               "{.arg patterns} has unknown pattern name(s): {.val {bad}}.",
               data = list(arg = "patterns", observed = bad), call = call)
  }
  invisible(NULL)
}

# Validate the named pattern-weight vector: uniquely named over the allowlist, each
# weight in [0, 1], summing to 1 (the cross-respondent pattern mixture). Returns it
# unchanged on success.
sim_validate_patterns <- function(patterns, call) {
  sim_check_pattern_names(patterns, call)
  if (anyNA(patterns) || any(patterns < -1e-8) || any(patterns > 1 + 1e-8)) {
    cier_abort("cier_error_input",
               "{.arg patterns} weights must lie in [0, 1].",
               data = list(arg = "patterns"), call = call)
  }
  if (abs(sum(patterns) - 1) > 1e-8) {
    cier_abort("cier_error_input",
               "{.arg patterns} weights must sum to 1.",
               data = list(arg = "patterns"), call = call)
  }
  patterns
}

# Resolve the onset window (lo, hi): NULL defaults to round(c(0.3, 0.8) * p)
# clamped into [1, p] (consistent with Welz & Alfons's 2023 sampled-onset design);
# an explicit window must be two whole numbers with 1 <= lo <= hi <= p.
sim_resolve_onset_window <- function(onset_window, p, call) {
  if (is.null(onset_window)) {
    return(as.integer(pmin(pmax(round(c(0.3, 0.8) * p), 1L), p)))
  }
  ok <- is_finite_whole(onset_window) && length(onset_window) == 2L &&
    all(onset_window >= 1) && all(onset_window <= p) &&
    onset_window[[1L]] <= onset_window[[2L]]
  if (!ok) {
    cier_abort(
      "cier_error_input",
      c("{.arg onset_window} must be two whole numbers 1 <= lo <= hi <= {p}.",
        "x" = "Got {.val {onset_window}}."),
      data = list(arg = "onset_window"), call = call
    )
  }
  as.integer(onset_window)
}

# ---- Extent / span draws ----------------------------------------------------

# Allocate each careless row an extent by EXACT proportional counts (not a sampled
# multinomial), the remainder full, then permute. Counts are clamped so the full
# bucket is never negative when prop_partial + prop_temporary rounds up to the
# whole sample.
sim_draw_extents <- function(n_careless, prop_partial, prop_temporary) {
  if (n_careless == 0L) {
    return(character(0L))
  }
  n_partial <- round(n_careless * prop_partial)
  n_temporary <- min(round(n_careless * prop_temporary), n_careless - n_partial)
  n_full <- n_careless - n_partial - n_temporary
  labels <- rep(c("full", "partial", "temporary"),
                c(n_full, n_partial, n_temporary))
  sample(labels)
}

# Per careless row, the inclusive careless span [onset, offset] given its extent:
# full = [1, p]; partial = [onset in window, p]; temporary = [onset, offset] with
# offset <= p - 1 so at least one recovery column follows.
sim_draw_spans <- function(extents, lo, hi, p) {
  n <- length(extents)
  onset <- integer(n)
  offset <- integer(n)
  is_full <- extents == "full"
  onset[is_full] <- 1L
  offset[is_full] <- p
  is_part <- extents == "partial"
  np <- sum(is_part)
  if (np > 0L) {
    onset[is_part] <- lo + sample.int(hi - lo + 1L, np, replace = TRUE) - 1L
    offset[is_part] <- p
  }
  is_temp <- extents == "temporary"
  nt <- sum(is_temp)
  if (nt > 0L) {
    thi <- min(hi, p - 2L)
    on_t <- lo + sample.int(thi - lo + 1L, nt, replace = TRUE) - 1L
    offset[is_temp] <- as.integer(on_t + floor(stats::runif(nt) * (p - on_t)))
    onset[is_temp] <- on_t
  }
  list(onset = onset, offset = offset)
}

# ---- Truth assembly ---------------------------------------------------------

# Build the n-row truth frame from the per-careless-row draws. Attentive rows carry
# the clean labels (pattern "attentive", extent "none", NA span, empty params).
sim_build_truth <- function(n, rows, pat, ext, onset, offset, pattern_params) {
  careless <- logical(n)
  pattern <- rep("attentive", n)
  extent <- rep("none", n)
  onset_item <- rep(NA_integer_, n)
  offset_item <- rep(NA_integer_, n)
  params <- rep(list(list()), n)
  if (length(rows) > 0L) {
    careless[rows] <- TRUE
    pattern[rows] <- pat
    extent[rows] <- ext
    onset_item[rows] <- onset
    offset_item[rows] <- offset
    for (i in seq_along(rows)) {
      params[[rows[[i]]]] <- pattern_params[[pat[[i]]]] %||% list()
    }
  }
  truth <- data.frame(careless, pattern, extent, onset_item, offset_item,
                      stringsAsFactors = FALSE)
  truth$params <- params
  truth
}

# ---- Plan -------------------------------------------------------------------

# Draw the careless plan for n respondents over p items. `prevalence` sets the
# exact careless count round(prevalence * n); `patterns` the pattern mixture;
# `prop_partial` / `prop_temporary` the extent shares (remainder full);
# `onset_window` the partial / temporary onset sampling range; `pattern_params` the
# per-pattern knobs recorded in the truth's params column. Returns the truth frame.
sim_build_plan <- function(n, p, prevalence, patterns = c(random = 1),
                           prop_partial = 0, prop_temporary = 0,
                           onset_window = NULL, pattern_params = list(),
                           call = rlang::caller_env()) {
  check_count(n, "n", call = call)
  check_count(p, "p", call = call)
  check_number(prevalence, "prevalence", lower = 0, upper = 1, call = call)
  patterns <- sim_validate_patterns(patterns, call)
  check_number(prop_partial, "prop_partial", lower = 0, upper = 1, call = call)
  check_number(prop_temporary, "prop_temporary", lower = 0, upper = 1, call = call)
  if (prop_partial + prop_temporary > 1 + 1e-8) {
    cier_abort("cier_error_input",
               "{.arg prop_partial} + {.arg prop_temporary} must not exceed 1.",
               data = list(arg = "prop_partial"), call = call)
  }
  win <- sim_resolve_onset_window(onset_window, p, call)
  if (prop_temporary > 0 && (p < 3L || win[[1L]] > p - 2L)) {
    cier_abort(
      "cier_error_input",
      c("Temporary carelessness needs at least 3 items and an onset window that \\
         leaves a recovery column.",
        "x" = "Got p = {p}, onset window [{win[[1L]]}, {win[[2L]]}]."),
      data = list(arg = "prop_temporary"), call = call
    )
  }
  n_careless <- round(n * prevalence)
  rows <- if (n_careless > 0L) sort(sample.int(n, n_careless)) else integer(0L)
  pat <- if (n_careless > 0L) {
    names(patterns)[sample.int(length(patterns), n_careless,
                               replace = TRUE, prob = patterns)]
  } else {
    character(0L)
  }
  ext <- sim_draw_extents(n_careless, prop_partial, prop_temporary)
  spans <- sim_draw_spans(ext, win[[1L]], win[[2L]], p)
  sim_build_truth(n, rows, pat, ext, spans$onset, spans$offset, pattern_params)
}

# ---- Content engine ---------------------------------------------------------

# Overwrite each careless non-speeder row's span with its pattern block. Each
# row's knobs come from the truth frame's own `params` column (the single
# provenance source sim_build_plan records), so what the truth records and what
# the engine applies cannot desync. Rows that share a pattern AND a span AND a
# params record are produced together (one producer call). Speeder and
# attentive rows -- and every cell outside a careless span -- are left untouched, so
# the attentive prefix / recovery tail is byte-identical to `attentive`. `items` is
# the validated check_items_simulate list; `truth` the sim_build_plan frame.
sim_apply_patterns <- function(attentive, items, truth,
                               call = rlang::caller_env()) {
  x <- attentive
  storage.mode(x) <- "integer"
  active <- which(truth$careless & truth$pattern != "speeder")
  if (length(active) == 0L) {
    return(x)
  }
  params_col <- truth$params %||% rep(list(list()), nrow(truth))
  # The params fingerprint extends the group key: sim_build_plan records one
  # params list per pattern, but a hand-built truth may vary knobs per row, and
  # grouping by pattern + span alone would silently hand such rows the first
  # row's knobs.
  params_key <- vapply(params_col[active], function(prm) {
    paste(deparse(prm), collapse = " ")
  }, character(1L))
  key <- paste(truth$pattern[active], truth$onset_item[active],
               truth$offset_item[active], params_key, sep = "\r")
  # Group rows that share a pattern AND a span AND params. Level order is
  # first-appearance (factor on `unique(key)`), NOT split()'s locale-sorted
  # default: the groups are visited in that order and each consumes RNG, so a
  # locale-dependent order would make the seeded output differ across machines
  # (the radix-sort lesson from the attentive layer).
  for (grp in split(active, factor(key, levels = unique(key)))) {
    first <- grp[[1L]]
    pat <- truth$pattern[[first]]
    cols <- truth$onset_item[[first]]:truth$offset_item[[first]]
    producer <- sim_block_fun(pat)
    block <- producer(length(grp), items$min[cols], items$max[cols],
                      params_col[[first]] %||% list(), call)
    x[grp, cols] <- block
  }
  x
}
