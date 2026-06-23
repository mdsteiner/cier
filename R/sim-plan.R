# Careless plan + content engine for cier_simulate(), pure internal kernels.
# sim_build_plan() draws which respondents are careless, with which pattern, extent
# (none / full / partial / temporary, drawn independently of pattern), and span,
# returning the per-respondent truth frame. sim_apply_patterns() overwrites each careless
# row's span with its pattern block (sim-patterns.R), leaving all other cells byte-identical
# to the attentive layer. speeded is added by a later timing layer, not here.

# ---- Validators -------------------------------------------------------------

# Check pattern-weight NAMES: uniquely-named numeric over the pattern allowlist.
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

# Validate the pattern-weight mixture: named over the allowlist, weights in [0, 1]
# summing to 1. Returns it unchanged.
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

# Resolve the onset window (lo, hi): NULL defaults to round(c(0.3, 0.8) * p) clamped
# into [1, p] (Welz & Alfons 2023 sampled-onset design); explicit windows must be two
# whole numbers with 1 <= lo <= hi <= p.
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

# Allocate extents by EXACT proportional counts (not a sampled multinomial), remainder
# full, then permute. n_temporary is clamped so the full bucket never goes negative when
# prop_partial + prop_temporary rounds up to the whole sample.
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

# Which careless rows carry fast response times. Speeder rows always speed (their
# carelessness IS the time shift); content patterns speed in an exact prop_speeded share,
# the rest at attentive pace ("slow careless"). prop_speeded >= 1 returns `careless`
# unchanged and draws nothing, keeping the default generator unchanged.
sim_draw_speeded <- function(careless, pattern, prop_speeded) {
  speeded <- careless
  if (prop_speeded >= 1) {
    return(speeded)
  }
  content <- which(careless & pattern != "speeder")
  n_content <- length(content)
  if (n_content == 0L) {
    return(speeded)
  }
  n_speed <- round(n_content * prop_speeded)
  speeded[content] <- sample(rep(c(TRUE, FALSE), c(n_speed, n_content - n_speed)))
  speeded
}

# Inclusive careless span [onset, offset] per row, by extent: full = [1, p];
# partial = [onset in window, p]; temporary = [onset, offset] with offset <= p - 1 so
# at least one recovery column follows.
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

# Build the n-row truth frame from the per-careless-row draws. Attentive rows get the
# clean labels (pattern "attentive", extent "none", NA span, empty params).
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

# Draw the careless plan for n respondents over p items, returning the truth frame.
# prevalence sets the exact careless count round(prevalence * n); patterns the mixture;
# prop_partial / prop_temporary the extent shares (remainder full); onset_window the
# partial/temporary onset range; pattern_params the per-pattern knobs recorded in the
# truth's params column.
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

# Overwrite each careless non-speeder row's span with its pattern block. Knobs come from
# the truth frame's own `params` column (the single provenance source), so truth and engine
# cannot desync. Rows sharing a pattern AND span AND params record are produced together.
# Speeder/attentive rows and every cell outside a careless span are left untouched.
# `items` is the validated check_items_simulate list.
sim_apply_patterns <- function(attentive, items, truth,
                               call = rlang::caller_env()) {
  x <- attentive
  storage.mode(x) <- "integer"
  active <- which(truth$careless & truth$pattern != "speeder")
  if (length(active) == 0L) {
    return(x)
  }
  params_col <- truth$params %||% rep(list(list()), nrow(truth))
  # params fingerprint extends the group key: a hand-built truth may vary knobs per row,
  # and grouping by pattern + span alone would silently hand such rows the first row's knobs.
  params_key <- vapply(params_col[active], function(prm) {
    paste(deparse(prm), collapse = " ")
  }, character(1L))
  key <- paste(truth$pattern[active], truth$onset_item[active],
               truth$offset_item[active], params_key, sep = "\r")
  # Level order is first-appearance (factor on `unique(key)`), NOT split()'s locale-sorted
  # default: each group consumes RNG, so a sorted order would make seeded output differ
  # across machines.
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
