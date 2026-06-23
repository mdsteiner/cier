# Public cier_simulate() orchestrator: wires the kernels (sim-attentive.R, sim-plan.R,
# sim-patterns.R, sim-times.R, sim-direct.R) into a validated cier_sim. No statistical
# logic here. The order attentive -> plan -> content -> times -> checks is load-bearing:
# reordering reorders every seeded draw.

# Equal weight on every pattern, per Welz & Alfons (2023). No study estimates a real-world
# mixture, so unequal weighting would be invented.
sim_default_patterns <- function() {
  stats::setNames(rep(1 / 8, 8L), sim_pattern_names())
}

# Validate the per-pattern knob list: a named list whose names subset the active patterns,
# each entry a list. Catches typos that would otherwise plant nothing silently.
sim_check_pattern_params <- function(pattern_params, patterns, call) {
  if (!is.list(pattern_params)) {
    cier_abort("cier_error_input",
               "{.arg pattern_params} must be a list.",
               data = list(arg = "pattern_params"), call = call)
  }
  if (length(pattern_params) == 0L) {
    return(invisible(NULL))
  }
  nms <- names(pattern_params)
  check_unique_names(nms, "pattern_params", call)
  reject_unknown_keys(nms, names(patterns), "pattern_params",
                      "pattern_params name", call)
  knob_names <- sim_pattern_knob_names()
  for (nm in nms) {
    knobs <- pattern_params[[nm]]
    if (!is.list(knobs)) {
      cier_abort(
        "cier_error_input",
        "{.arg pattern_params${nm}} must be a named list of knobs.",
        data = list(arg = "pattern_params", observed = nm), call = call
      )
    }
    if (length(knobs) > 0L) {
      check_unique_names(names(knobs), paste0("pattern_params$", nm), call)
      reject_unknown_keys(names(knobs), knob_names[[nm]],
                          paste0("pattern_params$", nm), "knob", call)
    }
  }
  invisible(NULL)
}

# Validate the numeric scalar knobs BEFORE the first RNG draw so a bad value aborts before
# any draws. Pure (no RNG), so hoisting leaves the stream unchanged. Redundant with the
# inner kernels' own copies (kept for direct callers). n_checks allows the default 0.
sim_validate_scalars <- function(prop_partial, prop_temporary, onset_window, p,
                                 n_checks, p_fail_careless, p_fail_attentive,
                                 call) {
  check_number(prop_partial, "prop_partial", lower = 0, upper = 1, call = call)
  check_number(prop_temporary, "prop_temporary", lower = 0, upper = 1,
               call = call)
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
  check_count_nonneg(n_checks, "n_checks", call = call)
  check_number(p_fail_careless, "p_fail_careless", lower = 0, upper = 1,
               call = call)
  check_number(p_fail_attentive, "p_fail_attentive", lower = 0, upper = 1,
               call = call)
  invisible(NULL)
}

#' Simulate survey responses with planted careless responding
#'
#' Generates a survey dataset -- raw responses, response times, and optional
#' attention checks -- in which a known share of respondents answer carelessly
#' with a known pattern, extent, and onset. Attentive responses come from a
#' graded response model over the scales in `items`; careless rows are overwritten
#' on their planted span and answer faster there. The returned object feeds the
#' shipped indices directly: `cier_screen(x$responses, x$items)`,
#' [cier_total_time()]`(x$seconds)`,
#' [cier_page_time()]`(x$page_seconds, x$items_per_page)`, and
#' [cier_attention()]`(x$checks, x$pass)`, each scored against `x$truth`.
#'
#' @details
#' **Scope.** For power analysis, study design, method comparison on planted patterns,
#' recovery / property tests, and replicating published designs. Not evidence of real-world
#' validity: simulated carelessness is more stylized than real, and a simulation-trained
#' advantage may not transfer (Schroeders et al. 2022). Anchor accuracy claims to labeled
#' real data.
#'
#' **Attentive model.** A graded response model: one factor per scale in `items`, default
#' loadings 0.7, identity factor correlations, latent traits from `trait_distribution`.
#' Reverse-keyed items load negatively, so the matrix is raw (as-clicked) and keying-sensitive
#' indices reverse-score internally from `items`. Thresholds come from `marginals` (presets or
#' probability vectors; default `"peaked"`, an as-clicked target) or raw `thresholds`
#' (mutually exclusive). Default loadings can leave synonym / antonym pairs below
#' [cier_psychsyn()] / [cier_psychant()]'s `critical_r`; raise `loadings` toward 0.8 or lower
#' `critical_r`.
#'
#' **Patterns.** `patterns` weights the eight styles: `random`, `straightline`, `midpoint`,
#' `extreme`, `diagonal` (with a bounce option), `alternating`, `markov`, and `speeder`
#' (attentive content but careless times); the default weights all equally. Per-style knobs
#' go in `pattern_params`, a named list of named knob lists (e.g.
#' `list(straightline = list(anchor = "value", value = 5))`). An unknown knob name is a typed
#' error; what is applied is recorded in `truth$params`. Accepted knobs:
#'
#' - `random`, `speeder`: none.
#' - `straightline`: `anchor` (`"position"` or `"value"`, default `"position"`),
#'   `q` (relative position 0..1), `value` (a whole number in every item's
#'   range), `switch_prob` (0..1, default 0).
#' - `midpoint`: `jitter` (non-negative whole number, default 0).
#' - `extreme`: `p_high` (0..1, default 0.5).
#' - `diagonal`: `step` (positive whole number, default 1), `start` (1..K_min),
#'   `bounce` (logical, default `FALSE`).
#' - `alternating`: `period` (whole `>= 2`, default 2), `values` (positions in
#'   1..K_min, **not** raw response values), `start_offset` (0..`period` - 1).
#' - `markov`: `transition` (a K_min x K_min row-stochastic matrix), `initial`
#'   (K_min probabilities).
#'
#' **Extent and onset.** Each careless row draws an extent: `full` (whole row), `partial`
#' (sampled onset to the end), or `temporary` (a bounded window, then recovered). The
#' inclusive careless span `[onset_item, offset_item]` is `[1, p]` for full; partial samples
#' its onset from `onset_window` (default `round(c(0.3, 0.8) * p)`) and runs to `p`; temporary
#' additionally samples an offset of at most `p - 1`, leaving a recovery tail of at least one
#' item (hence needs at least three items). Both content overwrite and time shift are confined
#' to the span.
#'
#' **Response times.** Lognormal per-item times with a per-respondent pace intercept
#' (`respondent_sd`), so attentive and careless distributions overlap. Careless rows shift to
#' the careless log-mean inside their span only. Defaults (attentive median 8 s/item; careless
#' median 1.5 s/item, below the 2 s/item page-screening floor of Bowling et al. 2023;
#' `sigma = 0.5`; `respondent_sd = 1.2`) are conventions, not estimates, overridable via
#' `timing`. The `speeded` truth column records rows carrying careless times: every `speeder`
#' row plus the `prop_speeded` share of content-careless rows (default 1). Set
#' `prop_speeded < 1` to plant "slow careless" rows (careless content at attentive pace).
#'
#' **Attention checks.** With `n_checks > 0`, separate five-point instructed-response checks
#' are generated (never written into `$responses`): careless rows fail each with
#' `p_fail_careless`, attentive rows with `p_fail_attentive` (defaults 0.75 / 0.05).
#' Failure is per-check independent, so one check overlaps the groups but the failed-check
#' count separates them sharply with several checks (AUC near 1.0 at 8). Lower
#' `p_fail_careless` for a realistic power analysis. Speeders are careless rows, so they also
#' fail at `p_fail_careless`.
#'
#' **Reproducibility.** With `seed`, the generator runs on a local RNG state (the session
#' stream is saved and restored); the same seed reproduces the object bytewise.
#'
#' @param n Positive whole number of respondents.
#' @param items Data frame describing the battery, one row per item: `scale`
#'   (character, defines the factor structure), `max` (largest response option,
#'   required), and optionally `min` (default 1) and `reverse_keyed` (logical).
#'   Heterogeneous ranges allowed. Returned untouched as `$items`.
#' @param prevalence Number in `[0, 1]` (default 0.2). Share made careless;
#'   exactly `round(n * prevalence)` rows (round half to even).
#' @param patterns Named numeric vector of pattern weights summing to 1, over
#'   `random`, `straightline`, `midpoint`, `extreme`, `diagonal`,
#'   `alternating`, `markov`, `speeder`. `NULL` (default) weights all equally.
#' @param pattern_params Named list of per-pattern knob lists, recorded in
#'   `truth$params`. Names must appear in `patterns`.
#' @param prop_partial,prop_temporary Numbers in `[0, 1]` summing to at most 1 (both default
#'   0). Shares of careless rows with `partial` / `temporary` extent, as exact counts like
#'   `prevalence`; `prop_partial` applied first, `prop_temporary` takes the remainder, the
#'   rest are `full`.
#' @param prop_speeded Number in `[0, 1]` (default 1). Share of content-careless rows (every
#'   style except `speeder`) that also answer at careless speed; the rest are "slow careless"
#'   (careless content, attentive pace). `speeder` rows always speed and attentive rows never,
#'   so default 1 means every careless row speeds (an exact count, no extra RNG).
#' @param onset_window Optional two whole numbers `c(lo, hi)` inside `[1, p]`: the sampling
#'   range for partial / temporary onsets. `NULL` (default) uses `round(c(0.3, 0.8) * p)`.
#' @param loadings,factor_cor,thresholds Optional GRM overrides: a `p x m` loading matrix
#'   (supplied **keyed-positive** -- magnitudes; reverse-keyed rows are negated internally, so
#'   all-negative reverse rows are a typed error, not double-negated forward data), an `m x m`
#'   factor correlation matrix, and a per-item list of strictly increasing threshold vectors.
#'   `thresholds` is mutually exclusive with `marginals`.
#' @param marginals Target observed distribution(s): a preset (`"peaked"`, `"uniform"`,
#'   `"skewed_right"`, `"skewed_left"`, `"bimodal"`), an explicit probability vector, or a
#'   per-scale / per-item spec. Default `"peaked"`. Met exactly (large-sample limit) only under
#'   `trait_distribution = "normal"`; otherwise the realised marginal drifts -- supply raw
#'   `thresholds` to pin it (mutually exclusive with `marginals`).
#' @param trait_distribution One of `"normal"`, `"skew_normal"`, `"t"`,
#'   `"bimodal"` for the latent traits.
#' @param trait_params Named list of distribution parameters (e.g.
#'   `list(df = 5)` for `"t"`).
#' @param timing Named list overriding response-time parameters `mu_att`, `mu_car` (log-mean
#'   seconds per item, attentive / careless), `sigma` (cell-level log-SD), `respondent_sd`
#'   (pace-intercept SD), `min_seconds` (strictly positive floor). Unknown names error.
#' @param items_per_page Page structure for `$page_seconds`: `NULL` (default) puts each
#'   contiguous run of `items$scale` on its own page; a single whole number chunks uniform
#'   pages of that size (last page takes the remainder); a vector of positive counts summing to
#'   the item count is used verbatim.
#' @param n_checks Non-negative whole number of attention checks (default 0:
#'   `$checks` and `$pass` are `NULL`).
#' @param p_fail_careless,p_fail_attentive Numbers in `[0, 1]`: per-check
#'   failure probabilities for careless and attentive rows.
#' @param seed Optional single whole number. Seeds the generator locally; the
#'   result reproduces bytewise.
#'
#' @return A `cier_sim` object: a list with `$responses` (n x p raw integer
#'   matrix), `$items` (untouched), `$seconds` (total time per respondent),
#'   `$page_seconds` with `$items_per_page`, `$checks` and `$pass` (or `NULL`),
#'   and `$truth` -- one row per respondent with `careless`, `pattern`
#'   (`"attentive"` on clean rows), `extent` (`"none"` / `"full"` / `"partial"` /
#'   `"temporary"`), `onset_item` and `offset_item` (inclusive careless span; `NA`
#'   on attentive rows), `speeded` (carries careless times -- see `prop_speeded`),
#'   and `params` (applied pattern knobs, a list-column: drop or serialize before
#'   `write.csv()`). Generator metadata (loadings, thresholds, residual variances,
#'   resolved timing, failure probabilities, seed) is attached as the
#'   `"generator"` attribute.
#'
#' @references
#' Bowling, N. A., Huang, J. L., Brower, C. K., & Bragg, C. B. (2023). The
#' quick and the careless: The construct validity of page time as a measure of
#' insufficient effort responding to surveys. *Organizational Research
#' Methods*, 26(2), 323-352. \doi{10.1177/10944281211056520}
#'
#' Bruhlmann, F., Petralito, S., Aeschbach, L. F., & Opwis, K. (2020). The
#' quality of data collected online: An investigation of careless responding
#' in a crowdsourced sample. *Methods in Psychology*, 2, 100022.
#' \doi{10.1016/j.metip.2020.100022}
#'
#' Curran, P. G. (2016). Methods for the detection of carelessly invalid
#' responses in survey data. *Journal of Experimental Social Psychology*, 66,
#' 4-19. \doi{10.1016/j.jesp.2015.07.006}
#'
#' Goldammer, P., Stöckli, P. L., Escher, Y. A., Annen, H., Jonas, K., &
#' Antonakis, J. (2024). Careless responding detection revisited: Accuracy of
#' direct and indirect measures. *Behavior Research Methods*, 56(8), 8422-8449.
#' \doi{10.3758/s13428-024-02484-3}
#'
#' Huang, J. L., Curran, P. G., Keeney, J., Poposki, E. M., & DeShon, R. P.
#' (2012). Detecting and deterring insufficient effort responding to surveys.
#' *Journal of Business and Psychology*, 27(1), 99-114.
#'
#' Schroeders, U., Schmidt, C., & Gnambs, T. (2022). Detecting careless
#' responding in survey data using stochastic gradient boosting.
#' *Educational and Psychological Measurement*, 82(1), 29-56.
#' \doi{10.1177/00131644211004708}
#'
#' Voss, N. M. (2024). The effects of careless responding on the fit of
#' confirmatory factor analysis and item response theory models. *Behavior
#' Research Methods*, 56(2), 577-599. \doi{10.3758/s13428-023-02074-9}
#'
#' Welz, M., & Alfons, A. (2023). When respondents don't care anymore:
#' Identifying the onset of careless responding. *arXiv:2303.07167*.
#' \doi{10.48550/arXiv.2303.07167}
#'
#' @seealso [cier_screen()] to score the responses; [cier_total_time()],
#'   [cier_page_time()], [cier_attention()] for the timing and direct slots.
#' @export
#' @examples
#' items <- data.frame(scale = rep(c("E", "A", "C"), each = 4L), max = 5)
#' sim <- cier_simulate(n = 120, items, prevalence = 0.2, seed = 2026)
#' sim
#'
#' # Score planted rows with a content index and the timing family.
#' table(flagged = cier_longstring(sim$responses)$flagged,
#'       truth = sim$truth$careless)
#' cier_total_time(sim$seconds)
cier_simulate <- function(n, items, prevalence = 0.2, patterns = NULL,
                          pattern_params = list(), prop_partial = 0,
                          prop_temporary = 0, prop_speeded = 1,
                          onset_window = NULL,
                          loadings = NULL, factor_cor = NULL,
                          thresholds = NULL, marginals = NULL,
                          trait_distribution = "normal",
                          trait_params = list(), timing = list(),
                          items_per_page = NULL, n_checks = 0L,
                          p_fail_careless = 0.75, p_fail_attentive = 0.05,
                          seed = NULL) {
  call <- rlang::caller_env()
  check_count(n, "n", call = call)
  it <- check_items_simulate(items, NROW(items), arg = "items", call = call)
  p <- length(it$scale)
  check_number(prevalence, "prevalence", lower = 0, upper = 1, call = call)
  check_number(prop_speeded, "prop_speeded", lower = 0, upper = 1, call = call)
  if (is.null(patterns)) {
    patterns <- sim_default_patterns()
  }
  patterns <- sim_validate_patterns(patterns, call)
  sim_check_pattern_params(pattern_params, patterns, call)
  sim_validate_scalars(prop_partial, prop_temporary, onset_window, p, n_checks,
                       p_fail_careless, p_fail_attentive, call)
  resolved_timing <- sim_resolve_timing(timing, call)
  pages <- sim_resolve_pages(items_per_page, it$scale, call)
  if (!is.null(seed)) {
    check_int(seed, "seed", call = call)
  }
  run <- function() {
    attentive <- sim_attentive_with_metadata(
      n, items, loadings = loadings, factor_cor = factor_cor,
      thresholds = thresholds, marginals = marginals,
      trait_distribution = trait_distribution, trait_params = trait_params,
      validated_items = it, call = call
    )
    truth <- sim_build_plan(n, p, prevalence, patterns,
                            prop_partial = prop_partial,
                            prop_temporary = prop_temporary,
                            onset_window = onset_window,
                            pattern_params = pattern_params, call = call)
    responses <- sim_apply_patterns(attentive$responses, it, truth,
                                    call = call)
    # speeded = every speeder row plus the prop_speeded share of content-careless rows.
    # At prop_speeded = 1 this draws nothing, so the default stream is unchanged.
    truth$speeded <- sim_draw_speeded(truth$careless, truth$pattern, prop_speeded)
    cells <- sim_times_lognormal(n, p, truth$speeded, truth$onset_item,
                                 truth$offset_item, timing = resolved_timing,
                                 call = call)
    direct <- sim_direct_checks(n, truth$careless, n_checks,
                                p_fail_careless, p_fail_attentive,
                                call = call)
    truth <- truth[, c("careless", "pattern", "extent", "onset_item",
                       "offset_item", "speeded", "params")]
    list(metadata = attentive$metadata, truth = truth,
         responses = responses, cells = cells, direct = direct)
  }
  # Non-NULL seed runs the draw under a local seed and restores the caller's RNG on exit;
  # NULL draws from the ambient stream. Nothing below draws RNG.
  result <- with_local_seed(seed, run)
  generator <- c(result$metadata,
                 list(timing = resolved_timing,
                      p_fail_careless = p_fail_careless,
                      p_fail_attentive = p_fail_attentive,
                      seed = seed))
  obj <- new_cier_sim(
    responses = result$responses, items = items,
    seconds = rowSums(result$cells),
    page_seconds = sim_page_totals(result$cells, pages),
    items_per_page = pages,
    checks = result$direct$checks, pass = result$direct$pass,
    truth = result$truth, generator = generator
  )
  validate_cier_sim(obj, call = call)
  obj
}
