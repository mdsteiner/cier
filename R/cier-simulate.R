# Purpose: The public cier_simulate() orchestrator. Wires the attentive GRM
#          layer (R/sim-attentive.R), the careless plan + content engine
#          (R/sim-plan.R, R/sim-patterns.R), the span-coupled response times
#          (R/sim-times.R), and the optional attention checks (R/sim-direct.R)
#          into a validated cier_sim object (R/cier-sim.R).
# Args:    See the public documentation below.
# Returns: A `cier_sim` passing validate_cier_sim().
# Invariants:
#   - No statistical logic here: every draw happens in a kernel. The frozen
#     orchestration order is attentive -> plan -> content -> times -> checks
#     (attentive first, so the clean layer is invariant to `prevalence` at a
#     fixed seed); reordering it reorders every seeded draw and must be a
#     deliberate, test-recorded act.
#   - `seed` applies locally (the caller's RNG state is saved and restored),
#     so a seeded call is reproducible without disturbing the session stream.

# The default pattern mixture: equal weight on every recognized pattern,
# following the equal-allocation simulation design of Welz & Alfons (2023),
# the one published precedent for mixing careless styles in one sample. No
# study estimates a real-world mixture (Voss 2024 and Schroeders et al. 2022
# evaluate one style per design cell), so any unequal weighting would be an
# invented number.
sim_default_patterns <- function() {
  stats::setNames(rep(1 / 8, 8L), sim_pattern_names())
}

# Validate the per-pattern knob list at the public boundary: a named list
# whose names are a subset of the active pattern weights, each entry itself a
# list. (sim_build_plan() records the knobs into the truth frame; the content
# engine reads them back from there, so a typo here would otherwise plant
# nothing and record nothing -- silently.)
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
  for (nm in nms) {
    if (!is.list(pattern_params[[nm]])) {
      cier_abort(
        "cier_error_input",
        "{.arg pattern_params${nm}} must be a named list of knobs.",
        data = list(arg = "pattern_params", observed = nm), call = call
      )
    }
  }
  invisible(NULL)
}

#' Simulate survey responses with planted careless responding
#'
#' Generates a complete survey dataset -- raw responses, response times, and
#' optional attention checks -- in which a known share of respondents answer
#' carelessly with a known pattern, extent, and onset. Attentive responses
#' come from a graded response model over the scales in `items`; careless
#' rows are then overwritten on their planted span and answer faster there.
#' The returned object feeds the shipped indices directly:
#' `cier_screen(x$responses, x$items)`, [cier_total_time()]`(x$seconds)`,
#' [cier_page_time()]`(x$page_seconds, x$items_per_page)`, and
#' [cier_attention()]`(x$checks, x$pass)`, each scored against `x$truth`.
#'
#' @details
#' **What this is for -- and what it is not.** `cier_simulate()` exists for
#' (a) power analysis and study design, (b) method comparison on *planted*
#' patterns, (c) recovery and property tests, and (d) replicating published
#' simulation designs. It is **not** evidence of real-world validity:
#' simulated carelessness is more stylized than even instructed carelessness,
#' and Schroeders et al. (2022) found that a detector's simulation-trained
#' advantage did not transfer to real data. Anchor every real-world accuracy
#' claim to labeled real data instead.
#'
#' **Attentive model.** A graded response model in the factor-analytic
#' parameterization: one factor per scale in `items`, default loadings 0.7,
#' identity factor correlations, latent traits from `trait_distribution`.
#' Reverse-keyed items load negatively, so the returned matrix is **raw**
#' (as-clicked) -- keying-sensitive indices reverse-score internally from
#' `items`, exactly as with field data. Item thresholds come from `marginals`
#' (named presets or explicit per-item / per-scale probability vectors;
#' default `"peaked"`, a symmetric triangular distribution) or from raw
#' `thresholds`; the two are mutually exclusive.
#'
#' **Patterns.** `patterns` weights the eight recognized styles: `random`,
#' `straightline` (anchor / switch knobs), `midpoint`, `extreme`, `diagonal`
#' (with a bounce option), `alternating`, `markov`, and `speeder` (attentive
#' *content* but careless *times* -- content indices should miss it, timing
#' should catch it). The default weights every style equally, following the
#' equal-allocation design of Welz & Alfons (2023); no published study
#' estimates a real-world style mixture, so the default makes no prevalence
#' claim. Per-style knobs go in `pattern_params` (e.g.
#' `list(straightline = list(anchor = "value", value = 5))`) and are recorded
#' verbatim in `truth$params` -- the truth always documents exactly what was
#' applied.
#'
#' **Extent and onset.** Each careless row independently draws an extent:
#' `full` (the whole row), `partial` (careless from a sampled onset to the
#' end), or `temporary` (careless on a bounded window, then recovered). The
#' careless span is the **inclusive** column range `[onset_item,
#' offset_item]`: full is `[1, p]`; partial samples its onset from
#' `onset_window` (default `round(c(0.3, 0.8) * p)`, the sampled-onset design
#' of Welz & Alfons 2023) and runs to `p`; temporary additionally samples an
#' offset of at most `p - 1`, so a recovery tail of at least one item always
#' follows (its onset is sampled below the window's upper bound where needed,
#' and the window may be as short as one item). `temporary` therefore needs
#' at least three items. Both the content overwrite and the careless time
#' shift are confined to that span.
#'
#' **Response times.** Lognormal per-item times with a per-respondent pace
#' intercept (`respondent_sd`), so the attentive and careless distributions
#' overlap -- by design, nothing separates perfectly. Careless rows shift to
#' the careless log-mean *inside their span only*. The defaults (attentive
#' median 8 s/item; careless median 1.5 s/item, below the 2 s/item per-page
#' screening floor of Bowling et al. 2023; `sigma = 0.5`;
#' `respondent_sd = 1.2`) are **documented conventions anchored to cited
#' thresholds, not estimates**, calibrated so the planted speeding is
#' detectable by [cier_page_time()]'s default rule while the total-time
#' rank-AUC stays inside the .66-.92 range reported across real datasets
#' (Goldammer et al. 2024; Huang et al. 2012). Every parameter is overridable
#' through `timing`. The truth's `speeded` column records which rows carry
#' careless times; in this version that is every careless row, including
#' speeders.
#'
#' **Attention checks.** With `n_checks > 0`, separate five-point
#' instructed-response checks are generated (never written into
#' `$responses`): careless rows fail each check with `p_fail_careless`,
#' attentive rows with `p_fail_attentive`. The defaults 0.75 / 0.05 are
#' documented conventions -- no published point estimates exist; the
#' bracketing evidence is Goldammer et al.'s (2024) bogus-item AUCs of
#' .66-.76 and Bruhlmann et al.'s (2020) 92 bogus-item failures among 394
#' respondents. Failure is probabilistic on both sides, so direct checks
#' never separate the groups perfectly either.
#'
#' **Reproducibility.** With `seed`, the generator runs on a local RNG state
#' (the session stream is saved and restored), and the same seed reproduces
#' the returned object bytewise.
#'
#' @param n Positive whole number. Number of respondents.
#' @param items A data frame describing the battery, one row per item:
#'   `scale` (character, defines the factor structure; a single scale is
#'   fine), `max` (largest response option, required on every item), and
#'   optionally `min` (smallest option, default 1) and `reverse_keyed`
#'   (logical). Heterogeneous ranges are allowed. Returned untouched as
#'   `$items`, ready for the index wrappers.
#' @param prevalence Number in `[0, 1]`. Share of respondents made careless;
#'   exactly `round(n * prevalence)` rows are drawn.
#' @param patterns Named numeric vector of pattern weights summing to 1, over
#'   `random`, `straightline`, `midpoint`, `extreme`, `diagonal`,
#'   `alternating`, `markov`, `speeder`. `NULL` (default) weights all eight
#'   equally.
#' @param pattern_params Named list of per-pattern knob lists, forwarded to
#'   the pattern generators and recorded in `truth$params`. Names must appear
#'   in `patterns`.
#' @param prop_partial,prop_temporary Numbers in `[0, 1]` summing to at most
#'   1. Shares of the careless rows whose extent is `partial` /
#'   `temporary`; the remainder are `full`. Both default to 0.
#' @param onset_window Optional two whole numbers `c(lo, hi)` inside
#'   `[1, p]`: the sampling range for partial / temporary onsets. `NULL`
#'   (default) uses `round(c(0.3, 0.8) * p)`.
#' @param loadings,factor_cor,thresholds Optional GRM structure overrides: a
#'   `p x m` loading matrix, an `m x m` factor correlation matrix, and a
#'   per-item list of strictly increasing threshold vectors. `thresholds` is
#'   mutually exclusive with `marginals`.
#' @param marginals Target observed response distribution(s): a preset name
#'   (`"peaked"`, `"uniform"`, `"skewed_right"`, `"skewed_left"`,
#'   `"bimodal"`), an explicit probability vector, or a per-scale / per-item
#'   specification. Default `"peaked"`.
#' @param trait_distribution One of `"normal"`, `"skew_normal"`, `"t"`,
#'   `"bimodal"` for the latent traits.
#' @param trait_params Named list of distribution parameters (e.g.
#'   `list(df = 5)` for `"t"`).
#' @param timing Named list overriding any of the response-time parameters
#'   `mu_att`, `mu_car` (log-mean seconds per item, attentive / careless),
#'   `sigma` (cell-level log-SD), `respondent_sd` (pace-intercept SD), and
#'   `min_seconds` (strictly positive floor). Unknown names are an error.
#' @param items_per_page Page structure for `$page_seconds`: `NULL` (default)
#'   puts each contiguous run of `items$scale` on its own page; a single
#'   whole number chunks uniform pages of that many items (the last page
#'   takes the remainder); a vector of whole positive counts summing to the
#'   item count is used verbatim.
#' @param n_checks Non-negative whole number of attention checks to generate
#'   (default 0: `$checks` and `$pass` are `NULL`).
#' @param p_fail_careless,p_fail_attentive Numbers in `[0, 1]`: per-check
#'   failure probabilities for careless and attentive rows.
#' @param seed Optional single whole number. When supplied, the generator is
#'   seeded locally and the result reproduces bytewise.
#'
#' @return A `cier_sim` object: a list with `$responses` (n x p raw integer
#'   matrix), `$items` (the input frame, untouched), `$seconds` (total
#'   completion time per respondent), `$page_seconds` (per-page totals) with
#'   `$items_per_page`, `$checks` and `$pass` (or `NULL`), and `$truth` -- one
#'   row per respondent with `careless`, `pattern` (`"attentive"` on clean
#'   rows), `extent` (`"none"` / `"full"` / `"partial"` / `"temporary"`),
#'   `onset_item` and `offset_item` (the inclusive careless span; `NA` on
#'   attentive rows), `speeded`, and `params` (the applied pattern knobs).
#'   Generator metadata (loadings, thresholds, residual variances, resolved
#'   timing parameters, failure probabilities, seed) is attached as the
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
#' Goldammer, P., Annen, H., Stockli, P. L., & Jonas, K. (2024). On the
#' utility of indirect methods for detecting faking and careless responding:
#' A comparison under experimentally induced and naturally occurring
#' conditions. *Behavior Research Methods*, 56(8), 8422-8449.
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
#' @seealso [cier_screen()] to score the simulated responses;
#'   [cier_total_time()], [cier_page_time()], and [cier_attention()] for the
#'   timing and direct slots.
#' @export
#' @examples
#' items <- data.frame(scale = rep(c("E", "A", "C"), each = 4L), max = 5)
#' sim <- cier_simulate(n = 120, items, prevalence = 0.2, seed = 2026)
#' sim
#'
#' # Score the planted rows with a content index and the timing family:
#' table(flagged = cier_longstring(sim$responses)$flagged,
#'       truth = sim$truth$careless)
#' cier_total_time(sim$seconds)
cier_simulate <- function(n, items, prevalence = 0.2, patterns = NULL,
                          pattern_params = list(), prop_partial = 0,
                          prop_temporary = 0, onset_window = NULL,
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
  if (is.null(patterns)) {
    patterns <- sim_default_patterns()
  }
  patterns <- sim_validate_patterns(patterns, call)
  sim_check_pattern_params(pattern_params, patterns, call)
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
    cells <- sim_times_lognormal(n, p, truth$careless, truth$onset_item,
                                 truth$offset_item, timing = resolved_timing,
                                 call = call)
    direct <- sim_direct_checks(n, truth$careless, n_checks,
                                p_fail_careless, p_fail_attentive,
                                call = call)
    # speeded records the timing fact (this row's span carries careless
    # times); in this version that is every careless row, kept as its own
    # column so the schema survives a future arm that decouples them.
    truth$speeded <- truth$careless
    truth <- truth[, c("careless", "pattern", "extent", "onset_item",
                       "offset_item", "speeded", "params")]
    list(metadata = attentive$metadata, truth = truth,
         responses = responses, cells = cells, direct = direct)
  }
  result <- if (is.null(seed)) {
    run()
  } else {
    saved <- globalenv()[[".Random.seed"]]  # NULL when no draw happened yet
    on.exit(restore_random_seed(saved), add = TRUE)
    set.seed(seed)
    run()
  }
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
