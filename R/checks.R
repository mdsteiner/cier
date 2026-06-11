# Purpose: Minimal user-input checks used by the public wrappers. Raise
#          `cier_error_input` (a bad user argument) rather than the
#          `cier_error_data` / `cier_error_state` raised by the schema
#          validators.
# Args:    The value under inspection plus its argument name.
# Returns: The input (invisibly) on success.
# Invariants: A failure always raises a typed cier condition.

# The check_* helpers wrap checkmate's predicate functions (which return TRUE or
# a message string) so the diagnostics stay typed cier_error_input conditions
# with the package's cli phrasing, rather than checkmate's plain assert errors.

check_string <- function(x, arg, call = rlang::caller_env()) {
  if (!isTRUE(checkmate::check_string(x, min.chars = 1L))) {
    cier_abort(
      "cier_error_input",
      "{.arg {arg}} must be a non-empty character string.",
      data = list(arg = arg), call = call
    )
  }
  invisible(x)
}

check_choice <- function(x, arg, choices, call = rlang::caller_env()) {
  if (!isTRUE(checkmate::check_choice(x, choices))) {
    cier_abort(
      "cier_error_input",
      "{.arg {arg}} must be one of {.val {choices}}.",
      data = list(arg = arg, observed = x), call = call
    )
  }
  invisible(x)
}

# Guard the two cutoff-override knobs an index exposes: a rate (`fpr` / `alpha` /
# `frac`) and a literal `cutoff`. They are mutually exclusive -- two ways to set
# the same cutoff -- so accepting both would be ambiguous. Abort when both are
# supplied. (Message stays generic so it reads correctly for every rate name.)
assert_single_override <- function(rate, rate_name, cutoff,
                                   call = rlang::caller_env()) {
  if (!is.null(rate) && !is.null(cutoff)) {
    cier_abort(
      "cier_error_input",
      c("Supply only one of {.arg {rate_name}} and {.arg cutoff}.",
        "i" = "They are two ways to set the same cutoff."),
      data = list(args = c(rate_name, "cutoff")), call = call
    )
  }
  invisible(NULL)
}

# Guard a set of mutually-exclusive cutoff-override knobs -- the n-way form of
# assert_single_override() for an index that exposes more than two. cier_total_time
# has three (`fpr`, `frac_median`, `cutoff`), each a different way to set the same
# cutoff. `args` is a named list of the knob values (NULL where unset); abort when
# more than one is supplied, naming the offenders so the message lists exactly the
# pair (or triple) given.
assert_single_cutoff <- function(args, call = rlang::caller_env()) {
  supplied <- names(args)[!vapply(args, is.null, logical(1L))]
  if (length(supplied) > 1L) {
    cier_abort(
      "cier_error_input",
      c("Supply only one of {.arg {supplied}}.",
        "i" = "They are different ways to set the same cutoff."),
      data = list(args = supplied), call = call
    )
  }
  invisible(NULL)
}

check_number <- function(x, arg, lower = -Inf, upper = Inf,
                         call = rlang::caller_env()) {
  if (!isTRUE(checkmate::check_number(x, lower = lower, upper = upper,
                                      finite = TRUE))) {
    cier_abort(
      "cier_error_input",
      "{.arg {arg}} must be a single finite number in [{lower}, {upper}].",
      data = list(arg = arg, observed = x), call = call
    )
  }
  invisible(x)
}

# A target tail probability / significance level: a single finite number
# strictly inside the open interval (0, 1) (e.g. `alpha`, `fpr`). checkmate
# guards the closed [0, 1]; `&& x > 0 && x < 1` opens both ends, so 0 (a cutoff
# that flags nobody) and 1 (a cutoff that flags everybody) are rejected.
check_open_unit <- function(x, arg, call = rlang::caller_env()) {
  is_open_unit <- isTRUE(
    checkmate::check_number(x, lower = 0, upper = 1, finite = TRUE)
  ) && x > 0 && x < 1
  if (!is_open_unit) {
    cier_abort(
      "cier_error_input",
      "{.arg {arg}} must be a single number in the open interval (0, 1).",
      data = list(arg = arg, observed = x), call = call
    )
  }
  invisible(x)
}

# A fraction of the item count: a single finite number in the half-open
# interval (0, 1] (0 would flag a zero-length run; values above 1 exceed the
# item count). checkmate guards the closed [0, 1]; `&& x > 0` opens the lower end.
check_fraction <- function(x, arg, call = rlang::caller_env()) {
  is_fraction <- isTRUE(
    checkmate::check_number(x, lower = 0, upper = 1, finite = TRUE)
  ) && x > 0
  if (!is_fraction) {
    cier_abort(
      "cier_error_input",
      "{.arg {arg}} must be a single number in the interval (0, 1].",
      data = list(arg = arg, observed = x), call = call
    )
  }
  invisible(x)
}

# A single logical flag: one non-missing TRUE/FALSE. Used for the boolean
# wrapper switches (e.g. `resample`).
check_flag <- function(x, arg, call = rlang::caller_env()) {
  if (!isTRUE(checkmate::check_flag(x))) {
    cier_abort(
      "cier_error_input",
      "{.arg {arg}} must be a single {.val TRUE} or {.val FALSE}.",
      data = list(arg = arg), call = call
    )
  }
  invisible(x)
}

# A count: a single positive whole number (>= 1, no NA). Used for the resampling
# iteration count (`n_resamples`); checkmate's integerish tolerance rejects a
# non-whole number such as 2.5.
check_count <- function(x, arg, call = rlang::caller_env()) {
  if (!isTRUE(checkmate::check_count(x, positive = TRUE))) {
    cier_abort(
      "cier_error_input",
      "{.arg {arg}} must be a single positive whole number.",
      data = list(arg = arg, observed = x), call = call
    )
  }
  invisible(x)
}

# A single integer (integerish, finite, any sign, no NA). Used for an RNG `seed`,
# where a fractional value would be silently truncated by `set.seed()`;
# checkmate's integerish tolerance rejects a non-whole number such as 1.9.
check_int <- function(x, arg, call = rlang::caller_env()) {
  if (!isTRUE(checkmate::check_int(x))) {
    cier_abort(
      "cier_error_input",
      "{.arg {arg}} must be a single integer.",
      data = list(arg = arg, observed = x), call = call
    )
  }
  invisible(x)
}

# Validate the two mutually-exclusive cutoff-override knobs every percentile
# index exposes: `fpr` (a target false-positive tail mass in the open unit
# interval) and a literal `cutoff` on the score (a finite number in
# `[lower, upper]` -- the score's natural range, e.g. `[-1, 1]` for a correlation
# or `[0, Inf)` for an SD). The shared front-end for cier_irv / cier_even_odd /
# cier_person_total / cier_personal_reliability, called before the kernel runs so
# a bad argument fails early.
check_percentile_overrides <- function(fpr, cutoff, lower = -Inf, upper = Inf,
                                       call = rlang::caller_env()) {
  if (!is.null(fpr)) check_open_unit(fpr, "fpr", call = call)
  if (!is.null(cutoff)) {
    check_number(cutoff, "cutoff", lower = lower, upper = upper, call = call)
  }
  assert_single_override(fpr, "fpr", cutoff, call = call)
  invisible(NULL)
}

# Coerce a user response payload to a validated numeric matrix. A data.frame or
# tibble is accepted and coerced (so users never call as.matrix()); `NA` is the
# only allowed missing marker -- NaN and infinities are rejected. This is the
# function-first input contract shared by every index wrapper.
check_responses <- function(responses, arg = "responses",
                            call = rlang::caller_env()) {
  if (!is.matrix(responses) && !is.data.frame(responses)) {
    cier_abort(
      "cier_error_input",
      c("{.arg {arg}} must be a matrix or data.frame (one row per respondent).",
        "x" = "Got {.cls {class(responses)}}; a bare vector or array is rejected."),
      data = list(arg = arg), call = call
    )
  }
  m <- as.matrix(responses)
  if (!is.numeric(m)) {
    cier_abort(
      "cier_error_input",
      "{.arg {arg}} must be numeric (a matrix or data.frame of numbers).",
      data = list(arg = arg), call = call
    )
  }
  if (nrow(m) == 0L || ncol(m) == 0L) {
    cier_abort(
      "cier_error_input",
      "{.arg {arg}} must have at least one row and one column.",
      data = list(arg = arg, nrow = nrow(m), ncol = ncol(m)), call = call
    )
  }
  if (any(is.nan(m) | is.infinite(m))) {
    cier_abort(
      "cier_error_input",
      c("{.arg {arg}} may contain only finite numbers or {.val NA}.",
        "x" = "Found {.val NaN} or infinite values."),
      data = list(arg = arg), call = call
    )
  }
  m
}

# Coerce and validate a per-respondent completion-time vector. cier_total_time
# takes a BARE numeric vector -- one total in seconds per respondent, what survey
# platforms export -- not a response matrix: a 2-D / data.frame input is ambiguous
# about which axis is respondents, so it is rejected with a fix-it hint. Every
# observed (non-NA) time must be finite and strictly positive (a duration cannot be
# zero or negative; NaN / infinite are rejected); NA is permitted and abstains.
check_seconds <- function(seconds, arg = "seconds", call = rlang::caller_env()) {
  if (!is.numeric(seconds) || !is.null(dim(seconds))) {
    cier_abort(
      "cier_error_input",
      c("{.arg {arg}} must be a numeric vector (one completion time per respondent).",
        "x" = "Got {.cls {class(seconds)}}; a 2-D or non-numeric input is rejected.",
        "i" = "Sum per-cell response times to one total per respondent first."),
      data = list(arg = arg), call = call
    )
  }
  if (length(seconds) == 0L) {
    cier_abort(
      "cier_error_input",
      "{.arg {arg}} must have at least one element.",
      data = list(arg = arg), call = call
    )
  }
  # is.nan() / is.infinite() distinguish NaN and +/-Inf from NA -- is.na() is TRUE
  # for NaN, so stripping NA first would hide a NaN. Check the raw vector before
  # dropping genuine NA (which abstains).
  if (any(is.nan(seconds) | is.infinite(seconds))) {
    cier_abort(
      "cier_error_input",
      c("{.arg {arg}} may contain only finite positive seconds or {.val NA}.",
        "x" = "Found {.val NaN} or infinite values."),
      data = list(arg = arg), call = call
    )
  }
  obs <- seconds[!is.na(seconds)]
  if (any(obs <= 0)) {
    cier_abort(
      "cier_error_input",
      c("{.arg {arg}} must be strictly positive completion times in seconds.",
        "x" = "Found zero or negative values; a completion time is greater than zero.",
        "i" = "Recode missing durations to {.val NA}; they abstain."),
      data = list(arg = arg), call = call
    )
  }
  as.numeric(seconds)
}

# Coerce and validate the per-page completion-time matrix cier_page_time takes:
# an n x pages numeric matrix (data.frame / tibble coerced) of per-page TOTAL
# times in seconds, one column per page. Reuses check_responses() for the
# matrix / numeric / finite / non-empty contract (a bare vector is rejected --
# which axis is the respondent is ambiguous), then adds the timing-specific rule
# that every observed (non-NA) cell is strictly positive (a page time cannot be
# zero or negative; NA marks an untimed page and contributes no evidence).
check_page_seconds <- function(page_seconds, arg = "page_seconds",
                               call = rlang::caller_env()) {
  m <- check_responses(page_seconds, arg = arg, call = call)
  # check_responses has already rejected NaN / infinite, so the only non-finite
  # entries left are NA (untimed pages); na.rm skips them. any() over no observed
  # cells (an all-NA matrix) is FALSE -- it passes here and abstains in the
  # kernel. Scans in place rather than materialising the observed-cell subset.
  if (any(m <= 0, na.rm = TRUE)) {
    cier_abort(
      "cier_error_input",
      c("{.arg {arg}} must be strictly positive page times in seconds.",
        "x" = "Found zero or negative values; a page time is greater than zero.",
        "i" = "Recode untimed pages to {.val NA}; they contribute no evidence."),
      data = list(arg = arg), call = call
    )
  }
  m
}

# Validate the per-page item-count vector cier_page_time takes: a plain numeric
# vector, one entry per column of `page_seconds`, every entry a positive whole
# number (the page total is divided by it to get the mean per-item time). NA,
# fractional, non-positive, non-numeric, wrong-length, and 2-D inputs all fail.
# Returns the vector as integer on success.
check_items_per_page <- function(items_per_page, n_pages,
                                 arg = "items_per_page",
                                 call = rlang::caller_env()) {
  ok <- is.numeric(items_per_page) && is.null(dim(items_per_page)) &&
    length(items_per_page) == n_pages && is_finite_whole(items_per_page) &&
    all(items_per_page >= 1)
  if (!ok) {
    cier_abort(
      "cier_error_input",
      c("{.arg {arg}} must give the item count for every page.",
        "x" = "Need a length-{n_pages} vector of positive whole numbers (one \\
               per column of {.arg page_seconds}).",
        "i" = "An untimed-but-present page still has its item count."),
      data = list(arg = arg, expected = n_pages), call = call
    )
  }
  as.integer(items_per_page)
}

# Validate the per-check pass-set spec cier_attention takes: a list with one
# element per column of `checks`, positional (pass[[j]] applies to column j). A
# respondent fails check j when their answered response is NOT in pass[[j]], so
# each pass-set must be a non-empty vector of finite passing values: an empty
# pass-set would fail everyone (a likely user mistake), and an NA / NaN / infinite
# passing value cannot be a real response code (the data side rejects them too).
# Returns the list unchanged on success.
check_pass <- function(pass, n_checks, arg = "pass", call = rlang::caller_env()) {
  if (!is.list(pass) || length(pass) != n_checks) {
    cier_abort(
      "cier_error_input",
      c("{.arg {arg}} must be a list with one pass-set per check column.",
        "x" = "Need a length-{n_checks} list (one element per column of \\
               {.arg checks})."),
      data = list(arg = arg, expected = n_checks), call = call
    )
  }
  ok <- all(vapply(
    pass,
    function(p) is.numeric(p) && length(p) >= 1L && all(is.finite(p)),
    logical(1L)
  ))
  if (!ok) {
    cier_abort(
      "cier_error_input",
      c("Each {.arg {arg}} element must be a non-empty vector of finite \\
         passing values.",
        "i" = "{.arg {arg}}[[j]] lists the value(s) that pass check j; an \\
               answered response outside it is a failure."),
      data = list(arg = arg), call = call
    )
  }
  pass
}

# Require integer-coded responses: every non-NA value must be a whole number.
# The Markov pattern index (cier_lazr) scores a transition matrix over discrete
# response anchors, so a continuous / averaged (POMP) score has no Markov chain;
# a typed error is clearer than silently dropping the non-integer transitions.
# Run AFTER check_responses (so NaN / infinite values are caught there first);
# the non-NA values are then finite numeric, so is_finite_whole() reduces to the
# wholeness test -- reused here to keep the integer definition in one place.
check_integer_responses <- function(m, arg = "responses",
                                    call = rlang::caller_env()) {
  obs <- m[!is.na(m)]
  if (length(obs) > 0L && !is_finite_whole(obs)) {
    cier_abort(
      "cier_error_input",
      c("{.arg {arg}} must be integer-coded responses (whole-number anchors).",
        "x" = "Found non-integer values; {.fun cier_lazr} scores a Markov chain \\
               over discrete response anchors.",
        "i" = "Recode to integer anchors, or use an index that accepts \\
               continuous scores such as {.fun cier_irv}."),
      data = list(arg = arg), call = call
    )
  }
  invisible(m)
}

# `scale`: required, character-coercible, every item labelled, >= min_scales
# distinct labels (the split-half indices correlate across scales).
check_items_scale <- function(items, min_scales, arg, call) {
  if (is.null(items$scale)) {
    cier_abort(
      "cier_error_input",
      "{.arg {arg}} must have a {.field scale} column (one label per item).",
      data = list(arg = arg), call = call
    )
  }
  scale <- as.character(items$scale)
  if (anyNA(scale) || any(!nzchar(scale))) {
    cier_abort(
      "cier_error_input",
      "{.field scale} must be a non-missing label for every item.",
      data = list(arg = arg), call = call
    )
  }
  n_scales <- length(unique(scale))
  if (n_scales < min_scales) {
    cier_abort(
      "cier_error_input",
      c("{.arg {arg}} must have at least {min_scales} distinct {.field scale}s.",
        "x" = "Found {n_scales}.",
        "i" = "The split-half consistency indices correlate across scales."),
      data = list(arg = arg, observed = n_scales, expected = min_scales),
      call = call
    )
  }
  scale
}

# `reverse_keyed`: optional logical, length n_items, no NA; defaults all-FALSE.
check_items_reverse <- function(items, n_items, arg, call) {
  rk <- items$reverse_keyed
  if (is.null(rk)) {
    return(rep(FALSE, n_items))
  }
  if (!is.logical(rk) || anyNA(rk)) {
    cier_abort(
      "cier_error_input",
      "{.field reverse_keyed} must be a logical value (no {.val NA}) per item.",
      data = list(arg = arg), call = call
    )
  }
  rk
}

# Predicate: a numeric vector of finite whole numbers (no NA / NaN / Inf, every
# element integer-valued). The shared core of every items-column validator --
# `max`, `min`, and their person-fit variants each layer their own bound or
# scope on top of it. `is.numeric()` is first so a non-numeric (or NULL) input
# short-circuits to FALSE before the finiteness / rounding predicates run.
is_finite_whole <- function(v) {
  is.numeric(v) && all(is.finite(v)) && all(v == round(v))
}

# `max`: the largest response option per item. Required -- a finite whole
# number of at least `min + 1` (two response options) -- only on items that are
# reverse-keyed (so they can be reverse-scored with the self-inverse reflection
# (min + max) - x; `mins` is the RESOLVED base from check_items_min, default 1);
# NA is permitted on forward items and the column may be absent entirely when
# nothing is reverse-keyed. Returns the column unchanged (or NULL when absent).
check_items_max <- function(items, reverse_keyed, mins, arg, call) {
  maxs <- items$max
  if (!any(reverse_keyed)) {
    return(maxs)
  }
  rev_maxs <- if (is.null(maxs)) NA_real_ else maxs[reverse_keyed]
  # is_finite_whole() rejects NA, NaN, and Inf (a non-finite maximum would
  # reflect to (min + Inf) - x and poison the reverse columns) and non-integers.
  # The bound is min + 1, NOT an absolute 2: a 0/1 item (min 0, max 1) is the
  # smallest valid scale.
  ok <- is_finite_whole(rev_maxs) && all(rev_maxs >= mins[reverse_keyed] + 1)
  if (!ok) {
    cier_abort(
      "cier_error_input",
      c("Reverse-keyed items need an integer {.field max} of at least \\
         {.field min} + 1.",
        "i" = "Set {.field max} (the largest response option) for every \\
               reverse-keyed item; {.field min} (the smallest, default 1) is \\
               read alongside it for the reflection (min + max) - x."),
      data = list(arg = arg), call = call
    )
  }
  maxs
}

# `min`: the response-scale minimum (base). Optional; defaults to 1 (the
# 1..max coding) when the column is absent. When supplied it generalises
# the reverse-keying reflection (min + max) - x, so a 0-based or bipolar scale
# reflects onto itself. Validated like max -- a finite whole number on every
# reverse-keyed item (any integer base: 0, negative, and bipolar are allowed;
# there is no lower bound) -- with NA permitted on forward items. Resolved
# BEFORE max, whose >= min + 1 bound reads it. Returns the column, or
# rep(1L, n_items) when absent.
check_items_min <- function(items, reverse_keyed, n_items, arg, call) {
  mins <- items$min
  if (is.null(mins)) {
    return(rep(1L, n_items))
  }
  rev_mins <- mins[reverse_keyed]
  ok <- length(rev_mins) == 0L || is_finite_whole(rev_mins)
  if (!ok) {
    cier_abort(
      "cier_error_input",
      c("Reverse-keyed items need a finite integer {.field min} (the scale base).",
        "i" = "{.field min} is the smallest response option (default 1); set it \\
               for every reverse-keyed item on a 0-based or bipolar scale."),
      data = list(arg = arg), call = call
    )
  }
  mins
}

# The shared item-frame precondition every items validator runs first: `items`
# must be a data.frame with exactly one row per item (column of `responses`).
# Returns invisibly NULL; aborts cier_error_input on a violation.
check_items_frame <- function(items, n_items, arg, call) {
  if (!is.data.frame(items)) {
    cier_abort(
      "cier_error_input",
      c("{.arg {arg}} must be a data.frame (one row per item).",
        "x" = "Got {.cls {class(items)}}."),
      data = list(arg = arg), call = call
    )
  }
  if (nrow(items) != n_items) {
    cier_abort(
      "cier_error_input",
      c("{.arg {arg}} must have one row per item (column of {.arg responses}).",
        "x" = "Got {nrow(items)} item row{?s} for {n_items} column{?s}."),
      data = list(arg = arg, observed = nrow(items), expected = n_items),
      call = call
    )
  }
  invisible(NULL)
}

# Validate the per-item `items` frame the split-half family uses (even-odd and
# personal reliability). `items` is a data.frame with one row per item, aligned
# to the columns of `responses`. Returns a normalized
# list(scale, reverse_keyed, min, max) on success; aborts cier_error_input
# on any malformed field. `max` is NULL when its column is absent
# (permitted only when nothing is reverse-keyed); `min` defaults to all-1.
check_items <- function(items, n_items, min_scales = 2L,
                        arg = "items", call = rlang::caller_env()) {
  check_items_frame(items, n_items, arg, call)
  scale <- check_items_scale(items, min_scales, arg, call)
  reverse_keyed <- check_items_reverse(items, n_items, arg, call)
  minimum <- check_items_min(items, reverse_keyed, n_items, arg, call)
  maximum <- check_items_max(items, reverse_keyed, minimum, arg, call)
  list(scale = scale, reverse_keyed = reverse_keyed,
       min = minimum, max = maximum)
}

# Validate the per-item `items` frame cier_ht() uses. Ht needs item metadata ONLY
# to reverse-score keyed items: mokken::coefH accepts a mix of category counts and
# the kernel never reads `max` directly, so -- unlike the Gnormed bridge --
# the response ranges need not be homogeneous and `max` is required only on
# reverse-keyed items (so they can be reverse-scored), exactly like the
# split-half family but without the `scale` requirement. Returns a normalized
# list(reverse_keyed, min, max); `max` is NULL when its column is absent
# (permitted when nothing is reverse-keyed) and `min` defaults to all-1.
check_items_ht <- function(items, n_items, arg = "items",
                           call = rlang::caller_env()) {
  check_items_frame(items, n_items, arg, call)
  reverse_keyed <- check_items_reverse(items, n_items, arg, call)
  minimum <- check_items_min(items, reverse_keyed, n_items, arg, call)
  maximum <- check_items_max(items, reverse_keyed, minimum, arg, call)
  list(reverse_keyed = reverse_keyed, min = minimum, max = maximum)
}

# `max` for the person-fit (Gnormed) bridge: required on EVERY item (not only
# reverse-keyed ones) -- a finite whole number of at least `min + 1` -- because
# reverse-keying, per-item zero-basing, and the category count
# Ncat = max - min + 1 all read it. Per-item validity is a PLAIN input error;
# the span homogeneity PerFit additionally needs is classified separately in
# check_items_span_homogeneous(). Returns the per-item vector on success.
check_items_max_personfit <- function(items, mins, arg, call) {
  maxs <- items$max
  # is_finite_whole() returns FALSE for an absent (NULL) column, so the missing
  # `max` case is caught here without a separate is.null() guard.
  ok <- is_finite_whole(maxs) && all(maxs >= mins + 1)
  if (!ok) {
    cier_abort(
      "cier_error_input",
      c("{.arg {arg}} needs an integer {.field max} of at least \\
         {.field min} + 1 on every item.",
        "i" = "{.field max} is the largest response option; {.field min} (the \\
               smallest, default 1) is read alongside it."),
      data = list(arg = arg), call = call
    )
  }
  maxs
}

# The span homogeneity PerFit's polytomous statistics require: ONE number of
# response categories, Ncat = max - min + 1, across all items (items may still
# differ in base: 1..5 and 0..4 both have five options). A heterogeneous span
# on otherwise-valid metadata is NOT a malformed frame -- it is accurate
# metadata for genuinely mixed-format data that the PerFit backend cannot
# score -- so the abort carries the cier_error_backend_limit subclass
# (mirroring mokken's 10-category ceiling in kernel_ht) and cier_screen()
# records the index as skipped-with-reason instead of aborting the battery.
# It takes the validated `maxs` returned by check_items_max_personfit (not the
# raw column), so a malformed max can never reach the span arithmetic and be
# misclassified as a backend limit. Returns the shared category count Ncat
# (max - min + 1), which the wrapper threads to the kernel -- the single place
# Ncat is derived.
check_items_span_homogeneous <- function(maxs, mins, arg, call) {
  spans <- maxs - mins
  if (length(unique(spans)) != 1L) {
    cier_abort(
      c("cier_error_backend_limit", "cier_error_input"),
      c("Gnormed cannot score items with different numbers of response \\
         categories.",
        "x" = "Observed numbers of options ({.field max} - {.field min} + 1): \\
               {.val {sort(unique(spans)) + 1}}.",
        "i" = "The PerFit backend scores one Ncat across all items; screen \\
               homogeneous item subsets separately, or use another index."),
      data = list(arg = arg, observed = sort(unique(spans)) + 1,
                  reason = "mixed response-category counts (PerFit scores a single Ncat)"),
      call = call
    )
  }
  as.integer(spans[[1L]] + 1L)
}

# `min` for the person-fit backends: the scale base, used to reverse-key AND to
# zero-base EVERY item (unlike the split-half family, which keys only reverse
# items), so it must be a finite whole number on every item. Optional; defaults
# to all-1 (the 1..max coding) when the column is absent.
check_items_min_personfit <- function(items, n_items, arg, call) {
  mins <- items$min
  if (is.null(mins)) {
    return(rep(1L, n_items))
  }
  ok <- is_finite_whole(mins)
  if (!ok) {
    cier_abort(
      "cier_error_input",
      c("{.field min} (the scale base) must be a finite integer on every item.",
        "i" = "{.field min} is the smallest response option (default 1); set it \\
               for a 0-based or bipolar scale."),
      data = list(arg = arg), call = call
    )
  }
  mins
}

# Validate the per-item `items` frame the Gnormed bridge uses. Unlike the
# split-half family this needs `max` on EVERY item, a homogeneous span
# (PerFit is single-Ncat), and does NOT use `scale`. Returns a normalized
# list(reverse_keyed, min, max, ncat); aborts cier_error_input on any malformed
# field, with the heterogeneous-span case additionally carrying
# cier_error_backend_limit (see check_items_span_homogeneous). `ncat` is the
# shared category count, computed once here and threaded to the kernel so it is
# never re-derived. The span check takes the VALIDATED `maximum`, so the
# per-item-validity-before-homogeneity order is enforced by the data flow.
check_items_personfit <- function(items, n_items, arg = "items",
                                  call = rlang::caller_env()) {
  check_items_frame(items, n_items, arg, call)
  reverse_keyed <- check_items_reverse(items, n_items, arg, call)
  minimum <- check_items_min_personfit(items, n_items, arg, call)
  maximum <- check_items_max_personfit(items, minimum, arg, call)
  ncat <- check_items_span_homogeneous(maximum, minimum, arg, call)
  list(reverse_keyed = reverse_keyed, min = minimum, max = maximum, ncat = ncat)
}

# Thin, mockable wrapper around requireNamespace() so tests can simulate an
# absent optional backend via testthat::local_mocked_bindings().
cier_namespace_present <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

# Abort with a typed input error when an optional backend package an index needs
# (a `Suggests` dependency) is not installed. Shared by the person-fit bridges
# (PerFit for Gnormed, mokken for Ht): those indices cannot compute their
# statistic without the backend, so absence is a precondition error the caller
# fixes by installing the package.
require_suggested <- function(pkg, fn, call = rlang::caller_env()) {
  if (!cier_namespace_present(pkg)) {
    cier_abort(
      "cier_error_input",
      c("{.pkg {pkg}} is required for {.fn {fn}} but is not installed.",
        "i" = "Install it with {.code install.packages(\"{pkg}\")}."),
      data = list(arg = pkg, fn = fn), call = call
    )
  }
  invisible(NULL)
}
