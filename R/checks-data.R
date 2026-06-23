# Coercers / validators for user response payloads (response matrix, completion-time
# vector and matrix, attention-check pass-sets). Abort with `cier_error_input` on bad
# input; return the coerced payload on success. Scalar / argument validators are in
# R/checks-args.R; `items`-frame validators in R/checks-items.R.

# Logical mask of data.frame columns to reject on (and, where named, to name). A logical
# column promotes silently to 0/1 in as.matrix(), so a non-numeric column with an observed
# value offends even when the matrix stays numeric. Non-data.frame: empty mask (no types).
nonnumeric_offenders <- function(responses, matrix_is_numeric) {
  if (!is.data.frame(responses)) {
    return(logical(0L))
  }
  nonnum <- !vapply(responses, is.numeric, logical(1L))
  if (matrix_is_numeric) {
    observed <- vapply(responses, function(col) any(!is.na(col)), logical(1L))
    nonnum & observed
  } else {
    nonnum
  }
}

# Coerce a response payload to a validated numeric matrix (data.frame / tibble coerced).
# `NA` is the only allowed missing marker (NaN / Inf rejected); an all-NA payload has no
# observed responses and abstains. Any observed non-numeric value is rejected. `arg` names
# the failing argument in the hint.
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
  # An all-NA payload (logical-typed) holds no observed responses: coerce to double so it
  # joins the abstention path instead of tripping the numeric gate.
  if (!is.numeric(m) && length(m) > 0L && all(is.na(m))) {
    storage.mode(m) <- "double"
  }
  numeric_m <- is.numeric(m)
  # A data.frame coerces to a character matrix once any column is non-numeric, and a
  # logical column coerces silently to 0/1 -- inspect source columns to catch both.
  # An offender with a blank / NA header is still rejected, only its printed name dropped.
  offenders <- nonnumeric_offenders(responses, numeric_m)
  if (!numeric_m || any(offenders)) {
    nm <- if (is.data.frame(responses)) names(responses)[offenders] else character(0L)
    bad_cols <- nm[!is.na(nm) & nzchar(nm)]
    head_line <- "{.arg {arg}} must be numeric (a matrix or data.frame of numbers)."
    name_line <- if (length(bad_cols) > 0L) {
      c("x" = "Non-numeric column{?s}: {.val {bad_cols}}.")
    }
    cier_abort(
      "cier_error_input",
      c(head_line, name_line,
        "i" = "An ID, label, or free-text column mixed in? Pass only the data \\
               columns (e.g. {.code {arg}[, keep_cols]})."),
      data = list(arg = arg, observed = bad_cols), call = call
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

# Coerce and validate a per-respondent completion-time vector (one total in seconds per
# respondent). cier_total_time rejects a 2-D input: which axis is respondents is ambiguous.
# Every observed time must be finite and strictly positive; NA is permitted and abstains.
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
  # is.na() is TRUE for NaN, so use is.nan() / is.infinite() on the raw vector to reject
  # NaN and +/-Inf before dropping genuine NA.
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
        "i" = "A zero total is almost always an export or merge artefact, not a \\
               real duration -- check the export. {.val NA} marks a genuinely \\
               missing time (it abstains)."),
      data = list(arg = arg), call = call
    )
  }
  as.numeric(seconds)
}

# Coerce and validate cier_page_time's n x pages matrix of per-page total times in seconds
# (data.frame / tibble coerced). Reuses check_responses() for the matrix / numeric / finite
# contract, then requires every observed cell non-negative (zero is valid speeding
# evidence; NA marks an untimed page).
check_page_seconds <- function(page_seconds, arg = "page_seconds",
                               call = rlang::caller_env()) {
  m <- check_responses(page_seconds, arg = arg, call = call)
  # check_responses already rejected NaN / Inf, so the only missing entries left are NA
  # (untimed pages), which na.rm skips. An all-NA matrix passes here and abstains.
  if (any(m < 0, na.rm = TRUE)) {
    cier_abort(
      "cier_error_input",
      c("{.arg {arg}} must be page times in seconds, zero or greater.",
        "x" = "Found negative values; a page time cannot be negative.",
        "i" = "A zero page time is kept as maximal speeding evidence (a per-item \\
               mean of 0 s is below {.arg min_seconds}); recode untimed pages to \\
               {.val NA}, which abstain."),
      data = list(arg = arg), call = call
    )
  }
  m
}

# Validate cier_page_time's per-page item-count vector: one positive whole number per
# column of `page_seconds` (the page total is divided by it for the mean per-item time).
# Returns the vector as integer.
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

# Validate cier_attention's per-check pass-set spec: a positional list, one element per
# column of `checks` (pass[[j]] applies to column j). Each pass-set must be a non-empty
# vector of finite passing values (an empty set would fail everyone). A named pass is
# cross-checked against check_names (see check_pass_alignment).
check_pass <- function(pass, n_checks, check_names = NULL, arg = "pass",
                       call = rlang::caller_env()) {
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
  check_pass_alignment(pass, check_names, arg, call)
  pass
}

# Cross-check a named `pass` against the columns of `checks`. `pass` binds by position
# (pass[[j]] -> column j), not by name, so names in a different order would silently score
# the wrong column. Engages only when a `pass` name matches a column; then the names must
# equal the column names position by position, else abort (no auto-reorder). Decorative
# labels referencing no column are ignored.
check_pass_alignment <- function(pass, check_names, arg, call) {
  pn <- names(pass)
  references_columns <- !is.null(pn) && !is.null(check_names) &&
    any(nzchar(pn)) && any(pn %in% check_names)
  if (!references_columns) {
    return(invisible(NULL))
  }
  mism <- which(is.na(pn) | is.na(check_names) | pn != check_names)
  if (length(mism) > 0L) {
    cier_abort(
      "cier_error_input",
      c("{.arg {arg}} names do not match the columns of {.arg checks}.",
        "x" = "Mismatch at {cli::qty(length(mism))}position{?s} {.val {mism}}: \\
               {.arg {arg}} {.val {pn[mism]}} vs {.arg checks} \\
               {.val {check_names[mism]}}.",
        "i" = "{.arg {arg}} binds to {.arg checks} by position, not by name. \\
               Order {.arg {arg}} to match the columns (same order), or drop \\
               the names."),
      data = list(arg = arg, observed = mism), call = call
    )
  }
  invisible(NULL)
}

# Require integer-coded responses: every non-NA value must be whole. cier_lazr scores a
# transition matrix over discrete anchors, so a continuous / averaged score has no Markov
# chain; a typed error beats silently dropping non-integer transitions. Run after
# check_responses, so non-NA values are already finite and the test reduces to wholeness.
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
