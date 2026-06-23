# Internal orchestration + vote-combiner helpers behind cier_screen(). Indices sharing
# a method-spec `vote_group` collapse into one vote firing when any member flagged the
# respondent (OR, NA -> FALSE). `flags` keeps NA, `votes` does not; agreement runs on
# the collapsed votes, never the raw flags.

# Abort on requested method ids the screen cannot run, distinguishing a
# registered-but-not-screenable standalone index (call it directly) from an unknown id.
# `all_methods` is the full internal method-spec vector.
screen_unknown_methods_abort <- function(bad, all_methods, screenable, call) {
  not_screenable <- intersect(bad, all_methods)
  unknown <- setdiff(bad, all_methods)
  msg <- character(0L)
  if (length(not_screenable) > 0L) {
    msg <- c(
      msg,
      "{cli::qty(not_screenable)}Method{?s} not run by {.fun cier_screen}: \\
       {.val {not_screenable}}.",
      "i" = "These are standalone indices -- call them directly, e.g. \\
             {.code {not_screenable[[1L]]}(responses)}."
    )
  }
  if (length(unknown) > 0L) {
    unknown_lines <- c(
      "{cli::qty(unknown)}Unknown method{?s}: {.val {unknown}}.",
      "i" = "Screenable methods: {.val {screenable}}."
    )
    # First line is the abort header when alone, a bullet otherwise.
    names(unknown_lines)[[1L]] <- if (length(not_screenable) > 0L) "x" else ""
    msg <- c(msg, unknown_lines)
  }
  cier_abort("cier_error_input", msg,
             data = list(arg = "methods", observed = bad), call = call)
}

# Resolve the run set: screenable method specs intersected with `methods` (NULL = all),
# in method-spec order. An unrunnable id aborts (see screen_unknown_methods_abort).
screen_resolve_methods <- function(methods, call) {
  screenable <- cier_screenable_methods()
  if (is.null(methods)) {
    return(screenable)
  }
  if (!is.character(methods) || anyNA(methods) || length(methods) == 0L) {
    cier_abort(
      "cier_error_input",
      "{.arg methods} must be {.val NULL} or a character vector of method ids.",
      data = list(arg = "methods"), call = call
    )
  }
  bad <- setdiff(methods, screenable)
  if (length(bad) > 0L) {
    screen_unknown_methods_abort(bad, cier_method_specs()$method, screenable, call)
  }
  screenable[screenable %in% methods]            # method-spec order
}

# Validate the per-index control list: a named list keyed by selected methods, each
# element an argument list. Names are checked against the selected set before skipping,
# so an override for an index that later skips is allowed.
screen_check_control <- function(control, selected, call) {
  if (!is.list(control)) {
    cier_abort(
      "cier_error_input",
      "{.arg control} must be a named list of per-method argument lists.",
      data = list(arg = "control"), call = call
    )
  }
  if (length(control) == 0L) {
    return(invisible(NULL))
  }
  nms <- names(control)
  check_unique_names(nms, "control", call = call)
  if (!all(vapply(control, is.list, logical(1L)))) {
    cier_abort(
      "cier_error_input",
      "Each {.arg control} entry must be a list of arguments for one index.",
      data = list(arg = "control"), call = call
    )
  }
  bad <- setdiff(nms, selected)
  if (length(bad) > 0L) {
    cier_abort(
      "cier_error_input",
      c("{.arg control} names method{?s} not selected to run: {.val {bad}}.",
        "i" = "Add it to {.arg methods}, or remove the override."),
      data = list(arg = "control", observed = bad), call = call
    )
  }
  for (nm in nms) {
    check_control_args(nm, control[[nm]], call)
  }
  invisible(NULL)
}

# Validate one control entry's arguments against the target index's formals: every
# element named, every name a real argument (excluding the screen-supplied `responses`
# / `items`). Turns a typo into a typed error rather than an opaque do.call failure.
check_control_args <- function(method, args, call) {
  if (length(args) == 0L) {
    return(invisible(NULL))
  }
  nms <- names(args)
  check_unique_names(nms, "control", call = call)
  allowed <- setdiff(names(formals(match.fun(method))), c("responses", "items"))
  bad <- setdiff(nms, allowed)
  if (length(bad) > 0L) {
    cier_abort(
      "cier_error_input",
      # qty(bad) pins the plural quantity: two vectors interpolate, so {?s} cannot
      # infer which drives it.
      c("Unknown {cli::qty(bad)}{.arg control} argument{?s} for \\
         {.val {method}}: {.val {bad}}.",
        "i" = "Valid arguments: {.val {allowed}}."),
      data = list(arg = "control", method = method, observed = bad), call = call
    )
  }
  invisible(NULL)
}

# The structural skip reason for one method, or NA when eligible to run. Missing `items`
# gates metadata indices; a malformed `items` is NOT a skip (the index's error propagates).
screen_skip_reason <- function(method, items) {
  if (method %in% cier_methods_requiring_items() && is.null(items)) {
    return("needs item metadata (items = NULL)")
  }
  NA_character_
}

# One-line note for a "0 / 0" pair-index line, from the captured cier_warning_no_pairs
# payload (`data`; the pairing kind comes from the payload, not the method id).
screen_no_pairs_note <- function(data) {
  antonym <- isTRUE(data$antonym)
  noun  <- if (antonym) "antonym" else "synonym"
  sweep <- if (antonym) {
    "cier_psychsyn_critval(antonym = TRUE)"
  } else {
    "cier_psychsyn_critval()"
  }
  sprintf(
    "no %s pairs at critical_r = %s (strongest in-tail r = %s); see %s",
    noun, format(data$critical_r), format(round(data$strongest_r, 3)), sweep
  )
}

# Forward a battery-wide `fpr` to percentile-default indices only; Mahalanobis `alpha`,
# longstring `frac`, and Gnormed nominal stay on their own scales. A per-index `control`
# `fpr` (or a literal `cutoff`, mutually exclusive in the wrapper) wins. Validated
# here, before the wrapper sees it.
screen_apply_fpr <- function(control, selected, fpr, call) {
  if (is.null(fpr)) {
    return(control)
  }
  check_open_unit(fpr, "fpr", call = call)
  for (m in intersect(selected, cier_percentile_methods())) {
    args <- control[[m]]
    if (is.null(args$fpr) && is.null(args$cutoff)) {
      args$fpr <- fpr
      control[[m]] <- args
    }
  }
  control
}

# Call one index wrapper with its base arguments plus the user's control overrides.
# Returns list(result, note): `result` is the cier_index, or a character skip reason
# when the index aborted on a typed backend limit; `note` is a one-line no-pairs
# annotation (or NULL) for print's "0 / 0" line.
screen_call_index <- function(method, responses, items, args) {
  base <- if (method %in% cier_methods_requiring_items()) {
    list(responses, items)
  } else {
    list(responses)
  }
  # Environment accumulator captures the no-pairs payload from the handler without `<<-`.
  cap <- new.env(parent = emptyenv())
  cap$note <- NULL
  # The screen owns the abstention narrative. Muffle ONLY the two abstention warnings
  # (no qualifying pairs -> insufficient_items; percentile cutoff on a tie mass ->
  # saturated_cutoff), both already visible in the printed "0 / 0" / flag rate. Targeted,
  # not a blanket suppressWarnings (a singular-covariance warning still propagates), and
  # only around the screen's call; a direct index call still warns.
  #
  # A typed backend limit (cier_error_backend_limit, e.g. Gnormed's single-Ncat contract
  # on mixed formats) must not abort the battery: catch it and return its reason so the
  # caller records a skip. Every other error -- malformed `items` included -- propagates.
  result <- tryCatch(
    withCallingHandlers(
      do.call(method, c(base, args)),
      # Capture the no-pairs abstention for the "0 / 0" annotation, then muffle.
      cier_warning_no_pairs = function(w) {
        cap$note <- screen_no_pairs_note(cier_condition_data(w))
        invokeRestart("muffleWarning")
      },
      cier_warning_insufficient_items = function(w) invokeRestart("muffleWarning"),
      cier_warning_saturated_cutoff = function(w) invokeRestart("muffleWarning"),
      # The missing-reverse_keyed message fires once at the screen level; muffle the
      # per-index repeats.
      cier_message_forward_keyed = function(m) invokeRestart("muffleMessage")
    ),
    cier_error_backend_limit = function(e) {
      reason <- cier_condition_data(e)$reason
      if (is.null(reason)) "backend limit" else reason
    }
  )
  list(result = result, note = cap$note)
}

# A data.frame with `n` rows and no columns: the empty flag / vote table when no index
# ran. Keeping the row count lets `rowSums()` yield one zero per respondent, and
# `ncol == 0` signals "no agreement".
empty_n_row_df <- function(n) {
  as.data.frame(matrix(nrow = n, ncol = 0L))
}

# Collapse a per-index flag table to per-vote-group votes by OR, with NA -> FALSE.
# Group order is the first-appearance order of `vote_group` across the flag columns.
collapse_votes <- function(flags, vote_group) {
  methods <- names(flags)
  if (length(methods) == 0L) {
    return(empty_n_row_df(nrow(flags)))
  }
  groups <- unique(unname(vote_group[methods]))
  votes <- lapply(groups, function(g) {
    members <- methods[vote_group[methods] == g]
    acc <- rep(FALSE, nrow(flags))
    for (m in members) {
      f <- flags[[m]]
      f[is.na(f)] <- FALSE                       # abstain -> did not flag
      acc <- acc | f
    }
    acc
  })
  names(votes) <- groups
  as.data.frame(votes, stringsAsFactors = FALSE, check.names = FALSE)
}

# The calibrated null nominal per vote group, for flag_agreement()'s per-vote excess.
# A singleton null-referenced vote carries that index's nominal (Mahalanobis `alpha`;
# Monte-Carlo `fpr` of Gnormed, or Ht when opted into "mc_null"), honouring a control
# override of method, rate, or a literal `cutoff`. Percentile votes are tautological ->
# NA. Null-referenced indices are singletons, so a multi-member group's NA never bites.
screen_null_rate <- function(groups, vote_group, control) {
  vapply(groups, function(g) {
    members <- names(vote_group)[vote_group == g]
    if (length(members) != 1L) {
      return(NA_real_)
    }
    row <- cier_method_spec(members)
    args <- control[[members]]
    # Honour a control `method` override (Ht's opt-in "mc_null"), not just the static
    # default -- else a control-calibrated Ht vote is mis-reported as percentile.
    method <- args$method %||% row$default_cutoff_method
    # A literal `cutoff` flags on an absolute threshold, unhooked from any calibrated
    # null -> NA, like the percentile votes.
    if (!is.null(args$cutoff)) {
      NA_real_
    } else if (identical(method, "chisq")) {
      if (is.null(args$alpha)) row$default_cutoff_value else args$alpha
    } else if (identical(method, "mc_null")) {
      if (is.null(args$fpr)) row$default_cutoff_value else args$fpr
    } else {
      NA_real_
    }
  }, numeric(1L))
}

# Assemble the cier_screen: build the per-index flag table, collapse to votes, run the
# agreement diagnostic on them (NULL when nothing ran), bind the skip + note rows.
build_cier_screen <- function(indices, skipped_rows, notes_rows, control, n) {
  ran <- names(indices)
  specs <- cier_method_specs()
  vote_group <- stats::setNames(specs$vote_group[match(ran, specs$method)], ran)
  flags <- if (length(indices) > 0L) {
    as.data.frame(lapply(indices, function(ix) ix$flagged), check.names = FALSE)
  } else {
    empty_n_row_df(n)
  }
  votes <- collapse_votes(flags, vote_group)
  agreement <- if (ncol(votes) > 0L) {
    null_rate <- screen_null_rate(colnames(votes), vote_group, control)
    flag_agreement(as.matrix(votes), null_rate = unname(null_rate))
  } else {
    NULL
  }
  skipped <- if (length(skipped_rows) > 0L) {
    do.call(rbind, skipped_rows)
  } else {
    data.frame(method = character(0L), reason = character(0L),
               stringsAsFactors = FALSE)
  }
  notes <- if (length(notes_rows) > 0L) {
    do.call(rbind, notes_rows)
  } else {
    data.frame(method = character(0L), note = character(0L),
               stringsAsFactors = FALSE)
  }
  new_cier_screen(indices, flags, vote_group, votes, agreement, skipped, notes, n)
}
