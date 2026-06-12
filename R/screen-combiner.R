# Purpose: The internal orchestration + vote-combiner helpers behind
#          cier_screen(). cier_screen() runs the registry's screenable indices
#          over one dataset, skips (with a recorded reason) the ones whose inputs
#          are absent, then COLLAPSES correlated indices to one vote before
#          counting cross-index agreement. The collapse is the only new logic in
#          the slice: indices sharing a registry `vote_group` fuse into a single
#          vote that fires when ANY member flagged the respondent (logical OR),
#          an abstaining (NA) member contributing FALSE. The screen never alters
#          a per-index result -- each `cier_index` it returns is byte-identical to
#          calling that index directly.
# Args:    See per-function documentation.
# Returns: helpers returning the run set, skip reasons, the collapsed votes, the
#          null-rate vector, and the assembled cier_screen.
# Invariants:
#   - The run set and every output column are in registry order.
#   - `votes` carries no NA (abstain collapses to FALSE); `flags` keeps the NA.
#   - The agreement runs on the COLLAPSED votes, never the raw per-index flags.

# Methods that need the optional `items` metadata frame; the screen skips these
# structurally when `items` is NULL, and supplies `items` as a base argument when
# present. The wrapper remains the single source of truth for WHICH columns each
# needs -- this set only decides screen eligibility. It is the SINGLE place that
# encodes "needs items" (both call sites below read it), so adding an items index
# means editing one function. A registry `requires_items` column would data-drive
# it, but the v0 battery is a frozen ten (a new index needs an ADR entry anyway),
# so a registry column is deferred rather than added speculatively.
screen_items_methods <- function() {
  c("cier_even_odd", "cier_personal_reliability", "cier_gnormed", "cier_ht")
}

# Resolve the run set: the screenable registry methods intersected with the
# user's `methods` request (NULL = all), in registry order. An unknown id aborts.
screen_resolve_methods <- function(methods, call) {
  reg <- load_method_registry()
  screenable <- reg$method[reg$screenable]
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
    cier_abort(
      "cier_error_input",
      c("Unknown screenable method{?s}: {.val {bad}}.",
        "i" = "See {.fun cier_methods} for the available set."),
      data = list(arg = "methods", observed = bad), call = call
    )
  }
  screenable[screenable %in% methods]            # keep registry order
}

# Validate the per-index control list: a named list whose names are among the
# selected methods and whose elements are themselves lists of wrapper arguments
# (spliced into the index call). The names are checked against the SELECTED set
# (before skipping), so an override for an index that later skips is allowed.
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
  # nzchar(NA) is TRUE, so an NA name must be caught explicitly here -- otherwise
  # it slips through to the "not selected" check with a confusing message.
  if (is.null(nms) || any(is.na(nms)) || any(!nzchar(nms))) {
    cier_abort(
      "cier_error_input",
      "Every {.arg control} entry must be named by its method id.",
      data = list(arg = "control"), call = call
    )
  }
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

# Validate one control entry's inner arguments against the target index's
# formals: every element must be named, and every name must be a real argument of
# that index (excluding the screen-supplied `responses` / `items` positionals). A
# typo (`fp` for `fpr`) or a collision with a managed positional (`responses` /
# `items`) is then a typed input error, not an opaque `do.call` failure or a
# silent positional bind.
check_control_args <- function(method, args, call) {
  if (length(args) == 0L) {
    return(invisible(NULL))
  }
  nms <- names(args)
  if (is.null(nms) || any(is.na(nms)) || any(!nzchar(nms))) {
    cier_abort(
      "cier_error_input",
      "Every {.arg control} argument for {.val {method}} must be named.",
      data = list(arg = "control", method = method), call = call
    )
  }
  allowed <- setdiff(names(formals(match.fun(method))), c("responses", "items"))
  bad <- setdiff(nms, allowed)
  if (length(bad) > 0L) {
    cier_abort(
      "cier_error_input",
      # qty(bad) pins the plural quantity: the message interpolates two vectors
      # ({method} and {bad}), so {?s} cannot infer which one drives it.
      c("Unknown {cli::qty(bad)}{.arg control} argument{?s} for \\
         {.val {method}}: {.val {bad}}.",
        "i" = "Valid arguments: {.val {allowed}}."),
      data = list(arg = "control", method = method, observed = bad), call = call
    )
  }
  invisible(NULL)
}

# The structural skip reason for one method, or NA when it is eligible to run.
# Missing `items` is checked first (it gates the metadata indices); then an
# absent Suggests backend (PerFit / mokken). A genuinely malformed `items` is
# NOT a skip -- the index's own typed error is left to propagate.
screen_skip_reason <- function(method, items, reg_row) {
  if (method %in% screen_items_methods() && is.null(items)) {
    return("needs item metadata (items = NULL)")
  }
  backend <- reg_row$backend
  if (!is.na(backend) && !cier_namespace_present(backend)) {
    return(paste0(backend, " not installed"))
  }
  NA_character_
}

# Call one index wrapper with the base arguments it needs plus the user's
# per-index control overrides. Matrix-only indices take just `responses`; the
# metadata indices also take `items`. Returns the `cier_index`, or a character
# skip REASON when the index aborted on a typed backend limit (see below).
screen_call_index <- function(method, responses, items, args) {
  base <- if (method %in% screen_items_methods()) {
    list(responses, items)
  } else {
    list(responses)
  }
  # The screen owns the abstention narrative. A fully-abstaining index (one that
  # scored no respondent -- e.g. psychsyn / psychant find no qualifying pairs on a
  # broad inventory) raises `cier_warning_insufficient_items` from its standalone
  # cutoff resolution, but the screen already reports that case transparently in
  # `$flags` / `print` as "0 / 0", so re-emitting it as a warning is redundant
  # noise that would fire on every broad-inventory screen. The same holds for
  # `cier_warning_saturated_cutoff` (a percentile cutoff sitting on a tie mass at
  # the score extreme, so the realised rate exceeds fpr): the screen already
  # prints the per-index flag rate, so the saturation note is redundant here --
  # and even-odd's point mass at +1 would otherwise fire it on essentially every
  # Likert screen. Muffle ONLY those two typed conditions (a targeted handler,
  # NOT a blanket suppressWarnings -- e.g. a singular-covariance warning still
  # propagates), and only around the screen's own call; a direct `cier_<index>()`
  # call still warns, which is where the ranking-convention caveat is wanted.
  #
  # A typed BACKEND LIMIT on otherwise-valid data (cier_error_backend_limit, e.g.
  # mokken's 10-category ceiling for cier_ht) must not abort the battery: catch
  # exactly that subclass and return its compact reason so the caller records the
  # index as skipped. Every other error -- a malformed `items` frame included --
  # still propagates so the user can fix the input.
  tryCatch(
    withCallingHandlers(
      do.call(method, c(base, args)),
      cier_warning_insufficient_items = function(w) invokeRestart("muffleWarning"),
      cier_warning_saturated_cutoff = function(w) invokeRestart("muffleWarning")
    ),
    cier_error_backend_limit = function(e) {
      reason <- cier_condition_data(e)$reason
      if (is.null(reason)) "backend limit" else reason
    }
  )
}

# A data.frame with `n` rows and no columns: the empty flag / vote table when no
# index ran. It preserves the respondent count so a downstream `rowSums()` still
# yields one zero per respondent (and `ncol == 0` signals "no agreement").
empty_n_row_df <- function(n) {
  as.data.frame(matrix(nrow = n, ncol = 0L))
}

# Collapse a per-index flag table to per-vote-group votes by OR, with NA ->
# FALSE. Group order is the first-appearance order of `vote_group` across the
# flag columns (registry order). Mirrors the independent oracle
# ref_collapse_votes() in tests/reference/ (kept in lock-step deliberately).
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
      acc <- acc | f                             # logical OR across members
    }
    acc
  })
  names(votes) <- groups
  as.data.frame(votes, stringsAsFactors = FALSE, check.names = FALSE)
}

# The calibrated null nominal per vote group, for flag_agreement()'s per-vote
# excess. A singleton vote whose index is null-referenced carries that index's
# nominal (Mahalanobis chi-square `alpha`; Gnormed's Monte-Carlo `fpr`),
# honouring a control override; every percentile vote -- including the
# all-percentile consistency collapse -- is tautological, so its null is NA.
# `vote_group` is the ran-methods map, so `names(vote_group)[vote_group == g]`
# recovers exactly the ran members of group `g` (the same set collapse_votes
# OR-fused). A collapsed multi-member group (length != 1) has no single nominal
# rate, so its null is NA -- which also means a (hypothetical) null-referenced
# index sharing a group would not surface its nominal; today the two
# null-referenced indices are both singleton groups, so this never bites.
screen_null_rate <- function(groups, vote_group, control, reg) {
  vapply(groups, function(g) {
    members <- names(vote_group)[vote_group == g]
    if (length(members) != 1L) {
      return(NA_real_)
    }
    row <- reg[reg$method == members, , drop = FALSE]
    args <- control[[members]]
    # A literal `cutoff` override unhooks the vote from any nominal rate (it flags
    # on an absolute threshold), so it is no longer referenced to a calibrated
    # null -- report NA, like the tautological percentile votes.
    if (!is.null(args$cutoff)) {
      NA_real_
    } else if (identical(row$default_cutoff_method, "chisq")) {
      if (is.null(args$alpha)) row$default_cutoff_value else args$alpha
    } else if (identical(row$default_cutoff_method, "perfit_null")) {
      if (is.null(args$fpr)) row$default_cutoff_value else args$fpr
    } else {
      NA_real_
    }
  }, numeric(1L))
}

# Assemble the cier_screen from the indices that ran and the recorded skips:
# build the per-index flag table, collapse it to votes, run the agreement
# diagnostic on the votes (NULL when nothing ran), and bind the skip rows.
build_cier_screen <- function(indices, skipped_rows, reg, control, n) {
  ran <- names(indices)
  vote_group <- stats::setNames(reg$vote_group[match(ran, reg$method)], ran)
  flags <- if (length(indices) > 0L) {
    as.data.frame(lapply(indices, function(ix) ix$flagged), check.names = FALSE)
  } else {
    empty_n_row_df(n)
  }
  votes <- collapse_votes(flags, vote_group)
  agreement <- if (ncol(votes) > 0L) {
    null_rate <- screen_null_rate(colnames(votes), vote_group, control, reg)
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
  new_cier_screen(indices, flags, vote_group, votes, agreement, skipped, n)
}
