# Internal (not a public registry): per-method cutoff defaults, flag direction,
# screen membership, vote grouping, and whether the index requires item metadata.

.cier_method_specs <- data.frame(
  method = c(
    "cier_longstring",
    "cier_irv",
    "cier_even_odd",
    "cier_personal_reliability",
    "cier_psychsyn",
    "cier_psychant",
    "cier_mahalanobis",
    "cier_person_total",
    "cier_gnormed",
    "cier_ht",
    "cier_autocorrelation",
    "cier_lazr",
    "cier_total_time",
    "cier_page_time",
    "cier_attention"
  ),
  default_cutoff_method = c(
    "fixed",
    "percentile",
    "percentile",
    "percentile",
    "percentile",
    "percentile",
    "chisq",
    "percentile",
    "mc_null",
    "percentile",
    "percentile",
    "percentile",
    "percentile",
    "fixed",
    "fixed"
  ),
  default_cutoff_value = c(
    0.5,
    0.05,
    0.05,
    0.05,
    0.05,
    0.05,
    0.001,
    0.05,
    0.05,
    0.05,
    0.05,
    0.05,
    0.05,
    1,
    1
  ),
  flag_direction = c(
    "upper",
    "lower",
    "upper",
    "upper",
    "lower",
    "upper",
    "upper",
    "lower",
    "upper",
    "lower",
    "upper",
    "upper",
    "lower",
    "upper",
    "upper"
  ),
  screenable = c(
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE
  ),
  vote_group = c(
    "cier_longstring",
    "cier_irv",
    "consistency",
    "consistency",
    "cier_psychsyn",
    "cier_psychant",
    "cier_mahalanobis",
    "cier_person_total",
    "cier_gnormed",
    "cier_ht",
    "cier_autocorrelation",
    "cier_lazr",
    "cier_total_time",
    "cier_page_time",
    "cier_attention"
  ),
  requires_items = c(
    FALSE,
    FALSE,
    TRUE,
    TRUE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    TRUE,
    TRUE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE
  ),
  stringsAsFactors = FALSE
)

cier_method_specs <- function() {
  validate_cier_method_specs(.cier_method_specs)
  .cier_method_specs
}

validate_cier_method_specs <- function(x) {
  validate_cier_method_spec_schema(x)
  validate_cier_method_spec_values(x)
  invisible(x)
}

validate_cier_method_spec_schema <- function(x) {
  required <- c("method", "default_cutoff_method", "default_cutoff_value",
                "flag_direction", "screenable", "vote_group", "requires_items")
  if (!identical(names(x), required)) {
    cier_abort(
      "cier_error_state",
      "{.cls cier_method_specs} columns do not match the internal schema."
    )
  }
  if (anyNA(x$method) || any(!nzchar(x$method)) ||
        anyDuplicated(x$method) > 0L) {
    cier_abort("cier_error_state", "Method ids must be unique and non-missing.")
  }
  invisible(x)
}

validate_cier_method_spec_values <- function(x) {
  if (anyNA(x$default_cutoff_method) ||
        length(setdiff(x$default_cutoff_method, cier_cutoff_methods())) > 0L) {
    cier_abort("cier_error_state", "Unknown method cutoff strategy.")
  }
  if (anyNA(x$flag_direction) ||
        length(setdiff(x$flag_direction, cier_flag_directions())) > 0L) {
    cier_abort("cier_error_state", "Unknown method flag direction.")
  }
  if (anyNA(x$default_cutoff_value) || any(x$default_cutoff_value <= 0)) {
    cier_abort("cier_error_state", "Method cutoff defaults must be positive.")
  }
  prob <- x$default_cutoff_method %in% c("percentile", "chisq", "mc_null")
  if (any(x$default_cutoff_value[prob] >= 1)) {
    cier_abort(
      "cier_error_state",
      "Probability cutoff defaults must be inside the open unit interval."
    )
  }
  if (anyNA(x$screenable) || anyNA(x$requires_items)) {
    cier_abort("cier_error_state", "Method flags must be non-missing logicals.")
  }
  if (anyNA(x$vote_group) || any(!nzchar(x$vote_group))) {
    cier_abort("cier_error_state", "Method vote groups must be non-missing.")
  }
  invisible(x)
}

cier_method_spec <- function(method, call = rlang::caller_env()) {
  check_string(method, "method", call = call)
  specs <- cier_method_specs()
  idx <- match(method, specs$method)
  if (is.na(idx)) {
    cier_abort(
      "cier_error_input",
      c("Unknown method {.val {method}}.",
        "i" = "Available methods: {.val {specs$method}}."),
      data = list(arg = "method", observed = method),
      call = call
    )
  }
  specs[idx, , drop = FALSE]
}

cier_screenable_methods <- function() {
  specs <- cier_method_specs()
  specs$method[specs$screenable]
}

cier_standalone_methods <- function() {
  specs <- cier_method_specs()
  specs$method[!specs$screenable]
}

cier_methods_requiring_items <- function() {
  specs <- cier_method_specs()
  specs$method[specs$requires_items]
}

cier_percentile_methods <- function() {
  specs <- cier_method_specs()
  specs$method[specs$default_cutoff_method == "percentile"]
}
