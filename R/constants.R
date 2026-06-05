# Purpose: Controlled vocabularies for cier, trimmed to what v0 reads. One
#          source of truth so the registry loader/validator pin against the
#          same level sets.
# Args:    None (accessors).
# Returns: Each accessor returns a fresh character vector; element [1] is the
#          documented default where the concept has one.
# Invariants:
#   - Level sets are class-stable; the registry validator relies on identity.
#   - Accessors are pure (return a copy each call, no side effects).

# Method families spanned by the v0 indices.
cier_family_levels <- function() {
  c("indirect", "personfit")
}

# Cutoff resolution strategies used by the registry.
cier_cutoff_methods <- function() {
  c("percentile", "fixed", "chisq")
}

# Which tail of a statistic flags carelessness.
cier_flag_directions <- function() {
  c("upper", "lower")
}
