# Purpose: Controlled vocabularies for cier, trimmed to what v0 reads. One
#          source of truth so the registry loader/validator pin against the
#          same level sets.
# Args:    None (accessors).
# Returns: Each accessor returns a fresh character vector; element [1] is the
#          documented default where the concept has one.
# Invariants:
#   - Level sets are class-stable; the registry validator relies on identity.
#   - Accessors are pure (return a copy each call, no side effects).

# Method families spanned by the indices. v0.2 adds the timing family (total /
# page time) and the direct family (attention checks); autocorrelation and Laz.R
# remain indirect (response-pattern) indices.
cier_family_levels <- function() {
  c("indirect", "personfit", "timing", "direct")
}

# Cutoff resolution strategies used by the registry. "perfit_null" is the
# PerFit Monte-Carlo null for the nonparametric person-fit indices (Gnormed):
# unlike the value-only strategies it is referenced to a simulated null and so is
# resolved at the bridge (it needs the fitted object), not in resolve_cutoff().
cier_cutoff_methods <- function() {
  c("percentile", "fixed", "chisq", "perfit_null")
}

# Which tail of a statistic flags carelessness.
cier_flag_directions <- function() {
  c("upper", "lower")
}
