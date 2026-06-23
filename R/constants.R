# Controlled vocabularies for cier internals. Each accessor returns a fresh
# character vector; element [1] is the documented default where one exists.

# Cutoff resolution strategies. "mc_null" (nonparametric Monte-Carlo null for
# person-fit indices) needs the scored response block, so it resolves at the
# bridge rather than from a value alone.
cier_cutoff_methods <- function() {
  c("percentile", "fixed", "chisq", "mc_null")
}

# Which tail of a statistic flags carelessness.
cier_flag_directions <- function() {
  c("upper", "lower")
}

# Cutoff provenance recorded on each cier_index; finer-grained than the method-spec
# `default_cutoff_method`: splits "fixed" into fraction / literal-count and adds the
# override-only resolvers (median-relative, Kneedle, literal). `NA` = no provenance.
cier_index_cutoff_provenance <- function() {
  c("percentile", "fixed_fraction", "fixed_count", "chisq", "mc_null",
    "median_relative", "kneedle", "literal")
}

# Latent trait distributions cier_simulate draws from.
sim_trait_distributions <- function() {
  c("normal", "skew_normal", "t", "bimodal")
}
