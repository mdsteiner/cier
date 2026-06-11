# Independent reference for cier_total_time() (timing family; Ward & Meade 2023,
# <doi:10.1146/annurev-psych-040422-045007>; Huang et al. 2012).
#
# total_time is not a statistic computed from a response matrix -- the input IS
# the per-respondent total completion time in seconds, so the per-respondent
# value is the validated seconds vector itself (an identity). What this slice
# pins is therefore the CONTRACT, not a derivation: that the wrapper applies no
# transform to the input, and that its three cutoff resolutions (empirical
# lower percentile, median-relative override, literal threshold) and the
# lower-tail flag rule are exactly the base-R quantities below.
#
# These helpers re-derive each quantity independently and NEVER call the
# production wrapper / resolvers, so any divergence is attributable to the
# wrapper. total_time has no CRAN parity partner (no package implements
# completion time as a C/IER index; verified 2026-06-10), so this oracle is the
# sole parity check -- oracle-only trust, like PR / RPR. Same base-R ops as the
# wrapper, so the target tolerance is 0 (exact); see tests/reference/TOLERANCES.md.

# Identity: the per-respondent value is the validated seconds vector, NA
# preserved (an all-NA respondent abstains). A mutant that sums, scales, ranks,
# or log-transforms the input diverges from this.
ref_total_time <- function(seconds) {
  as.numeric(seconds)
}

# Default cutoff: the empirical lower-tail percentile of the finite values at
# the target false-positive mass `fpr`. Lower direction takes the `fpr` quantile
# directly (NO flip to 1 - fpr). type = 7 is R's quantile default.
ref_total_time_percentile_cutoff <- function(value, fpr = 0.05) {
  finite <- value[is.finite(value)]
  if (length(finite) == 0L) {
    return(NA_real_)
  }
  as.numeric(stats::quantile(finite, probs = fpr, names = FALSE, type = 7L))
}

# Median-relative override: flag respondents faster than `frac` of the sample
# median (Leiner 2019 RSI; Greszki et al. 2015). The cutoff is `frac * median`
# of the finite values; an empty finite set yields NA (abstain).
ref_total_time_median_cutoff <- function(value, frac) {
  finite <- value[is.finite(value)]
  if (length(finite) == 0L) {
    return(NA_real_)
  }
  frac * stats::median(finite)
}

# Lower-tail flag rule: a respondent is flagged iff their (non-NA) value is at or
# below the cutoff (ties flag). An NA cutoff (abstaining resolver) flags nobody.
ref_total_time_flags <- function(value, cutoff) {
  if (is.na(cutoff)) {
    return(rep(FALSE, length(value)))
  }
  !is.na(value) & value <= cutoff
}
