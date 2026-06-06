# Purpose: Low-level numerical kernels shared by the indirect-index public
#          wrappers (R/cier-longstring.R, ...). Each kernel takes a numeric
#          matrix of responses and returns the per-respondent scores it is
#          responsible for. Single-kernel rule: one production implementation
#          per statistic.
# Args:    See per-kernel documentation below.
# Returns: Numeric vectors; never raises typed errors (the wrappers validate).
# Invariants:
#   - Kernels are pure (no I/O, no global state) and never mutate inputs.

# ---- Longstring -------------------------------------------------------------

# Maximum run length of consecutive identical responses per respondent, over
# the *raw* row (no scale blocking). Bytewise compatible with
# careless::longstring(): base::rle() treats each NA as its own run (NA == NA
# is NA, not TRUE), so identical values separated by NA do not merge and an
# all-NA row yields a max run length of 1. The wrapper applies NA-abstention for
# rows with no present responses; the kernel itself stays pure.
kernel_longstring <- function(responses) {
  vapply(
    seq_len(nrow(responses)),
    function(i) max(rle(responses[i, ])$lengths),
    numeric(1L)
  )
}
