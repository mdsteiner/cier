# Shared fixtures for the person-fit index tests (cier_gnormed, cier_ht). Both
# score a complete-case polytomous block, so they share the same generators.

# Reproducible polytomous matrix in 1..ncat coding (the wrappers zero-base). Uniform
# sampling spreads the person-fit values (needed for the direction / percentile tests).
poly_matrix <- function(n = 60L, p = 12L, ncat = 5L, seed = 21L) {
  withr::with_seed(seed, {
    m <- matrix(sample.int(ncat, n * p, replace = TRUE), nrow = n, ncol = p)
  })
  storage.mode(m) <- "double"
  m
}

# Item metadata: homogeneous span (1-based fixtures, so `max = ncat`), optional
# `reverse_keyed`. Doubles as each oracle's `data$items` (it reads items$reverse_keyed
# and items$max under its own 1-based contract).
poly_items <- function(p = 12L, reverse = FALSE, ncat = 5L) {
  rk <- if (length(reverse) == 1L) rep(reverse, p) else reverse
  data.frame(reverse_keyed = rk, max = ncat)
}
