# Shared fixtures for the split-half consistency index tests (cier_even_odd,
# cier_personal_reliability). Both score within-scale item halves, so they share the
# scale-blocked items builder, an independent block-index helper (the oracles need it
# but must not borrow scale_block_indices()), and a reverse-scorer. rand_matrix()
# stays local to each test file -- that name is reused with different defaults across
# other index tests, so a global definition would mislead.

# Scale-blocked `items`: `n_scales` scales of `per_scale` items. `reverse_keyed`
# defaults to alternating; pass FALSE for analytic / parity fixtures that isolate the
# correlation math from reverse-scoring.
blocked_items <- function(n_scales = 3L, per_scale = 4L, max = 5L,
                          reverse_keyed = NULL) {
  scale <- rep(LETTERS[seq_len(n_scales)], each = per_scale)
  if (is.null(reverse_keyed)) {
    reverse_keyed <- rep(c(FALSE, TRUE), length.out = length(scale))
  }
  data.frame(scale = scale, reverse_keyed = reverse_keyed,
             max = max, stringsAsFactors = FALSE)
}

# Build the scale blocks INDEPENDENTLY of production (the oracles need them but must
# not borrow scale_block_indices()).
blocks_from_scale <- function(items) {
  uniq <- unique(items$scale)
  lapply(uniq, function(s) which(items$scale == s))
}

# Independently reverse-score (min + max) - x on reverse items; fixtures here are
# 1-based, so the reflection is (1 + max) - x.
prescore <- function(x, items) {
  rk <- items$reverse_keyed
  x[, rk] <- (items$max[rk] + 1L) - x[, rk]
  x
}

# Shared input-validation contract for the split-half consistency indices
# (cier_even_odd, cier_personal_reliability): both reject the same bad scale structure
# / items frame / payload and the same malformed fpr / cutoff overrides (a value is a
# [-1, 1] correlation), the rejection firing before any scoring (so resample / seed do
# not matter and need not be passed). Each test passes its own wrapper; the fixtures
# are built inline (rand_matrix is file-local) and any valid-shaped numeric matrix
# triggers the error paths under test.
expect_split_half_input_rejected <- function(fn) {
  it  <- blocked_items(3L, 4L, reverse_keyed = FALSE)             # 3 scales, 12 items
  x12 <- matrix(as.double(rep_len(1:5, 120L)), nrow = 10L, ncol = 12L)
  # Structural: fewer than two scales, missing / misaligned items, non-numeric payload.
  it_one <- data.frame(scale = rep("A", 4L), reverse_keyed = FALSE, max = 5L)
  x_two  <- matrix(c(1, 2, 3, 4, 5, 4, 3, 2), nrow = 2L, byrow = TRUE)
  storage.mode(x_two) <- "double"
  expect_error(fn(x_two, it_one), class = "cier_error_input")
  it_noscale <- it
  it_noscale$scale <- NULL
  expect_error(fn(x12, it_noscale), class = "cier_error_input")
  expect_error(fn(matrix(as.double(rep_len(1:5, 50L)), 5L, 10L), it),  # 10 cols vs 12
               class = "cier_error_input")
  expect_error(fn(1:10, it), class = "cier_error_input")
  expect_error(fn(matrix(letters[1:12], nrow = 1L), it), class = "cier_error_input")
  # Overrides: fpr outside (0, 1) / non-scalar / non-numeric; cutoff outside [-1, 1] /
  # non-scalar / non-numeric; and fpr + cutoff supplied together.
  expect_error(fn(x12, it, fpr = 0), class = "cier_error_input")
  expect_error(fn(x12, it, fpr = 1), class = "cier_error_input")
  expect_error(fn(x12, it, fpr = -0.1), class = "cier_error_input")
  expect_error(fn(x12, it, fpr = NA_real_), class = "cier_error_input")
  expect_error(fn(x12, it, fpr = c(0.05, 0.1)), class = "cier_error_input")
  expect_error(fn(x12, it, fpr = "x"), class = "cier_error_input")
  expect_error(fn(x12, it, cutoff = c(0.5, 1)), class = "cier_error_input")
  expect_error(fn(x12, it, cutoff = NA_real_), class = "cier_error_input")
  expect_error(fn(x12, it, cutoff = "x"), class = "cier_error_input")
  expect_error(fn(x12, it, cutoff = -1.5), class = "cier_error_input")
  expect_error(fn(x12, it, cutoff = 1.5), class = "cier_error_input")
  expect_error(fn(x12, it, fpr = 0.1, cutoff = 0.5), class = "cier_error_input")
  invisible(NULL)
}
