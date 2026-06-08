# Paper-faithful reference implementation of the psychometric-synonyms C/IER
# index per Meade & Craig (2012, *Psychological Methods*). Independent of the
# production kernel: it re-derives the statistic from scratch with explicit
# nested loops over the inter-item correlation matrix and a per-respondent
# stacked-pair correlation, and never calls cier's kernel.
#
# Definition:
#   1. From the inter-item correlation matrix R (pairwise-complete), identify the
#      set of item pairs P = {(i, j): R_{ij} > critval, i > j}.
#   2. For each respondent, stack the values of the two columns of each pair into
#      two vectors (the larger item index first, the smaller second). Compute the
#      Pearson correlation across those stacked pair values.
#   3. Return that per-respondent correlation; low values flag carelessness.

ref_psychsyn_pairs <- function(x, critval) {
  cor_mat <- stats::cor(x, use = "pairwise.complete.obs")
  p <- ncol(cor_mat)
  pairs <- vector("list", 0L)
  for (j in seq_len(p - 1L)) {
    for (i in seq.int(j + 1L, p)) {
      if (!is.na(cor_mat[i, j]) && cor_mat[i, j] > critval) {
        pairs[[length(pairs) + 1L]] <- c(i, j)
      }
    }
  }
  if (length(pairs) == 0L) {
    return(matrix(integer(0L), ncol = 2L))
  }
  do.call(rbind, pairs)
}

ref_psychsyn_row <- function(row, pairs) {
  if (nrow(pairs) == 0L) return(NA_real_)
  a <- as.numeric(row[pairs[, 1L]])
  b <- as.numeric(row[pairs[, 2L]])
  complete <- !is.na(a) & !is.na(b)
  if (sum(complete) <= 2L) return(NA_real_)
  suppressWarnings(
    stats::cor(a[complete], b[complete], method = "pearson")
  )
}

ref_psychsyn <- function(x, critval = 0.60) {
  if (!is.matrix(x)) x <- as.matrix(x)
  pairs <- ref_psychsyn_pairs(x, critval)
  vapply(seq_len(nrow(x)),
         function(i) ref_psychsyn_row(x[i, ], pairs),
         numeric(1L))
}
