# Paper-faithful reference implementation of the psychometric-ANTONYMS C/IER
# index per Meade & Craig (2012, *Psychological Methods*). Independent of the
# production kernel: it re-derives the statistic from scratch with explicit
# nested loops over the inter-item correlation matrix and a per-respondent
# stacked-pair correlation, and never calls cier's kernel.
#
# The archive carried no dedicated antonyms oracle file (it was inline in the old
# test); this file is the standalone re-derivation.
#
# Definition (the synonyms definition with the threshold sign flipped):
#   1. From the inter-item correlation matrix R (pairwise-complete), identify the
#      set of item pairs P = {(i, j): R_{ij} < -|critval|, i > j} -- the strongly
#      NEGATIVELY correlated pairs (the "antonyms").
#   2. For each respondent, stack the values of the two columns of each pair into
#      two vectors (the larger item index first, the smaller second). Compute the
#      Pearson correlation across those stacked pair values.
#   3. Return that per-respondent correlation; HIGH values flag carelessness (a
#      careful respondent answers genuine antonyms in opposition -> strongly
#      negative; a careless respondent's decouples toward zero / positive).
#
# The per-row stacked-pair correlation is identical to the synonyms case, so this
# file REUSES `ref_psychsyn_row()` from ref-psychsyn-meade-craig-2012.R (the test
# file sources that oracle first). Only the pair-search sign differs -- a
# deliberate, sanctioned re-use.

ref_psychant_pairs <- function(x, critval) {
  cor_mat <- stats::cor(x, use = "pairwise.complete.obs")
  p <- ncol(cor_mat)
  pairs <- vector("list", 0L)
  for (j in seq_len(p - 1L)) {
    for (i in seq.int(j + 1L, p)) {
      if (!is.na(cor_mat[i, j]) && cor_mat[i, j] < -abs(critval)) {
        pairs[[length(pairs) + 1L]] <- c(i, j)
      }
    }
  }
  if (length(pairs) == 0L) {
    return(matrix(integer(0L), ncol = 2L))
  }
  do.call(rbind, pairs)
}

ref_psychant <- function(x, critval = 0.60) {
  if (!is.matrix(x)) x <- as.matrix(x)
  pairs <- ref_psychant_pairs(x, critval)
  vapply(seq_len(nrow(x)),
         function(i) ref_psychsyn_row(x[i, ], pairs),
         numeric(1L))
}

# Clean-reference variant: the antonym pairs are discovered on a
# SEPARATE reference sample `ref` (the strongly NEGATIVELY correlated pairs there),
# while every respondent is still scored on the analysis sample `x` -- the antonym
# mirror of ref_psychsyn_ref(). Independent of the production resolver; reuses the
# antonym pair-search sign and the shared ref_psychsyn_row() scorer, never calling
# cier's kernel.
ref_psychant_ref <- function(x, ref, critval = 0.60) {
  if (!is.matrix(x)) x <- as.matrix(x)
  if (!is.matrix(ref)) ref <- as.matrix(ref)
  pairs <- ref_psychant_pairs(ref, critval)
  vapply(seq_len(nrow(x)),
         function(i) ref_psychsyn_row(x[i, ], pairs),
         numeric(1L))
}
