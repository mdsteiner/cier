#' @param items A data.frame of item metadata, one row per item, aligned to the
#'   columns of `responses`. Must carry a `scale` column with at least two
#'   distinct labels. Optional `reverse_keyed` (logical, default none) marks
#'   reverse-keyed items; when any item is reverse-keyed an integer `max` column
#'   (the largest response option, at least `min + 1`) is required. Optional
#'   integer `min` gives the scale base (default `1`). When `responses` has
#'   column names, an optional `item` column (or matching row names) is
#'   cross-checked against them so a reordered frame is a typed error.
