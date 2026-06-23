#' Big-Five responses with independent careless-responding labels
#'
#' A trimmed copy of the Brühlmann et al. (2020) data-quality dataset: 394 crowdsourced
#' respondents who completed the 44-item Big Five Inventory (BFI-44) plus two
#' independent attention checks. The independent labels let index scores be compared
#' against a careless signal the indices did not produce.
#'
#' @format A data frame with 394 rows and 46 columns. The first 44 columns are
#'   the BFI-44 items, named `v_BFI_<scale><n>` where `<scale>` is one of the Big
#'   Five scales (`EX`, `AG`, `CON`, `NEU`, `OP`) and a trailing `_R` marks a
#'   reverse-keyed item; each is an integer 1–5 (5-point Likert). The last two
#'   columns are independent attention checks: `v_Bogus_Item` (integer 1–5, a
#'   bogus item with an obvious correct response) and `v_IRI` (integer 0–7, an
#'   instructed-response item directing the respondent to a specific option). A
#'   deviation on either attention check flags inattentive responding.
#'
#' @details
#' Bundled under CC BY 4.0 (see the package `LICENSE.note`). Only the 44 BFI items and the
#' two careless indicators are retained; no value was altered.
#'
#' Keying caveat: the authors' analysis script applies an extra transformation to one
#' openness item after the usual `_R` reverse-coding (`v_BFI_OP9_R <- 8 - v_BFI_OP9_R`,
#' "wrong item"). The bundled responses are as-is, so name-based reverse-keying
#' (`reverse_keyed = grepl("_R$", ...)`) follows the item names, not the authors'
#' corrected analysis; `vignette("published-results", package = "cier")` shows the
#' difference.
#'
#' @source Brühlmann, F., Petralito, S., Aeschbach, L. F., & Opwis, K. (2020).
#'   Open data on the OSF project \url{https://osf.io/9vjur/}.
#'
#' @references
#' Brühlmann, F., Petralito, S., Aeschbach, L. F., & Opwis, K. (2020). The
#' quality of data collected online: An investigation of careless responding in a
#' crowdsourced sample. *Methods in Psychology*, 2, 100022.
#' \doi{10.1016/j.metip.2020.100022}
#'
#' @examples
#' # The 44 BFI items are the first 44 columns; the last two are the labels.
#' str(bfi_careless[, 45:46])
"bfi_careless"
