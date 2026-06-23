# Shared kernels for the timing family. cier_total_time has no kernel (its value is the
# validated seconds vector itself), so this file opens with cier_page_time.

# Per-respondent count of pages whose mean per-item time is strictly below `min_seconds`
# (Bowling, Huang, Brower & Bragg 2023). Returns numeric length-n in [0, pages]; a
# respondent with no timed page abstains (NA). Untimed pages contribute no evidence.
kernel_page_time <- function(page_seconds, items_per_page, min_seconds) {
  # Column-major recycling matches sweep()'s per-column divide (incl. NA propagation)
  # without the generic overhead.
  per_item <- page_seconds / rep(items_per_page, each = nrow(page_seconds))
  observed <- !is.na(per_item)
  rapid <- observed & per_item < min_seconds
  value <- as.numeric(rowSums(rapid))
  value[rowSums(observed) == 0L] <- NA_real_
  unname(value)
}
