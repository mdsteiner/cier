# Build the bundled `bfi_careless` example dataset (data/bfi_careless.rda).
#
# Source: Bruhlmann, Petralito, Aeschbach & Opwis (2020), "The quality of data
# collected online", Methods in Psychology 2, 100022; data on OSF project 9vjur
# (CC BY 4.0). See data-raw/bruhlmann-2020-data-quality-SOURCE.md for the OSF
# download URL + MD5 to re-fetch the original. This reads the vetted archive
# copy of that data (semicolon-free, UTF-8, PII-removed) and trims it to the
# 44 BFI-44 items plus the two independent attention-check indicators.
#
# Re-run from the package root with:
#   Rscript data-raw/bfi_careless.R

src <- "archive/inst/extdata/bruhlmann-2020-data-quality.csv"
raw <- utils::read.csv(src, stringsAsFactors = FALSE)

bfi_items <- grep("^v_BFI_", names(raw), value = TRUE)   # 44 BFI-44 items
stopifnot(length(bfi_items) == 44L,
          all(c("v_Bogus_Item", "v_IRI") %in% names(raw)))

# BFI items first (source order; the `_R` suffix marks reverse-keyed items and
# EX/AG/CON/NEU/OP marks the Big-Five scale), then the two careless labels.
# The `_R` naming is kept as distributed even though the authors' analysis
# corrected v_BFI_OP9_R as mis-keyed -- see the keying note in
# data-raw/bruhlmann-2020-data-quality-SOURCE.md.
bfi_careless <- raw[, c(bfi_items, "v_Bogus_Item", "v_IRI")]
rownames(bfi_careless) <- NULL

dir.create("data", showWarnings = FALSE)
save(bfi_careless, file = "data/bfi_careless.rda", compress = "xz", version = 2L)
