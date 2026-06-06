# Source and licence: the `bfi_careless` example dataset

`data/bfi_careless.rda` is a trimmed copy of the Brühlmann et al. (2020)
data-quality dataset, bundled with `cier` as the example dataset used by the
help-page examples and the getting-started vignette. It carries the 44 BFI-44
personality items plus two **independent** careless-responding indicators (a
bogus item and an instructed-response item), so every v0 index can be
demonstrated and compared against an external label.

## Citation

Brühlmann, F., Petralito, S., Aeschbach, L. F., & Opwis, K. (2020). The quality
of data collected online: An investigation of careless responding in a
crowdsourced sample. *Methods in Psychology*, 2, 100022.
<https://doi.org/10.1016/j.metip.2020.100022>

## Licence

Creative Commons Attribution 4.0 International (CC BY 4.0) — the licence set on
the authors' OSF project, which permits redistribution and modification with
attribution.

- OSF project: <https://osf.io/9vjur/>
- Data component: `Data/data_anon.csv`, direct download
  <https://osf.io/download/ab2mk/> (semicolon-delimited; MD5
  `2a60f58261344a2559c7894366b0773f`; 394 rows x 163 columns).

## Changes made (per CC BY attribution terms)

This bundled object is **not** byte-identical to the original. Starting from the
original `data_anon.csv`, we:

1. removed two free-text / near-identifying columns (`v_lang_blob`,
   `v_Gender_other`) and re-saved as a UTF-8, comma-separated file (this is the
   `archive/inst/extdata/bruhlmann-2020-data-quality.csv` copy this script reads);
2. retained only the 44 `v_BFI_*` items and the two careless indicators
   `v_Bogus_Item` and `v_IRI`, dropping all other columns;
3. saved the result as a lazy-loaded R data object (`data/bfi_careless.rda`).

No response value, validity-indicator value, or timing value was altered. The
build is reproducible with `data-raw/bfi_careless.R`.
