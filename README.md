
<!-- README.md is generated from README.Rmd. Please edit that file. -->

# cier

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/mdsteiner/cier/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mdsteiner/cier/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/mdsteiner/cier/graph/badge.svg)](https://app.codecov.io/gh/mdsteiner/cier)
<!-- badges: end -->

`cier` detects careless and insufficient-effort responding (C/IER) in
survey data using **indirect, response-pattern indices**. It favours a
small, validated set of indices and an auditable,
false-positive-rate-controlled flagging layer over a large method count
– every number traces to a paper or a trusted package. The core battery
follows the indirect indices evaluated by Goldammer et al. (2024);
standalone sequence, timing, and attention-check indices complement it.

## The indices

| Index | Function | Construct family | Flags | Default cutoff |
|----|----|----|----|----|
| Longstring | `cier_longstring` | – | long identical runs | half the items |
| Response variability (IRV) | `cier_irv` | – | low variability | 5% percentile |
| Even-odd consistency | `cier_even_odd` | consistency | low consistency | 5% percentile |
| Personal reliability (PR/RPR) | `cier_personal_reliability` | consistency | low reliability | 5% percentile |
| Psychometric synonyms | `cier_psychsyn` | – | low agreement | 5% percentile |
| Psychometric antonyms | `cier_psychant` | – | high agreement | 5% percentile |
| Mahalanobis distance | `cier_mahalanobis` | – | multivariate outlier | chi-square, *p* \< .001 |
| Person-total correlation | `cier_person_total` | – | low correlation | 5% percentile |
| Gnormed (person-fit) | `cier_gnormed` | person-fit | aberrant pattern | Monte-Carlo null, *p* = .05 |
| Ht (person-fit) | `cier_ht` | person-fit | aberrant pattern | 5% percentile |

`cier_gnormed` / `cier_ht` use the suggested `PerFit` / `mokken`
packages and are skipped with a recorded reason when those are not
installed.

Five further indices are standalone functions – `cier_screen()` does not
run them (yet): they read inputs the response matrix does not carry
(timestamps, check keys) or overlap with the battery in ways still being
evaluated.

| Index | Function | Family | Flags | Default cutoff |
|----|----|----|----|----|
| Autocorrelation | `cier_autocorrelation` | sequence pattern | repetitive / periodic answers | 5% percentile |
| Laz.R (Markov predictability) | `cier_lazr` | sequence pattern | predictable transitions | 5% percentile |
| Total time | `cier_total_time` | timing | fast completion | 5% percentile (lower) |
| Page time | `cier_page_time` | timing | pages under 2 s/item | any rapid page |
| Attention checks | `cier_attention` | direct | failed checks | any failed check |

## Installation

`cier` is not yet on CRAN. Install the development version from GitHub –
with `dependencies = TRUE` so the suggested `PerFit` / `mokken` (which
back Gnormed and Ht) come along:

``` r
# pak is the current r-lib installer:
# install.packages("pak")
pak::pak("mdsteiner/cier", dependencies = TRUE)
```

## Quick start

`cier_screen()` runs the whole battery over one response set and prints
a transparent flag table. The four metadata-using indices need an
`items` frame aligned to the columns; the bundled `bfi_careless` makes
one easy to build from the column names:

``` r
library(cier)

resp <- bfi_careless[, 1:44]
nm   <- names(resp)
items <- data.frame(
  scale         = sub("^v_BFI_([A-Za-z]+)[0-9].*$", "\\1", nm),
  reverse_keyed = grepl("_R$", nm),
  max           = 5L
)

cier_screen(resp, items,
  control = list(cier_personal_reliability = list(seed = 1),
                 cier_gnormed = list(seed = 1)))
#> ── cier_screen ─────────────────────────────────────────────────────────────────
#> 10 indices on 394 respondents; 9 votes after collapsing shared constructs.
#> 
#> Per-index flags (transparent; not independent votes):
#>   cier_longstring            25 / 394 (6.3%)
#>   cier_irv                   20 / 394 (5.1%)
#>   cier_even_odd              37 / 375 (9.9%)  [consistency]
#>   cier_personal_reliability  19 / 377 (5.0%)  [consistency]
#>   cier_psychsyn              0 / 0 (--%)
#>   cier_psychant              0 / 0 (--%)
#>   cier_mahalanobis           33 / 394 (8.4%)
#>   cier_person_total          19 / 376 (5.1%)
#>   cier_gnormed               47 / 394 (11.9%)
#>   cier_ht                    19 / 377 (5.0%)
#> 
#> Cross-index agreement (observed vs independence baseline):
#>   flagged by >= 1 vote: 126 / 394 (32.0%); expected 41.3%
#>   flagged by >= 2 votes: 56 / 394 (14.2%); expected 8.6%  <- excess
#>   flagged by >= 3 votes: 15 / 394 (3.8%); expected 1.0%  <- excess
#>   flagged by >= 4 votes: 3 / 394 (0.8%); expected 0.1%  <- excess
#>   flagged by >= 5 votes: 1 / 394 (0.3%); expected 0.0%  <- excess
#>   flagged by >= 6 votes: 0 / 394 (0.0%); expected 0.0%
#>   flagged by >= 7 votes: 0 / 394 (0.0%); expected 0.0%
#>   flagged by >= 8 votes: 0 / 394 (0.0%); expected 0.0%
#>   flagged by >= 9 votes: 0 / 394 (0.0%); expected 0.0%
#> 
#> Skipped: 0
```

## How the flagging works

The detection signal lives in **multi-index agreement**, not any single
index. Even-odd, PR and RPR measure one construct, so they collapse to
**one vote** – never weighted as independent evidence. There is no
trained classifier in this release.

Three things to keep in mind:

- The percentile cutoff is a **ranking convention**, not a calibrated
  false-positive rate – it flags roughly `fpr`% of *this* sample by
  construction (more when coarse-scale responses tie at the cutoff), and
  is not Goldammer’s simulated-null Sen95 operating point. (Mahalanobis
  and Gnormed are the exceptions: they reference a calibrated null
  instead.)
- The bundled attention-check labels are a **partial, biased** signal,
  not ground truth – design checks in advance and do not treat them as a
  key.
- **Report results before and after exclusion**, and show the flag rate
  across `fpr` (for example 0.01, 0.05, 0.10).

See `vignette("cier")` for the full walk-through.

## Licence

Package code is MIT-licensed. The bundled `bfi_careless` dataset is
redistributed under CC BY 4.0 (Brühlmann et al., 2020); see
`LICENSE.note`.

## References

Brühlmann, F., Petralito, S., Aeschbach, L. F., & Opwis, K. (2020). The
quality of data collected online: An investigation of careless
responding in a crowdsourced sample. *Methods in Psychology*, 2, 100022.

Goldammer P., Stöckli, P. L., Escher, Y. A., Annen, H., Jonas, K., &
Antonakis, J. (2024) Careless responding detection revisited: Accuracy
of direct and indirect measures. *Behavior Research Methods*, 56,
8422-8449.
