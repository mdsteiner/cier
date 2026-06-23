# cier 0.1.0

First public release. `cier` detects careless and insufficient-effort responding
(C/IER) in survey data with a lean, auditable set of indices: every number traces
to a published method or a trusted package, and the flagging layer reports a
user-set target flag rate rather than a single take-it-or-leave-it label.

## Detection battery

`cier_screen()` runs ten indirect, response-pattern indices over one response set
and prints a transparent, per-index flag table. Shared constructs collapse to one
vote -- the detection signal lives in multi-index agreement, not any single index:

- Longstring and intra-individual response variability (IRV).
- Even-odd consistency and personal reliability (with the resampled RPR variant).
- Psychometric synonyms (`cier_psychsyn()`) and antonyms (`cier_psychant()`).
- Mahalanobis distance and the person-total correlation (`r_pbis`).
- The nonparametric person-fit statistics Gnormed and Ht, both computed in-package
  with no backend required (`PerFit` and `mokken` are retained only as the tests'
  independent parity oracles).

Most cutoffs default to an empirical percentile at a user-set target flag rate
(`fpr`, 5%), documented as a ranking convention -- not a calibrated false-positive
rate. Longstring keeps a fixed default (half the item count); Mahalanobis distance
(chi-square) and Gnormed (Monte-Carlo null) reference a proper null instead. The
battery follows the indirect indices evaluated by Goldammer et al. (2024).

## Standalone indices

Five further indices are exported as standalone functions; `cier_screen()` does not
run them, because they read inputs the response matrix does not carry, or overlap
the battery in ways still being evaluated:

- `cier_autocorrelation()` and `cier_lazr()` -- sequence-pattern family.
- `cier_total_time()` and `cier_page_time()` -- timing family.
- `cier_attention()` -- direct attention-check family.

## Data simulation

- `cier_simulate()` generates survey data with planted careless responding (a known
  pattern, extent, onset, and prevalence) for power analysis, method comparison, and
  recovery tests. It is a data generator, not an index: it has no cutoff, never
  enters `cier_screen()`, and its output is not evidence of real-world validity.

## Interface and behaviour

- `cier_screen(fpr = )` sweeps the percentile family's target flag rate in one call.
- `summary()` methods for `cier_index`, `cier_screen`, and `cier_sim`; the resolved
  cutoff and its provenance now print alongside each index.
- `cier_psychsyn()` / `cier_psychant()` gain an opt-in `reference` argument that
  discovers the synonym / antonym item pairs on a clean source while still scoring
  every respondent on the full sample.
- Percentile cutoffs abstain (rather than flag everyone) on constant scores or
  samples too small for the target rate, and warn when ties saturate the cutoff.
- `cier_autocorrelation()` defaults to `max_lag = min(n_items - 3, 10)` and abstains
  on respondents with too few answered pairs per lag.
- Input validation hardened: item metadata is cross-checked against the response
  columns; out-of-range values are caught on forward- and reverse-keyed items
  alike; a missing `reverse_keyed` column is reported once rather than silently
  treating every item as forward-keyed; and likely ID / non-item columns are named
  in the error.
- `cier_ht()` scores Ht with an in-package closed form instead of delegating to
  `mokken`. The result is value-identical (to machine precision), but the
  computation is **linear in the number of respondents** (one `O(n * p)` pass plus a
  per-respondent sort over the items, with linear memory) rather than building the
  quadratic person-pair matrix. Two consequences:
  the former 10-category ceiling is **lifted**, so wider scales (11-point, 0--100,
  and mixed-width batteries) now score; and `cier_screen()` runs Ht even when
  `mokken` is not installed. `mokken` stays a suggested package, now used only as the
  independent parity oracle in the tests.
- `cier_ht()` gains an **opt-in model-conforming Monte-Carlo null** cutoff, selectable
  with `method = "mc_null"` and reproducible through a new `seed` argument. It reuses
  the same in-package nonparametric null engine as `cier_gnormed()` -- a
  sum-score-conditional resample of the scored block -- scored by the linear Ht kernel
  on its lower tail, with constant (straightline) response vectors excluded. This
  polytomous Ht null is new (`PerFit` establishes a simulated null only for the
  dichotomous Ht) but is built from `PerFit`'s own mechanism, so it is consistent with
  Gnormed by construction. **The percentile cutoff remains the default**, so existing
  `cier_ht()` calls are unchanged and consume no random draws. On a straightliner-
  dominated sample the simulated null can degenerate (too few non-constant vectors to
  score); `cier_ht(method = "mc_null")` then abstains cleanly -- an `NA` cutoff with a
  typed warning, flagging no one -- rather than erroring.
- `cier_gnormed()` scores the normed Guttman-error statistic with an in-package
  closed form instead of delegating to `PerFit::Gnormed.poly()`, and resolves its
  default Monte-Carlo cutoff with an in-package nonparametric null that faithfully
  reproduces `PerFit::cutoff()`'s default mechanism (a sum-score-conditional
  resample). The score is value-identical to PerFit but **exact** where PerFit rounds
  its output to four decimals, so a Gnormed value can differ from a prior
  `PerFit`-scored run in the fourth decimal. `cier_screen()` now runs Gnormed even
  when `PerFit` is not installed; `PerFit` stays a suggested package, now used only
  as the independent parity oracle in the tests. The default cutoff is reproducible
  through the existing `seed` argument; with `seed = NULL` it varies run to run, as
  before.
