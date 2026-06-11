# Architecture decisions

The durable record of binding design decisions for cier, in plain terms. It is
deliberately short: one entry per genuinely binding choice — an object schema, a
default cutoff, a numerical tolerance, a cross-package parity contract, or a
scope boundary — not one per feature. Everything else belongs in a commit
message. New entries are added as such decisions are made during the build.

## Scope boundary

The v0 battery is the indirect, response-pattern C/IER indices evaluated by
Goldammer et al. (2024): personal reliability (including the resampled variant),
psychometric synonyms and antonyms, Mahalanobis distance, person-total
correlation, the nonparametric person-fit statistics Gnormed and Ht, and the
classic longstring and intra-individual response variability indices.

**v0.2 amendment (signed off 2026-06-10).** The boundary is extended to add five
published, cheap, oracle-able indices the bundle study needs on real data:

- `cier_autocorrelation` (indirect; Gottfried et al. 2022, PARE 27(2)) -- max
  absolute lag autocorrelation; targets the repetitive / periodic patterns the
  consistency family abstains on.
- `cier_lazr` (indirect; Biemann et al. 2025, ORM) -- first-order Markov response
  predictability; generalises longstring.
- `cier_total_time` (new **timing** family; Ward & Meade 2023; Huang et al. 2012)
  -- total completion time; the strongest single index on fully-careless cells.
- `cier_page_time` (timing; Bowling et al. 2023, ORM) -- count of pages faster
  than 2 s per item; catches within-survey speeding bursts total time misses.
- `cier_attention` (new **direct** family; Meade & Craig 2012; Goldammer et al.
  2024) -- count of failed instructed / bogus / infrequency checks.

Explicitly **out** (each evaluated study-side or excluded by prior evidence, not
shipped): the repetition index (the study calls `responsePatterns::rp.patterns()`
directly; archive evidence has it tracking Laz.R within ~.03 and it was the
archive's worst CPU hotspot), Stan time mixtures (archive measured AUC identical
to plain total time at MCMC cost), IRT person-fit / INDCHI / ML classifiers (weak
or below-chance on multi-scale real data, or label-dependent), self-report effort
**as an index** (a single user-thresholded column needs no statistic; it is a
label channel in the study), and the missingness / omission rate (one line of
base R, documented in the vignette not shipped).

Adding any further method family (IRT person-fit, model-based, machine-learning),
a learned combiner (including the post-study `cier_recommended()` bundle support),
a simulation engine, a specification-curve tool, or a report generator still
requires a new entry here and explicit sign-off before code.

## v0.2 index additions: families, cutoffs, deviations

The five v0.2 rows register metadata only; the wrappers land in later,
separately signed-off steps. The registry-encoded decisions and the deliberate
deviations from the source papers and the archived previous version:

- **Families.** Autocorrelation and Laz.R stay `indirect` (response-pattern
  indices); `timing` (total / page time) and `direct` (attention) are new family
  vocabulary levels (`cier_family_levels()`).
- **Screen wiring deferred.** All five ship `screenable = FALSE` with own-id vote
  groups, so the default `cier_screen()` run set is unchanged (the ten v0
  indices). The eventual "repetition" vote group (autocorrelation + Laz.R +
  longstring measure one construct) and the screen's `times` / `pages` / `checks`
  arguments are a post-study decision, not made here.
- **Corrected Gottfried DOI.** The archived registry cited the wrong DOI
  (`10.1177/00131644211046302`, an EPM DOI) for autocorrelation; the correct
  source is PARE 27(2), `10.7275/vyxb-gt24`.
- **Autocorrelation cutoff.** Empirical upper percentile -- this is the paper's
  *own* relative-ranking recommendation, so it is **not** a divergence.
- **Laz.R cutoff.** Empirical upper percentile, a **deliberate divergence** from
  Biemann et al.'s sample-specific Kneedle elbow: it is consistent with cier's
  single-`fpr`-knob ranking convention, and the score's documented
  sequence-length dependence makes an absolute cutoff indefensible. The study
  separately evaluates a paper-faithful Kneedle variant; the package may adapt.
- **Laz.R missing-data convention (resolved 2026-06-10 at the wrapper sign-off).**
  `cier_lazr` **drops NA transitions**: a transition is counted only when both
  endpoints are present, the denominator is the count of valid transitions (not
  the item count minus one), NA is never a Markov state, and a respondent with
  fewer than two valid transitions abstains (`NA`). This is the more defensible
  behaviour: the paper's footnote-2 reference code uses `useNA = "ifany"`, which
  makes a gap a *predictable* state and so scores a careful but incomplete
  respondent (a blank tail) as near-maximally careless -- and the `< 2`-observed
  guard cannot catch it. Drop-NA matches cier's abstain-on-missing house style
  and the autocorrelation sibling's pairwise handling. The study separately
  evaluates the paper-faithful NA-as-state variant; the registry row encodes none
  of this (it is a kernel decision).
- **Laz.R input shape: matrix-only, anchor-count-invariant.** `cier_lazr` takes
  only the response matrix -- no `items` -- exactly like `cier_autocorrelation`.
  The Laz.R value is **invariant to the assumed anchor count `s`** (it is
  determined purely by the observed transition counts, so an unobserved or higher
  anchor adds an all-zero row/column and changes nothing) and to the scale base
  (0-based and bipolar codings score identically). `items` / `s` would therefore
  buy only validation, not statistics, so they are omitted; non-integer responses
  are a typed error and the same-answer-scale / administration-order assumption is
  documented in the help page rather than enforced. This is a deliberate
  divergence from the plan's earlier optional-`items` proposal (written before the
  matrix-only autocorrelation sibling and against the superseded `categories`
  schema).
- **Total-time cutoff.** `percentile`, lower tail, fpr 0.05 -- the uniform knob,
  a divergence from the archive's stricter 0.01. A third mutually-exclusive
  cutoff override `frac_median` (flag respondents faster than a fraction of the
  sample median; Leiner 2019 RSI; Greszki et al. 2015) -- extending the two-knob
  override pattern to three for this index -- is recorded here as forthcoming and
  is built in the total-time wrapper, not now.
- **Attention citation.** The row cites Meade & Craig (2012), `10.1037/a0028085`
  (the foundational attention-check source, already a registered DOI); Goldammer
  et al. (2024) is named in the row notes, so the attention wrapper's help page
  cites it without a `\doi` (the references-DOI guard only constrains DOIs that
  appear in a references block).

## Timing family: a per-respondent seconds vector and a three-knob cutoff override

`cier_total_time` opens the timing family and fixes its input contract. Unlike
every response-pattern index it takes **no response matrix**: the argument
`seconds` is a bare numeric vector of one total completion time per respondent --
the shape survey platforms export. A two-dimensional input (matrix or data frame)
is a typed `cier_error_input`: which axis is the respondent is ambiguous, and
summing per-cell times to a per-respondent total is the user's one line of base R.
Each observed time must be **strictly positive** (a zero or negative duration
errors; `NaN` / infinite error -- and `is.na()` is `TRUE` for `NaN`, so the kernel
tests `is.nan()` / `is.infinite()` on the raw vector *before* dropping genuine
`NA`); a missing time (`NA`) abstains. This is a deliberate break from the archived
previous version, which summed a per-cell `times` matrix inside the heavy
`cier_data()` pipeline via a `kernel_total_time`. The lean input is pre-summed, so
there is **no kernel and no statistic**: the per-respondent value is the validated
`seconds` vector itself (an identity). The single-kernel rule governs statistics,
not an identity, so no `index-kernels-timing.R` is introduced here; the shared
timing kernel file lands with `cier_page_time`, the first timing index that
actually computes something.

The default cutoff is the empirical **lower** percentile at `fpr = 0.05` (low
totals flag speeders) -- the uniform knob, a divergence from the archive's stricter
0.01 (recorded under "v0.2 index additions"). Total time adds a **third**
mutually-exclusive cutoff override, `frac_median`: a median-relative rule that
flags respondents faster than a fraction of the **sample median** (Leiner's 2019
Relative Speed Index, `frac_median = 0.5`; Greszki et al. 2015 at 0.5 / 0.4 / 0.3).
It is anchored to the median, so it is robust to up to half the sample responding
carelessly -- exactly where the empirical percentile, which flags `fpr` by
construction, is least defensible. Its domain is `(0, 1]` (a fraction of the
median; above 1 would flag the slower-than-median half, nonsensical for a speeding
index), and it flags with the package's uniform lower comparator
(`value <= cutoff`) -- "faster than" is operationalised as `<=`, immaterial on
continuous times and keeping one comparator path. Because the median is
data-dependent, `frac_median` resolves through a dedicated **override** resolver,
`resolve_median_cutoff()` in `cutoff.R` (dispatched inline by the wrapper, as
`cier_longstring` dispatches its `frac`), not the value-only `resolve_cutoff()`; on
an all-missing vector it abstains to `NA` with the same
`cier_warning_insufficient_items` as the percentile path. The "every rate-based
**default** resolves through `resolve_cutoff`" rule is unaffected -- `frac_median`
is an override, not a default.

## Page time: a page-totals matrix, a count value, and a proportion override

`cier_page_time(page_seconds, items_per_page, min_seconds = 2, frac, cutoff)` is
the first timing index that actually computes something, so it opens the shared
timing kernel file `R/index-kernels-timing.R` (`cier_total_time` is an identity
with no kernel; the single-kernel rule governs statistics, not identities). Its
input contract and the deliberate deviations from the archived previous version:

- **Lean input: a page-totals matrix plus an explicit item-count vector.**
  `page_seconds` is an `n x pages` matrix of each respondent's **total** time on
  each page (one column per page -- the page-submit timer survey platforms
  export); `items_per_page` is the per-page item count. The per-item rate the
  Bowling et al. (2023) rule thresholds is the page total divided by
  `items_per_page[j]`. This replaces the archive's heavy path (a `cier_data`
  object carrying a per-cell `times` matrix, with page boundaries inferred from
  `items$page`) -- the same break `cier_total_time` made from the pipeline. The
  archive's page-boundary inference and its `cier_warning_page_fallback`
  (one-item-per-page fallback when `items$page` was absent) are **dropped**:
  pages are explicit columns now, so there is nothing to infer and no fallback to
  warn about.
- **Page-level NA, not cell-level.** A page with no recorded time (`NA`)
  contributes no evidence -- it is neither counted rapid nor counted toward the
  timed-page total -- and the per-item denominator is the **declared**
  `items_per_page`, not a count of answered cells (the archive averaged over
  answered items within a page). A respondent with at least one timed page scores
  a finite count; a respondent whose every page is `NA` abstains (`value = NA`).
  This is the deliberate simplification the one-number-per-page input buys.
- **The value is a count; the cutoff is fixed, with a proportion override.**
  Direction `upper`; the per-respondent value is the rapid-page **count** (lean
  bare-numeric kernel -- the archive's `fastest_page_z` / `n_pages_used` diag
  by-products are dropped, as the light `cier_index` has no diag field). The
  default cutoff is the cited `fixed = 1` (any rapid page flags). Because a single
  rapid page can be over-sensitive on a long survey, the wrapper exposes the
  **two-knob override pattern longstring already uses**: `frac`, a fraction of the
  **total page count** in `(0, 1]` resolving to `ceiling(frac * pages)` through
  the existing `resolve_fixed_cutoff(value, n_items)`, and a literal `cutoff`
  count in `[1, pages]` -- mutually exclusive via `assert_single_override`. The
  `frac` denominator is the survey's page count (`ncol`), a single sample-level
  cutoff, not a per-respondent count of answered pages. No `fpr` / percentile knob
  exists (the default is an absolute count, not an empirical tail), so an
  abstaining respondent never routes the cutoff through the percentile
  abstention. Strictly-positive page times are required (zero / negative / `NaN` /
  infinite are typed input errors; `NA` abstains), mirroring `check_seconds`.
- **Oracle-only trust.** No CRAN package implements page time as a C/IER index
  (verified 2026-06-10), so the parity check is the hand-rolled counting-rule
  oracle (`ref-page-time.R`) at tolerance 0 (exact integer counts), like
  `cier_total_time` / PR / RPR. Recorded in `tests/reference/TOLERANCES.md`.

## Cutoff philosophy

There is no ground truth in applied use, and no label-free rule can validate its
own cutoff. The default flags respondents at an empirical percentile for a
user-set target false-positive rate; this is a ranking convention, not a
calibrated rate, and is documented as such. A conforming-null reference is used
where it is cheap and principled (the chi-square distribution for Mahalanobis
distance; a Monte-Carlo null for the person-fit statistics). Fixed absolute
heuristics are not used as defaults. Detection relies on agreement across
indices rather than any single threshold.

## Cutoff resolution: one path, a single direction flip

All cutoffs resolve through one function, `resolve_cutoff()`, which branches on
the registry's `default_cutoff_method`: `percentile` (an empirical quantile at a
target false-positive rate `fpr`, default 0.05), `chisq` (`qchisq(1 - alpha, df)`
for Mahalanobis distance, `alpha` = 0.001), and `fixed` (a literal threshold, or
a fraction of the item count `ceiling(frac * ncol)` — longstring's default
`frac = 0.5`).

The percentile method applies the direction flip **exactly once**: `upper` flags
the high tail (the `1 - fpr` quantile), `lower` flags the low tail (the `fpr`
quantile). The registry stores the target false-positive tail mass `fpr`
(`0.05`) for **every** percentile index regardless of direction; the
`flag_direction` column selects the tail and `resolve_cutoff()` performs the
single flip. This is the deliberate fix for the v1 footgun where a stored
directional quantile (`0.95` for an upper index) was reused as `fpr` and
re-flipped via `1 - p`, landing on the wrong tail. A single `fpr` knob therefore
means the same target tail mass for every index, which is what lets the analyst
sweep `fpr` across {0.01, 0.05, 0.10}. The flag comparator (`>=` / `<=`) is
applied separately, by `apply_flag()`.

`resolve_cutoff()` returns a bare numeric scalar (`NA_real_`, with a typed
`cier_warning_insufficient_items`, when a percentile cutoff cannot be resolved
because no finite values remain). It builds no object: the index wrapper already
knows its method and direction from its registry row and attaches them to the
light `cier_index`.

## Agreement diagnostic: observed co-occurrence vs a Poisson-binomial null

Because an empirical-percentile cutoff flags its target rate by construction, the
per-index flag rate is tautological and is never presented as a false-positive
rate. The informative quantity is cross-index agreement: `flag_agreement()`
reports, for each level k, the observed share of respondents flagged by at least
k votes against the share expected if the votes fired independently. That
expectation is the exact upper tail of the **Poisson-binomial** distribution of
the per-vote flag rates (a sum of independent Bernoullis with *unequal*
probabilities — the correct null when, say, Mahalanobis flags ~0.1% and IRV ~5%;
a plain binomial assuming one shared rate would be wrong). Observed far above
expected makes a clustered careless subgroup visible. This is a descriptive
visibility diagnostic, not a calibrated test: the baseline is a null of "no
shared signal", not a claim the indices are independent (they are not, which is
why even-odd, PR, and RPR collapse to one vote upstream). A companion per-vote
table additionally reports each vote's excess over a supplied calibrated null;
that excess is informative only for the null-referenced indices (Mahalanobis
chi-square, the person-fit Monte-Carlo nulls) and is marked tautological (NA)
for the empirical-percentile votes.

The printed table's `<- excess` marker is **gated on chance**, not on a strict
point comparison. The expectation is an exact population quantity while the
observed share is a sample proportion, so under true independence the observed
exceeds the expected about half the time at each k — a bare `observed >
expected` marker fired on ~50% of clean, contamination-free screens (and on
rounding-invisible gaps like "0.3% vs expected 0.0%"). Each respondent is
independently flagged by >= k votes with probability `expected[k]` under the
null, so the observed count is Binomial(n, `expected[k]`); the marker now fires
only when the one-sided binomial tail P(count >= observed) is below 0.05. Still
descriptive (the per-k rows are not independent of each other; no multiplicity
correction), but it no longer advertises ordinary sampling noise as
contamination.

## Architecture: function-first indices

Each index is a documented function on a response matrix (a data.frame or
tibble is accepted and coerced internally, so users need not call
`as.matrix()`), not a method behind an input-object pipeline. The four indices
that need item metadata (even-odd, personal reliability, Gnormed, Ht) take a
single optional `items` data.frame with columns `scale`, `reverse_keyed`,
`categories`, and (optionally) `min`, the scale base for 0-based or bipolar
codings (one row per item); the other six need only the responses. An index
returns a light `cier_index` — a list-based S3 object (see "Object schema:
list-based `cier_index`" below) — assembled by one shared `new_cier_index()`
constructor. `cier_screen()` is a thin orchestrator that runs the selected
indices and returns a per-respondent flag table plus the count of flagged
*constructs* (even-odd and personal reliability collapse to one vote). There is
deliberately no `cier_data` / `cier_items` input class and no per-method or
heavy `validate_` layer: the matrix-plus-metadata convention keeps the surface
small and pushes validation to small per-function input checks. Cross-cutting
concerns stay shared and are the only retained foundation: typed conditions, the
method-properties registry, one cutoff resolver, and the one `cier_index`
constructor.

## Object schema: list-based `cier_index`

A `cier_index` is a list-based S3 object —
`structure(list(value, flagged, method, cutoff, direction), class = "cier_index")`
— not a `data.frame` carrying the metadata as attributes. Per-respondent `value`
and `flagged` are vectors (`NA` where the index abstains); the flag count and
rate are derived on `print` from `flagged` and are never stored, so they cannot
desynchronise from the data. One shared `new_cier_index()` is the single schema
definition for all indices and enforces the universal rule that `flagged` is
`NA` wherever `value` is `NA`. An `as.data.frame()` method returns the tidy
`data.frame(value, flagged)` for downstream analysis.

The list shape was chosen over data.frame-plus-attributes on the evidence that
custom attributes on a user-facing data.frame are silently dropped by
tidyverse/base operations and survive (stale) across row-subsetting — a real
defect for a research tool whose users pipe per-respondent scores (a subset
printed e.g. "2 of 1 respondent (200%)"). A list keeps the metadata robust and
discoverable (`out$method` works), is the canonical S3 record shape, and matches
the objects returned by `PerFit` / `mokken` and the prior version of the
package. `cier_screen()` (slice 12) follows the same robust shape.

## Cutoff overrides: a rate and a literal, mutually exclusive

Each index exposes its cutoff override as **two mutually-exclusive arguments**: a
rate (`fpr` for percentile indices, `alpha` for Mahalanobis; `frac`, a fraction
of the item count, for longstring) and a literal `cutoff` on the score. A single
overloaded argument was rejected: a value in `(0, 1]` cannot distinguish a rate
from a literal threshold for indices whose scores fall in that range (IRV SDs,
correlations, Gnormed), so e.g. an IRV `cutoff = 0.5` would be ambiguous.

**Validation lives at the public boundary; the resolver trusts its inputs.** Each
wrapper validates everything up front so a bad argument fails before the kernel
runs: `check_open_unit()` for a rate in the open interval `(0, 1)` (`fpr` /
`alpha`), `check_fraction()` for a fraction in `(0, 1]` (`frac`),
`check_number()` for a literal `cutoff` against index-specific bounds (`[1, p]`
for longstring, `[0, Inf)` for a non-negative score). Immediately after those
checks the wrapper calls `assert_single_override()` (reject both knobs), then the
kernel runs. The cutoff dispatch is **a literal `cutoff` used verbatim** (already
validated); otherwise the rate-based default resolved through the one resolver,
`resolve_cutoff(method = <registry method>, …)`. `resolve_cutoff()` and its
private helpers are internal and do **no** input re-checking (their `method` /
`direction` come from the registry and the rate/literal are wrapper-validated) —
they only do the math and signal the runtime percentile abstention
(`NA_real_` + `cier_warning_insufficient_items`).

Where that dispatch lives depends on the cutoff family. The **percentile**
indices (IRV, even-odd, person-total, personal reliability, Ht, and — via the
no-pairs-aware tail — psychsyn / psychant) share one helper,
`resolve_index_cutoff(value, row, fpr, cutoff, call)`, which composes the
dispatch with `apply_flag()` and `new_cier_index()` into the common
cutoff → flag → assemble tail; their wrappers end in that single call. The
non-percentile wrappers — longstring (`fixed`), Mahalanobis (`chisq`), and
Gnormed (`perfit_null`, resolved at the bridge from the fitted object) — inline
their own two-line dispatch, since each computes its cutoff differently. (An
earlier revision of this section described `resolve_index_cutoff()` as removed
in favour of fully-inlined dispatch; the shared helper was reinstated — without
the `rate_fn` closure it once took — when the percentile tail turned out to be
repeated verbatim across wrappers.) The one-cutoff-path rule holds throughout:
every rate-based **default** resolves through the single `resolve_cutoff()`.
Putting all validation in the public function gives the earliest possible
failure. The cost is that a future wrapper which forgets to validate would pass
bad input silently rather than get a typed error; the mitigation is the
wrapper-validation convention plus per-wrapper input-error tests.

`cier_total_time` is the one index with **three** mutually-exclusive knobs -- the
rate `fpr`, the median-relative `frac_median`, and the literal `cutoff` (see
"Timing family" above). They are guarded by `assert_single_cutoff()`, the n-way
generalisation of `assert_single_override()`: it takes a named list of the knob
values and aborts naming exactly the pair or triple supplied. The two-argument
`assert_single_override()` stays in use for the two-knob indices; the n-way form
is reached for only when an index exposes more than two ways to set its cutoff.

## Method-properties registry schema

The registry (`inst/extdata/method-properties.csv`) is the single source of truth
for cutoff defaults, flag direction, backend, screen membership, and the screen's
vote grouping. Its columns are `method`, `family`, `paper_year`,
`paper_citation_key`, `doi`, `default_cutoff_method`, `default_cutoff_value`,
`flag_direction`, `companion_methods`, `backend`, `screenable`, `vote_group`, and
`notes`. `vote_group` (added in the screen slice) labels the construct each index
votes for: indices sharing a label collapse to one vote in `cier_screen()`. Only
even-odd and personal reliability share a group (`consistency`) — they measure one
construct, so counting them as two votes would double-count; every other index's
`vote_group` is its own id. Keeping the grouping in the registry (rather than a
hard-coded map) makes vote membership data-driven and validated alongside
`screenable`. The boolean capability
flags (`requires_*`) and the `available` column are deliberately omitted: no v0
component reads them (each index's wrapper is the single source of truth for the
metadata it requires), the nonparametric battery makes an IRT-model flag
constant, and every shipped row is available. They will be re-added only if a
later method family needs them. `cier_personal_reliability` cites the
resampled-reliability paper (Goldammer et al. 2024); personal reliability's origin
(Jackson 1976) is recorded in the row's notes and the help-page references.

## Person-fit backends

The two nonparametric person-fit indices use different optional backends because
each is the only package providing the polytomous form survey data needs: normed
Guttman errors (Gnormed) via `PerFit`, person scalability (Ht) via `mokken`.
`PerFit::Ht` is dichotomous-only — it aborts on non-binary input — so it cannot
serve polytomous Likert responses; `mokken::coefH` on the transposed scale is the
polytomous Ht that matches the reference values. Both packages are optional
(`Suggests`), with a graceful skip when absent.

**Backend limits are typed and screen-survivable.** A backend's hard ceiling on
otherwise-valid data — concretely, `mokken`'s 10-category limit (`coefH` raw-
stops when the global zero-based range exceeds 9, so an 11-point or 0–100 item
cannot be scored) — is converted at the bridge into a typed `cier_error_input`
carrying an extra `cier_error_backend_limit` subclass, mirroring how
`personfit_zero_base()` converts `PerFit`'s terse aborts. The subclass exists
for `cier_screen()`: the screen catches **exactly** that class and records the
index as skipped-with-reason (the condition's `data$reason`), so one index's
backend ceiling cannot crash a ten-index battery; every other error — a
malformed `items` frame included — still propagates. A direct `cier_ht()` call
gets the typed error with the limit and remedy spelled out.

## Gnormed cutoff: PerFit Monte-Carlo null

Gnormed's default cutoff is **not** a sample percentile but the `PerFit`
Monte-Carlo null (`PerFit::cutoff`): `PerFit` resamples model-conforming response
vectors from the fitted object and takes the nominal-rate (`Blvl = fpr`) quantile
of the statistic. A real null was chosen over the percentile ranking convention
because, for a person-fit index, the null-referenced flag rate is **informative**
— it can exceed `fpr` under contamination, whereas a sample percentile flags `fpr`
by construction. The registry records this as a new `default_cutoff_method`,
`perfit_null` (value `0.05`, the nominal level), so the registry stays the honest
single source of the default.

Three consequences follow. (1) `perfit_null` is the package's first
**simulation-referenced** cutoff, so it is randomised; the wrapper exposes a
`seed` argument and applies it **locally** (saving and restoring the caller's
`.Random.seed`, exactly as `kernel_rpr` does), so a seeded call is reproducible
without disturbing the session RNG. (2) Unlike the value-only strategies
(`percentile` / `fixed` / `chisq`), the null needs the **fitted object**, not just
the score vector, so it is resolved at the bridge (`resolve_gnormed_cutoff` /
`resolve_perfit_null_cutoff`) rather than inside the value-only `resolve_cutoff()`
— the one documented exception to "all cutoffs resolve through `resolve_cutoff`".
The kernel therefore returns `list(value, fit)` so the simulation reuses the fit
rather than refitting. (3) `PerFit` is **required** for the standalone index
(scoring is impossible without it), so its absence is a typed
`cier_error_input`, not the architecture's "fall back to the percentile" path —
that fallback is unreachable here because there are no Gnormed values to rank when
`PerFit` is absent.

## Mahalanobis degenerate covariance: warn and abstain

The Mahalanobis index needs a sample covariance to compute any distance. When
none can be estimated — fewer than two respondents carry data, or the covariance
matrix is not invertible (more items than respondents or perfectly collinear
items, so `solve()` fails; or a pairwise covariance with undefined entries — two
items never co-answered, or an item answered by nobody — which yields a non-finite
inverse rather than a `solve()` error) — every respondent's `value` is `NA` and
no one is flagged, and
the wrapper raises a typed `cier_warning_singular_covariance` whose structured
`data$reason` names the cause (`"insufficient_responses"`,
`"singular_covariance"`, or `"indefinite_covariance"`). A silent all-`NA` result
was rejected: a singular
covariance is a substantive, actionable analyst event (a collinear or
over-wide item set), so it earns its own condition class rather than being folded
into the cutoff layer's `cier_warning_insufficient_items`. This is the only
wholesale-abstention case that warns; an individual all-`NA` respondent (alongside
others who answered) is ordinary per-row abstention and stays silent, as in every
other index. The kernel itself stays pure — it returns a status code, and the
wrapper raises the condition.

The third cause, an **indefinite** pairwise covariance, was added after review:
pairwise estimation assembles each covariance cell from a different subsample, so
under heavy or structured missingness the cells can be mutually inconsistent and
the matrix gains a negative eigenvalue while `solve()` still succeeds. The
bilinear form is then signed — a respondent can score a *negative* "squared
distance" that the upper-tail chi-square flag can never reach, and the ranking
among the positive rows is distorted too — so the distance is invalid for every
row and the kernel abstains wholesale (`chol()` is the test: it errors iff the
matrix is not positive definite; the inverse still comes from `solve()`, keeping
every positive-definite input byte-identical to the parity partners). Two
alternatives were evaluated and rejected. **Repairing the matrix** (projecting to
the nearest positive-definite matrix, e.g. `Matrix::nearPD()`) produces a
statistic that traces to no cited paper and no trusted package — it would diverge
*by construction* from `careless::mahad()` / `psych::outlier()` exactly on these
inputs, so the cross-package parity layer could no longer pin the kernel, and it
silently masks the data problem (the missingness pattern) the researcher needs to
see; the projection also carries its own tuning knobs (`eig.tol`, `conv.tol`)
with no reference value to validate against. **A hard abort** is inconsistent
with the established degenerate-covariance contract (singular Σ warns and
abstains) and would kill an entire `cier_screen()` run over one index's data
problem. Note that `careless::mahad()` itself returns the signed values on such
input, so the parity suite cannot guard this path — the dedicated
indefinite-fixture regression test is the only guard.

## Psychsyn/psychant kernel: vectorise (masked-sum), relax careless parity to 1e-12

The psychometric-synonyms / antonyms kernel originally scored each respondent with
a per-row `stats::cor()` call inside a `vapply` loop, which matched
`careless::psychsyn(resample_na = FALSE)` **bytewise** because `careless`'s own
`syn_for_one` scores `cor()` on the same full stacked pair vectors. Profiling
showed that per-row loop is ~85% of the call (one `cor()` per respondent); it
makes the index roughly 4–6× slower than necessary and scales linearly with the
respondent count.

The kernel now scores every respondent in one vectorised pass: a masked-sum
Pearson over the `n × K` stacked pair matrices (the same technique
`kernel_person_total` uses), with the two variance terms clamped at zero so
floating noise on a zero-variance pair side cannot send `sqrt()` to `NaN` with a
warning. The statistic is unchanged — it is the same Pearson correlation — but the
masked-sum sums in a different order, so the per-respondent scores move by
≤ 1.1e-13. The cross-package parity with `careless::psychsyn(resample_na = FALSE)`
is therefore **relaxed from bytewise (`0`) to `1e-12`** in `TOLERANCES.md`, and the
test assertion moves from `expect_equal(tolerance = 0)` to `tolerance = 1e-12`.
The independent-oracle parity (`ref_psychsyn`, which pre-filters complete cases) is
unaffected — the vectorised kernel matches it to ~5e-15, still well inside 1e-12.

The trade was accepted deliberately (Markus, this slice): the bytewise guarantee
buys exactness against one external package, whereas the vectorisation buys a 4–6×
speedup, sub-100 ms scoring at every realistic respondent count (so **no C++ /
`src/`** is warranted — adding a compiled backend would impose a Rtools/portability
burden against the "trace every number" mandate for no practical gain), and
uniformity with the already-vectorised `kernel_person_total` (itself `1e-12` vs
`PerFit`, for the same summation-order reason). The kernel stays a single pure
function shared by psychsyn (`pairing = "syn"`) and psychant (`"ant"`), so the
antonym index inherits both the speedup and the same tolerance.

A no-pairs result (no inter-item correlation clears `critical_r`, the common
case on broad inventories at the 0.60 default) warns with a **tailored** typed
condition, `cier_warning_no_pairs`, naming `critical_r`, the strongest in-tail
correlation, and the `cier_psychsyn_critval()` sweep — instead of the generic
percentile abstention ("no finite values remain"), which names neither cause
nor remedy. The subclass deliberately also carries
`cier_warning_insufficient_items` so `cier_screen()`'s targeted muffler keeps
covering it (the screen reports the case transparently as "0 / 0").

## Zero-variance abstention: exact, not cancellation-dependent

The two masked-sum Pearson kernels (`kernel_person_total`, `kernel_psychsyn`)
document that a constant (straightliner) respondent — or constant pair side —
has an undefined correlation and abstains (`NA`). For an integer constant the
deviation sum-of-squares cancels to exactly 0 and the abstention falls out of
the non-finite check; but for a **non-integer** constant (POMP / rescaled /
averaged scores are documented numeric input) floating cancellation lands the
term a few ulp on *either* side of zero. Tiny-negative sent `sqrt()` to `NaN`
with a leaked base-R locale-dependent warning (breaking the cli-only typed-
condition contract); tiny-positive leaked a spurious finite score — ~1e-7 for
person-total, ~1.0 (a perfect consistency score!) for psychsyn — that wrongly
entered the percentile pool and the flag denominator. Both kernels now (a)
clamp the variance terms at zero under the `sqrt` (silencing the warning) and
(b) detect a constant row/side **exactly** — masked `min == max` over the
answered values — and force `NA`, independent of which way the cancellation
fell. No tolerance is involved (exact float equality only), so no genuinely
varying respondent can be swept into abstention, and the cross-package parity
is unaffected (`careless` returns `NA` for these rows too, via `cor()`'s
zero-variance `NA`).

## Reverse-keying: the declared range is cross-checked against the data

`apply_split_half_keying()` validates more than the metadata's internal
type-consistency: before reflecting, it checks that every reverse-keyed
column's **observed** responses lie within the declared
`[min, min + categories - 1]`. A type-valid but wrong declaration — the classic
case is 0-based data declared `categories = 5` with the default `min = 1`,
reflecting `0 -> 6` / `4 -> 2` — previously produced off-scale reflected values
that silently corrupted the consistency score (flipping a substantial share of
flags with no signal to the user). The person-fit bridges already catch the
identical mistake in `personfit_zero_base()`; this gives the split-half family
(and Ht's reverse items) the equivalent typed `cier_error_input`, naming the
offending items. Only reverse-keyed columns are checked (forward items are
never reflected, and their `categories` is never read), and an all-NA reverse
column has no observed range to violate.

## cier_screen: a transparent flag-table combiner, no single label

`cier_screen()` runs the registry's screenable indices over one dataset and
returns a `cier_screen` — the per-index flag table (`$flags`), the per-construct
collapsed votes (`$votes`), and the cross-index agreement diagnostic
(`$agreement`). It is an orchestrator, **not a new statistic**: every `cier_index`
it returns is byte-identical to calling that index directly (the trust model is
internal parity at tolerance 0, since no external package computes this combined
screen). It deliberately produces **no single careless/not label** — the dropped
consensus/learned-combiner machinery is out of v0 scope; the screen reports the
count of flagged constructs and the agreement, and the researcher thresholds.

Three decisions (Markus, this slice):

- **Selectable indices.** `methods=` chooses which screenable indices run.
  Goldammer et al. (2024) report resampled personal reliability as the strongest
  single indirect indicator, so weaker indices must be off-selectable to avoid
  diluting it; the default runs all.
- **Redundancy collapse via the registry.** Indices sharing a `vote_group`
  collapse to one vote = logical OR of members' flags, an abstaining (`NA`) member
  counting as not flagged. Only even-odd and personal reliability share a group
  (`consistency`); the agreement diagnostic and the construct count run on the
  **collapsed** votes, never the raw per-index flags, so one construct is never
  weighted as two independent votes. The grouping is a registry column (above),
  not a hard-coded map.
- **Structural skips vs propagated errors.** A metadata index when `items` is
  `NULL`, or a `Suggests`-backed index (Gnormed/Ht) when its package is absent, is
  **skipped with a recorded reason** (`$skipped`, items reason taking precedence).
  A genuinely malformed `items` frame is **not** skipped — the index's own typed
  `cier_error_input` propagates, so the wrapper stays the single source of truth
  for the metadata it requires and the user is told to fix the data.

Tuning is a per-index `control` list (named by method id, each entry an argument
list spliced into that index's call), so any index's `fpr` / `alpha` / `cutoff` /
`seed` / `critical_r` / `n_resamples` is reachable without widening the screen's
own surface. The object is the robust list-based shape (not an attribute-laden
data.frame); the flag count and rate are derived on print, never stored.

A small accessor `cier_flagged_cases()` (an S3 generic over `cier_index` and
`cier_screen`) turns the agreement threshold into actionable respondent row
indices: `cier_flagged_cases(screen, min_votes = k)` returns `which(rowSums(votes)
>= k)` on the **collapsed** votes, so it equals the screen's "flagged by >= k
votes" count and never double-counts the consistency construct. It returns
positions only — no label, no exclusion — keeping the researcher in control of the
threshold.

## Continuous integration: the suite is the gate

CI (GitHub Actions, `r-lib/actions`) runs on push / pull-request to `main`:

- **`R-CMD-check`** on a lean cross-OS matrix — Ubuntu, Windows, macOS, all R
  `release`, `--no-manual`. This validates build / install / examples / vignette
  across platforms. The full devel/oldrel matrix is deliberately not run for this
  pure-R package; re-add it if a release needs the wider surface.
- **Gates live in the test suite, not in bespoke CI jobs.** The
  no-process-vocabulary, roxygen-up-to-date, and references-DOI guards are
  `testthat` tests. They inspect the development sources, which are absent under
  `R CMD check` (it tests the installed package), so a dedicated Ubuntu
  `tests-and-guards` job runs `devtools::test()` against the checked-out source
  with `NOT_CRAN=true` (so the `skip_on_cran()` guards and snapshot tests run) and
  `roxygen2` pinned to the recorded `RoxygenNote` (so the roxygen guard executes
  instead of self-skipping on a version mismatch).
- **Coverage** is reported to Codecov via `covr`; the floor is the 75% in the
  acceptance gate (project status, patch informational). Needs a `CODECOV_TOKEN`
  secret.
- **No automated README-render gate.** The README embeds a live `cier_screen()`
  chunk whose rendered output depends on the `PerFit` / `mokken` / `pandoc`
  versions, so a bytewise render-compare guard would false-fail on dependency
  drift (and rendering it in-process can destabilise the suite); README.md is
  regenerated from README.Rmd by the maintainer and checked at hand-off instead.

## Slow-test tier: `skip_if_slow()` opt-in

A `skip_if_slow()` test helper gates slow tests on the `CIER_SLOW_TESTS`
environment variable (default off): a normal local or CRAN run skips them; CI
opts in by setting `CIER_SLOW_TESTS=true`. Pair it with `skip_on_cran()` so the
tier is enforced in both directions. For the pure-R v0 battery nothing is yet
slow enough to tag (the `PerFit` / `mokken` cross-package parity already gates on
`skip_if_not_installed()` and runs sub-second), so the helper ships as the
documented convention with no test tagged; tag the first genuinely slow path
(e.g. a large-n recovery sweep) when it lands.

## Published-results reproduction vignette

The `published-results` vignette reproduces, index by index, the results of
four published studies (Bruhlmann et al. 2020; the two Goldammer et al. 2024
papers; Schroeders et al. 2022). Its binding choices:

- **Precomputed.** The shipped `vignettes/published-results.Rmd` is generated
  by the maintainers from `vignettes/published-results.Rmd.orig` (committed,
  `.Rbuildignore`d, knitted via `data-raw/knit-published-results.R` in a fresh
  `callr` session). Building the package evaluates nothing: no network, no
  `PerFit`/`mokken`/`mice` at build time. The code readers see in the vignette
  is the code that produced its output.
- **External data stays fetch-only.** The Goldammer polybox archives and the
  Schroeders OSF files carry no redistribution licence; they are downloaded
  into `tools::R_user_dir("cier", "cache")` at knit time and never bundled.
  The slow-tier tests pin each artifact's md5 and skip with a reason when an
  upstream file changes; the Bruhlmann cells need only the bundled
  `bfi_careless`.
- **The committed table is a frozen contract.** The knit writes
  `inst/extdata/published-results.csv`; the fast tier freezes the transcribed
  paper values (digest plus verbatim spot checks against the cited source
  tables) and recomputes every verdict from the stored values; the slow tier
  re-fetches and recomputes every cier cell at 5e-4 under the recorded seeds.
  Wording is honest by construction: a cell is called "matches" only at the
  paper's two-decimal precision (flagged counts: exactly), "close" within
  0.03 (counts: within 2), and anything beyond is "differs" with a mandatory
  written explanation and membership in a capped allowed list (see
  `tests/reference/TOLERANCES.md`, "Published-results reproduction").
- **Pre-registered divergences** (each capped and explained in place): the
  battery paper's Mahalanobis distance (whole-sample centroid vs the authors'
  careful-respondent reference), Ht on the within-subject contrasts,
  whole-sample synonym pairing on the two-wave studies, the careful-reference
  pairing companions, the RPR paper's SEN95 (their tie-corrected ROC
  regression vs the plain empirical cutoff), and the Bruhlmann resampled
  reliability / person-total conventions.
- **Two construction facts the reproduction settled** (now load-bearing in
  the vignette and its tests): the battery paper's resampled personal
  reliability splits within the 15 BFI-2 facets, not the 5 domains; and the
  RPR paper's "conventional PR" is the even-odd consistency, i.e. it maps to
  `cier_even_odd()`, where cier's value correlates -1.0 with the authors'
  stored column.
- **Multiple imputation was evaluated and not adopted.** Imputing abstained
  index values with `mice` moves the affected battery cells by at most a few
  thousandths and pulls the pairing cells away from the paper; the vignette
  documents the check, and cier's abstain-rather-than-impute behaviour
  stands. `mice` is needed only to re-knit, so it stays out of `DESCRIPTION`.

## Items schema: `min` / `max`, not `min` / `categories`

The per-item metadata frame now declares the response range directly — `max`
(the largest response option) replacing `categories` (the count of options) —
after applying the package to a real mixed-format survey showed the count
parameterisation forces users to maintain two interlocking columns and do the
`max = min + categories - 1` arithmetic the package can do itself (Markus,
2026-06-10; pre-release, clean break, no alias or deprecation shim). The two
carry identical information; everything derivable is derived: the reverse-key
reflection is `(min + max) - x` (the classic `(max + 1) - x` at the default
base `min = 1`), the keying cross-check reads `[min, max]`, and the person-fit
bridge's category count is `Ncat = max - min + 1`. The validity bound is
`max >= min + 1` (at least two response options) — deliberately **not** an
absolute `max >= 2`, so a 0/1 item is the smallest valid scale. `categories`
was never a homogeneity requirement: per-item heterogeneous ranges were and
are supported everywhere except the PerFit bridge (next entry).

## Gnormed on mixed response formats: a backend limit, not bad metadata

`PerFit`'s polytomous statistics score one `Ncat` across all items (the PerFit
manual: "The number of answer options, Ncat, is the same for all items"; cf.
Niessen et al. 2016, who recoded their battery to a uniform 4-point scale for
the same reason). cier expresses that contract as a **span homogeneity** check:
`max - min` must be constant — items may differ in base, so `1..5` and `0..4`
items (both five options) score together after per-item zero-basing. A
heterogeneous span on otherwise-valid metadata is *accurate* metadata for
genuinely mixed-format data, not a malformed frame, so the abort carries the
`cier_error_backend_limit` subclass and `cier_screen()` records Gnormed as
skipped-with-reason — the same line mokken's 10-category ceiling already drew —
instead of aborting the whole battery (Markus, 2026-06-10; this deliberately
flips the original "heterogeneous categories propagate" screen pin). Per-item
validity (absent / NA / fractional / non-finite / below-bound `max`) is checked
**before** homogeneity and stays a plain `cier_error_input`, so a genuine
metadata defect still propagates through the screen.

## Mixed response formats: documented caveats, no automatic remediation

The screen's documentation now lists the battery and adds a mixed-format
section; the guidance makes only sourced claims: Gnormed's single-`Ncat`
backend limit (PerFit manual; Niessen et al. 2016), Curran's (2016) "same
response format" condition for even-odd consistency, mokken's since-3.0.3
warning on mixed category counts, Gottfried et al.'s (2022) same-answer-scales
rule for autocorrelation (quoted in `cier_autocorrelation()`'s help), and the
mechanical fact that a fixed-position straightliner stops triggering the
zero-variance abstention once scale ranges differ between blocks (pair the
consistency indices with longstring / IRV). Published validation of the
battery is otherwise uniformly-formatted, and the docs say so rather than
extrapolate.

Two remediations were evaluated on simulated mixed-format data (a 76-item,
4/5/6/7-option design with uniform-random and fixed-position-straightliner
contamination, 8 replications) and **not adopted** (Markus, 2026-06-10):
splitting the battery into homogeneous-format subsets and combining per-subset
votes *degraded* the consistency family (even-odd AUC 0.93 -> 0.81, resampled
personal reliability 0.98 -> 0.90, IRV inverted) — these indices draw their
power from correlating across many scales — and per-item POMP rescaling, which
recovered the (small, ~0.01-0.02 AUC) heterogeneity dent, has no published
precedent in the careless-responding literature, so it stays out of the
documentation. The simulation itself: heterogeneity cost the consistency
family little (even-odd 0.956 -> 0.933, RPR 0.990 -> 0.975, person-total
0.997 -> 0.989 against random responders), left Mahalanobis and the synonym
indices exactly invariant, and Gnormed-per-subset was the one case where
splitting helps (it cannot run otherwise) — left as a documented user-side
option, not automated.
