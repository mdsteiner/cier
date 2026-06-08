# Architecture decisions

The durable record of binding design decisions for cier, in plain terms. It is
deliberately short: one entry per genuinely binding choice — an object schema, a
default cutoff, a numerical tolerance, a cross-package parity contract, or a
scope boundary — not one per feature. Everything else belongs in a commit
message. New entries are added as such decisions are made during the build.

## Scope boundary

cier ships the indirect, response-pattern C/IER indices evaluated by Goldammer
et al. (2024): personal reliability (including the resampled variant),
psychometric synonyms and antonyms, Mahalanobis distance, person-total
correlation, the nonparametric person-fit statistics Gnormed and Ht, and the
classic longstring and intra-individual response variability indices. Adding any
further method family (timing, IRT person-fit, model-based, machine-learning), a
learned combiner, a simulation engine, a specification-curve tool, or a report
generator requires a new entry here and explicit sign-off before code.

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

## Architecture: function-first indices

Each index is a documented function on a response matrix (a data.frame or
tibble is accepted and coerced internally, so users need not call
`as.matrix()`), not a method behind an input-object pipeline. The four indices
that need item metadata (even-odd, personal reliability, Gnormed, Ht) take a
single optional `items` data.frame with columns `scale`, `reverse_keyed`, and
`categories` (one row per item); the other six need only the responses. An index
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
kernel runs. The cutoff dispatch is then a plain conditional inline in the
wrapper: **a literal `cutoff` is used verbatim** (already validated); otherwise
the rate-based default is resolved through the one resolver,
`resolve_cutoff(method = <registry method>, …)`. `resolve_cutoff()` and its
private helpers are internal and do **no** input re-checking (their `method` /
`direction` come from the registry and the rate/literal are wrapper-validated) —
they only do the math and signal the runtime percentile abstention
(`NA_real_` + `cier_warning_insufficient_items`).

This keeps the one-cutoff-path rule — every rate-based **default** still resolves
through the single `resolve_cutoff()` — while removing the earlier
`resolve_index_cutoff()` indirection (a shared helper that took a per-index
`rate_fn` closure). Inlining the two-line dispatch in each wrapper is shorter and
more legible than the closure it replaced, and putting all validation in the
public function gives the earliest possible failure. The cost is that a future
wrapper which forgets to validate would pass bad input silently rather than get a
typed error; the mitigation is the wrapper-validation convention plus per-wrapper
input-error tests.

## Method-properties registry schema

The registry (`inst/extdata/method-properties.csv`) is the single source of truth
for cutoff defaults, flag direction, backend, and screen membership. Its columns
are `method`, `family`, `paper_year`, `paper_citation_key`, `doi`,
`default_cutoff_method`, `default_cutoff_value`, `flag_direction`,
`companion_methods`, `backend`, `screenable`, and `notes`. The boolean capability
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

## Mahalanobis degenerate covariance: warn and abstain

The Mahalanobis index needs a sample covariance to compute any distance. When
none can be estimated — fewer than two respondents carry data, or the covariance
matrix is not invertible (more items than respondents or perfectly collinear
items, so `solve()` fails; or a pairwise covariance with undefined entries — two
items never co-answered, or an item answered by nobody — which yields a non-finite
inverse rather than a `solve()` error) — every respondent's `value` is `NA` and
no one is flagged, and
the wrapper raises a typed `cier_warning_singular_covariance` whose structured
`data$reason` names the cause (`"insufficient_responses"` or
`"singular_covariance"`). A silent all-`NA` result was rejected: a singular
covariance is a substantive, actionable analyst event (a collinear or
over-wide item set), so it earns its own condition class rather than being folded
into the cutoff layer's `cier_warning_insufficient_items`. This is the only
wholesale-abstention case that warns; an individual all-`NA` respondent (alongside
others who answered) is ordinary per-row abstention and stays silent, as in every
other index. The kernel itself stays pure — it returns a status code, and the
wrapper raises the condition.

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
