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
for Mahalanobis distance, `alpha` = 0.001), and `fixed` (the longstring count,
resolved as `ceiling(0.5 * ncol)` in that wrapper and passed through verbatim).

The percentile method applies the direction flip **exactly once**: `upper` flags
the high tail (the `1 - fpr` quantile), `lower` flags the low tail (the `fpr`
quantile). The registry stores each method's *literal directional quantile*
(`cier_irv` 0.05/lower, `cier_even_odd` 0.95/upper); that stored value is the
documented result of `(fpr = 0.05, the row's flag direction)`, never an input
that is flipped again. This is the deliberate fix for the v1 footgun where a
stored directional quantile was re-flipped via `1 - p` and landed on the wrong
tail. A single `fpr` knob therefore means the same target tail mass for every
index, which is what lets the analyst sweep `fpr` across {0.01, 0.05, 0.10}. The
flag comparator (`>=` / `<=`) is applied separately, by `apply_flag()`.

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
returns a light `cier_index` — a `data.frame(value, flagged)` carrying its
method, cutoff, and flag direction as attributes, with a print method.
`cier_screen()` is a thin orchestrator that runs the selected indices and
returns a per-respondent flag table plus the count of flagged *constructs*
(even-odd and personal reliability collapse to one vote). There is deliberately
no `cier_data` / `cier_items` input class and no separate constructor/validator
layer: the matrix-plus-metadata convention matches the parity packages
(`careless`, `PerFit`, `mokken`), keeps the surface small, and pushes validation
to small per-function input checks. Cross-cutting concerns stay shared and are
the only retained foundation: typed conditions, the method-properties registry,
and one cutoff resolver.

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
