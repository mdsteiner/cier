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
