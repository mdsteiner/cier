# Numerical tolerances

The per-quantity tolerances the test suite holds the implementation to. These
are **binding**: loosening one to make a test pass is an anti-pattern. Either
fix the implementation or -- only with a recorded entry in `ADR.md` -- replace
the reference comparison.

The cutoff + flagging layer (slice 1) has no cross-package partner: `careless`,
`psych`, `PerFit`, and `mokken` neither resolve flagging cutoffs nor compute a
Poisson-binomial. Its references are therefore the base-R primitives it composes
(`stats::quantile`, `stats::qchisq`) and an independent `2^m`-enumeration oracle
for the agreement diagnostic. The per-index cross-package rows arrive with their
slices (2-11).

| Quantity | Target tolerance | Reference |
|---|---:|---|
| Percentile cutoff vs. `stats::quantile(x, p, type = 7)` | 0 (identical) | base R |
| Chi-square cutoff vs. `stats::qchisq(1 - alpha, df)` | 0 (identical) | base R |
| Fixed cutoff (verbatim passthrough of the resolved value) | 0 (identical) | -- |
| Agreement diagnostic Poisson-binomial tail vs. enumeration oracle | 1e-12 | `ref-poisson-binomial-enumeration.R` |

## Per-index parity (slices 2-11)

One row per statistic, added when its slice lands. Each index has an independent
definition oracle (re-derives the statistic by hand, never calls the production
kernel) and, where a partner package computes the same quantity, a cross-package
parity check.

| Index | Quantity | Target tolerance | Reference |
|---|---|---:|---|
| longstring | max run length vs. independent oracle (`ref_longstring$longest`) | 0 (exact) | `ref-longstring-johnson2005.R` |
| longstring | vs. `careless::longstring()` on complete data | 0 (bytewise) | `careless` (1.2.2) |
| irv | per-row sample SD vs. independent oracle (`ref_irv`) | 1e-10 | `ref-irv-marjanovic-2015.R` |
| irv | vs. `careless::irv()` (incl. `na.rm` rows) | 1e-10 | `careless` (1.2.2) |
| mahalanobis | per-row D² vs. independent oracle (`ref_mahalanobis`) | 1e-10 | `ref-mahalanobis-curran-2016.R` |
| mahalanobis | exact hand-computed 2-column fixture (`c(2, 0, 2, 2)`) | 1e-12 | worked by hand |
| mahalanobis | vs. `careless::mahad(flag = FALSE)` (incl. NA rows) | 1e-10 | `careless` (1.2.2) |
| mahalanobis | vs. `psych::outlier()` D² (incl. NA rows) | 1e-10 | `psych` (2.6.5) |
| person-total | per-row correlation vs. independent oracle (`ref_person_total`) | 1e-12 | `ref-person-total-donlon-fischer-1968.R` |
| person-total | exact hand-computed fixture (`c(1, 1, -1, 1)`, linear item-mean profile) | 1e-12 | worked by hand |
| person-total | vs. `PerFit::r.pbis()$PFscores` on complete data | 1e-4 | `PerFit` (4-dp output rounding) |
| even-odd | per-row `-SB(r)` vs. independent oracle (`ref_even_odd`) | 1e-12 | `ref-evenodd-curran-2016.R` |
| even-odd | analytic fixtures (consistent `-1`, inverse `+1`) | 1e-12 | worked by hand |
| even-odd | honouring `reverse_keyed` == independently pre-scored input | 1e-12 | property |
| even-odd | vs. `careless::evenodd(factors=)` on no-reverse-key data | 0 (bytewise) | `careless` (1.2.2) |
| personal-reliability (PR) | per-row `-SB(r)`, first/second-half vs independent oracle (`ref_pr`) | 1e-12 | `ref-pr-jackson-1976.R` |
| personal-reliability (PR) | analytic fixtures (consistent `-1`, inverse `+1`) | 1e-12 | worked by hand |
| personal-reliability | honouring `reverse_keyed` == independently pre-scored input | 1e-12 | property |
| personal-reliability (RPR) | per-row mean over 25 seeded random split-halves vs independent oracle (`ref_rpr`) | 1e-10 | `ref-rpr-goldammer-2024.R` |
| psychsyn | per-row stacked-pair correlation vs independent oracle (`ref_psychsyn`) | 1e-12 | `ref-psychsyn-meade-craig-2012.R` |
| psychsyn | orthogonal-contrast hand fixture (pair set `{(2,1),(4,3),(6,5)}`, values `c(NA, 1, 1, 1)`) | 1e-12 | worked by hand |
| psychsyn | vs. `careless::psychsyn(resample_na = FALSE)` on `careless_dataset` | 1e-12 | `careless` (1.2.2) |
| psychsyn (pairing) | `cier_synonym_pairs()` vs. `careless::psychsyn_critval()` (full pairing, both tails) | 1e-12 | `careless` (1.2.2) |
| psychant | per-row stacked-pair correlation vs independent oracle (`ref_psychant`) | 1e-12 | `ref-psychant-meade-craig-2012.R` |
| psychant | orthogonal-contrast hand fixture (antonym pair set `{(2,1),(4,3),(6,5)}`, values `c(NA, -1, -1, -1)`) | 1e-12 | worked by hand |
| psychant | vs. `careless::psychsyn(anto = TRUE, resample_na = FALSE)` on a planted antonym fixture | 1e-12 | `careless` (1.2.2) |
| gnormed | per-row value vs independent oracle `round(ref_personfit_gnormed_poly, 4)` | 0 | `ref-personfit-niessen-2016.R` |
| gnormed | vs `PerFit::Gnormed.poly` direct call (complete, non-keyed, n != p) | 0 (bytewise) | `PerFit` (1.4.7) |
| gnormed | reduces to dichotomous `PerFit::Gnormed` at `Ncat = 2` | 1e-9 (obs 0) | `PerFit` (1.4.7) |
| gnormed | Monte-Carlo null cutoff vs `PerFit::cutoff(fit, Blvl)$Cutoff`, same seed | 0 (bytewise) | `PerFit` (1.4.7) |
| ht | per-row value vs independent oracle `ref_personfit_ht_poly` | 1e-12 (obs ~5e-16) | `ref-personfit-niessen-2016.R` |
| ht | vs hand-built `mokken::coefH(t(zero_base(m)))$Hi` (complete, non-keyed, n != p) | 0 (bytewise) | `mokken` (3.1.2) |
| ht | reduces to dichotomous `PerFit::Ht` at `Ncat = 2` | 1e-4 (obs 4.8e-5) | `PerFit` (4-dp output rounding) |

Personal reliability (PR / RPR) has **no cross-package partner**: `careless`,
`psych`, `PerFit`, and `mokken` implement neither variant. The two independent
paper oracles above (`ref_pr`, `ref_rpr`) are therefore its parity checks. The
RPR oracle additionally coordinates its random-draw order with production so a
fixed-seed run matches bytewise -- a deliberate reproducibility constraint, not
a tautology (the statistic itself is re-derived from scratch).

Psychsyn's kernel scores each respondent with a **vectorised masked-sum** Pearson
over the stacked pair matrices (one `rowSums` pass, not a per-row `cor()` loop).
This is the same correlation as the per-row `cor()` that `careless::psychsyn` and
the independent oracle use, but it sums in a different order, so it matches both
at ~1e-13 — held at **1e-12**, not bytewise. (The earlier per-row-`cor()` kernel
matched `careless` exactly at `0`; it was vectorised for a 4–6× speedup, trading
the bytewise guarantee for scale and uniformity with the `person-total` kernel —
see `ADR.md`, "Psychsyn/psychant kernel: vectorise".) The oracle re-derives the
statistic by an independent path (it pre-filters complete cases before `cor()`),
and the kernel matches it to ~1e-13 — comfortably inside the 1e-12 row above.

Psychant shares that same vectorised kernel via the `pairing = "ant"` tail
(`r < -critical_r` instead of `r > critical_r`), so it carries the **identical
1e-12** against its own independent oracle (`ref_psychant`, which reuses the
psychsyn row correlation with the negated threshold) and against `careless`.
`careless::psychant()` does not surface `resample_na`, so the deterministic
parity comparison calls the underlying `psychsyn(anto = TRUE,
resample_na = FALSE)` directly. `careless_dataset` carries **no** antonym pairs
at `r < -0.60`, so the parity uses a constructed fixture with planted negative
structure (items mapped to a 1–5 Likert) plus injected `NA`s to pin NA
agreement.

Gnormed's production scorer **is** `PerFit::Gnormed.poly` (single-kernel rule), so
its genuine independent check is the closed-form oracle
`ref_personfit_gnormed_poly` (a popularity-rank numerator + max-plus-knapsack
denominator, re-derived from scratch and never calling the production bridge).
`PerFit` rounds its scores to 4 decimals and the oracle is exact, so
`round(oracle, 4)` matches **bytewise** (tolerance 0). The `PerFit::Gnormed.poly`
direct-call row is not a redundant scorer check: the bridge supplies the
preprocessing independently (zero-basing to `0..(Ncat - 1)`, persons-as-rows
orientation, reverse-keying, complete-casing), so an `n != p` fixture makes a
missing-transpose or raw-`1..k`-coding bridge diverge while the shared scorer
holds it bytewise. The cutoff is the PerFit Monte-Carlo null (`PerFit::cutoff`),
randomised but **reproducible under a seed**; the parity re-fits the same
zero-based block and seeds identically immediately before the call, so it matches
bytewise (the statistic is independent; only the RNG stream is coordinated -- the
same white-box reproducibility constraint as RPR, not a tautology).

Ht's production scorer **is** `mokken::coefH(t(z))$Hi` (single-kernel rule), so
its genuine independent check is the closed-form oracle `ref_personfit_ht_poly`
(the Frechet / rearrangement collapse of person scalability, re-derived from
scratch and never calling the production kernel). `coefH` returns full precision
(unlike `PerFit`'s 4-dp rounding), so the oracle holds to ~5e-16 -- recorded at
**1e-12**. The hand-built `mokken::coefH(t(.))` row is not a redundant scorer
check: it pins the bridge's persons-as-rows orientation (the transpose), so an
`n != p` fixture makes a missing-transpose mutant (`coefH(z)`, item scalability)
return the wrong length while the shared scorer holds it **bytewise**. (Global
zero-basing is translation-invariant for the covariance-based Ht, so it does not
need a dedicated parity row -- the oracle exercises the remaining preprocessing.)
The dichotomous-reduction row checks that the polytomous `coefH` path reduces to
the classic Ht at `Ncat = 2`, matching `PerFit::Ht` to that package's 4-dp output
rounding (**1e-4**, observed 4.8e-5). Unlike Gnormed, Ht has **no Monte-Carlo
null** row: no model-conforming null exists for the mokken-backed polytomous
statistic, so its cutoff is the empirical lower-tail percentile (`stats::quantile`,
already covered by the slice-1 cutoff rows).

## How to use this table

- Tests in `tests/testthat/test-cutoff.R` / `test-diagnostics.R` assert these
  with `expect_identical()` (the bytewise rows) or
  `expect_equal(..., tolerance = ...)`. The most safety-critical values are
  inlined at the assertion site so each test is self-documenting.
- A bytewise (`0`) tolerance means any non-zero difference (including
  floating-point noise) is a divergence and must be diagnosed before the test
  passes.

## When a tolerance is missed

1. **Diagnose first.** Is the difference definitional (e.g. a different
   `quantile` `type`, or a flip applied twice)? If so the comparison or the
   implementation is wrong, not the tolerance.
2. **Do not loosen the tolerance** to make the test pass. Fix the implementation
   or -- only with an `ADR.md` entry -- replace the reference comparison.

The agreement-diagnostic row is an **independent definition oracle**: the
enumeration re-derives the Poisson-binomial tail from first principles (sum over
all flag patterns), never by calling the production convolution kernel.
