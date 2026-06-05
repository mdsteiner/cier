# cier

`cier` detects careless / insufficient-effort responding (C/IER) in survey data
using a lean set of **indirect, response-pattern indices** from Goldammer et al.
(2024). S3 throughout, `cli` output, validated by reference oracles. Fast and
auditable: a non-statistician maintainer must be able to trace every number to a
paper or a trusted package.

## Scope (HARD boundary)

v0 ships exactly these indices: personal reliability (PR + resampled RPR),
psychometric synonyms, psychometric antonyms, Mahalanobis distance, person-total
(r_pbis), the nonparametric person-fit statistics Gnormed (`PerFit`) and Ht
(`mokken`), plus the classics longstring and IRV. Nothing else. A new method
family (timing, IRT person-fit, model-based, ML), a learned combiner, simulation,
spec-curve, or a report generator each needs an `ADR.md` entry + Markus's sign-off
**before** any code.

## Anatomy of a method

```
registry row (inst/extdata/method-properties.csv, cited cutoff + direction)
  -> pure kernel (R/kernels-*.R: math on matrices, no I/O, no state)
       -> thin wrapper R/cier-<verb>.R (<= 30 lines: validate -> kernel -> cutoff -> assemble)
            -> S3 cier_index (new_ / validate_ / print / summary)
```

One production implementation per statistic (single-kernel rule). Alternative /
paper re-derivations live only in `tests/reference/`. No statistical logic in
wrappers, formatters, or plotters.

## Build a method oracle-first (and stop for sign-off)

Anything numerically delicate — an index, a cutoff, a tolerance, an ordinal
correlation, an object schema — stops for Markus: write the paper reference
oracle and cite the value, name the parity package + tolerance, show the
recovery / parity delta, and get a go-ahead **before** writing the kernel. Then:
kernel -> property / invariant tests -> thin wrapper -> commit. The reference
value is the spec.

## Cutoffs (no ground truth)

Default = empirical **percentile at a user-set target false-positive rate** (5% =
Goldammer's Sen95), documented as a **ranking** convention, not a calibrated FPR
(a sample percentile always flags p% by construction). Use a proper null where it
is cheap: chi-square for Mahalanobis, the `PerFit` Monte-Carlo null for
Gnormed / Ht. Report a flag-rate / cross-index-agreement diagnostic so
contamination is visible. Never hard-code absolute heuristics ("RPR < .3",
"H >= .30"). The detection signal lives in multi-index agreement; even-odd, PR,
and RPR are one construct — never count them as independent votes.

## Trust model

- `tests/reference/` — each index reproduces its source paper's worked example,
  independent of the production kernel.
- Cross-package parity (`careless`, `psych`, `PerFit`, `mokken`) at a tolerance
  recorded in `tests/reference/TOLERANCES.md` (a binding contract; loosening it
  is a deliberate, recorded act).

## Acceptance gate (paste evidence, never assert)

`roxygen2::roxygenise()` · `devtools::test()` (fast tier green) ·
`devtools::check(document = FALSE)` (no new NOTE) · `lintr::lint_package()`
(0 new) · coverage >= 90% on changed files · report the reference / recovery
numbers. Tier tests: fast unit + oracle / parity (no skip) vs slow
(`skip_on_cran()` + `skip_if_slow()`). Per edit run only the affected file; run
the full fast tier pre-commit; `check()` pre-PR.

## Conventions

- `cli` only (typed `cier_abort` / `cier_warn` / `cier_inform`; `print` /
  `summary` to stdout). roxygen2 markdown on every export.
- **User-facing output is design-first:** propose a literal text mock-up, get
  Markus's approval, then implement and lock with a snapshot test. Never send a
  print method straight to code.
- **No labelled process vocabulary** (`Phase N`, `Step N`, `Slice N`, `Card`,
  `Decision N`, `ticket`, "see <internal doc>") in any **shipped** file — explain
  it plainly in place. A precise grep guard scans the packaged tree (the domain
  term "decision rule" is fine). `dev/` and `archive/` are `.Rbuildignore`d and
  exempt.
- One logical change per commit (Conventional Commits). No scratch files in the
  repo; `ADR.md` + `NEWS.md` are the durable record. Size budget: file <= 500,
  function <= 80, public wrapper <= 30. No `:::`, no `<<-`, no `library()` /
  `require()` in package code.
- Windows environment; use PowerShell for shell commands.

## Where things live

- `dev/restart/` — the design dossier that motivated this build (review,
  research, plan). Dev-only, `.Rbuildignore`d, **not** a user reference; binding
  facts graduate into `ADR.md` / the docs in plain terms.
- `archive/` — the previous exploratory version, kept on disk for porting
  validated kernels and reference oracles. `.gitignore`d + `.Rbuildignore`d.
