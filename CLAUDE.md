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

## Build a method test-first (and stop for sign-off)

Anything numerically delicate — an index, a cutoff, a tolerance, an ordinal
correlation, an object schema — stops for Markus first: state the definition,
cite the paper reference value, name the parity package + tolerance. Then build
in this order:

1. **Write the tests before the implementation** — the independent reference
   oracle (re-derives the statistic by hand, never calls the production code),
   the cross-package parity check at the recorded tolerance, the property /
   invariant tests, and the edge / degenerate cases.
2. **Adversarial test review (mandatory).** Run the `test-adversary` subagent on
   the new tests against the feature spec. It mutation-probes the suite (would a
   plausible-but-wrong implementation survive?). Close every surviving-mutant gap
   it flags **before any implementation exists**.
3. **Implement** the pure kernel, then the `<= 30`-line thin wrapper, until the
   tests pass — no more.
4. **Hand off for review.** Do NOT commit (see Conventions). Report the
   reference / parity / recovery numbers as evidence.

The reference value is the spec.

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
(0 new) · coverage: statistical kernels covered extensively (oracle / parity /
direction / reverse-keying / edges), overall >= 75% — do NOT pad with
defensive-branch tests on plumbing; keep the suite reviewable · report the
reference / recovery numbers. Tier tests: fast unit + oracle / parity (no skip) vs slow
(`skip_on_cran()` + `skip_if_slow()`). Per edit run only the affected file; run
the full fast tier + `check()` before each hand-off for review (no PRs in pre-release).

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
- **Never `git commit` or push.** Markus reviews and commits. Keep each change to
  one logical unit so a review maps to one commit, and suggest a
  Conventional-Commit message for him. `ADR.md` + `NEWS.md` are the durable
  record; no scratch files in the repo.
- **Work directly on `main`** — no feature branches or PRs during this pre-release.
  Hand off the working tree for review; Markus commits on `main`.
- Size budget: file <= 500, function <= 80, public wrapper <= 30; lintr
  `cyclocomp <= 20`, `object_length <= 40`, `line_length <= 100`. No `:::`, no
  `<<-`, no `library()` / `require()` in package code.
- Windows environment; use PowerShell for shell commands.

## Where things live

- **The build plan is `dev/restart/plan.md`** — the ordered sequence to implement,
  with `architecture.md`, `index-specs.md`, and `example-data.md` as references.
  Start there each session.
- `archive/dev/restart/` holds the original rationale (review, research, plan) if
  the "why" is ever needed. Dev-only, `.Rbuildignore`d, **not** a user reference;
  binding facts graduate into `ADR.md` / the docs in plain terms.
- `archive/` — the previous exploratory version, kept on disk for porting
  validated kernels and reference oracles. `.gitignore`d + `.Rbuildignore`d.
