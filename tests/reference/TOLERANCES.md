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
