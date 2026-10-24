# Changelog

## 0.4.0 release archive finalized — 24 October 2026

- Preserved the tested `0.4.0-joint1` estimation and postestimation code.
- Recorded completed Stata/MP 17 Windows testing in the README and full guide.
  The earlier release-test components passed across the supplied logs, including
  648 complete-fit comparisons. The author reported PASS for the separate
  upstream integration test.
- Included the corrected additional suites: 35 design/invariance cases; 24
  numerical designs with 48 fits, 960 joint candidate comparisons, and 96
  pretrend tests; four saved-results scenarios; and 17 actual-jwdid joint cases.
  The first three are verified by supplied logs; final PASS for the last was
  reported by the author. Together the four suites passed across these runs.
- Corrected test-only issues involving absent returned matrices, large Stata
  matrix expressions, and data-signature syntax/comparisons. No estimation or
  postestimation fix was needed for these additional runs.
- Described the independent reference as covariance-weighted minimum-distance
  J calculations, compared with Wald using the same compatible covariance.
  Robust inference requires an appropriate full robust joint covariance.
- Added a separately runnable installed-package smoke check. It has been
  reviewed but has not been executed here; the preparation environment has no
  Stata runtime.

The archive's installation manifest is complete, and production-file bytes
match the tested build. Validation is specific to the tested configurations;
it does not establish statistical coverage or global optimization. The record
does not claim a fresh simultaneous multi-session test, a clean installation
run, or runtime verification on Stata 16.1 or other operating systems.

## 0.4.0 — 23 October 2026

- Added `select(wbic) wbicscope(joint)`: joint Wald–BIC selection using the full covariance across cohorts and dates, followed by joint minimum-distance estimation. The default remains `wbicscope(cohort)`.
- Extended `jointmatrix()` to joint WBIC, with the existing Cartesian-product default and explicit `maxjoint()` cap. Candidate limits cause an error; no approximation or truncation is introduced.
- Added `estat pretrend` for joint testing of all eligible optional restrictions, a prespecified `prematrix()`, or a selected output's restriction mask. It returns Wald, the equivalent MD J value, degrees of freedom, and a chi-square reference p-value. Empty tests have zero degrees of freedom and no p-value.
- Added `e(wbic_scope)` and joint WBIC candidate diagnostics. Existing focused selection and numerical inference remain available.
- Revised both PDF guides, editable sources, README, and Stata help. The guides distinguish prespecified tests from selected-set diagnostics and explain joint-search feasibility.

The reported J uses the covariance-compatible Wald identity; it is not independently refitted. A selected-set chi-square p-value is a diagnostic and is not automatically valid after outcome-dependent selection. The maintained baseline is never tested.
