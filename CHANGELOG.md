# Changelog

## 0.5.0 - 26 October 2026

- Added an optional CS-DR adapter for the csdid 2.x interface, tested against
  version 2.0.1. It reads unaggregated group-time effects and their full
  analytical covariance, reconstructs normalization cells by name, and detects
  changes to the signed estimation data.
- Added nonlinear jwdid support for logit and Poisson through joint outcome-scale
  ATT/placebo margins, with conditional or unconditional covariance. The adapter
  validates convergence and known baseline normalizations before removing
  numerical derivative residue at structural zeros.
- Both new adapters retain the source cohort baseline and infer its anticipation
  convention when omitted. Unsupported rebasing, mismatched anticipation, and
  unsupported source specifications fail explicitly.
- Added source estimator/version/method/model and restriction-scale metadata.
  Updated display provenance and documented the canonical covariance-weighted
  moment adjustment. This is not an upstream nonlinear constrained refit or
  CS-DR re-estimation at alternative base periods.
- Preserved existing linear and matrix estimation behavior. Upstream packages
  remain optional and separately installed; no estimation command downloads or
  updates dependencies.
- Added controlled CS-DR fixtures, real-upstream integration suites, two examples,
  and an optional adapter-test driver. Updated both PDF manuals and help.
- Passed the Stata/MP 17 release gate, all actual adapter integrations, the new
  examples, an isolated local installation check, and six independent Python
  check scripts. See VALIDATION.md and validation/ for the executed record.

Initial support is limited to the documented balanced, unweighted panel designs
with never-treated controls and time-invariant numeric covariates. Numerical
validation does not establish statistical coverage or global optimization.

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
