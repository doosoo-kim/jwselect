# Version 0.5.0 validation - 26 October 2026

This release was built from the supplied jwselect-0.4.0.zip and tested using
Stata/MP 17 on Windows in separate batch sessions with synthetic data. The source
archive and the user's existing installed jwselect were not overwritten.

| Check | Result | Evidence |
| --- | --- | --- |
| Stata release gate, including CS-DR contract fixtures | PASS | validation/release_tests.log |
| Real csdid 2.0.1 adapter | PASS | validation/adapter_tests.log |
| Existing linear jwdid adapter | PASS | validation/adapter_tests.log |
| Nonlinear jwdid logit and Poisson adapters | PASS | validation/adapter_tests.log |
| New CS-DR and nonlinear examples | PASS | validation/examples.log |
| Isolated local net installation and installed-copy check | PASS | validation/installed_check.log |
| Six Python structural and independent numerical scripts | PASS | validation/python_checks.log |
| PDF guide compilation and visual review | PASS | Editable sources and PDFs in docs/ |

The actual upstream adapter driver ran all its suites; none was skipped.
The release gate includes worker/collector, cache/restart, complete-fit
optimization, joint selection/pretrend, and controlled CS-DR cases. The isolated
installation checked all 15 executable/help files outside the source folder.

## Upstream versions

- csdid 2.0.1, 23 October 2026, from commit
  `f98bcea0691b1670c409bb7e3cc81d60705b4f84` of
  https://github.com/pedrohcgs/csdid-stata. Downloaded into a project-local test
  directory; no global csdid installation was changed.
- Installed jwdid 2.201 and jwdid_estat 2.2.
- Installed reghdfe 6.13.1, hdfe 3.2.9, and ftools 2.50.0.
- Python 3.12 with project-local NumPy 2.5.4 and SciPy 1.18.1 for developer checks.
  Python is not a runtime dependency of jwselect.

## What the adapter tests verify

CS-DR tests verify named mapping of reordered cohort-time cells, reconstruction
of structural-zero reference cells, full covariance retention, baseline target
parity, joint-WBIC parity with explicit matrix input, upstream pretrend Wald
parity, anticipation, unit clustering, signed-data checks, and rejection of
unsupported versions, options, aggregates, missing cells, and bootstrap covariance.

Nonlinear tests verify independent inverse-link predictive contrasts, full
conditional and unconditional covariance parity with upstream margins, the
canonical ATT/placebo vector, independent closed-form MD projection, joint
WBIC, anticipation, convergence/index/method rejection, and restoration of data
and source results. Known structural-zero baseline margins are first validated
within 1e-10, then their numerical-derivative residue is set to exact zero.

The unchanged linear integration suite also passes. No runtime dependency is
required for the matrix interface or the controlled CS-DR fixture suite.

## Scope and limits

The first direct adapters deliberately support the documented unweighted,
balanced panels with observed never-treated controls, unit-level dependence,
and plain numeric covariates constant within unit. CS-DR requires analytical
covariance, universal base periods, and unit-spaced integer time. The nonlinear
adapter supports logit for binary outcomes and Poisson for an exponential mean.

The new adapters adjust baseline ATTs using average outcome-scale placebo
moments and their joint covariance. They do not refit constrained nonlinear
models or recompute CS-DR at alternative bases. The numerical Wald-J identity
belongs to this canonical covariance-weighted MD problem. Identification,
nuisance-model validity, joint limit theory, and post-selection inference
conditions remain substantive requirements. Adjusted nonlinear ATTs are not
clipped to natural bounds.

These tests do not establish finite-sample coverage or certify global numerical
optimization. A fresh simultaneous multi-session concurrency run, other Stata
versions/platforms, and unsupported panel/weight/covariate designs were not tested.

## Reproduction

From a fresh Stata session in the extracted package root, run
`do run_release_tests.do`. With the separately installed supported dependencies,
run `do run_adapter_tests.do`; inspect PASS/SKIP messages. Both drivers clear
in-memory data and estimates. The README documents independent Python checks and
the separately runnable installed-copy check.

SHA256SUMS.txt covers every distributed file other than itself. The companion
archive checksum identifies the final ZIP.
