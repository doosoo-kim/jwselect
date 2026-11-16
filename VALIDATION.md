# Release validation — 16 November 2026

Version 0.9.1 adds a postestimation comparison and updates the documentation. The existing selection, estimation, and inference algorithms are retained. The following checks were executed in fresh Stata/MP 17 sessions using this release; their source files and logs are included.

The subsequent November 16, 2026 terminology updates change documentation only and retain version 0.9.1 and build 0.9.1-20261116. They cover AMSE/MSE/RMSE wording, standard errors, estimator and procedure terminology, and explanations of selection, standardization, and numerical error controls. The validation runs below precede these wording updates; they are not new simulation runs. Estimation and inference code and stored diagnostic identifiers remain unchanged.

For the latest terminology update, SHA-256 comparisons verified that all 41 code, example, and test source files are byte-for-byte identical to the preceding package. Guide equations and executable examples are unchanged. Both guides were recompiled without unresolved references or overfull boxes, and every affected page passed visual inspection. The refreshed ZIP and checksum manifest were checked separately.

| Check | Evidence | Result |
|---|---|---|
| Existing MD covariance, WBIC, focused ties, and full index blocks | `tests/regression_090.do`, `.log` | Passed |
| Current build identity, incompatible modules, and prepared-job identity | `tests/build_identity.do`, `.log` | Passed |
| Matched baseline/own-cohort/joint estimates against an analytic two-cohort example; full covariance; empty and absent optional sets; manual and both WBIC scopes; all focused modes with distinct selected masks; unequal within-cohort weights; nonlinear blocks; state preservation; unsupported scalar-source rejection | `tests/comparison_091.do`, `.log` | Passed |
| Actual untrimmed `csdid` integration, trimming guard, source-state preservation, and matrix-route agreement | `tests/dr_untrimmed.do`, `.log` | Passed |
| Logit/Poisson full index blocks, raw/centered covariates, source covariance, and default dispatch | `tests/nonlinear_index_release.do`, `.log` | Passed |
| RCS transported inputs and scalar selection | `tests_rcs.do`, `.log` | Passed |
| RCS outcome-correct and propensity-correct branches; independent complete influence covariance checks under pooled and stratified sampling | `tests_rcs_branches.do`, `.log` | Passed |
| Scalar signatures and source design validation | `tests_rcs_schema.do`, `.log` | Passed |
| SA population-share covariance and downstream source workflow | `tests/sa_population.do`, `.log` | Passed |
| Separate `net install`, manifest/helper discovery, exact build identity, and all three runnable examples | `tests/isolated_install_091.do`, `.log` | Passed |

The isolated installation exercises `examples/comparison_demo.do`, `examples/matrix_demo.do` (including centered bootstrap and McCloskey inference), and `examples/rcs_composition_demo.do`. It uses a separate ado directory. Current version/build checks require `0.9.1` and `0.9.1-20261116`.

Both PDF guides compile without unresolved references or overfull boxes: the quick guide has 14 pages and the full guide has 43 pages. The distribution includes their editable LaTeX sources. Release assembly verifies every installation-manifest entry, supplies SHA-256 checksums, and checks the ZIP contents and integrity.

Run tests only in a separate Stata session without unsaved work; most clear the active data. Ordinary tests run from the extracted package root. The isolated installation test takes an absolute package path and an existing empty temporary directory. Upstream integration tests require their corresponding external estimator dependencies; no such test was skipped in the recorded run.

These checks verify implementation and packaging, not empirical identification, finite-sample efficiency dominance, or global numerical optimization. The complete historical Monte Carlo collection and empirical archive were not regenerated for this maintenance release. Their original numerical files remain unchanged. The full guide describes the remaining inference and numerical qualifications.
