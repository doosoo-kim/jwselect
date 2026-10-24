# jwselect

**Preperiod selection and minimum-distance estimation in Stata**  
Doosoo Kim ([doosoonet@gmail.com](mailto:doosoonet@gmail.com))  
Version 0.4.0 — 23 October 2026

Release archive finalized 24 October 2026. Tested on Stata/MP 17 for Windows;
the estimation code is unchanged from the tested `0.4.0-joint1` build.

**[Quick user guide (PDF)](docs/jwselect_quick_user_guide.pdf)** · [Full user guide (PDF)](docs/jwselect_user_guide.pdf)

`jwselect` selects optional preperiod restrictions and estimates cohort and aggregate treatment effects using cohort-time estimates and their full joint covariance. It works after a supported `jwdid` fit or with separately supplied coefficient and covariance matrices.

The command provides:

- **Manual selection:** use a prespecified set of preperiod restrictions, including baseline-only estimation.
- **Wald–BIC selection:** choose restrictions within each cohort (the default), or jointly across cohorts using the full covariance.
- **Joint pretrend tests:** test all eligible optional restrictions or a prespecified subset, with equivalent Wald and MD J statistics.
- **Focused selection:** balance estimated bias and variance for a cohort effect or aggregate, with numerical reselection inference.

Design options include cohort or common baselines, anticipation, treatment-effect horizons, and fixed aggregation weights. Start with the quick guide; the full guide and Stata's `help jwselect` provide the option reference and advanced workflows.

## Install

The main command requires **Stata 16.1 or later**. Prepared worker/collector commands and the release-test driver require **Stata 17 or later**. Python is only needed for developer checks, not for estimation.

1. On this repository's main page, choose **Code → Download ZIP** and extract the archive.
2. Locate the extracted folder containing `jwselect.pkg`.
3. In Stata, replace the example path below with that folder:

```stata
net install jwselect, from("C:/path/to/jwselect-main") replace
which jwselect
findfile _jwselect_mata.do
help jwselect
```

After updating an installation, save your work and restart Stata before using the new version.

For the direct `jwdid` workflow, install its dependencies as well:

```stata
ssc install require, replace
ssc install ftools, replace
ssc install reghdfe, replace
ssc install hdfe, replace
net install jwdid, from("https://raw.githubusercontent.com/friosavila/stpackages/main") replace
```

The matrix interface does not require `jwdid`. Upstream packages are installed separately and are not bundled here.

## Quick start

Replace the variable names with your own. Code never-treated units as `first_treat = 0`. The direct adapter supports a balanced, regularly spaced panel, unweighted linear estimation, unit-level clustering, and numeric covariates that are constant within each unit.

Run the following together from the Do-file Editor:

```stata
jwdid y x1 x2, ivar(id) tvar(year) gvar(first_treat) ///
    never method(regress) cluster(id)
estimates store JW_SOURCE

* Default: select within each cohort, then estimate by joint MD
jwselect, select(wbic) candidates(nested) postwindow(0 2)
estat selection
estat design
estimates store JW_WBIC

* Compare with the maintained-baseline estimate
estimates restore JW_SOURCE
jwselect, select(manual) postwindow(0 2)
estimates store JW_BASE
```

**Restore `JW_SOURCE` before every new direct `jwselect` specification.** Each run replaces Stata's active estimation results. A bare `jwselect` replays the current results.

Here, `postwindow(0 2)` averages the first three treatment periods within each cohort and requires those dates to exist for every cohort. The example uses the default nested candidate family. `candidates(all)` searches every local subset and can become large.

## Joint selection and pretrend testing

To select one restriction configuration jointly across cohorts:

```stata
estimates restore JW_SOURCE
jwselect, select(wbic) wbicscope(joint) ///
    candidates(nested) postwindow(0 2)
estat selection
estat scores
estimates store JW_JOINT

* Joint test of all eligible optional preperiod restrictions
estat pretrend
return list

* Diagnostic for the restrictions selected above
estat pretrend, selected
```

`wbicscope(joint)` minimizes **Wald − r × log(N)** over joint configurations, where `r` is the number of optional restrictions and `N` is the number of independent sampled units. It uses the full covariance across cohorts and dates, then performs joint MD estimation. The unchanged default, `wbicscope(cohort)`, selects separately within each cohort before joint MD estimation. The two selection rules can choose different restrictions.

By default, the joint family is the Cartesian product of the local families. With `K_g` optional dates per cohort, it contains `product(K_g + 1)` nested configurations or `2^(sum K_g)` all-subset configurations. `maxjoint(4096)` is a hard cap: exceeding it causes an error, with no approximation or truncation. A prespecified `jointmatrix(JM)` can supply a smaller joint family; see the guides for its mask format.

`estat pretrend` reports the equivalent joint Wald and MD J values, degrees of freedom, and a chi-square reference p-value. The J value is evaluated through the Wald identity under the compatible covariance-weighted MD model; it is not an independently refitted statistic. This tests optional restrictions, never the maintained baseline. For a prespecified subset, use `estat pretrend, prematrix(P)`, where `P` lists cohort and calendar-time pairs.

**A test of an outcome-selected set is a diagnostic:** its ordinary chi-square p-value is not automatically valid after selection. After focused selection, explicitly identify the output, for example `estat pretrend, selected output(1)`. `output()` is allowed only with `selected`; `selected` and `prematrix()` cannot be combined. With no restrictions, the command reports zero degrees of freedom and no p-value.

A self-contained two-cohort example is provided in [`examples/joint_demo.do`](examples/joint_demo.do).

## Focused selection

For cohort-focused selection and numerical confidence intervals:

```stata
estimates restore JW_SOURCE
jwselect, select(focused) focus(cohort) ///
    candidates(all) postwindow(0 2)
estat inference
```

Read the guide's **Search budget** discussion before interpreting focused results. Inspect the reported diagnostics and use the supplied intervals and p-values; `se_fixed` is a fixed-selection diagnostic. Budget exhaustion returns the best points found and calls for numerical stability checks. Search convergence does not certify global extrema.

## Matrix input

For separately validated cohort-time estimates and their joint covariance:

```stata
jwselect, bmat(B) vmat(V) cells(C) nunits(500) ///
    gap(1) select(wbic) candidates(all) postwindow(0 2)
```

`B` is the estimate vector; `V` is its full sampling covariance; `C` lists matching **cohort and calendar time** pairs. Replace `500` with the number of independent sampled units. The guide explains normalization cells, support, and covariance requirements. A self-contained matrix example is available in [`examples/matrix_demo.do`](examples/matrix_demo.do).

For robust pretrend inference, supply the appropriate heteroskedasticity- or
cluster-robust **full joint covariance**, including cross-cohort and cross-period
terms. The matrix interface uses the supplied covariance; it cannot infer or
repair its statistical construction. The direct adapter requires unit-level
clustering and supports its documented conditional/unconditional conventions.

## Citation and contact

When using the method, please cite the accompanying working paper:

Kim, Doosoo. 2026. *From Pretrend Tests to Preperiod Selection in Difference-in-Differences*. Working paper.

Please also report the software version used. Software citation metadata is provided in [`CITATION.cff`](CITATION.cff).

For questions or bug reports, use this repository's **Issues** tab or email [doosoonet@gmail.com](mailto:doosoonet@gmail.com). Include your Stata version, `which jwselect` output, the command, and a small reproducible example where possible.

## Development and checks

The maintained numerical sources are in [`mata/`](mata/). After editing them, run `python build.py` to regenerate `_jwselect_mata.do`. Both the modules and the generated loader belong in the repository. Editable LaTeX sources accompany both guide PDFs.

From the repository root, run the Python source and independent numerical checks:

```sh
python -m pip install -r requirements-dev.txt
python build.py --check
python tests/static_check.py
python tests/numerical_check.py
python tests/reselect_check.py
python tests/release030_reference.py
python tests/core_optimization_check.py
python tests/joint040_reference.py
```

These checks cover source structure and independent numerical calculations; they do not execute Stata or establish statistical coverage. See [CHANGELOG.md](CHANGELOG.md) for the release scope and validation limits.

For the Stata release checks, save your work, start a **fresh Stata 17 or later session**, change to the repository root, and run:

```stata
do run_release_tests.do
```

The driver clears data and estimates. Look for the final `jwselect 0.4.0 release tests: PASS` message. The separate `tests/integration_jwdid.do` checks the upstream adapter when its dependencies are installed; it does not install them.

The corrected additional tests are included in this release:

```stata
do run_extra_tests.do
```

They write a complete log and a four-suite summary in a new
`tests_extra/results_<timestamp>` folder. See
[`README_EXTRA_TESTS.txt`](README_EXTRA_TESTS.txt) for coverage and instructions.

### Validation record — 24 October 2026

Release-test components passed across the supplied Stata/MP 17 Windows logs,
including worker/collector, caching/streaming, restart/cleanup, 648 complete-fit
optimization comparisons, and the joint-selection/pretrend public API checks.
The separate upstream integration run was reported as PASS by the author.
All four additional suites also passed across the subsequent runs:

| Additional suite | Coverage | Evidence |
| --- | --- | --- |
| Actual `jwdid` fits | 17 joint-selection cases | Author-reported final PASS |
| Design/invariance | 35 cases | Supplied Stata log |
| Independent numerical reference | 24 designs, 48 fits, 960 joint candidate comparisons, 96 pretrend tests | Supplied Stata log |
| Saved results/postestimation | Four scenarios | Supplied Stata log |

The numerical reference independently minimizes a covariance-weighted MD
criterion and compares its J value with the Wald quadratic form. Synthetic
covariances check the numerical identity; the actual-fit tests exercise the
supported robust covariance workflow. These checks do not establish statistical
coverage, certify global numerical optimization, or verify every Stata version
and operating system. Actual simultaneous multi-session filesystem behavior
has its separate concurrency test. No fresh concurrent run or clean-install
run is claimed by this validation record.

The test corrections affected test logic, fixture loading, and diagnostics;
they required no change to estimation or postestimation code. The final archive
was checked for manifest completeness and byte equality of its production files
with the tested build.

### Check the installed copy

After installing with `net install` above, save your work and restart Stata.
Change to an existing writable folder outside the extracted package. Then run
the check with the **absolute source-folder path** as its argument, replacing
both example paths:

```stata
cd "C:/work/install-check"
do "C:/downloads/jwselect-0.4.0/run_installed_check.do" "C:/downloads/jwselect-0.4.0"
```

This clears data and estimates. It does not install packages or change
`adopath`. It checks that all installed command/help files resolve outside
the source folder, verifies the build, and runs a small self-contained
matrix example. It requires no `jwdid`. Look for
`jwselect 0.4.0 installed-package check: PASS`; the timestamped log is written
in the working directory. This new check was reviewed during final packaging
but has not been run in Stata here.
