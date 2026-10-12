# jwselect

**Preperiod selection and minimum-distance estimation in Stata**  
Doosoo Kim ([doosoonet@gmail.com](mailto:doosoonet@gmail.com))  
Version 0.3.0 (hotfix 6)

**[Read the quick user guide (PDF)](docs/jwselect_quick_user_guide.pdf)**

`jwselect` selects optional preperiod restrictions and estimates cohort and aggregate treatment effects using cohort-time estimates and their full joint covariance. It works after a supported `jwdid` fit or with separately supplied coefficient and covariance matrices.

The command provides:

- **Manual selection:** use a prespecified set of preperiod restrictions, including baseline-only estimation.
- **Wald–BIC selection:** choose restrictions using the pretrend Wald statistic and a BIC penalty.
- **Focused selection:** balance estimated bias and variance for a cohort effect or aggregate, with numerical reselection inference.

Design options include cohort or common baselines, anticipation, treatment-effect horizons, and fixed aggregation weights. The quick guide is the main user manual; Stata's `help jwselect` provides the full option reference.

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

* Select optional preperiod restrictions using Wald–BIC
jwselect, select(wbic) candidates(all) postwindow(0 2)
estat selection
estat design
estimates store JW_WBIC

* Compare with the maintained-baseline estimate
estimates restore JW_SOURCE
jwselect, select(manual) postwindow(0 2)
estimates store JW_BASE
```

**Restore `JW_SOURCE` before every new direct `jwselect` specification.** Each run replaces Stata's active estimation results. A bare `jwselect` replays the current results.

Here, `postwindow(0 2)` averages the first three treatment periods within each cohort and requires those dates to exist for every cohort. The example explicitly searches all preperiod subsets; the package default is `candidates(nested)`. All-subset searches can become large.

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

## Citation and contact

When using the method, please cite the accompanying working paper:

Kim, Doosoo. 2026. *From Pretrend Tests to Preperiod Selection in Difference-in-Differences*. Working paper.

Please also report the software version used. Software citation metadata is provided in [`CITATION.cff`](CITATION.cff).

For questions or bug reports, use this repository's **Issues** tab or email [doosoonet@gmail.com](mailto:doosoonet@gmail.com). Include your Stata version, `which jwselect` output, the command, and a small reproducible example where possible.

## Development and checks

The maintained numerical sources are in [`mata/`](mata/). After editing them, run `python build.py` to regenerate `_jwselect_mata.do`. Both the modules and the generated loader belong in the repository. The guide's editable [LaTeX source](docs/jwselect_quick_user_guide.tex) accompanies the PDF.

From the repository root, run the Python source and independent numerical checks:

```sh
python -m pip install -r requirements-dev.txt
python build.py --check
python tests/static_check.py
python tests/numerical_check.py
python tests/reselect_check.py
python tests/release030_reference.py
python tests/core_optimization_check.py
```

These checks passed when this repository package was prepared on 12 October 2026. They check source structure and independent numerical calculations; they do not execute Stata or establish statistical coverage.

For the Stata release checks, save your work, start a **fresh Stata 17 or later session**, change to the repository root, and run:

```stata
do run_release_tests.do
```

The driver clears data and estimates. Look for the final `jwselect 0.3.0 release tests: PASS` message. Stata was unavailable in the packaging environment, so these runtime checks were not rerun there. The separate `tests/integration_jwdid.do` checks the upstream adapter when its dependencies are installed; it does not install them.
