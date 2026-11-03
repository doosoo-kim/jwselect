# jwselect 0.8.8

`jwselect` selects optional preperiod restrictions and uses minimum distance to estimate cohort and aggregate effects in difference-in-differences designs. Ordinary Wald-BIC (WBIC) with a pointwise normal interval is the default. Manual and focused selection are also available.

Version 0.8.8 is prepared as the first public release in the 0.8 series. Its public predecessor is 0.7.0.

## Install in Stata

Extract the ZIP. In a fresh Stata session, use the path to the folder containing `jwselect.pkg`. If upgrading from 0.7.0, first run `ado uninstall jwselect` to remove its registered files.

```stata
net install jwselect, from("C:/path/to/jwselect-0.8.8") replace
which jwselect
help jwselect
```

The main command supports Stata 16.1 or later. Prepared focused-inference workers and their collector require Stata 17 or later. No external packages are needed for matrix input. For fitted-model input, install the upstream estimator you use. See the guides for supported designs and covariance requirements.

## First example

Run this from the extracted package folder after installation:

```stata
do examples/matrix_demo.do
```

The example supplies a cohort-time vector, its full sampling covariance, and a map of cohort and calendar dates. It fits a baseline model, runs ordinary WBIC, and compares two adjusted intervals using the same data. For a supported `jwdid`, `csdid`, or `jwselect_sa` fit, store the source estimates and restore them before each comparison:

```stata
estimates store sourcefit
jwselect, select(wbic)
estat selection
estimates restore sourcefit
jwselect, inference(wbic) wbicmethod(bootstrap)
estat inference
estimates restore sourcefit
jwselect, inference(wbic) wbicmethod(mc)
estat inference
```

Both adjusted methods keep the WBIC point estimate. `bootstrap` uses centered Gaussian reselection for pointwise inference. `mc` calibrates a marginal interval against a rejection-risk bound in the fitted Gaussian experiment. If calibration fails, McCloskey returns a labeled, protected Berger–Boos fallback by default. Use `mcfallback(error)` to stop instead. Read the endpoints from `e(ci)` or `estat inference`.

## Guides

- [Quick user guide](docs/jwselect_quick_user_guide.pdf) for installation, a first fit, and reading the results.
- [User guide](docs/jwselect_user_guide.pdf) for estimator inputs, selection, inference, and prepared focused-inference jobs.
- In Stata, use `help jwselect`, `help jwselect_sa`, `help jwselect_worker`, and `help jwselect_collect`.

## Citation

Kim, Doosoo. 2026. *From Pretrend Tests to Preperiod Selection in Difference-in-Differences*. Working paper. Please report the software version used. [Software citation](CITATION.cff)
