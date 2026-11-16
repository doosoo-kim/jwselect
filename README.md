# jwselect 0.9.1

`jwselect` selects optional preperiod restrictions and uses minimum distance to estimate cohort and aggregate effects in difference-in-differences designs. Ordinary Wald-BIC (WBIC) with a pointwise normal interval is the default. Manual and focused selection are also available.

Version 0.9.1 aligns the documentation with the manuscript dated November 16, 2026, and adds `estat compare` for matched baseline, own-cohort MD, and joint MD estimates. Cohort-wise Wald-BIC selection is followed by joint MD using the full covariance; the new comparison holds those selections and reporting targets fixed. The release retains full conditional-index blocks as the nonlinear default, doubly robust repeated-cross-section inputs with changing composition, Sun–Abraham population-share uncertainty, and the existing inference engines. See [revision notes](REVISION_NOTES.md) for scope and validation limits.

## Install in Stata

Extract the ZIP. Save your work and open a fresh Stata session, then use the path to the folder containing `jwselect.pkg`. This release reports `e(version) = "0.9.1"` and `e(build) = "0.9.1-20261116"` and rejects a mismatched loaded module. Recreate prepared inference jobs after upgrading. Confirm the path printed by `which jwselect` so an older installed copy does not take precedence.

```stata
net install jwselect, from("C:/path/to/jwselect-0.9.1") replace
which jwselect
help jwselect
```

The main command supports Stata 16.1 or later. Prepared focused-inference workers and their collector require Stata 17 or later. Matrix input and `jwselect_rcs` need no external estimator. Install the upstream estimator used by the other source routes. The installed help files describe current supported designs and covariance requirements.

## New source routes and restrictions

- After supported nonlinear `jwdid` logit or Poisson fits, `restrictions(index)` is the default. Each optional date selects the complete intercept-and-covariate coefficient block; Wald-BIC penalizes its identified rank. The outcome-scale ATT and index blocks use their full joint covariance. Explicit `restrictions(placebo)` retains the legacy average fitted-placebo option. A zero average placebo is a different restriction and does not imply the full conditional-index null. A weaker target-specific average formulation would require its own identifying map and estimation argument; it is not implemented by either option.
- The panel `csdid` adapter requires `pscoretrim(1)`, together with `method(dr) nevertreated base_period(universal) analytical pointwise`. The upstream default `.995` is rejected: active trimming is outside this untrimmed-ATT implementation.
- `jwselect_rcs` fits one DR repeated-cross-section target by standardizing cohort-date outcome means to the covariate distribution of its treated postperiod population. These standardized means are called transported means. Composition may change across dates. `sampling(pooled)` is the default for independent pooled observations; `sampling(stratified)` uses fixed date sample sizes and date-centered full influence functions. This route requires representative independent samples, overlap and the maintained target-standardized parallel-trends restriction. No anticipation must also hold after target standardization of every included treated preperiod; a zero mean only under the source cell's own covariate distribution need not suffice. It does not support overlapping panels, survey weights or arbitrary outcome-selective sampling.
- After `jwselect_sa`, `jwselect_sa_iw, postwindow(a b)` constructs a population IW target and joint covariance including estimated cohort-share uncertainty. Direct `jwselect` after `jwselect_sa` continues to use fixed aggregation weights.

The RCS and population-IW sources supply one prespecified scalar target. `jwselect` supports manual, WBIC and aggregate-focused selection for that target. Set its population, baseline and event window in the source command; it cannot be changed afterward. These routes do not average separately selected cohort effects.

```stata
jwselect_rcs y x1 x2, tvar(year) gvar(first_treat) cohort(4) target(5) baseline(3) sampling(stratified)
estimates store RCS_SOURCE
jwselect, select(wbic)

jwselect_sa y, ivar(id) tvar(year) gvar(first_treat)
jwselect_sa_iw, postwindow(0 2)
estimates store IW_SOURCE
jwselect, select(focused) focus(aggregate) inference(none)
```

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
estat compare
estimates restore sourcefit
jwselect, inference(wbic) wbicmethod(bootstrap)
estat inference
estimates restore sourcefit
jwselect, inference(wbic) wbicmethod(mc)
estat inference
```

Both adjusted methods keep the WBIC point estimate. `bootstrap` implements centered Gaussian bootstrap reselection: it repeats selection and estimation in each draw from the fitted Gaussian experiment and uses a centered, studentized error statistic. This numerically approximates a pointwise asymptotic construction. No higher-order accuracy improvement is asserted. Fixed simulation budgets leave quantile simulation error; increase draws and check sensitivity. `mc` calibrates a marginal interval against a rejection-probability bound in the fitted Gaussian experiment. If no search level passes, the default `mcfallback(terminal)` returns the analytic terminal interval in the same McCloskey family. `estat inference` identifies calibrated and terminal returns. `mcfallback(error)` stops when the search fails. Read the endpoints from `e(ci)` and inspect `e(mc_terminal_diagnostics)` for exceptional whole-line returns. The uniform coverage result requires the stated regularity conditions and rejection-probability bounds that remain conservative after numerical approximation and rounding. Current calculations use double precision with numerical padding, without a formal interval-arithmetic certificate. The retained `mcmode(certified)` compatibility label does not remove that limitation.

## Matched baseline, own-cohort, and joint MD comparison

After ordinary cohort-time or matrix-input estimation, run:

```stata
estat compare
matrix comparison = r(table)
estat compare, cells
```

The default table uses the reported cohort and aggregate targets in `e(output_map)`; `output(#)` selects one row, and `cells` instead compares the included cohort-time cells. `cells` and `output()` cannot be combined. The baseline uses no optional restrictions. Own-cohort MD applies each cohort's portion of the selected mask to its own effects; joint MD applies the full mask using the complete cross-cohort covariance. Both use identical selected restrictions and target weights, without reselection. Aggregate own-cohort variances retain cross-cohort covariance. Scalar population sources (`jwselect_rcs` and `jwselect_sa_iw`) do not supply this decomposition and are rejected.

The table contains `baseline se_baseline own se_own joint se_joint own_minus_base joint_minus_base joint_minus_own`. All standard errors are calculated by treating the selected restrictions as fixed. The table supplies neither selection-adjusted intervals nor tests of the differences; consult the original `e(ci)` and `e(pvalue)` only for the estimator and inference method that produced them. Point estimates can move in either direction, and asymptotic efficiency comparisons do not guarantee finite-sample improvement. `estat compare` preserves the active estimation results. See `help jwselect` for saved covariance, loading, mask, and target-map matrices and for focused-selection interpretation.

## Explicit all-subsets and centered-bootstrap workflow

The defaults are `candidates(nested) maxcand(4096)` and ordinary normal intervals treating the selected restrictions as fixed. The November 16, 2026 manuscript's empirical workflow instead uses all subsets within cohorts and explicit centered Gaussian reselection. The following template preserves one prespecified event window across comparison and inference; choose its dates for the application before inspecting results:

```stata
estimates restore sourcefit
jwselect, select(wbic) wbicscope(cohort) candidates(all) maxcand(16384) ///
    postwindow(0 2) inference(wbic) wbicmethod(bootstrap) reps(49999) seed(12345)
estat compare
matrix comparison = r(table)
estat inference
```

`maxcand(16384)` accommodates all subsets of up to 14 eligible optional dates per cohort, including the empty set. Larger families require a larger explicit limit; exceeding the limit stops rather than silently truncating the search. `wbicscope(cohort)` governs selection, while the resulting WBIC estimates remain joint MD. This is a workflow template, not a substitute for the manuscript's application data, maintained baseline, anticipation setting, target definitions, and covariance construction.

Focused inference retains numerical limitations: multistart search and independent simulation at discovered points do not certify global nuisance extrema. Reporting a target different from the focused selection target additionally needs the cross-target regularity/tie condition explained in the full guide; the software does not verify it automatically.

## Guides

Both user guides have been updated for this release. The full guide retains advanced explanations, saved-result schemas, numerical/parallel workflows, troubleshooting, and the legacy all-candidate envelope and cross-target extension. These advanced materials remain distinct from the default workflow.

- The [quick user guide](docs/jwselect_quick_user_guide.pdf) and [full user guide](docs/jwselect_user_guide.pdf) describe version 0.9.1, including the matched MD comparison, current defaults, RCS standardization, population-share SA targets, trimming contract, and inference qualifications. Editable LaTeX sources are included in `docs/`.
- Use `help jwselect`, `help jwselect_rcs`, `help jwselect_sa`, `help jwselect_sa_iw`, `help jwselect_worker`, and `help jwselect_collect`.

## Citation

Kim, Doosoo. 2026. *From Pretrend Tests to Preperiod Selection in Difference-in-Differences*. Working paper. Please report the software version used. [Software citation](CITATION.cff)
