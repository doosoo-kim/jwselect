# jwselect 0.9.0

Preperiod selection and minimum-distance estimation for difference-in-differences in Stata. Ordinary Wald-BIC (WBIC) with a pointwise normal interval is the default; manual and focused selection are also available.

## Installation

Extract this ZIP, save your work, and open a fresh Stata session. Install from the folder containing `jwselect.pkg`:

```stata
net install jwselect, from("C:/path/to/jwselect-0.9.0") replace
which jwselect
help jwselect
```

The main command requires Stata 16.1 or later. Prepared focused-inference workers and their collector require Stata 17 or later; the repeated-cross-section example also specifies Stata 17. Matrix input and `jwselect_rcs` need no external estimator. Other routes require the supported upstream packages described in the full guide. The bundled native acceleration plugin is for 64-bit Windows; the Mata engine is also available.

After estimation, the current build is reported as `e(build) = "0.9.0-20261104"`. Keep all modules on the same build and recreate prepared inference jobs after upgrading.

## Changes since 0.8.8

- **Repeated-cross-section DR:** New `jwselect_rcs` allows changing composition by transporting source cells to the treated postperiod covariate distribution. Its full covariance includes nuisance estimation and target-distribution uncertainty; pooled and date-stratified sampling are supported.
- **Sun-Abraham population aggregation:** New `jwselect_sa_iw` includes estimated cohort-share uncertainty. It supports the balanced, unweighted, no-covariate, never-treated design accepted by `jwselect_sa`. The direct route continues to use fixed reporting weights.
- **Nonlinear DiD:** The existing full conditional-index option, `restrictions(index)`, becomes the default after supported logit/Poisson fits. Explicit `restrictions(placebo)` remains available as a different restriction choice. A target-specific average-counterfactual extension is not implemented.
- **Panel DR:** The supported `csdid` adapter now requires `pscoretrim(1)`, together with `method(dr) nevertreated base_period(universal) analytical pointwise` and its other documented design conditions.
- **McCloskey inference:** An analytic terminal interval in the same McCloskey family replaces the former Berger-Boos fallback. New diagnostics identify finite terminal and exceptional whole-line returns. `mcfallback(error)` remains available.
- **Documentation:** Both PDF guides and installed help describe the current interfaces, assumptions, covariance conventions, numerical qualifications, and saved results.

Old do-files using `mcfallback(bb)`, `isreps()`, `bbcells()`, or `bbtol()` must be updated. Nonlinear results can change when `restrictions()` is omitted. Existing selection and inference engines otherwise retain their behavior for matched existing inputs.

The new RCS and population-IW sources supply one prespecified scalar target to manual, WBIC, or aggregate-focused selection. Set the target design in the source command; changing it requires refitting. RCS uses finite-dimensional OLS/logit nuisance models, independent representative samples, overlap, treatment-unaffected covariates, and target-standardized identifying restrictions. It does not support active trimming, overlapping panels, or survey weights.

## Guides and examples

- [Quick user guide](docs/jwselect_quick_user_guide.pdf): installation, first fits, interpretation, and short workflows.
- [Full user guide](docs/jwselect_user_guide.pdf): supported models, assumptions, option reference, saved results, inference, parallel workflows, and troubleshooting.
- In Stata: `help jwselect`, `help jwselect_rcs`, `help jwselect_sa`, `help jwselect_sa_iw`, `help jwselect_worker`, and `help jwselect_collect`.

Editable LaTeX sources for both guides are included in `docs/`.

After installation, run examples from the extracted package folder:

```stata
do examples/matrix_demo.do
do examples/rcs_composition_demo.do
```

Run examples in a separate Stata session without unsaved data; the RCS example clears the session. The PDFs and examples remain in the extracted folder rather than being installed by `net install`.

## Inference and validation scope

Read interval endpoints from `e(ci)` and inspect `estat inference`; for McCloskey terminal returns, also inspect `e(mc_terminal_diagnostics)`. Uniform coverage requires regular sampling inputs and outward numerical error control. Current McCloskey calculations use double precision and numerical padding, without a formal interval-arithmetic certificate. Focused numerical search does not certify global extrema. See the full guide for each inference method's interpretation.

Targeted implementation checks cover the nonlinear default, untrimmed panel DR, RCS transport and covariance, population-share uncertainty, scalar-input validation, and package installation. The manuscript's historical simulations and empirical tables have not been regenerated under this release.

## Citation

Kim, Doosoo. 2026. *From Pretrend Tests to Preperiod Selection in Difference-in-Differences*. Working paper. Please report the software version and build used. Machine-readable metadata are in [CITATION.cff](CITATION.cff).

Contact: doosoonet@gmail.com. File checksums are provided in `SHA256SUMS.txt`.
